# ------------------------------------------------------------------------------
# 2024/02/12: Ankur
# Python file to generate the mesh file for a capacitively coupled lambda/4 resonator
# ------------------------------------------------------------------------------
import gmsh
import sys
import os
import numpy as np
from itertools import chain

current_dir = os.path.dirname(__file__)
print(current_dir)
"""
 generate_coplanar_waveguide_lumped_mesh(;
        refinement::Integer = 0,
        trace_width_μm::Real = 30.0,
        gap_width_μm::Real = 18.0,
        separation_width_μm::Real = 200.0,
        ground_width_μm::Real = 800.0,
        substrate_height_μm::Real = 500.0,
        substrate_length_μm::Real = 4000.0,
        filename::AbstractString,
        verbose::Integer = 1,
    )
Generate a mesh for the coplanar waveguide with lumped ports using Gmsh

# Arguments

  - refinement - measure of how many elements to include, 0 is least
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - substrate_length_μm - length of the waveguides, in μm
  - filename - the filename to use for the generated mesh
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
"""

def generate_coplanar_waveguide_lumped_mesh(
    refinement: int = 1,
    trace_width_μm: float  = 30.0,
    gap_width_μm: float = 18.0,
    substrate_height_μm: float = 500.0,
    substrate_length_μm: float = 4000.0,
    filename: str = 'test.msh',
    verbose: int = 2
    ):

    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "cpw" in gmsh.model.list():
        gmsh.model.setCurrent("cpw")
        gmsh.model.remove()
    gmsh.model.add("cpw")
    
    sep_dz = 2000.0
    sep_dy = 0.5 * sep_dz

    # Mesh parameters
    l_trace = 1.5 * trace_width_μm * 2**(-refinement)
    l_farfield = 1.0 * substrate_height_μm * 2**(-refinement)

    substrate_width_μm = 3000
    substrate_length_μm = 7000

    ground_width_μm = (substrate_width_μm - 2 * gap_width_μm - trace_width_μm) // 2

    # Chip pattern
    dy = 0.0
    g1 = kernel.addRectangle(0.0, dy, 0.0, substrate_length_μm, ground_width_μm)
    dy += ground_width_μm
    n1 = kernel.addRectangle(0.0, dy, 0.0, substrate_length_μm, gap_width_μm)
    dy += gap_width_μm
    t1 = kernel.addRectangle(0.0, dy, 0.0, substrate_length_μm, trace_width_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, substrate_length_μm, gap_width_μm)
    dy += gap_width_μm
    g2 = kernel.addRectangle(0.0, dy, 0.0, substrate_length_μm, ground_width_μm)
    dy += ground_width_μm


    # Making the coupling capacitor symmetrically on each side
    resonance_freq = 10e9 #GHz
    quarter_length_μm = 1e6*(3e8/resonance_freq)/np.sqrt(np.mean(12+1)) //2 #in μm

    print(quarter_length_μm)

    cap_gap_μm = 5
    x0 = 0
    y0 = ground_width_μm + gap_width_μm
    cap_length_μm = substrate_length_μm/2 - quarter_length_μm/2 - cap_gap_μm

    # c1 = kernel.addRectangle(x0, y0, 0.0, cap_length_μm , trace_width_μm)
    x0 += cap_length_μm
    n3 = kernel.addRectangle(x0, y0, 0.0, cap_gap_μm, trace_width_μm)
    x0 += cap_gap_μm
    # t1 = kernel.addRectangle(x0, y0, 0.0, quarter_length_μm, trace_width_μm)
    x0 += quarter_length_μm
    n4 = kernel.addRectangle(x0, y0, 0.0, cap_gap_μm, trace_width_μm)
    x0 += cap_gap_μm
    # c2 = kernel.addRectangle(x0, y0, 0.0, cap_length_μm, trace_width_μm)

    print(n3, n4)
    t1 = kernel.cut([(2, t1)], [(2, n3), (2, n4)], removeObject=False, removeTool=False)[0]
    print(t1, n3, n4)

    # Substrate
    substrate = kernel.addBox(
        0.0, 
        0.0, 
        -substrate_height_μm,
        substrate_length_μm, 
        dy, 
        substrate_height_μm
    )

    # Exterior box
    domain = kernel.addBox(
        -0.5 * sep_dy,
        -sep_dy,
        -sep_dz,
        substrate_length_μm + sep_dy,
        dy + 2.0 * sep_dy,
        2.0 * sep_dz
    )

    _, domain_boundary = kernel.getSurfaceLoops(domain)
    domain_boundary = domain_boundary[0]

    # Ports
    dy = ground_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2a = kernel.addRectangle(substrate_length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2b = kernel.addRectangle(substrate_length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)

    # Embedding
    f = lambda x: (x[0]==2 or x[0]==3)

    geom_dimtags = list(filter(f,  kernel.getEntities()))

    _, geom_map = kernel.fragment(geom_dimtags, [])

    kernel.synchronize()


    # Add physical groups
    si_domain = [t[-1] for t in geom_map[geom_dimtags.index((3, substrate))]][0]

    air_domain = [t[-1] for t in geom_map[geom_dimtags.index((3, domain))] if t!=(3, si_domain)][0]

    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "si")
    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")

    m_list = [(2, g1), (2, g2)] + t1
    metal = []

    temp = [geom_map[i] for i, x in enumerate(geom_dimtags) if x in m_list]
    metal = [int(m[0][-1]) for m in temp]

    metal_group = gmsh.model.addPhysicalGroup(2, metal, -1, "metal")

    port1a = [g[-1] for g in geom_map[geom_dimtags.index((2, p1a))]]
    port2a = [g[-1] for g in geom_map[geom_dimtags.index((2, p2a))]]
    port1b = [g[-1] for g in geom_map[geom_dimtags.index((2, p1b))]]
    port2b = [g[-1] for g in geom_map[geom_dimtags.index((2, p2b))]]

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port2a_group = gmsh.model.addPhysicalGroup(2, port2a, -1, "port2a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")
    port2b_group = gmsh.model.addPhysicalGroup(2, port2b, -1, "port2b")

    ports_list = list(chain.from_iterable([port1a, port1b, port2a, port2b]))

    g_list = [(2, n1), (2, n2), (2, n3), (2, n4)]

    temp = [geom_map[i] for i, x in enumerate(geom_dimtags) if x in g_list]
    temp = list(chain.from_iterable(temp))

    gap = [m[-1] for m in temp if m[-1] not in ports_list]

    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")

    farfield = [geom_map[ii][0][-1] for ii, f in enumerate(geom_dimtags) if f[0]==2 and f[1] in domain_boundary]

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    f = lambda x: (x[0]==0)

    gap_points = list(filter(f, gmsh.model.getBoundary([(2, z) for z in gap], False, True, True)))

    f = lambda x: (x[0]==1)

    gap_curves = [int(l[-1]) for l in list(filter(f, gmsh.model.getBoundary([(2, z) for z in gap], False, False, False)))]

    gmsh.model.mesh.setSize(gap_points, l_trace)

    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", gap_curves)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", gap)
    gmsh.model.mesh.field.setNumber(1, "Power", 1.0)
    gmsh.model.mesh.field.setNumber(1, "DistMax", sep_dz)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", l_farfield)

    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "CurvesList", gap_curves)
    gmsh.model.mesh.field.setNumber(2, "Sampling", 10)

    gmsh.model.mesh.field.add("Threshold", 3)
    gmsh.model.mesh.field.setNumber(3, "InField", 2)
    gmsh.model.mesh.field.setNumber(3, "SizeMin", l_trace)
    gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(3, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(3, "DistMax", sep_dz)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [1, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(1)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)

    filename = current_dir + '/mesh/' + filename

    print(filename)

    gmsh.write(filename)

    """
    Add the point on the geometry where we want to compute the field and makse sure that it is meshed
    and then compute the field in post-processing
    """
    # Launch the GUI to see the results:
    if '-nopopup' not in sys.argv:
        gmsh.fltk.run()

    gmsh.finalize()

    # # Print some information
    if verbose > 0:
        print("\nFinished generating mesh. Physical group tags:")
        print("Si domain: ", si_domain_group)
        print("Air domain: ", air_domain_group)
        print("Farfield boundaries: ", farfield_group)
        print("Metal boundaries: ", metal_group)
        print("Trace negative boundaries: ", gap_group)

        print("\nMultielement lumped ports:")
        print("Port 1: ", port1a_group, ", ", port1b_group)
        print("Port 2: ", port2a_group, ", ", port2b_group)
        print()


filename = "lambda_resonator.msh"
generate_coplanar_waveguide_lumped_mesh(filename=filename, refinement=1)