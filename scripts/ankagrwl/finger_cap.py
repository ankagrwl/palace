import gmsh
import sys
import os

current_dir = os.path.dirname(__file__)
print(current_dir)

refinement: int = 1
trace_width_μm: float  = 30.0
gap_width_μm: float = 18.0
separation_width_μm: float = 200.0
ground_width_μm: float = 800.0
substrate_height_μm: float = 500.0
length_μm: float = 4000.0
filename: str = 'finger_cap.msh'
verbose: int = 2

kernel = gmsh.model.occ
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", verbose)

# Add model
if "cpw" in gmsh.model.list():
    gmsh.model.setCurrent("cpw")
    gmsh.model.remove()
gmsh.model.add("cpw")


num_fingers = 4
finger_l = 50 * 25.4
finger_w = 10 * 25.4
finger_gap = 10 * 25.4
l = finger_gap + finger_l + 2* finger_w
w = (2* num_fingers) * finger_w + (2*num_fingers-1) * finger_gap

def finger_block(starting_coord, num_fingers, finger_l, finger_w, finger_gap):
    """Eventually may want to pass dx and dy as an input"""
    dx, dy = starting_coord

    left_g1 = kernel.addRectangle(dx, dy, 0.0, finger_w, w)
    dy += w
    dx += finger_w
    temp_left = []
    temp_left.append(left_g1)

    for ii in range(num_fingers):
        t = kernel.addRectangle(dx, dy, 0.0, finger_l, -finger_w)
        dy -= 2*(finger_gap + finger_w)
        temp_left.append(t)


    dx += (finger_l + finger_gap)
    dy = substrate_w/2 - w/2

    right_g1 = kernel.addRectangle(dx, dy, 0.0, finger_w, w)
    temp_right = []
    temp_right.append(right_g1)

    for ii in range(num_fingers):
        t = kernel.addRectangle(dx, dy, 0.0, -finger_l, finger_w)
        dy += 2*(finger_gap + finger_w)
        temp_right.append(t)

    return temp_left, temp_right


def meander(starting_coord, height, gap_width, trace_width, num_periods):

    dx, dy = starting_coord

    period_length = 2 * gap_width + 3 * trace_width

    temp = []

    for ii in range(num_periods):
        x = dx + ii * period_length
        y = dy 
        t = kernel.addRectangle(x, y, 0, trace_width, height + trace_width)
        temp.append(t)

        x += trace_width
        y += height

        t=kernel.addRectangle(x, y, 0, gap_width, trace_width)
        temp.append(t)

        x += gap_width
        y += trace_width

        t=kernel.addRectangle(x, y, 0, trace_width, -2 * (trace_width + height))
        temp.append(t)

        x += trace_width
        y -= 2 * (height + trace_width)

        t=kernel.addRectangle(x, y, 0, gap_width, trace_width)
        temp.append(t)

        x += gap_width
        t=kernel.addRectangle(x, y, 0, trace_width, trace_width + height)
        temp.append(t)

        y += trace_width + height
        print(x, y)

    return temp


# Substrate/chip dimensions
# substrate_w = w + 2 * l
# substrate_l = 4 * l
substrate_w = 7000
substrate_l = 12000

substrate_h = 20 * 25.4

sep_dz = substrate_h
sep_dy = sep_dz

x0 = substrate_l/2 -l/2 
y0 = substrate_w/2 - w/2

left, right = finger_block([x0, y0], num_fingers, finger_l, finger_w, finger_gap)
left = [(2, m) for m in left]
right = [(2, m) for m in right]

height = 20 * 25.4
gap_width = 10 * 25.4
trace_width = 10 * 25.4
num_periods = 3

meander_tags = meander([x0 + l + l/2, substrate_w/2], height, gap_width, trace_width, num_periods)
meander_tags = [(2, m) for m in meander_tags]

print(left)
print(right)
print(meander_tags)

# Mesh parameters
l_trace = 1.0 * finger_w * 2**(-refinement)
l_farfield = 1.0 * substrate_h * 2**(-refinement)

# Substrate
substrate = kernel.addBox(
    0.0, 
    0.0, 
    -substrate_h,
    substrate_l, 
    substrate_w, 
    substrate_h
)

# Exterior box
domain = kernel.addBox(
    -sep_dy,
    -sep_dy,
    -substrate_h-sep_dz,
    substrate_l + 2.0 * sep_dy,
    substrate_w + 2.0 * sep_dy,
    substrate_h + 2.0 * sep_dz
)

_, domain_boundary = kernel.getSurfaceLoops(domain)
domain_boundary = domain_boundary[0]
print(domain_boundary)

# Ports
port_w = 2 * finger_w
dx = 0.0
dy = substrate_w/2- port_w/2

p1 = kernel.addRectangle(0.0, dy, 0.0, substrate_l/2 -l/2, port_w)
dx += (substrate_l/2 + l/2)
p2 = kernel.addRectangle(dx, dy, 0.0, substrate_l/2 -l/2, port_w)

#Gap
g = kernel.addRectangle(0, 0, 0, substrate_l, substrate_w)
gap = kernel.cut([(2, g)], left + right + [(2, p1)] + [(2, p2)], removeObject=True, removeTool=False)[0]
print(gap)

# Embedding
f = lambda x: (x[0]==2 or x[0]==3)

geom_dimtags = list(filter(f,  kernel.getEntities()))

_, geom_map = kernel.fragment(geom_dimtags, [])

kernel.synchronize()

print(geom_dimtags)
print(geom_map)

# Add physical groups
si_domain = [t[-1] for t in geom_map[geom_dimtags.index((3, substrate))]][0]

air_domain = [t[-1] for t in geom_map[geom_dimtags.index((3, domain))] if t!=(3, si_domain)][0]

# print(domain, si_domain, air_domain)

si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "si")
air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")


metal = []
temp_metal = left + right + meander_tags

temp = [geom_map[i] for i, x in enumerate(geom_dimtags) if x in temp_metal]
metal = [int(m[0][-1]) for m in temp]

metal_group = gmsh.model.addPhysicalGroup(2, metal, -1, "metal")


port1 = [g[-1] for g in geom_map[geom_dimtags.index((2, p1))]]
port2 = [g[-1] for g in geom_map[geom_dimtags.index((2, p2))]]


port1_group = gmsh.model.addPhysicalGroup(2, port1, -1, "port1")
port2_group = gmsh.model.addPhysicalGroup(2, port2, -1, "port2")


temp = [geom_map[i] for i, x in enumerate(geom_dimtags) if x in gap]
gap = [m[0][-1] for m in temp]

gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")

farfield = [geom_map[ii][0][-1] for ii, f in enumerate(geom_dimtags) if f[0]==2 and f[1] in domain_boundary]

print(farfield)

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

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()

if verbose > 0:
    print("\nFinished generating mesh. Physical group tags:")
    print("Si domain: ", si_domain_group)
    print("Air domain: ", air_domain_group)
    print("Farfield boundaries: ", farfield_group)
    print("Metal boundaries: ", metal_group)
    print("Trace negative boundaries: ", gap_group)

    print("\nMultielement lumped ports:")
    print("Port 1: ", port1_group)
    print("Port 2: ", port2_group)
    print()
