import gmsh
import sys
import os

current_dir = os.path.dirname(__file__)
print(current_dir)

refinement: int = 1
inductor_width_μm: float  = 30.0
inductor_gap_μm: float = 18.0
separation_width_μm: float = 200.0
ground_width_μm: float = 800.0
substrate_height_μm: float = 500.0
length_μm: float = 4000.0
filename: str = 'finger_cap.msh'
verbose: int = 2


class LCOscillator:
    def __init__(self, x0, y0):
        self.x = x0
        self.y = y0
        self.kernel = gmsh.model.occ
        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", verbose)

        # Add model
        if "cpw" in gmsh.model.list():
            gmsh.model.setCurrent("cpw")
        gmsh.model.remove()
        gmsh.model.add("cpw")

    def add_substrate(self, substrate_length, substrate_width, substrate_height):
        """
        Adds a substrate to the model.

        Args:
            substrate_length: length of the substrate
            substrate_width: width of the substrate
            substrate_height: height of the substrate

        Returns:
            substrate: gmsh.model.occ.addBox() object representing the substrate

        Example:
            substrate = add_substrate(100, 50, 10)
            Creates a substrate with initial vertex (0, 0, -10)
            and dimensions (100, 50, 10) along the positive (x, y, z) direction
        
        """
        self.substrate_length = substrate_length
        self.substrate_width = substrate_width
        self.substrate_height = substrate_height

        # Substrate
        substrate = self.kernel.addBox(
            0.0, 
            0.0, 
            -self.substrate_height,
            self.substrate_length, 
            self.substrate_width, 
            self.substrate_height
        )
        return substrate

    def add_domain(self, pad_x, pad_y, pad_z):
        """"
        Adds a domain (bounding box) to the simulation

        Args:
            pad_x: padding along the x-axis
            pad_y: padding along the y-axis
            pad_z: padding along the z-axis

        Returns:
            domain: tag representing the domain
            domain_boundary: tag representing the domain boundary

        Example:
            domain, domain_boundary = add_domain(10, 10, 10)
            Creates a domain with initial certex (0, 0, -10)
            and dimensions (100, 50, 10) along the positive (x, y, z) direction

        """

        # Exterior box
        domain = self.kernel.addBox(
            -pad_x,
            -pad_y,
            -self.substrate_height-pad_z,
            self.substrate_length + 2.0 * pad_x,
            self.substrate_width + 2.0 * pad_y,
            self.substrate_height + 2.0 * pad_z
        )

        _, domain_boundary = self.kernel.getSurfaceLoops(domain)
        domain_boundary = domain_boundary[0]

        return domain, domain_boundary

    def add_finger_capacitor(self, finger_l, finger_w, finger_gap, num_fingers):
        self.finger_l = finger_l
        self.finger_w = finger_w
        self.finger_gap = finger_gap
        self.num_fingers = num_fingers

        self.cap_block_length = finger_gap + finger_l + 2* finger_w
        self.cap_block_width = (2* num_fingers) * finger_w + (2*num_fingers-1) * finger_gap

        self.left_fingers, self.right_fingers = self.finger_cap_n((self.x, self.y - self.cap_block_width/2),
                                                                   self.num_fingers, 
                                                                   self.finger_l, 
                                                                   self.finger_w, 
                                                                   self.finger_gap)
        
    def finger_cap_n(self, 
                     starting_coord,
                     num_fingers: int = None, 
                     finger_l: float = None,
                     finger_w: float = None,
                     finger_gap: float = None):
        """
        Creates a set of interdigitated (finger) capacitors

        Args:
            starting_coord: (x, y) coordinates of the bottom left corner of the block
            num_fingers: number of fingers in the block
            finger_l: length of each finger
            finger_w: width of each finger
            finger_gap: gap between each finger

        Returns:
            left_fingers: list of left fingers
            right_fingers: list of right fingers

        Example:
            finger_block((0, 0), 2, 10, 1, 2)

            Creates a block with 2 fingers of length 10, width 1, and gap 2.
            The block bottom left corner is at (0, 0).
            It updates the global (x, y) --> (right most coorinate along x, center of the capcitor block along y)
            coordinates of the layout for further usage
            and returns the left and right fingers as lists.        
        
        """
        dx, dy = starting_coord

        left_g1 = self.kernel.addRectangle(dx, dy, 0.0, finger_w, self.cap_block_width)
        dy += self.cap_block_width
        dx += finger_w
        temp_left = []
        temp_left.append(left_g1)

        for _ in range(num_fingers):
            t = self.kernel.addRectangle(dx, dy, 0.0, finger_l, -finger_w)
            dy -= 2*(finger_gap + finger_w)
            temp_left.append(t)

        dx += (finger_l + finger_gap)
        dy -= (finger_gap + finger_w)

        right_g1 = self.kernel.addRectangle(dx, dy, 0.0, finger_w, self.cap_block_width)
        temp_right = []
        temp_right.append(right_g1)

        for _ in range(num_fingers):
            t = self.kernel.addRectangle(dx, dy, 0.0, -finger_l, finger_w)
            dy += 2*(finger_gap + finger_w)
            temp_right.append(t)

        dy += (finger_gap + 2*finger_w)

        self.x = dx
        self.y = dy - self.cap_block_width/2

        return temp_left, temp_right

    def add_inductor(self, inductor_w, inductor_gap, num_periods):
        self.inductor_w = inductor_w
        self.inductor_gap = inductor_gap
        self.num_periods = num_periods

        self.inductor = self.meander((self.x, self.y),
                                     self.inductor_w,
                                     self.inductor_gap,
                                     self.num_periods, 
                                    )

    def meander(self, starting_coord, inductor_width, inductor_gap, num_periods):
        """
        Creates a meander (inductor) and add leads on each side of the meander

        Args:
            starting_coord = (x, y) coordinates of the starting point of the meander
            inductor_width = width of the inductor (along the y-direction) 
            inductor_gap = thickness of the metal trace
            num_periods = number of periods in the meander

        Returns:
            meander: list of gmsh.model.occ.addRectangle() objects representing the meander

        Example:
            meander((0, 0), 3, 1, 2)

            Creates a meander with 3 width, 1 gap, and 2 periods. Add leads of length 3 and width 1
            on each side.
            The meander center left starts at (0, 0).
            It updates the global (x, y) --> (right most coorinate along x, center of the inductor along y)
            coordinates of the layout for further usage
            and returns the meander as a list of gmsh.model.occ.addRectangle() objects
        
        
        """

        dx, dy = starting_coord

        period_length = 2 * inductor_width + 3 * inductor_gap

        temp = []
        dy = dy - inductor_gap/2

        # Left connecting section
        t = self.kernel.addRectangle(dx, dy, 0, inductor_width, inductor_gap)
        temp.append(t)

        dx += inductor_width
        dy += inductor_gap/2

        for ii in range(num_periods):
            if ii==0:
                t = self.kernel.addRectangle(dx, dy + inductor_gap/2, 0, inductor_gap, inductor_width)
            else:
                t = self.kernel.addRectangle(dx, dy, 0, inductor_gap, inductor_gap + inductor_width)

            temp.append(t)

            dx += inductor_gap
            dy += inductor_width

            t = self.kernel.addRectangle(dx, dy, 0, inductor_width, inductor_gap)
            temp.append(t)

            dy += inductor_gap
            dx += inductor_width

            t = self.kernel.addRectangle(dx, dy, 0, inductor_gap, -2 * (inductor_width + inductor_gap))
            temp.append(t)

            dx += inductor_gap
            dy -= 2 * (inductor_width + inductor_gap)

            t = self.kernel.addRectangle(dx, dy, 0, inductor_width, inductor_gap)
            temp.append(t)

            dx += inductor_width

            if ii == (num_periods-1):
                t = self.kernel.addRectangle(dx, dy, 0, inductor_width, 0.5 * inductor_gap + inductor_width)
            else:    
                t = self.kernel.addRectangle(dx, dy, 0, inductor_width, inductor_gap + inductor_width)
            temp.append(t)

            dy += inductor_width + inductor_gap

        # Right connecting section
        dy -= inductor_gap/2
        t = self.kernel.addRectangle(dx, dy, 0, inductor_gap + inductor_width, inductor_gap)
        temp.append(t)

        self.dx = dx + inductor_gap + inductor_width
        self.dy = dy + inductor_gap/2

        return temp

    def add_port(self, starting_coord, port_w, port_l):
        """
        Add a port to the chip

        Args:
            port_w: width of the port
            port_l: length of the port
        
        Returns:
            port: gmsh.model.occ.addRectangle() object representing the port

        Example:
            add_port((0, 0), 2, 10)

            Creates a port with 2 width, 10 length.
            The port bottom left corner starts at (0, 0).
            It updates the global (x, y) --> (right most coorinate along x, center of the port along y)
            coordinates of the layout for further usage
            and returns the port as a gmsh.model.occ.addRectangle() object
        """

        p = self.kernel.addRectangle(starting_coord[0], starting_coord[1] - port_w/2, 0, port_l, port_w)
        self.x = starting_coord[0] + port_l
        self.y = starting_coord[1] + port_w/2
        return p
    

"""
substrate_w = 7000
substrate_l = 7000
substrate_h = 20 * 25.4

pad_z = substrate_h
pad_y = 0.1* pad_z
pad_x = pad_y


# Finger capacitor dimensions
num_fingers = 4
finger_l = 50 * 25.4
finger_w = 10 * 25.4
finger_gap = 10 * 25.4
l = finger_gap + finger_l + 2* finger_w
w = (2* num_fingers) * finger_w + (2*num_fingers-1) * finger_gap

# Meander dimensions
height = 20 * 25.4
inductor_gap = 10 * 25.4
inductor_width = 5 * 25.4
num_periods = 3

meander_length = num_periods * (2 * inductor_gap + 3 * inductor_width) + (2 * inductor_gap - 2 * inductor_width)

print(f'Length of inductor section is {meander_length} μm')
# Starting from left and adding each component

# Port 1
port_w = 2 * finger_w

dx = 0.0
dy = substrate_w/2- port_w/2

port_l = substrate_l/2 -l/2 - meander_length/2

p1 = kernel.addRectangle(0.0, dy, 0.0, port_l, port_w)

# Finger capacitors
dx += port_l
dy = substrate_w/2 - w/2

left, right = finger_block([dx, dy], num_fingers, finger_l, finger_w, finger_gap)
left = [(2, m) for m in left]
right = [(2, m) for m in right]


dx += l
dy = substrate_w/2
# Inductor
meander_tags = meander([dx, dy], height, inductor_gap, inductor_width, num_periods)
if not len(meander_tags) ==0:
    meander_tags = [(2, m) for m in meander_tags]
else:
    meander_tags = [()]

# Port 2
dx = substrate_l - port_l
dy -= port_w/2
p2 = kernel.addRectangle(dx, dy, 0.0, port_l, port_w)

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

#Gap
g = kernel.addRectangle(0, 0, 0, substrate_l, substrate_w)
gap = kernel.cut([(2, g)], left + right + meander_tags + [(2, p1)] + [(2, p2)], removeObject=True, removeTool=False)[0]

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

metal = []
temp_metal = left + right + meander_tags

temp = [geom_map[i] for i, x in enumerate(geom_dimtags) if x in temp_metal]
metal = [int(m[0][-1]) for m in temp]

metal_group = gmsh.model.addPhysicalGroup(2, metal, -1, "metal")
print(metal_group)

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
"""