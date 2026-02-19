#!/usr/bin/env nemesis
"""Generate a tri or quad mesh of a reverse fault with a splay fault using Gmsh, making
use of the built-in geometry engine.

We use the `gmsh_utils` module provided with PyLith. This module has helper functions
and classes for marking materials and boundaries compatible with PyLith. We also use the
`GenerateMesh` class from the module as a base class to our local `App` class. The
`GenerateMesh` class handles processing of command line options, initialization, and
finalizing the mesh.

Run `generate_gmsh.py --help` to see the command line options.

Run `generate_gmsh.py --write` to generate the mesh.

DOING LINEAR SLAB GEOMETRY - GOOD VERSION FOR ANY DEGREE SLAB - JHC on 7-Jul-2025

"""


###
# Changing the xpos to allow for slab to escape
###

# Import the math module for trigonometric functions
import math

# Import Gmsh Python interface
import gmsh

# Import the gmsh_utils Python module supplied with PyLith.
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh, group_exclude)

class App(GenerateMesh):
    """
    Application used to generate the mesh using Gmsh.

    App uses `GenerateMesh` from `gmsh_utils` for common functionality that we avoid
    duplicating in each of our examples.

    Domain is 200km by 100km.

    -100.0 km <= x <= 100.0 km
    -100.0 km <= y <= 0.0 km

    """
    # Define some constants that determine the geometry of the domain.
    DOMAIN_X = 20000.0e+3
    DOMAIN_Y = 10000.0e+3 
    LD = 40.0e+3 # locking depth (give positive value) [meters]
    FAULT_DIP = 10.0 # degrees (don't exceed yneg boundary, not too steep)
    SLAB_THK = 40e+3 # slab thickness [meters]

    DX_FAULT = 2000
    DX_BIAS = 1.08 

    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options
        # with the default cell type `tri` matching the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "./Mesh/mesh_tri_slab_10deg_40thc_final_v1.msh"

    def create_geometry(self):
        """Create geometry.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set local variables for domain size and corner of the domain.
        lx = self.DOMAIN_X
        ly = self.DOMAIN_Y
        x1 = -0.5 * lx
        y1 = -ly

        # Create points for domain perimeter (p1 = btm left, p2 = btm right, p3 = top right, p4 = top left)
        p1 = gmsh.model.geo.add_point(x1, y1, 0.0)
        p2 = gmsh.model.geo.add_point(x1+lx, y1, 0.0)
        p3 = gmsh.model.geo.add_point(x1+lx, y1+ly, 0.0)
        p4 = gmsh.model.geo.add_point(x1, y1+ly, 0.0)

        # Create points for main fault
        # w = self.FAULT_WIDTH
        fault_dip = self.FAULT_DIP / 180.0 * math.pi

        LDy = -self.LD # locking depth
        LDx = LDy / -math.tan(fault_dip) # locking depth
        fsx = 0.5 * lx # fault end
        fsy = -math.tan(fault_dip)*fsx # fault end
        p_SLABTOP_trench = gmsh.model.geo.add_point(0, 0, 0.0) # point of surface on fault
        p_SLABTOP_lock = gmsh.model.geo.add_point(LDx, LDy, 0.0) # point of locking depth on fault
        p_SLABTOP_east = gmsh.model.geo.add_point(fsx, fsy, 0.0)# point of fault connected to the side boundary

        # Create points for continental crust
        CRUST_east = gmsh.model.geo.add_point(fsx, 0-self.SLAB_THK, 0.0)
        # Use slabtop and slabtop lock points to connect the continental crust

        # Create points for slabbot fault
        p_SLABBOT_east = gmsh.model.geo.add_point(fsx, fsy-self.SLAB_THK, 0.0) # point exit east
        p_SLABBOT_trench = gmsh.model.geo.add_point(0, 0-self.SLAB_THK, 0.0) # point below trrench
        p_SLABBOT_west = gmsh.model.geo.add_point(x1, 0-self.SLAB_THK, 0.0) # outer (west)

        # Create curves. We store the curve tag as a data member
        # so that we can refer to them later.
        # We traverse the curves in a counter clock-wise direction.
        # If the curve is in the opposite direction, we use the negative tag.
       
        self.c_yneg = gmsh.model.geo.add_line(p1, p2)
        self.c_xpos_om = gmsh.model.geo.add_line(p2, p_SLABBOT_east) # oceanic mantle
        self.c_xpos_slab = gmsh.model.geo.add_line(p_SLABBOT_east, p_SLABTOP_east) # slab 
        self.c_xpos_cm = gmsh.model.geo.add_line(p_SLABTOP_east, CRUST_east) # continental mantle
        self.c_xpos_cont = gmsh.model.geo.add_line(CRUST_east, p3) # continental mantle
        self.c_ypos_cont = gmsh.model.geo.add_line(p3, p_SLABTOP_trench)
        # self.c_ypos_cm = gmsh.model.geo.add_line(CRUST_east, p_SLABTOP_lock) # continental mantle or bottom of the continental crust
        self.c_ypos_slab = gmsh.model.geo.add_line(p_SLABTOP_trench, p4)
        self.c_xneg_slab = gmsh.model.geo.add_line(p4, p_SLABBOT_west)
        self.c_xneg_om = gmsh.model.geo.add_line(p_SLABBOT_west, p1)

        # Slab top
        self.c_fault_u = gmsh.model.geo.add_line(p_SLABTOP_lock,p_SLABTOP_trench)
        self.c_fault_l = gmsh.model.geo.add_line(p_SLABTOP_east,p_SLABTOP_lock)
        self.p_fault_end = p_SLABTOP_east
        self.p_fault_lock = p_SLABTOP_lock
        self.p_fault_trench = p_SLABTOP_trench

        # Slab bottom
        self.c_slabbot_u = gmsh.model.geo.add_line(p_SLABBOT_trench,p_SLABBOT_west)
        self.c_slabbot_l = gmsh.model.geo.add_line(p_SLABBOT_east,p_SLABBOT_trench)
        self.p_slabbot_end = p_SLABBOT_east

        ### Curve loops ###
        # Continent
        # c0 = gmsh.model.geo.add_curve_loop([self.c_xpos_cont, self.c_ypos_cont, -self.c_fault_u, -self.c_ypos_cm]) # footwall loop
        # self.s_concrust = gmsh.model.geo.add_plane_surface([c0])        
        # Slab
        c0 = gmsh.model.geo.add_curve_loop([self.c_xpos_slab, self.c_fault_l, self.c_fault_u, self.c_ypos_slab, self.c_xneg_slab, -self.c_slabbot_l,-self.c_slabbot_u]) # footwall loop
        self.s_oceancrust = gmsh.model.geo.add_plane_surface([c0])
        # Oceanic mantle
        c0 = gmsh.model.geo.add_curve_loop([self.c_xpos_om, self.c_slabbot_l, self.c_slabbot_u, self.c_xneg_om, self.c_yneg]) # footwall loop
        self.s_ocemantle = gmsh.model.geo.add_plane_surface([c0])
        # Continental mantle
        # c0 = gmsh.model.geo.add_curve_loop([self.c_xpos_cm,  self.c_ypos_cm, -self.c_fault_l]) # footwall loop
        c0 = gmsh.model.geo.add_curve_loop([self.c_xpos_cm, self.c_xpos_cont, self.c_ypos_cont, -self.c_fault_u, -self.c_fault_l]) # footwall loop
        self.s_conmantle = gmsh.model.geo.add_plane_surface([c0])

        gmsh.model.geo.synchronize()


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Create three materials (slab, crust, and wedge).
        # We use the `MaterialGroup` data class defined in `gmsh_utils.`
        # The tag argument specifies the integer tag for the physical group.
        # The entities argument specifies the array of surfaces for the material.
        materials = (
            # MaterialGroup(tag=1, entities=[self.s_concrust]),
            MaterialGroup(tag=2, entities=[self.s_oceancrust]),
            MaterialGroup(tag=3, entities=[self.s_conmantle]),
            MaterialGroup(tag=4, entities=[self.s_ocemantle]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the faults.
        # We use the `VertexGroup` data class defined in `gmsh_utils`.
        # The name and tag specify the name and tag assigned to the physical group.
        # The dimension and entities specify the geometric entities to include in the physical
        # group.
        vertex_groups = (
            VertexGroup(name="boundary_xneg_tmp", tag=110, dim=1, entities=[self.c_xneg_om]),
            VertexGroup(name="boundary_xneg_slab", tag=14, dim=1, entities=[self.c_xneg_slab]),
            VertexGroup(name="boundary_xpos_top", tag=15, dim=1, entities=[self.c_xpos_cm,self.c_xpos_cont]),
            VertexGroup(name="boundary_xpos_bot", tag=16, dim=1, entities=[self.c_xpos_om]),
            VertexGroup(name="boundary_xpos", tag=11, dim=1, entities=[self.c_xpos_om, self.c_xpos_slab, self.c_xpos_cm,self.c_xpos_cont]),
            # VertexGroup(name="boundary_xpos_fw", tag=14, dim=1, entities=[self.c_xpos_fw, self.c_xpos_hw]),
            VertexGroup(name="boundary_yneg", tag=12, dim=1, entities=[self.c_yneg]),
            VertexGroup(name="boundary_ypos", tag=13, dim=1, entities=[self.c_ypos_cont, self.c_ypos_slab]),
            VertexGroup(name="groundsurf", tag=100, dim=1, entities=[self.c_ypos_cont]),
            VertexGroup(name="fault_slabtop", tag=30, dim=1, entities=[self.c_fault_l, self.c_fault_u]),
            # VertexGroup(name="fault_slabtop_end", tag=31, dim=0, entities=[self.p_fault_end]),
            VertexGroup(name="fault_slabtop_lower", tag=20, dim=1, entities=[self.c_fault_l]),
            VertexGroup(name="fault_slabtop_lower_end", tag=21, dim=0, entities=[self.p_fault_end]),
            VertexGroup(name="fault_slabbot", tag=22, dim=1, entities=[self.c_slabbot_l,self.c_slabbot_u]),
            VertexGroup(name="fault_slabbot_end", tag=23, dim=0, entities=[self.p_slabbot_end]),
            VertexGroup(name="fault_slabtop_lock", tag=24, dim=1, entities=[self.c_fault_u]),
            VertexGroup(name="fault_slabtop_lock_end", tag=25, dim=0, entities=[self.p_fault_lock]),

        )
        for group in vertex_groups:
            group.create_physical_group()
        # group_exclude("bndry_east_mantle_tmp", "fault_slabbot", new_name="bndry_east_mantle", new_tag=13)
        group_exclude("boundary_xneg_tmp", "fault_slabbot", new_name="boundary_xneg", new_tag=10)

            

    def generate_mesh(self, cell):
        """Generate the mesh.

        This method is abstract in the base class and must be implemented
        in our local App class.
        """
        # Set discretization size with geometric progression from distance to the fault.
        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 1)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the fault.
        fault_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", [self.c_fault_u])
        # gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", self.c_fault_u.tolist() + [self.c_fault_l])
        # gmsh.model.mesh.field.setNumbers(fault_distance, "PointsList", [self.p_fault_lock])

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.
        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(fault_distance, min_dx=self.DX_FAULT, bias=self.DX_BIAS)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        # Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        if cell == "quad":
            # Generate a tri mesh and then recombine cells to form quadrilaterals.
            # We use the Frontal-Delaunay for Quads algorithm.
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.recombine()
            gmsh.model.mesh.generate(2)
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


# If script is called from the command line, run the application.
if __name__ == "__main__":
    App().main()


# End of file
