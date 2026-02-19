#!/usr/bin/env nemesis
"""Generate a tri or quad mesh of a subduction zone vertical profile using Gmsh, making
use of the built-in geometry engine.

Points have been projected from longitude/latitude into a local
transverse Mercator projection. PyLith uses the Proj.4 library
for geographic projections. The proj parameters are:

+proj=tmerc +datum=WGS84 +lon_0=142.0 +lat_0=38.0 +k=0.9996

so that the local origin is at a longitude of 142.0 degrees (WGS84)
and a latitude of 38.0 degrees (WGS84).

Run `generate_gmsh.py --help` to see the command line options.
"""
import gmsh
import numpy as np
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh, group_exclude)

class App(GenerateMesh):
    """
    Application for generating the mesh.
    """
    X_WEST = -1500.0e+3
    X_EAST = 2500.0e+3
    Y_BOT = -1000.0e+3
    Y_MOHO = -40.0e+3

    DX_FAULT = 2.0e+3
    DX_BIAS = 1.05

    TOPO_POINTS = (
        (   X_EAST, 100.0),
        (839.1e+3, 100.0),
        (739.1e+3, 1000.0),
        (639.1e+3, 1000.0),
        (539.1e+3, 1000.0),
        (439.1e+3, 100.0),
        (351.2e+3, 0.0),
        (263.4e+3,    -300.0),
        (175.6e+3,  -700.0),
        ( 87.7e+3,    -1500.0),
        (   0.0e+3,  -3000.0),
        (  -87.7e+3, -3000.0),
        ( -165.6e+3, -3000.0),
        ( -263.4e+3, -3000.0),
        ( -351.2e+3, -3000.0),
        ( -439.1e+3, -3000.0),
        (-539.1e+3, -3000.0),
        (-639.1e+3, -3000.0),
        (-739.1e+3, -3000.0),
        (-839.1e+3, -3000.0),
        (   X_WEST, -3000.0),
    )


    TOPO_EAST = 0
    TOPO_TRENCH = 10
    TOPO_WEST = len(TOPO_POINTS)-1


    # Top of slab from Slab 2 - see MATLAB script (Extract_Slab2.m)
    ## CHANGE INPUT OF SLAB ACCORDINGLY HERE
    # Flipping and reading the slab top points
    SLABTOP_POINTS = np.flipud(np.loadtxt(fname='./Geometries/Polynomial/slab_top_geom_poly4_v1.txt')) 

    # Flipping and reading the slab top points
    SLABBOT_POINTS = np.flipud(np.loadtxt(fname='./Geometries/Polynomial/slab_bot_geom_poly4_v1.txt')) 

    # Get the difference between the trench and the slab (to be adjusted later on)
    differenceX = TOPO_POINTS[TOPO_TRENCH][0] - SLABTOP_POINTS[-1,0]

    # Adjusting the position of the slab closer to the trench
    SLABTOP_POINTS[:,0] += differenceX
    SLABBOT_POINTS[:,0] += differenceX

    SLABTOP_EAST = 0
    SLABTOP_COSEISMIC = 73 # this cannot be higher (shallower) than moho, index number must be smaller than moho 
    SLABTOP_MOHO = 74    

    SLABTOP_WEST = len(SLABTOP_POINTS)-1 

    Y_BOT = SLABBOT_POINTS[0,1]

    # Add a bottom plate of the oceanic crust here (connects from bottom plate)
    SLABBOT_ADD = (
        (-205.6e+3, Y_MOHO),
        (  X_WEST, Y_MOHO),
    )

    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options
        # with the default cell type `tri` matching the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "mesh_tri_poly_test_v1.msh"

    def create_geometry(self):
        """Create geometry.
        """
        # Create curve for topography/bathymetry
        pts_topo = []
        for x, y in self.TOPO_POINTS:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_topo.append(pt)
        # c_topo = gmsh.model.geo.add_spline(pts_topo)
        c_topo = gmsh.model.geo.add_polyline(pts_topo)
        p_topo_west = pts_topo[self.TOPO_WEST]
        p_topo_east = pts_topo[self.TOPO_EAST]
        p_topo_trench = pts_topo[self.TOPO_TRENCH]

        pts_slabtop = []
        for x, y in self.SLABTOP_POINTS:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_slabtop.append(pt)
        c_slabtop = gmsh.model.geo.add_bspline(pts_slabtop + [p_topo_trench])
        self.p_slabtop_east = pts_slabtop[self.SLABTOP_EAST]
        self.p_slabtop_west = pts_slabtop[self.SLABTOP_WEST]
        self.p_slabtop_coseismic = pts_slabtop[self.SLABTOP_COSEISMIC]
        p_slabtop_moho = pts_slabtop[self.SLABTOP_MOHO]

        # # Create b-spline curve for the bottom of the slab.
        # We translate points to the east. A better approach would be to
        # move points normal to the slab to preserve uniform thickness.

        pts_slabbot = []
        for x, y in self.SLABBOT_POINTS[0:-1]:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_slabbot.append(pt)
        for x, y in self.SLABBOT_ADD:
            pt = gmsh.model.geo.add_point(x, y, 0.0)
            pts_slabbot.append(pt)
        self.c_slabbot = gmsh.model.geo.add_bspline(pts_slabbot)
        self.p_slabbot_west = pts_slabbot[-1]
        self.p_slabbot_east = pts_slabbot[0]

        # Create curve for bottom edge of slab
        self.c_slab_end = gmsh.model.geo.add_polyline([pts_slabtop[0], pts_slabbot[0]])
        # c_slab_end = gmsh.model.geo.add_polyline([self.p_slabtop_east, self.p_slabbot_east])

        # # Create top of mantle (Ignore mantle wedge if don't want to force con_moho)
        p_moho_east = gmsh.model.geo.add_point(self.X_EAST, self.Y_MOHO, 0.0)
        self.c_conmoho = gmsh.model.geo.add_line(p_moho_east, p_slabtop_moho) # original
        # p_mantlewedge = gmsh.model.geo.add_point(self.TOPO_POINTS[8][0], self.Y_MOHO, 0.0) # for adding a point in con_moho (wedge size testing)
        # self.c_conmoho = gmsh.model.geo.add_polyline([p_moho_east, p_mantlewedge, p_slabtop_moho]) # for adding a point in con_moho (wedge size testing)

        # # Create lateral edges and bottom
        p_bot_west = gmsh.model.geo.add_point(self.X_WEST, self.Y_BOT, 0.0)
        p_bot_east = gmsh.model.geo.add_point(self.X_EAST, self.Y_BOT, 0.0)
        c_east = gmsh.model.geo.add_polyline([p_bot_east, p_moho_east, p_topo_east])
        c_west = gmsh.model.geo.add_polyline([p_bot_west, self.p_slabbot_west, p_topo_west])

        # self.c_bot = gmsh.model.geo.add_polyline([p_bot_west, p_bot_east])
        self.c_bot_oce = gmsh.model.geo.add_polyline([p_bot_west, pts_slabbot[0]])
        self.c_bot_con = gmsh.model.geo.add_polyline([pts_slabtop[0], p_bot_east])

        
        # # ----------------------------------------------------------------------
        # # Split curves to form bounding curves for each material
        # # Constructing the entire boundary curves as splines and then breaking them into
        # # pieces preserves C1 continuity in the curves.
        # # ----------------------------------------------------------------------

        curves = gmsh.model.geo.split_curve(c_topo, [p_topo_trench])
        self.c_topo_west = curves[1]
        self.c_topo_east = curves[0]

        curves = gmsh.model.geo.split_curve(c_slabtop, [self.p_slabtop_coseismic, p_slabtop_moho])
        self.c_slabtop_mantle_lower = curves[0]
        self.c_slabtop_mantle_upper = curves[1]
        self.c_slabtop_crust = curves[2]

        curves = gmsh.model.geo.split_curve(c_west, [self.p_slabbot_west])
        self.c_west_mantle = curves[0]
        self.c_west_crust = curves[1]

        curves = gmsh.model.geo.split_curve(c_east, [p_moho_east])
        self.c_east_mantle = curves[0]
        self.c_east_crust = curves[1]

        # ----------------------------------------------------------------------
        # Create surfaces from bounding curves
        # They are mix in direction of each of the loops
        # ----------------------------------------------------------------------
        # below is counter clockwise (continental crust)
        loop = gmsh.model.geo.add_curve_loop([
            self.c_east_crust,
            -self.c_conmoho,
            -self.c_slabtop_crust,
            self.c_topo_east,            
            ])
        self.s_concrust = gmsh.model.geo.add_plane_surface([loop])

        # Below is counter clockwise (oceanic crust)
        loop = gmsh.model.geo.add_curve_loop([
            -self.c_slab_end,
            -self.c_slabbot,
            -self.c_west_crust,
            self.c_topo_west,
            self.c_slabtop_crust,
            self.c_slabtop_mantle_upper,
            self.c_slabtop_mantle_lower,
            ])
        self.s_oceancrust = gmsh.model.geo.add_plane_surface([loop])

        # Below is continental mantle counter clockwise 
        loop = gmsh.model.geo.add_curve_loop([
            self.c_east_mantle,
            self.c_bot_con,
            -self.c_slabtop_mantle_lower,
            -self.c_slabtop_mantle_upper,
            self.c_conmoho,
            ])
        self.s_conmantle = gmsh.model.geo.add_plane_surface([loop])

        # Below is oceanic mantle counter clockwise
        loop = gmsh.model.geo.add_curve_loop([
            self.c_bot_oce,
            -self.c_west_mantle,
            self.c_slabbot,
            ])
        self.s_ocemantle = gmsh.model.geo.add_plane_surface([loop])

        gmsh.model.geo.synchronize()


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented.
        """
        # Create materials matching surfaces.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_concrust]),
            MaterialGroup(tag=2, entities=[self.s_oceancrust]),
            MaterialGroup(tag=3, entities=[self.s_conmantle]),
            MaterialGroup(tag=4, entities=[self.s_ocemantle]),
        )
        for material in materials:
            material.create_physical_group()

        # Create physical groups for the boundaries and the fault.
        vertex_groups = (
            VertexGroup(name="groundsurf", tag=10, dim=1, entities=[self.c_topo_east, self.c_topo_west]),
            VertexGroup(name="bndry_east", tag=11, dim=1, entities=[self.c_east_mantle, self.c_east_crust]),
            VertexGroup(name="bndry_west_crust", tag=12, dim=1, entities=[self.c_west_crust]),
            VertexGroup(name="bndry_west_mantle_tmp", tag=113, dim=1, entities=[self.c_west_mantle]),
            VertexGroup(name="bndry_bot_con", tag=14, dim=1, entities=[self.c_bot_con]),
            VertexGroup(name="bndry_bot_oce", tag=15, dim=1, entities=[self.c_bot_oce]),
            VertexGroup(name="bndry_bot_slab_end", tag=16, dim=1, entities=[self.c_slab_end]),
            VertexGroup(name="fault_coseismic", tag=20, dim=1, entities=[self.c_slabtop_mantle_upper, self.c_slabtop_crust]),
            VertexGroup(name="fault_coseismic_edge", tag=30, dim=0, entities=[self.p_slabtop_coseismic]),
            VertexGroup(name="fault_slabtop", tag=21, dim=1, entities=[self.c_slabtop_mantle_lower, self.c_slabtop_mantle_upper, self.c_slabtop_crust]),
            VertexGroup(name="fault_slabtop_edge", tag=31, dim=0, entities=[self.p_slabtop_east]),
            VertexGroup(name="fault_slabbot", tag=22, dim=1, entities=[self.c_slabbot]),
            VertexGroup(name="fault_slabbot_edge", tag=32, dim=0, entities=[self.p_slabbot_east]),
        )


        for group in vertex_groups:
            group.create_physical_group()
        # group_exclude("bndry_east_mantle_tmp", "fault_slabbot", new_name="bndry_east_mantle", new_tag=13)
        group_exclude("bndry_west_mantle_tmp", "fault_slabbot", new_name="bndry_west_mantle", new_tag=13)

    def generate_mesh(self, cell):
        """Generate the mesh.
        """
        # Set discretization size with geometric progression from distance to the fault.
        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the fault.
        fault_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", [self.c_slabtop_mantle_upper, self.c_slabtop_crust,self.c_topo_east])
        # gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", [self.c_slabtop_mantle_lower, self.c_slabtop_mantle_upper, self.c_slabtop_crust])

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
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")


if __name__ == "__main__":
    App().main()


# End of file
