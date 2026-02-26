#!/usr/bin/env python

import xarray as xr
import numpy as np
from argparse import ArgumentParser
import os
import sys
import metpy

if 'esmpyMKFILE' not in os.environ:
    os.environ['esmpyMKFILE'] = '/opt/conda/analysis3-25.10/lib/esmpy.mk'

import esmpy
esmpy.Manager(debug=True)

def add_options(parser):
    """
    Process the command-line arguments
    """
    parser.description = 'Regrid a MOM6 ocean mesh onto a regular grid with the same bounds.\n The user specifies the output resolution'

    parser.add_argument('-f','--input-ocean-mesh-file', 
        dest='mesh_file', 
        type=str, 
        help="The rMOM6 esmpy mesh file, e.g. access-rom3-esmpymesh.nc")
    parser.add_argument('-r','--resolution', 
        dest='atm_res', 
        type=str, 
        help="The desired atmospheric resolution in degrees, e.g. 0.025")


def create_center_based_grid(src_grid):
    # 1. Extract the center coordinates from the source
    # These will have shape (ni, nj)
    src_lon_centers = src_grid.get_coords(coord_dim=0, staggerloc=esmpy.StaggerLoc.CENTER)
    src_lat_centers = src_grid.get_coords(coord_dim=1, staggerloc=esmpy.StaggerLoc.CENTER)
    
    # 2. Define dimensions for the new grid
    # To treat centers as corners, the new grid 'cells' are the gaps between centers
    ni, nj = src_lon_centers.shape
    new_grid_shape = np.array([ni - 1, nj - 1])
    
    # 3. Create the new Grid object
    # We use the same coordinate system (usually Spherical)
    new_grid = esmpy.Grid(max_index=new_grid_shape, 
                         coord_sys=src_grid.coord_sys, 
                         staggerloc=[esmpy.StaggerLoc.CENTER, esmpy.StaggerLoc.CORNER])
    
    # 4. Assign the original centers to the new grid's CORNERS
    new_lon_corners = new_grid.get_coords(coord_dim=0, staggerloc=esmpy.StaggerLoc.CORNER)
    new_lat_corners = new_grid.get_coords(coord_dim=1, staggerloc=esmpy.StaggerLoc.CORNER)
    
    new_lon_corners[...] = src_lon_centers
    new_lat_corners[...] = src_lat_centers
    
    # 5. (Optional) Calculate new centers for the new grid
    # Often done via simple averaging of the new corners
    # new_grid.add_coords(staggerloc=esmpy.StaggerLoc.CENTER)
    
    return new_grid


def build_grid_cons(bounds, res):
    """
    Build an esmpy (ESMF) grid object in memory using a modification 
    of Kieran's method at https://gist.github.com/kieranricardo/eb98f76235255efff800d28c2442e5c3
    """
    
    lon0, lon1 = bounds['lon_min'],bounds['lon_max']
    lat0, lat1 = bounds['lat_min'],bounds['lat_max']

    # Let's enforce the distance b/w the bounds is a multiple of the input resolution
    # to prevent rounding errors

    # TO DO - we need some robust checking here to ensure
    # the grid bounds are a multiple of the required
    # resolution

    delta_lat = np.ptp([bounds['lat_min'],bounds['lat_max']])
    delta_lon = np.ptp([bounds['lat_min'],bounds['lat_max']])

    nlat = int(np.round(delta_lat/res))
    nlon = int(np.round(delta_lon/res))

    max_index = np.array([nlon, nlat])
    grid = esmpy.Grid(max_index, 
                    coord_sys=esmpy.CoordSys.SPH_DEG,
                    staggerloc=[esmpy.StaggerLoc.CORNER])


    # Get the coordinate arrays for the cell vertices (or corners)
    grid_lon_v = grid.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    grid_lat_v = grid.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)

    # Define the coordinates for the vertices 
    # The coordinate array size is now (nlat+1) x (nlog+1)
    lat_vertices = np.linspace(lat0, lat1, max_index[0] + 1)
    lon_vertices = np.linspace(lon0, lon1, max_index[1] + 1)

    # Assign coordinates to the grid object
    for j in range(max_index[0] + 1):
        for i in range(max_index[1] + 1):
            grid_lat_v[j, i] = lat_vertices[j]
            grid_lon_v[j, i] = lon_vertices[i]

    # Let's create a dual grid where we create add the 
    # co-ordinates of the grid centers
    grid.add_coords(staggerloc=esmpy.StaggerLoc.CENTER)

    grid_lon_v = grid.get_coords(coord_dim=0, staggerloc=esmpy.StaggerLoc.CENTER)
    grid_lat_c = grid.get_coords(coord_dim=1, staggerloc=esmpy.StaggerLoc.CENTER)

    # Define the coordinates for the centers
    lat_centres = np.linspace(lat0+res/2, lat1-res/2, max_index[0])
    lon_centres = np.linspace(lon0+res/2, lon1-res/2, max_index[1])

    # Assign coordinates to the grid object
    for j in range(max_index[0]):
        for i in range(max_index[1]):
            grid_lat_c[j, i] = lat_centres[j]
            grid_lon_c[j, i] = lon_centres[i]

    return grid, lat_vertices, lon_vertices


if __name__ == "__main__":

    parser = ArgumentParser()
    add_options(parser)

    if len(sys.argv) == 1:
        print (" INFO : No input data supplied.")
        print (" INFO : usage: python %prog [options].")
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    MOM6_mesh_file = args.mesh_file
    atm_res = float(args.atm_res)

    MOM6_mesh, MOM6_mask, bounds = load_ocn_data(ocean_file)

    # Create the MOM6 mask field
    MOM6_field = esmpy.Field(MOM6_mesh,meshloc=esmpy.MeshLoc.ELEMENT)

    # Transfer the mask data into the field
    MOM6_field.data[:] = MOM6_mask

    # Create the target UM Grid
    UM_grid, lats, lons = build_grid_cons(bounds, 
                                          atm_res)

    # Create a UM dummy field
    UM_field = esmpy.Field(UM_grid,staggerloc=esmpy.StaggerLoc.CENTER)

    # To use the conservative regridder, we are going to define a dual grid composted of the UM grid's cell centres

    MOM6_to_UM_cons = esmpy.api.regrid.Regrid(
      MOM6_field, 
      UM_field, 
      unmapped_action=esmpy.api.constants.UnmappedAction.IGNORE,
      regrid_method=esmpy.api.constants.RegridMethod.CONSERVE,
      norm_type=esmpy.api.constants.NormType.DSTAREA, 
      factors=True
    )

    # Now we interpolate the UM field (defined at the grid centre) back to the vertices