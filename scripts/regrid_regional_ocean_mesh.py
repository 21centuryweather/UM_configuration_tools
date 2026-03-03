#!/usr/bin/env python

import xarray as xr
import numpy as np
from argparse import ArgumentParser
import os
import sys
import metpy

if 'ESMFMKFILE' not in os.environ:
    os.environ['ESMFMKFILE'] = '/opt/conda/analysis3-25.10/lib/esmf.mk'

import esmpy
esmpy.Manager(debug=True)

def add_options(parser):
    """
    Process the command-line arguments
    """
    parser.description = 'Regrid a MOM6 ocean mesh onto a regular grid with the same bounds.\n The user specifies the output resolution'

    parser.add_argument('-f','--input-ocean-mesh-file', 
        dest='ocean_file', 
        type=str, 
        help="The rMOM6 ESMF mesh file, e.g. access-rom3-ESMFmesh.nc")
    parser.add_argument('-r','--resolution', 
        dest='atm_res', 
        type=str, 
        help="The desired atmospheric resolution in degrees, e.g. 0.025")


def load_ocn_data(ocn_mesh_fp):
    """
    Load the rMOM6 ocean mesh file and return at ESMPY mesh object
    """

    mesh_ds = xr.load_dataset(
        ocn_mesh_fp
    )

    bounds = get_bounds(mesh_ds)

    ocn_mask = mesh_ds.elementMask.data

    # We need to logically invert this because MOM6 and the UM tread land/ocean 
    # points differently. i.e. MOM6 has ocean points as '1' (True), the UM has them
    # as '0' (False)
    #ocn_mask = np.logical_not(ocn_mask).astype('int')

    ocn_mesh = esmpy.api.mesh.Mesh(
        filename=ocn_mesh_fp, 
        filetype=esmpy.api.constants.FileFormat.ESMFMESH
    )

    return ocn_mesh, ocn_mask, bounds


def get_bounds(ocn_ds):
    """
    Return the bounds of the the ocean mesh dataset
    """
    coords = ocn_ds.nodeCoords.data

    # Find the min and max 
    lons = coords[:,0]
    lats = coords[:,1]

    bounds = { 'lat_min' : lats.min(), 
               'lat_max' : lats.max(), 
               'lon_min' : lons.min(), 
               'lon_max' : lons.max() }

    return bounds


def build_grid_cons(bounds, res):
    """
    Build an esmpy (ESMF) grid object in memory using a modification 
    of Kieran's method at https://gist.github.com/kieranricardo/eb98f76235255efff800d28c2442e5c3
    """
    
    lon0, lon1 = bounds['lon_min'],bounds['lon_max']
    lat0, lat1 = bounds['lat_min'],bounds['lat_max']

    # Let's enforce the distance b/w the bounds is a multiple of the input resolution
    # to prevent rounding errors

    delta_lat = np.ptp([bounds['lat_min'],bounds['lat_max']])
    delta_lon = np.ptp([bounds['lat_min'],bounds['lat_max']])

    nlat = int(np.round(delta_lat/res))
    nlon = int(np.round(delta_lon/res))

    max_index = np.array([nlon, nlat])
    grid = esmpy.Grid(max_index, 
                    coord_sys=esmpy.CoordSys.SPH_DEG,
                    staggerloc=[esmpy.StaggerLoc.CORNER])


    # Get the coordinate arrays for the cell vertices 
    grid_lon_v = grid.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    grid_lat_v = grid.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)

    # Define the coordinates for the vertices 
    # The coordinate array size is now (nlat+1) x (nlog+1)
    lat_vertices = np.linspace(lat0, lat1, max_index[0] + 1)
    lon_vertices = np.linspace(lon0, lon1, max_index[1] + 1)

    # Assign coordinates to the grid object
    for j in range(max_index[0] + 1):
        for i in range(max_index[1] + 1):
            grid_lat_c[j, i] = lat_vertices[j]
            grid_lon_c[j, i] = lon_vertices[i]

    return grid, lat_vertices, lon_vertices


def regrid_cons(grid,
                ocn_mesh,
                ocn_mask,
                bounds, 
                lat_vertices,
                lon_vertices):
    """
    Regrid my conservative grid object
    """
    dst_field_conserve = esmpy.Field(grid, name='source_data_conserve')

    # Define the output data on the edges of the cells (which is N_y x N_x)
    dst_lon_edges = ocn_mesh.coords[0][0]
    dst_lat_edges = ocn_mesh.coords[0][1]
    dst_field_conserve.data[:] = 100.0 + np.outer(src_lat_centers, np.cos(np.deg2rad(src_lon_centers)))

    # Invert the mask
    ocn_mask = np.logical_not(ocn_mask).astype('int')

    ocn_field = esmpy.Field(ocn_mesh, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    ocn_field.data[:] = ocn_mask

    ocn_to_atm_cons = esmpy.api.regrid.Regrid(
      ocn_field, 
      src_field_conserve, 
      unmapped_action=esmpy.api.constants.UnmappedAction.IGNORE,
      regrid_method=esmpy.api.constants.RegridMethod.CONSERVE,
      norm_type=esmpy.api.constants.NormType.DSTAREA, 
      factors=True
    )

    new_ocn_frac = src_field_conserve.data.reshape((len(src_lat_centers),len(src_lon_centers)))

    # Output to array
    lat_coord = xr.DataArray(
      dims=['latitude'],
      coords=dict(latitude=src_lat_centers),
      data=src_lat_centers,
      attrs=dict(standard_name='latitude', units='degrees_north')
    )

    lon_coord = xr.DataArray(
      dims=['longitude'],
      coords=dict(longitude=src_lon_centers),
      data=src_lon_centers,
      attrs=dict(standard_name='longitude', units='degrees_east')
     )

    da = xr.DataArray(
    data=new_ocn_frac,
    dims=["latitude", "longitude"],
    attrs=dict(um_stash_source='m01s00i505',
               grid_mapping='geog_cs',
               earth_radius=6371229.0)
    )
    
    # --- 2. Create the GeogCS Grid Mapping Variable ---
    # This variable is crucial for defining the coordinate system (CS).
    # It must be a scalar variable (value does not matter)
    grid_mapping_var = xr.Variable(
        dims=(),
        data=0,
        attrs={
            # Standard CF attribute indicating a geographical coordinate system
            'grid_mapping_name': 'latitude_longitude',

            # CF attributes defining the spherical earth datum (radius)
            # 6371229.0 m is the radius for Met Office standard sphere
            'earth_radius': 6371229.0,
        }
    )

    ds = xr.Dataset(
            data_vars=dict(
                land_binary_mask=da,
                geog_cs=grid_mapping_var
                ),
            coords=dict(
                latitude=lat_coord,
                longitude=lon_coord,
            ),
         attrs=dict(Conventions='CF-1.7')
         )


    return ds

    
if __name__ == "__main__":

    parser = ArgumentParser()
    add_options(parser)

    if len(sys.argv) == 1:
        print (" INFO : No input data supplied.")
        print (" INFO : usage: python %prog [options].")
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    ocean_file = args.ocean_file
    atm_res = float(args.atm_res)

    ocn_mesh, ocn_mask, bounds = load_ocn_data(ocean_file)

    atm_grid_cons, lats, lons = build_grid_cons(bounds, 
                                                atm_res)

    con_ds = regrid_cons(atm_grid_cons, 
                         ocn_mesh,
                         ocn_mask,
                         bounds, 
                         lats, 
                         lons)

    out_fp = 'target_grid.nc'
    con_ds.to_netcdf(out_fp)