#!/usr/bin/env python
import esmpy
import xarray as xr
import numpy as np
from argparse import ArgumentParser
import os
import sys

def add_options(parser):
    """
    Process the command-line arguments
    """
    parser.description = 'Plot the memory usage of the Unified Model I/O server'

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

def create_atmos_mask(bounds, resolution):
    """
    Create a atmospheric land/sea dataarray within the given bounds at the desired resolution
    The mesh must have an even number of points.
    """
    epsilon = 1e-9 

    lat_coords = np.arange(bounds['lat_min'], bounds['lat_max'] + resolution - epsilon, resolution)

    N_LAT_NODES = len(lat_coords)

    if N_LAT_NODES % 2 != 0:
        print (f'ERROR : Number of atmospheric latitude points is {N_LAT_NODES}')
        print ('ERROR : Change resolution to ensure an even number of points within your ocean mesh bounds' )
        sys.exit(1)

    lon_coords = np.arange(bounds['lon_min'], bounds['lon_max'] + resolution - epsilon, resolution)

    N_LON_NODES = len(lon_coords)

    if N_LON_NODES % 2 != 0:
        print (f'ERROR : Number of atmospheric longitude points is {N_LON_NODES}')
        print ('ERROR : Change resolution to ensure an even number of points within your ocean mesh bounds' )
        sys.exit(1)

    # Build 
    atm_grid = build_grid(bounds, atm_res)
    # Temp array - will regrid the ocean land/sea mask in due course
    

    da =  xr.DataArray(
        data,
        coords={
            'lat': lat_coords,
            'lon': lon_coords
        },
        dims=['lat', 'lon'],
        name='land_binary_mask',
        attrs={
            'description': 'MOM6 land/sea mask regridded to atmospheric resolution'
        }
    )

    return da


def build_grid(bounds, nlat):
    """
    Build an esmpy (ESMF) grid object in memory using Kieran's method at https://gist.github.com/kieranricardo/eb98f76235255efff800d28c2442e5c3
    """
    nlon = nlat
    lon0, lon1 = bounds['lon_min'],bounds['lon_max']
    lat0, lat1 = bounds['lat_min'],bounds['lat_max']

    max_index = np.array([nlon, nlat])
    grid = esmpy.Grid(max_index, num_peri_dims=1, staggerloc=[esmpy.StaggerLoc.CENTER])

    dlat = (lat1 - lat0) / nlat
    lat = (np.arange(nlat) * dlat) + lat0 + 0.5 * dlat

    dlon = (lon1 - lon0) / nlon
    lon = (np.arange(nlon) * dlon) + lon0 + 0.5 * dlon

    grid.get_coords(0)[:] = lon[:, None]
    grid.get_coords(1)[:] = lat[None, :]

    # corner 
    grid.add_coords([esmpy.StaggerLoc.CORNER])

    grid_lon_corner = grid.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    grid_lat_corner = grid.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)

    grid_lon_corner[:] = lon[:, None] - 0.5 * dlon
    grid_lat_corner[:, :-1] = lat[None, :] - 0.5 * dlat
    grid_lat_corner[:, -1] = 90.0

    return grid, lat, lon



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

    out_fp = 'dummy.nc'
    nlat=26
    nlon=nlat
    ocn_mesh, ocn_mask, bounds = load_ocn_data(ocean_file)
    atm_grid, lats, lons = build_grid(bounds, nlat)

    atm_field = esmpy.Field(atm_grid, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    atm_field.data[:] = 0.0    

    ocn_field = esmpy.Field(ocn_mesh, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    ocn_field.data[:] = ocn_mask

    ocn_to_atm_cons = esmpy.api.regrid.Regrid(
      ocn_field, atm_field, 
      unmapped_action=esmpy.api.constants.UnmappedAction.IGNORE,
      regrid_method=esmpy.api.constants.RegridMethod.CONSERVE,
      norm_type=esmpy.api.constants.NormType.DSTAREA, factors=True
    )

    new_ocn_frac = atm_field.data.reshape((nlat, nlon))

    new_land_frac = 1.0 - new_ocn_frac

    new_land_frac[new_land_frac < 0.01] = 0.0
    new_land_frac[new_land_frac > 1.0] = 1.0

    lat_coord = xr.DataArray(
      dims=['lat'],
      coords=dict(lat=lats),
      data=lats,
      attrs=dict(standard_name='latitude', units='degrees_north')
    )

    lon_coord = xr.DataArray(
      dims=['lon'],
      coords=dict(lon=lons),
      data=lons,
      attrs=dict(standard_name='longitude', units='degrees_east')
     )

    da = xr.DataArray(
    data=new_land_frac,
    dims=["lat", "lon"],
    coords=dict(
        lat=lats,
        lon=lons,
    ),
    attrs=dict(um_stash_source='m01s00i505')
    )

    ds = xr.Dataset(dict(landfrac=da))

    ds.to_netcdf(out_fp)