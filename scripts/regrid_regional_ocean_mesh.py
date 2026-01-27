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


def build_my_grid(bounds, resolution):
    """
    Build an esmpy (ESMF) grid object in memory co-ords which align with the edge of the MOM6 mesh

    """
    epsilon = 1e-9 

    lat = np.arange(bounds['lat_min'], bounds['lat_max'] + resolution - epsilon, resolution)

    dlat = resolution   
    dlon = dlat

    nlat = len(lat)

    if nlat % 2 != 0:
        print (f'ERROR : Number of atmospheric latitude points is {nlat}')
        print ('ERROR : Change resolution to ensure an even number of points within your ocean mesh bounds' )
        sys.exit(1)

    lon = np.arange(bounds['lon_min'], bounds['lon_max'] + resolution - epsilon, resolution)

    nlon = len(lon)

    if nlon % 2 != 0:
        print (f'ERROR : Number of atmospheric longitude points is {nlon}')
        print ('ERROR : Change resolution to ensure an even number of points within your ocean mesh bounds' )
        sys.exit(1)

    max_index = np.array([nlon, nlat])
    grid = esmpy.Grid(max_index, num_peri_dims=1, staggerloc=[esmpy.StaggerLoc.CENTER])

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


def build_grid_cons(bounds, nlat, nlon):
    """
    Build an esmpy (ESMF) grid object in memory using a modification 
    of Kieran's method at https://gist.github.com/kieranricardo/eb98f76235255efff800d28c2442e5c3
    """
    
    lon0, lon1 = bounds['lon_min'],bounds['lon_max']
    lat0, lat1 = bounds['lat_min'],bounds['lat_max']

    max_index = np.array([nlon, nlat])
    grid = esmpy.Grid(max_index, 
                    coord_sys=esmpy.CoordSys.SPH_DEG,
                    staggerloc=[esmpy.StaggerLoc.CORNER])


    # Get the coordinate arrays for the cell corners (vertices)
    grid_lon_c = grid.get_coords(0, staggerloc=esmpy.StaggerLoc.CORNER)
    grid_lat_c = grid.get_coords(1, staggerloc=esmpy.StaggerLoc.CORNER)

    # Define the coordinates for the corners (vertices)
    # The coordinate array size is now (nlat+1) x (nlog+1)
    lat_corners = np.linspace(lat0, lat1, max_index[0] + 1)
    lon_corners = np.linspace(lon0, lon1, max_index[1] + 1)

    # Assign coordinates to the grid object
    for j in range(max_index[0] + 1):
        for i in range(max_index[1] + 1):
            grid_lat_c[j, i] = lat_corners[j]
            grid_lon_c[j, i] = lon_corners[i]

    return grid, lat_corners, lon_corners


def regrid(atm_grid, ocn_mesh, ocn_mask, nlat, nlon, lats, lons):
    """
    Use the ESMF regridder to grid the ocn_mask onto the 
    atm_grid
    """
    ocn_mesh, ocn_mask, bounds = load_ocn_data(ocean_file)

    atm_field = esmpy.Field(atm_grid, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    atm_field.data[:] = 0.0    

    ocn_field = esmpy.Field(ocn_mesh, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    ocn_field.data[:] = ocn_mask

    ocn_to_atm_cons = esmpy.api.regrid.Regrid(
      ocn_field, atm_field, 
      unmapped_action=esmpy.api.constants.UnmappedAction.IGNORE,
      regrid_method=esmpy.api.constants.RegridMethod.CONSERVE,
      norm_type=esmpy.api.constants.NormType.DSTAREA, 
      factors=True
    )

    new_ocn_frac = atm_field.data.reshape((nlat, nlon))

    # Output to array
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
    data=new_ocn_frac,
    dims=["lat", "lon"],
    coords=dict(
        lat=lats,
        lon=lons,
    ),
    attrs=dict(um_stash_source='m01s00i505')
    )

    ds = xr.Dataset(dict(landfrac=da))

    return ds


def regrid_cons(grid,ocn_mesh, nlat, nlon, lat_corners,lon_corners):
    """
    Regrid my conservative grid object
    """
    src_field_conserve = esmpy.Field(grid, name='source_data_conserve')

    # Define the data on the centers of the cells (which is N_y x N_x)
    src_lat_centers = (lat_corners[:-1] + lat_corners[1:]) / 2.0
    src_lon_centers = (lon_corners[:-1] + lon_corners[1:]) / 2.0
    src_field_conserve.data[:] = 100.0 + np.outer(src_lat_centers, np.cos(np.deg2rad(src_lon_centers)))

    ocn_mesh, ocn_mask, bounds = load_ocn_data(ocean_file)

    ocn_field = esmpy.Field(ocn_mesh, meshloc=esmpy.api.constants.MeshLoc.ELEMENT)
    ocn_field.data[:] = ocn_mask

    ocn_to_atm_cons = esmpy.api.regrid.Regrid(
      ocn_field, src_field_conserve, 
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

    out_fp = 'dummy_kieran.nc'
    nlat=24
    nlon=24

    ocn_mesh, ocn_mask, bounds = load_ocn_data(ocean_file)
    atm_grid, lats, lons = build_grid(bounds, nlat)
    ds = regrid(atm_grid, ocn_mesh, ocn_mask, nlat, nlon, lats, lons)
    ds.to_netcdf(out_fp)

    my_atm_grid, lats, lons = build_my_grid(bounds, atm_res)
    nlon = 36
    my_ds = regrid(my_atm_grid, ocn_mesh, ocn_mask, nlat, nlon, lats,
    lons)

    out_fp = 'dummy_paul.nc'
    my_ds.to_netcdf(out_fp)

    atm_grid_cons, lats, lons = build_grid_cons(bounds, nlat, nlon)
    con_ds = regrid_cons(atm_grid_cons, ocn_mesh, nlat, nlon, lats, lons)

    out_fp = 'dummy_cons.nc'
    con_ds.to_netcdf(out_fp)