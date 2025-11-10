#!/usr/bin/env python
import esmpy
import xarray as xr
import numpy as np
from argparse import ArgumentParser
import os
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

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


def create_atmos_mask(bounds, resolution, mygrid=True):
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
    if mygrid:
        atm_grid = build_my_grid(lat_coords, lon_coords)
    else:
        atm_grid = build_grid(bounds, atm_res)
   
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


def build_grid(lat_coords, lon_coords):
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


def regrid(atm_grid, ocn_mesh, ocn_mask, nlat, nlon, lats, lons):
    """
    Use the ESMF regridder to grid the ocn_mask onto the 
    atm_grid
    """
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
    nlat=200
    nlon=nlat

    ocn_mesh, ocn_mask, bounds = load_ocn_data(ocean_file)
    atm_grid, lats, lons = build_grid(bounds, nlat)
    ds = regrid(atm_grid, ocn_mesh, ocn_mask, nlat, nlon, lats, lons)
    ds.to_netcdf(out_fp)

    my_atm_grid, lats, lons = build_my_grid(bounds, nlat)
    ds = regrid(my_atm_grid, ocn_mesh, ocn_mask, nlat, nlon, lats,
    lons)

    out_fp = 'dummy_paul.nc'
    ds.to_netcdf(out_fp)

    # Load in for visualisation
    ESMF_mesh = xr.load_dataset(ocean_file)
    ocean_mask = xr.open_dataset('/g/data/gb02/pag548/ocean-config-for-rcm3-dev/input/tiny-ocean-rcm3-dev/ocean_mask.nc')
   
    nx = len(ocean_mask.nx)
    ny = len(ocean_mask.ny)

    ESMF_mask = xr.DataArray(ESMF_mesh.elementMask.values.reshape(ny,nx), 
                dims=['ny','nx'],
                coords={'ny':ocean_mask.ny, 
                        'nx':ocean_mask.nx})
                        
    bathym = xr.open_dataset('/g/data/gb02/pag548/ocean-config-for-rcm3-dev/input/tiny-ocean-rcm3-dev/bathymetry.nc')

    hgrid = xr.open_dataset('/g/data/gb02/pag548/ocean-config-for-rcm3-dev/input/tiny-ocean-rcm3-dev/hgrid.nc')

    # To extract lat/lons from the hgrid (Madi's method)
    # --- MOM6 native mask 
    # Mask is defined at T-cell centers. Lat/lon at these points can be extracted from hgrid:
    geo_lon_t = hgrid.x[1::2,1::2]
    geo_lat_t = hgrid.y[1::2,1::2]

    ESMF_lat = xr.DataArray(ESMF_mesh.centerCoords[:, 1].values.reshape(ny,nx), 
                            dims=['ny','nx'],
                            coords={'ny':ocean_mask.ny, 
                                    'nx':ocean_mask.nx})

    ESMF_lon = xr.DataArray(ESMF_mesh.centerCoords[:, 0].values.reshape(ny,nx), 
                            dims=['ny','nx'],
                            coords={'ny':ocean_mask.ny, 
                                    'nx':ocean_mask.nx})

    ESMF_ds = xr.Dataset({
        'mask': ESMF_mask,
        'geo_lat_t':  ESMF_lat,
        'geo_lon_t':  ESMF_lon
    })

    # Plot these T cell centres using pcolormesh
    fig,ax=plt.subplots(1,1)
    p1 = ax.pcolormesh(ESMF_ds.geo_lon_t, ESMF_ds.geo_lat_t, ESMF_ds.mask, shading='auto', vmin = 0, vmax = 1, edgecolors='black')
    fig.colorbar(p1, ax=ax, orientation='vertical', label='Ocean mask ESMF')

    # Now for gemini solution
    # https://gemini.google.com/share/aa910f049570

    nodeCoords = ESMF_mesh.nodeCoords.data

    # Find the min and max 
    lon = nodeCoords[:,0]
    lat = nodeCoords[:,1]

    elementConn = ESMF_mesh.elementConn.values
    elementMask = ESMF_mesh.elementMask.values

    elementConn_0based = elementConn - 1

    tri1 = elementConn_0based[:, [0, 1, 2]]
    tri2 = elementConn_0based[:, [0, 2, 3]]

    # Stack the new triangles
    triangles = np.vstack([tri1, tri2])

    # Duplicate the mask data for the new set of triangles
    facecolors = np.hstack([elementMask, elementMask])

    # 2. Set up the Cartopy GeoAxes
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree()) # Or another projection

    plot = ax.tripcolor(
    lon,                       # X-coordinates (Longitude) of the nodes
    lat,                       # Y-coordinates (Latitude) of the nodes
    triangles=triangles,       # Connectivity matrix for the triangles
    facecolors=facecolors,     # Data defined at the element center (mask values)
    cmap='viridis',            # Custom land/sea colormap
    vmin=0, vmax=1,            # Set min/max for the colormap to ensure 0 and 1 are covered
    edgecolors='lightgray',         # 'k' or 'lightgray' to show grid lines
    transform=ccrs.PlateCarree() # Specify the coordinate system of the input data
    )

    # Plot lat/lon axes
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
 
