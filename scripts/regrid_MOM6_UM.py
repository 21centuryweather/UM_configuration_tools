#!/usr/bin/env python

import xarray as xr
import numpy as np
from argparse import ArgumentParser
import os
import sys
import esmpy
esmpy.Manager(debug=True)

GRID_STAGGER_DICT = { 0 : 'CENTER',
                      1 : 'EDGE1',
                      2 : 'EDGE2',
                      3 : 'CORNER'}

def add_options(parser):
    """
    Process the command-line arguments
    """
    parser.description = 'Regrid a MOM6 ocean mesh onto a regular grid with the same bounds.\n The user specifies the output resolution and output file'

    parser.add_argument('-i','--input-ocean-mesh-file', 
        dest='mesh_file', 
        type=str, 
        help="The rMOM6 esmpy mesh file, e.g. access-rom3-esmpymesh.nc")
    parser.add_argument('-o','--output-land-sea-mask-file', 
        dest='outfile', 
        type=str, 
        help="The output land-sea mask at the desired resolutio, e.g. target_mask.nc")
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


def build_grid_cons(bounds, 
                    res):
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

    grid_lon_c = grid.get_coords(coord_dim=0, staggerloc=esmpy.StaggerLoc.CENTER)
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


def create_dataarray(field,
                     index):
    """
    Create a datarray for a specified field
    """

    coords = field.grid.coords
    lons = coords[index][0]
    lats = coords[index][1]

    data = field.data

    # Output to DataArray
    lat_coord = xr.DataArray(
      dims=['latitude'],
      coords=dict(latitude=lats[:,0]),
      data=lats[:,0],
      attrs=dict(standard_name='latitude', units='degrees_north')
    )

    lon_coord = xr.DataArray(
      dims=['longitude'],
      coords=dict(longitude=lons[0]),
      data=lons[0],
      attrs=dict(standard_name='longitude', units='degrees_east')
     )

    da = xr.DataArray(
    data=data,
    dims=["latitude", "longitude"],
    coords=dict(latitude=lat_coord,
                longitude=lon_coord),
    attrs=dict(um_stash_source='m01s00i505',
               grid_mapping='geog_cs',
               earth_radius=6371229.0)
    )

    return da


def create_grid_mapping_variable():
    """
    Defines the co-ordinate system required for UM ancillary
    generation
    """
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

    return grid_mapping_var


def convert_gridded_field(field,
                         outfilename):
    """
    Outputs a gridded ESMF 2-Dfield object to xarray dataset
    """
    # Get data shape
    shape = field.data.shape

    # Get the available grid staggers. This returns a boolean list of length 4
    # The entries correspond to 
    # [ CENTER, EDGE1, EDGE2, CORNER ]
    # See https://earthsystemmodeling.org/esmpy_doc/release/ESMF_8_0_1/html/StaggerLoc.html#ESMF.api.constants.StaggerLoc

    stagger = field.grid.staggerloc
    coords = field.grid.coords
    for i,stag in enumerate(stagger):

        if stag:
            lons = coords[i][0]
            lats = coords[i][1]

            if (lons.shape == shape) and (lats.shape == shape):
                print (f"INFO : Field has data corresponding to {GRID_STAGGER_DICT[i]}")

                field_da = create_dataarray(field,
                                            i)
                grid_var = create_grid_mapping_variable()

                field_ds = xr.Dataset(
                    data_vars=dict(
                    land_binary_mask=field_da,
                    geog_cs=grid_var
                    ),
                    coords=dict(
                    latitude=field_da.latitude,
                    longitude=field_da.longitude,
                    ),
                    attrs=dict(Conventions='CF-1.7')
                    )

                field_ds.to_netcdf(outfilename)
                print(f"INFO : Writing {outfilename}")
    
    return field_ds


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
    target_mask_file = args.outfile
    atm_res = float(args.atm_res)

    #Load the Mesh file
    print(f"INFO : Opening {MOM6_mesh_file}")
    MOM6_mesh, MOM6_mask, bounds = load_ocn_data(MOM6_mesh_file)

    # Create the MOM6 mask field
    MOM6_field = esmpy.Field(MOM6_mesh,meshloc=esmpy.MeshLoc.ELEMENT)

    # Transfer the mask data into the field
    MOM6_field.data[:] = MOM6_mask

    # Create the target UM Grid
    UM_grid, lats, lons = build_grid_cons(bounds, 
                                          atm_res)

    # Create a UM dummy field based on the UM corners, defined at 
    # the grid centres
    UM_field = esmpy.Field(UM_grid,staggerloc=esmpy.StaggerLoc.CENTER)

    # To use the conservative regridder, we are going to define a dual grid composed of the UM grid's cell centres

    MOM6_to_UM_cons = esmpy.api.regrid.Regrid(
      MOM6_field, 
      UM_field, 
      unmapped_action=esmpy.api.constants.UnmappedAction.IGNORE,
      regrid_method=esmpy.api.constants.RegridMethod.CONSERVE,
      norm_type=esmpy.api.constants.NormType.DSTAREA, 
      factors=True
    )

    # Now we interpolate the UM field (defined at the grid centre) back to the vertices
    UM_field_v = esmpy.Field(UM_grid,staggerloc=esmpy.StaggerLoc.CORNER)
    
    UM_c_to_v = esmpy.api.regrid.Regrid(
      UM_field, 
      UM_field_v, 
      unmapped_action=esmpy.api.constants.UnmappedAction.IGNORE,
      regrid_method=esmpy.api.constants.RegridMethod.BILINEAR,
      norm_type=esmpy.api.constants.NormType.DSTAREA, 
      factors=True
    )

    # To do -output both for testing sake
    ds_c = convert_gridded_field(UM_field,
                        'target_centres.nc')

    ds_v = convert_gridded_field(UM_field_v,
                        'target_corners.nc')