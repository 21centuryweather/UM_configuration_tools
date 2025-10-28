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

    # Create connectivity list
    # The mesh has (N_LON_NODES - 1) * (N_LAT_NODES - 1) elements
    N_ELEMENTS_LON = N_LON_NODES - 1
    N_ELEMENTS_LAT = N_LAT_NODES - 1
    N_ELEMENTS = N_ELEMENTS_LON * N_ELEMENTS_LAT

    connectivity = []
    
    # Loop over all elements in the mesh (rows 'r' and columns 'c')
    for r in range(N_ELEMENTS_LAT): # Element row index (0 to N_LAT_NODES - 2)
        for c in range(N_ELEMENTS_LON): # Element column index (0 to N_LON_NODES - 2)
            
            # The node index of the bottom-left corner (BL) of the current element (c, r)
            # Node indices are calculated in Row-Major order (row * width + col)
            BL = r * N_LON_NODES + c
            
            # Define the four nodes of the element in Counter-Clockwise (CCW) order:
            # BL: Bottom-Left (r, c)
            # BR: Bottom-Right (r, c+1)
            # TR: Top-Right (r+1, c+1)
            # TL: Top-Left (r+1, c)
            
            BR = BL + 1
            TR = BL + N_LON_NODES + 1
            TL = BL + N_LON_NODES
            
            # Append the indices for the current element
            connectivity.extend([BL, BR, TR, TL])

    element_to_node_conn = np.array(connectivity, dtype=np.int32)


    # Temp array - will regrid the ocean land/sea mask in due course
    data = np.ones([len(lat_coords),len(lon_coords)])

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


def build_mesh(lat_coords, lon_coords):
    """
    Build an esmpy (ESMF) mesh object in memory
    """
    N_LAT_NODES = len(lat_coords)
    N_LON_NODES = len(lon_coords)

    nx = N_LON_NODES-1 #No of cells
    ny = N_LAT_NODES-1

    x = lon_coords
    y = lat_coords

    # Total number of nodes
    num_nodes = (nx + 1) * (ny + 1)
    # Total number of quads
    num_elems = nx * ny

    # Create 2D coordinate arrays
    xv, yv = np.meshgrid(x, y, indexing='xy')

    # Flatten to 1D for ESMF
    nodeCoord = np.vstack([xv.flatten(), yv.flatten()]).T

    # Node IDs (1-based indexing)
    nodeId = np.arange(1, num_nodes + 1, dtype=np.int32)
    # Single process ownership (for parallel runs you'd set per-rank)
    nodeOwner = np.zeros(num_nodes, dtype=np.int32)

    elemId = np.arange(1, num_elems + 1, dtype=np.int32)
    elemType = np.full(num_elems, esmpy.MeshElemType.QUAD, dtype=np.int32)

    elemConn = []
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i
            n1 = n0 + 1
            n2 = n0 + nx + 2
            n3 = n0 + nx + 1
            elemConn.extend([n0, n1, n2, n3])
    elemConn = np.array(elemConn, dtype=np.int32)


    # Initialize ESMF
    mgr = esmpy.Manager(debug=True)

    # Create a Mesh with 2D parametric and spatial dimensions
    mesh = esmpy.Mesh(parametric_dim=2, spatial_dim=2)

    # Add nodes
    mesh.add_nodes(
        node_count=num_nodes,
        node_ids=nodeId,
        node_coords=nodeCoord,
        node_owners=nodeOwner
    )

    # Add elements (cells)
    mesh.add_elements(
        num_elem=num_elems,
        elemId=elemId,
        elemType=elemType,
        elemConn=elemConn
    )

    print(f"Mesh created: {mesh.size[esmpy.MeshLoc.NODE]} nodes, {mesh.size[esmmpy.MeshLoc.ELEM]} elements.")


def regrid_land_sea_mask(ocn_mesh, ocn_mask, atm_da):
    """
    Read in the ocean mesh, ocean land/sea mask and the atmospheric mesh defined as a dataarray
    Create an atmospheric mesh object and then regrid the ocean land/sea mask onto the atmospheric mesh.
    """


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
    atm_da = create_atmos_mask(bounds, atm_res)
