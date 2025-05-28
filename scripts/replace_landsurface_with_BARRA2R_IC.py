# Copyright 2024 ACCESS-NRI (https://www.access-nri.org.au/)
# See the top-level COPYRIGHT.txt file for details.
#
# SPDX-License-Identifier: Apache-2.0
#
# Created by: Chermelle Engel <Chermelle.Engel@anu.edu.au>

import os
import sys
from glob import glob
from pathlib import Path

import iris
import mule
import numpy as np
import xarray as xr

ROSE_DATA = os.environ.get('ROSE_DATA', "")
# Base directory of the ERA5-land archive on NCI
BARRA_DIR = os.path.join(ROSE_DATA, 'etc', 'barra_r2')


class ReplaceOperator(mule.DataOperator):
    """ Mule operator for replacing the data"""
    def __init__(self):
        pass
    def new_field(self, sources):
        #print('new_field')
        return sources[0]
    def transform(self, sources, result):
        #print('transform')
        return sources[1]


class bounding_box(): 
    """ Container class to hold spatial extent information."""
    def __init__(self, ncfname, maskfname, var):
        """
        Initialization function for bounding_box class

        Parameters
        ----------
        ncfname : POSIX string
            POSIX path to input data (from the NetCDF archive)
        maskfname : POSIX string
            POSIX path to mask information to define the data to be cut out
        var : string
            The name of the mask variable that defines the spatial extent

        Returns
        -------
        None.  The variables describing the spatial extent are definied within bounding box object.

        """
        # Read in the mask and get the minimum/maximum latitude and longitude information
        if Path(maskfname).exists():
            d = iris.load(maskfname, var)
            d = d[0]
            d = xr.DataArray.from_iris(d)

            lons = d['longitude'].data
            lonmin = np.min(lons)
            lonmax = np.max(lons)
            lats = d['latitude'].data
            latmin = np.min(lats)
            latmax = np.max(lats)
            d.close()
        else:
            print(f'ERROR: File {maskfname} not found', file=sys.stderr)
            raise

        # Read in the file from the high-res netcdf archive
        if Path(ncfname).exists():
            d = xr.open_dataset(ncfname)
        else:
            print(f'ERROR: File {ncfname} not found', file=sys.stderr)
            sys.exit(1)


        # Get the longitude information
        lons = d['lon'].data

        # Coping with numerical inaccuracy
        adj = (lons[1] - lons[0])/2.

        # Work out which longitudes define the minimum/maximum extents of the grid of interest
        lonmin_index = np.argwhere( (lons > lonmin - adj) & (lons < lonmin + adj))[0][0]
        lonmax_index = np.argwhere( (lons > lonmax - adj) & (lons < lonmax + adj))[0][0]

        # Get the latitude information
        lats = d['lat'].data

        # Work out which longitudes define the minimum/maximum extents of the grid of interest
        # use the same adjustment as for longitude
        latmin_index = np.argwhere( (lats > latmin - adj) & (lats < latmin + adj))[0][0]
        latmax_index = np.argwhere( (lats > latmax - adj) & (lats < latmax + adj))[0][0]

        # Set the boundaries
        self.lonmin = lonmin_index
        self.lonmax = lonmax_index
        self.latmin = latmin_index
        self.latmax = latmax_index


def get_BARRA_nc_data(ncfname, FIELDN, wanted_dt, NLAYERS, bounds):
    """
    Function to get the BARA2-R data for a single land/surface variable.

    Parameters
    ----------
    ncfname : string
        The name of the file to read
    FIELDN : string
        The name of the variable in the file to read
    wanted_dt : string
        The date-time required in "%Y%m%d%H%M" format
    NLAYERS : int
        The number of layers in the multi-resolution grid (1 or more)
    bounds : bounding_box object
        A bounding box object defining the spatial extent to keep

    Returns
    -------
    2d numpy array 
        A 2-D numpy array containg the field data for the date/time and spatial extent
    """

    # Open the file containing the data
    if Path(ncfname).exists():
        d = xr.open_dataset(ncfname)
    else:
        print(f'ERROR: File {ncfname} not found', file=sys.stderr)
        sys.exit(1)
    
    # Find the array index for the date/time of interest
    times = d['time'].dt.strftime("%Y%m%d%H%M").data
    TM = times.tolist().index(wanted_dt)

    # retrieve the spatial extends of interest
    lonmin_index, lonmax_index = bounds.lonmin, bounds.lonmax
    latmin_index, latmax_index = bounds.latmin, bounds.latmax

    # Read the data
    try:
        if NLAYERS>1:
          data = d[FIELDN][TM, :, latmin_index:latmax_index+1, lonmin_index:lonmax_index+1]
        else:
          data = d[FIELDN][TM, latmin_index:latmax_index+1, lonmin_index:lonmax_index+1]
    except KeyError:
        print(fname)
        print(f'ERROR: Variable temp not found in file', file=sys.stderr)
        sys.exit(1)

    d.close()

    return data.data


def swap_land_barra(mask_fullpath, ec_cb_file_fullpath, ic_date):
    """
    Function to get the BARRA2-R data for all land/surface variables.

    Parameters
    ----------
    mask_fullpath : Path
        Path to the mask defining the spatial extent
    ic_file_fullpath : string
        Path to file with the coarser resolution data to be replaced with ".tmp" appended at end
    ic_date : string
        The date-time required in "%Y%m%d%H%M" format

    Returns
    -------
    None.
        The file is replaced with a version of itself holding the higher-resolution data.
    """


    # create name of file to be replaced
    ec_cb_file = ec_cb_file_fullpath.parts[-1].replace('.tmp', '')
   
    # create date/time useful information
    yyyy = ic_date[0:4]
    mm = ic_date[4:6]
    ic_z_date = ic_date.replace('T', '').replace('Z', '')
    
    # Path to input file
    ff_in = ec_cb_file_fullpath.as_posix().replace('.tmp', '')

    # Path to output file
    ff_out = ec_cb_file_fullpath.as_posix()
    #print(ff_in, ff_out)
    
    # Read input file
    mf_in = mule.load_umfile(ff_in)
    
    # Create Mule Replacement Operator
    replace = ReplaceOperator() 

    # Read in the surface temperature data from the archive
    BARRA_FIELDN = 'ts'
    indir = os.path.join(BARRA_DIR, '1hr',BARRA_FIELDN, 'latest')
    barra_files = glob(os.path.join(indir, BARRA_FIELDN + '*' + yyyy + mm + '*nc'))
    barra_fname = indir + '/' + barra_files[0].split('/')[-1]

    # Work out the grid bounds using the surface temperature file
    bounds = bounding_box(barra_fname, mask_fullpath.as_posix(), "land_binary_mask")

    # Read in the surface temperature data (and keep to use for replacement)
    data = get_BARRA_nc_data(barra_fname, BARRA_FIELDN, ic_z_date ,-1 ,bounds)
    surface_temp = data.copy()
    
    # Read in the soil moisture data (and keep to use for replacement)
    BARRA_FIELDN = 'mrsol'
    indir = os.path.join(BARRA_DIR, '3hr',BARRA_FIELDN, 'latest')
    barra_files = glob(os.path.join(indir, BARRA_FIELDN + '*' + yyyy + mm + '*nc'))
    barra_fname = indir + '/' + barra_files[0].split('/')[-1]
    data = get_BARRA_nc_data(barra_fname, BARRA_FIELDN, ic_date.replace('T', '').replace('Z', ''), 4, bounds)
    mrsol = data.copy()

    # Read in the soil temperature data (and keep to use for replacement)
    BARRA_FIELDN = 'tsl'
    indir = os.path.join(BARRA_DIR, '3hr',BARRA_FIELDN, 'latest')
    barra_files = glob(os.path.join(indir, BARRA_FIELDN + '*' + yyyy + mm + '*nc'))
    barra_fname = indir + '/' + barra_files[0].split('/')[-1]
    data = get_BARRA_nc_data(barra_fname, BARRA_FIELDN, ic_date.replace('T', '').replace('Z', ''), 4, bounds)
    tsl = data.copy()
    
    # Set up the output file
    mf_out = mf_in.copy()

    # For each field in the input write to the output file (but modify as required)
    for f in mf_in.fields:
    
      #print(f.lbuser4, f.lblev, f.lblrec, f.lbhr, f.lbcode)
      if f.lbuser4 == 9:
        # replace coarse soil moisture with high-res information
        current_data = f.get_data()
        data = mrsol[f.lblev-1, :, :]
        data = np.where(np.isnan(data), current_data, data)
        mf_out.fields.append(replace([f, data]))
      elif f.lbuser4 == 20:
        # replace coarse soil temperature with high-res information
        current_data = f.get_data()
        data = tsl[f.lblev-1, :, :]
        data = np.where(np.isnan(data), current_data, data)
        mf_out.fields.append(replace([f, data]))
      elif f.lbuser4 == 24:
        # replace surface temperature with high-res information
        current_data = f.get_data()
        data = surface_temp
        data = np.where(np.isnan(data), current_data, data)
        mf_out.fields.append(replace([f, data]))
      else:
        mf_out.fields.append(f)
    
    mf_out.validate = lambda *args, **kwargs: True
    mf_out.to_file(ff_out)
