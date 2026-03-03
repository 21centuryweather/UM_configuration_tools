#!/usr/bin/env python
# Loads vegetation fraction ancillary and checks factions sm to one
# per https://code.metoffice.gov.uk/doc/ancil/ants/2.1/lib/ants.analysis.html#ants.analysis.cover_mapping.normalise_fractions
# works for ants 2.1.0 or on current xp65 conda/analysis3-25.08

import sys
import ants
import xarray as xr
import argparse

stashid = {'vegfrac'  : 'm01s00i216'}

parser = argparse.ArgumentParser(description='Fix vegetation fraction ancillary')
parser.add_argument('source_fpath', help='Path to source vegetation fraction file')
parser.add_argument('target_fpath', help='Path to target output file')
parser.add_argument('lnd_mask_fpath', help='Path to land sea mask file')
args = parser.parse_args()

source_fpath = args.source_fpath
target_fpath = args.target_fpath
lnd_mask = ants.load_cube(args.lnd_mask_fpath)

cube = ants.load_cube(source_fpath, constraint=stashid['vegfrac'])
cube_updated = cube.copy()

ants.analysis.cover_mapping.normalise_fractions(cube_updated)

# cube_updated = xr.DataArray.from_iris(cube_updated)
nans = cube_updated.data.mask & lnd_mask.data
cube_updated.data[nans.astype(bool)] = 0



# deal with different i/o in ANTS v0/v1 and v2
print('ANTS version:', ants.__version__)
if ants.__version__[0] != '2':
    print('saving with ants < 2 API')
    ants.save(cube_updated, target_fpath, zlib=True)
else:
    print('saving with ants v2 API')
    ants.io.save.ancil(cube_updated, target_fpath)
    ants.io.save.netcdf(cube_updated, f'{target_fpath}.nc', zlib=True)