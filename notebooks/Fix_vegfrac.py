#!/usr/bin/env python
# Loads vegetation fraction ancillary and checks factions sm to one
# per https://code.metoffice.gov.uk/doc/ancil/ants/2.1/lib/ants.analysis.html#ants.analysis.cover_mapping.normalise_fractions
# works for ants 2.1.0 or on current xp65 conda/analysis3-25.08

import sys
import ants
import xarray as xr

stashid = {'vegfrac'  : 'm01s00i216'}

source_fpath='/scratch/gb02/sl5165/u-dr651/n320e/mom025_20250315/cable_veg/cable_vegfrac_1850.anc'
target_fpath = '/scratch/gb02/pc2687/cable_vegfrac_1850.anc'

cube = ants.load_cube(source_fpath, constraint=stashid['vegfrac'])
cube_updated = cube.copy()
lnd_mask = ants.load_cube('/g/data/gb02/public/AM3/ancils/u-dj813_om3-025deg/n320e/mom025_20250315/land_sea_mask/etop01/qrparm.mask')

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