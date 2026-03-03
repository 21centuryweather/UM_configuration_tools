# Create a target land-sea mask using subset of the standard
# ERA5 land-sea mask from the default ACCESS-rAM3 Lismore domain

import xarray as xr
from pathlib import Path

SOURCE_DIR = Path('/home/548/pag548/cylc-run/u-dg767/share/data/ancils/Lismore/era5')
TARGET_DIR = Path('/home/548/pag548/code/UM_config_tools/scripts')

source_mask_file = SOURCE_DIR / 'qrparm.mask_sea.nc'
target_mask_file = TARGET_DIR / 'mask_test.nc'

# Open the source mask

sm = xr.open_dataset(source_mask_file).land_binary_mask

# Take a subset
# Without land
subset = sm[:50,:50]
# With land
#subset = sm[50:100,50:100]

da = xr.DataArray(
    data=subset.data,
    dims=["latitude", "longitude"],
    attrs=dict(um_stash_source='m01s00i505',
               grid_mapping='geog_cs',
               earth_radius=6371229.0)
    )
    

# Createa GeogCS Grid Mapping Variable
grid_mapping_var = xr.Variable(
        dims=(),
        data=0,
        attrs={
            'grid_mapping_name': 'latitude_longitude',
            'earth_radius': 6371229.0,
        }
    )

# Write the dataset
ds = xr.Dataset(
        data_vars=dict(
            land_binary_mask=da,
            geog_cs=grid_mapping_var
            ),
        coords=dict(
            latitude=subset.latitude,
            longitude=subset.longitude,
        ),
        attrs=dict(Conventions='CF-1.7')
        )

ds.to_netcdf('mask_test.nc')