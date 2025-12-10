import numpy as np
import xarray as xr
import os
import sys
import matplotlib.pyplot as plt

if 'ESMFMKFILE' not in os.environ:
    os.environ['ESMFMKFILE'] = '/opt/conda/analysis3-25.09/lib/esmf.mk'

import esmpy
# Initialize ESMF Manager (assuming this is done)
esmpy.Manager(debug=True)

# --- Source Grid Definition for Conservative Regridding ---
src_size = np.array([10, 20])  # Number of cell centers [Lat, Lon]
# Define grid asking for CORNER stagger location
src_grid_conserve = esmpy.Grid(src_size,
                               coord_sys=esmpy.CoordSys.SPH_DEG,
                               staggerloc=esmpy.StaggerLoc.CORNER)

# Get the coordinate arrays for the cell corners
src_lat_c = src_grid_conserve.get_coords(1, esmpy.StaggerLoc.CORNER) # Index 1 for Lat
src_lon_c = src_grid_conserve.get_coords(0, esmpy.StaggerLoc.CORNER) # Index 0 for Lon

# Define the coordinates for the CORNERS (from 20N-40N, 10W-30E)
# The coordinate array size is now (10+1) x (20+1)
lat_corners = np.linspace(20, 40, src_size[0] + 1)
lon_corners = np.linspace(-10, 30, src_size[1] + 1)

# Assign coordinates to the grid object
for j in range(src_size[0] + 1):
    for i in range(src_size[1] + 1):
        src_lat_c[j, i] = lat_corners[j]
        src_lon_c[j, i] = lon_corners[i]

# Create the Source Field
# The data is always placed at the CENTER stagger location by default,
# and ESMPy calculates the cell areas based on the corner coordinates.
src_field_conserve = esmpy.Field(src_grid_conserve, name='source_data_conserve')

# Define the data on the centers of the cells (which is N_y x N_x)
src_lat_centers = (lat_corners[:-1] + lat_corners[1:]) / 2.0
src_lon_centers = (lon_corners[:-1] + lon_corners[1:]) / 2.0
src_field_conserve.data[:] = 100.0 + np.outer(src_lat_centers, np.cos(np.deg2rad(src_lon_centers)))

# --- Destination Grid Definition for Conservative Regridding ---
dst_size = np.array([30, 60])
dst_grid_conserve = esmpy.Grid(dst_size,
                               coord_sys=esmpy.CoordSys.SPH_DEG,
                               staggerloc=esmpy.StaggerLoc.CORNER)

# Define coordinates for the CORNERS
dst_lat_c = dst_grid_conserve.get_coords(1, esmpy.StaggerLoc.CORNER)
dst_lon_c = dst_grid_conserve.get_coords(0, esmpy.StaggerLoc.CORNER)

dst_lat_corners = np.linspace(20, 40, dst_size[0] + 1)
dst_lon_corners = np.linspace(-10, 30, dst_size[1] + 1)

for j in range(dst_size[0] + 1):
    for i in range(dst_size[1] + 1):
        dst_lat_c[j, i] = dst_lat_corners[j]
        dst_lon_c[j, i] = dst_lon_corners[i]

# Define the data on the centers of the cells (which is N_y x N_x)
dst_lat_centers = (dst_lat_corners[:-1] + dst_lat_corners[1:]) / 2.0
dst_lon_centers = (dst_lon_corners[:-1] + dst_lon_corners[1:]) / 2.0

# Create the Destination Field
dst_field_conserve = esmpy.Field(dst_grid_conserve, name='destination_data_conserve')
dst_field_conserve.data[:] = 1e20

# --- Create the Regridder ---
regridder_conserve = esmpy.Regrid(src_field_conserve,
                                   dst_field_conserve,
                                   # === KEY CHANGE: Use CONSERVE method ===
                                   regrid_method=esmpy.RegridMethod.CONSERVE,
                                   unmapped_action=esmpy.UnmappedAction.IGNORE)

# --- Apply the Regridding ---
# The result is placed directly into dst_field_conserve.data
dst_field_conserve = regridder_conserve(src_field_conserve, dst_field_conserve)

# Get the regridded data
regridded_data_conserve = dst_field_conserve.data

# --- 5. Create Xarray DataArrays (src_da and dst_da) ---
print("Creating Xarray DataArrays...")

# Extract 1D coordinate arrays (assuming uniform grids)

# Create Source DataArray
src_da = xr.DataArray(
    src_field_conserve.data,
    coords={'lat': src_lat_centers, 'lon': src_lon_centers},
    dims=['lat', 'lon'],
    name='Source_Data'
)

# Create Destination/Regridded DataArray
dst_da = xr.DataArray(
    dst_field_conserve.data, 
    coords={'lat': dst_lat_centers, 'lon': dst_lon_centers},
    dims=['lat', 'lon'],
    name='Regridded_Data'
)
print("DataArrays created successfully.")


# --- 6. Plotting the Results ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('ESMPy Conservative Regridding: Source vs. Destination', fontsize=16)

# Subplot 1: Source Grid
src_da.plot.pcolormesh(
    ax=axes[0],
    x='lon',
    y='lat',
    vmin=110,
    cmap='viridis',
    add_colorbar=True
)
axes[0].set_title(f'Source Grid Data ({src_da.shape[0]}x{src_da.shape[1]})')
axes[0].set_xlabel('Longitude (°)')
axes[0].set_ylabel('Latitude (°)')

# Subplot 2: Destination/Regridded Grid
dst_da.plot.pcolormesh(
    ax=axes[1],
    x='lon',
    y='lat',
    cmap='viridis',
    vmin=110,
    add_colorbar=True
)
axes[1].set_title(f'Regridded Data ({dst_da.shape[0]}x{dst_da.shape[1]})')
axes[1].set_xlabel('Longitude (°)')
axes[1].set_ylabel('Latitude (°)')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# --- Verification ---
print("\nSource Data Min/Max:", src_field_conserve.data.min(), src_field_conserve.data.max())
print("Regridded Data Min/Max:", dst_field_conserve.data.min(), dst_field_conserve.data.max())
print("Shape of Regridded Data:", dst_field_conserve.data.shape)


# --- Cleanup ---
src_grid_conserve.destroy()
dst_grid_conserve.destroy()
src_field_conserve.destroy()
dst_field_conserve.destroy()
regridder_conserve.destroy()