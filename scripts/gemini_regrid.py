import numpy as np
import esmpy
import xarray as xr
import matplotlib.pyplot as plt

# 1. Initialize ESMF Manager
print("Initializing ESMPy Manager...")
esmpy.Manager(debug=True)

# --- 2. Define the Source Grid and Field (Coarse) ---
src_size = np.array([10, 20])  # [Lat, Lon] dimensions (y, x)
src_grid = esmpy.Grid(src_size,
                      coord_sys=esmpy.CoordSys.SPH_DEG, 
                      staggerloc=esmpy.StaggerLoc.CENTER) 

# Define coordinates
src_lat = src_grid.get_coords(1)
src_lon = src_grid.get_coords(0)

lat_bounds = np.linspace(20, 40, src_size[0])
lon_bounds = np.linspace(-10, 30, src_size[1])

for j in range(src_size[0]):
    for i in range(src_size[1]):
        src_lat[j, i] = lat_bounds[j]
        src_lon[j, i] = lon_bounds[i]

# Create Source Field and data
src_field = esmpy.Field(src_grid, name='source_data')
# Simple data function
src_field.data[:] = 100.0 + np.outer(lat_bounds, np.cos(np.deg2rad(lon_bounds)))

# --- 3. Define the Destination Grid and Field (Fine) ---
dst_size = np.array([30, 60])
dst_grid = esmpy.Grid(dst_size,
                      coord_sys=esmpy.CoordSys.SPH_DEG,
                      staggerloc=esmpy.StaggerLoc.CENTER)

# Define coordinates
dst_lat = dst_grid.get_coords(1)
dst_lon = dst_grid.get_coords(0)

dst_lat_bounds = np.linspace(20, 40, dst_size[0])
dst_lon_bounds = np.linspace(-10, 30, dst_size[1])

for j in range(dst_size[0]):
    for i in range(dst_size[1]):
        dst_lat[j, i] = dst_lat_bounds[j]
        dst_lon[j, i] = dst_lon_bounds[i]

# Create Destination Field
dst_field = esmpy.Field(dst_grid, name='destination_data')
dst_field.data[:] = 1e20

# --- 4. Create and Apply Regridder ---
regridder = esmpy.Regrid(src_field,
                         dst_field,
                         regrid_method=esmpy.RegridMethod.BILINEAR, 
                         unmapped_action=esmpy.UnmappedAction.IGNORE,
                         extrap_method=esmpy.ExtrapMethod.NEAREST_IDAVG)

# Perform the regridding
print("Performing regridding...")
dst_field = regridder(src_field, dst_field)


# --- 5. Create Xarray DataArrays (src_da and dst_da) ---
print("Creating Xarray DataArrays...")

# Extract 1D coordinate arrays (assuming uniform grids)
src_lat_1D = src_lat[:, 0]
src_lon_1D = src_lon[0, :]
dst_lat_1D = dst_lat[:, 0]
dst_lon_1D = dst_lon[0, :]

# Create Source DataArray
src_da = xr.DataArray(
    src_field.data,
    coords={'lat': src_lat_1D, 'lon': src_lon_1D},
    dims=['lat', 'lon'],
    name='Source_Data'
)

# Create Destination/Regridded DataArray
dst_da = xr.DataArray(
    dst_field.data, 
    coords={'lat': dst_lat_1D, 'lon': dst_lon_1D},
    dims=['lat', 'lon'],
    name='Regridded_Data'
)
print("DataArrays created successfully.")


# --- 6. Plotting the Results ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle('ESMPy Bilinear Regridding: Source vs. Destination', fontsize=16)

# Subplot 1: Source Grid
src_da.plot.pcolormesh(
    ax=axes[0],
    x='lon',
    y='lat',
    vmin=100,
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
    vmin=100,
    add_colorbar=True
)
axes[1].set_title(f'Regridded Data ({dst_da.shape[0]}x{dst_da.shape[1]})')
axes[1].set_xlabel('Longitude (°)')
axes[1].set_ylabel('Latitude (°)')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# --- Verification ---
print("\nSource Data Min/Max:", src_field.data.min(), src_field.data.max())
print("Regridded Data Min/Max:", dst_field.data.min(), dst_field.data.max())
print("Shape of Regridded Data:", dst_field.data.shape)

# --- Cleanup ---
# Destroy the objects to release memory
src_grid.destroy()
dst_grid.destroy()
src_field.destroy()
dst_field.destroy()
regridder.destroy()

print("\nCleanup complete. ESMF objects destroyed.")