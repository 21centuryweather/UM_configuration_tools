from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import PowerNorm
from PIL import Image

# Plot OLR over map using AM3 outputs
# Author: Shixue Li
# Updated by: Pao Corrales

year = "2007"
month = "feb"

file_dict = {
    "n96":  f"/g/data/gb02/public/AM3/output/control/n96/am3a.pa{year}{month}.nc",
    "n512": f"/g/data/av17/access-nri/AM3/n512/test2a.pa{year}{month}.nc"
}

target_names = ["n96", "n512"]

img        = Image.open("/g/data/gb02/sl5165/am_view/world.topo.bathy.200408.3x5400x2700.jpg")
img_array  = np.array(img.convert("RGB"))
img_extent = [-180, 180, -90, 90]

def read_data(file_dict, target_names):
    data = {}
    for name, path in file_dict.items():
        if name not in target_names:
            continue

        ds = xr.open_dataset(path)
        data[name] = {
            "lon": ds["lon"].values,
            "lat": ds["lat"].values,
            "time": ds["time_0"].values,
            "olr": ds["fld_s02i205"].values,
        }
        ds.close()
    return data

def save_one_time(t, data, img_array, img_extent, out_dir="/home/565/pc2687/plots"):
    names = list(data.keys())
    n_panels = len(names)

    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(14, 5),
        subplot_kw={"projection": ccrs.Robinson(central_longitude=180)}
    )
    if n_panels == 1:
        axes = [axes]

    for ax, name in zip(axes, names):
        olr = data[name]["olr"][t, :, :]
        olr_plot = np.ma.masked_where(olr > 240, olr)  # optional: hide high values

        lons = data[name]["lon"]
        lats = data[name]["lat"]

        ax.set_global()
        ax.imshow(
            img_array,
            origin="upper",
            extent=img_extent,
            transform=ccrs.PlateCarree(),
            zorder=0
        )

        mesh = ax.pcolormesh(
            lons, lats, olr_plot,
            transform=ccrs.PlateCarree(),
            cmap="Greys",
            norm=PowerNorm(gamma=2.0, vmin=0, vmax=240),
            shading="auto",
            zorder=1
        )
        ax.set_title(f"{name.upper()}")

    cbar_ax = fig.add_axes([0.25, 0.08, 0.5, 0.03])
    fig.colorbar(mesh, cax=cbar_ax, orientation="horizontal").set_label("OLR (W/m²)")

    ts = pd.to_datetime(data["n96"]["time"][t]) + pd.DateOffset(months=6, days=11)
    this_time = ts.strftime("%Y%m%d_%H%M%S")

    out_path = Path(out_dir)
    out_path.mkdir(exist_ok=True)

    plt.tight_layout(rect=[0, 0.12, 1, 1])
    plt.savefig(out_path / f"olr_n96_n512_-{this_time}.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    data = read_data(file_dict, target_names)

    n_time = min(data[name]["olr"].shape[0] for name in target_names)

    with ProcessPoolExecutor(max_workers=8) as pool:
        futures = [
            pool.submit(
                save_one_time,
                t,
                data,
                img_array,
                img_extent,
                "/scratch/gb02/pc2687/am3-plots/"
            )
            for t in range(n_time)
        ]

        for f in futures:
            f.result()