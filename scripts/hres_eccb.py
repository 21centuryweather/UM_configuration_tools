# Copyright 2024 ACCESS-NRI (https://www.access-nri.org.au/)
# See the top-level COPYRIGHT.txt file for details.
#
# SPDX-License-Identifier: Apache-2.0
#
# Created by: Chermelle Engel <Chermelle.Engel@anu.edu.au>

"""
Replace the land/surface fields in the ec_cb000 file with higher-resolution
era5-land or BARRA2-R data (if requested).
"""

import argparse
import shutil
from pathlib import Path

import pandas

from replace_landsurface import replace_landsurface_with_BARRA2R_IC, replace_landsurface_with_ERA5land_IC

def main():
    """
    The main function that creates a worker pool and generates single GRIB files 
    for requested date/times in parallel.

    Parameters
    ----------
    None.  The arguments are given via the command-line

    Returns
    -------
    None.  The ec_cb000 file is updated and overwritten
    """

    # Parse the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--mask', required=True, type=Path)
    parser.add_argument('--file', required=True, type=Path)
    parser.add_argument('--start', required=True, type=pandas.to_datetime)
    parser.add_argument('--type', default="era5land")
    parser.add_argument('--hres_ic', type=Path)
    args = parser.parse_args()
    print(args)

    # Convert the date/time to a formatted string
    t = args.start.strftime("%Y%m%dT%H%MZ")
    print(args.mask, args.file, t)

    # If necessary replace ERA5 land/surface fields with higher-resolution options
    if "era5land" in args.type:
        replace_landsurface_with_ERA5land_IC.swap_land_era5land(args.mask, args.file, t)
        shutil.move(args.file.as_posix(), args.file.as_posix().replace('.tmp', ''))
    elif "barra" in args.type:
        replace_landsurface_with_BARRA2R_IC.swap_land_barra(args.mask, args.file, t)
        shutil.move(args.file.as_posix(), args.file.as_posix().replace('.tmp', ''))
    elif "astart" in args.type:
        print("Fields not swapped out for ECCB files when using start dump as replacement option.")
    else:
        print("No need to swap out IC")

if __name__ == '__main__':
    main()

