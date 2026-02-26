# Build a script using very low level iris routines to
# correctly read a 0.1 degree land-sea mask at the correct
# precision

import iris
from iris.fileformats._ff import FF2PP
from pathlib import Path
from iris.fileformats.pp_load_rules import convert
import numpy as np

def regular_points(zeroth, step, count):
    """Make an array of regular points.

    Create an array of `count` points from `zeroth` + `step`, adding `step` each
    time. In float32 if this gives a sufficiently regular array (tested with
    points_step) and float64 if not.

    Parameters
    ----------
    zeroth : number
        The value *prior* to the first point value.
    step : number
        The numeric difference between successive point values.
    count : number
        The number of point values.

    Notes
    -----
    This function does maintain laziness when called; it doesn't realise data.
    See more at :doc:`/userguide/real_and_lazy_data`.
    """

    def make_steps(dtype: np.dtype):
        start = np.add(zeroth, step, dtype=dtype)
        steps = np.multiply(step, np.arange(count), dtype=dtype)
        return np.add(start, steps, dtype=dtype)

    points = make_steps(np.float64)
    
    return points

def replace_lat_lon(cube,filename:
    """
    Replace the 32-bit lat and lon points of input cube with 64 bit points
    contructed from the PP field of the input cube
    """

    field, = FF2PP(filename), read_data=False)

    bdx = field.bdx
    bzx = field.bzx
    lbnpt= field.lbnpt

    bdy = field.bdy
    bzy = field.bzy
    lbrow = field.lbrow

    lon_points = regular_points(bzx, bdx, lbnpt)
    lat_points = regular_points(bzy, bdy, lbrow)

    cube.coord('longitude').points = lon_points
    cube.coord('latitude').points = lat_points

    return cube


ANCIL_DIR = Path('/scratch/gb02/pag548/cylc-run/rCM3-test-UM-ancil/share/data/ancils/Lismore/d1100')
MASK = ANCIL_DIR / 'qrparm.mask'
#breakpoint()
mask, = iris.load(MASK)

ff2pp, = FF2PP(MASK, read_data=False)

bdx = ff2pp.bdx
bzx = ff2pp.bzx
lbnpt= ff2pp.lbnpt

bdy = ff2pp.bdy
bzy = ff2pp.bzy
lbrow = ff2pp.lbrow

# Note these are read as 64 bit floats

lon_points = regular_points(bzx, bdx, lbnpt)
lat_points = regular_points(bzy, bdy, lbrow)

mask = replace_lat_lon(mask,ff2pp)
