# Build a script using very low level iris routines to
# correctly read a 0.1 degree land-sea mask at the correct
# precision

import iris
from iris.fileformats._ff import FF2PP
from pathlib import Path


ANCIL_DIR = Path('/scratch/gb02/pag548/cylc-run/rCM3-test-UM-ancil/share/data/ancils/Lismore/d1100')
MASK = ANCIL_DIR / 'qrparm.mask'
breakpoint()
mask = iris.load(MASK)

ff2pp = FF2PP(MASK, read_data=False)

fields = list(ff2pp)

bdx = fields[0].bdx
bzx = fields[0].bzx
lbnpt= fields[0].lbnpt

# Note these are read as 64 bit floats