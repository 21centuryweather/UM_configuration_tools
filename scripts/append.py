#!/usr/bin/env python
"""
The CAP reads specifically 10 fields from the input file when generating the
snow ancillaries (ancilSmcsnow).  There is more than 10 fields present in the
soils ancillary which means that ancilSmcsnow is dependent on field order...
We move those fields it doesn't need to the end of the file to avoid a
segmentation fault.

"""
import warnings
import sys
sys.path.insert(0,'/home/548/pag548/cylc-run/rCM3-ancil-suite/share/fcm_make_ants/build/lib/')

import ants
import ants.io.save as save


def main(filenames, output):
    cubes = ants.load(filenames)

    # Sort by the following stash items (lbuser4).
    items = [40, 41, 43, 207, 47, 44, 46, 48, 220, 223, 8, 418, 419, 420]
    end = len(items) + 1
    
    def sortme(cube):
        stash = cube.attributes['STASH']
        try:
            order = items.index(stash.item)
        except ValueError:
            warnings.warn('Was this stash expected?? {}'.format(stash))
        return order
    cubes = sorted(cubes, key=sortme)

    # Check staggering
    grids = [
        cube.attributes["grid_staggering"]
        for cube in cubes
        if "grid_staggering" in cube.attributes
        ]
    unique_grids = set(grids)
    if len(unique_grids)> 1:
        # We have different staggeging
         most_common = max(set(grids), key=grids.count)

         # Alter cubes to enforce this stagger value
         for cube in cubes:
            if cube.attributes["grid_staggering"] != most_common:
                cube.attributes["grid_staggering"] = most_common

    save.ancil(cubes, output)


if __name__ == '__main__':
    parser = ants.AntsArgParser()
    args = parser.parse_args()
    main(args.sources, args.output)
