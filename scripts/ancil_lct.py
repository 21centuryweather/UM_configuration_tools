#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
"""
Land cover type fraction application
************************************

Derive the land cover type fraction, land cover fraction and optionally the
landsea mask ancillaries by utilising the provided landcover types source with
its corresponding transform (crosswalk table), mapping the source classes to
the output JULES target classes.
If providing a pre-existing landsea mask, the land cover type fraction will be
made consistent with it, otherwise, the land cover type fraction field itself
will be used to derive and output a landsea mask.

Here are the steps taken:

- Transform source land classifications to the target classifications
  (ordinarily this target classification will be JULES classes).
- Regrid these target classes to the provided target grid.
- Extract the land cover fraction field ('m01s00i505').
- Set ice to whole fractions by applying a >= 50% threshold of ice to 1,
  otherwise set to 0 and re-distributing differences amongst non-zero
  land fractions (the UM defines the requirement for whole fraction ice).
- Remove non-glacial ice (driven by UM behaviour with isolated ice points).
    - Identify ice seed points by those locations of ice which are surrounded
      by 5x5 cells of ice.
    - Find contiguous ice regions using these seed points by search
      neighbouring cells (including diagonals).
    - Replace ice with soils for those locations not amongst these contiguous
      regions of ice found.
- Ensure land fractions add up to 1 and ocean fraction add up to 1 by:
    - If an optionally supplied landsea mask (``--target-lsm``) is provided, ensure
      fields are consistent with the supplied mask.
    - If the mask is to be derived by the land cover type fraction field
      (``--landseamask-out``), a threshold of >= 50% ocean is to mean ocean and
      differences are redistributed amongst non-zero land type fraction fields.

Fields returned:
- JULES land cover type fraction ancillary fields ('m01s00i216').
- Land sea mask ('m01s00i030')
- Land area fractions ('m01s00i505')

"""
import json
import os
import warnings

import sys
sys.path.insert(0,'/home/548/pag548/cylc-run/u-dq487/share/fcm_make_ants/build/lib/')

import ants
import ants.decomposition as decomp
import ants.fileformats.cover_mapping as cover_mapping
import ants.io.save as save
import ants.utils
import iris
import numpy as np
from proc_ants import lct


def _load_ostia(ostia):
    # DEV: This function is here as part of development but isn't currently used.
    #      The OSTIA mask cannot be loaded without processing due to bad metadata.
    # Note that there may be at some point some cause to interpret the OSTIA
    # mask in some way and make some conditional on the lakes included in the
    # mask based on their presence in OSTIA.
    #
    # For example, we can fish out lat-lon points for each lake geometry in
    # Natural Earth as follows:
    # >>> records = proc_ants.lakes.get_lake_geoms()
    # >>> for record in records:
    #    geom = record.geometry.buffer(0)
    #    x, y = geom.representative_point().xy
    #    x = x[0]
    #    y = y[0]
    #
    # However, we should note that not all lakes have identified geometries in
    # Natural Earth or necessarily with one that is aligned to data derived
    # from the CCI.  This means that relating lakes between OSTIA and the CCI
    # might be rather difficult to achieve programatically.
    #
    fields = iris.fileformats.pp.load(ostia)

    def mutate(field):
        field.lbproc = 0
        return field

    tweaked_fields = (mutate(field) for field in fields)
    ostia_cubes = ants.fileformats.pp.pp2cubes(tweaked_fields)
    assert len(ostia_cubes) == 1, "More than one OSTIA cube?"
    ostia_cube = ostia_cubes[0]
    ostia_lsm = ostia_cube.copy(np.ma.getmaskarray(ostia_cube.data))
    return ostia_lsm


def _prepare_mask_cube(cube):
    """
    Sets the attributes needed for mask files, and converts datatype to integer.

    Operates in place.

    Attributes 'valid_min' and 'valid_max' are set to ensure netCDF result is
    correct for boolean data.

    STASH code is set to m01s00i030.

    Data type is set to integer.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Mask cube to prepare.

    Returns
    -------
    : None
    Operates in place.

    """
    cube.attributes["valid_min"] = 0
    cube.attributes["valid_max"] = 1
    cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i030")
    cube.data = cube.data.astype("int8")


def _validate_transform_against_source(cube, transform_path):
    """
    Checks whether a transform file is appropriate for a given cube and raises
    an error if not.

    Parameters
    ----------
    cube : :class:`iris.cube.Cube`
        Cube to check against
    transform_path : str
        Path to the transform json file to be used with the cube
    """

    with open(transform_path) as mappings_file:
        mappings = json.load(mappings_file)

    # The transformer is tolerant of different cases for using the markings
    # convert all to lower case here though for ease of comparison.
    labels = [flag.lower() for flag in mappings["source"]]
    types_in_cube = [
        flag.lower() for flag in cube.attributes["flag_meanings"].split(" ")
    ]

    # Iterate over all labels in the cube attributes and warn if they are
    # from the transforms file
    missing_details = []
    for source_type in types_in_cube:
        if source_type not in labels:
            missing_details.append(source_type)

    if missing_details:
        details = ("\n").join(missing_details)
        raise ValueError(
            "One or more source types not present in transform file:\n" + details
        )


def load_data(
    source_path,
    transform_path,
    target_grid=None,
    target_landseamask=None,
    land_fraction_threshold=None,
):
    source = ants.io.load.load_cube(source_path)
    if target_grid:
        target_cube = ants.io.load.load_grid(target_grid)
    else:
        target_cube = ants.io.load.load_landsea_mask(
            target_landseamask, land_fraction_threshold
        )

    source = source.extract(ants.ExtractConstraint(target_cube))
    target_cube = target_cube.extract(ants.ExtractConstraint(source, fix_period=True))

    _validate_transform_against_source(source, transform_path)

    trans = cover_mapping.load_cover_mapper(transform_path)
    return source, target_cube, trans


def gen_lct(src_cube, grid_cube, transform, min_frac=0.5):
    operation = ants.analysis.SCTTransformer(transform)
    lct_cube = decomp.decompose(operation, src_cube, grid_cube)
    lct_cube.units = "1"
    lct_cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi("m01s00i216")
    lct_cube.rename("vegetation_area_fraction")

    # Any masked data is nan in the land cover type fraction as the ocean level
    # is used to derive the landsea mask.
    if np.ma.is_masked(lct_cube.data):
        lct_cube.data[lct_cube.data.mask] = np.nan
        lct_cube.data = lct_cube.data.data

    # Derive the land area fraction before we remove the ocean.
    # This field is not currently used directly by the model.
    land_fraction_cube = lct_cube.extract(iris.Constraint(pseudo_level=0))
    land_fraction_cube.remove_coord("pseudo_level")
    land_fraction_cube.data = 1 - land_fraction_cube.data
    land_fraction_cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(
        "m01s00i505"
    )
    land_fraction_cube.rename("land_area_fraction")

    lct.set_whole_fraction_ice(lct_cube)
    lct.remove_non_glacial_ice(lct_cube)
    # Exclude the ocean level - All land fractions must add up to 1 which means
    # that points must be 100% land or 100% ocean.  The array is masked as a
    # result (representing this ocean).
    lct_cube, land_mask_cube = lct.remove_ocean_level(lct_cube, min_frac=min_frac)

    # Extract landsea mask.
    land_mask_cube = lct_cube.slices_over("pseudo_level").next()
    land_mask_cube.remove_coord("pseudo_level")
    # Ensure land_mask_cube is appropriately masked
    ants.utils.cube.fix_mask(land_mask_cube)
    land_mask_cube.data = np.logical_not(land_mask_cube.data.mask).astype(
        "int8", copy=False
    )
    nan_values = np.isnan(lct_cube.data[0])
    # Inherit unknown values from the lct in the form of a mask.
    if nan_values.any():
        land_mask_cube.data = np.ma.array(land_mask_cube.data, mask=nan_values)
        warnings.warn(
            "There are nan values in the vegetation fraction. These will translate "
            "to masked data in the land sea mask and may cause unexpected results "
            "if the mask is used with ants.analysis.make_consistent_with_lsm."
        )

    _prepare_mask_cube(land_mask_cube)
    land_mask_cube.rename("land_binary_mask")
    lsm_cubes = iris.cube.CubeList([land_mask_cube, land_fraction_cube])
    return lct_cube, lsm_cubes


def main(
    source_path,
    transform_path,
    out_filepath,
    target_grid,
    landseamask_in,
    landseamask_out_root,
    land_fraction_threshold,
    netcdf_only,
):
    (
        source,
        grid,
        src_trans,
    ) = load_data(
        source_path,
        transform_path,
        target_grid,
        landseamask_in,
        land_fraction_threshold,
    )

    # Select land fraction to mask out
    min_frac = 0.5
    if landseamask_in:
        min_frac = 0.0

    lct_cube, lsm_cubes = gen_lct(source, grid, src_trans, min_frac=min_frac)

    landseamask_out_root=True
    out_dir = '/home/548/pag548/cylc-run/u-dg767/share/data/ancils/Lismore/era5/'
    if landseamask_out_root:
        ants.config.dirpath_writeable(out_dir)
        _prepare_mask_cube(lsm_cubes[0])
        land_mask = lsm_cubes[0]
        land_fraction = lsm_cubes[1]
        ocean_data = np.ones_like(land_fraction.data)
        ocean_data[land_fraction.data < 1] = 0
        ocean_mask = land_fraction.copy(ocean_data)
        ocean_mask.standard_name = "land_binary_mask"
        _prepare_mask_cube(ocean_mask)

        cubes_to_save = {
            # "filename": (cube, fill_value),
            "qrparm.mask": (land_mask, -1),
            "qrparm.mask_sea": (ocean_mask, -1),
            "qrparm.landfrac": (land_fraction, None),
        }
        for filename, (cube, fill_value) in cubes_to_save.items():
            filepath = os.path.join(out_dir, filename)
            save_ants(cube, filepath, fill_value, netcdf_only)

    if landseamask_in:
        ants.analysis.make_consistent_with_lsm(lct_cube, grid, True)
    save_ants(lct_cube, out_filepath, None, netcdf_only)
    return lct_cube, lsm_cubes[0], lsm_cubes[1]


def save_ants(cube, filename, fill_value, netcdf_only):
    if not netcdf_only:
        save.ancil(cube, filename)
    save.netcdf(cube, filename, fill_value=fill_value)


def _get_parser():
    parser = ants.AntsArgParser(target_lsm=True, target_grid=True)
    parser.add_argument(
        "--transform-path", type=str, help="Filepath to crosswalk table.", required=True
    )
    msg = (
        "Output directory for the landseamask as derived by the the land cover type "
        "fraction field."
    )
    parser.add_argument("--landseamask-output-root", type=str, help=msg, required=False)
    return parser


if __name__ == "__main__":
    args = _get_parser().parse_args()
    if args.landseamask_output_root and args.target_lsm:
        raise ValueError(
            "Output landseamask filepath (landseamask_filepath) "
            "provided, along with an input landsea mask "
            "(target-lsm).  These are mutually exclusive "
            "arguments."
        )
    main(
        args.sources,
        args.transform_path,
        args.output,
        args.target_grid,
        args.target_lsm,
        args.landseamask_output_root,
        args.land_threshold,
        args.netcdf_only
    )

    