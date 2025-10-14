#!/bin/bash

module use /g/data/hr22/modulefiles
module load cylc7/24.03
module load python2-as-python

module use ~access/modules
module use /g/data/access/ngm/modules
module load cap/9.2

export CYLC_FLOW_VERSION=7.9.8
export ROSE_VERSION=2019.01.8
export CYLC_VERSION=$CYLC_FLOW_VERSION
export PATH=$CYLC_DIR/../24.03/bin:$PATH

export CYLC_SUITE_RUN_DIR="$HOME/cylc-run/u-dq487"
export CYLC_SUITE_DEF_PATH="${HOME}/cylc-run/u-dq487"
export CYLC_SUITE_DEF_PATH_ON_SUITE_HOST="/home/548/pag548/cylc-run/u-dq487"
export CYLC_SUITE_UUID="2f3b4304-7093-4c0f-b35a-c10faa698b7d"

# CYLC TASK ENVIRONMENT:
export CYLC_TASK_JOB="1/Flagship_ERA5to1km_12km-halo_ancil_lct/01"
export CYLC_TASK_NAMESPACE_HIERARCHY="root HOST_HPC HOST_ANTS_MPP HOST_ANTS_SERIAL ANCIL_ANTS ANCIL_LCT Flagship_ERA5to1km_12km-halo_ancil Flagship_ERA5to1km_12km-halo_ancil_lct"
export CYLC_TASK_DEPENDENCIES="Flagship_ERA5to1km_12km-halo_ancil_top.1 install_ants.1 install_ugants.1"
export CYLC_TASK_TRY_NUMBER=1
    
ULIMIT_SETTING="unlimited"
TMPDIR="$CYLC_TASK_WORK_DIR"
TIDS="/g/data/access/TIDS"
ROSE_TASK_N_JOBS="${PBS_NCPUS:-1}"
ANTS_NPROCESSES="1"
ANTS_VERSION="2.1.0"
ANTS_MODULE="ants/ug-2.1.0"
CONTRIB_APPS="$CYLC_SUITE_RUN_DIR/share/contrib_apps"
PATH_PREPEND="$CYLC_SUITE_RUN_DIR/share/fcm_make_ants/build/bin"
PYTHONPATH_PREPEND="$CYLC_SUITE_RUN_DIR/share/fcm_make_ants/build/lib"
ANCIL_MASTER=~/"$ROSE_SUITE_DIR_REL/share/data/etc/ancil_master_ants/"
ANCIL_PREPROC_PATH="$ROSE_DATA/etc/ants_preproc"
TRANSFORM_DIR="/g/data/access/TIDS/UM/ancil/data/transforms"
PATH="$ROSE_SUITE_DIR/site/nci-gadi:$PATH"
CONTRIB_PATH="${CYLC_SUITE_DEF_PATH}/src/contrib"
ANCIL_TARGET_PATH="$ROSE_DATA/ancils/Flagship_ERA5to1km/12km-halo"
HORIZ_GRID="$ROSE_DATA/ancils/Flagship_ERA5to1km/12km-halo/grid.nl"
VERT_LEV="$ROSE_DATA/ancils/Flagship_ERA5to1km/12km-halo/L70_80km"
CAPGRID="$ROSE_DATA/ancils/Flagship_ERA5to1km/12km-halo/grid.nl"
CAPHORIZGRID=""
VARIABLE="F"
ROSE_TASK_APP="ancil_lct"
export ULIMIT_SETTING TMPDIR TIDS ROSE_TASK_N_JOBS ANTS_NPROCESSES ANTS_VERSION ANTS_MODULE CONTRIB_APPS PATH_PREPEND PYTHONPATH_PREPEND ANCIL_MASTER ANCIL_PREPROC_PATH TRANSFORM_DIR PATH CONTRIB_PATH ANCIL_TARGET_PATH HORIZ_GRID VERT_LEV CAPGRID CAPHORIZGRID VARIABLE ROSE_TASK_APP

source=${ANCIL_MASTER}/vegetation/cover/cci/v3/vegetation_fraction.nc
transformpath=${TRANSFORM_DIR}/cci2jules_ra1.json
output_vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci_pre_c4
output_lsm=${ANCIL_TARGET_PATH}
target_grid=${HORIZ_GRID}
ANTS_CONFIG=${PWD}/rose-app-run.conf

~/cylc-run/u-dq487/share/fcm_make_ants/build/bin/ants-launch

#ants-launch ${CONTRIB_APPS}/LCT/ancil_lct.py ${source} --target-grid ${target_grid} --transform-path ${transformpath} -o ${output_vegfrac} --landseamask-output ${output_lsm} --ants-config ${ANTS_CONFIG}

# From rose-app-run.conf
#[ants_decomposition]
#x_split=20
#y_split=20