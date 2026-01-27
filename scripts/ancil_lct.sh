#!/bin/bash

#module use /g/data/hr22/modulefiles
#module load cylc7/24.03
#module load python2-as-python

module use ~access/modules
module use /g/data/access/ngm/modules
module load cap/9.2
#module load /g/data/access/ngm/modules/analysis3/23.07
module use /g/data/xp65/public/modules
module load conda/analysis3-25.10
set -x

export HOME="/home/548/pag548"

export CYLC_SUITE_RUN_DIR="/home/548/pag548/cylc-run/u-dg767"
export CYLC_SUITE_DEF_PATH="${HOME}/cylc-run/u-dg767"
export CYLC_SUITE_DEF_PATH_ON_SUITE_HOST="/home/548/pag548/cylc-run/u-dg767"
export CYLC_SUITE_UUID="192cf8a1-510c-4992-b2d8-5df71e94b1b9"
 
# CYLC TASK ENVIRONMENT:
export CYLC_TASK_JOB="1/Lismore_era5_ancil_lct/01"
export CYLC_TASK_NAMESPACE_HIERARCHY="root HOST_HPC HOST_ANTS_MPP HOST_ANTS_SERIAL ANCIL_ANTS ANCIL_LCT Lismore_era5_ancil Lismore_era5_ancil_lc
t"
export CYLC_TASK_DEPENDENCIES="Lismore_era5_ancil_top.1 ants_package_build.1"
export CYLC_TASK_TRY_NUMBER=1

    
ULIMIT_SETTING="unlimited"
TMPDIR="$CYLC_TASK_WORK_DIR"
TIDS="/g/data/access/TIDS"
ROSE_TASK_N_JOBS="${PBS_NCPUS:-1}"
ANTS_NPROCESSES="1"
ANTS_VERSION="2.1.0"
ANTS_MODULE="ants/ug-2.1.0"

#Manually add ROSE DIRS
ROSE_DATA="$CYLC_SUITE_RUN_DIR/share/data"

CONTRIB_APPS="$CYLC_SUITE_RUN_DIR/share/contrib_apps"
PATH_PREPEND="${HOME}/cylc-run/u-dq487/share/fcm_make_ants/build/bin"
PYTHONPATH_PREPEND="${HOME}/cylc-run/u-dq487/share/fcm_make_ants/build/lib"
#PATH_PREPEND=/g/data/gb02/public/ants_build/bin
#PYTHONPATH_PREPEND=/g/data/gb02/public/ants_build/lib/
ANCIL_MASTER="$CYLC_SUITE_RUN_DIR/share/data/etc/ancil_master_ants/"
ANCIL_PREPROC_PATH="$ROSE_DATA/etc/ants_preproc"
TRANSFORM_DIR="/g/data/access/TIDS/UM/ancil/data/transforms"
PATH="$ROSE_SUITE_DIR/site/nci-gadi:$PATH"
CONTRIB_PATH="${CYLC_SUITE_DEF_PATH}/src/contrib"
TRANSFORM_DIR="/g/data/access/TIDS/UM/ancil/data/transforms"
ANCIL_TARGET_PATH="$ROSE_DATA/ancils/Lismore/era5"
HORIZ_GRID="$ROSE_DATA/ancils/Lismore/era5/grid.nl"
VERT_LEV="$ROSE_DATA/ancils/Lismore/era5/L70_80km"
CAPGRID="$ROSE_DATA/ancils/Lismore/era5/grid.nl"
CAPHORIZGRID=""
VARIABLE="F"
ROSE_TASK_APP="ancil_lct"
export ULIMIT_SETTING TMPDIR TIDS ROSE_TASK_N_JOBS ANTS_NPROCESSES ANTS_VERSION ANTS_MODULE CONTRIB_APPS PATH_PREPEND PYTHONPATH_PREPEND ANCIL_MASTER ANCIL_PREPROC_PATH TRANSFORM_DIR PATH CONTRIB_PATH ANCIL_TARGET_PATH HORIZ_GRID VERT_LEV CAPGRID CAPHORIZGRID VARIABLE ROSE_TASK_APP

source=${ANCIL_MASTER}/vegetation/cover/cci/v3/vegetation_fraction.nc
transformpath=${TRANSFORM_DIR}/cci2jules_ra1.json
output_vegfrac=${ANCIL_TARGET_PATH}/qrparm.veg.frac_cci_pre_c4
output_lsm=${ANCIL_TARGET_PATH}
target_grid=${HORIZ_GRID}
target_lsm="${HOME}/code/UM_config_tools/scripts/dummy_cons.nc"
ANTS_CONFIG=${CYLC_SUITE_RUN_DIR}/work/1/Lismore_era5_${ROSE_TASK_APP}/rose-app-run.conf

# To load the executable
ANTS_LAUNCH="${HOME}/cylc-run/u-dq487/share/fcm_make_ants/build/bin/ants-launch"
#ANTS_LAUNCH=/g/data/gb02/public/ants_build/bin/ants-launch

export PATH=$PATH_PREPEND:$PATH
export PYTHONPATH=$PYTHONPATH_PREPEND:$PYTHONPATH_PREPEND/ants

echo "ANCIL_MASTER=$ANCIL_MASTER"
echo "PYTHONPATH=$PYTHONPATH"
# To launch the executable
#echo "ants-launch ${CONTRIB_APPS}/LCT/ancil_lct.py ${source} --target-lsm ${target_lsm} --transform-path ${transformpath} -o ${output_vegfrac} --landseamask-output ${output_lsm} --ants-config ${ANTS_CONFIG}"

${ANTS_LAUNCH} ${HOME}/code/UM_config_tools/scripts/ancil_lct.py ${source} --target-lsm ${target_lsm} --transform-path ${transformpath} -o ${output_vegfrac}  --ants-config ${ANTS_CONFIG}

# From rose-app-run.conf
#[ants_decomposition]
#x_split=20
#y_split=20