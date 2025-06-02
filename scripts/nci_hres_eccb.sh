#!/bin/bash

module use /g/data/hr22/modulefiles
module load cylc7/24.03
module load python2-as-python

export CYLC_VERSION=$CYLC_FLOW_VERSION
export CYLC_FLOW_VERSION=7.9.8
export ROSE_VERSION=2019.01.8
export PATH=$CYLC_DIR/../24.03/bin:$PATH

export CYLC_SUITE_RUN_DIR="/home/548/pag548/cylc-run/u-dq126"
export CYLC_SUITE_DEF_PATH="${HOME}/cylc-run/u-dq126"
export CYLC_SUITE_DEF_PATH_ON_SUITE_HOST="/home/548/pag548/cylc-run/u-dq126"
export CYLC_SUITE_UUID="d6eb1444-9256-4ede-8879-67901637acea"

# CYLC TASK ENVIRONMENT:
export CYLC_TASK_JOB="20220226T0000Z/nci_hres_eccb/02"
export CYLC_TASK_NAMESPACE_HIERARCHY="root HOST_HPC nci_hres_eccb"
export CYLC_TASK_DEPENDENCIES="ec_um_recon_000.20220226T0000Z install_cold_hpc.20220226T0000Z"
export CYLC_TASK_TRY_NUMBER=2

TMPDIR="$CYLC_TASK_WORK_DIR"
UMDIR="/g/data/access/projects/access/umdir"
TIDS="/g/data/access/TIDS"
ROSE_TASK_N_JOBS="${PBS_NCPUS:-4}"
ROSE_DATAC=/scratch/gb02/pag548/cylc-run/u-dq126/share/cycle/20220226T0000Z
ATMOS_LAUNCHER="mpirun"
RECON_LAUNCHER="mpirun"
PROJECT="gb02"
ECCBFILE="$ROSE_DATAC/ec/um/ec_cb000"
MASK="/scratch/gb02/pag548/cylc-run/u-dq124/share/data/ancils/Flagship_ERA5to1km/BARRA-R2-halo/qrparm.mask"
#START="$(rose date -c -f %Y%m%d%H%M)"
START=20220226
TYPE="barra"
ROSE_TASK_APP="nci_hres_eccb"
ROSE_DATA="/scratch/gb02/pag548/cylc-run/u-dq126/share/data/"
export TMPDIR UMDIR TIDS ROSE_TASK_N_JOBS ATMOS_LAUNCHER RECON_LAUNCHER PROJECT ECCBFILE MASK START TYPE ROSE_TASK_APP ROSE_DATA

module use /projects/access/modules
ulimit -s unlimited
module load openmpi/4.1.4
module load intel-mkl/2019.3.199
module load python2/2.7.16

module use /g/data/vk83/modules
module load conda/access-ram/2025.03.0
export BARRA_DIR=$ROSE_DATA/etc/barra_r2
export ERA_DIR=$ROSE_DATA/etc/era5_land

export PATH=/g/data/hr22/apps/cylc7/rose_2019.01.8/bin:/g/data/vk83/apps/conda/access-ram/2025.03.0/bin:/home/548/pag548/cylc-run/u-dq126/bin:/home/548/pag548/cylc-run/u-dq126/bin:/home/548/pag548/cylc-run/u-dq126/share/fcm_make_surf/build/bin:/apps/python2/2.7.16/bin:/apps/openmpi/wrapper/fortran:/apps/openmpi/wrapper:/apps/openmpi/4.1.4/bin:/g/data/hr22/apps/cylc7/cylc_7.9.9/../24.03/bin:/g/data/hr22/apps/cylc7/24.03/bin:/g/data/hr22/apps/mosrs-setup/2.0.1/bin:/g/data/hr22/apps/cylc7/cylc_7.9.9/bin:/opt/pbs/default/bin:/opt/nci/bin:/opt/bin:/opt/Modules/v4.3.0/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/local/pbs/bin:/local/pbs/bin


#hres_eccb --mask $MASK --file $ECCBFILE --start $START --type $TYPE
#python hres_eccb.py --mask $MASK --file $ECCBFILE --start $START --type $TYPE
