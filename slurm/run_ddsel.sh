#!/bin/bash 
#SBATCH -N 1                  # number of nodes
#SBATCH -c 1
#SBATCH --ntasks=1            # number of tasks
#SBATCH -t 48:00:00           # maximum run time in [HH:MM:SS] or [MM:SS] or [minutes]

## submit the job with OBSID as an argument
OBSID=${1}

#################################################################################
## Cluster specific directories to change
## PLEASE SEE slurm/add_these_to_bashrc.txt 

export LINCDIR=${SOFTWAREDIR}/LINC
export FLOCSDIR=${SOFTWAREDIR}/flocs
BINDPATHS=${SOFTWAREDIR},${DATA_DIR},${SCRATCH_DIR}

## FOR TOIL
export TOIL_SLURM_ARGS="${CLUSTER_OPTS} --export=ALL -t 24:00:00"
export CWL_SINGULARITY_CACHE=${SOFTWAREDIR}/singularity
export TOIL_CHECK_ENV=True

#################################################################################
## IN GENERAL DO NOT TOUCH ANYTHING BELOW HERE

## define the data directories
DATADIR=${DATA_DIR}/${OBSID}
PROCDIR=${DATA_DIR}/processing
OUTDIR=${PROCDIR}/${OBSID}
WORKDIR=${SCRATCH_DIR}/${OBSID}/workdir
OUTPUT=${OUTDIR}
JOBSTORE=${OUTDIR}/jobstore
TMPD=${WORKDIR}/tmp
LOGSDIR=${OUTDIR}/logs
mkdir -p ${WORKDIR}
mkdir -p ${TMPD}
mkdir -p ${TMPD}_interim
mkdir -p ${LOGSDIR}

## location of LINC
LINC_DATA_ROOT=${LINCDIR}

# Pass along necessary variables to the container.
CONTAINERSTR=$(singularity --version)
if [[ "$CONTAINERSTR" == *"apptainer"* ]]; then
    export APPTAINER_CACHEDIR=${SOFTWAREDIR}/singularity
    export APPTAINER_TMPDIR=${APPTAINER_CACHEDIR}/tmp
    export APPTAINER_PULLDIR=${APPTAINER_CACHEDIR}/pull
    export APPTAINER_BIND=${BINDPATHS}
    export APPTAINERENV_LINC_DATA_ROOT=${LINC_DATA_ROOT}
    export APPTAINERENV_LOGSDIR=${LOGSDIR}
    export APPTAINERENV_TMPDIR=${TMPDIR}
    export SINGULARITYENV_PREPEND_PATH=${LINC_DATA_ROOT}/scripts
    #export APPTAINERENV_PYTHONPATH=${LINC_DATA_ROOT}/scripts:\$PYTHONPATH
else
    export SINGULARITYENV_LINC_DATA_ROOT=${LINC_DATA_ROOT}
    export SINGULARITYENV_LOGSDIR=${LOGSDIR}
    export SINGULARITYENV_TMPDIR=${TMPDIR}
    export SINGULARITYENV_PREPEND_PATH=${LINC_DATA_ROOT}/scripts
fi

## go to working directory
cd ${OUTDIR}

## pipeline input 
python3 ${SOFTWAREDIR}/lotss-hba-survey/lotsshr_ddsel_json.py ${DATADIR}/selfcal

ulimit -n 8192

toil-cwl-runner --no-read-only --singularity --bypass-file-store --moveExports --jobStore=${JOBSTORE} --logFile=${OUTDIR}/job_output.txt --workDir=${WORKDIR} --outdir=${OUTPUT} --retryCount=0 --writeLogsFromAllJobs=True --writeLogs=${LOGSDIR} --tmp-outdir-prefix=${TMPD}/ --tmpdir-prefix=${TMPD}_interim/ --disableAutoDeployment=True --preserve-environment ${APPTAINERENV_PYTHONPATH} ${SINGULARITYENV_PREPEND_PATH} ${APPTAINERENV_LINC_DATA_ROOT} ${APPTAINER_BIND} ${APPTAINER_PULLDIR} ${APPTAINER_TMPDIR} ${APPTAINER_CACHEDIR} --batchSystem=slurm ${SOFTWAREDIR}/VLBI-cwl/workflows/subworkflows/ddcal_pre_selection.cwl config.json 

