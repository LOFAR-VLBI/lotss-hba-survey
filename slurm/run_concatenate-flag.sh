#!/bin/bash 
#SBATCH -N 1                  # number of nodes
#SBATCH -c 1                 # number of cores 
#SBATCH --ntasks=1            # number of tasks
#SBATCH -t 72:00:00           # maximum run time in [HH:MM:SS] or [MM:SS] or [minutes]

## submit the job with OBSID as an argument
OBSID=${1}

#################################################################################
## Cluster specific directories to change
## PLEASE SEE slurm/add_these_to_bashrc.txt 

export VLBIDIR=${SOFTWAREDIR}/VLBI-cwl
export LINCDIR=${SOFTWAREDIR}/LINC
export FLOCSDIR=${SOFTWAREDIR}/flocs
export LOFARHELPERS=${SOFTWAREDIR}/lofar_helpers
export FACETSELFCAL=${SOFTWAREDIR}/lofar_facet_selfcal
BINDPATHS=${SOFTWAREDIR},${DATA_DIR},${WORK_DIR}

## FOR TOIL
export TOIL_SLURM_ARGS="${CLUSTER_OPTS} --export=ALL -t 24:00:00 -N 1 --ntasks=1"
export CWL_SINGULARITY_CACHE=${SOFTWAREDIR}/singularity
export TOIL_CHECK_ENV=True

#################################################################################
## IN GENERAL DO NOT TOUCH ANYTHING BELOW HERE

DATADIR=${DATA_DIR}/${OBSID}/setup
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
export LINC_DATA_ROOT=${LINCDIR}

# Pass along necessary variables to the container.
export APPTAINER_CACHEDIR=${SOFTWAREDIR}/singularity
export APPTAINER_TMPDIR=${APPTAINER_CACHEDIR}/tmp
export APPTAINER_PULLDIR=${APPTAINER_CACHEDIR}/pull
export APPTAINER_BIND=${BINDPATHS}
export APPTAINERENV_LINC_DATA_ROOT=${LINC_DATA_ROOT}
#### PATH: note that apptainer has a bug and does not use APPTAINERENV_PREPEND_PATH correctly
export SINGULARITYENV_PREPEND_PATH=${VLBIDIR}/scripts:${LINCDIR}/scripts
export APPTAINERENV_PYTHONPATH=${VLBIDIR}/scripts:${LINCDIR}/scripts:\$PYTHONPATH

## go to working directory
cd ${OUTDIR}

## list of measurement sets - THIS WILL NEED TO BE CHECKED
apptainer exec -B ${PWD},${BINDPATHS} --no-home ${LOFAR_SINGULARITY} python3 ${FLOCSDIR}/runners/create_ms_list.py VLBI concatenate-flag --linc ${LINCDIR} --aoflagger_memory_fraction 20 ${DATADIR}/ >> create_ms_list.log 2>&1

echo LINC starting
TMPID=`echo ${OBSID} | cut -d'/' -f 1`

ulimit -n 8192

toil-cwl-runner --no-read-only --singularity --bypass-file-store --jobStore=${JOBSTORE} --logFile=${OUTDIR}/job_output.txt --workDir=${WORKDIR} --outdir=${OUTPUT} --retryCount 0 --writeLogsFromAllJobs TRUE --writeLogs=${LOGSDIR} --tmp-outdir-prefix=${TMPD}/ --coordinationDir=${OUTPUT} --tmpdir-prefix=${TMPD}_interim/ --disableAutoDeployment True --preserve-environment ${APPTAINERENV_PYTHONPATH} ${SINGULARITYENV_PREPEND_PATH} ${APPTAINERENV_LINC_DATA_ROOT} ${APPTAINER_BIND} ${APPTAINER_PULLDIR} ${APPTAINER_TMPDIR} ${APPTAINER_CACHEDIR} --batchSystem slurm ${VLBIDIR}/workflows/concatenate-flag.cwl mslist_VLBI_concatenate-flag.json

if grep 'CWL run complete' ${OUTDIR}/job_output.txt
then 
	echo 'SUCCESS: Pipeline finished successfully' > ${OUTDIR}/finished.txt
else
	echo "**FAILURE**: Pipeline failed with exit status: ${?}" > ${OUTDIR}/finished.txt
fi

