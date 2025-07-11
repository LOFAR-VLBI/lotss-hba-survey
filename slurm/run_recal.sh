#!/bin/bash 
#SBATCH -N 1                  # number of nodes
#SBATCH -c 32                 # number of cores
#SBATCH --ntasks=1            # number of tasks
#SBATCH -t 6:00:00           # maximum run time in [HH:MM:SS] or [MM:SS] or [minutes]

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
export BINDPATHS=${SOFTWAREDIR},${DATA_DIR}

#################################################################################
## IN GENERAL DO NOT TOUCH ANYTHING BELOW HERE

## define the data directories
DATADIR=${DATA_DIR}/${OBSID}/selfcal
PROCDIR=${DATA_DIR}/processing
OUTDIR=${PROCDIR}/${OBSID}_${SLURM_ARRAY_TASK_ID}
LOGSDIR=${OUTDIR}/logs

LINENUM=$((${SLURM_ARRAY_TASK_ID}+1))
TMP=`sed -n "${LINENUM}p" ${DATA_DIR}/${OBSID}/../recalibration_list.csv`
TARGET=`echo ${TMP} | cut -d',' -f 1`
DDCAL=`echo ${TMP} | cut -d',' -f 2`

TARGETINMS=`ls -d ${DATADIR}/${TARGET}/ILTJ*`
DDCALH5=`ls -d ${DATADIR}/${DDCAL}/merged*h5`

mkdir -p ${OUTDIR}
cp -r ${TARGETINMS} ${OUTDIR}
cp -r ${DDCALH5} ${OUTDIR}

cd ${OUTDIR}
TARGETMS=`ls -d ILTJ*`
DDCALH5=`ls -d merged*h5`

apptainer exec -B ${PWD},${BINDPATHS} --no-home ${LOFAR_SINGULARITY} python3 ${FACETSELFCAL}/facetselfcal.py ${TARGETMS} --helperscriptspath ${FACETSELFCAL} --helperscriptspathh5merge ${LOFARHELPERS} --configpath ${VLBIDIR}/target_selfcal_config.txt --targetcalILT=tec --ncpu-max-DP3solve=56 --preapplyH5-list '["'${DDCALH5}'"]' > facet_selfcal.log 2>&1

: << END
## check if it finishes 
if compgen -G "merged_selfcalcycle*h5" > /dev/null; then
    ## find the last cycle to run to copy over
    TMP=`find -iname "merged_selfcalcycle*h5" | sort -n | tail -1`
    CYCLE=`echo ${TMP} | cut -d'_' -f 2`
	## clean up and write a finished.txt
	mkdir tmp
	mv * tmp/
	mv tmp/${TARGETMS} .
	#mv tmp/${TARGETMS}.copy .
	mv tmp/merged*${CYCLE}*.h5 .
	mv tmp/plotlosoto* .
	mv tmp/*png .
	mv tmp/*MFS-image.fits .
	mv tmp/*MFS-image-pb.fits .
	mv tmp/facet_selfcal.log .
	mv tmp/selfcal.log .
	rm -r tmp
	ILTJ=`ls -d ILTJ*ms | cut -d'_' -f 1`
	mkdir ${ILTJ}
	mv * ${ILTJ}/
	echo 'SUCCESS: Pipeline finished successfully' > ${OUTDIR}/finished.txt
	echo 'Resolved /fake/workflow/selfcal.cwl' > ${OUTDIR}/job_output.txt
elif cat facet_selfcal.log | grep -q "FAILED"; then
    ## it's crashed
    echo "**FAILURE**: Pipeline failed" > ${OUTDIR}/finished.txt
	echo 'Resolved /fake/workflow/selfcal.cwl' > ${OUTDIR}/job_output.txt    
else
    echo "**FAILURE**: Pipeline failed" > ${OUTDIR}/finished.txt
	echo 'Resolved /fake/workflow/selfcal.cwl' > ${OUTDIR}/job_output.txt
fi
END
