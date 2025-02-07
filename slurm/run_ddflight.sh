#!/bin/bash 
#SBATCH --job-name=test_ddf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=54
#SBATCH --exclusive 
#SBATCH --time=72:00:00

## submit the job with OBSID as an argument
OBSID=${1}

echo "Starting up, field is " ${OBSID}
hostname
if ! test -d ${DATA_DIR}/${OBSID}/ddfpipeline
then
	mkdir ${DATA_DIR}/${OBSID}/ddfpipeline
fi
cd ${DATA_DIR}/${OBSID}/ddfpipeline
mkdir template
mv ${DATA_DIR}/${OBSID}/ddfsolutions/* template/
cp ${DDF_PIPELINE_INSTALL}/examples/tier1-rerun.cfg .
sed -i "s~\\\$\\\$~${BOOTSTRAP_DIR}~g" tier1-rerun.cfg
## move results here if measurement sets dont already exist
if compgen -G "L*MS" > /dev/null
then
	echo "files already present"
else
	mv ${DATA_DIR}/${OBSID}/HBA_target/results/* .
fi

DICOMODEL=template/image*DicoModel
singularity exec -B ${PWD},${SOFTWAREDIR} --no-home ${DDFPIPELINE_SINGULARITY} python3 ${SOFTWAREDIR}/lotss-hba-survey/fix_dicomodel.py ${DICOMODEL}

singularity exec -B ${PWD},${BOOTSTRAP_DIR} --no-home ${DDFPIPELINE_SINGULARITY} CleanSHM.py
singularity exec -B ${PWD},${BOOTSTRAP_DIR} --no-home ${DDFPIPELINE_SINGULARITY} make_mslists.py
singularity exec -B ${PWD},${BOOTSTRAP_DIR} --no-home ${DDFPIPELINE_SINGULARITY} pipeline.py tier1-rerun.cfg

if test -f image_full_ampphase_di_m.NS.app.restored.fits
then
	echo "SUCCESS: Pipeline finished successfully" > finished.txt
	## move solutions etc to ../ddfsolutions
	if ! test -d ${DATA_DIR}/${OBSID}/ddfsolutions
	then
		mkdir ${DATA_DIR}/${OBSID}/ddfsolutions
	fi
	mv SOLSDIR ${DATA_DIR}/${OBSID}/ddfsolutions
	## check whether the bootstrap will need to be applied
	FILES=(logs/KillMS-L*DIS2_full.log)
	TMP=`grep ' InCol' ${FILES[1]}`
	INCOL=`echo ${TMP} | cut -d' ' -f 7`
	if [ "${INCOL}" = "DATA" ]
	then
		echo "used bootstrapped model, no corrections necessary"
	else
		## will also need     
		#L*frequencies.txt (can be reconstructed if it doesn't exist)
		#L*crossmatch-results-2.npy
		cp L*frequencies.txt ${DATA_DIR}/${OBSID}/ddfsolutions
		cp L*crossmatch-results-2.npy ${DATA_DIR}/${OBSID}/ddfsolutions
	fi
else
	echo "Pipeline did not report finishing successfully. Please check processing" > finished.txt
fi

## copy files over to the processing directory where the monitor script checks for them
## and make sure that there will be the right things to check
if ! test -d ${DATA_DIR}/processing/${OBSID}
then
	mkdir ${DATA_DIR}/processing/${OBSID}
fi
cp finished.txt ${DATA_DIR}/processing/${OBSID}/
echo -e "ddflight.cwl Resolved \n\n\n\n\n\n\n\n\n\n" > ${DATA_DIR}/processing/${OBSID}/job_output.txt
cp big-mslist.txt ${DATA_DIR}/processing/${OBSID}/mslist-ddflight.json
sed -i 's~ms~ path~g' ${DATA_DIR}/processing/${OBSID}/mslist-ddflight.json
