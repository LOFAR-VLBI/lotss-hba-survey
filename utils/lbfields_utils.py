#!/usr/bin/python

# Utilities for monitor_lbfields and friends

from time import sleep
import datetime
from surveys_db import SurveysDB, tag_field, get_cluster

import os
import threading
import glob
import requests
import stager_access
from rclone import RClone   ## DO NOT pip3 install --user python-rclone -- use https://raw.githubusercontent.com/mhardcastle/ddf-pipeline/master/utils/rclone.py
from download_file import download_file ## in ddf-pipeline/utils
#import progress_bar
from sdr_wrapper import SDR
from reprocessing_utils import do_sdr_and_rclone_download, do_rclone_download
from tasklist import set_task_list
from calibrator_utils import get_linc, download_ddfpipeline_solutions, download_field_calibrators, unpack_calibrator_sols, compare_solutions
import numpy as np


def update_status(name,status,stage_id=None,time=None,workdir=None,av=None,survey=None):
    # adapted from surveys_db
    # utility function to just update the status of a field
    # name can be None (work it out from cwd), or string (field name)
    with SurveysDB(survey=survey) as sdb:
        idd=sdb.db_get('lb_fields',name)
        if idd is None:
          raise RuntimeError('Unable to find database entry for field "%s".' % name)
        idd['status']=status
        tag_field(sdb,idd,workdir=workdir)
        if time is not None and idd[time] is None:
            idd[time]=datetime.datetime.now()
        if stage_id is not None:
            idd['staging_id']=stage_id
        sdb.db_set('lb_fields',idd)

##############################
## job management

def run_apptainer( command ):
    singularity = os.getenv('LOFAR_SNGULARITY')
    bindpaths = ','.join([os.getenv('SOFTWAREDIR'),os.getenv('DATA_DIR')])
    os.system( 'apptainer exec -B {:s} --no-home {:s} {:s}'.format( bindpaths, singularity, command ) )

def restart_toil_job( field, obsid, workflow ):
    softwaredir = os.getenv('SOFTWAREDIR')
    slurmscript = os.path.join( softwaredir, 'autoPILOT/slurm', 'run_{:s}.sh'.format( workflow.replace('HBA_','' ) ) )
    ## read the file
    with open( slurmscript, 'r' ) as f:
        lines = f.readlines()
    for i in np.arange(0,len(lines)):
        if 'toil-cwl-runner' in lines[i]:
            lines[i] = lines[i].replace('--singularity','--singularity --restart') 
        if 'create_ms_list.py' in lines[i]:
            lines[i] = '## '+lines[i]
    with open( 'restart_slurm.sh','w' ) as f:
        for line in lines:
            f.write(line)    
    cluster_opts = os.getenv('CLUSTER_OPTS')
    if workflow == 'split-directions':
        command = chunk_imagecat( fieldobsid )
    else:
        command = "sbatch -J {:s} {:s} restart_slurm.sh {:s}/{:s}".format(field, cluster_opts, field, obsid)
    update_status(field,'Queued')
    os.system( 'rm {:s}'.format( os.path.join( os.getenv('DATA_DIR'),'processing',field,obsid,'finished.txt') ) )
    os.system( 'rm {:s}'.format( os.path.join( os.getenv('DATA_DIR'),'processing',field,obsid,'job_output.txt') ) )
    if os.system(command):
        update_status(field,"Submission failed")


##############################
## do things by obsid

def get_obsids( name, survey=None ):
    with SurveysDB(survey=survey) as sdb:
        sdb.execute('select * from observations where field="'+name+'"')
        fld = sdb.cur.fetchall()
    obsids = [ val['id'] for val in fld ]
    return(obsids)

def get_local_obsid( name ):
    basedir = os.getenv('DATA_DIR')
    fielddir = os.path.join(basedir, name)
    fieldfiles = glob.glob(os.path.join(fielddir,'*'))
    fielddirs = []
    for ff in fieldfiles:
        tmp = os.path.splitext(ff)
        if tmp[1] == '':
            tmp2 = os.path.basename(tmp[0])
            if tmp2.isdigit():
                fielddirs.append(tmp2)
    return(fielddirs)

##############################
## finding and checking solutions 

def collect_solutions_lhr( caldir ):
    caldir = caldir.rstrip('/')
    survey=None
    obsid = os.path.basename(caldir)
    namedir = os.path.dirname(caldir)
    name = os.path.basename(namedir)
    tasklist = []

    update_status( name,'Solutions' )
    success = True

    with SurveysDB(survey=survey) as sdb:
        sdb.execute('select * from observations where id="'+obsid+'"')
        fld = sdb.cur.fetchall()
    calibrator_id = fld[0]['calibrator_id']

    ## get the calibrator solutions
    result = download_field_calibrators(obsid,caldir)
    if len(result[int(obsid)]) == 0:
        ## need to re-run calibrator .... shouldn't ever be in this situation but here fore completeness
        tasklist.append('calibrator')
        tasklist.append('target_VLBI')
        tasklist.append('delay-calibration')
        tasklist.append('delay')
        tasklist.append('split-directions')
        tasklist.append('selfcal')
    else:
        solutions = unpack_calibrator_sols(caldir,result)
        if len(solutions) >= 1:
            print('One or more calibrator found, comparing solutions ...')
            best_sols = compare_solutions(solutions)
            print('Best solutions are {:s}, cleaning up others.'.format(best_sols[0]))
            os.system('cp {:s} {:s}/LINC-cal_solutions.h5'.format(best_sols[0],os.path.dirname(best_sols[0])))
            for sol in solutions:
                os.system('rm -r {:s}/{:s}*'.format(os.path.dirname(best_sols[0]),os.path.basename(sol).split('_')[0]))
            tasklist.append('target_VLBI')
            tasklist.append('delay-calibration')
            tasklist.append('delay')
            tasklist.append('split-directions')
            tasklist.append('selfcal')
    if success:
        ## set the task list in the lb_operations table
        set_task_list(obsid,tasklist)
        ## set the status back to staging if no runtime errors
        update_status( name, 'Staging' )

def get_calibrators( field ):
    ## get obsids
    obsids = get_obsids( field )
    ## stage the data
    for obsid in obsids:
        ss = 'flocs-lta search-id --sasid {:s} --freq_end=168. --get-surls --stage-products calibrator'.format(str(obsid))
        os.system( ss )
        os.system('rm 20*log')

    all_uris = []
    for obsid in obsids:
        with open( 'srms_{:s}.txt'.format(obsid) ) as f:
            lines = f.readlines()
        uris = [ line.rstrip('\n') for line in lines ]
        all_uris = all_uris + uris

    ## submit as a single request to get the stage id

    stage_status = stager_access.get_status(s)



    ## download the data

    ## run LINC calibrator
    flocs_common_options = "--record-toil-stats --scheduler slurm --slurm-queue {:s} --slurm-account {:s} --runner toil --rundir {:s} --outdir {:s} ".format(os.getenv('SLURM_QUEUES'),os.getenv('SLURM_ACCOUNT'), rundir,outdir)
    calibrator_directory = os.path.join(os.getenv('DATA_DIR'),fieldobsid,'calibrator')  ## doesn't actually exist for lotss-hr because we always just have calibrator solutions already
    calibrator_options = '--slurm-time 24:00:00 --save-raw-solutions {:s} {:s}'.format(calibrator_directory)
    command = 'flocs-run linc calibrator '+flocs_common_options+calibrator_options





def run_task( fieldobsid, task ):

    field = fieldobsid.split('/')[0]
    if os.getenv("SCRATCH_DIR") is not None: #This should work, untested, fingers crossed
        rundir = os.path.join(os.getenv('SCRATCH_DIR'),fieldobsid,'rundir')
    else:
        rundir = os.path.join(os.getenv('DATA_DIR'),'processing',fieldobsid,'rundir')
    outdir = os.path.join(os.getenv('DATA_DIR'),'processing',fieldobsid)
    os.makedirs(rundir)
    os.makedirs(outdir)

    fielddir = os.path.join(os.getenv('DATA_DIR'),fieldobsid)

    flocs_common_options = "--record-toil-stats --scheduler slurm --slurm-queue {:s} --slurm-account {:s} --runner toil --rundir {:s} --outdir {:s} ".format(os.getenv('SLURM_QUEUES'),os.getenv('SLURM_ACCOUNT'), rundir,outdir)
    ## Frits to look at this and check that this is sensible (as well as stuff below)

    if task == 'calibrator':
        ## need to stage and download calibrators
        calibrator_directory = os.path.join(os.getenv('DATA_DIR'),fieldobsid,'calibrator')  ## doesn't actually exist for lotss-hr because we always just have calibrator solutions already
        os.mkdirs( calibrator_directory )

        ## use obsid and flocs-lta to stage and download calibrators - Frits

        calibrator_dirs = glob.glob( calibrator_directory + '/*' )
        cal_success = 0
        for calibrator_dir in calibrator_dirs:
            calibrator_options = '--slurm-time 24:00:00 --save-raw-solutions {:s}'.format(calibrator_directory)
            command = 'flocs-run linc calibrator '+flocs_common_options+calibrator_options
            ## run linc calibrator
            ## now tar the results
            ## and upload to spider
            ## and move the results locally
            tarfiles = glob.glob( 'somewhere' ) ### UPDATE!
            success = upload_to_spider( tarfiles, spider_location='disk/surveys/' )
            cal_success = cal_success + success 

    elif task == 'target_VLBI':
        ## flocs-run linc target
        cal_solutions = os.path.join( fielddir, 'LINC-cal_solutions.h5' )
        ## get skymodel -- update to get LoTSS skymodel
        target_skymodel = os.path.join( fielddir, 'target.skymodel' )
        mslist = glob.glob( os.path.join( fielddir, '*.MS' ) )
        ss = "python3 {:s}/autoPILOT/download_skymodel_target.py --Radius 5. --Source LOTSS --DoDownload True --targetname={:s} --fluxlimit 0.01 {:s} {:s}".format( softwaredir, field, mslist[0], target_skymodel )
        run_apptainer( ss )
        target_options = "--slurm-time 48:00:00 --output-fullres-data --min-unflagged-fraction 0.05 --offline-workers --target_skymodel {:s} --cal-solutions {:s} {:s}".format(cal_solutions,fielddir)
        commmand = 'flocs-run linc target '+flocs_common_options+target_options
    elif task == 'delay-calibration':
        datadir = os.path.join( fielddir, 'HBA_target_VLBI', 'results' )
        delaycal_options = "--slurm-time 48:00:00 --delay-calibrator {:s} --ms-suffix dp3concat {:s}".format(delay_catalogue,datadir)
        command = 'flocs-run vlbi delay-calibration '+flocs_common_options+delaycal_options
    elif task == 'delay':
        update_status('DelayCheck')
    elif task == 'split-directions':
        ## NEED TO SORT OUT SELFCAL
        datadir = os.path.join( fielddir, 'HBA_target_VLBI', 'results' )
        delaycal_solutions = glob.glob( os.path.join( fielddir, delay-calibration, '*verified.h5' ) )[0]
        nchunks = chunk_imagecat( fieldobsid )
        commands = []
        for i in range(nchunks):
            image_catalogue = os.path.join(os.getenv('DATA_DIR'),field,'image_catalogue_{:s}.csv'.format(str(i+1)))
            split_options = "--slurm-time 48:00:00 --delay-solset {:s} --source-catalogue {:s} --no-do-selfcal --ms-suffix dp3concat {:s}".format(delaycal_solutions,image_catalogue,datadir)
            command = 'flocs-run vlbi split-directions '+flocs_common_options+split_options
            commands.append(command)
            ## sort out how to do multiprocessing
    elif task == 'selfcal':
        datadir = os.path.join( os.getenv('DATA_DIR'), fieldobsid, 'split-directions' )
        msfiles = glob.glob( os.path.join( datadir, 'ILTJ*' ) )
        with open( os.path.join( os.path.dirname(msfiles[0]), 'targetlist.txt' ), 'w' ) as f:
            for msfile in msfiles:
                f.write('{:s}\n'.format(msfile) )
            with open('selfcal_{:s}.sh'.format(field),'w') as f:
                f.write('#!/bin/bash -l\n\n')
                f.write('#SBATCH --ntasks=1\n')
                f.write('#SBATCH --cpus-per-task=32\n')
                f.write('#SBATCH --job-name=selfcal\n')
                f.write('#SBATCH -t 6:00:00\n\n')
                f.write('DATADIR={:s}\n'.format(datadir))
                f.write('OUTDIR={:s}/selfcal_${SLURM_ARRAY_TASK_ID}\n'.format(outdir))
                f.write('TARGETINMS=`sed -n "${SLURM_ARRAY_TASK_ID}p" ${DATADIR}/targetlist.txt`\n')
                f.write('mkdir -p ${OUTDIR}\n')
                f.write('mv ${TARGETINMS ${OUTDIR}\n')
                f.write('cd ${OUTDIR}\n')
                f.write('TARGETMS=`ls -d ILTJ*`\n')
                f.write("apptainer exec -B {:s},{:s} --no-home {:s} facetselfcal --configpath ${VLBIDIR}/target_selfcal_config.txt --targetcalILT=tec --ncpu-max-DP3solve=32 > facet_selfcal.log 2>&1") 
            os.system('sbatch {:s} --array=1-{:s}%10 selfcal_{:s}.sh'.format(os.getenv('CLUSTER_OPTS'),str(len(msfiles))) )


    os.system(command)

    
        ## an example of the outdir, which here is J090630+693030
        #(pyenv) [lofarvlbi-fsweijen@ui-01 J090630+693030]$ ls LINC_calibrator_L749280_2025_10_31-08_34_26/
        #jobstore  logs_LINC_calibrator.tar  results_LINC_calibrator


def collect_solutions( caldir ):
    survey=None
    obsid = os.path.basename(caldir)
    namedir = os.path.dirname(caldir)
    name = os.path.basename(namedir)
    tasklist = []

    update_status( name,'Solutions' )
    success = True

    with SurveysDB(survey=survey) as sdb:
        sdb.execute('select * from observations where id="'+obsid+'"')
        fld = sdb.cur.fetchall()
    calibrator_id = fld[0]['calibrator_id']

    ## check if linc/prefactor 3 has been run - for calibrator solutions - which are contained in LINC-target_solutions.h5
    linc_check, macname = get_linc( obsid, caldir )

    ## download previous ddfpipeline for re-running
    templatedir = os.path.join(caldir,'ddfpipeline/template')
    os.makedirs(templatedir)
    try:
        download_ddfpipeline_solutions(name,templatedir,ddflight=True)
    except RuntimeError:
        success = False

    if linc_check: 
        ## rename linc solutions so target_VLBI will pick them up -- but what will happen because the target solutions already exist?
        os.system('cp {:s}/LINC-target_solutions.h5 {:s}/LINC-cal_solutions.h5'.format(caldir,caldir))
        tasklist.append('target_VLBI')
        tasklist.append('ddflight')
        tasklist.append('concatenate-flag')
        tasklist.append('process-ddf')
        tasklist.append('phaseup-concat')
        tasklist.append('delay')
        tasklist.append('split')
        tasklist.append('selfcal') 
        ## some more widefield things .... 
    else:
        print('valid LINC solutions not found. Checking lb_calibrators.')
        ## linc is not good
        result = download_field_calibrators(obsid,caldir)
        solutions = unpack_calibrator_sols(caldir,result)
        if len(solutions) >= 1:
            print('One or more calibrator found, comparing solutions ...')
            best_sols = compare_solutions(solutions)
            print('Best solutions are {:s}, cleaning up others.'.format(best_sols[0]))
            os.system('cp {:s} {:s}/LINC-cal_solutions.h5'.format(best_sols[0],os.path.dirname(best_sols[0])))
            for sol in solutions:
                os.system('rm -r {:s}/{:s}*'.format(os.path.dirname(best_sols[0]),os.path.basename(sol).split('_')[0]))
            tasklist.append('target_VLBI')
            tasklist.append('ddflight')
            tasklist.append('concatenate-flag')
            tasklist.append('process-ddf')
            tasklist.append('phaseup-concat')
            tasklist.append('delay')
            tasklist.append('split')
            tasklist.append('selfcal')
            ## some more widefield things ....
        else:
            ## need to re-run calibrator .... shouldn't ever be in this situation!
            success = False
            tasklist.append('calibrator')
            tasklist.append('target_VLBI')
            tasklist.append('ddflight')
            tasklist.append('concatenate-flag')
            tasklist.append('process-ddf')
            tasklist.append('phaseup-concat')
            tasklist.append('delay')
            tasklist.append('split')
            tasklist.append('selfcal')
    if success:
        ## set the task list in the lb_operations table
        set_task_list(obsid,tasklist)
        ## set the status back to staging if no runtime errors
        update_status( name, 'Staging' )

##############################
## staging

def stage_field( name, survey=None ):
    ## get the observation ids 
    with SurveysDB(survey=survey) as sdb:
        sdb.cur.execute('select * from observations where field=%s',(name,))
        obs = sdb.cur.fetchall()
    ## get the observation ids which are not failed
    obsids = []
    for nobs in range(len(obs)):
        if obs[nobs]['status'] in ['Archived','DI_Processed']:
            obsids.append(obs[nobs]['id'])
    ## the field directory
    fielddir = os.path.join(str(os.getenv('DATA_DIR')),str(name))
    caldirs = []
    csrmfiles = []
    for obsid in obsids:
        caldir = os.path.join(fielddir, str(obsid))
        srmfile = 'srms_{:s}.txt'.format(str(obsid))
        csrmfile = os.path.join(caldir,srmfile)
        csrmfiles.append(csrmfile)
        caldirs.append(caldir)
        if not os.path.exists(caldir):
            os.makedirs(caldir)
            ss = 'flocs-lta search-id --sasid {:s} --freq_end=168. --get-surls'.format(str(obsid))
            os.system( ss )
            os.system('mv {:s} {:s}'.format(srmfile, csrmfile))
            os.system('rm 20*log')
        else:
            if os.path.exists(csrmfile):
                with open(csrmfile,'r') as f:
                    lines = f.readlines()
                uris = [ line.rstrip('\n') for line in lines ]
                tarfiles = glob.glob(os.path.join(caldir,'*.tar'))
                trfs = [ val.split('/')[-1] for val in tarfiles ]
                ## remove the already downloaded files from uris_to_stage ...
                idxs = [ i for i,val in enumerate(uris) if val.split('/')[-1] in trfs ]
                for idx in idxs:
                    uris.pop(idx)
                with open(csrmfile,'w') as f:
                    for uri in uris:
                        f.write(uri)
                        f.write('\n')
            else:
                print('Something bad happened for {:s}'.format(str(obsid)))
    ## now collect everything that needs to be staged
    uris = []
    for csrmfile in csrmfiles:
        with open(csrmfile,'r') as f:
            lines = f.readlines()
        for line in lines:
            uris.append(line.rstrip('\n'))

    stage_id = stager_access.stage(uris)
    update_status(name,'Staging', stage_id = stage_id )
    return(caldirs)

##############################
## downloading

def make_macaroon_file( obsid, site, url, macdict ):
    for mac in macdict:
        try: 
            bearer_token = mac[site]
        except:
            print('Token not found for {:s}'.format(site))
    macname = os.path.join( os.getenv('PWD'), obsid + '.conf' )
    with open( macname, 'w' ) as f:
        f.write( 'type = webdav\n' )
        f.write( 'vendor = other\n' )
        f.write( 'url = {:s}\n'.format(url) )
        f.write( 'bearer_token = {:s}'.format(bearer_token) )
    return(macname)

def do_download( name ):
    basedir = os.getenv('DATA_DIR')
    update_status(name,'Downloading')
    ## get the staging id from the surveys database
    with SurveysDB(readonly=True) as sdb:
        idd=sdb.db_get('lb_fields',name)
        stage_id = idd['staging_id']
    surls = stager_access.get_surls_online(stage_id)
    macaroons = stager_access.get_macaroons(stage_id)
    
    if len(surls) == 0:
        print('Something went wrong with the download for {:s} (staging id {:s})'.format(name,str(stage_id)))
        update_status(name,'Download failed')
    else:
        tmp_obsids = [ val.split('/')[-2] for val in surls ]
        obsids = np.unique(tmp_obsids)
        all_tarfiles = []
        for obsid in obsids:
            tmp = os.path.join(str(os.getenv('DATA_DIR')),str(name))
            caldir = os.path.join(tmp,obsid)
            obsid_surls = [ surl for surl in surls if obsid in surl ]
            ## get project
            project = surls[0].split('projects/')[-1].split('/')[0]
            url = ':'.join(surls[0].split(':')[0:2]) + ':' + surls[0].split(':')[-1].split('/')[0]

            if 'juelich' in obsid_surls[0]:
                mac_name = make_macaroon_file( obsid, 'Juelich', url, macaroons )
            elif 'psnc' in obsid_surls[0]:
                mac_name = make_macaroon_file( obsid, 'Poznan', url, macaroons )
            if 'surf' in obsid_surls[0]:
                mac_name = make_macaroon_file( obsid, 'SURF', url, macaroons )

            if os.path.exists(mac_name):
                files = [ project+val.split(project)[-1] for val in obsid_surls ]
                rc = RClone( mac_name, debug=True)
                rc.get_remote()
                for f in files:
                    d = rc.execute(['-P','--no-check-certificate=true','copy',rc.remote + f]+[caldir]) 
                if d['err'] or d['code']!=0:
                    update_status(field,'rclone failed')
                    print('Rclone failed for field {:s}'.format(field))
                os.system('rm {:s}'.format(mac_name))
                ## check that everything was downloaded
                tarfiles = check_tarfiles( caldir )
                all_tarfiles = all_tarfiles + tarfiles
                if len(all_tarfiles) == len(surls):
                    print('Download successful for {:s}'.format(name) )
                    update_status(name,'Downloaded',stage_id=0)
            else:
                print('Failed to create macaroon')
                update_status(name,'Macaroon failed')

def check_tarfiles( caldir ):
    trfs = glob.glob(os.path.join(caldir,'*tar'))
    for trf in trfs:
        os.system( 'tar -tvf {:s} > tmp.txt 2>&1'.format(trf) )
        with open( 'tmp.txt', 'r' ) as f:
            lines = f.readlines()
        if 'tar: Unexpected EOF in archive\n' in lines:
            os.system( 'rm {:s}'.format(trf) )
        if 'tar: Exiting with failure status due to previous errors\n' in lines:
            os.system( 'rm {:s}'.format(trf) )
    os.system('rm tmp.txt')    
    trfs = glob.glob(os.path.join(caldir,'*tar'))
    return(trfs)

##############################
## unpacking

def dysco_compress(caldir,msfile):
    msfile = os.path.join(caldir,msfile)
    success=True
    with open(os.path.join(caldir,'dysco_compress.parset'),'w') as f:
        f.write('msin={:s}\n'.format(msfile))
        f.write('msin.datacolumn=DATA\n')
        f.write('msout={:s}.tmp\n'.format(msfile))
        f.write('msout.datacolumn=DATA\n')
        f.write('msout.storagemanager=dysco\n')
        f.write('steps=[count]')
    sing_img = os.getenv('LOFAR_SINGULARITY')
    os.system('singularity exec -B {:s} {:s} DP3 {:s}'.format(os.getcwd(),sing_img,os.path.join(caldir,'dysco_compress.parset')))
    if os.path.exists('{:s}.tmp'.format(msfile)):
        os.system('rm -r {:s}'.format(msfile))
        os.system('mv {:s}.tmp {:s}'.format(msfile,msfile))
    else:
        print('something went wrong with dysco compression for {:s}'.format(msfile))
        success=False
    return(success)

def dysco_compress_job(caldir):
    success=True
    os.system('ls -d {:s}/*.MS > {:s}/myfiles.txt'.format(caldir,caldir))
    file_number = len(open("{:s}/myfiles.txt".format(caldir), "r").readlines())
    command = 'sbatch -W --array=1-{:n}%5 {:s} {:s}/autoPILOT/slurm/dysco.sh {:s}'.format(file_number,os.getenv('CLUSTER_OPTS'),os.getenv('SOFTWAREDIR'),caldir)
    if os.system(command):
        print("Something went wrong with the dysco compression job!")
        success = False
    os.system('rm {:s}/myfiles.txt'.format(caldir))
    return success

def do_unpack(field):
    update_status(field,'Unpacking')
    caldir = os.path.join(str(os.getenv('DATA_DIR')),field)
    obsdirs = glob.glob(os.path.join(caldir,'*'))
    tmp = [ val for val in obsdirs if os.path.isdir(val) ]
    obsdirs = tmp
    for obsdir in obsdirs:
        ## get the tarfiles
        tarfiles = glob.glob(os.path.join(obsdir,'*tar'))
        for trf in tarfiles:
            with open('unpack_{:s}.sh'.format(field),'w') as f:
                f.write('#!/bin/bash -l\n\n')
                f.write('#SBATCH --ntasks=1\n')
                f.write('#SBATCH --cpus-per-task=1\n')
                f.write('#SBATCH --job-name=untar\n')
                f.write('#SBATCH -t 4:00:00\n\n')
                f.write('OUTDIR={:s}\n'.format(os.path.dirname(trf)))
                f.write('cd ${OUTDIR}\n')
                f.write("apptainer exec -B {:s},{:s} --no-home {:s} python3 {:s}/autoPILOT/unpack_and_dysco_compress.py {:s}".format(os.getenv('SOFTWAREDIR'),os.getenv('DATA_DIR'),os.getenv('LOFAR_SINGULARITY'),os.getenv('SOFTWAREDIR'),trf))
            os.system('sbatch {:s} -W unpack_{:s}.sh'.format(os.getenv('CLUSTER_OPTS'),field) )
    ## check that everything unpacked
    success = 0    
    for obsdir in obsdirs:
        failed = glob.glob(os.path.join(obsdir,'failed*txt'))
        if len(failed) == 0:
            success += 1
    if success == len(obsdirs):
        update_status(field,'Unpacked')
    else:
        update_status(field,'Unpack failed')    

##############################
## get LINC input

def get_linc_inputs( field, obsid ):
    datadir = os.path.join( os.getenv('DATA_DIR'), field, obsid )
    softwaredir = os.getenv('SOFTWAREDIR')
    mslist = glob.glob( os.path.join( datadir, '*.MS' ) )
    singularity_img = os.getenv('LOFAR_SINGULARITY')
    ## download TGSS skymodel
    skymodel = os.path.join( datadir, 'target.skymodel' )
    cmd = "apptainer exec -B {:s},{:s} --no-home {:s} python3 {:s}/LINC/scripts/download_skymodel_target.py --targetname={:s} {:s} {:s}".format( os.getcwd(), softwaredir, singularity_img, softwaredir, field, mslist[0], skymodel )
    if os.system(cmd):
        update_status(field,"TGSS failed")
    #Download IONEX
    ionexpath = datadir
    cal_solutions = os.path.join( datadir, 'LINC-cal_solutions.h5' )
    cmd = "apptainer exec -B {:s},{:s} --no-home {:s} spinifex get_rm_h5parm_from_ms {:s} -o {:s} --solset-name target --soltab-name spinifex".format( os.getcwd(), softwaredir, singularity_img, mslist[0 ], cal_solutions)
    cc = os.system(cmd)
    if cc == 256:
        os.system(cmd.replace('http://ftp.aiub.unibe.ch/CODE/','http://chapman.upc.es/'))
    #if os.system(cmd):
    #    update_status(field,"IONEX failed")

##############################
## split directions

def chunk_imagecat( fieldobsid, numdirs=10, catname='image_catalogue.csv', nchunkspertime=2 ):
    basedir = os.getenv('DATA_DIR')
    name = fieldobsid.split('/')[0]
    obsid = fieldobsid.split('/')[1]
    catfile = os.path.join( basedir, '{:s}/{:s}'.format(name,catname) )

    with open( catfile, 'r' ) as f:
        lines = f.readlines()
    header = lines[0]
    lines = lines[1:]

    nchunks = int(np.ceil(len(lines)/numdirs))
    chunk = 1
    for i in np.arange(nchunks):
        j = 10*i
        if i < nchunks - 1:
            mylines = lines[j:j+numdirs]
        else:
            mylines = lines[j:]
        with open( catfile.replace('.csv','_{:s}.csv'.format(str(chunk))), 'w' ) as f:
            f.write(header)
            for myline in mylines:
                f.write(myline)
        chunk = chunk + 1
    
    #cmd = 'sbatch -J split {:s} --array=1-{:s}%{:s} {:s}/autoPILOT/slurm/run_split-directions.sh {:s}'.format(os.getenv('CLUSTER_OPTS'),str(nchunks),str(nchunkspertime),os.getenv('SOFTWAREDIR'),fieldobsid)
    return( nchunks )



##############################
## verifying

def get_workflow_obsid(outdir):
    ## get the workflow that was run
    with open(os.path.join(outdir,'job_output.txt'),'r') as f:
        for line in f:
            if 'Resolved' in line:
                break
            elif 'Workflow' in line:
                if '.cwl' in line:
                    break
    tmp = line.split('.cwl')
    workflow = os.path.basename(tmp[0])
    obsid = os.path.basename(outdir).split("_")[0]
    return(workflow,obsid)

def check_field(field):
    field_obsids = get_local_obsid(field)
    for tmp_obsid in field_obsids:
        tmp_fieldobsid = '{:s}/{:s}'.format(field,tmp_obsid)
        ## check for presence of processing directories
        tmp_outdir = glob.glob(os.path.join(os.getenv('DATA_DIR'),'processing','{:s}*'.format(tmp_fieldobsid)))
        if len(tmp_outdir) > 0:
            ## check if directory is not empty
            contents = os.listdir(tmp_outdir[0])
            if len(contents) > 0:
                obsid = tmp_obsid
    fieldobsid = '{:s}/{:s}'.format(field,obsid)
    procdir = os.path.join(str(os.getenv('DATA_DIR')),'processing')
    outdirs = glob.glob(os.path.join(procdir,'{:s}*'.format(fieldobsid)))
    ## frits will change flocs-run to put log in output directory
    finished = glob.glob(os.path.join(procdir,'{:s}*'.format(fieldobsid),'log*.txt') )
    success = []
    if len(outdirs) == len(finished) and len(finished) > 0:
        for outdir in outdirs:
            with open(os.path.join(outdir,'finished.txt'),'r') as f:
                lines = f.readlines()
            ## need to update based on flocs output
            if 'SUCCESS: Pipeline finished successfully' in lines[0]:
                success.append(1)
            else:
                success.append(0)
        if sum(success) == len(outdirs):
            ## everything finished successfully
            success = 'Finished'
            workflow, obsid = get_workflow_obsid(outdirs[0])        
        else:
            idx = np.where(np.asarray(success) == 0)
            print('The following pipelines did not finish successfully, please check processing:')
            print( np.asarray(outdirs)[idx] )    
            success = 'Failed'
            ## will also need to update this
            workflow, obsid = get_workflow_obsid(outdirs[0])
    else:
        ## the number of finished != number of directories - the process is still running
        success = 'Running'
        workflow, obsid = None, None
    return success, workflow, obsid

def cleanup_step(field, fieldobsid):
    basedir = os.getenv('DATA_DIR')
    procdir = os.path.join(str(os.getenv('DATA_DIR')),'processing')
    field_procdirs = glob.glob( os.path.join(procdir,field, fieldobsid+'*') )
    workflow, tmpid = get_workflow_obsid(field_procdirs[0])
    field_datadir = os.path.join(basedir,field,fieldobsid)
    workflowdir = os.path.join(field_datadir,workflow)
    os.makedirs(workflowdir,exist_ok=True)
    ## clean the scratch directory if it's different than procdir
    if os.getenv('SCRATCH_DIR') != procdir:
        os.system('rm -rf {:s}/*'.format(os.path.join(os.getenv('SCRATCH_DIR'),fieldobsid)) )
    for field_procdir in field_procdirs:
        ## remove logs directory (run was successful)
        os.system('rm -rf {:s}'.format(os.path.join(field_procdir,'logs')))
        ## same for tmp directory
        os.system('rm -rf {:s}'.format(os.path.join(field_procdir,'tmp*')))
        ## and workdir
        os.system('rm -rf {:s}'.format(os.path.join(field_procdir,'workdir')))
        ## move everything else to the data directory and rename MSs
        remaining_files = glob.glob(os.path.join(field_procdir,'*'))
        for ff in remaining_files:           
            dest = os.path.join(workflowdir,os.path.basename(ff).replace('out_',''))
            os.system('mv {:s} {:s}'.format(ff,dest))
        ## remove data from previous step if required
        if workflow in ['setup','delay-calibration','target_VLBI']:
            os.system('rm -r {:s}'.format(os.path.join(field_datadir, '*.MS')))
        if workflow in ['HBA_target','target_VLBI']:
            os.system('cp {:s} {:s}'.format(os.path.join(workflowdir,'LINC-cal_solutions.h5'),os.path.join(field_datadir,'LINC-target_solutions.h5')))
            ## a results sub-directory in the results directory in 
        if workflow in ['concatenate-flag']:
            os.system('rm -r {:s}'.format(os.path.join(field_datadir,'setup/L*MS')))
        os.system('rmdir {:s}'.format(field_procdir) )

def upload_to_spider( tarfiles, spider_location='disk/surveys/' ):
    macaroon_dir = os.getenv('MACAROON_DIR')
    macaroon = glob.glob(os.path.join(macaroon_dir,'*maca_lofarvlbi.conf'))[0]
    rc = RClone( macaroon, debug=True )
    rc.get_remote()
    spider_location = spider_location.lstrip('/').rstrip('/') + '/'
    success = 0
    for tarfile in tarfiles:
        d = rc.execute_live(['-P', 'copy', tarfile]+[rc.remote + '/' + spider_location ])
        if d['err'] or d['code']!=0:
            success += 1
            print('Rclone failed for file {:s}'.format(tarfile))
    return(success)

def do_verify(field):
    tarfiles = glob.glob(field+'*tgz')
    success = upload_to_spider( tarfiles, spider_location='disk/surveys' )
    if success == 0:
        print('Tidying uploaded directory for',field)
        update_status(field,'Complete')
        ## delete the directory
        os.system( 'rm -r {:s}'.format(os.path.join(procdir,field)))
        ## delete the initial data
        os.system( 'rm -r {:s}'.format(os.path.join(basedir,field)))
        ## delete the tarfile
        os.system( 'rm {:s}.tgz'.format(field))
    else:
        print('Rclone failed for field {:s}'.format(field) )
        update_status(field,'rclone failed')

def archive_lbfield( field, operation='mv' ):
    if operation == 'copy':
        cmd = 'cp -r '
    else:
        cmd = 'mv '
    ## make an output directory to put things in
    fielddir = os.path.join( os.getenv('DATA_DIR'), field )
    resultsdir = os.path.join( fielddir, 'archive' )
    finaldir = os.path.join( fielddir, 'final' )
    if not os.path.exists(resultsdir):
        os.makedirs(resultsdir)
        os.makedirs(finaldir)

    ## catalogue files
    catdir = os.path.join(fielddir, 'catalogue' )
    catfiles = glob.glob( os.path.join( catdir, '*_catalogue.csv' ) )
    catpngs = glob.glob( os.path.join( catdir, '*png' ) )
    fits = glob.glob( os.path.join( catdir, 'fits' ) )
    inspection = glob.glob( os.path.join( catdir, 'inspection' ) )
    allfiles = catfiles + catpngs + fits + inspection 
    for af in allfiles:
        os.system( cmd + af + ' ' + finaldir + '/' )

    obsids = get_local_obsid( field )
    for obsid in obsids:
        obsdir = os.path.join( fielddir, obsid )
        results_obsdir = os.path.join( resultsdir, obsid )
        os.makedirs( results_obsdir )
        ## general solutions
        lincsols = glob.glob( os.path.join( obsdir, 'LINC*h5' ) )
        if os.path.exists( os.path.join( obsdir, 'delay-calibration' ) ):
            deldir = 'delay-calibration'
        else:
            deldir = 'phaseup-concat'
        tmp = glob.glob( os.path.join( obsdir, deldir, 'merged*verified.h5' ) )
        tmp1 = glob.glob( os.path.join( obsdir, deldir, 'verified_solution_plots' ) )
        tmp2 = glob.glob( os.path.join( obsdir, deldir, 'pipelinesols/ILTJ*dp3-concat' ) )
        delcal = tmp + tmp1 + tmp2
        allfiles = lincsols + delcal
        for af in allfiles:
            os.system( cmd + af + ' ' + results_obsdir + '/' )

        ## selfcal solutions
        sources1 = glob.glob( os.path.join( obsdir, 'selfcal', 'ILTJ*' ) )
        sources2 = glob.glob( os.path.join( obsdir, 'selfcal', 'ILTJ*', 'first_selfcal' ) )
        sources = sources1 + sources2
        os.makedirs( os.path.join( results_obsdir, 'selfcal' ) )
        for source in sources:
            src = os.path.basename( source )
            if src == 'first_selfcal':
                src = os.path.join( os.path.basename( os.path.dirname( source ) ), 'first_selfcal')
            src_outdir = os.path.join( results_obsdir, 'selfcal', src )
            os.makedirs( src_outdir )
            mslist = glob.glob( os.path.join( source, 'ILTJ*ms' ) )
            h5parms = glob.glob( os.path.join( source, 'merged*h5' ) )
            plots = glob.glob( os.path.join( source, 'plotlosoto*' ) )
            mfs = glob.glob( os.path.join( source, '*MFS*fits' ) )
            pngs = glob.glob( os.path.join( source, '*png' ) )
            logs = glob.glob( os.path.join( source, '*log' ) )
            srcfiles = mslist + h5parms + plots + mfs + pngs + logs
            for sf in srcfiles:
                os.system( cmd + sf + ' ' + src_outdir + '/' )
    pwd = os.getenv('PWD')
    os.chdir(os.getenv('DATA_DIR'))
    os.system( 'tar cvzf {:s}_archive.tgz {:s}'.format(field, os.path.join(field,'archive')) )
    os.system( 'tar cvzf {:s}_final.tgz {:s}'.format(field, os.path.join(field,'final')) )
    os.chdir(pwd)

    tarfiles = glob.glob(field+'*tgz')
    success = upload_to_spider( tarfiles, spider_location='disk/fields/' )
    if success == 0:
        print('Tidying uploaded directory for',field)
        update_status(field,'Complete')
        ## delete the directory
        os.system( 'rm -r {:s}'.format(os.path.join(procdir,field)))
        ## delete the initial data
        os.system( 'rm -r {:s}'.format(os.path.join(basedir,field)))
    if success == 4:
        ## delete the tarfile
        os.system( 'rm {:s}*.tgz'.format(field))

