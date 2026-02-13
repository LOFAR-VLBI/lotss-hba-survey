#!/usr/bin/python

# clone of monitor.py to do LB target processing
from __future__ import print_function
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
from tasklist import *
from calibrator_utils import *
from plot_field import *
import numpy as np
from lbfields_utils import *


#################################
## CLUSTER SPECIFICS - use environment variables

'''
PLEASE SEE THE slurm/add_these_to_bashrc.txt file
'''

user = os.getenv('USER')
if len(user) > 20:
    user = user[0:20]
cluster = os.getenv('DDF_PIPELINE_CLUSTER')
if cluster is None:
    raise RuntimeError('DDF_PIPELINE_CLUSTER must be set')
softwaredir = os.getenv('SOFTWAREDIR')
basedir = os.getenv('DATA_DIR')
if basedir is None:
    raise RuntimeError('DATA_DIR must point to your data directory')
if not os.getenv('LOFAR_SINGULARITY'):
    raise RuntimeError('LOFAR_SINGULARITY must point to a singularity image')
procdir = os.path.join(basedir,'processing')

solutions_thread=None
solutions_name=None
download_thread=None
download_name=None
stage_thread=None
stage_name=None
unpack_thread=None
unpack_name=None
verify_thread=None
verify_name=None

totallimit=20
staginglimit=0 #2
maxstaged=6

## cluster specific queuing limits
if cluster == 'spider':
    maxqueue = 10
elif cluster == 'cosma':
    maxqueue = 3
else:
    #default
    maxqueue = 10
    
'''
updated in MySQL_utils.py:
update_status
get_lbcalibrator
set_lbcalibrator

Logic is as follows:

1. if there is a not started dataset, first operation is always to stage a dataset (NB a different operation if it's on rclone or on SDR -- can do SDR first). At most one staging thread. Set status to Staged on sucessful complete. From this point on we only look for datasets that are associated with the local cluster.

2. any Staged dataset can be downloaded. At most one download thread: set status to Downloaded on successful complete

3. any Downloaded dataset can be unpacked. Set status to Unpacked when done.

4. any Unpacked dataset can have processing run on it. Set status to Started on start. (status set to Verified on upload)

5. upload data products back
'''

while True:

    ## collect information on what's started 
    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select * from lb_fields where clustername="'+cluster+'" and username="'+user+'" order by priority,id')
        result=sdb.cur.fetchall()
        sdb.cur.execute('select * from lb_fields where status="Not started" and priority>0 order by priority,id')
        result2=sdb.cur.fetchall()
        if len(result2)>0:
            nextfield=result2[0]['id']
        else:
            nextfield=None

    d={} ## len(fd)
    fd={}  ## dictionary of fields of a given status type
    for r in result:
        status=r['status']
        if status in d:
            d[status]=d[status]+1
            fd[status].append(r['id'])
        else:
            d[status]=1
            fd[status]=[r['id']]
    d['Not started']=len(result2)
    print('\n\n-----------------------------------------------\n\n')
    print('LB target reprocessing status on cluster %s' % (cluster))
    print(datetime.datetime.now())
    print()
    failed=0
    for k in sorted(d.keys()):
        print("%-20s : %i" % (k,d[k]))
        if 'ailed' in k:
            failed+=d[k]

    ## print out the information
    print()
    ksum=(len(glob.glob(basedir+'/*'))-4)-failed
    if ksum<0: ksum=0
    print(ksum,'live directories out of',totallimit)
    print('Next field to work on is',nextfield)

    ## check what's already running
    if solutions_thread is not None:
        print('Solutions thread is running (%s)' % solutions_name)
    if download_thread is not None:
        print('Download thread is running (%s)' % download_name)
    if unpack_thread is not None:
        print('Unpack thread is running (%s)' % unpack_name)
    if stage_thread is not None:
        print('Stage thread is running (%s)' % stage_name)
    if verify_thread is not None:
        print('Verify thread is running (%s)' % verify_name)

    ## check if anything has terminated
    if solutions_thread is not None and not solutions_thread.is_alive():
        print('Solutions thread seems to have terminated')
        solutions_thread=None
    if download_thread is not None and not download_thread.is_alive():
        print('Download thread seems to have terminated')
        download_thread=None
    if unpack_thread is not None and not unpack_thread.is_alive():
        print('Unpack thread seems to have terminated')
        unpack_thread=None
    if stage_thread is not None and not stage_thread.is_alive():
        print('Stage thread seems to have terminated')
        stage_thread=None
    if verify_thread is not None and not verify_thread.is_alive():
        print('Verify thread seems to have terminated')
        verify_thread=None

    ## Check if staging should be started
    check_stage = 0
    if 'Staging' in d.keys():
        check_stage += d['Staging']
    if 'Solutions' in d.keys():
        check_stage += d['Solutions']
    if 'Staged' in d.keys():
        check_stage += d['Staged']
    if check_stage  <=  staginglimit:
        do_stage = True
    else:
        do_stage = False

    if do_stage and nextfield is not None:
        stage_name=nextfield
        print('We need to stage a new field (%s)' % stage_name)
        solutions_names = stage_field(stage_name)
        ## while staging, collect the solutions
        for solutions_name in solutions_names:
            collect_solutions_lhr(solutions_name)
        ## and download the catalogue
        with SurveysDB(survey=None) as sdb:
            idd=sdb.db_get('lb_fields',stage_name)
        generate_catalogues( float(idd['ra']), float(idd['decl']), targRA = float(idd['ra']), targDEC = float(idd['decl']),
             im_radius=1.24, bright_limit_Jy=5., lotss_result_file='image_catalogue.csv', delay_cals_file='delay_calibrators.csv', 
             match_tolerance=5., image_limit_Jy=0.01, vlass=True, html=False, outdir=os.path.dirname(solutions_names[0]) )

    ## for things that are staging, check the status 
    if 'Staging' in d.keys():
        ## get the staging ids and then check if they are complete
        for field in fd['Staging']:
            ## get the stage id
            r = [ item for item in result if item['id'] == field ][0]
            s = r['staging_id']
            stage_status = stager_access.get_status(s)
            #    "new", "scheduled", "in progress", "aborted", "failed", "partial success", "success", "on hold" 
            if stage_status == 'success' or stage_status == 'completed':
                print('Staging for {:s} is complete, updating status'.format(str(r['staging_id'])))
                update_status(r['id'],'Staged') ## don't reset the staging id till download happens
            else:
                print('Staging for {:s} is {:s} (staging id {:s})'.format(field,stage_status,str(s)))

    ## for things that are staged, start downloading 
    ## if there is not already a download going (one downnload at time allowed)
    if ksum<totallimit and 'Staged' in d and download_thread is None:
        download_name=fd['Staged'][0]
        print('We need to download a new file (%s)!' % download_name)
        download_thread=threading.Thread(target=do_download, args=(download_name,))
        download_thread.start()

    ## unpack the files if downloaded
    if 'Downloaded' in d and unpack_thread is None:
        unpack_name=fd['Downloaded'][0]
        print('We need to unpack a new file (%s)!' % unpack_name)
        unpack_thread=threading.Thread(target=do_unpack, args=(unpack_name,))
        unpack_thread.start()

    ## get number of things that are queued
    if 'Queued' in d:
        nq = d['Queued']
    else:
        nq = 0

    ## if unpacked, start the processing
    if 'Unpacked' in d:
        for field in fd['Unpacked']:
            if nq > maxqueue:
                print( 'Queue is full, {:s} waiting for submission'.format(field) )
            else: 
                nq = nq + 1
                print('Running a new job',field)
                field_obsids = get_local_obsid(field)
                for obsid in field_obsids:
                    next_task = get_task_list(obsid)[0]
                    print('next task is',next_task)
                    fieldobsid = '{:s}/{:s}'.format(field,obsid)
                    ## multiprocessing to run the task
                    run_task( fieldobsid, next_task )
                    ## NEED TO FIX
                    update_status(field,"Submission failed")

    if 'Queued' in d:
        for field in fd['Queued']:
            print('Verifying processing for',field)
            ## will need to fix this
            result, workflow, obsid = check_field(field)
            if result == 'Running':
                pass
            elif result == 'Finished':
                ## get task that was run from the finished.txt to mark done
                mark_done(obsid,workflow.replace('HBA_',''))
                ## and cleanup after the step
                cleanup_step(field,obsid)
                ## start next step
                remaining_tasks = get_task_list(obsid)
                if len(remaining_tasks) > 0:
                    next_task = remaining_tasks[0]
                    fieldobsid = '{:s}/{:s}'.format(field,obsid)
                    run_task( fieldobsid, next_task )
                else:
                    update_status(field,'Verified')  
            else:
                update_status(field,'Workflow failed')

    ## this will also need to be changed to use macaroons to copy back to spider
    if 'Verified' in d and verify_thread is None:
        verify_name = fd['Verified'][0]
        verify_thread=threading.Thread(target=do_verify, args=(verify_name,))
        verify_thread.start()


    print('\n\n-----------------------------------------------\n\n')
    
    sleep(60)
