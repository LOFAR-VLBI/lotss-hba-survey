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
    
with SurveysDB(readonly=True) as sdb:
    #### CHANGE QUERIES TO USE TARGET TABLE
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

print()
ksum=(len(glob.glob(basedir+'/*'))-4)-failed
if ksum<0: ksum=0
print(ksum,'live directories out of',totallimit)
print('Next field to work on is',nextfield)

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

