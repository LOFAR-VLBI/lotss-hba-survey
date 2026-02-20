#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Search LoTSS cutout for a given RA and Dec
Written by Roland Timmerman
"""


#Imports
import os
import sys
import subprocess
import argparse
import glob
import numpy as np
import bdsf


def makeBBSmodel(filename, outfile):
    img = bdsf.process_image(filename, mean_map='zero', adaptive_rms_box=True, rms_map=True, rms_box = (100,10), adaptive_thresh=50, rms_box_bright=(40,10), frequency=144e6)
    img.write_catalog(format='bbs', bbs_patches='single', outfile=outfile, clobber=True)
    return outfile 

def dl_lotss(source_ids, ras, decs, scale=0.1, overwrite=False):
    
    cwd = os.getcwd()
    
    fits_files = []
    
    for idx, source_id in enumerate(source_ids):
        
        #Gather source info
        ra = ras[idx]
        dec = decs[idx]
        
        #Define output filename for this cutout
        out_filename = f"{source_id}_LoTSS.fits"

        #Obtain list of all files present in output directory before download
        prior_files = glob.glob(os.path.join(cwd,'*'))
        
        #Check if output file already exists, if so, skip this file
        if not os.path.exists(os.path.join(cwd,out_filename)) or overwrite:
            
            #Download cutout
            subprocess.call(f"everystamp download --ra {ra} --dec {dec} --survey lotss --lotss_release dr2 --mode fits --size {scale}", shell=True)
            post_files = glob.glob(os.path.join(cwd,'*'))

            #Rename cutout to desired filename
            new_filename = os.path.basename(list(set(post_files)-set(prior_files))[0])
            subprocess.call(f"mv {new_filename} {out_filename}", shell=True)
        
        fits_files.append(out_filename)
    
    return fits_files
        
        
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        '''
        Searches LoTSS for all sources in the LoTSS-HR catalogue per field
        User needs to supply image_catalogue.csv-like file containing Source ID, RA (in decimal degrees), Dec (in decimal degrees), etc...
        ''', formatter_class=argparse.RawTextHelpFormatter)
    
    #Required parameters
    parser.add_argument('--input', help="(Required) image_catalogue.csv-like file", type=str)

    #Optional parameters
    parser.add_argument('--outdir', help="Output directory for all files", default=".", type=str)
    parser.add_argument('--scale', help="Maximum image size in degrees", default=0.1, type=float)
    parser.add_argument('--overwrite', help="Overwrite existing fits files", default=False, type=bool)
    parser.add_argument('--run_bdsf', help="If true: run bdsf on postage stamps", default=False, type=bool)

    options = parser.parse_args()
    args = vars(options)
        
    catalogue = np.loadtxt(args['input'], dtype=str, skiprows=1, usecols=(0,1,2), delimiter=',')
    source_ids = catalogue[:,0].T
    ras = catalogue[:,1].T.astype(float)
    decs = catalogue[:,2].T.astype(float)
    
    subprocess.call(f"mkdir -p {args['outdir']}", shell=True)
    os.chdir(args['outdir'])
    
    lotss_fits = dl_lotss(source_ids, ras, decs, scale=args['scale'], overwrite=args['overwrite']) 
    
    if args['run_bdsf']:
        for fits_file in lotss_fits:
            base_name = fits_file.removesuffix(".fits")
            skymodel = f"{base_name}_BDSF.skymodel"
            makeBBSmodel(fits_file, skymodel)
