# import bdsf
import os, glob, re
import argparse
import numpy as np
import subprocess
from astropy.table import Table, vstack
from astropy.io import fits
import astropy.units as u

import aplpy
import pylab as plt


################################################
## borrowed from facetselfcal and adjusted to pass back rms

def findrms(mIn, maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m = mIn[np.abs(mIn) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for i in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs((rms - rmsold) / rmsold) < diff: break
        rmsold = rms
    return rms

def get_image_dynamicrange(image):
    """
    Get dynamic range of an image (peak over rms)

    Args:
        image (str): FITS image file name .
     
    Returns:
        DR (float): Dynamic range vale.
    """

    # print('Compute image dynamic range (peak over rms): ', image)
    hdul = fits.open(image)
    image_rms = findrms(np.ndarray.flatten(hdul[0].data))
    DR = np.nanmax(np.ndarray.flatten(hdul[0].data)) / image_rms
    nDR = -np.nanmin(np.ndarray.flatten(hdul[0].data)) / image_rms
    hdul.close()
    return DR, nDR, image_rms

################################################


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def plot_source(source, image, outdir, lotss_dir, dr, ndr, rms):
    font_scale=16
    fits_output_png = os.path.join(outdir, f"{source}.png")

    f = aplpy.FITSFigure(image, hdu = 0)

    #make the image squared
    plt.rcParams['figure.figsize'] = [10,10]
    
    hdul = fits.open(image)
    image_data = hdul[0].data
    
    rms_factor = 2.5
    vmin = -rms_factor*rms
    vmax = np.max((75*rms, np.max(image_data)))
    power_scaling = np.log(0.2)/np.log(2*rms_factor*rms/(vmax-vmin))
    
    f.show_colorscale(cmap='cubehelix_r', stretch='power', exponent=power_scaling, vmin=vmin, vmax=vmax, interpolation = 'none')
    
    n=4
    lotss_fits = os.path.join(lotss_dir, f"{source}_LoTSS.fits")
    _, _, lotss_rms = get_image_dynamicrange(lotss_fits)
    levels_rms = [-lotss_rms*n,lotss_rms*n,lotss_rms*n*2,lotss_rms*n*4,lotss_rms*n*8,lotss_rms*n*16,lotss_rms*n*32,lotss_rms*n*64,lotss_rms*n*128, lotss_rms*n*256, lotss_rms*n*512, lotss_rms*n*1024, lotss_rms*n*2048]
    f.show_contour(lotss_fits, levels=levels_rms, colors='k', linewidths=1, zorder=5)
    
    #colorbar commands
    f.add_colorbar()
    f.colorbar.set_axis_label_font(size=0.9*font_scale, family='sans-serif', variant='small-caps')
    f.colorbar.set_axis_label_text('Jy/beam')
    f.colorbar.set_font(size=0.8*font_scale)
    f.colorbar.set_location('right')
    
    #axis labels and ticks setting
    f.tick_labels.set_font(size=font_scale, family='sans-serif', variant='small-caps')
    f.axis_labels.set_font(size=font_scale, family='sans-serif', variant='small-caps')
    f.axis_labels.set_xtext("Right Ascension (J2000)")
    f.axis_labels.set_ytext("Declination (J2000)")
    f.ticks.set_xspacing(30/3600)
    f.ticks.set_yspacing(25/3600)
    f.ax.coords[1].ticklabels.set_rotation(90)
    
    f.add_label(0.03, 0.96, source, relative=True, horizontalalignment='left', size=0.8*font_scale, weight='bold')
    f.add_label(0.03, 0.92, f"DR: {dr:.1f}", relative=True, horizontalalignment='left', size=0.8*font_scale)
    f.add_label(0.03, 0.88, rf"$\sigma_\mathrm{{rms}}:$ {1e6*rms:.2f} uJy/b", relative=True, horizontalalignment='left', size=0.8*font_scale)
    
    f.save(fits_output_png, dpi=200)
    f.close()

def main(pointing):
    
    pointing = pointing.rstrip('/')

    catfile = f"{pointing}_inspection.csv".format(pointing)

    #Pointing input file
    imcat_file = os.path.join(os.getenv('DATA_DIR'),pointing,'image_catalogue.csv')
    imcat = Table.read(imcat_file,format='csv')

    #The output catalogue
    outcat = os.path.join(os.getenv('DATA_DIR'),pointing,'inspection',catfile)
    
    #Run LoTSS download
    lotss_dir = os.path.join(os.getenv('DATA_DIR'),pointing,'LoTSS')
    lotss_dl_py = os.path.join(os.getenv('SOFTWAREDIR'),'lotss-hba-survey','lotss_dl_lhr.py')
    subprocess.call(f"python3 {lotss_dl_py} --input {imcat_file} --outdir {lotss_dir}", shell=True)

    #Get list of sources - second glob is for the case of multi fields
    selfcaldirs = glob.glob(os.path.join(os.getenv('DATA_DIR'),pointing,'*/selfcal')) + glob.glob(os.path.join(os.getenv('DATA_DIR'),pointing,'selfcal'))
    sources = []
    for selfcaldir in selfcaldirs:
        tmp_sources = glob.glob(os.path.join(selfcaldir,'ILTJ*'))
        sources = sources + tmp_sources
    
    inspection_cat = Table()

    for source_file in sources:
        ## get lotss info on source
        source = os.path.basename(source_file)

        runs = ['']
        alternate_dir = 'first_selfcal'
        if os.path.exists(os.path.join(source_file, alternate_dir)):
            runs.append(alternate_dir)

        drs = []
        ndrs = []
        rmss = []
        for run in runs:
            run_drs = []
            run_ndrs = []
            run_rmss = []
            if os.path.exists( os.path.join( source_file, run, 'image_000-MFS-image-pb.fits' ) ):
                pb = True
                imfiles = natural_sort( glob.glob( os.path.join( source_file, run, 'image_*-MFS-image-pb.fits' ) ) )
            else:
                pb = False
                imfiles = natural_sort( glob.glob( os.path.join( source_file, run, 'image_*-MFS-image.fits' ) ) )
            for imfile in imfiles:
                dr, ndr, rms = get_image_dynamicrange( imfile )
                run_drs.append(dr)
                run_ndrs.append(ndr)
                run_rmss.append(rms)

            drs.append(run_drs)
            ndrs.append(run_ndrs)
            rmss.append(run_rmss)
        
        try:
            regular = True
            drs = np.array(drs)
            ndrs = np.array(ndrs)
            rmss = np.array(rmss)
        except ValueError:
            regular = False

        DR_minimum = 20
        
        if len(runs) > 1 and regular: #Target has multiple runs
            #Case: Barely any bright emission
            if drs[1,0] < DR_minimum: #(check first image of original imaging run)
                #Select first image of re-calibration run
                best_run = runs[0]
                iter_idx = (0,0)
                iteration = '000'
                dr=drs[iter_idx]
                ndr=ndrs[iter_idx]
                rms=rmss[iter_idx]
            else:
                check_iter = drs / np.max(drs) + np.min(rmss) / rmss + np.min(ndrs) / ndrs
                iter_idx = np.unravel_index(np.argmax(check_iter, axis=None), check_iter.shape)
                
                best_run = runs[iter_idx[0]]
                iteration = f'{iter_idx[1]:03}'
                dr=drs[iter_idx]
                ndr=ndrs[iter_idx]
                rms=rmss[iter_idx]
        elif regular: #Target was selected as dd-calibrator
            check_iter = drs / np.max(drs) + np.min(rmss) / rmss + np.min(ndrs) / ndrs
            iter_idx = np.unravel_index(np.argmax(check_iter, axis=None), check_iter.shape)
                
            best_run = runs[iter_idx[0]]
            iteration = f'{iter_idx[1]:03}'
            
            dr=drs[iter_idx]
            ndr=ndrs[iter_idx]
            rms=rmss[iter_idx]
        elif len(runs) > 1:
            if len(drs[0]) == 0:
                good_idx = 1
            else:
                good_idx = 0
            drs = drs[good_idx]
            ndrs = ndrs[good_idx]
            rmss = rmss[good_idx]
            
            check_iter = drs / np.max(drs) + np.min(rmss) / rmss + np.min(ndrs) / ndrs
            iter_idx = np.unravel_index(np.argmax(check_iter, axis=None), check_iter.shape)
            
            best_run = runs[good_idx]
            iteration = f'{iter_idx[0]:03}'
            
            dr=drs[iter_idx[0]]
            ndr=ndrs[iter_idx[0]]
            rms=rmss[iter_idx[0]]
        else:
            best_run = runs[0]
            iteration = '000'
            dr=drs[0][0]
            ndr=ndrs[0][0]
            rms=rmss[0][0]
        
        if pb:
            print(f"Best image for {source_file}: {best_run}/image_{iteration}-MFS-image-pb.fits")
            best_image = os.path.join(source_file, best_run, f"image_{iteration}-MFS-image-pb.fits")
        else:
            print(f"Best image for {source_file}: {best_run}/image_{iteration}-MFS-image.fits")
            best_image = os.path.join(source_file, best_run, f"image_{iteration}-MFS-image.fits")
            
        #Make cutout inspection image
        inspection_dir = os.path.join(os.getenv('DATA_DIR'),pointing,'inspection','images')
        try:
            os.makedirs(inspection_dir)
        except FileExistsError:
            pass
        
        plot_source(source, best_image, inspection_dir, lotss_dir, dr=dr, ndr=ndr, rms=rms)

        tmp_cat = Table()
        tmp_cat.add_column( [source], name='Source_id' )
        # tmp_cat.add_column( [dr], name='DR' )
        # tmp_cat.add_column( [ndr], name='NDR' )
        # tmp_cat.add_column( [rms], name='rms' )
        tmp_cat.add_column( [best_image], name='image' )
        tmp_cat.add_column( [False], name='flagged' )
        inspection_cat = vstack([inspection_cat,tmp_cat])
    inspection_cat.sort('Source_id')
    inspection_cat.write(outcat,format='csv',overwrite=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( dest='pointing', type=str )
    args = parser.parse_args()
    main(args.pointing)
