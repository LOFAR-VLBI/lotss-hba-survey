import bdsf
from astropy.table import Table, vstack
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import os, glob, re
import argparse

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def fix_imhead( fitsfile ):
    data, header = fits.getdata( fitsfile, header=True )
    if header['CRVAL1'] < 0:
        print( 'RA is negative, updating header.' )
        header['CRVAL1'] = header['CRVAL1'] + 360.
        fits.writeto( fitsfile, data, header, overwrite=True )
        print( 'done.' )
    else:
        print( 'RA is already positive, exiting.' )

def collapse_source( components ):
    
def get_component( src ):
    tmp = Table()
    tmp.add_column( [src.source_id], name='component' )
    tmp.add_column( [src.posn_sky_centroid[0]], name='RA' )
    tmp.add_column( [src.posn_sky_centroidE[0]], name='E_RA' )
    tmp.add_column( [src.posn_sky_centroid[1]], name='DEC' )
    tmp.add_column( [src.posn_sky_centroidE[1]], name='E_DEC' )
    tmp.add_column( [src.total_flux], name='Total_flux' )
    tmp.add_column( [src.total_fluxE], name='E_Total_flux' )
    tmp.add_column( [src.peak_flux_max], name='Peak_flux' )
    tmp.add_column( [src.peak_flux_maxE], name='E_Peak_flux' )
    tmp.add_column( [src.size_sky[0]], name='Maj' )
    tmp.add_column( [src.size_skyE[0]], name='E_Maj' )
    tmp.add_column( [src.size_sky[1]], name='Min' )
    tmp.add_column( [src.size_skyE[1]], name='E_Min' )
    tmp.add_column( [src.size_sky[2]], name='PA' )
    tmp.add_column( [src.size_skyE[2]], name='E_PA' )
    tmp.add_column( [src.deconv_size_sky[0]], name='DC_Maj' )
    tmp.add_column( [src.deconv_size_skyE[0]], name='E_DC_Maj' )
    tmp.add_column( [src.deconv_size_sky[1]], name='DC_Min' )
    tmp.add_column( [src.deconv_size_skyE[1]], name='E_DC_Min' )
    tmp.add_column( [src.deconv_size_sky[2]], name='DC_PA' )
    tmp.add_column( [src.deconv_size_skyE[2]], name='E_DC_PA' )
    tmp.add_column( [src.rms_isl], name='Isl_rms' )
    tmp.add_column( [src.code], name='S_Code' )
    return(tmp)

def get_source_info( source_file, sc_round='000' ):
    ## source name
    source_name = os.path.basename(source_file)
    ## source detection 
    imf = os.path.join( source_file, 'image_{:s}-MFS-image-pb.fits'.format(sc_round) ) ## true flux
    appf = os.path.join( source_file, 'image_{:s}-MFS-image.fits'.format(sc_round) ) ## apparent flux
    ## fix image header 
    fix_imhead(appf)
    ## run pybdsf - use pb-corrected if it exists
    if os.path.exists(imf):
        fix_imhead(imf)
        img = bdsf.process_image(imf, detection_image=appf, thresh_isl=thresh_isl, thresh_pix=thresh_pix, rms_box=(100,10), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(40,10), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)
    else:
        img = bdsf.process_image(appf, thresh_isl=thresh_isl, thresh_pix=thresh_pix, rms_box=(100,10), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(40,10), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=False, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5, advanced_opts=True, blank_limit=None, frequency=restfrq)
    ## get a table of source information
    sources = img.sources
    source_table = Table()
    for source in sources:
        source_table = vstack([source_table,get_component(source)])
    source_table.add_column(source_name, index=0, name='Source_id' )
    return( source_table ) 

def main( pointing, outdir='catalogue', catfile='pointing_catalogue.fits', update=False ):

    ## bdsf settings
    thresh_pix = 5
    thresh_isl = 4
    restfrq = 144000000.0

    ## pointing input file
    imcat_file = os.path.join(os.getenv('DATA_DIR'),pointing,'image_catalogue.csv')
    imcat = Table.read(imcat_file,format='csv')

    ## the output catalogue
    outcat = os.path.join(os.getenv('DATA_DIR'),pointing,catfile)

    ## output directory
    bdsf_outdir = os.path.join(os.getenv('DATA_DIR'),pointing,outdir)
    if not os.path.exists(bdsf_outdir):
        os.mkdir(bdsf_outdir)
    ## go to this directory
    os.chdir(bdsf_outdir)

    ## get list of sources - second glob is for the case of multi fields
    selfcaldirs = glob.glob(os.path.join(os.getenv('DATA_DIR'),pointing,'*/selfcal')) + glob.glob(os.path.join(os.getenv('DATA_DIR'),pointing,'selfcal'))
    sources = []
    for selfcaldir in selfcaldirs:
        tmp_sources = glob.glob(os.path.join(selfcaldir,'ILTJ*'))
        sources = sources + tmp_sources

    ## If catalogue already exists then open it
    if os.path.exists(outcat) and not update:
        ## danger!
        print('catalogue already exists and update is not specified, stopping.')
        return
    elif os.path.exists(outcat) and update:
        print('catalogue exists and will be updated.')
        pointing_cat = Table.read(outcat,format='fits')
    else:
        print('catalogue does not exist, creating.')
        pointing_cat = Table()

    for source_file in sources:
        ## get lotss info on source
        source = os.path.basename(source_file)
        idx = np.where(imcat['Source_id'] == source)[0]
        lotss_info = imcat[idx]
        majax = lotss_info['Majax']*u.arcsec
        radius = lotss_info['Radius']
        src_coords = SkyCoord( lotss_info['RA'], lotss_info['DEC'], unit='deg' )


        ''' for use when early stopping is in play 
        ## find last selfcal iteration
        last_sc = os.path.basename( natural_sort( glob.glob(os.path.join( source_file, 'image_*-MFS-image.fits')) )[-1] ).split('-')[0].replace('image_','') 
        if int(last_sc) > 0:
            source_info = get_source_info( source_file, sc_round=last_sc )
        else:
            ## there is no source!
            source_info = get_source_info( source_file, sc_round='000' )  
        '''

        ## until then get both first and last selfcal
        source_info_000 = get_source_info( source_file, sc_round='000' )
        for cn in source_info_000.colnames:
            source_info_000.rename_column( cn, '{:s}_0'.format(cn) )
        last_sc = os.path.basename( natural_sort( glob.glob(os.path.join( source_file, 'image_*-MFS-image.fits')) )[-1] ).split('-')[0].replace('image_','') 
        source_info = get_source_info( source_file, sc_round=last_sc )

        if os.path.exists('Default-'+source+'.srl.fits'):
            #If no sources are found, no output file is written
            srl = Table.read('Default-'+source+'.srl.fits', format='fits')
            srl_coords = SkyCoord(srl['RA'], srl['DEC'], unit='deg' )
            seps = srl_coords.separation(src_coords)
            tmp_srl = srl[np.where(seps <= majax)]
            ## add radius from phase centre
            if len(tmp_srl) > 0:
                #Catch situation where no sources lie within the selection
                tmp_srl['Radius'] = radius
                pointing_cat = vstack([pointing_cat,tmp_srl])

    pointing_cat.write(outcat,format='fits')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( dest='pointing', type=str )
    parser.add_argument( '--outdir',type=str,default='catalogue',help='directory to store outputs')
    parser.add_argument( '--catalogue-name',type=str,default='pointing_catalogue.fits',help='name of output catalogue file' )
    parser.add_argument( '--update', action='store_true', default=False
    args = parser.parse_args()
    main( args.pointing, outdir=args.outdir, catfile=args.catalogue_name, update=args.update )

    

