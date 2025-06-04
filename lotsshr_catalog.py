import bdsf
from astropy.table import Table, vstack, join
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import os, glob, re
import argparse

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

    print('Compute image dynamic range (peak over rms): ', image)
    hdul = fits.open(image)
    image_rms = findrms(np.ndarray.flatten(hdul[0].data))
    DR = np.nanmax(np.ndarray.flatten(hdul[0].data)) / image_rms
    hdul.close()
    return DR, image_rms

################################################


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
    ## sum total fluxes
    total_flux = np.sum( components['Total_flux'] )
    e_total_flux = np.sqrt( np.sum( np.power( components['E_Total_flux'], 2. ) ) )
    ## directly use peak flux
    peak_flux = np.max( components['Peak_flux'] )
    e_peak_flux = components['E_Peak_flux'][np.where( components['Peak_flux'] == peak_flux )[0]]
    ## flux-weighted average of RA and DEC
    ra_avg = np.average( components['RA'], weights=components['Total_flux'] )
    dec_avg = np.average( components['DEC'], weights=components['Total_flux'] )
    ## errors on the RA and DEC
    e_ra_avg = np.sqrt( np.sum( np.power( components['E_RA'], 2. ) ) )
    e_dec_avg = np.sqrt( np.sum( np.power( components['E_DEC'], 2. ) ) )
    ## rms
    isl_rms = np.mean( components['Isl_rms'] )
    if len(components) > 1:
        ncomp = len(components)
        scode = 'M'
    else:
        ncomp = 1
        scode = components['S_Code'][0]

    ## size and PA ?
    #PA_avg = np.arctan2( np.sum( np.sin( components['PA']*np.pi/180.), np.sum( np.cos( components['PA']*np.pi/180. ) ) ) * 180. / np.pi

    tmp = Table()
    tmp.add_column( [ra_avg], name='RA' )
    tmp.add_column( [e_ra_avg], name='E_RA' )
    tmp.add_column( [dec_avg], name='DEC' )
    tmp.add_column( [e_dec_avg], name='E_DEC' )
    tmp.add_column( [total_flux], name='Total_flux' )
    tmp.add_column( [e_total_flux], name='E_Total_flux' )
    tmp.add_column( [peak_flux], name='Peak_flux' )
    tmp.add_column( [e_peak_flux], name='E_Peak_flux' )
    tmp.add_column( [isl_rms], name='Isl_rms' )
    tmp.add_column( [ncomp], name='Ncomponents' )
    tmp.add_column( [scode], name='S_Code' )
    return( tmp )

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

    ## bdsf settings
    thresh_pix = 5
    thresh_isl = 5

    ## source name
    source_name = os.path.basename(source_file)
    ## source detection images
    imf = os.path.join( source_file, 'image_{:s}-MFS-image-pb.fits'.format(sc_round) ) ## true flux
    appf = os.path.join( source_file, 'image_{:s}-MFS-image.fits'.format(sc_round) ) ## apparent flux
    ## fix image header 
    fix_imhead(appf)
    ## run pybdsf - use pb-corrected if it exists
    if os.path.exists(imf):
        fix_imhead(imf)
        bdsf_im = imf
    else:
        bdsf_im = appf

    print(bdsf_im)

    ## get the noise
    with fits.open( bdsf_im ) as hdul:
        imnoise = findrms(np.ndarray.flatten(hdul[0].data))
        try:
            reffreq = hdul[0].header['CRVAL3']
        except:
            reffreq = hdul[0].header['REF_FREQ']

    img = bdsf.process_image(bdsf_im, 
                        adaptive_rms_box=True, 
                        adaptive_thresh=100, 
                        rms_box=(150,10), 
                        rms_box_bright=(40,10),
                        rms_map=False, 
                        rms_value = imnoise,
                        mean_map='zero', 
                        thresh_isl=thresh_isl, 
                        thresh_pix=thresh_pix, 
                        group_by_isl=False,
                        group_tol=10.0, 
                        flagging_opts=True, 
                        flag_maxsize_fwhm=0.5,
                        atrous_do=True, 
                        atrous_jmax=4,
                        frequency=reffreq)

    ## get a table of source information
    sources = img.sources
    source_table = Table()
    for source in sources:
        source_table = vstack([source_table,get_component(source)])
    source_table.add_column(source_name, index=0, name='Source_id' )
    return( source_table ) 

def main( pointing, outdir='catalogue', catfile='catalogue.fits', update=False ):

    catfile = '{:s}_{:s}'.format(pointing,catfile)

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
        combined_cat = Table()

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

        drs = []
        rmss = []
        if os.path.exists( os.path.join( source_file, 'image_000-MFS-image-pb.fits' ) ):
            imfiles = natural_sort( glob.glob( os.path.join( source_file, 'image_*-MFS-image-pb.fits' ) ) )
        else:
            imfiles = natural_sort( glob.glob( os.path.join( source_file, 'image_*-MFS-image.fits' ) ) )
        for imfile in imfiles:
            dr, rms = get_image_dynamicrange( imfile )
            drs.append(dr)
            rmss.append(rms)

        ## decide which iteration to use with normalised combination
        check_iter = drs / np.max(drs) + rmss / np.max(rmss)
        iter_idx = np.where( check_iter == np.max(check_iter) )[0][0]
        iteration = f'{iter_idx:03}'
        
        ## get source information - add to pointing catalogue
        source_info = get_source_info( source_file, sc_round=iteration )
        pointing_cat = vstack([pointing_cat,source_info])

        ## collapse the information - add to image catalogue
        source_coll = collapse_source( source_info )
        source_coll.add_column( [source], name='Source_id', index=0 )
        source_coll.add_column( [pointing], name='Pointing' )
        source_coll.add_column( [radius], name='Radius' )
        source_coll.add_column( [iteration], name='SelfcalIteration' )
        combined_cat = vstack( [combined_cat, source_coll] )


    ## for testing
    t_dr = []
    t_rms = []
    for source_file in sources:
        drs = []
        rmss = []
        if os.path.exists( os.path.join( source_file, 'image_000-MFS-image-pb.fits' ) ):
            imfiles = natural_sort( glob.glob( os.path.join( source_file, 'image_*-MFS-image-pb.fits' ) ) )
        else:
            imfiles = natural_sort( glob.glob( os.path.join( source_file, 'image_*-MFS-image.fits' ) ) )
        for imfile in imfiles:
            dr, rms = get_image_dynamicrange( imfile )
            drs.append(dr)
            rmss.append(rms)
        t_dr.append(drs)
        t_rms.append(rmss)

    with open( 'iteration_info.csv', 'w' ) as f:
        f.write('source,dr0,dr1,dr2,dr3,dr4,dr5,dr6,dr7,dr8,dr9,rms0,rms1,rms2,rms3,rms4,rms5,rms6,rms7,rms8,rms9\n')
        for i in np.arange(0,len(sources)):
            src = os.path.basename(sources[i])
            a = [ str(val) for val in t_dr[i] ]
            b = [ str(val) for val in t_rms[i] ]
            f.write('{:s},{:s},{:s}\n'.format(src,','.join(a),','.join(b)) )


    
    ## join combined with lotss
    for cn in combined_cat.colnames:
        if cn in imcat.colnames:
            combined_cat.rename_column( cn, '{:s}_HR'.format(cn) )
    combined_cat.rename_column('Source_id_HR','Source_id')
    
    lotss_combined = join(left=imcat, right=combined_cat, keys='Source_id')
    
    ## write out catalogues
    pointing_cat.write(outcat,format='fits')
    lotss_combined.write(outcat.replace('{:s}_'.format(pointing),'combined_'),format='fits')



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( dest='pointing', type=str )
    parser.add_argument( '--outdir',type=str,default='catalogue',help='directory to store outputs')
    parser.add_argument( '--catalogue-filestem',type=str,default='catalogue.fits',help='suffix of catalogue file' )
    parser.add_argument( '--update', action='store_true', default=False )
    args = parser.parse_args()
    main( args.pointing, outdir=args.outdir, catfile=args.catalogue_filestem, update=args.update )

    

