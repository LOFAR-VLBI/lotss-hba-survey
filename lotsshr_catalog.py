import os, glob
import numpy as np
import bdsf
import argparse
from astropy.table import Table, vstack, join
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats

import aplpy
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath


################################################
## borrowed from facetselfcal and adjusted to pass back rms

def sign(x):
    if x<0:
        return -1
    else:
        return 1

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

def mask_outliers(array):
    array = np.array(array)
    mask = np.ones(len(array), dtype=bool)
    for idx in range(4):
        temp_array = array[mask]
        med = np.median(temp_array)
        std = np.std(temp_array)
        mask = np.array(array < med+3*std) & np.array(array > med-3*std)
    return mask
        
def wcs_2d_from_header(header):
    wcs = WCS(header)
    if header['NAXIS'] == 4:
        wcs = wcs.dropaxis(2).dropaxis(2)
    return wcs

def identify_contour(x, y, segs, lotss_idx, grab_nearest=False):
    starting_level = 1
    for layer in range(starting_level, len(segs)):
        seg_layer = segs[layer]
        for idx, polygon in enumerate(seg_layer):
            path = mpltPath.Path(polygon)
            inside = path.contains_points(np.array([x,y]).T)
            if inside[lotss_idx] == True and np.sum(inside.astype(bool)) == 1:
                return layer, idx
    if grab_nearest:
        seg_layer = segs[starting_level]
        min_sep = 1e6
        min_idx = 0
        for idx, polygon in enumerate(seg_layer):
            x_poly, y_poly = polygon.T
            sep = np.sqrt((x_poly-x[lotss_idx])**2+(y_poly-y[lotss_idx])**2)
            if np.min(sep) < min_sep:
                min_sep = np.min(sep)
                min_idx = idx
        return starting_level, min_idx
    return None, None

def fix_imhead(fitsfile):
    data, header = fits.getdata(fitsfile, header=True)
    if header['CRVAL1'] < 0:
        header['CRVAL1'] = header['CRVAL1'] + 360.
        fits.writeto(fitsfile, data, header, overwrite=True)
            
def reorder_cat(srcs, comps):
    src_map = {}
    for idx, src_id in enumerate(srcs['source_id']):
        src_map[src_id] = idx
    
    Nsrcs = len(srcs)
    Ncomps = len(comps)
    for idx in range(Ncomps):
        comps['source_id'][idx] = src_map[comps['source_id'][idx]]
    srcs['source_id'] = range(Nsrcs)
    comps['gaus_num'] = range(Ncomps)
    return srcs, comps

def collapse_source( components ):
    ## sum total fluxes
    total_flux = np.sum(components['Total_flux'])
    if total_flux > 0:
        # Calculate error on total flux density
        e_total_flux = np.sqrt(np.sum(components['E_Total_flux']**2))
        
        # Grab peak flux -RT: I don't think this would actually generally correspond to the brightest pixel if components overlap, but okay for now.
        peak_flux_idx = np.argmax(components['Peak_flux'])
        peak_flux = components['Peak_flux'][peak_flux_idx]
        e_peak_flux = components['E_Peak_flux'][peak_flux_idx]
        
        # Flux-weighted average of RA and DEC
        ra_avg = np.average(components['RA'], weights=components['Total_flux'])
        dec_avg = np.average(components['DEC'], weights=components['Total_flux'])
        
        # Errors on the RA and DEC
        e_ra_avg = np.sqrt(np.sum(components['E_RA']**2))
        e_dec_avg = np.sqrt(np.sum(components['E_DEC']**2))
        
        # rms
        isl_rms = np.mean(components['Isl_rms'])
        
        ncomp = len(components)
        if ncomp > 1:
            scode = 'M'
        else:
            scode = components['S_Code'][0]
    else:
        ra_avg = components['RA'][0]
        e_ra_avg = 0
        dec_avg = components['DEC'][0]
        e_dec_avg = 0
        e_total_flux = 0
        peak_flux = 0
        e_peak_flux = 0
        isl_rms = 0
        ncomp = 0
        scode = 'X'

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

def get_source(src):
    tmp = Table()
    tmp.add_column( [src.source_id], name='source_id' )
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

def get_component(comp):
    tmp = Table()
    tmp.add_column( [comp.source_id], name='source_id' )
    tmp.add_column( [comp.gaus_num], name='gaus_num' )
    tmp.add_column( [comp.centre_sky[0]], name='RA' )
    tmp.add_column( [comp.centre_skyE[0]], name='E_RA' )
    tmp.add_column( [comp.centre_sky[1]], name='DEC' )
    tmp.add_column( [comp.centre_skyE[1]], name='E_DEC' )
    tmp.add_column( [comp.total_flux], name='Total_flux' )
    tmp.add_column( [comp.total_fluxE], name='E_Total_flux' )
    tmp.add_column( [comp.peak_flux], name='Peak_flux' )
    tmp.add_column( [comp.peak_fluxE], name='E_Peak_flux' )
    tmp.add_column( [comp.size_sky[0]], name='Maj' )
    tmp.add_column( [comp.size_skyE[0]], name='E_Maj' )
    tmp.add_column( [comp.size_sky[1]], name='Min' )
    tmp.add_column( [comp.size_skyE[1]], name='E_Min' )
    tmp.add_column( [comp.size_sky[2]], name='PA' )
    tmp.add_column( [comp.size_skyE[2]], name='E_PA' )
    tmp.add_column( [comp.deconv_size_sky[0]], name='DC_Maj' )
    tmp.add_column( [comp.deconv_size_skyE[0]], name='E_DC_Maj' )
    tmp.add_column( [comp.deconv_size_sky[1]], name='DC_Min' )
    tmp.add_column( [comp.deconv_size_skyE[1]], name='E_DC_Min' )
    tmp.add_column( [comp.deconv_size_sky[2]], name='DC_PA' )
    tmp.add_column( [comp.deconv_size_skyE[2]], name='E_DC_PA' )
    tmp.add_column( [comp.rms], name='rms' )
    tmp.add_column( [comp.code], name='S_Code' )
    return(tmp)

def add_zero_source(tmp, ra, dec, component=False):
    if component:
        tmp.add_column([0], name='gaus_num')
    tmp.add_column( [0],   name='source_id')
    tmp.add_column( [ra],  name='RA' )
    tmp.add_column( [0.0], name='E_RA' )
    tmp.add_column( [dec], name='DEC' )
    tmp.add_column( [0.0], name='E_DEC' )
    tmp.add_column( [0.0], name='Total_flux' )
    tmp.add_column( [0.0], name='E_Total_flux' )
    tmp.add_column( [0.0], name='Peak_flux' )
    tmp.add_column( [0.0], name='E_Peak_flux' )
    tmp.add_column( [0.0], name='Maj' )
    tmp.add_column( [0.0], name='E_Maj' )
    tmp.add_column( [0.0], name='Min' )
    tmp.add_column( [0.0], name='E_Min' )
    tmp.add_column( [0.0], name='PA' )
    tmp.add_column( [0.0], name='E_PA' )
    tmp.add_column( [0.0], name='DC_Maj' )
    tmp.add_column( [0.0], name='E_DC_Maj' )
    tmp.add_column( [0.0], name='DC_Min' )
    tmp.add_column( [0.0], name='E_DC_Min' )
    tmp.add_column( [0.0], name='DC_PA' )
    tmp.add_column( [0.0], name='E_DC_PA' )
    if component:
        tmp.add_column( [0.0], name='rms' )
    else:
        tmp.add_column( [0.0], name='Isl_rms' )
    tmp.add_column( ['X'], name='S_Code' )
    return(tmp)

def gaia_quasar_match(ra, dec, max_off):
    from astroquery.vizier import Vizier

    t = Table()
    t.add_column( [ra * u.deg], name="_RAJ2000")
    t.add_column( [dec * u.deg], name="_DEJ2000")
    CatNorth = "J/ApJS/271/54/table4"

    v = Vizier(catalog=CatNorth, columns=["Gaia", "RA_ICRS", "DE_ICRS"])

    # Match against table
    # print(f"Querying GAIA quasars for {ra}, {dec}")
    out = v.query_region(t, radius=f"{max_off}s", cache=False)
    if out:
        return out[0]
    else:    
        return None

def match_lotss(ra, dec, max_off):
    from astroquery.vizier import Vizier

    t = Table()
    t.add_column( [ra * u.deg], name="_RAJ2000")
    t.add_column( [dec * u.deg], name="_DEJ2000")
    LoTSSDR2 = "J/A+A/659/A1/catalog"

    v = Vizier(catalog=LoTSSDR2, columns=["Source", "RAJ2000", "DEJ2000", "Speak", "e_Speak", "SpeakTot", "e_SpeakTot", "Islrms", "SCode"])

    # Match against table
    # print(f"Querying GAIA quasars for {ra}, {dec}")
    out = v.query_region(t, radius=f"{max_off}s", cache=False)
    if out:
        return out[0]
    else:    
        return None

def get_source_info(source_file, source_name, lotss_info):

    ## bdsf settings
    thresh_pix = 5
    thresh_isl = 5

    ## fix image header 
    fix_imhead(source_file)

    ## get the noise
    with fits.open(source_file) as hdul:
        imnoise = findrms(np.ndarray.flatten(hdul[0].data))
        try:
            reffreq = hdul[0].header['CRVAL3']
        except:
            reffreq = hdul[0].header['REF_FREQ']

    img = bdsf.process_image(source_file, 
                        adaptive_rms_box=True, 
                        adaptive_thresh=100, 
                        rms_box=(100,10), 
                        rms_box_bright=(40,10),
                        rms_map=True, 
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
                        atrous_orig_isl=True,
                        frequency=reffreq,
                        quiet=True)

    ## get a table of source information
    sources = img.sources
    source_table = Table()
    component_table = Table()
    if len(sources) > 0:
        for source in sources:
            source_table = vstack([source_table,get_source(source)])
        gaussians = img.gaussians
        for gauss in gaussians:
            component_table = vstack([component_table,get_component(gauss)])
    else:
        source_table = add_zero_source(source_table, float(lotss_info['RA'][0]), float(lotss_info['DEC'][0]))
        component_table = add_zero_source(component_table, float(lotss_info['RA'][0]), float(lotss_info['DEC'][0]), component=True)
    # source_table.add_column(source_name, index=0, name='Source_id' )
    return source_table, component_table 

def plot_source(source, image, outdir, lotss_dir, component_info, imcat, lotss_idx, flagged=False):
    font_scale=16
    fits_output_png = os.path.join(outdir, f"{source}.png")

    _, _, rms = get_image_dynamicrange(image)
    
    f = aplpy.FITSFigure(image, hdu = 0)

    #make the image squared
    plt.rcParams['figure.figsize'] = [10,10]
    # plt.rcParams['text.usetex'] = False
    #using Helvetica font family
    # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    
    hdul = fits.open(image)
    image_data = hdul[0].data
    image_header = hdul[0].header
    
    rms_factor = 2.5
    vmin = -rms_factor*rms
    vmax = np.max((75*rms, np.max(image_data)))
    power_scaling = np.log(0.2)/np.log(2*rms_factor*rms/(vmax-vmin))
    
    f.show_colorscale(cmap='cubehelix_r', stretch='power', exponent=power_scaling, vmin=vmin, vmax=vmax, interpolation='none')
    
    n=4
    lotss_fits = os.path.join(lotss_dir, f"{source}_LoTSS.fits")
    _, _, lotss_rms = get_image_dynamicrange(lotss_fits)
    levels_rms = [-lotss_rms*n,lotss_rms*n,lotss_rms*n*2,lotss_rms*n*4,lotss_rms*n*8,lotss_rms*n*16,lotss_rms*n*32,lotss_rms*n*64,lotss_rms*n*128, lotss_rms*n*256, lotss_rms*n*512, lotss_rms*n*1024, lotss_rms*n*2048]
    f.show_contour(lotss_fits, levels=levels_rms, colors='k', linewidths=1, zorder=5)
    
    #Find main source contour and determine flux density
    wcs = wcs_2d_from_header(image_header)
    x, y = wcs.wcs_world2pix(imcat['RA'], imcat['DEC'], 0)
    segs = f.get_layer('contour_set_1').allsegs
    contour_idx, seg_idx = identify_contour(x, y, segs, lotss_idx, grab_nearest=True)
    base_contour = segs[contour_idx][seg_idx]
    
    x_comps, y_comps = wcs.wcs_world2pix(component_info['RA'], component_info['DEC'], 0)
    path = mpltPath.Path(base_contour)
    inside = np.array(path.contains_points(np.array([x_comps, y_comps]).T)).astype(bool)
    fluxdensity = np.sum(component_info['Total_flux'][inside].astype(float))
    
    #Find nearby LoTSS detections
    margin = 1.25
    image_size = image_header['NAXIS1']*np.abs(image_header['CDELT2'])/2*margin
    seren_mask = np.zeros(len(component_info), dtype=bool)
    for idx, lotss_source in enumerate(imcat):
        if idx != lotss_idx and np.abs(lotss_source['RA'] - imcat['RA'][lotss_idx])/np.cos(np.deg2rad(imcat['DEC'][lotss_idx])) < image_size and np.abs(lotss_source['DEC'] - imcat['DEC'][lotss_idx]) < image_size:
            contour_idx, seg_idx = identify_contour(x, y, segs, idx)
            if contour_idx is None:
                continue
            base_contour = segs[contour_idx][seg_idx]
            path = mpltPath.Path(base_contour)
            seren_inside = np.array(path.contains_points(np.array([x_comps, y_comps]).T)).astype(bool)
            seren_mask = seren_mask | seren_inside
    
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
    f.add_label(0.03, 0.92, rf"$I_\mathrm{{144~MHz}}:$ {fluxdensity:.3e} Jy", relative=True, horizontalalignment='left', size=0.8*font_scale)
    f.add_label(0.03, 0.88, rf"$\sigma_\mathrm{{rms}}:$ {1e6*rms:.2f} $\mu$Jy/b", relative=True, horizontalalignment='left', size=0.8*font_scale)
    if flagged:
            f.add_label(0.97, 0.03, "FLAGGED", relative=True, horizontalalignment='right', size=1.2*font_scale, weight='bold', c='r')

    f.save(fits_output_png, dpi=200)
    f.close()
    
    return inside, seren_mask

def main( pointing, outdir='catalogue', processingdir='postprocessing', inspection_file='default', cat_format='csv', update=True ):

    pointing = pointing.rstrip('/')
    
    ## output directories
    postprocessingdir = os.path.join(os.getenv('DATA_DIR'),pointing,processingdir)
    outputdir = os.path.join(os.getenv('DATA_DIR'),pointing,outdir)
    if not os.path.exists(postprocessingdir):
        os.mkdir(postprocessingdir)
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
        
    pointing_catfile = os.path.join(outputdir,f"{pointing}_source_catalogue.{cat_format}")
    gaussian_catfile = os.path.join(outputdir,f"{pointing}_gauss_catalogue.{cat_format}")
    seren_src_catfile = os.path.join(outputdir,f"{pointing}_serendip_src_catalogue.{cat_format}")
    seren_gauss_catfile = os.path.join(outputdir,f"{pointing}_serendip_gauss_catalogue.{cat_format}")

    ## pointing input file
    imcat_file = os.path.join(os.getenv('DATA_DIR'),pointing,'image_catalogue.csv')
    imcat = Table.read(imcat_file,format='csv')
    
    if inspection_file == "default":
        inspection_file = f"{pointing}_inspection.csv"
        inspection_file = os.path.join(os.getenv('DATA_DIR'),pointing,'inspection',inspection_file)
        inspection_cat = Table.read(inspection_file,format='csv')

    ## First work in the postprocessing directory
    os.chdir(postprocessingdir)

    ## get list of sources - second glob is for the case of multi fields
    selfcaldirs = glob.glob(os.path.join(os.getenv('DATA_DIR'),pointing,'*/selfcal')) + glob.glob(os.path.join(os.getenv('DATA_DIR'),pointing,'selfcal'))
    sources = []
    for selfcaldir in selfcaldirs:
        tmp_sources = glob.glob(os.path.join(selfcaldir,'ILTJ*'))
        sources = sources + tmp_sources

    pointing_cat = Table()
    gaussian_cat = Table()
    seren_src_cat = Table()
    seren_gauss_cat = Table()
    
    #Prepare lists
    ra_diff = []
    dec_diff = []
    snr_diff = []
    
    lotss_flux = []
    lotss_flux_e = []
    lotsshr_flux = []
    lotsshr_flux_e = []
    
    #Maximum permittable offset in arcseconds between LoTSS-HR components and GAIA
    max_off = 20
    
    #Minimum significance for astrometry
    min_rms_level = 25

    #Iterate over the sources for astrometry and flux scale
    for source_file in sources:
        #Get LoTSS info on source
        source = os.path.basename(source_file)
        
        #Obtain entry in inspection catalogue
        idx = np.where(inspection_cat['Source_id'] == source)[0]
        inspection_info = inspection_cat[idx]
        
        #Obtain LoTSS info
        lotss_idx = np.where(imcat['Source_id'] == source)[0]
        lotss_info = imcat[lotss_idx]
        
        #Define source and component lists
        source_bdsf = os.path.join(postprocessingdir, f"{source}_src_BDSF.csv")
        comp_bdsf = os.path.join(postprocessingdir, f"{source}_comp_BDSF.csv")
        
        #Try to open existing source and component tables - else generate these
        if os.path.exists(source_bdsf) and os.path.exists(comp_bdsf):
            source_info = Table.read(source_bdsf,format='csv')
            component_info = Table.read(comp_bdsf,format='csv')
        else:
            source_info, component_info = get_source_info(inspection_info['image'][0], source, lotss_info)
            source_info.write(source_bdsf,format='csv',overwrite=True)
            component_info.write(comp_bdsf,format='csv',overwrite=True)
        
        #If source is flagged, don't use it for astrometry or flux scale
        if inspection_info['flagged'][0] != "False":
            continue
        
        #Check if compact
        if len(component_info) == 2:
            #A substantial amount of compact sources get identified as two components in LoTSS-HR. Can probably include these as "compact"
            ra_sep = (((component_info['RA'][0] - component_info['RA'][1]) * u.deg).to(u.arcsec).value)*np.cos(np.deg2rad(component_info['DEC'][0]))
            dec_sep = ((component_info['DEC'][0] - component_info['DEC'][1]) * u.deg).to(u.arcsec).value
            total_sep = np.sqrt(ra_sep**2 + dec_sep**2)
            if total_sep < 1: # and component_info['S_Code'][0] == 'M' and component_info['S_Code'][1] == 'M':
                # print(f"{source} consists of two nearby compact components in LoTSS-HR")
                two_nearby = True
            else:
                two_nearby = False
        else:
            two_nearby = False
        compact_enough = (len(component_info) == 1 and component_info['S_Code'][0] == 'S') or two_nearby
        if compact_enough:
            lotss_match = match_lotss(source_info['RA'][0], source_info['DEC'][0], max_off)
            if lotss_match is not None:
                idx = np.where(lotss_match['Source'] == source)[0]
                if len(idx) > 0:
                    match = lotss_match[idx]
                    lotss_peak_ratio = float(match['Speak'][0])/float(match['SpeakTot'][0])
                    if lotss_peak_ratio > 0.7 or match['SCode'][0] == 'S':
                        #Obtain flux ratio measurement
                        lotss_flux.append(float(match['SpeakTot'][0])/1e3)
                        lotss_flux_e.append(float(match['e_SpeakTot'][0])/1e3)
                        lotss_hr = np.sum(component_info['Total_flux'].astype(float))
                        lotss_hr_e = np.sqrt(np.sum(component_info['E_Total_flux'].astype(float))**2)
                        lotsshr_flux.append(lotss_hr)
                        lotsshr_flux_e.append(lotss_hr_e)
                        
                        #Obtain astrometry measurement
                        lotss_hr_ra = np.average(component_info['RA'], weights=component_info['Total_flux'])
                        lotss_hr_dec = np.average(component_info['DEC'], weights=component_info['Total_flux'])
                        ra_sep = (((lotss_hr_ra-match['RAJ2000'][0]) * u.deg).to(u.arcsec).value)*np.cos(np.deg2rad(match['DEJ2000'][0]))
                        dec_sep = ((lotss_hr_dec-match['DEJ2000'][0]) * u.deg).to(u.arcsec).value
                        snr_sep = lotss_hr/lotss_hr_e
                        ra_diff.append(ra_sep)
                        dec_diff.append(dec_sep)
                        snr_diff.append(snr_sep)
        
    #Move to final output folder
    os.chdir(outputdir)
    
    #Plot astrometry results
    ra_diff = np.array(ra_diff)
    dec_diff = np.array(dec_diff)
    snr_diff = np.array(snr_diff)
    
    if len(ra_diff)>0:
        ra_mask = mask_outliers(ra_diff)
        dec_mask = mask_outliers(dec_diff)
        mask = ra_mask & dec_mask
                
        ra_offset = np.average(ra_diff[mask], weights=snr_diff[mask])
        dec_offset = np.average(dec_diff[mask], weights=snr_diff[mask])
        
        total_offset = np.sqrt((ra_diff[mask]-ra_offset)**2+(dec_diff[mask]-dec_offset)**2)
        _, astro_err = stats.rayleigh.fit(total_offset)
        
        print(f"RA offset: {ra_offset:.5f} arcseconds")
        print(f"Dec offset: {dec_offset:.5f} arcseconds")
        print(f"Astrometric error: {astro_err:.5f} arcseconds")

        fig = plt.figure(figsize=(8,7))
        ax = plt.gca()
        plt.scatter(ra_diff[mask], dec_diff[mask], c='royalblue', alpha=0.6, s=2e2/np.sqrt(snr_diff[mask]),zorder=10)
        plt.scatter(ra_diff[~mask], dec_diff[~mask], c='r', alpha=0.75, s=2e2/np.sqrt(snr_diff[~mask]),zorder=8)
        plt.axvline(ra_offset, c='royalblue', linestyle='dashed', alpha=0.5, linewidth=1)
        plt.axhline(dec_offset, c='royalblue', linestyle='dashed', alpha=0.5, linewidth=1)
        plt.axvline(0, c='k', linestyle='dotted', linewidth=0.5)
        plt.axhline(0, c='k', linestyle='dotted', linewidth=0.5)
        plt.xlabel("Offset in R.A. [arcseconds]")
        plt.ylabel("Offset in Dec. [arcseconds]")
        plt.xlim((max_off, -max_off))
        plt.ylim((-max_off, max_off))
        plt.text(.02, .98, f"Field: {pointing}\nR.A. offset: {ra_offset:.4f}''\nDec. offset: {dec_offset:.4f}''\nScatter: {astro_err:.4f}''\nN={np.sum(mask.astype(int))}", ha='left', va='top', transform=ax.transAxes)
        plt.tight_layout()
        plt.savefig('astrometry.png', dpi=300)
        plt.close(fig)

    #Plot flux scale results
    lotss_flux = np.array(lotss_flux)
    lotss_flux_e = np.array(lotss_flux_e)
    lotsshr_flux = np.array(lotsshr_flux)
    lotsshr_flux_e = np.array(lotsshr_flux_e)
    
    if len(lotss_flux)>0:
        flux_ratio = lotss_flux/lotsshr_flux
        flux_ratio_e = flux_ratio*np.sqrt((lotss_flux_e/lotss_flux)**2+(lotsshr_flux_e/lotsshr_flux)**2)
        
        flux_correction = np.average(flux_ratio, weights=flux_ratio_e**-2)
        print(f"Flux scale factor: {flux_correction:.5f}")
    
        fig = plt.figure(figsize=(8,7))
        ax = plt.gca()
        plt.errorbar(lotsshr_flux, lotss_flux, yerr=lotss_flux_e, xerr=lotsshr_flux_e, ls='none', c='royalblue', fmt='o')
        plt.xscale('log')
        plt.yscale('log')
        xlim = plt.xlim()
        ylim = plt.ylim()
        new_low = min(xlim[0],ylim[0])
        new_high = max(xlim[1],ylim[1])
        plt.plot([1e-6, 1e3], [1e-6, 1e3], c='k', linestyle='dotted', linewidth=0.5, label='1:1')
        plt.plot([1e-6, 1e3], [flux_correction*1e-6, flux_correction*1e3], c='royalblue', alpha=0.5, linestyle='dashed', label='Flux scale fit')
        plt.xlim((new_low, new_high))
        plt.ylim((new_low, new_high))
        plt.xlabel("Flux density in LoTSS-HR [Jy]")
        plt.ylabel("Flux density in LoTSS [Jy]")
        plt.text(.02, .98, f"Field: {pointing}\nFlux scale correction: {flux_correction:.4f}\nN={len(lotss_flux)}", ha='left', va='top', transform=ax.transAxes)
        plt.legend(loc='lower right')
        plt.tight_layout()
        plt.savefig('flux_scale.png', dpi=300)
        plt.close(fig)
        
    #Test values, do not use for real catalogues
    # ra_offset = 0.000
    # dec_offset = 0.000
    # flux_correction = 1.000
    
    #Output folder for corrected fits images
    corrected_fitsfolder = os.path.join(outputdir, 'fits')
    if not os.path.exists(corrected_fitsfolder):
        os.mkdir(corrected_fitsfolder)
        
    final_plotfolder = os.path.join(outputdir, 'inspection')
    if not os.path.exists(final_plotfolder):
        os.mkdir(final_plotfolder)
        
    lotss_dir = os.path.join(os.getenv('DATA_DIR'),pointing,'LoTSS')
        
    #Iterate over the sources for astrometry and flux scale
    for source_file in sources:
        source = os.path.basename(source_file)
        
        lotss_idx = np.where(imcat['Source_id'] == source)[0]
        lotss_info = imcat[lotss_idx]
                        
        idx = np.where(inspection_cat['Source_id'] == source)[0]
        inspection_info = inspection_cat[idx]
        
        source_bdsf = os.path.join(postprocessingdir, f"{source}_src_BDSF.csv")
        comp_bdsf = os.path.join(postprocessingdir, f"{source}_comp_BDSF.csv")
        source_info = Table.read(source_bdsf,format='csv')
        component_info = Table.read(comp_bdsf,format='csv')
        
        flagged = inspection_info['flagged'][0] != "False"
        
        #Apply corrections to fits file
        corrected_fitsfile = os.path.join(corrected_fitsfolder, f"{source}_LoTSS-HR.fits")
        if not os.path.exists(corrected_fitsfile):
                
            #Open up the fits image
            hdul = fits.open(inspection_info['image'][0])
            data = hdul[0].data
            header = hdul[0].header
            
            #Apply corrections - flux correction
            data *= flux_correction
            
            #Apply corrections - astrometry
            if 'RA' in header['CTYPE1'] and 'DEC' in header['CTYPE2'] and header['CUNIT1'] == 'deg' and header['CUNIT2'] == 'deg':
                #Header is as expected
                delta_ra_deg = (ra_offset * u.arcsec).to(u.deg).value
                delta_dec_deg = ((dec_offset * u.arcsec).to(u.deg).value) / np.cos(np.deg2rad(header['CRVAL2']))
                
                header['CRVAL1'] -= delta_ra_deg
                header['CRVAL2'] -= delta_dec_deg
            else:
                raise KeyError('FITS headers are not as expected. Cataloguing code may need an update')
            
            newhdu = fits.PrimaryHDU(data=data, header=header)
            newhdu.writeto(corrected_fitsfile, overwrite=True)
        
        #Apply corrections to preliminary catalogues
        component_info['Total_flux'] *= flux_correction
        component_info['E_Total_flux'] *= flux_correction
        component_info['Peak_flux'] *= flux_correction
        component_info['E_Peak_flux'] *= flux_correction
        component_info['rms'] *= flux_correction
        component_info['RA'] -= ((ra_offset * u.arcsec).to(u.deg).value) / np.cos(np.deg2rad(component_info['DEC']))
        component_info['DEC'] -= (dec_offset * u.arcsec).to(u.deg).value
        
        source_info['Total_flux'] *= flux_correction
        source_info['E_Total_flux'] *= flux_correction
        source_info['Peak_flux'] *= flux_correction
        source_info['E_Peak_flux'] *= flux_correction
        source_info['Isl_rms'] *= flux_correction
        source_info['RA'] -= ((ra_offset * u.arcsec).to(u.deg).value) / np.cos(np.deg2rad(source_info['DEC']))
        source_info['DEC'] -= (dec_offset * u.arcsec).to(u.deg).value
        
        #Plot current field and find sources associated with LoTSS source
        inside_mask, seren_mask = plot_source(source, corrected_fitsfile, final_plotfolder, lotss_dir, component_info, imcat, lotss_idx, flagged)
        inside_mask = inside_mask | np.array(component_info['S_Code'] == 'X')
        
        #Split out catalogues into source-associated components and serendipitous detections
        component_info.sort('source_id')
        assoc_comps = component_info[inside_mask]
        assoc_src_idx = np.array(list(set(assoc_comps['source_id'])))
        
        assoc_idx = np.array([idx for idx, source_id in enumerate(source_info['source_id']) if source_id in assoc_src_idx])
        assoc_srcs = source_info[assoc_idx]
                
        if len(assoc_srcs) == 0:
            assoc_srcs = add_zero_source(Table(), float(lotss_info['RA'][0]), float(lotss_info['DEC'][0]))
        else:
            #Modify and complete source-associated components
            assoc_srcs, assoc_comps = reorder_cat(assoc_srcs, assoc_comps)
                
        assoc_comps.add_column([pointing]*len(assoc_comps), index=0, name="Pointing")
        assoc_comps.add_column([source]*len(assoc_comps), index=0, name="LoTSS_ID")
        assoc_comps.add_column([flagged]*len(assoc_comps), name="Flagged")
            
        #Combine into source-based catalog
        src_collapsed = collapse_source(assoc_srcs)
        src_collapsed.add_column([pointing], index=0, name="Pointing")
        src_collapsed.add_column([source], index=0, name="LoTSS_ID")
        src_collapsed.add_column([flagged], name="Flagged")
        
        #Apply filters
        mask = assoc_comps['S_Code'] != 'X'
        assoc_comps = assoc_comps[mask]
    
        #Append source to pointing catalogue
        pointing_cat = vstack([pointing_cat, src_collapsed])

        #Append gaussians to catalogue
        if len(assoc_comps) > 0:
            gaussian_cat = vstack([gaussian_cat, assoc_comps])
            
        #Now the serendipitous detections
        if not flagged:
            mask = ~inside_mask & ~seren_mask
            seren_comps = component_info[mask]
            seren_src_idx = np.array(list(set(seren_comps['source_id'])))
            
            seren_idx = np.array([idx for idx, source_id in enumerate(source_info['source_id']) if source_id in seren_src_idx])
            seren_srcs = source_info[seren_idx]
            
            seren_srcs, seren_comps = reorder_cat(seren_srcs, seren_comps)

            try:
                max_src_idx = np.max(seren_src_cat['source_id'])
            except KeyError:
                #Needed for first iteration
                max_src_idx = 0
            seren_srcs['source_id'] += max_src_idx + 1
            seren_comps['source_id'] += max_src_idx + 1
            
            #Append source to serendipitous source catalogue
            if len(seren_srcs) > 0:
                seren_src_cat = vstack([seren_src_cat, seren_srcs])

            #Append gaussians to serendipitous catalogue
            if len(seren_comps) > 0:
                seren_gauss_cat = vstack([seren_gauss_cat, seren_comps])
    
    #Sort final catalogues
    pointing_cat.sort('LoTSS_ID')
    gaussian_cat.sort(['LoTSS_ID','source_id','gaus_num'])
    
    seren_src_cat.add_column([pointing]*len(seren_src_cat), index=0, name="Pointing")
    seren_gauss_cat.add_column([pointing]*len(seren_gauss_cat), index=0, name="Pointing")
    
    #Write out catalogues
    pointing_cat.write(pointing_catfile, format=cat_format, overwrite=update)
    gaussian_cat.write(gaussian_catfile, format=cat_format, overwrite=update)
    seren_src_cat.write(seren_src_catfile, format=cat_format, overwrite=update)
    seren_gauss_cat.write(seren_gauss_catfile, format=cat_format, overwrite=update)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument( dest='pointing', type=str )
    parser.add_argument( '--outdir',type=str,default='catalogue',help='directory to store outputs')
    parser.add_argument( '--postprocessdir',type=str,default='postprocessing',help='directory to store intermediate products')
    # parser.add_argument( '--catalogue-filestem',type=str,default='catalogue.fits',help='suffix of catalogue file' )
    parser.add_argument( '--catalogue-format',type=str,default='csv',help='file format for catalogue (csv, fits)' )
    parser.add_argument( '--inspection',type=str,default='default',help='inspection catalogue file' )
    parser.add_argument( '--update', action='store_true', default=True )
    args = parser.parse_args()
    main( args.pointing, outdir=args.outdir, processingdir=args.postprocessdir, inspection_file=args.inspection, cat_format=args.catalogue_format, update=args.update)
