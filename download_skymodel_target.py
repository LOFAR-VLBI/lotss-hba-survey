#!/usr/bin/env python
import astropy.units as u
import os,sys
import glob
import casacore.tables as pt
import mocpy
import numpy as np
import time
import lsmtool
import lsmtool.skymodel
import subprocess

########################################################################
def grab_coord_MS(MS):
    """
    Read the coordinates of a field from one MS corresponding to the selection given in the parameters

    Parameters
    ----------
    MS : str
        Full name (with path) to one MS of the field

    Returns
    -------
    RA, Dec : "tuple"
        coordinates of the field (RA, Dec in deg , J2000)
    """

    # reading the coordinates ("position") from the MS
    # NB: they are given in rad,rad (J2000)
    [[[ra,dec]]] = pt.table(MS+'::FIELD', readonly=True, ack=False).getcol('PHASE_DIR')

    # RA is stocked in the MS in [-pi;pi]
    # => shift for the negative angles before the conversion to deg (so that RA in [0;2pi])
    if ra<0:
        ra=ra+2*np.pi

    # convert radians to degrees
    ra_deg =  ra/np.pi*180.
    dec_deg = dec/np.pi*180.

    # and sending the coordinates in deg
    return(ra_deg,dec_deg)


########################################################################
def input2strlist_nomapfile(invar):
    """ from bin/download_IONEX.py
    give the list of MSs from the list provided as a string
    """

    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            str_list = [invar.strip(' \'\"')]
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return(str_list)


########################################################################
def main(ms_input, SkymodelPath, Radius=5., DoDownload="True", Source="TGSS", targetname = "pointing", fluxlimit = None):
    """
    Download the skymodel for the target field

    Parameters
    ----------
    ms_input : str
        String from the list (map) of the target MSs
    SkymodelPath : str
        Full name (with path) to the skymodel; if YES is true, the skymodel will be downloaded here
    Radius : string with float (default = "5.")
        Radius for the skymodel cone search in degrees
    DoDownload : str ("Force" or "True" or "False")
        Do or do not download the skymodel
        "Force": download skymodel and delete existing skymodel if needed.
        "True" or "Yes": use existing skymodel file if it exists, download skymodel
                         if it does not.
        "False" or "No": Do not download skymodel, raise an exception if skymodel
                         file does not exist.
    targetname : str
        Give the patch a certain name, default: "pointing"
    """

    FileExists = os.path.isfile(SkymodelPath)
    if (not FileExists and os.path.exists(SkymodelPath)):
        raise ValueError("download_skymodel_target: Path: \"%s\" exists but is not a file!"%(SkymodelPath))
    download_flag = False
    if not os.path.exists(os.path.dirname(SkymodelPath)):
        os.makedirs(os.path.dirname(SkymodelPath))
    if DoDownload.upper() == "FORCE":
        if FileExists:
            os.remove(SkymodelPath)
        download_flag = True
    elif DoDownload.upper() == "TRUE" or DoDownload.upper() == "YES":
        if FileExists:
            print("USING the exising skymodel in "+ SkymodelPath)
            return(0)
        else:
            download_flag = True
    elif DoDownload.upper() == "FALSE" or DoDownload.upper() == "NO":
         if FileExists:
            print("USING the exising skymodel in "+ SkymodelPath)
            return(0)
         else:
            raise ValueError("download_skymodel_target: Path: \"%s\" does not exist and skymodel download is disabled!"%(SkymodelPath))

    # If we got here, then we are supposed to download the skymodel.
    assert download_flag is True # Jaja, belts and suspenders...
    print("DOWNLOADING skymodel for the target into "+ SkymodelPath)

    # Reading a MS to find the coordinate (casacore)
    [RATar,DECTar]=grab_coord_MS(input2strlist_nomapfile(ms_input)[0])
    RATar *= u.deg
    DECTar *= u.deg
    Radius *= u.deg

    if Source == 'LOTSS':
        print('Checking LoTSS coverage for the requested centre and radius.')
        mocpath = os.path.join(os.path.dirname(SkymodelPath), 'dr3-moc.moc')
        subprocess.run(['wget', 'https://lofar-surveys.org/public/DR3/catalogues/dr3-moc.moc', '-O', mocpath], capture_output=True, check=True)
        moc = mocpy.MOC.from_fits(mocpath)
        print(RATar, DECTar, Radius)
        covers_centre = moc.contains(RATar, DECTar)
        # Checking single coordinates, so get rid of the arRATary.
        covers_left = moc.contains(RATar - Radius, DECTar)[0]
        covers_right = moc.contains(RATar + Radius, DECTar)[0]
        covers_bottom = moc.contains(RATar, DECTar - Radius)[0]
        covers_top = moc.contains(RATar, DECTar + Radius)[0]
        if covers_centre and not (covers_left and covers_right and covers_bottom and covers_top):
            raise ValueError('Incomplete LoTSS coverage for the requested centre and radius!')
        elif not covers_centre:
            raise ValueError('Requested centre not covered by LoTSS!')
        elif not covers_centre and not (covers_left and covers_right and covers_bottom and covers_top):
            raise ValueError('No LoTSS coverage for the requested centre and radius!')
        else:
            print('Complete LoTSS coverage for the requested centre and radius.')

    # Downloading the skymodel, skip after five tries
    errorcode = 1
    tries     = 0
    while errorcode != 0 and tries < 5:
        if Source == 'TGSS':
            errorcode = os.system("wget -O "+SkymodelPath+ " \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv5.cgi?coord="+str(RATar.value)+","+str(DECTar.value)+"&radius="+str(Radius.value)+"&unit=deg&deconv=y\' ")
        elif Source == 'GSM':
            errorcode = os.system("wget -O "+SkymodelPath+ " \'https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord="+str(RATar.value)+","+str(DECTar.value)+"&radius="+str(Radius.value)+"&unit=deg&deconv=y\' ")
        elif Source == 'LOTSS':
            try:
                lotssmodel = lsmtool.skymodel.SkyModel('lotss', VOPosition=[RATar.value, DECTar.value], VORadius=Radius.value)
                lotssmodel.write(SkymodelPath)
                if len(lotssmodel) > 0:
                    errorcode = 0
            except ConnectionError:
                if tries == 5:
                    raise IOError('Download of LoTSS sky model failed after {} attempts.'.format(5))
                else:
                    print('Attempt #{0:d} to download LoTSS sky model failed. Attempting '
                                 '{1:d} more times.'.format(tries, 5 - tries))
        time.sleep(5)
        tries += 1

    if not os.path.isfile(SkymodelPath):
        raise IOError("download_skymodel_target: Path: \"%s\" does not exist after trying to download the skymodel."%(SkymodelPath))

    # Treat all sources as one group (direction)
    skymodel = lsmtool.load(SkymodelPath)
    if fluxlimit:
        skymodel.remove('I<' + str(fluxlimit))
    skymodel.group('single', root = targetname)
    skymodel.write(clobber=True)
    
    return(0)


########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=' Download the skymodel for the target field')

    parser.add_argument('MSfile', type=str, nargs='+',
                        help='One (or more MSs) for which a skymodel will be download.')
    parser.add_argument('SkyTar', type=str,
                        help='Full name (with path) to the skymodel; the skymodel will be downloaded here')
    parser.add_argument('--Radius', type=float, default=5.,
                        help='Radius for the skymodel cone search in degrees')
    parser.add_argument('--Source', type=str, default='TGSS',
                        help='Choose source for skymodel: TGSS, GSM or LOTSS')
    parser.add_argument('--DoDownload', type=str, default="True",
                        help='Do or do not download the skymodel ("Force" or "True" or "False").')
    parser.add_argument('--targetname', type=str, default='pointing',
                        help='Name of the patch of the skymodel')
    parser.add_argument('--fluxlimit', type=float, default=None,
                        help='Remove sources from the skymodel below the given fluxlimit [Jy],')

    args = parser.parse_args()

    main(args.MSfile, args.SkyTar, args.Radius, args.DoDownload, args.Source, args.targetname, args.fluxlimit)
