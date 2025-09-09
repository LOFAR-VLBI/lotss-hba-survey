from astropy.table import Table, vstack, join
from astropy.io import fits
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import os, glob, re
import matplotlib
from matplotlib import pyplot as plt
import argparse

##############################################################
## Plotting housekeeping

matplotlib.rcParams['legend.frameon'] = False
matplotlib.rcParams['axes.labelsize'] = 'large'
## set up some colours
n = 255
mycols = plt.cm.viridis(np.linspace(0, 1,n))
mycols_m = plt.cm.inferno(np.linspace(0, 1,n))

##############################################################

def main( pointing, phasediff_file ):
    phasediff = Table.read(phasediff_file,format='csv')

    ## find a good cutoff
    test_cutoffs = np.arange(0.4,2.4,0.1)[::-1]
    for tc in test_cutoffs:
        ddcal_idx = np.where(phasediff['spd_score'] < tc )[0]
        recal_idx = np.where(phasediff['spd_score'] >= tc )[0]
        ddcal = phasediff[ddcal_idx]
        recal = phasediff[recal_idx]
        ddcal_coords = SkyCoord( ddcal['RA'], ddcal['DEC'], unit='deg' )
        recal_coords = SkyCoord( recal['RA'], recal['DEC'], unit='deg' )
        test_separations = []
        for i in np.arange(0,len(recal)):
            coords = recal_coords[i]
            seps = coords.separation(ddcal_coords)
            test_separations.append(np.min(seps).value)
        print(np.max(test_separations))
        if np.max(test_separations) < 0.5:
            cutoff = tc

    ## plotting
    fig = plt.figure(figsize=(20,10))
    p1 = plt.axes([0.05,0.1,sbsizex*5/10,sbsizey])
    p1.hist( phasediff['spd_score'], bins=20 )
    p1.axvline(x=cutoff,color='black',linewidth=3)
    p1.text( 0.1,0.9,'Cutoff={:.2}'.format(cutoff), transform=p1.transAxes, size=16)
    p1.xaxis.set_tick_params(labelsize=30)
    p1.yaxis.set_tick_params(labelsize=30)
    p2 = plt.axes([0.12+sbsizex*5/10,0.1,sbsizex*5/10,sbsizey])
    p2.scatter( phasediff['RA'], phasediff['DEC'], color='gray' )
    p2.scatter( phasediff['RA'][np.where(phasediff['spd_score'] < cutoff)[0]], phasediff['DEC'][np.where(phasediff['spd_score'] < cutoff)[0]], marker='X', color='red', s=80 )
    p2.xaxis.set_tick_params(labelsize=30)
    p2.yaxis.set_tick_params(labelsize=30)
    fig.savefig('phasediff_distribution.png')
    plt.close()
    fig.clear()

    ## now generate the final list
    ddcal_idx = np.where(phasediff['spd_score'] < cutoff )[0]
    recal_idx = np.where(phasediff['spd_score'] >= cutoff )[0]

    ddcal = phasediff[ddcal_idx]
    recal = phasediff[recal_idx]

    ddcal_coords = SkyCoord( ddcal['RA'], ddcal['DEC'], unit='deg' )
    recal_coords = SkyCoord( recal['RA'], recal['DEC'], unit='deg' )

    nearest = []
    separations = []
    for i in np.arange(0,len(recal)):
        source = recal['source'][i].split('_')[0]
        coords = recal_coords[i]
        seps = coords.separation(ddcal_coords)
        nearest_idx = np.where( seps == np.min(seps) )[0][0]
        nearest.append(ddcal['source'][nearest_idx].split('_')[0])
        separations.append(np.min(seps).value)

    outfile = os.path.join(os.getenv('DATA_DIR'),pointing,'recalibration_list.csv')
    with open(outfile,'w') as f:
        f.write('Source,Nearest,separation\n')
        for i in np.arange(0,len(recal)):
            f.write('{:s},{:s},{:.4}\n'.format(recal['source'][i].split('_')[0],nearest[i],separations[i]))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('pointing', type=str, help='name of the pointing')
    parser.add_argument('phasediff', type=str, help='phasediff file path')

    args = parser.parse_args()
    main( args.pointing, args.phasediff )



