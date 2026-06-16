import glob
import os
import argparse
import re

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def get_selfcal_files( directory ):
    image_files = glob.glob( os.path.join( directory, 'image*MFS-image.fits' ) )
    image_pb_files = glob.glob( os.path.join( directory, 'image*MFS-image-pb.fits' ) )
    image_pngs = glob.glob( os.path.join( directory, '*png' ) )
    losoto = glob.glob( os.path.join( directory, 'plotlosoto*' ) )
    src = os.path.basename(losoto[0]).replace('plotlosoto','').split('_')[0]
    h5_files = natural_sort( glob.glob( os.path.join( directory, 'merged*{:s}*h5'.format( src ) ) ) )
    logfiles = glob.glob( os.path.join( directory, '*.log') )
    keep_files = image_files + image_pb_files + image_pngs + [h5_files[-1]] + losoto + logfiles
    return( keep_files )


def main( selfcaldir ):
    sources = glob.glob( os.path.join( selfcaldir, 'ILTJ*' ) )
    pointing = sources[0].split('/')[-4]
    recal_sources = glob.glob( os.path.join( os.getenv('DATA_DIR'), 'processing',pointing,'*', 'ILTJ*ms' ) )
    for recal_source in recal_sources:
        src = os.path.basename(recal_source).split('_')[0]
        src_idx = [ i for i, val in enumerate(sources) if src in val ]
        src_dir = sources[src_idx[0]]
        src_files = glob.glob( os.path.join( src_dir, '*' ) )
        os.makedirs(os.path.join(src_dir,'first_selfcal'))
        for srcf in src_files:
            if 'ILTJ' not in os.path.basename(srcf):
                os.system('mv {:s} {:s}/'.format(srcf,os.path.join(src_dir,'first_selfcal')) )
            if 'plotlosoto' in srcf or 'merged' in srcf:
                os.system('mv {:s} {:s}/'.format(srcf,os.path.join(src_dir,'first_selfcal')) )
        new_files = get_selfcal_files( os.path.dirname( recal_source ) )
        for nf in new_files:
            os.system( 'mv {:s} {:s}/'.format( nf, src_dir ) )
        os.system( 'rm -r {:s}'.format( os.path.dirname( recal_source ) ) )
        
        
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( dest='selfcaldir', type=str )
    args = parser.parse_args()
    main( args.selfcaldir )

