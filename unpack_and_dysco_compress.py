#!/usr/bin/python3

import argparse
import os
import casacore.tables as ct

def main( tarfile ):
    basedir = os.path.dirname(tarfile)
    trf = os.path.basename(tarfile)
    os.chdir(basedir)
    os.system( 'tar -xvf {:s}'.format(os.path.basename(tarfile)) )
    msfile = '_'.join(trf.split('_')[0:-1])
    try: 
        with ct.table(msfile) as tab:
            has_dysco = ( tab.getdesc()["DATA"]["dataManagerGroup"] == "DyscoData" )
        if not has_dysco:
            os.rename(msfile, msfile + ".nodysco")
            os.system("DP3 numthreads=1 msin={:s}.nodysco msout={:s} msout.storagemanager=dysco steps=[]".format(msfile, msfile))
            os.system('rm -r {:s}.nodysco'.format(msfile))

        os.remove(trf)
    except:
        with open('failed_{:s}.txt'.format(msfile),'w') as f:
            f.write("{:s} is not a valid MeasurementSet".format(msfile))

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('tarfile',help='tarfile')
    args = parser.parse_args()
    main(args.tarfile)


