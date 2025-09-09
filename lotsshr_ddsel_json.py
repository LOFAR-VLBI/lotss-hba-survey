import glob
import argparse
import os

def main( selfcaldir ):
    softwaredir = os.getenv('SOFTWAREDIR')
    msins = glob.glob( os.path.join( selfcaldir, 'ILTJ*/ILTJ*' ) )
    with open( 'config.json', 'w' ) as f:
        f.write('{\n')
        f.write('  "msin": [\n')
        for msin in msins:
            f.write('    {\n')
            f.write('      "class": "Directory",\n')
            f.write('      "path": "{:s}"\n'.format(msin))
            ss = '    },\n'
            if msin == msins[-1]:
                f.write(ss.replace(',',''))
            else:
                f.write(ss)
        f.write('  ]\n')
#        f.write('  ],\n')
#        f.write('  "selfcal": {\n')
#        f.write('    "class": "Directory",\n')
#        f.write('    "path": "{:s}"\n'.format(os.path.join(softwaredir,'lofar_facet_selfcal')))
#        f.write('  },\n')
#        f.write('  "lofar_helpers": {\n')
#        f.write('    "class": "Directory",\n')
#        f.write('    "path": "{:s}"\n'.format(os.path.join(softwaredir,'lofar_helpers')))
#        f.write('  }\n')
        f.write('}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument( dest='selfcaldir', type=str )
    args = parser.parse_args()
    main( args.selfcaldir )


    

