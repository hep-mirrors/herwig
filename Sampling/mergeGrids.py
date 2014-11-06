#! /usr/bin/env python

import os
import sys
import glob

"""\
%prog [--setupfile=FILE] RUNNAME

Combine Herwig++ grid files
"""

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage=__doc__)

    parser.add_option('-x', '--setupfile', type='string',
                      help='Specify the setup file which has been used.',
                      default='',
                      dest='setupFile')

    opts, args = parser.parse_args()

    if len(args) < 1:
        sys.stderr.write('Please specify a run name\n')
        sys.exit(1)

    gridId = 'Matchbox/' + args[0] + '/' + args[0]
    setupName=opts.setupFile
    if setupName:
        gridId = gridId + '-' + setupName

    gridFiles=glob.glob(gridId + '-integrationJob*-grids.xml')

    if not gridFiles:
        sys.stderr.write('No grid files have been found to combine\n')
        sys.exit(1)        

    gridCombined = open(gridId + '-grids.xml','w')
    gridCombined.write('<Grids>\n')

    for gridFile in gridFiles:
        grid = open(gridFile,'r')
        gridContent = grid.read()
        gridContent = gridContent.replace('<Grids>','')
        gridContent = gridContent.replace('</Grids>','')
        gridCombined.write(gridContent)
        grid.close()

    gridCombined.write('</Grids>\n')

