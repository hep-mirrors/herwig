#! /usr/bin/env python

import os
import sys
import glob

"""\
%prog [--setupfile=FILE] [--tag=TAG] RUNNAME

Combine Herwig grid files
"""

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage=__doc__)

    parser.add_option('-x', '--setupfile', type='string',
                      help='Specify the setup file which has been used.',
                      default='',
                      dest='setupFile')

    parser.add_option('-t', '--tag', type='string',
                      help='Specify the tag name which has been used.',
                      default='',
                      dest='tagName')

    opts, args = parser.parse_args()

    if len(args) < 1:
        sys.stderr.write('Please specify a run name\n')
        sys.exit(1)

    runName=args[0]
    if runName.endswith('.run'):
        runName = runName[:-4]
    setupName=opts.setupFile
    tagName=opts.tagName
    if setupName:
        runName = runName + '/' + setupName
    if tagName:
        runName = runName + '/' + tagName

    gridId = 'Herwig-scratch/' + runName

    # print 'Looking in ' + gridId

    gridFiles=sorted(glob.glob(gridId + '/integrationJob*/HerwigGrids.xml'))

    if not gridFiles:
        sys.stderr.write('No grid files have been found to combine\n')
        sys.exit(1)

    # print gridFiles

    gridCombined = open(gridId + '/HerwigGrids.xml','w')
    gridCombined.write('<Grids>\n')


    for gridFile in gridFiles:
        grid = open(gridFile,'r')
	print(gridFile)
        gridContent = grid.read()
        gridContent = gridContent.replace('<Grids>\n','')
        gridContent = gridContent.replace('</Grids>\n','')
        gridCombined.write(gridContent)
        grid.close()

    gridCombined.write('</Grids>\n')

