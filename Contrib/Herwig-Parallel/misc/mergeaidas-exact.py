#! /usr/bin/env python

import sys, os, copy
from math import sqrt
from subprocess import Popen
import lighthisto

## Try to load faster but non-standard cElementTree module
try:
    import xml.etree.cElementTree as ET
except ImportError:
    try:
        import cElementTree as ET
    except ImportError:
        try:
            import xml.etree.ElementTree as ET
        except:
            sys.stderr.write("Can't load the ElementTree XML parser: please install it!\n")
            sys.exit(1)

from optparse import OptionParser
parser = OptionParser(usage="%prog aidafile [aidafile2 ...]")
parser.add_option("-o", "--outfile", dest="OUTFILE",
                  default="mergedrew.aida", help="file for merged aida output.")
parser.add_option("-s", "--sum",
                  action="store_true", dest="performSum", default=False,
                  help="sum the bin values instead of averaging")
opts, args = parser.parse_args()
headerprefix = ""

if len(args) < 1:
    sys.stderr.write("Must specify at least one AIDA histogram file\n")
    sys.exit(1)


try:
    outaida = open(opts.OUTFILE, "w")
except:
    sys.stderr.write("Couldn't open outfile %s for writing." % opts.OUTFILE)

try:
    outdat = open(opts.OUTFILE.replace(".aida", ".dat"), "w")
except:
    sys.stderr.write("Couldn't open outfile %s for writing." % opts.OUTFILE.replace(".aida", ".dat"))

# Get attempts for all jobs
attempts = []
ATTEMPTS = 0
f_job_results = open('job-results.out','r')
for line in f_job_results:
  n = int(line.split()[1])
  attempts.append(n)
  ATTEMPTS += n
f_job_results.close()

## Get histos
inhistos = {}
weights = {}
for aidafile in args:
    tree = ET.parse(aidafile)
    for dps in tree.findall("dataPointSet"):
        h = lighthisto.Histo.fromDPS(dps)
        if not inhistos.has_key(h.fullPath()):
            inhistos[h.fullPath()] = {}
        inhistos[h.fullPath()][aidafile] = h

## Merge histos
outhistos = {}
for path, hs in inhistos.iteritems(): # loop over the observables
    #print path, hs
    outhistos[path] = copy.deepcopy(hs.values()[0])
    for i, b in enumerate(outhistos[path].getBins()): # loop over the bins
        j = 0
        ns = 0.0
        ns2 = 0.0
        n2d2 = 0.0
        for infile, h in hs.iteritems(): # loop over the jobs
          ns += attempts[j]*h.getBin(i).val
          ns2 += attempts[j]*h.getBin(i).val*h.getBin(i).val
          n2d2 += attempts[j]*(attempts[j]-1)*h.getBin(i).getErr()*h.getBin(i).getErr()
          j += 1
        outhistos[path].getBin(i).val = ns / float(ATTEMPTS)
        outhistos[path].getBin(i).setErr(sqrt(n2d2 + ns2 - ns*ns/float(ATTEMPTS)) / sqrt(float(ATTEMPTS*(ATTEMPTS-1))))

## Write out merged histos
#print sorted(outhistos.values())
outdat.write("\n\n".join([h.asFlat() for h in sorted(outhistos.values())]))
outdat.write("\n")
outdat.close()

Popen(["flat2aida", opts.OUTFILE.replace(".aida", ".dat")], stdout=outaida).wait()

os.unlink(opts.OUTFILE.replace(".aida", ".dat"))
