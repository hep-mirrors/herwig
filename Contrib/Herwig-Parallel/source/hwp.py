#! /usr/bin/env python

import string
import sys
import os
from ConfigParser import SafeConfigParser

"""read a file returning the lines in reverse order for each call of readline()
This actually just reads blocks (4096 bytes by default) of data from the end of
the file and returns last line in an internal buffer.  I believe all the corner
cases are handled, but never can be sure..."""

class BackwardsReader:
  def readline(self):
    while len(self.data) == 1 and ((self.blkcount * self.blksize) < self.size):
      self.blkcount = self.blkcount + 1
      line = self.data[0]
      try:
        self.f.seek(-self.blksize * self.blkcount, 2) # read from end of file
        self.data = string.split(self.f.read(self.blksize) + line, '\n')
      except IOError:  # can't seek before the beginning of the file
        self.f.seek(0)
        self.data = string.split(self.f.read(self.size - (self.blksize * (self.blkcount-1))) + line, '\n')

    if len(self.data) == 0:
      return ""

    # self.data.pop()
    # make it compatible with python <= 1.5.1
    line = self.data[-1]
    self.data = self.data[:-1]
    return line + '\n'

  def close(self):
    self.f.close()

  def __init__(self, file, blksize=4096):
    """initialize the internal structures"""
    # get the file size
    self.size = os.stat(file)[6]
    # how big of a block to read from the file...
    self.blksize = blksize
    # how many blocks we've read
    self.blkcount = 1
    self.f = open(file, 'rb')
    # if the file is smaller than the blocksize, read a block,
    # otherwise, read the whole thing...
    if self.size > self.blksize:
      self.f.seek(-self.blksize * self.blkcount, 2) # read from end of file
    self.data = string.split(self.f.read(self.blksize), '\n')
    # strip the last item if it's empty...  a byproduct of the last line having
    # a newline at the end of it
    if not self.data[-1]:
      # self.data.pop()
      self.data = self.data[:-1]


def assertRunIsNotCompressed(runName):
  if os.path.isfile(runName+'/Jobs.tar.gz'):
    sys.stderr.write("The specified run "+runName+" is already compressed!\n")
    sys.exit(1)

def checkConfig(configParserClusters, configParserQueues):
  clusters = []
  for cluster in configParserClusters.sections():
    clusters.append(cluster)
    options = configParserClusters.options(cluster)
    if not 'joblist' in options:
      sys.stderr.write("The cluster "+cluster+" does not have the option 'joblist'! Please correct the configuration file Herwig-Parallel/config/clusters.conf accordingly.")
      sys.exit(1)
    if not 'status' in options:
      sys.stderr.write("The cluster "+cluster+" does not have the option 'status'! Please correct the configuration file Herwig-Parallel/config/clusters.conf accordingly.")
      sys.exit(1)
    if not 'statusqueued' in options:
      sys.stderr.write("The cluster "+cluster+" does not have the option 'statusQueued'! Please correct the configuration file Herwig-Parallel/config/clusters.conf accordingly.")
      sys.exit(1)
    if not 'statusrunning' in options:
      sys.stderr.write("The cluster "+cluster+" does not have the option 'statusRunning'! Please correct the configuration file Herwig-Parallel/config/clusters.conf accordingly.")
      sys.exit(1)
    if not 'abort' in options:
      sys.stderr.write("The cluster "+cluster+" does not have the option 'abort'! Please correct the configuration file Herwig-Parallel/config/clusters.conf accordingly.")
      sys.exit(1)
  if len(clusters) == 0:
    sys.stderr.write("Please specify at least one cluster in Herwig-Parallel/config/clusters.conf!")
    sys.exit(1)

  queues = []
  for queue in configParserQueues.sections():
    queues.append(queue)
    options = configParserQueues.options(queue)
    if not 'cluster' in options:
      sys.stderr.write("The queue "+queue+" does not have the option 'cluster'! Please correct the configuration file Herwig-Parallel/config/queues.conf accordingly.")
      sys.exit(1)
    if not configParserQueues.get(queue,'cluster') in clusters:
      sys.stderr.write("The cluster "+configParserQueues.get(queue,'cluster')+" specified for the queue "+queue+" was not configured in the configuration file Herwig-Parallel/config/clusters.conf! Please correct either the cluster or the queue configuration file.")
      sys.exit(1)
    if not 'submit' in options:
      sys.stderr.write("The queue "+queue+" does not have the command 'submit'! Please correct the configuration file Herwig-Parallel/config/queues.conf accordingly.")
      sys.exit(1)
  if len(queues) == 0:
    sys.stderr.write("Please specify at least one queue in Herwig-Parallel/config/queues.conf!")
    sys.exit(1)
  return(queues)

def addToIndex(JobID, Path):
  home = os.path.expanduser('~')
  index = home + '/.Herwig-Parallel/index'
  if not os.path.exists(home+'/.Herwig-Parallel'):
    os.makedirs(home+'/.Herwig-Parallel')
  f = open(index,'a')
  f.write('{} {}\n'.format(JobID, Path))
  f.close()
  return(True)

def findInIndex(JobIDs):
  if len(JobIDs) == 0: return(None)
  home = os.path.expanduser('~')
  index = home + '/.Herwig-Parallel/index'

  if not os.path.isfile(index):
    sys.stderr.write('! Could not find the run index which contains the job ids.\n')
    sys.stderr.write('! The run index will be filled once you start runs with Herwig-Parallel.\n')
    sys.stderr.write("! It is located at '~/.Herwig-Parallel/index'\n")
    sys.stderr.write('\n')
    return(None)

  f = open(index,'r')
  pairs = []
  for line in f:
    job = line.split()[0]
    if job in JobIDs:
      pair = [job, line.replace(job,'').replace('\n','').strip()]
      pairs.insert(0, pair)
  f.close()

  paths = []
  for JobID in JobIDs:
    item = []
    for pair in pairs:
      if pair[0] == JobID:
        item.append(pair[1])
    paths.append([JobID, item])

  return(paths)

def checkIndex():

  nLines = 25000 # keep the latest 25000 jobs in index
  home = os.path.expanduser('~')
  index = home + '/.Herwig-Parallel/index'

  if not os.path.isfile(index):
    sys.stderr.write('! Could not find the run index.\n')
    sys.stderr.write('! The run index will be filled once you start runs with Herwig-Parallel.\n')
    sys.stderr.write("! It is located at '~/.Herwig-Parallel/index'\n")
    return(False)

  size = os.path.getsize(index)
  if size > 10000000: # = 10 MB
    print('! The size of the run index exceeded 10 MB.')
    print("! It's current size is {:.2f} MB.".format(float(size)/(1000000)))
    print('! Therefore the run index will now be trimmed down to {} lines.'.format(nLines))
    print('')
    result = trimIndex(nLines)
    if result == True:
      print('> Trimmed down run index successfully.')
      print('')
      return(True)
    else:
      return(False)

def trimIndex(nLines):

  home = os.path.expanduser('~')
  index = home + '/.Herwig-Parallel/index'

  try:
    # copy old index file to have a backup
    os.rename(index, index+'.orig')

    # read the last nLines lines from the old index file
    lines = []
    fo = BackwardsReader(index+'.orig')
    for l in range(nLines):
      lines.insert(0, fo.readline())
    fo.close()

    # write new (trimmed) index file
    fn = open(index, 'w')
    for l in lines:
      fn.write(l)
    fn.close()

    # remove old index file and return True if everything has worked out alright
    os.remove(index+'.orig')
    return(True)

  except:
    return(False)