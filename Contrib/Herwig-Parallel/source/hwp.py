#! /usr/bin/env python

import sys
import os
from ConfigParser import SafeConfigParser

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
