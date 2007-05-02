// -*- C++ -*-
#ifndef HERWIG_CluHadConfig_H
#define HERWIG_CluHadConfig_H

#include "Herwig++/Config/Herwig.h"
#include "Cluster.fh"
#include "ClusterDecayer.fh"
#include "ClusterFissioner.fh"

namespace Herwig { 
using namespace ThePEG;

/** \ingroup Hadronization
 *  \brief This file contains the typedef declarations used in Hadronization
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This is the declaration of the CluHadConfig.h header file.
 *
 *  Handy header file to be included in all Hadronization classes. <BR>
 *  It contains only some useful typedefs.
 *
 *  See also:
 *  Herwig.h
 */ 

  /**
   * Typedef for a vector of ClusterPtr
   */
  typedef vector<ClusterPtr> ClusterVector;

  /**
   * Typedef for a vector of tClusterPtr
   */
  typedef vector<tClusterPtr> tClusterVector;

}

#endif // HERWIG_CluHadConfig_H 



