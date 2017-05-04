// -*- C++ -*-
//
// CluHadConfig.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_CluHadConfig_H
#define HERWIG_CluHadConfig_H

#include "ThePEG/Config/ThePEG.h"
#include "Cluster.fh"
#include "ClusterDecayer.fh"
#include "ClusterFissioner.fh"
#include "HadronSelector.fh"

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



