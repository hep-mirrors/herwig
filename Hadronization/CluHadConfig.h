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

  typedef vector<ClusterPtr> ClusterVector;
  typedef vector<tClusterPtr> tClusterVector;

  class PartonSplitter;
  typedef Ptr<PartonSplitter>::pointer PartonSplitterPtr;
  typedef Ptr<PartonSplitter>::transient_pointer tPartonSplitterPtr;

  class ClusterFinder;
  typedef Ptr<ClusterFinder>::pointer ClusterFinderPtr;
  typedef Ptr<ClusterFinder>::transient_pointer tClusterFinderPtr;

  class ColourReconnector;
  typedef Ptr<ColourReconnector>::pointer ColourReconnectorPtr;
  typedef Ptr<ColourReconnector>::transient_pointer tColourReconnectorPtr; 

  class LightClusterDecayer;
  typedef Ptr<LightClusterDecayer>::pointer LightClusterDecayerPtr;
  typedef Ptr<LightClusterDecayer>::transient_pointer tLightClusterDecayerPtr;

  class HadronSelector;
  typedef Ptr<HadronSelector>::pointer HadronSelectorPtr;
  typedef Ptr<HadronSelector>::transient_pointer tHadronSelectorPtr;

} // end Herwig namespace

#endif // HERWIG_CluHadConfig_H 



