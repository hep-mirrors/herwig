// -*- C++ -*-
#ifndef HERWIG_CluHadConfig_H
#define HERWIG_CluHadConfig_H
/*! \file CluHadConfig.h 
 *  \brief This file contains the typedef declarations used in Hadronization
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *  \ingroup Hadronization
 *
 * This is the declaration of the CluHadConfig.h header file.
 *
 * Handy header file to be included in all Hadronization classes. <BR>
 * It contains only some useful typedefs.
 *
 * See also:
 * Herwig.h
 */ 

#include "Herwig++/Config/Herwig.h"

namespace Herwig { 

  using namespace ThePEG;

  class Cluster;
  typedef Ptr<Cluster>::pointer ClusterPtr;
  typedef Ptr<Cluster>::transient_pointer tClusterPtr;
  typedef vector<ClusterPtr> ClusterVector;
  typedef vector<tClusterPtr> tClusterVector;

  class GlobalParameters;
  typedef Ptr<GlobalParameters>::pointer GlobParamPtr;
  typedef Ptr<GlobalParameters>::transient_pointer tGlobParamPtr;

  class PartonSplitter;
  typedef Ptr<PartonSplitter>::pointer PartonSplitterPtr;
  typedef Ptr<PartonSplitter>::transient_pointer tPartonSplitterPtr;

  class ClusterFinder;
  typedef Ptr<ClusterFinder>::pointer ClusterFinderPtr;
  typedef Ptr<ClusterFinder>::transient_pointer tClusterFinderPtr;

  class ColourReconnector;
  typedef Ptr<ColourReconnector>::pointer ColourReconnectorPtr;
  typedef Ptr<ColourReconnector>::transient_pointer tColourReconnectorPtr; 
 
  class ClusterFissioner;
  typedef Ptr<ClusterFissioner>::pointer ClusterFissionerPtr;
  typedef Ptr<ClusterFissioner>::transient_pointer tClusterFissionerPtr;

  class LightClusterDecayer;
  typedef Ptr<LightClusterDecayer>::pointer LightClusterDecayerPtr;
  typedef Ptr<LightClusterDecayer>::transient_pointer tLightClusterDecayerPtr;

  class ClusterDecayer;
  typedef Ptr<ClusterDecayer>::pointer ClusterDecayerPtr;
  typedef Ptr<ClusterDecayer>::transient_pointer tClusterDecayerPtr;

  class HadronSelector;
  typedef Ptr<HadronSelector>::pointer HadronSelectorPtr;
  typedef Ptr<HadronSelector>::transient_pointer tHadronSelectorPtr;
} // end Herwig namespace

/*! \defgroup Hadronization */

#endif // HERWIG_CluHadConfig_H 



