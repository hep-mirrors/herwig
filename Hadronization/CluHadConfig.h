// -*- C++ -*-
#ifndef HERWIG_CluHadConfig_H
#define HERWIG_CluHadConfig_H
//
// This is the declaration of the <!id>CluHadConfig.h<!!id> header file.
//
// CLASSDOC SUBSECTION Description:
//
// Handy header file to be included in all Hadronization classes. <BR>
// It contains only some useful typedefs.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:Herwig.html">Herwig.h</a>,
// 

#include "Herwig++/Config/Herwig.h"

namespace Herwig { 

  using namespace Pythia7;

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

  class HadronsSelector;
  typedef Ptr<HadronsSelector>::pointer HadronsSelectorPtr;
  typedef Ptr<HadronsSelector>::transient_pointer tHadronsSelectorPtr;
} // end Herwig namespace


#endif // HERWIG_CluHadConfig_H 



