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

  class Component;
  typedef Ptr<Component>::pointer CompPtr;
  typedef Ptr<Component>::transient_pointer tCompPtr;
  typedef vector<CompPtr> CollecCompPtr;

  class Cluster;
  typedef Ptr<Cluster>::pointer CluPtr;
  typedef Ptr<Cluster>::transient_pointer tCluPtr;
  typedef vector<CluPtr> CollecCluPtr;

} // end Herwig namespace


#endif // HERWIG_CluHadConfig_H 



