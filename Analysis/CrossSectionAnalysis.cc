// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CrossSectionAnalysis class.
//

#include "CrossSectionAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "Herwig/Sampling/GeneralSampler.h"

#include "Herwig/Utilities/XML/ElementIO.h"

using namespace Herwig;

CrossSectionAnalysis::CrossSectionAnalysis() {}

CrossSectionAnalysis::~CrossSectionAnalysis() {}

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

void CrossSectionAnalysis::dofinish() {

  AnalysisHandler::dofinish();

  Ptr<StandardEventHandler>::tptr seh =
    dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(generator()->eventHandler());
  Ptr<GeneralSampler>::tptr sampler =
    dynamic_ptr_cast<Ptr<GeneralSampler>::tptr>(seh->sampler());

  unsigned long attemptedPoints = sampler->attempts();
  double sumOfWeights = sampler->sumWeights();
  double sumOfSquaredWeights = sampler->sumWeights2();
  CrossSection maxXSection = sampler->maxXSec();

  XML::Element elem(XML::ElementTypes::Element,"Run");

  elem.appendAttribute("name",generator()->runName());
  elem.appendAttribute("attemptedPoints",attemptedPoints);
  elem.appendAttribute("sumOfWeights",sumOfWeights*maxXSection/picobarn);
  elem.appendAttribute("sumOfSquaredWeights",sumOfSquaredWeights*sqr(maxXSection/picobarn));

  XML::Element xhistos(XML::ElementTypes::Element,"Histograms");

  elem.append(xhistos);

  string fname = name() + ".xml";
  ofstream runXML(fname.c_str());
  runXML << setprecision(16);
  XML::ElementIO::put(elem,runXML);

}

IBPtr CrossSectionAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr CrossSectionAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void CrossSectionAnalysis::persistentOutput(PersistentOStream &) const {}

void CrossSectionAnalysis::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<CrossSectionAnalysis,AnalysisHandler>
  describeHerwigCrossSectionAnalysis("Herwig::CrossSectionAnalysis", "HwJetsAnalysis.so");

void CrossSectionAnalysis::Init() {

  static ClassDocumentation<CrossSectionAnalysis> documentation
    ("There is no documentation for the CrossSectionAnalysis class");

}

