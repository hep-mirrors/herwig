// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeightAnalyzer class.
//

#include "WeightAnalyzer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

WeightAnalyzer::WeightAnalyzer() 
  : sumWeights(0.0), sumPositiveWeights(0.0),
    sumNegativeWeights(0.0),
    sumGroupWeights(0.0), sumPositiveGroupWeights(0.0),
    sumNegativeGroupWeights(0.0),
    maxDeviationGroupWeight(0.0),
    maxDeviationEventWeight(0.0) {}

WeightAnalyzer::~WeightAnalyzer() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

void WeightAnalyzer::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);

  sumWeights += event->weight();
  if ( event->weight() > 0.0 )
    sumPositiveWeights += event->weight();
  if ( event->weight() < 0.0 )
    sumNegativeWeights += event->weight();

  tSubProPtr sub = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr grp = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);

  double sumEvents = 0.0;
  double sumGroups = 0.0;

  sumGroupWeights += event->weight()*sub->groupWeight();
  if ( event->weight()*sub->groupWeight() > 0.0 )
    sumPositiveGroupWeights += event->weight()*sub->groupWeight();
  if ( event->weight()*sub->groupWeight() < 0.0 )
    sumNegativeGroupWeights += event->weight()*sub->groupWeight();

  sumEvents += event->weight()*sub->groupWeight();
  sumGroups += sub->groupWeight();

  if ( grp ) {
    
    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {

      sumGroupWeights += event->weight()*(**s).groupWeight();

      if ( event->weight()*(**s).groupWeight() > 0.0 )
	sumPositiveGroupWeights += event->weight()*(**s).groupWeight();
      if ( event->weight()*(**s).groupWeight() < 0.0 )
	sumNegativeGroupWeights += event->weight()*(**s).groupWeight();

      sumEvents += event->weight()*(**s).groupWeight();
      sumGroups += (**s).groupWeight();

    }

  }

  maxDeviationGroupWeight = max(maxDeviationGroupWeight,abs(sumGroups-1));
  maxDeviationEventWeight = max(maxDeviationEventWeight,abs(sumEvents/event->weight()-1));

}

void WeightAnalyzer::dofinish() {
  AnalysisHandler::dofinish();

  string dataName = generator()->filename() + "-weights.dat";

  ofstream data(dataName.c_str());

  data << setprecision(20)
       << "--------------------------------------------------------------------------------\n"
       << "WeightAnalyzer information\n"
       << "--------------------------------------------------------------------------------\n"
       << "sum of weights                        : " << sumWeights << "\n"
       << "sum of positive weights               : " << sumPositiveWeights << "\n"
       << "sum of negative weights               : " << sumNegativeWeights << "\n" 
       << "sum of weights (from groups)          : " << sumGroupWeights << "\n"
       << "sum of positive weights (from groups) : " << sumPositiveGroupWeights << "\n"
       << "sum of negative weights (from groups) : " << sumNegativeGroupWeights << "\n"
       << "maximum devation group weights        : " << maxDeviationGroupWeight << "\n"
       << "maximum devation event weights        : " << maxDeviationEventWeight << "\n"
       << "--------------------------------------------------------------------------------\n"
       << flush;

}

void WeightAnalyzer::doinitrun() {
  AnalysisHandler::doinitrun();
}


IBPtr WeightAnalyzer::clone() const {
  return new_ptr(*this);
}

IBPtr WeightAnalyzer::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void WeightAnalyzer::persistentOutput(PersistentOStream &) const {}

void WeightAnalyzer::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<WeightAnalyzer,AnalysisHandler>
  describeHerwigWeightAnalyzer("Herwig::WeightAnalyzer", "HwMatchbox.so");

void WeightAnalyzer::Init() {

  static ClassDocumentation<WeightAnalyzer> documentation
    ("There is no documentation for the WeightAnalyzer class");

}

