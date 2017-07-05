// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanPT class.
//

#include "DrellYanPT.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DrellYanPT::DrellYanPT()
 : _Zpt(0.,250.,250), _Wppt(0.,250.,250), _Wmpt(0.,250.,250) {}

void DrellYanPT::dofinish() {
  AnalysisHandler::dofinish();
  ofstream outZ ("pt_Z.dat");
  _Zpt.normaliseToCrossSection();
  _Zpt.simpleOutput(outZ,true);
  ofstream outWm ("pt_Wm.dat");
  _Wmpt.normaliseToCrossSection();
  _Wmpt.simpleOutput(outWm,true);
  ofstream outWp ("pt_Wp.dat");
  _Wppt.normaliseToCrossSection();
  _Wppt.simpleOutput(outWp,true);
}

void DrellYanPT::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator send=event->primaryCollision()->steps().end();
  for(;sit!=send;++sit) {
    ParticleSet part=(**sit).all();
    ParticleSet::const_iterator iter=part.begin();
    ParticleSet::const_iterator end =part.end();
    for( ;iter!=end;++iter) {
      if(((**iter).id()==ParticleID::Z0||(**iter).id()==ParticleID::gamma)
	 && (**iter).children().size()==2) {
	_Zpt.addWeighted((**iter).momentum().perp()/GeV,event->weight());
      } else if ((**iter).id()==ParticleID::Wplus && (**iter).children().size()==2) {
	_Wppt.addWeighted((**iter).momentum().perp()/GeV,event->weight());
	
      } else if ((**iter).id()==ParticleID::Wminus && (**iter).children().size()==2) {
	_Wmpt.addWeighted((**iter).momentum().perp()/GeV,event->weight());
      }
    }
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<DrellYanPT,AnalysisHandler>
describeHerwigDrellYanPT("Herwig::DrellYanPT", "HwAnalysis.so");

void DrellYanPT::Init() {

  static ClassDocumentation<DrellYanPT> documentation
    ("Analyses the pt of weak bosons produces in Drell-Yan processes.");


}

