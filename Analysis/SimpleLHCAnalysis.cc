// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleLHCAnalysis class.
//

#include "SimpleLHCAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SimpleLHCAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SimpleLHCAnalysis::~SimpleLHCAnalysis() {}

namespace {
  inline Lorentz5Momentum getMomentum(tcPPtr particle) {
    //return particle->momentum();
    Lorentz5Momentum tmp = particle->children()[0]->next()->momentum();
    tmp += particle->children()[1]->next()->momentum();
    tmp.rescaleMass();
    return tmp;

  }
}


void SimpleLHCAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum pz;
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator send=event->primaryCollision()->steps().end();
  for(;sit!=send;++sit)
    {
      ParticleSet part=(**sit).all();
      ParticleSet::const_iterator iter=part.begin();
      ParticleSet::const_iterator end =part.end();
      for( ;iter!=end;++iter)
	{
	  if((**iter).id()==ParticleID::Z0||
	     ((**iter).id()==ParticleID::gamma&&!(**iter).children().empty()))
	    {
	      pz=getMomentum(*iter);
	      Energy pt = pz.perp()/GeV;
	      Energy mz = pz.mass()/GeV;
	      if(mz>20.&&mz<80.) (_ptZ[1])+=(pt);
	      else if (mz>80.&&mz<100.) (_ptZ[2])+=(pt);
	      else if (mz>100.) (_ptZ[3])+=(pt);
	      (_ptZ[0])+=(pt);
	      (_mZ)+=(mz);
	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
	      (_rapZ)+=(rap);
	      (_phiZ)+=pz.phi();
	    } else if ((**iter).id()==ParticleID::Wplus) {
	      pz=getMomentum(*iter);
	      Energy pt = pz.perp()/GeV;
	      Energy mz = pz.mass()/GeV;
	      if(mz>20.&&mz<80.) (_ptWp[1])+=(pt);
	      else if (mz>80.&&mz<100.) (_ptWp[2])+=(pt);
	      else if (mz>100.) (_ptWp[3])+=(pt);
	      (_ptWp[0])+=(pt);
	      (_mWp)+=(mz);
	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
	      (_rapWp)+=(rap);
	      (_phiWp)+=pz.phi();
	    } else if ((**iter).id()==ParticleID::Wminus) {
	      pz=getMomentum(*iter);
	      Energy pt = pz.perp()/GeV;
	      Energy mz = pz.mass()/GeV;
	      if(mz>20.&&mz<80.) (_ptWm[1])+=(pt);
	      else if (mz>80.&&mz<100.) (_ptWm[2])+=(pt);
	      else if (mz>100.) (_ptWm[3])+=(pt);
	      (_ptWm[0])+=(pt);
	      (_mWm)+=(mz);
	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
	      (_rapWm)+=(rap);
	      (_phiWm)+=pz.phi();
	    }
	}
    }
}

LorentzRotation SimpleLHCAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void SimpleLHCAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
}

void SimpleLHCAnalysis::analyze(tPPtr) {}

void SimpleLHCAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void SimpleLHCAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<SimpleLHCAnalysis> SimpleLHCAnalysis::initSimpleLHCAnalysis;
// Definition of the static class description member.

void SimpleLHCAnalysis::Init() {

  static ClassDocumentation<SimpleLHCAnalysis> documentation
    ("There is no documentation for the SimpleLHCAnalysis class");

}

