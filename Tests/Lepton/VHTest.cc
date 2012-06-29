// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VHTest class.
//

#include "VHTest.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void VHTest::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator stest =event->primaryCollision()->steps().end();
  StepVector::const_iterator send=sit;
  ++send;
  if(send==stest) --send;
  ++send;
  if(send==stest) --send;
  ++send;
  Lorentz5Momentum pz;
  for(;sit!=send;++sit) {
    ParticleSet part=(**sit).all();
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      if((**iter).id()==ParticleID::h0) {
	if((**iter).momentum().m()>390.*GeV) cerr << (**iter) << "\n" << *event << "\n";
	*_mH     += (**iter).momentum().m()/GeV;
	*_phiH   += (**iter).momentum().phi()+Constants::pi;
	*_thetaH += (**iter).momentum().cosTheta();
      }
      else if((**iter).id()==ParticleID::Z0) {
	bool fermion=true;
	for(unsigned int ix=0;ix<(**iter).children().size();++ix) {
	  if(abs((**iter).children()[ix]->id())>16) {
	    fermion=false;
	    break;
	  }
	}
	if(fermion) {
	  *_mZ     += (**iter).momentum().m()/GeV;
	  *_phiZ   += (**iter).momentum().phi()+Constants::pi;
	  *_thetaZ += (**iter).momentum().cosTheta();	  
	  for(unsigned int ix=0;ix<(**iter).children().size();++ix) {
	    if((**iter).children()[ix]->id()==11) {
	      *_phil[0]   += (**iter).children()[ix]->momentum().phi()+Constants::pi;
	      *_thetal[0] += (**iter).children()[ix]->momentum().cosTheta();
	    }
	    else if((**iter).children()[ix]->id()==-11) {
	      *_phil[1]   += (**iter).children()[ix]->momentum().phi()+Constants::pi;
	      *_thetal[1] += (**iter).children()[ix]->momentum().cosTheta();
	    }
	  }
	}
      }
    }
  }
}

NoPIOClassDescription<VHTest> VHTest::initVHTest;
// Definition of the static class description member.

void VHTest::Init() {

  static ClassDocumentation<VHTest> documentation
    ("There is no documentation for the VHTest class");

}

void VHTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  string title,species;
  title = "mass of H";
  _mH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of H";
  _thetaH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of H";
  _phiH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "mass of Z";
  _mZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of Z";
  _thetaZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of Z";
  _phiZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of e-";
  _thetal[0]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "phi of e-";
  _phil[0]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of e+";
  _thetal[1]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "phi of e+";
  _phil[1]->topdrawOutput(outfile,Frame,"BLACK",title);
}

void VHTest::doinitrun() {
  AnalysisHandler::doinitrun();
  if(getParticleData(ParticleID::h0)->mass()>200.*GeV) 
    _mH     = new_ptr(Histogram(200.,            400.,200));
  else
    _mH     = new_ptr(Histogram(114.,            116.0,200));
  _mZ     = new_ptr(Histogram(  0.0,            200.0,400));
  _thetaH = new_ptr(Histogram( -1.0,              1.0,200));
  _thetaZ = new_ptr(Histogram( -1.0,              1.0,200));
  _phiH   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _phiZ   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  for(unsigned int ix=0;ix<2;++ix) {
    _thetal[ix] = new_ptr(Histogram( -1.0,              1.0,200));
    _phil[ix]   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  }
}
