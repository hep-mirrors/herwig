// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZHTest class.
//

#include "ZHTest.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"


using namespace Herwig;

void ZHTest::analyze(tEventPtr event, long ieve, int loop, int state) {
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
	*_mH     += (**iter).momentum().m()/GeV;
	*_phiH   += (**iter).momentum().phi()+Constants::pi;
	*_yH += (**iter).momentum().rapidity();
	*_ptH += (**iter).momentum().perp()/GeV;
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
	  *_yZ += (**iter).momentum().rapidity();
	  *_ptZ += (**iter).momentum().perp()/GeV;
	  for(unsigned int ix=0;ix<(**iter).children().size();++ix) {
	    if((**iter).children()[ix]->id()==11) {
	      *_phil[0]   += (**iter).children()[ix]->momentum().phi()+Constants::pi;
	      *_yl[0]     += (**iter).children()[ix]->momentum().rapidity();
	      *_ptl[0]    += (**iter).children()[ix]->momentum().perp()/GeV;
	    }
	    else if((**iter).children()[ix]->id()==-11) {
	      *_phil[1]   += (**iter).children()[ix]->momentum().phi()+Constants::pi;
	      *_yl[1]     += (**iter).children()[ix]->momentum().rapidity();
	      *_ptl[1]    += (**iter).children()[ix]->momentum().perp()/GeV;
	    }
	  }
	}
      }
    }
  }
}

NoPIOClassDescription<ZHTest> ZHTest::initZHTest;
// Definition of the static class description member.

void ZHTest::Init() {

  static ClassDocumentation<ZHTest> documentation
    ("There is no documentation for the ZHTest class");

}

void ZHTest::doinitrun() {
  AnalysisHandler::doinitrun();
  if(getParticleData(ParticleID::h0)->mass()>200.*GeV) 
    _mH     = new_ptr(Histogram(200.,            400.,200));
  else
    _mH     = new_ptr(Histogram(125.,            127.0,200));
  _mZ     = new_ptr(Histogram(  0.0,            200.0,400));
  _yH     = new_ptr(Histogram( -10.0,            10.0,200));
  _yZ     = new_ptr(Histogram( -10.0,            10.0,200));
  _ptH    = new_ptr(Histogram( 0.,              7000.0,1000));
  _ptZ    = new_ptr(Histogram( 0.,              7000.0,1000));
  _phiH   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _phiZ   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  for(unsigned int ix=0;ix<2;++ix) {
    _ptl[ix]  =  new_ptr(Histogram( 0.,              7000.0,1000));
    _yl[ix]   = new_ptr(Histogram( -10.0,            10.0,200));
    _phil[ix] = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  }
}

void ZHTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  string title,species;
  title = "mass of H";
  _mH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pt of H";
  _ptH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of H";
  _yH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of H";
  _phiH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "mass of Z";
  _mZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pt of Z";
  _ptZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of Z";
  _yZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of Z";
  _phiZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pt of e-";
  _ptl[0]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of e-";
  _yl[0]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "phi of e-";
  _phil[0]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "pt of e+";
  _ptl[1]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of e+";
  _yl[1]->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "phi of e+";
  _phil[1]->topdrawOutput(outfile,Frame,"BLACK",title);
}
