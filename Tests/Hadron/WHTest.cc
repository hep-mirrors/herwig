// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WHTest class.
//

#include "WHTest.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"


using namespace Herwig;

void WHTest::analyze(tEventPtr event, long ieve, int loop, int state) {
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
      else if(abs((**iter).id())==ParticleID::Wplus) {
	bool fermion=true;
	for(unsigned int ix=0;ix<(**iter).children().size();++ix) {
	  if(abs((**iter).children()[ix]->id())>16) {
	    fermion=false;
	    break;
	  }
	}
	if(fermion) {
	  *_mW[0]     += (**iter).momentum().m()/GeV;
	  *_phiW[0]   += (**iter).momentum().phi()+Constants::pi;
	  *_yW[0]     += (**iter).momentum().rapidity();
	  *_ptW[0]    += (**iter).momentum().perp()/GeV;
	  if((**iter).id()>0) {
	    *_mW[1]     += (**iter).momentum().m()/GeV;
	    *_phiW[1]   += (**iter).momentum().phi()+Constants::pi;
	    *_yW[1]     += (**iter).momentum().rapidity();
	    *_ptW[1]    += (**iter).momentum().perp()/GeV;
	  }
	  else {
	    *_mW[2]     += (**iter).momentum().m()/GeV;
	    *_phiW[2]   += (**iter).momentum().phi()+Constants::pi;
	    *_yW[2]     += (**iter).momentum().rapidity();
	    *_ptW[2]    += (**iter).momentum().perp()/GeV;
	  }
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
	    else if((**iter).children()[ix]->id()==12) {
	      *_phil[2]   += (**iter).children()[ix]->momentum().phi()+Constants::pi;
	      *_yl[2]     += (**iter).children()[ix]->momentum().rapidity();
	      *_ptl[2]    += (**iter).children()[ix]->momentum().perp()/GeV;
	    }
	    else if((**iter).children()[ix]->id()==-12) {
	      *_phil[3]   += (**iter).children()[ix]->momentum().phi()+Constants::pi;
	      *_yl[3]     += (**iter).children()[ix]->momentum().rapidity();
	      *_ptl[3]    += (**iter).children()[ix]->momentum().perp()/GeV;
	    }
	  }
	}
      }
    }
  }
}

NoPIOClassDescription<WHTest> WHTest::initWHTest;
// Definition of the static class description member.

void WHTest::Init() {

  static ClassDocumentation<WHTest> documentation
    ("There is no documentation for the WHTest class");

}

void WHTest::doinitrun() {
  AnalysisHandler::doinitrun();
  if(getParticleData(ParticleID::h0)->mass()>200.*GeV) 
    _mH     = new_ptr(Histogram(200.,            400.,200));
  else
    _mH     = new_ptr(Histogram(125.,            127.0,200));
  _yH     = new_ptr(Histogram( -10.0,            10.0,200));
  _ptH    = new_ptr(Histogram( 0.,              7000.0,1000));
  _phiH   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  for(unsigned int ix=0;ix<3;++ix) {
    _mW[ix]     = new_ptr(Histogram(  0.0,            200.0,400));
    _yW[ix]     = new_ptr(Histogram( -10.0,            10.0,200));
    _ptW[ix]    = new_ptr(Histogram( 0.,              7000.0,1000));
    _phiW[ix]   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  }
  for(unsigned int ix=0;ix<4;++ix) {
    _ptl[ix]  =  new_ptr(Histogram( 0.,              7000.0,1000));
    _yl[ix]   = new_ptr(Histogram( -10.0,            10.0,200));
    _phil[ix] = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  }
}

void WHTest::dofinish() {
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
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==0) species = "W all";
    else if(ix==1) species = "W+";
    else if(ix==2) species = "W-";
    title = "mass of "+species;
    _mW[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "pt of "+species;
    _ptW[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "rapidity of "+species;
    _yW[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "azimuth of "+species;
    _phiW[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
  }
  for(unsigned int ix=0;ix<4;++ix) {
    if(ix==0) species="e-";
    else if(ix==1) species="e+";
    else if(ix==2) species="nu_e";
    else if(ix==3) species="nu_ebar";
    title = "pt of "+species;
    _ptl[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "rapidity of "+species;
    _yl[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "phi of "+species;
    _phil[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
  }
}
