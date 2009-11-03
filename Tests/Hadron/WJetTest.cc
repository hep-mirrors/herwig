// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WJetTest class.
//

#include "WJetTest.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;


void WJetTest::analyze(tEventPtr event, long ieve, int loop, int state) {
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
      if (abs((**iter).id())==ParticleID::Wplus) {
	pz=(*iter)->momentum();
	double pt = pz.perp()/GeV;
	double y  = pz.rapidity();
	double phi=pz.phi()+Constants::pi;
	double m  = pz.m()/GeV;
	*_ptW  [0] += pt;
	*_mW   [0] += m;
	*_yW   [0] += y;
	*_phiW [0] += phi;
	if((**iter).id()>0) {
	  *_ptW  [1] += pt;
	  *_mW   [1] += m;
	  *_yW   [1] += y;
	  *_phiW [1] += phi;
	}
	else {
	  *_ptW  [2] += pt;
	  *_mW   [2] += m;
	  *_yW   [2] += y;
	  *_phiW [2] += phi;
	}
      }
      else if((**iter).id()==ParticleID::eminus) {
	*_ptl[0] += (*iter)->momentum().perp()/GeV;
	*_yl [0] += (*iter)->momentum().rapidity();
	*_phil[0] += (*iter)->momentum().phi()+Constants::pi;
      }
      else if((**iter).id()==ParticleID::eplus) {
	*_ptl[1] += (*iter)->momentum().perp()/GeV;
	*_yl [1] += (*iter)->momentum().rapidity();
	*_phil[1] += (*iter)->momentum().phi()+Constants::pi;
      }
      else if((**iter).id()==ParticleID::nu_e) {
	*_ptl[2] += (*iter)->momentum().perp()/GeV;
	*_yl [2] += (*iter)->momentum().rapidity();
	*_phil[2] += (*iter)->momentum().phi()+Constants::pi;
      }
      else if((**iter).id()==ParticleID::nu_ebar) {
	*_ptl[3] += (*iter)->momentum().perp()/GeV;
	*_yl [3] += (*iter)->momentum().rapidity();
	*_phil[3] += (*iter)->momentum().phi()+Constants::pi;
      }
      
    }
  }
}

NoPIOClassDescription<WJetTest> WJetTest::initWJetTest;
// Definition of the static class description member.

void WJetTest::Init() {

  static ClassDocumentation<WJetTest> documentation
    ("There is no documentation for the WJetTest class");

}

void WJetTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  string title,species;
  for(unsigned int ix=0;ix<3;++ix) {
    if     (ix==0) species = "all W";
    else if(ix==1) species = "W+";
    else if(ix==2) species = "W-";
    title = "pT of " +species;
    _ptW[ix]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
    title = "mass of " + species;
    _mW[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "rapidity of " + species;
    _yW[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "azimuth of " + species;
    _phiW[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
  }
  for(unsigned int ix=0;ix<4;++ix) {
    if     (ix==0) species = "e-";
    else if(ix==1) species = "e+";
    else if(ix==2) species = "nu_e";
    else if(ix==3) species = "nu_ebar";
    title = "pT of " +species;
    _ptl[ix]->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
    title = "rapidity of " + species;
    _yl[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
    title = "azimuth of " + species;
    _phil[ix]->topdrawOutput(outfile,Frame,"BLACK",title);
  }
}

void WJetTest::doinitrun() {
  AnalysisHandler::doinitrun();
  for(unsigned int ix=0;ix<3;++ix) {
    _ptW  [ix] = new_ptr(Histogram(  0.,7000.,700));
    _mW   [ix] = new_ptr(Histogram( 50.,110. ,600));
    _yW   [ix] = new_ptr(Histogram(-10.,10   ,200));
    _phiW [ix] = new_ptr(Histogram(0.,2.*Constants::pi,200 ));
  }
  for(unsigned int ix=0;ix<4;++ix) {
    _ptl  [ix] = new_ptr(Histogram(  0.,7000.,700));
    _yl   [ix] = new_ptr(Histogram(-10.,  10.,200 ));
    _phil [ix] = new_ptr(Histogram(0.,2.*Constants::pi,200 ));
  }
}
