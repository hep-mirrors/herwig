// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZJetTest class.
//

#include "ZJetTest.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;


void ZJetTest::analyze(tEventPtr event, long ieve, int loop, int state) {
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
      if (abs((**iter).id())==ParticleID::Z0) {
	pz=(*iter)->momentum();
	double pt = pz.perp()/GeV;
	double y  = pz.rapidity();
	double m  = pz.m()/GeV;
	*_ptZ += pt;
	*_mZ  += m;
	*_yZ  += y;
	*_phiZ+= pz.phi()+Constants::pi;
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

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<ZJetTest,AnalysisHandler>
describeHerwigZJetTest("Herwig::ZJetTest", "HadronTest.so");

void ZJetTest::Init() {

  static ClassDocumentation<ZJetTest> documentation
    ("There is no documentation for the ZJetTest class");

}

void ZJetTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  string title,species;
  title = "pT of Z";
  _ptZ->topdrawOutput(outfile,Frame|Ylog,"BLACK",title);
  title = "mass of Z";
  _mZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "rapidity of Z";
  _yZ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of Z";
  _phiZ->topdrawOutput(outfile,Frame,"BLACK",title);
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

void ZJetTest::doinitrun() {
  AnalysisHandler::doinitrun();
  _ptZ = new_ptr(Histogram(  0.,7000.,700));
  _mZ  = new_ptr(Histogram( 60.,120. ,600));
  _yZ  = new_ptr(Histogram(-10.,10   ,200));
  _phiZ= new_ptr(Histogram(0.,2.*Constants::pi,200 ));
  for(unsigned int ix=0;ix<4;++ix) {
    _ptl  [ix] = new_ptr(Histogram(  0.,7000.,700));
    _yl   [ix] = new_ptr(Histogram(-10.,  10.,200 ));
    _phil [ix] = new_ptr(Histogram(0.,2.*Constants::pi,200 ));
  }
}
