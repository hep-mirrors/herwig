// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFTest class.
//

#include "VBFTest.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

void VBFTest::analyze(tEventPtr event, long ieve, int loop, int state) {
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
    ParticleSet part;
    (**sit).selectFinalState(inserter(part));
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      if((**iter).id()==ParticleID::h0) {
	*_mH   += (**iter).momentum().m()/GeV;
	*_cosH += (**iter).momentum().cosTheta();
	*_phiH += (**iter).momentum().phi()+Constants::pi;
	*_eH   += (**iter).momentum().t()/GeV;
      }
      else if((**iter).id()==ParticleID::nu_e) {
	*_cosnu += (**iter).momentum().cosTheta();
	*_phinu += (**iter).momentum().phi()+Constants::pi;
	*_enu   += (**iter).momentum().t()/GeV;
      }
      else if((**iter).id()==ParticleID::nu_ebar) {
	*_cosnub += (**iter).momentum().cosTheta();
	*_phinub += (**iter).momentum().phi()+Constants::pi;
	*_enub   += (**iter).momentum().t()/GeV;
      }
      else if((**iter).id()==ParticleID::eminus) {
	*_cosem += (**iter).momentum().cosTheta();
	*_phiem += (**iter).momentum().phi()+Constants::pi;
	*_eem   += (**iter).momentum().t()/GeV;
      }
      else if((**iter).id()==ParticleID::eplus) {
	*_cosep  += (**iter).momentum().cosTheta();
	*_phiep  += (**iter).momentum().phi()+Constants::pi;
	*_eep    += (**iter).momentum().t()/GeV;
      }
    }
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<VBFTest,AnalysisHandler>
describeHerwigVBFTest("Herwig::VBFTest", "LeptonTest.so");

void VBFTest::Init() {

  static ClassDocumentation<VBFTest> documentation
    ("There is no documentation for the VBFTest class");

}

void VBFTest::doinitrun() {
  AnalysisHandler::doinitrun();
  if(getParticleData(ParticleID::h0)->mass()>200.*GeV) 
    _mH     = new_ptr(Histogram(200.,            400.,200));
  else
    _mH     = new_ptr(Histogram(114.,            116.0,200));
  _cosH     = new_ptr(Histogram( -1.0,              1.0,200));
  _phiH     = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _eH       = new_ptr(Histogram(  0.0,1000.,1000));
  _cosnu    = new_ptr(Histogram( -1.0,              1.0,200));
  _phinu    = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _enu      = new_ptr(Histogram(  0.0,1000.,1000));
  _cosnub   = new_ptr(Histogram( -1.0,              1.0,200));
  _phinub   = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _enub     = new_ptr(Histogram(  0.0,1000.,1000));
  _cosem    = new_ptr(Histogram( -1.0,              1.0,200));
  _phiem    = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _eem      = new_ptr(Histogram(  0.0,1000.,1000));
  _cosep    = new_ptr(Histogram( -1.0,              1.0,200));
  _phiep    = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  _eep      = new_ptr(Histogram(  0.0,1000.,1000));
}

void VBFTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  string title,species;
  title = "mass of H";
  _mH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of H";
  _cosH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "Energy of H";
  _eH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of H";
  _phiH->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of nu";
  _cosnu->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "Energy of nu";
  _enu->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of nu";
  _phinu->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of nub";
  _cosnub->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "Energy of nub";
  _enub->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of nub";
  _phinub->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of em";
  _cosem->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "Energy of em";
  _eem->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of em";
  _phiem->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "theta of ep ";
  _cosep ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "Energy of ep ";
  _eep ->topdrawOutput(outfile,Frame,"BLACK",title);
  title = "azimuth of ep ";
  _phiep ->topdrawOutput(outfile,Frame,"BLACK",title);
}
