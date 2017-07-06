// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QQHTest class.
//

#include "QQHTest.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

QQHTest::QQHTest() : quarkFlavour_(6) {}

void QQHTest::analyze(tEventPtr event, long ieve, int loop, int state) {
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
  Lorentz5Momentum ptotal;
  double y1(0),y2(0),y3(0);
  for(;sit!=send;++sit) {
    ParticleSet part;
    (**sit).selectFinalState(inserter(part));
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      if((**iter).id()==quarkFlavour_) {
	ptotal+=(**iter).momentum();
	QpT_ ->addWeighted((**iter).momentum().perp()/GeV,event->weight());
	Qrap_->addWeighted((**iter).momentum().rapidity(),event->weight());
	Qphi_->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
	y1 = (**iter).momentum().rapidity();
      }
      else if((**iter).id()==-quarkFlavour_) {
	ptotal+=(**iter).momentum();
	QBpT_ ->addWeighted((**iter).momentum().perp()/GeV,event->weight());
	QBrap_->addWeighted((**iter).momentum().rapidity(),event->weight());
	QBphi_->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
	y2 = (**iter).momentum().rapidity();
      }
      else if((**iter).id()==ParticleID::h0) {
	HpT_ ->addWeighted((**iter).momentum().perp()/GeV,event->weight());
	Hrap_->addWeighted((**iter).momentum().rapidity(),event->weight());
	Hphi_->addWeighted((**iter).momentum().phi()+Constants::pi,event->weight());
	y3 = (**iter).momentum().rapidity();
	ptotal+=(**iter).momentum();
      }
    }
  }
  mass_->addWeighted(ptotal.m()/GeV,event->weight());
  y12_->addWeighted(y1-y2,event->weight());
  y13_->addWeighted(y1-y3,event->weight());
  y23_->addWeighted(y2-y3,event->weight());
}

IBPtr QQHTest::clone() const {
  return new_ptr(*this);
}

IBPtr QQHTest::fullclone() const {
  return new_ptr(*this);
}

void QQHTest::persistentOutput(PersistentOStream & os) const {
  os << quarkFlavour_;
}

void QQHTest::persistentInput(PersistentIStream & is, int) {
  is >> quarkFlavour_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QQHTest,AnalysisHandler>
describeHerwigQQHTest("Herwig::QQHTest", "HadronTest.so");

void QQHTest::Init() {

  static ClassDocumentation<QQHTest> documentation
    ("There is no documentation for the QQHTest class");

  static Switch<QQHTest,int> interfaceQuarkFlavour
    ("QuarkFlavour",
     "The flavour of the heavy quark",
     &QQHTest::quarkFlavour_, 6, false, false);
  static SwitchOption interfaceQuarkFlavourBottom
    (interfaceQuarkFlavour,
     "Bottom",
     "bottom quarks",
     5);
  static SwitchOption interfaceQuarkFlavourTop
    (interfaceQuarkFlavour,
     "Top",
     "top quarks",
     6);
}

void QQHTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  QpT_->topdrawOutput(outfile,Frame,"BLACK","Pt of Q");
  Qrap_->topdrawOutput(outfile,Frame,"BLACK","Rapidity of Q");
  Qphi_->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for Q");
  QBpT_->topdrawOutput(outfile,Frame,"BLACK","Pt of Qbar");
  QBrap_->topdrawOutput(outfile,Frame,"BLACK","Rapidity of Qbar");
  QBphi_->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for Qbar");
  HpT_ ->topdrawOutput(outfile,Frame,"BLACK","Pt of H ");
  Hrap_ ->topdrawOutput(outfile,Frame,"BLACK","Rapidity of H ");
  Hphi_ ->topdrawOutput(outfile,Frame,"BLACK","Azimuthal angle for H ");
  mass_ ->topdrawOutput(outfile,Frame,"BLACK","system mass");
  QpT_->normaliseToCrossSection();
  QpT_->topdrawOutput(outfile,Frame,"BLACK","Pt of Q");
  Qrap_->normaliseToCrossSection();
  Qrap_->topdrawOutput(outfile,Frame,"BLACK","Rapidity of Q");
  QBpT_->normaliseToCrossSection();
  QBpT_->topdrawOutput(outfile,Frame,"BLACK","Pt of Qbar");
  QBrap_->normaliseToCrossSection();
  QBrap_->topdrawOutput(outfile,Frame,"BLACK","Rapidity of Qbar");
  HpT_ ->normaliseToCrossSection();
  HpT_ ->topdrawOutput(outfile,Frame,"BLACK","Pt of H ");
  Hrap_ ->normaliseToCrossSection();
  Hrap_ ->topdrawOutput(outfile,Frame,"BLACK","Rapidity of H ");
  mass_ ->normaliseToCrossSection();
  mass_ ->topdrawOutput(outfile,Frame,"BLACK","system mass");
  y12_->topdrawOutput(outfile,Frame,"BLACK","yQQbar");
  y13_->topdrawOutput(outfile,Frame,"BLACK","yQH");
  y23_->topdrawOutput(outfile,Frame,"BLACK","yQbarH");
  y12_->normaliseToCrossSection();
  y13_->normaliseToCrossSection();
  y23_->normaliseToCrossSection();
  y12_->topdrawOutput(outfile,Frame,"BLACK","yQQbar");
  y13_->topdrawOutput(outfile,Frame,"BLACK","yQH");
  y23_->topdrawOutput(outfile,Frame,"BLACK","yQbarH");
}


void QQHTest::doinitrun() {
  AnalysisHandler::doinitrun();
  QpT_   = new_ptr(Histogram(0.,1500.,300));
  Qrap_  = new_ptr(Histogram(-15,15,300));
  Qphi_  = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  QBpT_  = new_ptr(Histogram(0.,1500.,300));
  QBrap_ = new_ptr(Histogram(-15,15,300));
  QBphi_ = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  HpT_   = new_ptr(Histogram(0.,1500.,300));
  Hrap_  = new_ptr(Histogram(-15,15,300));
  Hphi_  = new_ptr(Histogram(  0.0,2.0*Constants::pi,200));
  mass_  = new_ptr(Histogram(0.,3000.,300));
  y12_ = new_ptr(Histogram(-15,15,300));
  y13_ = new_ptr(Histogram(-15,15,300));
  y23_ = new_ptr(Histogram(-15,15,300));
}
