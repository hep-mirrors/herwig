// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HJetsAnalysis class.
//

#include "HJetsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HJetsAnalysis::HJetsAnalysis() {}

HJetsAnalysis::~HJetsAnalysis() {}

IBPtr HJetsAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr HJetsAnalysis::fullclone() const {
  return new_ptr(*this);
}

void HJetsAnalysis::reconstructHardObjects(ParticleVector& parts) {
  ParticleVector::iterator higgs = parts.begin();
  for ( ; higgs != parts.end(); ++higgs ) {
    if ( (**higgs).id() == ParticleID::h0 )
      break;
  }
  if ( higgs == parts.end() )
    throw Exception() << "No Higgs found in HJetsAnalysis"
		      << Exception::abortnow;
  hardObjectMomentum("h") = (**higgs).momentum();
  parts.erase(higgs);
}

void HJetsAnalysis::analyzeSpecial(long id, double weight) {
  if ( nJets() < 2 )
    return;
  higgsYStar().count(Statistics::EventContribution(yStar(jetMomentum(1),jetMomentum(2),hardObjectMomentum("h")),weight,0.1),id);
  if ( nJets() > 2 )
    thirdJetYStar().count(Statistics::EventContribution(yStar(jetMomentum(1),jetMomentum(2),jetMomentum(3)),weight,0.1),id);
  if ( nJets() > 3 )
    fourthJetYStar().count(Statistics::EventContribution(yStar(jetMomentum(1),jetMomentum(2),jetMomentum(4)),weight,0.1),id);
  if ( nJets() > 1 ) {
    LorentzMomentum p12 = jetMomentum(1) + jetMomentum(2);
    double dphi = JetsPlusAnalysis::PairProperties::dPhi(p12,hardObjectMomentum("h"));
    jet12HiggsDeltaPhi().count(Statistics::EventContribution(dphi,weight,0.1),id);
    // Jeppe delta phi
    double minRap = jetMomentum(1).rapidity();
    double maxRap = minRap;
    for ( unsigned int k = 2; k <= nJets(); ++k ) {
      double yjet = jetMomentum(k).rapidity();
      minRap = min(minRap,yjet);
      maxRap = max(maxRap,yjet);
    }
    double hrap = hardObjectMomentum("h").rapidity();
    if ( minRap < hrap && hrap < maxRap ) {
      LorentzMomentum pLeft, pRight;
      for ( unsigned int k = 1; k <= nJets(); ++k ) {
	const LorentzMomentum& pjet = jetMomentum(k);
	if ( pjet.rapidity() < hrap )
	  pLeft += pjet;
	else
	  pRight += pjet;
      }
      double jeppedphi = JetsPlusAnalysis::PairProperties::dPhi(pLeft,pRight);
      jeppeDeltaPhi().count(Statistics::EventContribution(jeppedphi,weight,0.1),id);
    }
  }
}

void HJetsAnalysis::finalize(XML::Element& xhistos) {
  if ( !theHiggsYStar.bins().empty() ) {
    theHiggsYStar.finalize();
    xhistos.append(theHiggsYStar.toXML());
  }
  if ( !theThirdJetYStar.bins().empty() ) {
    theThirdJetYStar.finalize();
    xhistos.append(theThirdJetYStar.toXML());
  }
  if ( !theFourthJetYStar.bins().empty() ) {
    theFourthJetYStar.finalize();
    xhistos.append(theFourthJetYStar.toXML());
  }
  if ( !theJet12HiggsDeltaPhi.bins().empty() ) {
    theJet12HiggsDeltaPhi.finalize();
    xhistos.append(theJet12HiggsDeltaPhi.toXML());
  }
  if ( !theJeppeDeltaPhi.bins().empty() ) {
    theJeppeDeltaPhi.finalize();
    xhistos.append(theJeppeDeltaPhi.toXML());
  }
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HJetsAnalysis::persistentOutput(PersistentOStream &) const {}

void HJetsAnalysis::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HJetsAnalysis,Herwig::JetsPlusAnalysis>
  describeHerwigHJetsAnalysis("Herwig::HJetsAnalysis", "JetCuts.so HwJetsAnalysis.so");

void HJetsAnalysis::Init() {

  static ClassDocumentation<HJetsAnalysis> documentation
    ("There is no documentation for the HJetsAnalysis class");

}

