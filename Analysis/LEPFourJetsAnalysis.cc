// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LepFourJetsAnalysis class.
//

#include "LEPFourJetsAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LEPFourJetsAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LEPFourJetsAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
}

LorentzRotation LEPFourJetsAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LEPFourJetsAnalysis::analyze(const tPVector & particles) {
  _kint.clearMap();
  KtJet::KtEvent ev = KtJet::KtEvent(_kint.convert(particles), 1, 1, 1);
  // four jet distributions
  ev.findJetsY(0.008); 
  vector<KtJet::KtLorentzVector> ktjets = ev.getJetsE();
  vector<Lorentz5Momentum> jets;

  if (ktjets.size() == 4) {
    for (int j=0; j<4; ++j) {
      if (! ktjets[j].isJet()) {
	throw Exception() << "LEPFourJetsAnalysis: Trying to extract jet " 
			  << "momenta from a single particle." 
			  << Exception::warning;
      }
      jets.push_back(KtJetInterface::convert(ktjets[j]));
    }
    *_cchiBZ += abs(cosChiBZ(jets));
    *_cphiKSW += cosPhiKSW(jets);
    *_cthNR += abs(cosThetaNR(jets));
    *_ca34 += cosAlpha34(jets);
  }
}

void LEPFourJetsAnalysis::analyze(tPPtr) {}

void LEPFourJetsAnalysis::persistentOutput(PersistentOStream & ) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void LEPFourJetsAnalysis::persistentInput(PersistentIStream & , int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<LEPFourJetsAnalysis> LEPFourJetsAnalysis::initLEPFourJetsAnalysis;
// Definition of the static class description member.

void LEPFourJetsAnalysis::Init() {

  static ClassDocumentation<LEPFourJetsAnalysis> documentation
    ("There is no documentation for the LEPFourJetsAnalysis class");

}

