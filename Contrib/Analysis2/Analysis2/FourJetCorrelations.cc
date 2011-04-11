// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "FourJetCorrelations.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FourJetCorrelations.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

FourJetCorrelations::~FourJetCorrelations() {}

void FourJetCorrelations::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _cosAlpha34 << _cosChiBZ << _cosPhiKSW << _cosThetaNR;
}

void FourJetCorrelations::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _cosAlpha34 >> _cosChiBZ >> _cosPhiKSW >> _cosThetaNR;
}

ClassDescription<FourJetCorrelations> FourJetCorrelations::initFourJetCorrelations;
// Definition of the static class description member.

void FourJetCorrelations::Init() {

  static ClassDocumentation<FourJetCorrelations> documentation
    ("Four jet correlations at e+e- colliders (preferably LEP)");

  static Parameter<FourJetCorrelations,string> interfaceCosAlpha34
    ("CosAlpha34",
     "Options for CosAlpha34",
     &FourJetCorrelations::_cosAlpha34, "",
     false, false);

  static Parameter<FourJetCorrelations,string> interfaceCosChiBZ
    ("CosChiBZ",
     "Options for CosChiBZ",
     &FourJetCorrelations::_cosChiBZ, "",
     false, false);

  static Parameter<FourJetCorrelations,string> interfaceCosPhiKSW
    ("CosPhiKSW",
     "Options for CosPhiKSW",
     &FourJetCorrelations::_cosPhiKSW, "",
     false, false);

  static Parameter<FourJetCorrelations,string> interfaceCosThetaNR
    ("CosThetaNR",
     "Options for CosThetaNR",
     &FourJetCorrelations::_cosThetaNR, "",
     false, false);

}

void FourJetCorrelations::doinit() throw(InitException) {
  Analysis2Base::doinit();

  if (!jetFinder())
    Throw<InitException>() << "YMerge : No JetFinder has been set, giving up.";

  int plotFlags = HistogramOutput::Frame | HistogramOutput::Errorbars;
  Histogram2Options options (plotFlags);

  insert("CosAlpha34",_cosAlpha34,options);
  insert("CosChiBZ",_cosChiBZ,options);
  insert("CosPhiKSW",_cosPhiKSW,options);
  insert("CosThetaNR",_cosThetaNR,options);

}

void FourJetCorrelations::analyze(const tPVector &) {

  pair<vector<Lorentz5Momentum>,double> ev;

  while (*eventExtractor() >> ev) {

    jetFinder()->use(ev.first);
    jetFinder()->findJetsY();

    if (jetFinder()->getNJets() != 4) continue;

    jetFinder()->sortEnergy();

    book(cosAlpha34(jetFinder()->jets()),"CosAlpha34",ev.second);
    book(abs(cosChiBZ(jetFinder()->jets())),"CosChiBZ",ev.second);
    book(cosPhiKSW(jetFinder()->jets()),"CosPhiKSW",ev.second);
    book(abs(cosThetaNR(jetFinder()->jets())),"CosThetaNR",ev.second);

  }

}

void FourJetCorrelations::dofinish() {

  finish("CosAlpha34");
  finish("CosChiBZ");
  finish("CosPhiKSW");
  finish("CosThetaNR");

  Analysis2Base::dofinish();
}
