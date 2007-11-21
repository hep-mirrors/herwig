// -*- C++ -*-
//
// KtJetInterface.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "Herwig++/Interfaces/KtJetInterface.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace Herwig;
using namespace ThePEG;
using namespace KtJet;

vector<KtLorentzVector> KtJetInterface::convert(const tPVector &pv) {
  vector<KtLorentzVector> rval;
  for(tPVector::const_iterator it = pv.begin(); it != pv.end(); it++) {
    rval.push_back(convert(*it));
    Kt2PythiaMap[rval.back().getID()] = (*it)->number();
  }
  return rval;
}
   
KtLorentzVector KtJetInterface::convert(tcPPtr particle) {
  Lorentz5Momentum p = particle->momentum();
  return KtLorentzVector(p.x()/MeV, p.y()/MeV, p.z()/MeV, p.e()/MeV);
}

int KtJetInterface::getThePEGID(KtLorentzVector &kv) {
  return Kt2PythiaMap[kv.getID()];
}

LorentzMomentum KtJetInterface::convert(const KtLorentzVector & kt) {
  return LorentzMomentum(kt.x()*MeV, kt.y()*MeV, kt.z()*MeV, kt.e()*MeV);
}

vector<LorentzMomentum> KtJetInterface::convert(const vector<KtLorentzVector> & kt) {
  vector<LorentzMomentum> ret;
  ret.reserve(kt.size());
  for (vector<KtLorentzVector>::const_iterator it = kt.begin();
       it != kt.end(); ++it) {
    ret.push_back(convert(*it));
  }
  return ret;
}
