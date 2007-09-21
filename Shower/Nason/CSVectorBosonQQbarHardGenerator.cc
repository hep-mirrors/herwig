// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CSVectorBosonQQbarHardGenerator class.
//

#include "CSVectorBosonQQbarHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace Herwig;

void CSVectorBosonQQbarHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << ounit(_ptmin,GeV);
}

void CSVectorBosonQQbarHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> iunit(_ptmin,GeV);
}

ClassDescription<CSVectorBosonQQbarHardGenerator> 
CSVectorBosonQQbarHardGenerator::initCSVectorBosonQQbarHardGenerator;
// Definition of the static class description member.

void CSVectorBosonQQbarHardGenerator::Init() {

  static ClassDocumentation<CSVectorBosonQQbarHardGenerator> documentation
    ("There is no documentation for the CSVectorBosonQQbarHardGenerator class");

  static Reference<CSVectorBosonQQbarHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &CSVectorBosonQQbarHardGenerator::_alphaS, false, false, true, false, false);

  static Parameter<CSVectorBosonQQbarHardGenerator,Energy> interfacepTMin
    ("pTMin",
     "The minimum pt of the hardest emission",
     &CSVectorBosonQQbarHardGenerator::_ptmin, GeV, 0.5*GeV, 0.2*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

bool CSVectorBosonQQbarHardGenerator::canHandle(ShowerTreePtr tree) {
  if(tree->incomingLines().size()!=1) return false;    
  if((tree->incomingLines().begin()->first->id()==22)&&
  (tree->incomingLines().begin()->first->progenitor()->id()==23)) return false;
  map<ShowerProgenitorPtr,tShowerParticlePtr> outgoing=tree->outgoingLines();
  if(outgoing.size()!=2) return false;
  if(abs(outgoing.begin()->first->progenitor()->id())>6)  return false;
  if(outgoing.begin()->first->progenitor()->id()!=
     -1*outgoing.rbegin()->first->progenitor()->id())     return false;
  return true;
}

NasonTreePtr CSVectorBosonQQbarHardGenerator::generateHardest(ShowerTreePtr tree) {
  // Get the progenitors: Q and Qbar.
  vector<tcPDPtr> partons(2);
  ShowerProgenitorPtr 
    QProgenitor   =tree->outgoingLines().begin()->first,
    QbarProgenitor=tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap(QProgenitor   ,QbarProgenitor);
  partons[0]=QProgenitor->progenitor()->dataPtr();
  partons[1]=QbarProgenitor->progenitor()->dataPtr();
}

void CSVectorBosonQQbarHardGenerator::generate(Energy2 s) {
  double c = 2.*_alphaS->overestimateValue()/3./Constants::pi;
  Energy2 pt2(s);
  do {
    pt2 = pt2*exp(sqrt(log(UseRandom::rnd())/c));
  }
  while(pt2>sqr(_ptmin));
}
