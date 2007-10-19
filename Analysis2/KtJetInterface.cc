// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KtJetInterface class.
//

#include "KtJetInterface.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KtJetInterface.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KtJetInterface::~KtJetInterface() {}

void KtJetInterface::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _collisionType << _distanceScheme << _recombinationScheme;
}

void KtJetInterface::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _collisionType >> _distanceScheme >> _recombinationScheme;
}

ClassDescription<KtJetInterface> KtJetInterface::initKtJetInterface;
// Definition of the static class description member.

void KtJetInterface::Init() {

  static ClassDocumentation<KtJetInterface> documentation
    ("Interface to KtJet");


  static Switch<KtJetInterface,int> interfaceCollisionType
    ("CollisionType",
     "Set  the collision type used by KtJet",
     &KtJetInterface::_collisionType, 1, false, false);
  static SwitchOption interfaceCollisionTypeee
    (interfaceCollisionType,
     "ee",
     "ee",
     1);
  static SwitchOption interfaceCollisionTypeeP
    (interfaceCollisionType,
     "eP",
     "eP",
     2);
  static SwitchOption interfaceCollisionTypePe
    (interfaceCollisionType,
     "Pe",
     "Pe",
     3);
  static SwitchOption interfaceCollisionTypePP
    (interfaceCollisionType,
     "PP",
     "PP",
     4);


  static Switch<KtJetInterface,int> interfaceDistanceScheme
    ("DistanceScheme",
     "Set the distance scheme (\"angle\") used by KtJet",
     &KtJetInterface::_distanceScheme, 1, false, false);
  static SwitchOption interfaceDistanceSchemeangular
    (interfaceDistanceScheme,
     "angular",
     "angular",
     1);
  static SwitchOption interfaceDistanceSchemedeltaR
    (interfaceDistanceScheme,
     "deltaR",
     "deltaR",
     2);
  static SwitchOption interfaceDistanceSchemeQCD
    (interfaceDistanceScheme,
     "QCD",
     "QCD",
     3);



  static Switch<KtJetInterface,int> interfaceRecombinationScheme
    ("RecombinationScheme",
     "Set the recombination scheme used by KtJet",
     &KtJetInterface::_recombinationScheme, 1, false, false);
  static SwitchOption interfaceRecombinationSchemeE
    (interfaceRecombinationScheme,
     "E",
     "E",
     1);
  static SwitchOption interfaceRecombinationSchemePt
    (interfaceRecombinationScheme,
     "Pt",
     "Pt",
     2);
  static SwitchOption interfaceRecombinationSchemePt2
    (interfaceRecombinationScheme,
     "Pt2",
     "Pt2",
     3);
  static SwitchOption interfaceRecombinationSchemeEt
    (interfaceRecombinationScheme,
     "Et",
     "Et",
     4);
  static SwitchOption interfaceRecombinationSchemeEt2
    (interfaceRecombinationScheme,
     "Et2",
     "Et2",
     5);

}

void KtJetInterface::use (tcEventPtr evt, bool inclusive) {
  JetFinder::use(evt, inclusive);
  convert();
  if (inclusive)
    clusterInclusive();
  else
    cluster();
}

void KtJetInterface::convert () {
  _lastMomenta.clear();
  tPVector final = lastEvent()->getFinalState();
  for(tPVector::iterator p = final.begin();
      p != final.end(); ++p) {
    _lastMomenta.push_back(KtJet::KtLorentzVector((**p).momentum().x()/GeV,
						  (**p).momentum().y()/GeV,
						  (**p).momentum().z()/GeV,
						  (**p).momentum().t()/GeV
						  ));
  }
}

void KtJetInterface::convert (const vector<KtJet::KtLorentzVector>& ktvectors) {
  list<Lorentz5Momentum> tmp;
  for(vector<KtJet::KtLorentzVector>::const_iterator p =ktvectors.begin(); p != ktvectors.end(); ++p) {
    tmp.push_back(Lorentz5Momentum((*p).px()*GeV,(*p).py()*GeV,(*p).pz()*GeV,(*p).e()*GeV));
  }
  jets(tmp);
}

void KtJetInterface::cluster () {
  _lastKtEvent =
    std::auto_ptr<KtJet::KtEvent>
    (new KtJet::KtEvent(_lastMomenta,_collisionType,_distanceScheme,_recombinationScheme));
}

void KtJetInterface::clusterInclusive () {
  _lastKtEvent =
    std::auto_ptr<KtJet::KtEvent>
    (new KtJet::KtEvent(_lastMomenta,_collisionType,_distanceScheme,_recombinationScheme));
}

