// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "KtJetFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KtJetFinder.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

KtJetFinder::~KtJetFinder() {}

void KtJetFinder::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _collisionType << _distanceScheme << _recombinationScheme;
}

void KtJetFinder::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _collisionType >> _distanceScheme >> _recombinationScheme;
}

ClassDescription<KtJetFinder> KtJetFinder::initKtJetFinder;
// Definition of the static class description member.

void KtJetFinder::Init() {

  static ClassDocumentation<KtJetFinder> documentation
    ("Interface to KtJet");


  static Switch<KtJetFinder,int> interfaceCollisionType
    ("CollisionType",
     "Set  the collision type used by KtJet",
     &KtJetFinder::_collisionType, 1, false, false);
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


  static Switch<KtJetFinder,int> interfaceDistanceScheme
    ("DistanceScheme",
     "Set the distance scheme (\"angle\") used by KtJet",
     &KtJetFinder::_distanceScheme, 1, false, false);
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



  static Switch<KtJetFinder,int> interfaceRecombinationScheme
    ("RecombinationScheme",
     "Set the recombination scheme used by KtJet",
     &KtJetFinder::_recombinationScheme, 1, false, false);
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

void KtJetFinder::use (const vector<Lorentz5Momentum>& evt, bool inclusive) {
  JetFinder::use(evt, inclusive);
  convert();
  if (inclusive)
    clusterInclusive();
  else
    cluster();
}

void KtJetFinder::convert () {
  _lastMomenta.clear();
  for(vector<Lorentz5Momentum>::const_iterator p = lastEvent().begin();
      p != lastEvent().end(); ++p) {
    _lastMomenta.push_back(KtJet::KtLorentzVector(p->x()/GeV,
						  p->y()/GeV,
						  p->z()/GeV,
						  p->t()/GeV
						  ));
  }
}

void KtJetFinder::convert (const vector<KtJet::KtLorentzVector>& ktvectors) {
  vector<Lorentz5Momentum> tmp;
  for(vector<KtJet::KtLorentzVector>::const_iterator p =ktvectors.begin(); p != ktvectors.end(); ++p) {
    tmp.push_back(Lorentz5Momentum((*p).px()*GeV,(*p).py()*GeV,(*p).pz()*GeV,(*p).e()*GeV));
  }
  jets(tmp);
}

void KtJetFinder::cluster () {
  _lastKtEvent =
    std::auto_ptr<KtJet::KtEvent>
    (new KtJet::KtEvent(_lastMomenta,_collisionType,_distanceScheme,_recombinationScheme));
}

void KtJetFinder::clusterInclusive () {
  _lastKtEvent =
    std::auto_ptr<KtJet::KtEvent>
    (new KtJet::KtEvent(_lastMomenta,_collisionType,_distanceScheme,_recombinationScheme));
}

