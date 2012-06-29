// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "FastJetFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FastJetFinder.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

FastJetFinder::~FastJetFinder() {}

void FastJetFinder::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _jetFinder << _strategy << _recombinationScheme;
}

void FastJetFinder::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _jetFinder >> _strategy >> _recombinationScheme;
}

ClassDescription<FastJetFinder> FastJetFinder::initFastJetFinder;
// Definition of the static class description member.

void FastJetFinder::Init() {

  static ClassDocumentation<FastJetFinder> documentation
    ("JetFinder using the FastJetFinder library.");


  static Switch<FastJetFinder,int> interfaceJetFinder
    ("JetFinder",
     "Set the jet finder type",
     &FastJetFinder::_jetFinder, 0, false, false);
  static SwitchOption interfaceJetFinderkt_algorithm
    (interfaceJetFinder,
     "kt_algorithm",
     "kt_algorithm",
     0);
  static SwitchOption interfaceJetFindercambridge_algorithm
    (interfaceJetFinder,
     "cambridge_algorithm",
     "cambridge_algorithm",
     1);


  static Switch<FastJetFinder,int> interfaceStrategy
    ("Strategy",
     "The strategy used by FastJetFinder",
     &FastJetFinder::_strategy, 1, false, false);
  static SwitchOption interfaceStrategyN2MinHeapTiled
    (interfaceStrategy,
     "N2MinHeapTiled",
     "N2MinHeapTiled",
     -4);
  static SwitchOption interfaceStrategyN2Tiled
    (interfaceStrategy,
     "N2Tiled",
     "N2Tiled",
     -3);
  static SwitchOption interfaceStrategyN2PoorTiled
    (interfaceStrategy,
     "N2PoorTiled",
     "N2PoorTiled",
     -2);  
  static SwitchOption interfaceStrategyN2Plain
    (interfaceStrategy,
     "N2Plain",
     "N2Plain",
     -1);  
  static SwitchOption interfaceStrategyN3Dumb
    (interfaceStrategy,
     "N3Dumb",
     "N3Dumb",
     0);  
  static SwitchOption interfaceStrategyBest
    (interfaceStrategy,
     "Best",
     "Best",
     1);  
  static SwitchOption interfaceStrategyNlnN
    (interfaceStrategy,
     "NlnN",
     "NlnN",
     2);  
  static SwitchOption interfaceStrategyNlnN3pi
    (interfaceStrategy,
     "NlnN3pi",
     "NlnN3pi",
     3);  
  static SwitchOption interfaceStrategyNlnN4pi
    (interfaceStrategy,
     "NlnN4pi",
     "NlnN4pi",
     4);  


  static Switch<FastJetFinder,int> interfaceRecombinationScheme
    ("RecombinationScheme",
     "The recombination scheme to be used",
     &FastJetFinder::_recombinationScheme, 0, false, false);
  static SwitchOption interfaceRecombinationSchemeE_scheme
    (interfaceRecombinationScheme,
     "E_scheme",
     "E_scheme",
     0);
  static SwitchOption interfaceRecombinationSchemept_scheme
    (interfaceRecombinationScheme,
     "pt_scheme",
     "pt_scheme",
     1);
  static SwitchOption interfaceRecombinationSchemept2_scheme
    (interfaceRecombinationScheme,
     "pt2_scheme",
     "pt2_scheme",
     2);
  static SwitchOption interfaceRecombinationSchemeEt_scheme
    (interfaceRecombinationScheme,
     "Et_scheme",
     "Et_scheme",
     3);
  static SwitchOption interfaceRecombinationSchemeEt2_scheme
    (interfaceRecombinationScheme,
     "Et2_scheme",
     "Et2_scheme",
     4);
  static SwitchOption interfaceRecombinationSchemeBIpt_scheme
    (interfaceRecombinationScheme,
     "BIpt_scheme",
     "BIpt_scheme",
     5);
  static SwitchOption interfaceRecombinationSchemeBIpt2_scheme
    (interfaceRecombinationScheme,
     "BIpt2_scheme",
     "BIpt2_scheme",
     6);


}

void FastJetFinder::convert () {
  _lastPseudojets.clear();
  Lorentz5Momentum sum = Lorentz5Momentum (0.*GeV,0.*GeV,0.*GeV,0.*GeV);
  for(vector<Lorentz5Momentum>::const_iterator p = lastEvent().begin();
      p != lastEvent().end(); ++p) {
    sum = sum + (*p).momentum();
    _lastPseudojets.push_back(fastjet::PseudoJet((*p).momentum().x()/GeV,
						 (*p).momentum().y()/GeV,
						 (*p).momentum().z()/GeV,
						 (*p).momentum().t()/GeV
						 ));
  }
  _E2vis = sum.m2();
}

void FastJetFinder::convert (const vector<fastjet::PseudoJet>& pv) {
  list<Lorentz5Momentum> tmp;
  for(vector<fastjet::PseudoJet>::const_iterator p =pv.begin(); p != pv.end(); ++p) {
    tmp.push_back(Lorentz5Momentum((*p).px()*GeV,(*p).py()*GeV,(*p).pz()*GeV,(*p).e()*GeV));
  }
  jets(tmp);
}

void FastJetFinder::cluster () {
  _lastClusterSequence =
    std::auto_ptr<fastjet::ClusterSequence>
    (new fastjet::ClusterSequence(_lastPseudojets,_jetDefinition));
}

void FastJetFinder::use (const vector<Lorentz5Momentum>& evt, bool) {
  JetFinder::use(evt);
  _lastPseudojets.clear();
  convert();
  cluster();
}

