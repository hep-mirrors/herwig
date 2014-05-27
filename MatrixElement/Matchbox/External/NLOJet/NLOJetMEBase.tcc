// -*- C++ -*-
//
// NLOJetMEBase.tcc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEBase class.
//

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"



#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {

template<unsigned int N, unsigned int I, unsigned int F>
NLOJetMEBase<N,I,F>::NLOJetMEBase() 
  : MatchboxMEBase() {}

template<unsigned int N, unsigned int I, unsigned int F>
NLOJetMEBase<N,I,F>::~NLOJetMEBase() {}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetMEBase<N,I,F>::getDiagrams() const {

  useMe();

  set<PDPtr> inQuarks;

  for ( PDVector::const_iterator p = subProcess().legs.begin();
	p != subProcess().legs.end(); ++p ) {
    if ( abs((**p).id()) < 6 ) {
      if ( (**p).id() > 0 )
	inQuarks.insert(*p);
      if ( (**p).id() < 0 )
	inQuarks.insert((**p).CC());
    }
  }

  quark.clear();
  copy(inQuarks.begin(),inQuarks.end(),back_inserter(quark));

  antiquark.clear();
  for ( PDVector::const_iterator q = quark.begin();
	q != quark.end(); ++q )
    antiquark.push_back((**q).CC());

  doGetDiagrams();

}

struct EqualPID {
  inline bool operator()(cPDPtr a, cPDPtr b) const {
    return a->id() == b->id();
  }
};

struct SortPID {
  inline bool operator()(cPDPtr a, cPDPtr b) const {
    return a->id() < b->id();
  }
};

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetMEBase<N,I,F>::addSafe(DiagPtr diag) const {

  assert(subProcess().legs.size() == diag->partons().size());
  if ( subProcess().legs[0]->id() != diag->partons()[0]->id() ||
       subProcess().legs[1]->id() != diag->partons()[1]->id() )
    return;
  multiset<cPDPtr,SortPID> sublegs;
  copy(subProcess().legs.begin() + 2,subProcess().legs.end(),inserter(sublegs));
  multiset<cPDPtr,SortPID> diaglegs;
  copy(diag->partons().begin() + 2,diag->partons().end(),inserter(diaglegs));
  if ( equal(sublegs.begin(),sublegs.end(),diaglegs.begin(),EqualPID()) )
    add(diag);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetMEBase<N,I,F>::doinit() {
  MatchboxMEBase::doinit();
  theNLOJetAmplitude = 
    dynamic_ptr_cast<typename Ptr<NLOJetAmplitude<N,I,F> >::ptr>(matchboxAmplitude());
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetMEBase<N,I,F>::persistentOutput(PersistentOStream & os) const {
  os << theNLOJetAmplitude << quark << antiquark;
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetMEBase<N,I,F>::persistentInput(PersistentIStream & is, int) {
  is >> theNLOJetAmplitude >> quark >> antiquark;
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetMEBase<N,I,F>::Init() {

  static ClassDocumentation<NLOJetMEBase<N,I,F> > documentation
    ("NLOJetMEBase provides the base class for all NLO jet processes.",
     "Matrix elements have been calculated using nlojet++");

  static RefVector<NLOJetMEBase<N,I,F>,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &NLOJetMEBase<N,I,F>::quark, -1, false, false, true, true, false);

}

}
