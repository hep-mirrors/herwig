// -*- C++ -*-
//
// NLOJetPhasespace.tcc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetPhasespace class.
//

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "nlo++/bits/psg-phasespace_n0i0f0.h"
#include "nlo++/bits/psg-phasespace_n0i2f0.h"
//#include "nlo++/bits/psg-phasespace_n1i1f0.h"
#include "nlo++/bits/psg-phasespace_n2i2f0.h"

namespace Herwig {

template<unsigned int N, unsigned int I, unsigned int F>
NLOJetPhasespace<N,I,F>::NLOJetPhasespace() 
  : MatchboxPhasespace(), theNLOPhasespace(0), nloEvent(0), nloHardEvent(0) {}

template<unsigned int N, unsigned int I, unsigned int F>
NLOJetPhasespace<N,I,F>::~NLOJetPhasespace() {
  if ( theNLOPhasespace ) {
    delete theNLOPhasespace;
    theNLOPhasespace = 0;
  }
  
  for(typename  map <  XCPtr, NLOEvent* >::iterator it= theXCtoEvents.begin();
      it != theXCtoEvents.end(); ++it )
    delete (*it).second;
  nloEvent = 0;
  for(typename  map <  XCPtr, NLOEvent* >::iterator it= theXCtoHardEvents.begin();
      it != theXCtoHardEvents.end(); ++it )
    delete (*it).second;
  nloHardEvent = 0;
}

template<unsigned int N, unsigned int I, unsigned int F>
IBPtr NLOJetPhasespace<N,I,F>::clone() const {
  return new_ptr(*this);
}

template<unsigned int N, unsigned int I, unsigned int F>
IBPtr NLOJetPhasespace<N,I,F>::fullclone() const {
  return new_ptr(*this);
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetPhasespace<N,I,F>::doinit() {
  MatchboxPhasespace::doinit();
  double s = sqr(generator()->maximumCMEnergy()/GeV);
  if ( theNLOPhasespace ) 
    delete theNLOPhasespace;
  theNLOPhasespace = new NLOPhasespace(&nloRnd,s);
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetPhasespace<N,I,F>::doinitrun() {
  MatchboxPhasespace::doinitrun();
  double s = sqr(generator()->maximumCMEnergy()/GeV);
  if ( theNLOPhasespace ) 
    delete theNLOPhasespace;
  theNLOPhasespace = new NLOPhasespace(&nloRnd,s);
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetPhasespace<N,I,F>::prepare(tStdXCombPtr xc, bool) {
  theLastXComb = xc;

  if( theXCtoEvents.find( lastXCombPtr() ) != theXCtoEvents.end() ){
    nloEvent=theXCtoEvents[lastXCombPtr()];
  }
  else{
    nloEvent= new NLOEvent(nHard-F);
    theXCtoEvents[lastXCombPtr()]=nloEvent;
  }
  if( theXCtoHardEvents.find(lastXCombPtr()) != theXCtoHardEvents.end() ){
    nloHardEvent=theXCtoHardEvents[lastXCombPtr()];
  }
  else{
    nloHardEvent= new NLOEvent(mePartonData().size()-2-N-F);
    theXCtoHardEvents[lastXCombPtr()]=nloHardEvent;
  }
}

template<unsigned int N, unsigned int I, unsigned int F>
int NLOJetPhasespace<N,I,F>::nDim(int nFinal) const { 
  // rambo + x's + dipoles
  return
    4*(N+nHard) + 2*I + 8*(nFinal-N-nHard);
}

template<unsigned int N, unsigned int I, unsigned int F>
double NLOJetPhasespace<N,I,F>::generateTwoToNKinematics(const double* r,
							 vector<Lorentz5Momentum>& momenta) {
  nloRnd.numbers(r);
  double weight = (*nloPhasespace())(*nloHardEvent);
  if ( weight == 0. )
    return 0.;
  weight *= (*nloPhasespace())(*nloHardEvent,*nloEvent);
  if ( weight == 0. )
    return 0.;

  NLOMomentumConverter converter;

  momenta[0] = converter((*nloEvent)[-1]);
  momenta[1] = converter((*nloEvent)[0]);

  for ( unsigned int k = 1; k <= F; ++k )
    momenta[k+1] = converter((*nloEvent)[-1-k]);

  for ( unsigned int k = 1; k <= momenta.size()-2-F; ++k )
    momenta[k+F+1] = converter((*nloEvent)[k]);

  lastXCombPtr()->meta(NLOMetaKeys::HadronicEvent,*nloEvent);

  double x1 = momenta[0].plus()/lastParticles().first->momentum().plus();
  double x2 = momenta[1].minus()/lastParticles().second->momentum().minus();

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat((momenta[0]+momenta[1]).m2());

  return pow(lastSHat()/GeV2,-(double)(momenta.size()-4))*weight/(x1*x2);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetPhasespace<N,I,F>::persistentOutput(PersistentOStream & os) const {
  os << nHard;
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetPhasespace<N,I,F>::persistentInput(PersistentIStream & is, int) {
  is >> nHard;
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetPhasespace<N,I,F>::Init() {

  static ClassDocumentation<NLOJetPhasespace<N,I,F> > documentation
    ("NLOJetPhasespace provides an interface to "
     "nlojet's phasespace generator classes.");


  static Parameter<NLOJetPhasespace,unsigned int> interfaceNHard
    ("NHard",
     "The number of outgoing partons from which on dipole generation is used.",
     &NLOJetPhasespace::nHard, 2, 1, 0,
     false, false, Interface::lowerlim);

}

}
