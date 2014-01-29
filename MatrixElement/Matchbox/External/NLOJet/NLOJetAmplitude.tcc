// -*- C++ -*-
//
// NLOJetAmplitude.tcc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetAmplitude class.
//

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SpinorHelicity.h"

#include "NLOJetPhasespace.h"

namespace Herwig {

template<unsigned int N, unsigned int I, unsigned int F>
NLOJetAmplitude<N,I,F>::NLOJetAmplitude() 
  : MatchboxAmplitude(), calculateInvariants(true), theIdenticalQuarks(1) {}

template<unsigned int N, unsigned int I, unsigned int F>
NLOJetAmplitude<N,I,F>::~NLOJetAmplitude() {}

template<unsigned int N, unsigned int I, unsigned int F>
const typename NLOJetAmplitude<N,I,F>::NLOEvent& 
NLOJetAmplitude<N,I,F>::lastNLOEvent() const {
  if ( lastXComb().hasMeta(NLOMetaKeys::HadronicEvent) ) {
    return lastXComb().template meta<NLOEvent>(NLOMetaKeys::HadronicEvent);
  }

  NLOMomentumConverter converter;

  theLastNLOEvent.resize(meMomenta().size()-2-F);

  theLastNLOEvent[-1] = converter(meMomenta()[0]);
  theLastNLOEvent[0] = converter(meMomenta()[1]);

  for ( unsigned int k = 1; k <= F; ++k )
    theLastNLOEvent[-1-k] = converter(meMomenta()[k+1]);

  for ( unsigned int k = 1; k <= meMomenta().size()-2-F; ++k )
    theLastNLOEvent[k] = converter(meMomenta()[k+F+1]);

  return theLastNLOEvent;
}

struct matchQuark {
  long checkid;
  explicit matchQuark(int ci) : checkid(ci) {}
  inline bool operator()(tcPDPtr data) const {
    return abs(data->id()) == checkid;
  }
};

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetAmplitude<N,I,F>::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr) {
  if ( !calculateInvariants )
    return;

  theLastInvariants.calculate(lastNLOEvent());
  theIdenticalQuarks = 0;
  theOutgoingIdenticalQuarks = 0;
  for ( int q = 1; q < 6; ++q ) {
    theIdenticalQuarks = 
      max(theIdenticalQuarks,count_if(mePartonData().begin(),
				      mePartonData().end(),
				      matchQuark(q)));
    theOutgoingIdenticalQuarks = 
      max(theOutgoingIdenticalQuarks,count_if(mePartonData().begin()+2,
					      mePartonData().end(),
					      matchQuark(q)));
  }
  theIdenticalQuarks /= 2;
  calculateInvariants = false;
}

template<unsigned int N, unsigned int I, unsigned int F>
double NLOJetAmplitude<N,I,F>::colourCorrelatedME2(pair<int,int> ij) const {
  double Nc = generator()->standardModel()->Nc();
  double cfac = mePartonData()[ij.first]->id() == ParticleID::g ? Nc : (sqr(Nc)-1.)/(2.*Nc);
  if ( !calculateColourCorrelator(ij) )
    return lastColourCorrelator(ij)/cfac;
  pair<int,int> ijp(ij.first-1,ij.second-1);
  double res = 
    pow(4.*Constants::pi*SM().alphaS(),double(orderInAlphaS()))*
    pow(lastSHat()/GeV2,double(mePartonData().size()-4))*
    treeLevelCC(ijp,crossingMap());
  lastColourCorrelator(ij,res);
  return res/cfac;
}

template<unsigned int N, unsigned int I, unsigned int F>
double NLOJetAmplitude<N,I,F>::spinColourCorrelatedME2(pair<int,int> ij,
						       const SpinCorrelationTensor& c) const {

  using namespace SpinorHelicity;

  double Nc = generator()->standardModel()->Nc();
  double cfac = mePartonData()[ij.first]->id() == ParticleID::g ? Nc : (sqr(Nc)-1.)/(2.*Nc);

  pair<int,int> ijp(ij.first-1,ij.second-1);
  pair<double,Complex> dc = treeLevelSCC(ijp,crossingMap());

  Lorentz5Momentum p = meMomenta()[ij.first];
  Lorentz5Momentum n = meMomenta()[ij.second];

  LorentzVector<complex<Energy> > num =
    PlusSpinorCurrent(PlusConjugateSpinor(n),MinusSpinor(p)).eval();

  complex<Energy> den =
    sqrt(2.)*PlusSpinorProduct(PlusConjugateSpinor(n),PlusSpinor(p)).eval();

  LorentzVector<Complex> polarization(num.x()/den,num.y()/den,num.z()/den,num.t()/den);

  Complex pFactor = (polarization*c.momentum())/sqrt(abs(c.scale()));

  double res = dc.first*(-c.diagonal() + (c.scale() > ZERO ? 1. : -1.)*norm(pFactor));

  res += 2.*(c.scale() > ZERO ? 1. : -1.)*real(sqr(pFactor)*dc.second);

  return 
    pow(4.*Constants::pi*SM().alphaS(),double(orderInAlphaS()))*
    pow(lastSHat()/GeV2,double(mePartonData().size()-4))*
    res/cfac;

}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetAmplitude<N,I,F>::flushCaches() {
  MatchboxAmplitude::flushCaches();
  calculateInvariants = true;
}

static long pows10[6] = {

  1,10,100,1000,10000,100000

};

template<unsigned int N, unsigned int I, unsigned int F>
pair<unsigned int, unsigned int> NLOJetAmplitude<N,I,F>::countColoured(const PDVector& pdIn) const {
  PDVector pd = pdIn;
  if ( pd.size() < 2 )
    return make_pair(0,0);
  if ( pd[0]->CC() )
    pd[0] = pd[0]->CC();
  if ( pd[1]->CC() )
    pd[1] = pd[1]->CC();
  long quarkSum = 0;
  unsigned int ng = 0;
  unsigned int nq = 0;
  for ( PDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p ) {
    if ( (**p).id() == 21 )
      ng += 1;
    if ( abs((**p).id()) < 6 ) {
      nq += 1;
      if ( (**p).id() > 0 )
	quarkSum += pows10[abs((**p).id())];
      else
	quarkSum -= pows10[abs((**p).id())];
    }
  }
  if ( nq % 2 != 0 || quarkSum != 0 )
    return make_pair(0,0);
  return make_pair(nq,ng);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetAmplitude<N,I,F>::persistentOutput(PersistentOStream &) const {
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetAmplitude<N,I,F>::persistentInput(PersistentIStream &, int) {
}

template<unsigned int N, unsigned int I, unsigned int F>
void NLOJetAmplitude<N,I,F>::Init() {

  static ClassDocumentation<NLOJetAmplitude<N,I,F> > documentation
    ("NLOJetAmplitude provides an interface to nlojet's amplitude classes.");

}

}
