// -*- C++ -*-
//
// AmplitudeCache.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

namespace Herwig {
namespace SpinorHelicity {

template<typename AmplitudeKey>
void AmplitudeCache<AmplitudeKey>::nPoints(int n) {

  assert( n <= MAX_N );

  theNPoints = n;

  theMasses.fill({});
  theMomenta.fill({});
  theCrossingSigns.fill({});
  thePlusSpinors.fill(PlusSpinor());
  thePlusConjugateSpinors.fill(PlusConjugateSpinor());  
  theInvariants.fill({});
  thePlusProducts.fill({});
  thePlusCurrents.fill({});

  reset();    
}

template<typename AmplitudeKey>
void AmplitudeCache<AmplitudeKey>::amplitudeScale(Energy s) const { 
  theScale = s;
  reset();
}

template<typename AmplitudeKey>
void AmplitudeCache<AmplitudeKey>::momentum(int k, const LorentzMomentum& p, 
			      bool getSpinors,
			      Energy mass) const {
  theMasses[k] = mass/theScale;
  theMomenta[k] = p;
  if ( getSpinors ) {
    theCrossingSigns[k] = p.t() > ZERO ? 1 : -1;
    thePlusSpinors[k] = PlusSpinor(p);
    thePlusConjugateSpinors[k] = PlusConjugateSpinor(p);
  }
}

template<typename AmplitudeKey>
void AmplitudeCache<AmplitudeKey>::reset() const {
  getInvariant.fill(true);
  getPlusProduct.fill(true);
  getPlusCurrent.fill(true);
  for_each(theCachedAmplitudes.begin(),theCachedAmplitudes.end(),boolResetter());
  for_each(theCachedCurrents.begin(),theCachedCurrents.end(),boolResetter());
}

}}
