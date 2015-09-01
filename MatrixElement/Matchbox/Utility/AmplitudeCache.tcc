// -*- C++ -*-
//
// AmplitudeCache.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

namespace Herwig {
namespace SpinorHelicity {

template<class AmplitudeKey>
void AmplitudeCache<AmplitudeKey>::nPoints(int n) {

  theNPoints = n;

  theMasses.clear();
  theMomenta.clear();
  theCrossingSigns.clear();
  thePlusSpinors.clear();
  thePlusConjugateSpinors.clear();  
  theInvariants.clear();
  thePlusProducts.clear();
  thePlusCurrents.clear();
  getInvariant.clear();
  getPlusProduct.clear();
  getPlusCurrent.clear();

  theMasses.resize(n);
  theMomenta.resize(n);
  theCrossingSigns.resize(n);
  thePlusSpinors.resize(n);
  thePlusConjugateSpinors.resize(n);  
  theInvariants.resize(n,vector<double>(n));
  thePlusProducts.resize(n,vector<Complex>(n));
  thePlusCurrents.resize(n,vector<LorentzVector<Complex> >(n));
  getInvariant.resize(n,vector<bool>(n));
  getPlusProduct.resize(n,vector<bool>(n));
  getPlusCurrent.resize(n,vector<bool>(n));

  reset();    
}

template<class AmplitudeKey>
void AmplitudeCache<AmplitudeKey>::amplitudeScale(Energy s) const { 
  theScale = s;
  reset();
}

template<class AmplitudeKey>
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

template<class AmplitudeKey>
void AmplitudeCache<AmplitudeKey>::reset() const {
  for_each(getInvariant.begin(),getInvariant.end(),boolVectorResetter());
  for_each(getPlusProduct.begin(),getPlusProduct.end(),boolVectorResetter());
  for_each(getPlusCurrent.begin(),getPlusCurrent.end(),boolVectorResetter());
  for_each(theCachedAmplitudes.begin(),theCachedAmplitudes.end(),boolResetter());
  for_each(theCachedCurrents.begin(),theCachedCurrents.end(),boolResetter());
}

}}
