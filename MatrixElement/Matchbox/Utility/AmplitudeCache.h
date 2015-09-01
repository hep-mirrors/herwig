// -*- C++ -*-
//
// AmplitudeCache.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_AmplitudeCache_H
#define HERWIG_AmplitudeCache_H

#include "Herwig/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "ThePEG/Config/algorithm.h"

namespace Herwig {

using namespace ThePEG;

namespace SpinorHelicity {

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Caching for amplitudes using spinor helicity techniques.
 *
 */
template<class AmplitudeKey>
class AmplitudeCache {

  typedef map<AmplitudeKey,pair<bool,Complex> > AmplitudeCacheMap;
  typedef map<AmplitudeKey,pair<bool,LorentzVector<Complex> > > CurrentCacheMap;

  /**
   * The number of points
   */
  int theNPoints;

  /**
   * The energy scale to obtain dimensionless
   * quantities.
   */
  mutable Energy theScale;

  /**
   * Masses indexed by partons
   */
  mutable vector<double> theMasses;

  /**
   * Momenta indexed by partons
   */
  mutable vector<LorentzMomentum> theMomenta;

  /**
   * Crossing signs indexed by partons
   */
  mutable vector<int> theCrossingSigns;

  /**
   * Plus spinors indexed by partons
   */
  mutable vector<PlusSpinor> thePlusSpinors;

  /**
   * Plus conjugate spinors indexed by partons
   */
  mutable vector<PlusConjugateSpinor> thePlusConjugateSpinors;

  /**
   * Invariants indexed by partons
   */
  mutable vector<vector<double> > theInvariants;

  /**
   * Flag products to be recalculated
   */
  mutable vector<vector<bool> > getInvariant;

  /**
   * Spinor products indexed by partons
   */
  mutable vector<vector<Complex> > thePlusProducts;

  /**
   * Flag products to be recalculated
   */
  mutable vector<vector<bool> > getPlusProduct;

  /**
   * Spinor currents indexed by partons
   */
  mutable vector<vector<LorentzVector<Complex> > > thePlusCurrents;

  /**
   * Flag currents to be recalculated
   */
  mutable vector<vector<bool> > getPlusCurrent;

  /**
   * Cache intermediate amplitudes by index
   */
  mutable AmplitudeCacheMap theCachedAmplitudes;

  /**
   * The last query for a cached amplitude
   */
  mutable typename AmplitudeCacheMap::iterator theLastAmplitude;

  /**
   * Cache intermediate currents by index
   */
  mutable CurrentCacheMap theCachedCurrents;

  /**
   * The last query for a cached current
   */
  mutable typename CurrentCacheMap::iterator theLastCurrent;

  /**
   * Helper to reset flags
   */
  struct boolResetter {
    void operator()(vector<bool>::reference flag) const {
      flag = true;
    }
    void operator()(pair<const AmplitudeKey,pair<bool,Complex> >& flag) const {
      flag.second.first = true;
    }
    void operator()(pair<const AmplitudeKey,pair<bool,LorentzVector<Complex> > >& flag) const {
      flag.second.first = true;
    }
  };

  /**
   * Helper to reset flags
   */
  struct boolVectorResetter {
    void operator()(vector<bool>& flags) const {
      for_each(flags.begin(),flags.end(),boolResetter());
    }
  };

public:

  /**
   * Constructor
   */
  AmplitudeCache()
    : theNPoints(0) {}

  /**
   * Prepare for n-point amplitude
   */
  void nPoints(int n);

  /**
   * Return the number of points
   */
  int nPoints() const {
    return theNPoints;
  }

  /**
   * Set the energy scale to obtain dimensionless
   * quantities and flag all quantities to be recalculated.
   */
  void amplitudeScale(Energy s) const;

  /**
   * Set the momentum for the k'th parton
   * and its associated mass.
   */
  void momentum(int k, const LorentzMomentum& p,
		bool getSpinors = true,
		Energy mass = ZERO) const;

  /**
   * Reset flags
   */
  void reset() const;

public:

  /**
   * Return the momentum for the k'th parton
   */
  LorentzVector<double> momentum(int k) const { return theMomenta[k]/theScale; }

  /**
   * Get the energy scale to obtain dimensionless
   * quantities and flag all quantities to be recalculated.
   */
  Energy amplitudeScale() const { return theScale; }

  /**
   * Return the mass associated to the k'th parton
   */
  double mass(int k) const { return theMasses[k]; }

  /**
   * Return the crossing sign for the
   * i'th parton
   */
  int crossingSign(int i) const { return theCrossingSigns[i]; }

  /**
   * Return the crossing sign for the
   * i'th and j'th parton
   */
  double crossingSign(int i, int j) const { return theCrossingSigns[i]*theCrossingSigns[j]; }

  /**
   * Return (ij)
   */
  double invariant(int i, int j) const {
    if ( i== j )
      return 0.;
    if ( i > j )
      swap(i,j);
    if ( getInvariant[i][j] ) {
      getInvariant[i][j] = false;
      theInvariants[i][j] = 2.*(momentum(i)*momentum(j));
    }
    return theInvariants[i][j];
  }

  /**
   * Return <ij>
   */
  Complex plusProduct(int i, int j) const {
    if ( i== j )
      return 0.;
    bool swapij = (i > j);
    if ( swapij )
      swap(i,j);
    if ( getPlusProduct[i][j] ) {
      getPlusProduct[i][j] = false;
      thePlusProducts[i][j] = 
	PlusSpinorProduct(thePlusConjugateSpinors[i],
			  thePlusSpinors[j]).eval() / theScale;
    }
    return swapij ? -thePlusProducts[i][j] : thePlusProducts[i][j];
  }

  /**
   * Return [ij]
   */
  Complex minusProduct(int i, int j) const {
    if ( i== j )
      return 0.;
    return -crossingSign(i,j)*conj(plusProduct(i,j));
  }

  /**
   * Return <i|\gamma^\mu|j]
   */
  LorentzVector<Complex> plusCurrent(int i, int j) const {
    bool swapij = (i > j);
    if ( swapij )
      swap(i,j);
    if ( getPlusCurrent[i][j] ) {
      getPlusCurrent[i][j] = false;
      if ( i != j ) {
	thePlusCurrents[i][j] = 
	  PlusSpinorCurrent(thePlusConjugateSpinors[i],
			    MinusSpinor(theMomenta[j])).eval() / theScale;
      } else {
	thePlusCurrents[i][j] = 2.*momentum(i);
      }
    }
    return swapij ? crossingSign(i,j)*thePlusCurrents[i][j].conjugate() : thePlusCurrents[i][j];
  }

  /**
   * Return [i|\gamma^\mu|j>
   */
  LorentzVector<Complex> minusCurrent(int i, int j) const {
    return plusCurrent(j,i);
  }

public:

  /**
   * Return true, if the given amplitude
   * needs to be recalculated.
   */
  bool getAmplitude(const AmplitudeKey& key) const {
    static Complex czero;
    if ( ( theLastAmplitude = theCachedAmplitudes.find(key) )
	 == theCachedAmplitudes.end() ) {
      theLastAmplitude = theCachedAmplitudes.insert(make_pair(key,make_pair(true,czero))).first;
    }
    return theLastAmplitude->second.first;
  }

  /**
   * Cache an amplitude
   */
  void cacheAmplitude(Complex amp) const {
    theLastAmplitude->second = make_pair(false,amp);
  }

  /**
   * Return a cached amplitude
   */
  const Complex& cachedAmplitude() const {
    return theLastAmplitude->second.second;
  }

  /**
   * Return true, if the given current
   * needs to be recalculated.
   */
  bool getCurrent(const AmplitudeKey& key) const {
    static LorentzVector<Complex> czero;
    if ( ( theLastCurrent = theCachedCurrents.find(key) )
	 == theCachedCurrents.end() ) {
      theLastCurrent = theCachedCurrents.insert(make_pair(key,make_pair(true,czero))).first;
    }
    return theLastCurrent->second.first;
  }

  /**
   * Cache an current
   */
  void cacheCurrent(const LorentzVector<Complex>& curr) const {
    theLastCurrent->second = make_pair(false,curr);
  }

  /**
   * Return a cached current
   */
  const LorentzVector<Complex>& cachedCurrent() const {
    return theLastCurrent->second.second;
  }

};

}

}

#include "AmplitudeCache.tcc"

#endif // HERWIG_AmplitudeCache_H
