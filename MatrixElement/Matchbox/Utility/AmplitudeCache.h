// -*- C++ -*-
//
// AmplitudeCache.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_AmplitudeCache_H
#define HERWIG_AmplitudeCache_H

#include "Herwig/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "ThePEG/Config/algorithm.h"
#include <array>

namespace Herwig {

using namespace ThePEG;
using std::array;

namespace SpinorHelicity {

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Caching for amplitudes using spinor helicity techniques.
 *
 */
template<typename AmplitudeKey>
class AmplitudeCache {

  typedef map<AmplitudeKey,pair<bool,Complex> > AmplitudeCacheMap;
  typedef map<AmplitudeKey,pair<bool,LorentzVector<Complex> > > CurrentCacheMap;

  /**
   * Maximum N we can handle, SYM_N is storage size for a symmetric matrix of N * N elements
   */
  enum { MAX_N = 7, SYM_N = MAX_N*(MAX_N+1)/2 };

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
  mutable array<double,MAX_N> theMasses;

  /**
   * Momenta indexed by partons
   */
  mutable array<LorentzMomentum,MAX_N> theMomenta;

  /**
   * Crossing signs indexed by partons
   */
  mutable array<int,MAX_N> theCrossingSigns;

  /**
   * Plus spinors indexed by partons
   */
  mutable array<PlusSpinor,MAX_N> thePlusSpinors;

  /**
   * Plus conjugate spinors indexed by partons
   */
  mutable array<PlusConjugateSpinor,MAX_N> thePlusConjugateSpinors;

  /**
   * Invariants indexed by partons
   */
  mutable array<double,SYM_N> theInvariants;

  /**
   * Flag products to be recalculated
   */
  mutable array<bool,SYM_N> getInvariant;

  /**
   * Spinor products indexed by partons
   */
  mutable array<Complex,SYM_N> thePlusProducts;

  /**
   * Flag products to be recalculated
   */
  mutable array<bool,SYM_N> getPlusProduct;

  /**
   * Spinor currents indexed by partons
   */
  mutable array<LorentzVector<Complex>,SYM_N> thePlusCurrents;

  /**
   * Flag currents to be recalculated
   */
  mutable array<bool,SYM_N> getPlusCurrent;

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
   * Helper function to index symmetric arrays, assumes i <= j.
   * Usual indexing function (N*i + j) corrected by triangle number for i-th row.
   */
  inline size_t idx(size_t i, size_t j) const {
    return MAX_N * i - (i + 1) * i / 2 + j;
  }

  /**
   * Helper to reset flags
   */
  struct boolResetter {
    void operator()(pair<const AmplitudeKey,pair<bool,Complex> >& flag) const {
      flag.second.first = true;
    }
    void operator()(pair<const AmplitudeKey,pair<bool,LorentzVector<Complex> > >& flag) const {
      flag.second.first = true;
    }
  };

public:

  /**
   * Constructor
   */
  AmplitudeCache() : theNPoints(0) {}

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
    if ( i == j ) return 0.;
    if ( i > j  ) swap(i,j);
    if ( getInvariant[idx(i,j)] ) {
      getInvariant[idx(i,j)] = false;
      theInvariants[idx(i,j)] = 2.*(momentum(i)*momentum(j));
    }
    return theInvariants[idx(i,j)];
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
    if ( getPlusProduct[idx(i,j)] ) {
      getPlusProduct[idx(i,j)] = false;
      thePlusProducts[idx(i,j)] = 
	Complex(PlusSpinorProduct(thePlusConjugateSpinors[i],
				  thePlusSpinors[j]).eval() / theScale);
    }
    return swapij ? -thePlusProducts[idx(i,j)] : thePlusProducts[idx(i,j)];
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
    if ( getPlusCurrent[idx(i,j)] ) {
      getPlusCurrent[idx(i,j)] = false;
      if ( i != j ) {
	thePlusCurrents[idx(i,j)] = 
	  PlusSpinorCurrent(thePlusConjugateSpinors[i],
			    MinusSpinor(theMomenta[j])).eval() / theScale;
      } else {
	thePlusCurrents[idx(i,j)] = 2.*momentum(i);
      }
    }
    return swapij ? crossingSign(i,j)*thePlusCurrents[idx(i,j)].conjugate() : thePlusCurrents[idx(i,j)];
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
