// -*- C++ -*-
//
// DensityOperator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_DensityOperator_H
#define Herwig_DensityOperator_H
//
// This is the declaration of the DensityOperator class.
//

#include "ThePEG/Handlers/HandlerBase.h"

#include "Herwig/MatrixElement/Matchbox/Utility/ColourBasis.h"

#include <tuple>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace Herwig {

using namespace ThePEG;

typedef boost::numeric::ublas::vector<Complex> CVector;



/**
 * Here is the documentation of the DensityOperator class.
 *
 * @see \ref DensityOperatorInterfaces "The interfaces"
 * defined for DensityOperator.
 */
class DensityOperator: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DensityOperator();

  /**
   * The destructor.
   */
  virtual ~DensityOperator();
  //@}

public:
  
  /**
   * Clears theDensityOperatorMap.
   */
  void clear();

  /**
   * Prepare for the given sub process.
   */
  void prepare(const cPDVector&);

  /**
   * Fill the density operator for the given hard subprocess, summing over all
   * helicity configurations.
   */
  void fill(const Ptr<MatchboxXComb>::ptr,
	    const cPDVector&, const vector<Lorentz5Momentum>& momenta);
  
  /**
   * Evolve the density operator, by 
   * M_{n+1} = -\sum_{i,k}{-4*pi*alpha_s/Ti2*V_{ij,k} T_{i,n}M_nT_{k,n}^\dag},
   * see arXiv:1206.0180 eq. (5), note that the pi*pj factor is assumed to be
   * included in V_{ij,k}.
   */
  void evolve(const map<pair<size_t,size_t>,Complex>& Vijk, 
	      const cPDVector& before, 
	      const cPDVector& after,
	      const map<std::tuple<size_t,size_t,size_t>,map<size_t,size_t> >& emissionsMap,
	      const bool splitAGluon,
	      const bool initialGluonSplitting);
  
  /**
   * Calculate the colour matrix element correction.
   * -(1+delta(i,gluon))/Ti^2 Tr(Sn+1 Ti Mn Tk^dagger)/Tr(Sn Mn)
   * where the bracket in front compensates for the gluon symmetry factor,
   * Ti^2 is C_f or C_a, Sn+1 is the matrix of scalar products, and
   * Ti is the radiation matrix.
   * The first arg contains (emitter index, spectator index, emission pid)
   *
   */
  double colourMatrixElementCorrection(const std::tuple<size_t,size_t,long>& ikemission,
				       const cPDVector& particles);

  /**
   * Checking colour conservation for the colour matrix element corrections.
   */
  void colourConservation(const cPDVector& particles);

  /**
   * Get the colour basis.
   */
  Ptr<ColourBasis>::tptr colourBasis() { return theColourBasis; }
  
  /**
   * Get the colour basis.
   */
  const Ptr<ColourBasis>::tptr colourBasis() const { return theColourBasis; }
  
  /**
   * Set the colour basis.
   */
  void colourBasis(Ptr<ColourBasis>::ptr ptr) { theColourBasis = ptr; }
  
  /**
   * Get the correlator map.
   */
  const map<pair<vector<PDT::Colour>,pair<size_t,size_t> >,double>& correlatorMap() const {
    return theCorrelatorMap;
  }

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
  * Number of colours used in colourNorm.
  */
  double Nc;

  /**
   * QCD vertex normalization.
   */
  double TR;

  /**
   * Normalization of colour charges \mathbf{T}_{ij}^2.
   */
  double colourNorm(const cPDPtr particle);

  /**
   * Fast evaluation of Tij*Mn, where a Tij is the matrix from ColourBasis::charge,
   * which is a sparse matrix, and Mn is the density operator, a dense matrix.
   *
   */
  matrix<Complex> prodSparseDense(const compressed_matrix<double>&,
				  const matrix<Complex>&);
  /**
   * Fast evaluation of TijMn*Tkdagger, where a TijMn is the result from the method
   * prodSparseDense, a dense matrix, and Tkdagger is the transponse conjugate of
   * the matrix from ColourBasis::charge, a sparse matrix.
   *
   */
  matrix<Complex> prodDenseSparse(const matrix<Complex>&,
				  const compressed_matrix<double>&);

  /**
   * Boosts a vector of momenta to the rest frame of the initial pair
   * of particles (the first 2 elements of the argument vector). Returns
   * the boosted vectors
   */
  vector<Lorentz5Momentum> boostToRestFrame(const vector<Lorentz5Momentum>& momenta);

  /**
   * Boosts a vector of momenta to the rest frame of the initial pair
   */
  bool compareMomentum(const Lorentz5Momentum& p, const Lorentz5Momentum& q);
  
  /**
   * Mapping of colour structures to density operator matrices.
   *
   */
  map<vector<PDT::Colour>,matrix<Complex> > theDensityOperatorMap;
  
  /**
   * Mapping of colour structures and legs to colour correlators. 
   */
  map<pair<vector<PDT::Colour>,pair<size_t,size_t> >,double> theCorrelatorMap;

  /**
   * A map from the hard subprocess particles to a map of amplitude colour
   * basis order to the normal ordered colour basis. 
   */
  map<cPDVector, map<size_t,size_t> > theColourBasisToColourBasisMap;
  
  /**
   * Colour basis used.
   */
  Ptr<ColourBasis>::ptr theColourBasis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DensityOperator & operator=(const DensityOperator &) = delete;

};

}

#endif /* Herwig_DensityOperator_H */
