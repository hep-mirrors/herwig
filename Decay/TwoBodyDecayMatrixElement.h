// -*- C++ -*-
//
// TwoBodyMatrixElement.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoBodyDecayMatrixElement_H
#define HERWIG_TwoBodyDecayMatrixElement_H
//
// This is the declaration of the TwoBodyDecayMatrixElement class.

#include "DecayMatrixElement.h"

namespace Herwig {

using namespace ThePEG;


/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  Implementation of the complex matrix element for a two-body decay
 * a decay.
 *  An arbitary number of external particles are supported.
 *  @see RhoDMatrix
 *  @see ProductionMatrixElement
 *  @see DecayVertex
 */
class TwoBodyDecayMatrixElement: public DecayMatrixElement {
  
public:
      
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  TwoBodyDecayMatrixElement() : DecayMatrixElement(2) {}

  /** 
   * Constructor for two body decay.
   * @param spinin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   */
  TwoBodyDecayMatrixElement(PDT::Spin spinin,PDT::Spin outspin1,PDT::Spin outspin2)
    : DecayMatrixElement(2,spinin) {
    outspin().push_back(outspin1);
    outspin().push_back(outspin2);
  }

  /** 
   * Constructor for arbitary body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param spinout \f$2S+1\f$ for the decay products.
   */
  TwoBodyDecayMatrixElement(PDT::Spin inspin,vector<PDT::Spin> spinout)
    : DecayMatrixElement(2,inspin) {
    assert(spinout.size()==2);
    outspin() = spinout;
  }

  /** 
   * Constructor for arbitary body decay.
   * @param extspin  \f$2S+1\f$ external particles.
   */
  TwoBodyDecayMatrixElement(vector<PDT::Spin> extspin)
    : DecayMatrixElement(2,extspin[0]) {
    assert(extspin.size()==3);
    outspin() = vector<PDT::Spin>(extspin.begin()+1,extspin.end());
  }
  //@}

public:

  /**
   * Spin Density matrices
   */
  //@{
  /** 
   * Calculate the decay matrix for this decay.
   * @param rhoout The \f$D\f$ matrix for this decay.
   */
  RhoDMatrix calculateDMatrix(const vector<RhoDMatrix> & rhoout) const;

  /** 
   * Calculate the \f$\rho\f$ matrix for a given outgoing particle.
   * @param ipart The outgoing particle the \f$\rho\f$ matrix is needed for
   * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
   * @param rhoout he \f$D\f$ matrices for the other decay products.
   */
  RhoDMatrix calculateRhoMatrix(int ipart,const RhoDMatrix & rhoin,
				const vector<RhoDMatrix> & rhoout) const;

  /** 
   * Contract the matrix element with the \f$\rho\f$ matrix of the 
   * incoming particle. The spins of the decay products are summed over.
   * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
   */
  Complex contract(const RhoDMatrix & rhoin) const;

  /** 
   * Contract the matrix element with the \f$\rho\f$ matrix of the 
   * incoming particle. The spins of the decay products are summed over.
   * @param con The conjugate matrix elemetn for the contraction
   * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
   */
  Complex contract(const TwoBodyDecayMatrixElement & con, 
  		   const RhoDMatrix & rhoin);
  //@}

public:
  
  /** 
   * Access to the individual helicity components. 
   */
  //@{
  /** 
   * Get the helicity components for a two body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   */
  Complex   operator () (unsigned int inhel, unsigned int outhel1,
			 unsigned int outhel2) const {
    return matrixElement_[inhel][outhel1][outhel2];
  }

  /** 
   * Set the helicity components for a two body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   */
  Complex & operator () (unsigned int inhel, unsigned int outhel1,
			 unsigned int outhel2) {
    return matrixElement_[inhel][outhel1][outhel2];
  }

  /** 
   * Get the helicity components for a three body decay
   */
  Complex   operator () (unsigned int , unsigned int ,
			 unsigned int , unsigned int ) const {
    assert(false);
    static const Complex out = 0.;
    return out;
  }

  /** 
   * Set the helicity components for a three body decay
   */
  Complex & operator () (unsigned int , unsigned int ,
			 unsigned int , unsigned int ) {
    assert(false);
    static Complex out = 0.;
    return out;
  }

  /** 
   * Get the helicity components for a four body decay
   */
  Complex   operator () (unsigned int , unsigned int ,
			 unsigned int , unsigned int ,
			 unsigned int ) const {
    assert(false);
    static const Complex out = 0.;
    return out;
  }

  /** 
   * Set the helicity components for a four body decay
   */
  Complex & operator () (unsigned int , unsigned int ,
			 unsigned int , unsigned int ,
			 unsigned int ) {
    assert(false);
    static Complex out = 0.;
    return out;
  }

  /** 
   * Get the helicity components for a five body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   */
  Complex   operator () (unsigned int , unsigned int ,
			 unsigned int , unsigned int ,
			 unsigned int , unsigned int ) const {
    assert(false);
    static const Complex out = 0.;
    return out;
  }

  /** 
   * Set the helicity components for a five body decay
   */
  Complex & operator () (unsigned int , unsigned int ,
			 unsigned int , unsigned int ,
			 unsigned int , unsigned int ) {
    assert(false);
    static Complex out = 0.;
    return out;
  }

  /** 
   * Get the helicity components for a six body decay
   */
  Complex   operator () (unsigned int ,unsigned int ,
			 unsigned int ,unsigned int ,
			 unsigned int ,unsigned int ,
			 unsigned int ) const {
    assert(false);
    static const Complex out = 0.;
    return out;
  }

  /** 
   * Set the helicity components for a six body decay
   */
  Complex & operator () (unsigned int , unsigned int ,
			 unsigned int , unsigned int ,
			 unsigned int , unsigned int ,
			 unsigned int ) {
    assert(false);
    static Complex out = 0.;
    return out;
  }

  /**
   * Get the helicity components for an \f$n\f$-body decay.
   * @param in The helicities of the external particles.
   */
  Complex   operator () (const vector<unsigned int> & in) const {
    assert(in.size()==3);
    return matrixElement_[in[0]][in[1]][in[2]];
  }

  /**
   * Set the helicity components for an \f$n\f$-body decay.
   * @param in The helicities of the external particles.
   */
  Complex & operator () (const vector<unsigned int> & in) {
    assert(in.size()==3);
    return matrixElement_[in[0]][in[1]][in[2]];
  }
  //@}

  /**
   *  Member to zero all the elements for the matrix element
   */
  void zero() {
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy) {
	for(unsigned int iz=0;iz<5;++iz) {
	  matrixElement_[ix][iy][iz]=0.;
	}
      }
    }
  }

private:

  /**
   * Storage of the matrix element, a vector is better for memory usage.
   */
  Complex matrixElement_[5][5][5];

};

}

#endif /* HERWIG_TwoBodyDecayMatrixElement_H */
