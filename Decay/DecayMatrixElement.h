// -*- C++ -*-
#ifndef Herwig_DecayMatrixElement_H
#define Herwig_DecayMatrixElement_H
//
// This is the declaration of the DecayMatrixElement class.
//

#include "ThePEG/Pointer/ReferenceCounted.h"
#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>
#include "DecayMatrixElement.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DecayMatrixElement class.
 */
class DecayMatrixElement: public Pointer::ReferenceCounted {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DecayMatrixElement(unsigned int ntemp=999, PDT::Spin spin = PDT::SpinUndefined) 
    : nOut_(ntemp), inSpin_(spin) {}

  /**
   * The destructor.
   */
  virtual ~DecayMatrixElement();
  //@}

public: 
      
  /**
   * Access to the spins of the particles
   */
  //@{
  /** 
   * Get the spin of the incoming particle.
   */
  PDT::Spin inspin() const {return inSpin_;}

  /** 
   * Get the spins of the outgoing particles.
   */
  const vector<PDT::Spin> & outspin() const {return outSpin_;}
  //@}

  /**
   *  Get the number of outgoing particles
   */
  unsigned int nOut() const {return nOut_;}

  /**
   * Spin Density matrices
   */
  //@{
  /** 
   * Calculate the decay matrix for this decay.
   * @param rhoout The \f$D\f$ matrix for this decay.
   */
  virtual RhoDMatrix calculateDMatrix(const vector<RhoDMatrix> & rhoout) const = 0;

  /** 
   * Calculate the \f$\rho\f$ matrix for a given outgoing particle.
   * @param ipart The outgoing particle the \f$\rho\f$ matrix is needed for
   * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
   * @param rhoout he \f$D\f$ matrices for the other decay products.
   */
  virtual RhoDMatrix calculateRhoMatrix(int ipart,const RhoDMatrix & rhoin,
					const vector<RhoDMatrix> & rhoout) const = 0;

  /** 
   * Contract the matrix element with the \f$\rho\f$ matrix of the 
   * incoming particle. The spins of the decay products are summed over.
   * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
   */
  virtual Complex contract(const RhoDMatrix & rhoin) const = 0;

  // /** 
  //  * Contract the matrix element with the \f$\rho\f$ matrix of the 
  //  * incoming particle. The spins of the decay products are summed over.
  //  * @param con The conjugate matrix elemetn for the contraction
  //  * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
  //  */
  // Complex contract(const DecayMatrixElement & con, 
  // 		   const RhoDMatrix & rhoin);
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
  virtual Complex   operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2) const = 0;
  
  /** 
   * Set the helicity components for a two body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   */
  virtual Complex & operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2) = 0;
  
  /** 
   * Get the helicity components for a three body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   */
  virtual Complex   operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3) const = 0;
  
  /** 
   * Set the helicity components for a three body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   */
  virtual Complex & operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3) = 0;
  
  /** 
   * Get the helicity components for a four body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   */
  virtual Complex   operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3,
				 unsigned int outhel4) const = 0;
  
  /** 
   * Set the helicity components for a four body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   */
  virtual Complex & operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3,
				 unsigned int outhel4) = 0;
  
  /** 
   * Get the helicity components for a five body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   */
  virtual Complex   operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3,
				 unsigned int outhel4,unsigned int outhel5) const = 0;
  
  /** 
   * Set the helicity components for a five body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   */
  virtual Complex & operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3,
				 unsigned int outhel4,unsigned int outhel5) = 0;
  
  /** 
   * Get the helicity components for a six body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   * @param outhel6 The helicity of the sixth  decay product.
   */
  virtual Complex   operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3,
				 unsigned int outhel4,unsigned int outhel5,
				 unsigned int outhel6) const = 0;
  
  /** 
   * Set the helicity components for a six body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   * @param outhel6 The helicity of the sixth  decay product.
   */
  virtual Complex & operator () (unsigned int inhel,unsigned int outhel1,
				 unsigned int outhel2,unsigned int outhel3,
				 unsigned int outhel4,unsigned int outhel5,
				 unsigned int outhel6) = 0;
  
  /**
   * Get the helicity components for an \f$n\f$-body decay.
   * @param in The helicities of the external particles.
   */
  virtual Complex   operator () (const vector<unsigned int> & in) const = 0;
  
  /**
   * Set the helicity components for an \f$n\f$-body decay.
   * @param in The helicities of the external particles.
   */
  virtual Complex & operator () (const vector<unsigned int> & in) = 0;
  //@}
  
  /**
   *  Member to zero all the elements for the matrix element
   */
  virtual void zero() = 0;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DecayMatrixElement & operator=(const DecayMatrixElement &) = delete;

protected:

  /** 
   * Get the spins of the outgoing particles.
   */
  vector<PDT::Spin> & outspin() {return outSpin_;}
  
private:
  
  /**
   * Number of outgoing particles.
   */
  unsigned int nOut_;

  /**
   * Spin of the incoming particle as 2s+1.
   */
  PDT::Spin inSpin_;

  /**
   * Spins of the outgoing particles.
   */
  vector<PDT::Spin> outSpin_;

};

}

#endif /* Herwig_DecayMatrixElement_H */
