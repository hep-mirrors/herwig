// -*- C++ -*-
//
// DecayMatrixElement.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DecayMatrixElement_H
#define HERWIG_DecayMatrixElement_H
//
// This is the declaration of the DecayMatrixElement class.

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>

namespace Herwig {

using namespace ThePEG;


/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  Implementation of the complex matrix element for a decay.
 *  An arbitary number of external particles are supported.
 *
 *  @see RhoDMatrix
 *  @see ProductionMatrixElement
 *  @see DecayVertex
 */
class DecayMatrixElement: public Base {
  
public:
      
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  DecayMatrixElement() : _nout(999) {}

  /** 
   * Constructor for two body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   */
  DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2)
    : _nout(2), _inspin(inspin), _outspin(2) {
    _outspin[0] = outspin1;
    _outspin[1] = outspin2;
    setMESize();
  }

  /** 
   * Constructor for three body decay. 
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   */
  DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,
		     PDT::Spin outspin2,PDT::Spin outspin3) 
    : _nout(3), _inspin(inspin), _outspin(3) {
    _outspin[0] = outspin1;
    _outspin[1] = outspin2;
    _outspin[2] = outspin3;
    setMESize();
  }

  /** 
   * Constructor for four body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   */
  DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2,
		     PDT::Spin outspin3,PDT::Spin outspin4)
    : _nout(4), _inspin(inspin), _outspin(4) {
    _outspin[0] = outspin1;
    _outspin[1] = outspin2;
    _outspin[2] = outspin3;
    _outspin[3] = outspin4;
    setMESize();
  }

  /**
   * Constructor for five body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   * @param outspin5 \f$2S+1\f$ for the fifth  decay product.
   */
  DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2,
		     PDT::Spin outspin3,PDT::Spin outspin4,PDT::Spin outspin5)
    : _nout(5), _inspin(inspin), _outspin(5) {
    _outspin[0] = outspin1;
    _outspin[1] = outspin2;
    _outspin[2] = outspin3;
    _outspin[3] = outspin4;
    _outspin[4] = outspin5;
    setMESize();
  }

  /** 
   * Constructor for six body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   * @param outspin5 \f$2S+1\f$ for the fifth  decay product.
   * @param outspin6 \f$2S+1\f$ for the sixth  decay product.
   */
  DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2,
		     PDT::Spin outspin3,PDT::Spin outspin4,PDT::Spin outspin5,
		     PDT::Spin outspin6) 
    : _nout(6), _inspin(inspin), _outspin(6) {
    _outspin[0] = outspin1;
    _outspin[1] = outspin2;
    _outspin[2] = outspin3;
    _outspin[3] = outspin4;
    _outspin[4] = outspin5;
    _outspin[5] = outspin6;
    setMESize();
  }

  /** 
   * Constructor for arbitary body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin \f$2S+1\f$ for the decay products.
   */
  DecayMatrixElement(PDT::Spin inspin,vector<PDT::Spin> outspin)
    : _nout(outspin.size()), _inspin(inspin), _outspin(outspin) {
    setMESize();
  }

  /** 
   * Constructor for arbitary body decay.
   * @param extspin  \f$2S+1\f$ external particles.
   */
  DecayMatrixElement(vector<PDT::Spin> extspin)
    : _nout(extspin.size()-1), _inspin(extspin[0]),
      _outspin(extspin.begin()+1,extspin.end()) {
    setMESize();
  }
  //@}  

public: 
      
  /**
   * Access to the spins of the particles
   */
  //@{
  /** 
   * Get the spin of the incoming particle.
   */
  PDT::Spin inspin() {return _inspin;}

  /** 
   * Get the spins of the outgoing particles.
   */
  vector<PDT::Spin> outspin() {return _outspin;}
  //@}

  /**
   *  Get the number of outgoing particles
   */
  unsigned int nOut() const {return _nout;}

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
  Complex contract(const DecayMatrixElement & con, 
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
  Complex   operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2) const {
    unsigned int iloc = inhel*_constants[1]+
      outhel1*_constants[2]+outhel2*_constants[3];
    assert(_outspin.size()==2&&iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /** 
   * Set the helicity components for a two body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   */
  Complex & operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2) {
    unsigned int iloc = inhel*_constants[1]+
      outhel1*_constants[2]+outhel2*_constants[3];
    assert(_outspin.size()==2&&iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /** 
   * Get the helicity components for a three body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   */
  Complex   operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3) const {
    unsigned int iloc = inhel*_constants[1]+outhel1*_constants[2]+
      outhel2*_constants[3]+outhel3*_constants[4];
    assert(_outspin.size()==3&&iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /** 
   * Set the helicity components for a three body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   */
  Complex & operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3) {
    unsigned int iloc = inhel*_constants[1]+outhel1*_constants[2]+
      outhel2*_constants[3]+outhel3*_constants[4];
    assert(_outspin.size()==3&&iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /** 
   * Get the helicity components for a four body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   */
  Complex   operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3,
			 unsigned int outhel4) const {
    vector<unsigned int> itemp(5);
    itemp[0]=inhel  ; itemp[1]=outhel1;
    itemp[2]=outhel2; itemp[3]=outhel3;
    itemp[4]=outhel4 ;
    return (*this)(itemp);
  }

  /** 
   * Set the helicity components for a four body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   */
  Complex & operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3,
			 unsigned int outhel4) {
    vector<unsigned int> itemp(5);
    itemp[0]=inhel  ; itemp[1]=outhel1;
    itemp[2]=outhel2; itemp[3]=outhel3;
    itemp[4]=outhel4; 
    return (*this)(itemp);
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
  Complex   operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3,
			 unsigned int outhel4,unsigned int outhel5) const {
    vector<unsigned int> itemp(6);
    itemp[0]=inhel  ; itemp[1]=outhel1;
    itemp[2]=outhel2; itemp[3]=outhel3;
    itemp[4]=outhel4 ;itemp[5]=outhel5;
    return (*this)(itemp);
  }

  /** 
   * Set the helicity components for a five body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   */
  Complex & operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3,
			 unsigned int outhel4,unsigned int outhel5) {
    vector<unsigned int> itemp(6);
    itemp[0]=inhel  ; itemp[1]=outhel1;
    itemp[2]=outhel2; itemp[3]=outhel3;
    itemp[4]=outhel4; itemp[5]=outhel5;
    return (*this)(itemp);
  }

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
  Complex   operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3,
			 unsigned int outhel4,unsigned int outhel5,
			 unsigned int outhel6) const {
    vector<unsigned int> itemp(7);
    itemp[0]=inhel  ; itemp[1]=outhel1;
    itemp[2]=outhel2; itemp[3]=outhel3;
    itemp[4]=outhel4 ;itemp[5]=outhel5;
    itemp[6]=outhel6;
    return (*this)(itemp);
  }

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
  Complex & operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2,unsigned int outhel3,
			 unsigned int outhel4,unsigned int outhel5,
			 unsigned int outhel6) {
    vector<unsigned int> itemp(7);
    itemp[0]=inhel  ; itemp[1]=outhel1;
    itemp[2]=outhel2; itemp[3]=outhel3;
    itemp[4]=outhel4; itemp[5]=outhel5;
    itemp[6]=outhel6;
    return (*this)(itemp);
  }

  /**
   * Get the helicity components for an \f$n\f$-body decay.
   * @param in The helicities of the external particles.
   */
  Complex   operator () (const vector<unsigned int> & in) const {
    assert(in.size()==_outspin.size()+1);
    // contribution  for the incoming particle
    unsigned int iloc(in[0]*_constants[1]);
    // contributions for the outgoing particles
    for(unsigned int ix=1;ix<in.size();++ix){iloc+=in[ix]*_constants[ix+1];}
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /**
   * Set the helicity components for an \f$n\f$-body decay.
   * @param in The helicities of the external particles.
   */
  Complex & operator () (const vector<unsigned int> & in) {
    assert(in.size()==_outspin.size()+1);
    // contribution  for the incoming particle
    unsigned int iloc(in[0]*_constants[1]);
    // contributions for the outgoing particles
    for(unsigned int ix=1;ix<in.size();++ix){iloc+=in[ix]*_constants[ix+1];}
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }
  //@}

  /**
   *  Member to zero all the elements for the matrix element
   */
  void zero() {
    for(unsigned int ix=0;ix<_matrixelement.size();++ix)
      _matrixelement[ix]=0.;
  }
  
public:
  
  /**
   * Reset the matrix element. 
   */
  void reset(const DecayMatrixElement & x) const {
    _nout = x._nout;
    _inspin = x._inspin;
    _outspin =x._outspin;
    _matrixelement=x._matrixelement;
    _constants=x._constants;
  }
  
public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
private:
  
  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<DecayMatrixElement> initDecayMatrixElement;
  
  /**
   * Private and non-existent assignment operator.
   */
  DecayMatrixElement & operator=(const DecayMatrixElement &);
  
private:
  
  /**
   * Set the size of the vector containing the matrix element.
   */
  void setMESize() {
    int isize=_inspin;
    for(unsigned int ix=0;ix<_outspin.size();++ix) isize*=_outspin[ix];
    _matrixelement.resize(isize,0.);
    // set up the constants for the mapping of helicity to vector index
    _constants.resize(_outspin.size()+2);
    int temp=1;
    for(unsigned int ix=_outspin.size();ix>0;--ix) {
      temp*=_outspin[ix-1];_constants[ix]=temp;
    }
    temp*=_inspin;_constants[0]=temp;
    _constants[_outspin.size()+1]=1;
  }
  
private:
  
  /**
   * Number of outgoing particles.
   */
  mutable unsigned int _nout;

  /**
   * Spin of the incoming particle as 2s+1.
   */
  mutable PDT::Spin _inspin;

  /**
   * Spins of the outgonig particles.
   */
  mutable vector<PDT::Spin> _outspin;

  /**
   * Storage of the matrix element, a vector is better for memory usage.
   */
  mutable vector<Complex> _matrixelement;

  /**
   * Constants needed to map the index of the vector to a helicity structure.
   */
  mutable vector<unsigned int> _constants;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of DecayMatrixElement.
 */
template <>
struct BaseClassTrait<Herwig::DecayMatrixElement,1> {
    /** Typedef of the base class of DecayMatrixElement. */
  typedef Base NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::DecayMatrixElement>
  : public ClassTraitsBase<Herwig::DecayMatrixElement> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::DecayMatrixElement"; }
};

/** @endcond */

}

#endif /* HERWIG_DecayMatrixElement_H */
