// -*- C++ -*-
//
// GeneralDecayMatrixElement.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GeneralDecayMatrixElement_H
#define HERWIG_GeneralDecayMatrixElement_H
//
// This is the declaration of the GeneralDecayMatrixElement class.

#include "DecayMatrixElement.h"
#include "GeneralDecayMatrixElement.fh"

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
class GeneralDecayMatrixElement: public DecayMatrixElement {
  
public:
      
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  GeneralDecayMatrixElement() : DecayMatrixElement(999) {}

  /** 
   * Constructor for two body decay.
   * @param spinin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   */
  GeneralDecayMatrixElement(PDT::Spin spinin,
			    PDT::Spin outspin1,PDT::Spin outspin2)
    : DecayMatrixElement(2,spinin) {
    outspin().push_back(outspin1);
    outspin().push_back(outspin2);
    setMESize();
  }

  /** 
   * Constructor for three body decay. 
   * @param spinin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   */
  GeneralDecayMatrixElement(PDT::Spin spinin,PDT::Spin outspin1,
			    PDT::Spin outspin2,PDT::Spin outspin3) 
    : DecayMatrixElement(3,spinin) {
    outspin().push_back(outspin1);
    outspin().push_back(outspin2);
    outspin().push_back(outspin3);
    setMESize();
  }

  /** 
   * Constructor for four body decay.
   * @param spinin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   */
  GeneralDecayMatrixElement(PDT::Spin spinin,PDT::Spin outspin1,PDT::Spin outspin2,
			    PDT::Spin outspin3,PDT::Spin outspin4)
    : DecayMatrixElement(4,spinin) {
    outspin().push_back(outspin1);
    outspin().push_back(outspin2);
    outspin().push_back(outspin3);
    outspin().push_back(outspin4);
    setMESize();
  }

  /**
   * Constructor for five body decay.
   * @param spinin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   * @param outspin5 \f$2S+1\f$ for the fifth  decay product.
   */
  GeneralDecayMatrixElement(PDT::Spin spinin,PDT::Spin outspin1,PDT::Spin outspin2,
			    PDT::Spin outspin3,PDT::Spin outspin4,PDT::Spin outspin5)
    : DecayMatrixElement(5,spinin) {
    outspin().push_back(outspin1);
    outspin().push_back(outspin2);
    outspin().push_back(outspin3);
    outspin().push_back(outspin4);
    outspin().push_back(outspin5);
    setMESize();
  }

  /** 
   * Constructor for six body decay.
   * @param spinin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   * @param outspin5 \f$2S+1\f$ for the fifth  decay product.
   * @param outspin6 \f$2S+1\f$ for the sixth  decay product.
   */
  GeneralDecayMatrixElement(PDT::Spin spinin,PDT::Spin outspin1,PDT::Spin outspin2,
			    PDT::Spin outspin3,PDT::Spin outspin4,PDT::Spin outspin5,
			    PDT::Spin outspin6) 
    : DecayMatrixElement(6,spinin)  {
    outspin().push_back(outspin1);
    outspin().push_back(outspin2);
    outspin().push_back(outspin3);
    outspin().push_back(outspin4);
    outspin().push_back(outspin5);
    outspin().push_back(outspin6);
    setMESize();
  }

  /** 
   * Constructor for arbitary body decay.
   * @param spinin \f$2S+1\f$ for the decaying particle
   * @param spinout \f$2S+1\f$ for the decay products.
   */
  GeneralDecayMatrixElement(PDT::Spin spinin,vector<PDT::Spin> spinout)
    : DecayMatrixElement(spinout.size(),spinin) {
    outspin() = spinout;
    setMESize();
  }
  
  /** 
   * Constructor for arbitary body decay.
   * @param extspin  \f$2S+1\f$ external particles.
   */
  GeneralDecayMatrixElement(vector<PDT::Spin> extspin)
    : DecayMatrixElement(int(extspin.size())-1,extspin[0]) {
    outspin() = vector<PDT::Spin>(extspin.begin()+1,extspin.end());
    setMESize();
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
  Complex contract(const GeneralDecayMatrixElement & con, 
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
    unsigned int iloc = inhel*constants_[1]+
      outhel1*constants_[2]+outhel2*constants_[3];
    assert(outspin().size()==2&&iloc<matrixElement_.size());
    return matrixElement_[iloc];
  }

  /** 
   * Set the helicity components for a two body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   */
  Complex & operator () (unsigned int inhel,unsigned int outhel1,
			 unsigned int outhel2) {
    unsigned int iloc = inhel*constants_[1]+
      outhel1*constants_[2]+outhel2*constants_[3];
    assert(outspin().size()==2&&iloc<matrixElement_.size());
    return matrixElement_[iloc];
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
    unsigned int iloc = inhel*constants_[1]+outhel1*constants_[2]+
      outhel2*constants_[3]+outhel3*constants_[4];
    assert(outspin().size()==3&&iloc<matrixElement_.size());
    return matrixElement_[iloc];
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
    unsigned int iloc = inhel*constants_[1]+outhel1*constants_[2]+
      outhel2*constants_[3]+outhel3*constants_[4];
    assert(outspin().size()==3&&iloc<matrixElement_.size());
    return matrixElement_[iloc];
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
    assert(in.size()==outspin().size()+1);
    // contribution  for the incoming particle
    unsigned int iloc(in[0]*constants_[1]);
    // contributions for the outgoing particles
    for(unsigned int ix=1;ix<in.size();++ix){iloc+=in[ix]*constants_[ix+1];}
    assert(iloc<matrixElement_.size());
    return matrixElement_[iloc];
  }

  /**
   * Set the helicity components for an \f$n\f$-body decay.
   * @param in The helicities of the external particles.
   */
  Complex & operator () (const vector<unsigned int> & in) {
    assert(in.size()==outspin().size()+1);
    // contribution  for the incoming particle
    unsigned int iloc(in[0]*constants_[1]);
    // contributions for the outgoing particles
    for(unsigned int ix=1;ix<in.size();++ix){iloc+=in[ix]*constants_[ix+1];}
    assert(iloc<matrixElement_.size());
    return matrixElement_[iloc];
  }
  //@}

  /**
   *  Member to zero all the elements for the matrix element
   */
  void zero() {
    for(unsigned int ix=0;ix<matrixElement_.size();++ix)
      matrixElement_[ix]=0.;
  }
  
  /**
   * Set the size of the vector containing the matrix element.
   */
  void setMESize() {
    int isize = inspin();
    for(unsigned int ix=0;ix<outspin().size();++ix) isize*=outspin()[ix];
    matrixElement_.resize(isize,0.);
    // set up the constants for the mapping of helicity to vector index
    constants_.resize(outspin().size()+2);
    int temp=1;
    for(unsigned int ix=outspin().size();ix>0;--ix) {
      temp*=outspin()[ix-1];constants_[ix]=temp;
    }
    temp *= inspin();
    constants_[0]=temp;
    constants_[outspin().size()+1]=1;
  }
  
private:

  /**
   * Storage of the matrix element, a vector is better for memory usage.
   */
  mutable vector<Complex> matrixElement_;

  /**
   * Constants needed to map the index of the vector to a helicity structure.
   */
  mutable vector<unsigned int> constants_;

};

}

#endif /* HERWIG_GeneralDecayMatrixElement_H */
