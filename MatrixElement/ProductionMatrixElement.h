// -*- C++ -*-
//
// ProductionMatrixElement.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ProductionMatrixElement_H
#define HERWIG_ProductionMatrixElement_H
//
// This is the declaration of the ProductionMatrixElement class.

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>

namespace Herwig {


using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The storage of the helicity amplitude expression for the matrix element 
 *  of a hard process. Two incoming particles and an arbitary number of 
 *  external particles are supported.
 *
 *  @see DecayMatrixElement
 *  @see RhoDMatrix
 *  @see HardVertex
 *
 *  \author Peter Richardson
 */

class ProductionMatrixElement {  
      
public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor for 2-1 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out \f$2S+1\f$ for the outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out) {
    _nout=2;
    _inspin.resize(2);
    _inspin[0]=in1;
    _inspin[1]=in2;
    _outspin.push_back(out);
    setMESize();
  }

  /**
   * Constructor for 2-2 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
			  PDT::Spin out2) {
    _nout=2;
    _inspin.resize(2);
    _inspin[0]=in1; 
    _inspin[1]=in2;
    _outspin.push_back(out1);
    _outspin.push_back(out2);
    setMESize();
  }
  
  /**
   * Constructor for 2-3 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
			  PDT::Spin out2,PDT::Spin out3) {
    _inspin.resize(2);
    _nout=3;
    _inspin[0]=in1;
    _inspin[1]=in2;
    _outspin.push_back(out1);
    _outspin.push_back(out2);
    _outspin.push_back(out3);
    setMESize();
  }

  /**
   * Constructor for 2-4 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   * @param out4 \f$2S+1\f$ for the fourth outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
			  PDT::Spin out2,PDT::Spin out3, PDT::Spin out4) {
    _nout=4;
    _inspin.resize(2);
    _inspin[0]=in1;
    _inspin[1]=in2;
    _outspin.push_back(out1);
    _outspin.push_back(out2);
    _outspin.push_back(out3);
    _outspin.push_back(out4);
    setMESize();
  }

  /**
   * Constructor for 2-5 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   * @param out4 \f$2S+1\f$ for the fourth outgoing particle.
   * @param out5 \f$2S+1\f$ for the fifth outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
			  PDT::Spin out2,PDT::Spin out3, PDT::Spin out4,
			  PDT::Spin out5) {
    _nout=5;
    _inspin.resize(2);
    _inspin[0]=in1;
    _inspin[1]=in2;
    _outspin.push_back(out1);
    _outspin.push_back(out2);
    _outspin.push_back(out3);
    _outspin.push_back(out4);
    _outspin.push_back(out5);
    setMESize();
  }

  /**
   * Constructor for 2-6 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   * @param out4 \f$2S+1\f$ for the fourth outgoing particle.
   * @param out5 \f$2S+1\f$ for the fifth outgoing particle.
   * @param out6 \f$2S+1\f$ for the sixth outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
			  PDT::Spin out2,PDT::Spin out3, PDT::Spin out4,
			  PDT::Spin out5, PDT::Spin out6) {
    _nout=6;
    _inspin.resize(2);
    _inspin[0]=in1;
    _inspin[1]=in2;
    _outspin.push_back(out1);
    _outspin.push_back(out2);
    _outspin.push_back(out3);
    _outspin.push_back(out4);
    _outspin.push_back(out5);
    _outspin.push_back(out6);
    setMESize();
  }
  
  /**
   * Constructor for 2-n scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out A vector containing \f$2S+1\f$ for the outgoing particles.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,vector<PDT::Spin> out) {
    _inspin.resize(2);
    _nout=out.size(); 
    _inspin[0]=in1;
    _inspin[1]=in2;
    _outspin=out;
    setMESize();
  }
  
  /**
   * Default constructor.
   */
  ProductionMatrixElement() {};
  //@}

public:
     
  /** @name Access to the spins of the particles. */
  //@{
  /**
   * Get the spins of the incoming particles particle
   * @return A vector containing \f$2S+1\f$ for the two incoming particles.
   */
  vector<PDT::Spin> inspin() {return _inspin;}

  /**
   * Get the spins of the outgoing particles.
   * @return A vector containing \f$2S+1\f$ for the outgoing particles.
   */
  vector<PDT::Spin> outspin() {return _outspin;}
  //@}  

public:

  /** @name Access to the individual helicity components. */
  //@{
  /**
   * Access the helicity components for a 2-1 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out The helicity of the outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (unsigned int in1,unsigned int in2,
			 unsigned int out) const {
    assert(_outspin.size()==1);
    unsigned int iloc = in1*_constants[1] + in2*_constants[2] + out*_constants[3];
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }
  
  /**
   * Access the helicity components for a 2-1 scattering. This method supplies
   * the component and allows it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out The helicity of the outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (unsigned int in1,unsigned int in2,
			 unsigned int out) {
    assert(_outspin.size()==1);
    unsigned int iloc = in1*_constants[1] + in2*_constants[2] + out*_constants[3];
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /**
   * Access the helicity components for a 2-2 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2) const {
    assert(_outspin.size()==2);
    unsigned int iloc = in1*_constants[1] + in2*_constants[2] +
      out1*_constants[3] + out2*_constants[4];
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /**
   * Access the helicity components for a 2-2 scattering. This method supplies
   * the component and allows it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2) {
    assert(_outspin.size()==2);
    unsigned int iloc = in1*_constants[1] + in2*_constants[2] +
      out1*_constants[3] + out2*_constants[4];
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }

  /**
   * Access the helicity components for a 2-3 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3) const {
    assert(_outspin.size()==3);
    vector<unsigned int> ivec(5);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    return (*this)(ivec);
  }

  /**
   * Access the helicity components for a 2-3 scattering. This method supplies
   * the component and allows it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3) {
    assert(_outspin.size()==3);
    vector<unsigned int> ivec(5);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    return (*this)(ivec);
  }

  /**
   * Access the helicity components for a 2-4 scattering.  This method supplies
   * the component but does not allow it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @param out4 The helicity of the fourth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3,unsigned int out4) const {
    assert(_outspin.size()==4);
    vector<unsigned int> ivec(6);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    ivec[5]=out4;
    return (*this)(ivec);
  }
  
  /**
   * Access the helicity components for a 2-4 scattering. This method supplies
   * the component and allows it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @param out4 The helicity of the fourth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3, unsigned int out4) {
    assert(_outspin.size()==4);
    vector<unsigned int> ivec(6);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    ivec[5]=out4;
    return (*this)(ivec);
  }

  /**
   * Access the helicity components for a 2-5 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @param out4 The helicity of the fourth outgoing particle.
   * @param out5 The helicity of the fifth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3,unsigned int out4,
			 unsigned int out5) const {
    assert(_outspin.size()==5);
    vector<unsigned int> ivec(7);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    ivec[5]=out4;
    ivec[6]=out5;
    return (*this)(ivec);
  }
  
  /**
   * Access the helicity components for a 2-5 scattering. This method supplies
   * the component and allows it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @param out4 The helicity of the fourth outgoing particle.
   * @param out5 The helicity of the fifth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3, unsigned int out4,
			 unsigned int out5) {
    assert(_outspin.size()==5);
    vector<unsigned int> ivec(7);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    ivec[5]=out4;
    ivec[6]=out5;
    return (*this)(ivec);
  }

  /**
   * Access the helicity components for a 2-6 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @param out4 The helicity of the fourth outgoing particle.
   * @param out5 The helicity of the fifth outgoing particle.
   * @param out6 The helicity of the sixth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3,unsigned int out4,
			 unsigned int out5,unsigned int out6) const {
    assert(_outspin.size()==6);
    vector<unsigned int> ivec(8);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    ivec[5]=out4;
    ivec[6]=out5;
    ivec[7]=out6;
    return (*this)(ivec);
  }

  /**
   * Access the helicity components for a 2-6 scattering. This method supplies
   * the component and allows it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out1 The helicity of the first outgoing particle.
   * @param out2 The helicity of the second outgoing particle.
   * @param out3 The helicity of the third outgoing particle.
   * @param out4 The helicity of the fourth outgoing particle.
   * @param out5 The helicity of the fifth outgoing particle.
   * @param out6 The helicity of the sixth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (unsigned int in1,unsigned int in2,
			 unsigned int out1,unsigned int out2,
			 unsigned int out3, unsigned int out4,
			 unsigned int out5, unsigned int out6) {
    assert(_outspin.size()==6);
    vector<unsigned int> ivec(8);
    ivec[0]=in1;
    ivec[1]=in2;
    ivec[2]=out1;
    ivec[3]=out2;
    ivec[4]=out3;
    ivec[5]=out4;
    ivec[6]=out5;
    ivec[7]=out6;
    return (*this)(ivec);
  }

  /**
   * Access the helicity components for a 2-n scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param hel The helicities of the incoming and outgoing particles
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (vector<unsigned int> hel) const {
    assert(_outspin.size() == hel.size()-2);
    unsigned int iloc(0),ix;
    // incoming and outgoing particles
    for(ix=0;ix<hel.size();++ix)
      iloc += hel[ix]*_constants[ix+1];
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }
  
  /**
   * Access the helicity components for a 2-n scattering. This method supplies
   * the component and allows it to be changed.
   * @param hel The helicities of the incoming and outgoing particles
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (vector<unsigned int> hel) {
    assert(_outspin.size() == hel.size()-2);
    unsigned int iloc=0,ix;
    // incoming particles
    for(ix=0;ix<hel.size();++ix)
      iloc += hel[ix]*_constants[ix+1];
    assert(iloc<_matrixelement.size());
    return _matrixelement[iloc];
  }
  //@}

public:

  /**
   * Calculate the decay matrix for an incoming particle.
   */
  RhoDMatrix calculateDMatrix(int,const RhoDMatrix &,
			      const vector<RhoDMatrix> &) const;

  /**
   * Calculate the rho matrix for a given outgoing particle.
   */
  RhoDMatrix calculateRhoMatrix(int,const RhoDMatrix &,
                                const RhoDMatrix &,
				const vector<RhoDMatrix> &) const;

  /**
   *  Compute the spin averaged matrix element
   */
  double average() const;

  /**
   *  Compute the spin average matrix element
   */
  double average(const RhoDMatrix & in1, 
		 const RhoDMatrix & in2) const;

  /**
   *  Compute the spin average matrix element
   */
  Complex average(const ProductionMatrixElement & me2,
		  const RhoDMatrix & in1, 
		  const RhoDMatrix & in2) const;
  
public:

  /**
   * Reset the matrix element.
   */
  void reset(const ProductionMatrixElement & x) const {
    _nout = x._nout;
    _inspin = x._inspin;
    _outspin = x._outspin;
    _matrixelement = x._matrixelement;
    _constants     = x._constants;
  }
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
private:
  
  /**
   * Set the size of the vector containing the matrix element.
   */
  void setMESize() {
    unsigned int ix;
    int isize=_inspin[0]*_inspin[1];
    for(ix=0;ix<_outspin.size();++ix)
      isize*=_outspin[ix];
    // zero the matrix element
    _matrixelement.resize(isize,0.);
    // set up the constants for the mapping of helicity to vectro index
    _constants.resize(_outspin.size()+3);
    unsigned int temp=1;
    for(ix=_outspin.size()+1;ix>1;--ix) {
      temp*=_outspin[ix-2];
      _constants[ix]=temp;
    }
    temp*=_inspin[1];_constants[1]=temp;
    temp*=_inspin[0];_constants[0]=temp;
    _constants[_outspin.size()+2]=1;
  }
  
private:
  
  /**
   * Number of outgoing particles.
   */
  mutable unsigned int _nout;

  /**
   * Spin of the incoming particles as 2s+1.
   */
  mutable vector<PDT::Spin> _inspin;

  /**
   * Spins of the outgoing particles.
   */
  mutable vector<PDT::Spin> _outspin;

  /**
   * Storage of the matrix element, a vector is better for memory usage.
   */
  mutable vector<Complex> _matrixelement;

  /**
   * Constants needed to map the index of the vector to a helicity structure.
   */
  mutable vector<int> _constants;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ProductionMatrixElement.
 */
template <>
struct BaseClassTrait<Herwig::ProductionMatrixElement,1> {
  /** Typedef of the base class of ProductionMatrixElement. */
  typedef Base NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ProductionMatrixElement>
  : public ClassTraitsBase<Herwig::ProductionMatrixElement> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::ProductionMatrixElement"; }
};

/** @endcond */
  
}

#endif /* HERWIG_ProductionMatrixElement_H */
