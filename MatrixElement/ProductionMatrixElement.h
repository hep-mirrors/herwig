// -*- C++ -*-
//
// ProductionMatrixElement.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ProductionMatrixElement_H
#define HERWIG_ProductionMatrixElement_H
//
// This is the declaration of the ProductionMatrixElement class.

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
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out);

  /**
   * Constructor for 2-2 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
			  PDT::Spin out2);
  
  /**
   * Constructor for 2-3 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,PDT::Spin out1,
			  PDT::Spin out2,PDT::Spin out3);

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
			  PDT::Spin out2,PDT::Spin out3, PDT::Spin out4);

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
			  PDT::Spin out5);

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
			  PDT::Spin out5, PDT::Spin out6);
  
  /**
   * Constructor for 2-n scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out A vector containing \f$2S+1\f$ for the outgoing particles.
   */
  ProductionMatrixElement(PDT::Spin in1,PDT::Spin in2,vector<PDT::Spin> out);
  
  /**
   * Default constructor.
   */
  ProductionMatrixElement() {}
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
			 unsigned int out) const;
  
  /**
   * Access the helicity components for a 2-1 scattering. This method supplies
   * the component and allows it to be changed.
   * @param in1 The helicity of the first incoming particle.
   * @param in2 The helicity of the second incoming particle.
   * @param out The helicity of the outgoing particle.
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (unsigned int in1,unsigned int in2,
			 unsigned int out);

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
			 unsigned int out1,unsigned int out2) const;

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
			 unsigned int out1,unsigned int out2);

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
			 unsigned int out3) const;

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
			 unsigned int out3);

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
			 unsigned int out3,unsigned int out4) const;
  
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
			 unsigned int out3, unsigned int out4);

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
			 unsigned int out5) const;
  
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
			 unsigned int out5);

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
			 unsigned int out5,unsigned int out6) const;

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
			 unsigned int out5, unsigned int out6);

  /**
   * Access the helicity components for a 2-n scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param hel The helicities of the incoming and outgoing particles
   * @return The matrix element for the given helicities.
   */
  Complex   operator () (vector<unsigned int> hel) const;
  
  /**
   * Access the helicity components for a 2-n scattering. This method supplies
   * the component and allows it to be changed.
   * @param hel The helicities of the incoming and outgoing particles
   * @return The matrix element for the given helicities.
   */
  Complex & operator () (vector<unsigned int> hel);
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
  void reset(const ProductionMatrixElement & x) const;
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
private:
  
  /**
   * Set the size of the vector containing the matrix element.
   */
  void setMESize();
  
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

#endif /* HERWIG_ProductionMatrixElement_H */
