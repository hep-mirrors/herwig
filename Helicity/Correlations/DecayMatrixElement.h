// -*- C++ -*-
#ifndef HERWIG_DecayMatrixElement_H
#define HERWIG_DecayMatrixElement_H
//
// This is the declaration of the DecayMatrixElement class.

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>
#include <ThePEG/Helicity/RhoDMatrix.h>
// #include "DecayMatrixElement.fh"
// #include "DecayMatrixElement.xh"

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;

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
  inline DecayMatrixElement();

  /** 
   * Constructor for two body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   */
  inline DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2);

  /** 
   * Constructor for three body decay. 
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   */
  inline DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,
			    PDT::Spin outspin2,PDT::Spin outspin3);

  /** 
   * Constructor for four body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   */
  inline DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2,
			    PDT::Spin outspin3,PDT::Spin outspin4);

  /**
   * Constructor for five body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin1 \f$2S+1\f$ for the first  decay product.
   * @param outspin2 \f$2S+1\f$ for the second decay product.
   * @param outspin3 \f$2S+1\f$ for the third  decay product.
   * @param outspin4 \f$2S+1\f$ for the fourth decay product.
   * @param outspin5 \f$2S+1\f$ for the fifth  decay product.
   */
  inline DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2,
			    PDT::Spin outspin3,PDT::Spin outspin4,PDT::Spin outspin5);

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
  inline DecayMatrixElement(PDT::Spin inspin,PDT::Spin outspin1,PDT::Spin outspin2,
			    PDT::Spin outspin3,PDT::Spin outspin4,PDT::Spin outspin5,
			    PDT::Spin outspin6);

  /** 
   * Constructor for arbitray body decay.
   * @param inspin \f$2S+1\f$ for the decaying particle
   * @param outspin \f$2S+1\f$ for the decay products.
   */
  inline DecayMatrixElement(PDT::Spin inspin,vector<PDT::Spin> outspin);

  /** 
   * Constructor for arbitray body decay.
   * @param extspin  \f$2S+1\f$ external particles.
   */
  inline DecayMatrixElement(vector<PDT::Spin> extspin);

  /**
   * Copy-constructor.
   */
  inline DecayMatrixElement(const DecayMatrixElement &);

  /**
   * Destructor.
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
  inline PDT::Spin inspin();

  /** 
   * Get the spins of the outgoing particles.
   */
  inline vector<PDT::Spin> outspin();
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
  RhoDMatrix calculateDMatrix(vector<RhoDMatrix> rhoout);

  /** 
   * Calculate the \f$\rho\f$ matrix for a given outgoing particle.
   * @param ipart The outgoing particle the \f$\rho\f$ matrix is needed for
   * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
   * @param rhoout he \f$D\f$ matrices for the other decay products.
   */
  RhoDMatrix calculateRhoMatrix(int ipart,RhoDMatrix rhoin,vector<RhoDMatrix> rhoout);

  /** 
   * Contract the matrix element with the \f$\rho\f$ matrix of the 
   * incoming particle. The spins of the decay products are summed over.
   * @param rhoin The \f$\rho\f$ matrix for the decaying particle.
   */
  Complex contract(RhoDMatrix & rhoin);
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
  inline Complex   operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2) const;

  /** 
   * Set the helicity components for a two body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   */
  inline Complex & operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2);

  /** 
   * Get the helicity components for a three body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   */
  inline Complex   operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3) const;

  /** 
   * Set the helicity components for a three body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   */
  inline Complex & operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3);

  /** 
   * Get the helicity components for a four body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   */
  inline Complex   operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3,
				unsigned int outhel4) const;

  /** 
   * Set the helicity components for a four body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   */
  inline Complex & operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3,
				unsigned int outhel4);

  /** 
   * Get the helicity components for a five body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   */
  inline Complex   operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3,
				unsigned int outhel4,unsigned int outhel5) const;

  /** 
   * Set the helicity components for a five body decay
   * @param inhel The helicity of the decaying particle.
   * @param outhel1 The helicity of the first  decay product.
   * @param outhel2 The helicity of the second decay product.
   * @param outhel3 The helicity of the third  decay product.
   * @param outhel4 The helicity of the fourth decay product.
   * @param outhel5 The helicity of the fifth  decay product.
   */
  inline Complex & operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3,
				unsigned int outhel4,unsigned int outhel5);

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
  inline Complex   operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3,
				unsigned int outhel4,unsigned int outhel5,
				unsigned int outhel6) const;

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
  inline Complex & operator () (unsigned int inhel,unsigned int outhel1,
				unsigned int outhel2,unsigned int outhel3,
				unsigned int outhel4,unsigned int outhel5,
				unsigned int outhel6);

  /**
   * Get the helicity components for an \f$n\f$-body decay.
   * @param exthel The helicities of the external particles.
   */
  inline Complex   operator () (vector<unsigned int> exthel) const;

  /**
   * Set the helicity components for an \f$n\f$-body decay.
   * @param exthel The helicities of the external particles.
   */
  inline Complex & operator () (vector<unsigned int> exthel);
  //@}
  
public:
  
  /**
   * Reset the matrix element. 
   */
  inline void reset(const DecayMatrixElement &) const;
  
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
  inline void setMESize();
  
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
}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of DecayMatrixElement.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::DecayMatrixElement,1> {
    /** Typedef of the base class of DecayMatrixElement. */
  typedef Base NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::DecayMatrixElement>
  : public ClassTraitsBase<Herwig::Helicity::DecayMatrixElement> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::Helicity::DecayMatrixElement"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwCorrelations.so"; }

};

}

#include "DecayMatrixElement.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayMatrixElement.tcc"
#endif

#endif /* HERWIG_DecayMatrixElement_H */
