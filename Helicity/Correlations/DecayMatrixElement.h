// -*- C++ -*-
#ifndef HERWIG_DecayMatrixElement_H
#define HERWIG_DecayMatrixElement_H
//
// This is the declaration of the DecayMatrixElement class.

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
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
 *  Implementation of the complex matrix element for a decay
 *  an arbitary number of external particles are supported.
 *
 *  @see RhoDMatrix
 *  @see ProductionMatrixElement
 *  @see DecayVertex
 */
class DecayMatrixElement: public Base {
  
public:
      
  inline DecayMatrixElement();

  /** 
   * Constructor for two body decay. 
   */
  inline DecayMatrixElement(int,int,int);

  /** 
   * Constructor for three body decay. 
   */
  inline DecayMatrixElement(int,int,int,int);

  /** 
   * Constructor for four body decay.
   */
  inline DecayMatrixElement(int,int,int,int,int);

  /**
   * Constructor for five body decay.
   */
  inline DecayMatrixElement(int,int,int,int,int,int);

  /** 
   * Constructor for six body decay.
   */
  inline DecayMatrixElement(int,int,int,int,int,int,int);

  /** 
   * Constructors for arbitray body decay.
   */
  inline DecayMatrixElement(int,vector<int>);
  inline DecayMatrixElement(vector<int>);

  /** 
   * Standard ctors and dtor.
   */
  inline DecayMatrixElement(const DecayMatrixElement &);
  virtual ~DecayMatrixElement();
  
public: 
      
  /** 
   * Get the spin of the incoming particle.
   */
  inline int inspin();

  /** 
   * Get the spins of the outgoing particles.
   */
  inline vector<int> outspin();

public:

  /** 
   * Calculate the decay matrix for this decay.
   */
  RhoDMatrix calculateDMatrix(vector<RhoDMatrix>);

  /** 
   * Calculate the rho matrix for a given outgoing particle.
   */
  RhoDMatrix calculateRhoMatrix(int,RhoDMatrix,vector<RhoDMatrix>);

  /** 
   * Contract the matrix element with the rho matrix of the 
   * incoming particle.
   */
  inline Complex contract(RhoDMatrix &);

public:
  
  /** 
   * Access to the individual helicity components. 
   */

  /** 
   * 2 body decay.
   */
  inline Complex   operator () (int,int,int) const;
  inline Complex & operator () (int,int,int);

  /** 
   * 3 body decay. 
   */
  inline Complex   operator () (int,int,int,int) const;
  inline Complex & operator () (int,int,int,int);

  /**
   * 4 body decay.
   */
  inline Complex   operator () (int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int);

  /**
   * 5 body decay.
   */
  inline Complex   operator () (int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int);

  /**
   * 6 body decay.
   */
  inline Complex   operator () (int,int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int,int);

  /**
   * n body decay.
   */
  inline Complex   operator () (vector<int>) const;
  inline Complex & operator () (vector<int>);
  
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
  mutable int _nout;

  /**
   * Spin of the incoming particle as 2s+1.
   */
  mutable int _inspin;

  /**
   * Spins of the outgonig particles.
   */
  mutable vector<int> _outspin;

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
}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of DecayMatrixElement.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::DecayMatrixElement,1> {
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
  static string className() { return "/Herwig++/Helicity/DecayMatrixElement"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "hwCorrelations.so"; }

};

}

#include "DecayMatrixElement.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayMatrixElement.tcc"
#endif

#endif /* HERWIG_DecayMatrixElement_H */
