// -*- C++ -*-
#ifndef HERWIG_DecayMatrixElement_H
#define HERWIG_DecayMatrixElement_H
//
// This is the declaration of the <!id>DecayMatrixElement<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  Implementation of the complex matrix element for a decay
//  an arbitary number of external particles are supported.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="RhoDMatrix.html">RhoDMatrix.h</a>,
// <a href="ProductionMatrixElement.html">ProductionMatrixElement.h</a>,
// <a href="DecayVertex.html">DecayVertex.h</a>.
// 
// Author: Peter Richardson
//

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
#include <ThePEG/Helicity/RhoDMatrix.h>
// #include "DecayMatrixElement.fh"
// #include "DecayMatrixElement.xh"

namespace Herwig {
namespace Helicity {

using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;

class DecayMatrixElement: public Base {
  
public:
      
  inline DecayMatrixElement();

  inline DecayMatrixElement(int,int,int);
  // constructor for two body decay

  inline DecayMatrixElement(int,int,int,int);
  // constructor for three body decay

  inline DecayMatrixElement(int,int,int,int,int);
  // constructor for four body decay

  inline DecayMatrixElement(int,int,int,int,int,int);
  // constructor for five body decay

  inline DecayMatrixElement(int,int,int,int,int,int,int);
  // constructor for six body decay

  inline DecayMatrixElement(int,vector<int>);
  inline DecayMatrixElement(vector<int>);
  // constructors for arbitray body decay

  inline DecayMatrixElement(const DecayMatrixElement &);
  virtual ~DecayMatrixElement();
  // Standard ctors and dtor.
  
public: 
      
  inline int inspin();
  // get the spin of the incoming particle

  inline vector<int> outspin();
  // get the spins of the outgoing particles

public:
  RhoDMatrix calculateDMatrix(vector<RhoDMatrix>);
  // calculate the decay matrix for this decay

  RhoDMatrix calculateRhoMatrix(int,RhoDMatrix,vector<RhoDMatrix>);
  // calculate the rho matrix for a given outgoing particle

  inline Complex contract(RhoDMatrix &);
  // contract the matrix element with the rho matrix of the incoming particle

public:
  
  // access to the individual helicity components
  inline Complex   operator () (int,int,int) const;
  inline Complex & operator () (int,int,int);
  // 2 body decay
  inline Complex   operator () (int,int,int,int) const;
  inline Complex & operator () (int,int,int,int);
  // 3 body decay
  inline Complex   operator () (int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int);
  // 4 body decay
  inline Complex   operator () (int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int);
  // 5 body decay
  inline Complex   operator () (int,int,int,int,int,int,int) const;
  inline Complex & operator () (int,int,int,int,int,int,int);
  // 6 body decay
  inline Complex   operator () (vector<int>) const;
  inline Complex & operator () (vector<int>);
  // n body decay
  
public:
  
  inline void reset(const DecayMatrixElement &) const;
  // reset the matrix element
  
public:
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
private:
  
  static NoPIOClassDescription<DecayMatrixElement> initDecayMatrixElement;
  // Describe a concrete class without persistent data.
  
  DecayMatrixElement & operator=(const DecayMatrixElement &);
  // Private and non-existent assignment operator.
  
private:
  
  inline void setMESize();
  // set the size of the vector containing the matrix element
  
private:
  
  mutable int _nout;
  // number of outgoing particles
  mutable int _inspin;
  // spin of the incoming particle as 2s+1
  mutable vector<int> _outspin;
  // spins of the outgonig particles
  mutable vector<Complex> _matrixelement;
  // storage of the matrix element, a vector is better for memory usage
  mutable vector<int> _constants;
  // constants needed to map the index of the vector to a helicity structure
};

}
}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of DecayMatrixElement.
template <>
struct BaseClassTrait<Herwig::Helicity::DecayMatrixElement,1> {
  typedef Base NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::Helicity::DecayMatrixElement>
  : public ClassTraitsBase<Herwig::Helicity::DecayMatrixElement> {
  static string className() { return "/Herwig++/Helicity/DecayMatrixElement"; }
  // Return the class name.
  static string library() { return "hwCorrelations.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "DecayMatrixElement.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayMatrixElement.tcc"
#endif

#endif /* HERWIG_DecayMatrixElement_H */
