// -*- C++ -*-
#ifndef HERWIG_StandardCKM_H
#define HERWIG_StandardCKM_H
//
// This is the declaration of the StandardCKM class.

#include <ThePEG/Config/Complex.h>
#include <ThePEG/StandardModel/CKMBase.h>
// #include "StandardCKM.fh"
// #include "StandardCKM.xh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Models
 * 
 *  StandardCKM inherits from CKMBase and implements the standard 
 *  parameterization of the CKM matrix in terms of three angles and 
 *  a phase. It provides access to the unsquared matrix from helicity 
 *  amplitude calculations.
 *
 * @see CKMBase
 */
class StandardCKM: public CKMBase {

public:

  /**
   * Standard ctors and dtor.
   */
  inline StandardCKM();
  inline StandardCKM(const StandardCKM &);
  virtual ~StandardCKM();
  
  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  
public:
  
  /**
   * Return the matrix of squared matrix elements.
   */
  virtual vector< vector<double> >  getMatrix(unsigned int nFamilies) const;

  /**
   * Return the matrix of matrix elements.
   */
  virtual vector< vector<Complex> >
  getUnsquaredMatrix(unsigned int nFamilies) const;
  
public:
  
  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  
  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();
  
protected:
  
  /**
   * Standard clone methods.
   */
  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;
  
  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void
  rebind(const TranslationMap & trans) throw(RebindException);
  
  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();
  
private:
  
  double theta12;
  double theta13;
  double theta23;
  double delta;
  
private:
  
  static ClassDescription<StandardCKM> initStandardCKM;
  
  /**
   * Private and non-existent assignment operator.
   */
  StandardCKM & operator=(const StandardCKM &);
  
};

}

#include "StandardCKM.icc"

namespace ThePEG {
template <>
struct BaseClassTrait<Herwig::StandardCKM,1> {
  typedef CKMBase NthBase;
};

template <>
struct ClassTraits<Herwig::StandardCKM>: public ClassTraitsBase<Herwig::StandardCKM> {
  static string className() { return "/Herwig++/StandardCKM"; }
  static string library() { return "libHwStandardModel.so"; }
};

}


#endif /* HERWIG_StandardCKM_H */
