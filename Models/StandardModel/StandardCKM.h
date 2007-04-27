// -*- C++ -*-
#ifndef HERWIG_StandardCKM_H
#define HERWIG_StandardCKM_H
//
// This is the declaration of the StandardCKM class.
//
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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline StandardCKM();

  /**
   * Copy-constructor.
   */
  inline StandardCKM(const StandardCKM &);

  /**
   * Destructor.
   */
  virtual ~StandardCKM();
  //@}

public:
  
  /**
   * Return the matrix of squared CKM matrix elements. The returned
   * matrix should be for \a nf families.
   */
  virtual vector< vector<double> >  getMatrix(unsigned int nf) const;

  /**
   * Return the matrix of CKM matrix elements. The returned
   * matrix should be for \a nf families.
   */
  virtual vector< vector<Complex> > getUnsquaredMatrix(unsigned int nf) const;
  
public:
  
  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}    
  
private:
  
  /**
   * The \f$\theta_{12}\f$ angle.
   */
  double theta12;

  /**
   * The \f$\theta_{13}\f$ angle.
   */
  double theta13;

  /**
   * The \f$\theta_{23}\f$ angle.
   */
  double theta23;

  /**
   * The \f$\delta\f$ angle describing the phase.
   */
  double delta;
  
private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StandardCKM> initStandardCKM;
  
  /**
   * Private and non-existent assignment operator.
   */
  StandardCKM & operator=(const StandardCKM &);
  
};

}

#include "StandardCKM.icc"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of StandardCKM. */
template <>
struct BaseClassTrait<Herwig::StandardCKM,1> {
  /** Typedef of the first base class of StandardCKM. */
  typedef CKMBase NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  StandardCKM class and the shared object where it is
 *  defined. */
template <>
struct ClassTraits<Herwig::StandardCKM>: public ClassTraitsBase<Herwig::StandardCKM> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::StandardCKM"; }
};

/** @endcond */

}


#endif /* HERWIG_StandardCKM_H */
