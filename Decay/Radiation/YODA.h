// -*- C++ -*-
#ifndef HERWIG_YODA_H
#define HERWIG_YODA_H
//
// This is the declaration of the YODA class.
//

#include "DecayRadiationGenerator.h"
#include "GeneralDipole.fh"
#include "FFDipole.fh"
#include "IFDipole.fh"
#include "YODA.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 * Here is the documentation of the YODA class.
 *
 * @see \ref YODAInterfaces "The interfaces"
 * defined for YODA.
 */
class YODA: public DecayRadiationGenerator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline YODA();

  /**
   * The copy constructor.
   */
  inline YODA(const YODA &);

  /**
   * The destructor.
   */
  virtual ~YODA();
  //@}

public:

  /**
   *  Member to generate the photons in the decay. This must be implemented
   *  in classes inheriting from this one to produce the radiation.
   * @param p The decaying particle
   * @param children The decay products
   * @return The decay products with additional radiation
   */
  virtual ParticleVector generatePhotons(const Particle & p,ParticleVector children);

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

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
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<YODA> initYODA;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  YODA & operator=(const YODA &);

private:

  /**
   *  The final-final dipole
   */
  FFDipolePtr _ffdipole;

  /**
   *  The initial-final dipole
   */
  IFDipolePtr _ifdipole;

  /**
   *  The general dipole
   */
  GeneralDipolePtr _generaldipole;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of YODA. */
template <>
struct BaseClassTrait<Herwig::YODA,1> {
  /** Typedef of the first base class of YODA. */
  typedef Herwig::DecayRadiationGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the YODA class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::YODA>
  : public ClassTraitsBase<Herwig::YODA> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::YODA"; }
  /** Return the name of the shared library be loaded to get
   *  access to the YODA class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwDecRad.so"; }
};

/** @endcond */

}

#include "YODA.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "YODA.tcc"
#endif

#endif /* HERWIG_YODA_H */
