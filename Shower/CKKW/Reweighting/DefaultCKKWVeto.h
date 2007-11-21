// -*- C++ -*-
//
// DefaultCKKWVeto.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DefaultCKKWVeto_H
#define HERWIG_DefaultCKKWVeto_H
//
// This is the declaration of the DefaultCKKWVeto class.
//

#include "Herwig++/Shower/Base/ShowerVeto.h"
#include "DefaultCKKWVeto.fh"

#include "DefaultJetMeasure.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * DefaultCKKWVeto is a ShowerVeto wrapper 
 * around DefaultJetMeasure.
 *
 *@author Simon Plaetzer
 *
 * @see \ref DefaultCKKWVetoInterfaces "The interfaces"
 * defined for DefaultCKKWVeto.
 */
class DefaultCKKWVeto: public ShowerVeto {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DefaultCKKWVeto();

  /**
   * Construct giving a DefaultJetMeasure object.
   */
  inline explicit DefaultCKKWVeto (DefaultJetMeasurePtr);

  /**
   * The destructor.
   */
  virtual ~DefaultCKKWVeto();
  //@}

public:

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  inline virtual bool vetoTimeLike (tcShowerProgenitorPtr, tcShowerParticlePtr, const Branching&);

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  inline virtual bool vetoSpaceLike (tcShowerProgenitorPtr, tcShowerParticlePtr, const Branching&);

  /**
   * Enable the veto.
   */
  inline void enable (bool en = true);

  /**
   * Disable the veto.
   */
  inline void disable ();

public:

  /**
   * Set the event generator
   */
  inline void eventGenerator(tEGPtr);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


protected:

  /** @name Standard Interfaced functions. */
  //@{
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
   * The default jet resolution used
   */
  DefaultJetMeasurePtr  _resolution;

  /**
   * True, if the veto is enabled
   */
  bool _enabled;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DefaultCKKWVeto> initDefaultCKKWVeto;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DefaultCKKWVeto & operator=(const DefaultCKKWVeto &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DefaultCKKWVeto. */
template <>
struct BaseClassTrait<Herwig::DefaultCKKWVeto,1> {
  /** Typedef of the first base class of DefaultCKKWVeto. */
  typedef Herwig::ShowerVeto NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DefaultCKKWVeto class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DefaultCKKWVeto>
  : public ClassTraitsBase<Herwig::DefaultCKKWVeto> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DefaultCKKWVeto"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DefaultCKKWVeto is implemented. It may also include several, space-separated,
   * libraries if the class DefaultCKKWVeto depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "DefaultCKKWVeto.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultCKKWVeto.tcc"
#endif

#endif /* HERWIG_DefaultCKKWVeto_H */
