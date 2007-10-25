// -*- C++ -*-
#ifndef HERWIG_NBodyDecayConstructorBase_H
#define HERWIG_NBodyDecayConstructorBase_H
//
// This is the declaration of the NBodyDecayConstructorBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Utilities/Exception.h"
#include "NBodyDecayConstructorBase.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * This is the base class for NBodyDecayConstructors. An N-body 
 * decay constructor should inherit from this and implement the 
 * DecayList virtual funtcion to create the decays and decayers.  
 *
 * @see \ref NBodyDecayConstructorInterfaces "The interfaces"
 * defined for NBodyDecayConstructor. 
 */
class NBodyDecayConstructorBase: public Interfaced {

public:

  /**
   * The default constructor.
   */
  inline NBodyDecayConstructorBase();

  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const vector<PDPtr> & part)=0;

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<NBodyDecayConstructorBase> initNBodyDecayConstructorBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NBodyDecayConstructorBase & operator=(const NBodyDecayConstructorBase &);

};

  /** An Exception class that can be used by all inheriting classes to
   * indicate a setup problem. */
  class NBodyDecayConstructorError : public Exception {};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NBodyDecayConstructorBase. */
template <>
struct BaseClassTrait<Herwig::NBodyDecayConstructorBase,1> {
  /** Typedef of the first base class of NBodyDecayConstructorBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NBodyDecayConstructorBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NBodyDecayConstructorBase>
  : public ClassTraitsBase<Herwig::NBodyDecayConstructorBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NBodyDecayConstructorBase"; }      /** Return the name of the shared library be loaded to get
   *  access to the NBodyDecayConstructorBase class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwModelGenerator.so"; }
};

/** @endcond */

}

#include "NBodyDecayConstructorBase.icc"

#endif /* HERWIG_NBodyDecayConstructorBase_H */
