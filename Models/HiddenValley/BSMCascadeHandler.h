// -*- C++ -*-
#ifndef HERWIG_BSMCascadeHandler_H
#define HERWIG_BSMCascadeHandler_H
//
// This is the declaration of the BSMCascadeHandler class.
//

#include "Herwig++/Shower/ShowerHandler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the BSMCascadeHandler class.
 *
 * @see \ref BSMCascadeHandlerInterfaces "The interfaces"
 * defined for BSMCascadeHandler.
 */
class BSMCascadeHandler: public ShowerHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BSMCascadeHandler();

  /**
   * The destructor.
   */
  virtual ~BSMCascadeHandler();
  //@}

  /**
   * The main method which manages the multiple interactions and starts
   * the shower by calling cascade(sub, lastXC).
   */
  virtual void cascade();

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

//   /**
//    * Identify the particles in the hard process and decayed particles
//    * which need to be showered
//    */
//   virtual void findShoweringParticles();

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BSMCascadeHandler> initBSMCascadeHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BSMCascadeHandler & operator=(const BSMCascadeHandler &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BSMCascadeHandler. */
template <>
struct BaseClassTrait<Herwig::BSMCascadeHandler,1> {
  /** Typedef of the first base class of BSMCascadeHandler. */
  typedef Herwig::ShowerHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BSMCascadeHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BSMCascadeHandler>
  : public ClassTraitsBase<Herwig::BSMCascadeHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BSMCascadeHandler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * BSMCascadeHandler is implemented. It may also include several, space-separated,
   * libraries if the class BSMCascadeHandler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so HwHiddenValleyModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BSMCascadeHandler_H */
