// -*- C++ -*-
#ifndef HERWIG_FFMqx2qgxDipoleKernel_H
#define HERWIG_FFMqx2qgxDipoleKernel_H
//
// This is the declaration of the FFMqx2qgxDipoleKernel class.
//

#include "DipoleSplittingKernel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer, Martin Stoll, Stephen Webster
 *
 * \brief FFMqx2qgxDipoleKernel implements the q -> qg
 * splitting off a final-final dipole
 *
 */
class FFMqx2qgxDipoleKernel: public DipoleSplittingKernel {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FFMqx2qgxDipoleKernel();

  /**
   * The destructor.
   */
  virtual ~FFMqx2qgxDipoleKernel();
  //@}

public:

  /**
   * Return true, if this splitting kernel
   * applies to the given dipole index.
   */
  virtual bool canHandle(const DipoleIndex&) const;

  /**
   * Return true, if this splitting kernel is
   * the same for the given index a, as the given
   * splitting kernel for index b.
   */
  virtual bool canHandleEquivalent(const DipoleIndex& a,
				   const DipoleSplittingKernel& sk,
				   const DipoleIndex& b) const;

  /**
   * Return the emitter data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr emitter(const DipoleIndex&) const;

  /**
   * Return the emission data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr emission(const DipoleIndex&) const;

  /**
   * Return the spectator data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr spectator(const DipoleIndex&) const;

  /**
   * Evaluate this splitting kernel for the given
   * dipole splitting.
   */
  virtual double evaluate(const DipoleSplittingInfo&) const;

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
  static ClassDescription<FFMqx2qgxDipoleKernel> initFFMqx2qgxDipoleKernel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FFMqx2qgxDipoleKernel & operator=(const FFMqx2qgxDipoleKernel &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FFMqx2qgxDipoleKernel. */
template <>
struct BaseClassTrait<Herwig::FFMqx2qgxDipoleKernel,1> {
  /** Typedef of the first base class of FFMqx2qgxDipoleKernel. */
  typedef Herwig::DipoleSplittingKernel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FFMqx2qgxDipoleKernel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FFMqx2qgxDipoleKernel>
  : public ClassTraitsBase<Herwig::FFMqx2qgxDipoleKernel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::FFMqx2qgxDipoleKernel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FFMqx2qgxDipoleKernel is implemented. It may also include several, space-separated,
   * libraries if the class FFMqx2qgxDipoleKernel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_FFMqx2qgxDipoleKernel_H */
