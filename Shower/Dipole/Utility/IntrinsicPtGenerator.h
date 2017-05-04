// -*- C++ -*-
//
// IntrinsicPtGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_IntrinsicPtGenerator_H
#define HERWIG_IntrinsicPtGenerator_H
//
// This is the declaration of the IntrinsicPtGenerator class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Vectors/SpinOneLorentzRotation.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 * 
 * \brief IntrinsicPtGenerator generates intrinsic pt for massless
 * incoming partons in a shower independent way.
 *
 * @see \ref IntrinsicPtGeneratorInterfaces "The interfaces"
 * defined for IntrinsicPtGenerator.
 */
class IntrinsicPtGenerator: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  IntrinsicPtGenerator();

  /**
   * The destructor.
   */
  virtual ~IntrinsicPtGenerator();
  //@}

public:

  /**
   * Generate intrinsic pt for the given incoming
   * partons and return the transformation to be
   * applied on the final state particles. Add the
   * old incoming partons to the given list.
   */
  SpinOneLorentzRotation kick(PPair& in,
			      PList& intermediates);

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
   * The mean of the Gaussian distribution for
   * the intrinsic pt of valence partons.
   */
  Energy theValenceIntrinsicPtScale;

  /**
   * The mean of the Gaussian distribution for
   * the intrinsic pt of sea partons.
   */
  Energy theSeaIntrinsicPtScale;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<IntrinsicPtGenerator> initIntrinsicPtGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IntrinsicPtGenerator & operator=(const IntrinsicPtGenerator &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of IntrinsicPtGenerator. */
template <>
struct BaseClassTrait<Herwig::IntrinsicPtGenerator,1> {
  /** Typedef of the first base class of IntrinsicPtGenerator. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the IntrinsicPtGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::IntrinsicPtGenerator>
  : public ClassTraitsBase<Herwig::IntrinsicPtGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::IntrinsicPtGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * IntrinsicPtGenerator is implemented. It may also include several, space-separated,
   * libraries if the class IntrinsicPtGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_IntrinsicPtGenerator_H */
