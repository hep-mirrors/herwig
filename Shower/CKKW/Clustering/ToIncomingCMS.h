// -*- C++ -*-
//
// ToIncomingCMS.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ToIncomingCMS_H
#define HERWIG_ToIncomingCMS_H
//
// This is the declaration of the ToIncomingCMS class.
//

#include "PostClustering.h"
#include "ToIncomingCMS.fh"

#include "ThePEG/Vectors/LorentzRotation.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * A boost to the CMS of the incoming particles emerging after
 * clustering.
 *
 *@author Simon Plaetzer
 *
 * @see \ref ToIncomingCMSInterfaces "The interfaces"
 * defined for ToIncomingCMS.
 */
class ToIncomingCMS: public PostClustering {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~ToIncomingCMS();
  //@}

public:

  /**
   * Initialize the PostClustering object given the
   * particles after clustering has been performed.
   */
  virtual void initialize (const vector<tClusteringParticlePtr>&);

  /**
   * Apply the transformation
   */
  virtual void doTransform (tClusteringParticlePtr);

public:

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


private:

  /**
   * The boost to be applied
   */
  LorentzRotation _boost;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<ToIncomingCMS> initToIncomingCMS;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ToIncomingCMS & operator=(const ToIncomingCMS &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ToIncomingCMS. */
template <>
struct BaseClassTrait<Herwig::ToIncomingCMS,1> {
  /** Typedef of the first base class of ToIncomingCMS. */
  typedef Herwig::PostClustering NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ToIncomingCMS class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ToIncomingCMS>
  : public ClassTraitsBase<Herwig::ToIncomingCMS> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ToIncomingCMS"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ToIncomingCMS is implemented. It may also include several, space-separated,
   * libraries if the class ToIncomingCMS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ToIncomingCMS.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ToIncomingCMS.tcc"
#endif

#endif /* HERWIG_ToIncomingCMS_H */
