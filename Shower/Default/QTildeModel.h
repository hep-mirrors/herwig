// -*- C++ -*-
//
// QTildeModel.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_QTildeModel_H
#define HERWIG_QTildeModel_H
//
// This is the declaration of the QTildeModel class.
//

#include "Herwig++/Shower/Base/ShowerModel.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * The QTildeModel class inherits from the ShowerModel class and implements the
 * checkConsistency member for the default Herwig++ Shower.
 *
 * @see \ref QTildeModelInterfaces "The interfaces"
 * defined for QTildeModel.
 */
class QTildeModel: public ShowerModel {

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  The implementation of the virtual member from the base class to
   *  check that the correct objects are loaded
   */
  virtual void checkConsistency() throw(InitException);

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<QTildeModel> initQTildeModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeModel & operator=(const QTildeModel &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QTildeModel. */
template <>
struct BaseClassTrait<Herwig::QTildeModel,1> {
  /** Typedef of the first base class of QTildeModel. */
  typedef Herwig::ShowerModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QTildeModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QTildeModel>
  : public ClassTraitsBase<Herwig::QTildeModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QTildeModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QTildeModel is implemented. It may also include several, space-separated,
   * libraries if the class QTildeModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_QTildeModel_H */
