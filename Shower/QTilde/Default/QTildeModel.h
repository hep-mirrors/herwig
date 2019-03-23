// -*- C++ -*-
//
// QTildeModel.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_QTildeModel_H
#define HERWIG_QTildeModel_H
//
// This is the declaration of the QTildeModel class.
//

#include "Herwig/Shower/QTilde/Base/ShowerModel.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * The QTildeModel class inherits from the ShowerModel class and implements the
 * checkConsistency member for the default Herwig Shower.
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
  virtual void checkConsistency();

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeModel & operator=(const QTildeModel &) = delete;

};

}

#endif /* HERWIG_QTildeModel_H */
