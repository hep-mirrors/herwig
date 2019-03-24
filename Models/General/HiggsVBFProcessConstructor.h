// -*- C++ -*-
#ifndef HERWIG_HiggsVBFProcessConstructor_H
#define HERWIG_HiggsVBFProcessConstructor_H
//
// This is the declaration of the HiggsVBFProcessConstructor class.
//

#include "HardProcessConstructor.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiggsVBFProcessConstructor class.
 *
 * @see \ref HiggsVBFProcessConstructorInterfaces "The interfaces"
 * defined for HiggsVBFProcessConstructor.
 */
class HiggsVBFProcessConstructor: public HardProcessConstructor {

public:

  /**
   * The default constructor.
   */
  HiggsVBFProcessConstructor();

  /**
   * Main function called to start constructing the diagrams for 
   * the 2->2 process
   */
  void constructDiagrams();

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiggsVBFProcessConstructor & 
  operator=(const HiggsVBFProcessConstructor &) = delete;

private:

  /**
   *  The outgoing higgs bosons
   */
  PDVector _higgs;

  /**
   *  Collision Type
   */
  bool _type;

  /**
   *  Treatment of the Higgs width
   */
  unsigned int _shapeOpt;

  /**
   *  which intermediates to include
   */
  unsigned int _intermediates;
};

}

#endif /* HERWIG_HiggsVBFProcessConstructor_H */
