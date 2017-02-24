// -*- C++ -*-
#ifndef HERWIG_HiggsVectorBosonProcessConstructor_H
#define HERWIG_HiggsVectorBosonProcessConstructor_H
//
// This is the declaration of the HiggsVectorBosonProcessConstructor class.
//

#include "HardProcessConstructor.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiggsVectorBosonProcessConstructor class.
 *
 * @see \ref HiggsVectorBosonProcessConstructorInterfaces "The interfaces"
 * defined for HiggsVectorBosonProcessConstructor.
 */
class HiggsVectorBosonProcessConstructor: public HardProcessConstructor {

public:

  /**
   * The default constructor.
   */
  HiggsVectorBosonProcessConstructor();

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HiggsVectorBosonProcessConstructor> 
  initHiggsVectorBosonProcessConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiggsVectorBosonProcessConstructor & 
  operator=(const HiggsVectorBosonProcessConstructor &);

private:

  /**
   *  The allowed outgoing vector bosons
   */
  PDVector _vector;

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
   *  The shower coupling for the Matrix Element corrections
   */
  ShowerAlphaPtr _alpha;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiggsVectorBosonProcessConstructor. */
template <>
struct BaseClassTrait<Herwig::HiggsVectorBosonProcessConstructor,1> {
  /** Typedef of the first base class of HiggsVectorBosonProcessConstructor. */
  typedef Herwig::HardProcessConstructor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiggsVectorBosonProcessConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiggsVectorBosonProcessConstructor>
  : public ClassTraitsBase<Herwig::HiggsVectorBosonProcessConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiggsVectorBosonProcessConstructor"; }
};

/** @endcond */

}

#endif /* HERWIG_HiggsVectorBosonProcessConstructor_H */
