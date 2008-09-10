// -*- C++ -*-
//
// PTVeto.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PTVeto_H
#define HERWIG_PTVeto_H
//
// This is the declaration of the PTVeto class.
//

#include "ShowerVeto.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 * 
 * Veto shower emissions according to transverse
 * momentum of a branching.
 *
 * Emissions may be vetoed, if the pt is above and/or
 * below given values.
 *
 * @see \ref PTVetoInterfaces "The interfaces"
 * defined for PTVeto.
 */
class PTVeto: public ShowerVeto {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline PTVeto() : ShowerVeto(ShowerVeto::Emission) {}
  //@}

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

public:

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoTimeLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			     const Branching&b ) {
    return _vetoTimelike && b.kinematics->pT() > _maxPT &&
      b.kinematics->pT() < _minPT;
  }

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoSpaceLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			      const Branching & b) {
    return _vetoSpacelike && b.kinematics->pT() > _maxPT &&
      b.kinematics->pT() < _minPT;
  }



protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PTVeto> initPTVeto;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PTVeto & operator=(const PTVeto &);

private:

  /**
   * The pT above which emissions are to be vetoed
   */
  Energy _maxPT;

  /**
   * The pT below which emissions are to be vetoed
   */
  Energy _minPT;

  /**
   * Apply the veto to timelike showering
   */
  bool _vetoTimelike;

  /**
   * Apply the veto to spacelike showering
   */
  bool _vetoSpacelike;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PTVeto. */
template <>
struct BaseClassTrait<Herwig::PTVeto,1> {
  /** Typedef of the first base class of PTVeto. */
  typedef Herwig::ShowerVeto NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PTVeto class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PTVeto>
  : public ClassTraitsBase<Herwig::PTVeto> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PTVeto"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PTVeto is implemented. It may also include several, space-separated,
   * libraries if the class PTVeto depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_PTVeto_H */
