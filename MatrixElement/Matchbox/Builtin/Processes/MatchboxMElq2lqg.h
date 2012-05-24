// -*- C++ -*-
//
// MatchboxMElq2lqg.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMElq2lqg_H
#define HERWIG_MatchboxMElq2lqg_H
//
// This is the declaration of the MatchboxMElq2lqg class.
//

#include "Herwig++/MatrixElement/Matchbox/Builtin/Processes/MatchboxMElP2lJetJet.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxMElq2lqg implements the matrix element
 * for charged lepton + quark -> charged lepton + quark + gluon
 *
 */
class MatchboxMElq2lqg: public MatchboxMElP2lJetJet {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxMElq2lqg();

  /**
   * The destructor.
   */
  virtual ~MatchboxMElq2lqg();
  //@}

public:

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

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
  static ClassDescription<MatchboxMElq2lqg> initMatchboxMElq2lqg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMElq2lqg & operator=(const MatchboxMElq2lqg &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MatchboxMElq2lqg. */
template <>
struct BaseClassTrait<Herwig::MatchboxMElq2lqg,1> {
  /** Typedef of the first base class of MatchboxMElq2lqg. */
  typedef Herwig::MatchboxMElP2lJetJet NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MatchboxMElq2lqg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MatchboxMElq2lqg>
  : public ClassTraitsBase<Herwig::MatchboxMElq2lqg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MatchboxMElq2lqg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MatchboxMElq2lqg is implemented. It may also include several, space-separated,
   * libraries if the class MatchboxMElq2lqg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MatchboxMElq2lqg_H */
