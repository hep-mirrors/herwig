// -*- C++ -*-
//
// MatchboxMEqg2llbarq.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMEqg2llbarq_H
#define HERWIG_MatchboxMEqg2llbarq_H
//
// This is the declaration of the MatchboxMEqg2llbarq class.
//

#include "Herwig++/MatrixElement/Matchbox/Builtin/Processes/MatchboxMEPP2llbarJet.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxMEqg2llbarq implements the matrix element
 * for quark + gluon -> quark + charged lepton pair
 *
 */
class MatchboxMEqg2llbarq: public MatchboxMEPP2llbarJet {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxMEqg2llbarq();

  /**
   * The destructor.
   */
  virtual ~MatchboxMEqg2llbarq();
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

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;

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
   * Set which of the incoming partons
   * is the incoming gluon.
   */
  int theWhichGluon;
  
  /**
   * A user defined scale; if zero, the center of mass
   * energy is used.
   */
  Energy theUserScale;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MatchboxMEqg2llbarq> initMatchboxMEqg2llbarq;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMEqg2llbarq & operator=(const MatchboxMEqg2llbarq &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MatchboxMEqg2llbarq. */
template <>
struct BaseClassTrait<Herwig::MatchboxMEqg2llbarq,1> {
  /** Typedef of the first base class of MatchboxMEqg2llbarq. */
  typedef Herwig::MatchboxMEPP2llbarJet NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MatchboxMEqg2llbarq class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MatchboxMEqg2llbarq>
  : public ClassTraitsBase<Herwig::MatchboxMEqg2llbarq> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MatchboxMEqg2llbarq"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MatchboxMEqg2llbarq is implemented. It may also include several, space-separated,
   * libraries if the class MatchboxMEqg2llbarq depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MatchboxMEqg2llbarq_H */
