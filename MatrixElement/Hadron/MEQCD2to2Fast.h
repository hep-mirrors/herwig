// -*- C++ -*-
//
// MEQCD2to2Fast.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEQCD2to2Fast_H
#define HERWIG_MEQCD2to2Fast_H
//
// This is the declaration of the MEQCD2to2Fast class.
//

#include "Herwig++/MatrixElement/HwME2to2Base.h"
#include "ThePEG/Repository/UseRandom.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEQCD2to2Fast class implements the matrix elements for
 * QCD \f$2\to2\f$ scattering processes using hard coded formulae and
 * as such can not include spin correlations. It is designed to be a faster
 * replacement for MEQCD2to2 for use in the underlying event.
 *
 * @see \ref MEQCD2to2FastInterfaces "The interfaces"
 * defined for MEQCD2to2Fast.
 */
class MEQCD2to2Fast: public HwME2to2Base {

public:

  /**
   * The default constructor.
   */
  inline MEQCD2to2Fast();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
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

protected:

  /**
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   * Matrix element for \f$gg\to gg\f$.
   */
  inline double gg2ggME() const;

  /**
   * Matrix element for \f$gg\to q\bar{q}\f$
   */
  inline double gg2qqbarME() const;

  /**
   * Matrix element for \f$q\bar{q}\to gg\f$
   */
  inline double qqbar2ggME() const;

  /**
   * Matrix element for \f$qg\to qg\f$
   */
  inline double qg2qgME() const;

  /**
   * Matrix elements for \f$\bar{q}g\to \bar{q}g\f$.
   */
  inline double qbarg2qbargME() const;

  /**
   * Matrix element for \f$qq\to qq\f$
   */
  inline double qq2qqME() const;

  /**
   * Matrix element for \f$\bar{q}\bar{q}\to \bar{q}\bar{q}\f$
   */
  inline double qbarqbar2qbarqbarME() const;

  /**
   * Matrix element for \f$q\bar{q}\to q\bar{q}\f$
   */
  inline double qqbar2qqbarME() const;
  //@}

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
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEQCD2to2Fast> initMEQCD2to2Fast;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEQCD2to2Fast & operator=(const MEQCD2to2Fast &);

private:

  /**
   *  Maximum numbere of quark flavours to include
   */
  unsigned int _maxflavour;

  /**
   *  Processes to include
   */
  unsigned int _process;

  /**
   *  Colour flow
   */
  mutable unsigned int _flow;

  /**
   *  Diagram
   */
  mutable unsigned int _diagram;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEQCD2to2Fast. */
template <>
struct BaseClassTrait<Herwig::MEQCD2to2Fast,1> {
  /** Typedef of the first base class of MEQCD2to2Fast. */
  typedef Herwig::HwME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEQCD2to2Fast class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEQCD2to2Fast>
  : public ClassTraitsBase<Herwig::MEQCD2to2Fast> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEQCD2to2Fast"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEQCD2to2Fast is implemented. It may also include several, space-separated,
   * libraries if the class MEQCD2to2Fast depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadronFast.so"; }
};

/** @endcond */

}

#include "MEQCD2to2Fast.icc"

#endif /* HERWIG_MEQCD2to2Fast_H */
