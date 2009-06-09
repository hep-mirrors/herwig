// -*- C++ -*-
//
// METRP2to2.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_METRP2to2_H
#define HERWIG_METRP2to2_H
//
// This is the declaration of the METRP2to2 class.
//

#include "Herwig++/MatrixElement/HwME2to2Base.h"
#include "ThePEG/Repository/UseRandom.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The METRP2to2 class implements the matrix elements for
 * Transplanckian \f$2\to2\f$ scattering process
 */
class METRP2to2: public HwME2to2Base {

public:

  /**
   * The default constructor.
   */
  inline METRP2to2();

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
   * Matrix element for \f$q_i q_j \to q_i q_j\f$.
   */
  inline double ME() const;

  
  //@}
  
  inline double Any(double s, double t) const;
  inline double fny(double n, double bc, double q) const;
  inline double METRP2to2::fnyasympt(double n, double y) const;
  inline double METRP2to2::ffile(string inputf, double x) const;
  inline double METRP2to2::interp(double y, double f0, double f1, double y0, double y1) const ;
  inline double METRP2to2::fpoint(int n, double x) const;
  inline double METRP2to2::bccalc(double s) const;

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
  static ClassDescription<METRP2to2> initMETRP2to2;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  METRP2to2 & operator=(const METRP2to2 &);

private:

  /**
   *  Maximum numbere of quark flavours to include
   */
  unsigned int _maxflavour;

  unsigned int _ndim;

  double _planckmass;

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
 *  base classes of METRP2to2. */
template <>
struct BaseClassTrait<Herwig::METRP2to2,1> {
  /** Typedef of the first base class of METRP2to2. */
  typedef Herwig::HwME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the METRP2to2 class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::METRP2to2>
  : public ClassTraitsBase<Herwig::METRP2to2> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::METRP2to2"; }
  /**
   * The name of a file containing the dynamic library where the class
   * METRP2to2 is implemented. It may also include several, space-separated,
   * libraries if the class MEQCD2to2Fast depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "METRP2to2.so"; }
};

/** @endcond */

}

#include "METRP2to2.icc"

#endif /* HERWIG_METRP2to2_H */
