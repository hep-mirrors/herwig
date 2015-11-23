// -*- C++ -*-
//
// METRP2to2.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_METRP2to2_H
#define HERWIG_METRP2to2_H
//
// This is the declaration of the METRP2to2 class.
//
// The implementation of this process is based upon hep-ph/0112161 by G.F. Giudice, R. Rattazzi, J.D. Wells.

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig/Utilities/Interpolator.h"


namespace Herwig {
using namespace ThePEG;

/**
 * The METRP2to2 class implements the matrix elements for
 * Transplanckian \f$2\to2\f$ scattering process
 */
class METRP2to2: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  METRP2to2();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const { return 0; }

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const { return 0; }

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

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;
  //@}



protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
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

protected:

  /** @name Helper functions for me2. */
  //@{
  /**
   * The function which calculates b_c according to hep-ph/0112161, eq.(18)
   */
  InvEnergy bccalc(Energy2 s) const;

  /**
   * A_ny calculates part of the matrix element squared (divided by 16 * pi^2) according to hep-ph/0112161. eq.(19)
   */ 
  double A_ny(Energy2 s, Energy2 t) const;

  /**
   * fpoint initializes the matrix of pre-calculated points for the functions F_n(y) (hep-ph/0112161, eq.(20)
   */
  double fpoint(double x) const;

  /**
   * The asymptotic form of the F_n(y) functions, used for y>20, according to hep-ph/0112161, eq. (25)
   */ 
  double fnyasympt(double y) const;
  //@}

  /**
   * Setup the interpolator 
   */
  void setup_interpolator();
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  METRP2to2 & operator=(const METRP2to2 &);

private:

  /**
   * Interpolator
   */ 
  Interpolator<double, double>::Ptr _interpol;
  
  /**
   *  Maximum number of quark flavours to include
   */
  unsigned int _maxflavour;

  /**
   *  Number of Extra dimensions (>=2)
   */
  unsigned int _ndim;

  /**
   *  The Extra-dimensional Planck mass
   */
  Energy _planckmass;


  /**
   *  Processes to include
   */
  unsigned int _process;
};

}

#endif /* HERWIG_METRP2to2_H */
