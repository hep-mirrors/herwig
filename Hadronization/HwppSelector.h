// -*- C++ -*-
//
// HwppSelector.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HwppSelector_H
#define HERWIG_HwppSelector_H
//
// This is the declaration of the HwppSelector class.
//

#include "HadronSelector.h"
#include "HwppSelector.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup hadronization
 * The HwppSelector class selects the hadrons produced in cluster decay using
 * the Herwig variant of the cluster model.
 *
 * @see \ref HwppSelectorInterfaces "The interfaces"
 * defined for HwppSelector.
 */
class HwppSelector: public HadronSelector {

public:

  /**
   * The default constructor.
   */
  HwppSelector() : HadronSelector(1), _mode(1)
  {}

  /**
   *
   * This method is used to choose a pair of hadrons.
   *
   * Given the mass of a cluster and the particle pointers of its 
   * two (or three) constituents, this returns the pair of particle pointers of
   * the two hadrons with proper flavour numbers. 
   * Furthermore, the first of the two hadron must have the 
   * constituent with par1, and the second must have the constituent with par2.
   * At the moment it does *nothing* in the case that also par3 is present.
   *
   * Kupco's method is used, rather than one used in FORTRAN HERWIG
   * The idea is to build on the fly a table of all possible pairs
   * of hadrons (Had1,Had2) (that we can call "cluster decay channels")
   * which are kinematically above threshold  and have flavour 
   * Had1=(par1,quarktopick->CC()), Had2=(quarktopick,par2), where quarktopick
   * is the poniter of: 
   *    ---  d, u, s, c, b  
   *                        if either par1 or par2 is a diquark;      
   *    ---  d, u, s, c, b, dd, ud, uu, sd, su, ss, 
   *                        cd, cu, cs, cc, bd, bu, bs, bc, bb
   *                        if both par1 and par2  are quarks.
   * The weight associated with each channel is given by the product
   * of: the phase space available including the spin factor 2*J+1, 
   *     the constant weight factor for chosen idQ, 
   *     the octet-singlet isoscalar mixing factor, and finally 
   *     the singlet-decuplet weight factor.
   */
  pair<tcPDPtr,tcPDPtr> chooseHadronPair(const Energy cluMass,tcPDPtr par1, 
						   tcPDPtr par2,tcPDPtr par3 = PDPtr()) const
   ;

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HwppSelector & operator=(const HwppSelector &);

private:

  /**
   *  Which algorithm to use
   */
  unsigned int _mode;
};

}

#endif /* HERWIG_HwppSelector_H */
