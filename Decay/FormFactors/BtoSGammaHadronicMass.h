// -*- C++ -*-
//
// BtoSGammaHadronicMass.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_BtoSGammaHadronicMass_H
#define HERWIG_BtoSGammaHadronicMass_H
//
// This is the declaration of the BtoSGammaHadronicMass class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "BtoSGammaHadronicMass.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 * The BtoSGammaHadronicMass class is the base class for the implementation
 * of models of the hadronic mass spectrum in \f$B\to s\gamma\f$ decays.
 * Classes inheriting from this class should implement the hadronicMass()
 * member which should return a value of the hadronic mass selected from
 * the distribution.
 *
 * The parameters relating to the minimum and maximum values of the mass are stored
 * in this class.
 *
 */
class BtoSGammaHadronicMass: public Interfaced {

public:

  /**
   * The default constructor.
   */
  BtoSGammaHadronicMass() : _minMass(825*MeV),_maxMass(5300*MeV) {}

  /**
   * Virtual member which must be implemented in classes inheriting from this
   * class to return the hadronic mass.
   * @param mb The mass of the decaying B meson
   * @param mquark The minimum mass of the hadronic system based on the consistuent quark
   * masses.
   * @return The hadronic mass
   */
  virtual Energy hadronicMass(Energy mb,Energy mquark) =0;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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

  /** @name Access to the limits on the mass. */
  //@{
  /**
   *  Minimum mass
   */
  Energy minMass() const {return _minMass;}

  /**
   *  Maximum mass
   */
  Energy maxMass() const {return _maxMass;}
  //@}

  /** @name Functions for the fermi motion needed in classes inheriting from this */
  //@{
  /**
   * Exponential function of the form, \f$N(1-x)^ae^{-3\bar{\Lambda}^2x/\lambda_1}\f$,
   * where 
   * \f$x=k_+/\bar{\Lambda}\f$ taken from hep-ph/9805303
   * @param scale The energy scale, \f$k_+\f$, at which to evaluate the function.
   * @param lambda The hadronic scale, \f$\bar{\Lambda}\f$
   * @param a The shape parameter, \f$a\f$.
   * @param norm The normalisation, \f$N\f$.
   * @param lambda1 Scale related to kinetic energy of b quark, \f$\lambda_1\f$.
   */
  InvEnergy exponentialFermiFunction(Energy scale,Energy lambda, double a,
				     InvEnergy norm,Energy2 lambda1 ) const {
    double x(scale/lambda);
    return norm*pow(1.-x,a)*exp(-3.*sqr(lambda)/lambda1*x);
  }
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<BtoSGammaHadronicMass> initBtoSGammaHadronicMass;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BtoSGammaHadronicMass & operator=(const BtoSGammaHadronicMass &) = delete;

private:

  /**
   * The minimum value of the hadronic mass
   */
  Energy _minMass;

  /**
   * The maximum value of the hadronic mass
   */
  Energy _maxMass;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BtoSGammaHadronicMass. */
template <>
struct BaseClassTrait<Herwig::BtoSGammaHadronicMass,1> {
  /** Typedef of the first base class of BtoSGammaHadronicMass. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BtoSGammaHadronicMass class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BtoSGammaHadronicMass>
  : public ClassTraitsBase<Herwig::BtoSGammaHadronicMass> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BtoSGammaHadronicMass"; }
};

/** @endcond */

}

#endif /* HERWIG_BtoSGammaHadronicMass_H */
