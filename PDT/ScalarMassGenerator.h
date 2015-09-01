// -*- C++ -*-
//
// ScalarMassGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ScalarMassGenerator_H
#define HERWIG_ScalarMassGenerator_H
// This is the declaration of the ScalarMassGenerator class.

#include "GenericMassGenerator.h"
#include "ScalarMassGenerator.fh"
#include "ThePEG/Config/Complex.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  The <code>ScalarMassGenerator</code> class is designed for the generation
 *  of the masses of the \f$a_0\f$ and \f$f_0\f$ mesons which have \f$K\bar{K}\f$
 *  modes close to the on-shell mass of the particle. 
 *
 *  The form based on the Flatte parameterisation of PLB63, 224, we use a weight
 * \f[\frac{1}{\pi}\frac{m\Gamma(m)}{|M^2-m^2-i\sum_ig^2_i\rho_i|^2}\f],
 * where
 * -  \f$g_i\f$  is the coupling for a given decay mode
 * -  \f$\rho_i=2p_i/m\f$ is Lorentz-invariant phase-space where \f$p_i\f$ is the 
 *     momentum release in the decay, this analytically continued
 *     below the threshold.
 * In this case the running width given by the sum of the running partial widths
 * \f[\Gamma_i(m) = 2g^2_i\frac{p_i}{m^2}\f],
 * and we differ from the Flatte approach in not analytically 
 * continuing below the threshold for the numerator.
 *
 * @see MassGenerator
 * @see GenericMassGenerator
 * 
 */
class ScalarMassGenerator: public GenericMassGenerator {

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:

  /**
   * Weight for the factor for an off-shell mass
   * @param mass The off-shell mass
   * @param shape The type of shape to use as for the BreitWignerShape interface
   * @return The weight.
   */
  virtual double weight(Energy mass,int shape) const;

  /**
   * output for the database
   */
  virtual void dataBaseOutput(ofstream &,bool);

protected:

  /**
   *  Return the full weight
   */
  virtual InvEnergy2 BreitWignerWeight(Energy q, int shape) const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ScalarMassGenerator> initScalarMassGenerator;

  /**
   * Private and non-existent assignment operator.
   */
  ScalarMassGenerator & operator=(const ScalarMassGenerator &);

private:

  /**
   * couplings for the decay channels
   */
  vector<Energy> _coupling;

  /**
   * The first outgoing mass for the channels
   */
  vector<Energy> _mass1;

  /**
   * The second outgoing mass for the channels
   */
  vector<Energy> _mass2;

  /**
   * calculated values to speed things up
   */
  //@{
  /**
   *  Maximum mass squared
   */
  vector<Energy2> _m2plus;

  /**
   *  Minimum mass squared
   */
  vector<Energy2> _m2minus;
  //@}

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ScalarMassGenerator.
 */
template <>
 struct BaseClassTrait<Herwig::ScalarMassGenerator,1> {
  /** Typedef of the base class of ScalarMassGenerator. */
  typedef Herwig::GenericMassGenerator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ScalarMassGenerator>
  : public ClassTraitsBase<Herwig::ScalarMassGenerator> {
  /** Return the class name. */
  static string className() { return "Herwig::ScalarMassGenerator"; }
};

/** @endcond */

}

#endif /* HERWIG_ScalarMassGenerator_H */
