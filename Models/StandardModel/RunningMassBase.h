// -*- C++ -*-
//
// RunningMassBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RunningMassBase_H
#define HERWIG_RunningMassBase_H
//
// This is the declaration of the RunningMassBase class.

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Models
 * 
 *  Base class for running mass calculations.
 */
class RunningMassBase: public Interfaced {
  
public:
  
  /**
   * Return the running mass for a given scale \f$q^2\f$ and particle type.
   * @param q2 The scale \f$q^2\f$.
   * @param part The ParticleData pointer
   */
  virtual Energy value(Energy2 q2,tcPDPtr part) const = 0;
 
  /**
   * Return the masses used.
   */
  virtual vector<Energy> mass() const = 0;

  /**
   * Return the \f$i\f$ th element of the mass array.
   * @param i The element to return
   */
  Energy massElement(unsigned int i) const {return _theMass[i];}

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
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<RunningMassBase> initRunningMassBase;
  
  /**
   * Private and non-existent assignment operator.
   */
  RunningMassBase & operator=(const RunningMassBase &);
  
private:
  
  /**
   * Flavour thresholds and the masses, set at initialization.
   */
  vector<Energy> _theMass;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of RunningMassBase.
 */
template <>
struct BaseClassTrait<Herwig::RunningMassBase,1> {
  /** Typedef of the base class of RunningMassBase. */
  typedef Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::RunningMassBase>
  : public ClassTraitsBase<Herwig::RunningMassBase> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::RunningMassBase"; }
};

/** @endcond */
  
}

#endif /* HERWIG_RunningMassBase_H */
