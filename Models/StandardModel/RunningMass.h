// -*- C++ -*-
//
// RunningMass.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RunningMass_H
#define HERWIG_RunningMass_H
//
// This is the declaration of the RunningMass class.

#include "RunningMassBase.h"
#include "StandardModel.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Models
 * 
 *  Implementation of the 1 or 2 loop QCD running mass.
 *
 *  @see RunningMassBase
 */
class RunningMass: public RunningMassBase {
  
public:
  
  /**
   * Default constructor.
   */
  RunningMass()  : _theQCDOrder(1), _theMaxFlav(6), _lightOption(1), _heavyOption(0) {}

public:
  
  /**
   * Return the running mass for a given scale \f$q^2\f$ and particle type.
   * @param q2 The scale \f$q^2\f$.
   * @param part The ParticleData pointer
   */
  virtual Energy value(Energy2 q2,tcPDPtr part) const;
  
  /**
   * Return the masses used.
   */
  virtual vector<Energy> mass() const;

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
  static ClassDescription<RunningMass> initRunningMass;
  
  /**
   * Private and non-existent assignment operator.
   */
  RunningMass & operator=(const RunningMass &);
  
private:

  /**
   * Order in alphaS.
   */
  unsigned int _theQCDOrder;

  /**
   * The maximum number of flavours.
   */
  unsigned int _theMaxFlav;

  /**
   * The power for the running mass calculation.
   */
  vector<double> _thePower;

  /**
   * The coefficients for the running mass calculation.
   */
  vector<double> _theCoefficient;

  /**
   * Pointer to the StandardModel object.
   */
  tcSMPtr _theStandardModel;

  /**
   *  Option to use pole masses for u,d,s
   */
  unsigned int _lightOption;

  /**
   *  Option to use pole masses for c,b
   */
  unsigned int _heavyOption;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of RunningMass.
 */
template <>
struct BaseClassTrait<Herwig::RunningMass,1> {
  /** Typedef of the base class of RunningMass. */
  typedef Herwig::RunningMassBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::RunningMass>
  : public ClassTraitsBase<Herwig::RunningMass> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::RunningMass"; }
};

/** @endcond */
  
}

#endif /* HERWIG_RunningMass_H */
