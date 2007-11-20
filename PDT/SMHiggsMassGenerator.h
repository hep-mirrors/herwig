// -*- C++ -*-
//
// SMHiggsMassGenerator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHiggsMassGenerator_H
#define HERWIG_SMHiggsMassGenerator_H
//
// This is the declaration of the SMHiggsMassGenerator class.
//

#include "GenericMassGenerator.h"
#include "SMHiggsWidthGenerator.h"
#include "SMHiggsMassGenerator.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The SMHiggsMassGenerator class implements the generation
 * of the Higgs boson mass according to the prescription of hep-ph/9505211
 *
 * @see \ref SMHiggsMassGeneratorInterfaces "The interfaces"
 * defined for SMHiggsMassGenerator.
 */
class SMHiggsMassGenerator: public GenericMassGenerator {

public:

  /**
   * The default constructor.
   */
  inline SMHiggsMassGenerator();

  /**
   * Weight for the factor for an off-shell mass
   * @param mass The off-shell mass
   * @param shape The type of shape to use as for the BreitWignerShape interface
   * @return The weight.
   */
  inline virtual double weight(Energy mass,int shape) const;

  /**
   * Weight for the factor for an off-shell mass
   * @param mass The off-shell mass
   * @param shape The type of shape to use as for the BreitWignerShape interface
   * @return The weight.
   */
  inline InvEnergy2 BreitWignerWeight(Energy mass,int shape) const;

  /**
   * Return true if this mass generator can handle the given particle type.
   * @param part The particle data pointer of the particle.
   * @return True ig this class can handle the particle and false otherwise
   */
  bool accept(const ParticleData & part) const;

  /**
   * output for the database
   */
  virtual void dataBaseOutput(ofstream &,bool);

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMHiggsMassGenerator> initSMHiggsMassGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsMassGenerator & operator=(const SMHiggsMassGenerator &);

private:

  /**
   *  Option for the line-shape
   */
  unsigned int _shape;

  /**
   *  The width generator
   */
  SMHiggsWidthGeneratorPtr _hwidth;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHiggsMassGenerator. */
template <>
struct BaseClassTrait<Herwig::SMHiggsMassGenerator,1> {
  /** Typedef of the first base class of SMHiggsMassGenerator. */
  typedef Herwig::GenericMassGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHiggsMassGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHiggsMassGenerator>
  : public ClassTraitsBase<Herwig::SMHiggsMassGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHiggsMassGenerator"; }
};

/** @endcond */

}

#include "SMHiggsMassGenerator.icc"

#endif /* HERWIG_SMHiggsMassGenerator_H */
