// -*- C++ -*-
//
// SMHiggsMassGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHiggsMassGenerator_H
#define HERWIG_SMHiggsMassGenerator_H
//
// This is the declaration of the SMHiggsMassGenerator class.
//

#include "GenericMassGenerator.h"
#include "GenericWidthGenerator.h"

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
  SMHiggsMassGenerator() : _shape(1) {}

  /**
   * Weight for the factor for an off-shell mass
   * @param q The off-shell mass
   * @param shape The type of shape to use as for the BreitWignerShape interface
   * @return The weight.
   */
  virtual double weight(Energy q, int shape) const {
    Energy2 q2    = sqr(q);
    Energy2 mass2 = sqr(nominalMass());
    Energy2 mwidth= nominalMass()*nominalWidth();
    return BreitWignerWeight(q,shape)*(sqr(mass2-q2)+sqr(mwidth))/mwidth;
  }

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

  /**
   * Weight for the factor for an off-shell mass
   * @param q The off-shell mass
   * @param shape The type of shape to use as for the BreitWignerShape interface
   * @return The weight.
   */
  virtual InvEnergy2 BreitWignerWeight(Energy q,int shape) const {
    useMe();
    pair<Energy,Energy> widths = shape!=2 ? _hwidth->width(q,*particle()) :
      make_pair(nominalWidth(),nominalWidth());
    Energy2 q2 = sqr(q);
    Energy4 sq=sqr(q2-sqr(nominalMass()));
    Energy2 num = widths.first*q;
    double fact = 1.;
    if(_shape==1) fact *= pow<4,1>(nominalMass()/q);
    if( shape==3) num=GeV2;
    return num*fact/Constants::pi/(sq+sqr(widths.second*q)*fact);
  }

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMHiggsMassGenerator> initSMHiggsMassGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsMassGenerator & operator=(const SMHiggsMassGenerator &) = delete;

private:

  /**
   *  Option for the line-shape
   */
  unsigned int _shape;

  /**
   *  The width generator
   */
  GenericWidthGeneratorPtr _hwidth;

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

#endif /* HERWIG_SMHiggsMassGenerator_H */
