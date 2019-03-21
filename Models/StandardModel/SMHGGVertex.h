// -*- C++ -*-
//
// SMHGGVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHGGVertex_H
#define HERWIG_SMHGGVertex_H
//
// This is the declaration of the SMHGGVertex class.
//

#include "Herwig/Models/General/VVSLoopVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The <code>SMHGGVertex</code> class implements the 
 * setCoupling member for the Standard Model effective 
 * vertex Higgs-gluon-gluon. 
 */
class SMHGGVertex: public VVSLoopVertex {
  
public:

  /**
   * The default constructor.
   */
  SMHGGVertex();

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

  /** 
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to first particle
   *@param part3 ParticleData pointer to first particle
   */
  virtual void setCoupling (Energy2 q2, tcPDPtr part1, tcPDPtr part2, tcPDPtr part3);

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

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMHGGVertex> initSMHGGVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHGGVertex & operator=(const SMHGGVertex &) = delete;

  /**
   * Storage of couplings
   */
  //@{
  /**
   * Last value of the coupling calculated
   */
  Complex _couplast;

  /**
   * The scale \f$q^2\f$ at which coupling was last evaluated
   */
  Energy2 _q2last;
  //@}

  /**
   * Pointer to Standard Model object
   */
  tcHwSMPtr _theSM;

  /**
   * Mass of W-boson for higgs coupling
   */
  Energy _mw;

  /**
   * define quark mass scheme (fixed/running)
   */
  unsigned int massopt;

  /**
   * The minimum flavour number in quark loops
   */
  int _minloop;
  
  /**
   * The maximum flavour number in quark loops
   */
  int _maxloop;

  /**
   * Loop calculations: A1 for spin-1/2 particles (see details in ``Higgs Hunter's Guide'')
   */
  Complex Af(double labmda) const;

  /**
   * Loop calculations: W2 function (see details in NPB297,221)
   */
  Complex W2(double lambda) const;

  /**
   * Switch between two representations of coefficients (_a00,_a11,_a12,_a21,_a22,_aEp):
   * suitable for the simplified H-g-g and H-gamma-gamma vertices and 
   * suitable for the Passarino-Veltman tensor reduction scheme
   */
  unsigned int _CoefRepresentation;

};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHGGVertex. */
template <>
struct BaseClassTrait<Herwig::SMHGGVertex,1> {
  /** Typedef of the first base class of SMHGGVertex. */
  typedef Herwig::VVSLoopVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHGGVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHGGVertex>
  : public ClassTraitsBase<Herwig::SMHGGVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHGGVertex"; }
};

/** @endcond */

}

#endif /* HERWIG_SMHGGVertex_H */
