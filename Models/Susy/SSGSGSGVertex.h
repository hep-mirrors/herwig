// -*- C++ -*-
//
// SSGSGSGVertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSGSGSGVertex_H
#define HERWIG_SSGSGSGVertex_H
//
// This is the declaration of the SSGSGSGVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "SSGSGSGVertex.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the g-\f$\tilde{g}\f$-\f$\tilde{g}\f$ vertex. It 
 * inherits from FFVVertex and implements the setCoupling virtual method.
 *
 * @see FFVVertex
 */
class SSGSGSGVertex: public FFVVertex {

public:

  /**
   * The default constructor.
   */
  inline SSGSGSGVertex();

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */  
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
			   tcPDPtr part2,tcPDPtr part3);


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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<SSGSGSGVertex> initSSGSGSGVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSGSGSGVertex & operator=(const SSGSGSGVertex &);

private:

  /**
   * The value of the coupling when last evaluated
   */
  Complex _couplast;

  /**
   * The scale at which the coupling was last evaluated
   */
  Energy2 _q2last;
};
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSGSGSGVertex. */
template <>
struct BaseClassTrait<Herwig::SSGSGSGVertex,1> {
  /** Typedef of the first base class of SSGSGSGVertex. */
  typedef ThePEG::Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSGSGSGVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSGSGSGVertex>
  : public ClassTraitsBase<Herwig::SSGSGSGVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSGSGSGVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSGSGSGVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSGSGSGVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#include "SSGSGSGVertex.icc"

#endif /* HERWIG_SSGSGSGVertex_H */
