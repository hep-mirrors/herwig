// -*- C++ -*-
//
// SSWSSVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSWSSVertex_H
#define HERWIG_SSWSSVertex_H
//
// This is the declaration of the SSWSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "MSSM.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This is the implementation of the coupling of an SM gauge boson to
 * a pair of sfermions.
 * It inherits from VSSVertex and implements the setCoupling() method.
 *
 * @see VSSVertex
 */
class SSWSSVertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  SSWSSVertex();

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
   * Initialize this object after the setup phase before saving an
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
  static ClassDescription<SSWSSVertex> initSSWSSVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSWSSVertex & operator=(const SSWSSVertex &);

  /**
   * Value of \f$sin(\theta_w)\f$
   */
  double _sw;

  /**
   * Value of \f$cos(\theta_w)\f$
   */
  double _cw;
  
  /**
   * Scale at which the coupling was last evaluated
   */
  Energy2 _q2last;

  /**
   * Value of coupling when last evaluated
   */
  Complex _couplast;

  /**
   * Stau mixing matrix
   */
  tMixingMatrixPtr _stau;
  
  /**
   * Stop mixing matrix
   */
  tMixingMatrixPtr _stop;

  /**
   * Sbottom mixing matrix
   */
  tMixingMatrixPtr _sbottom;
  
  /**
   * The up type sfermion present when the vertex was evaluated. 
   */
  long _ulast;

  /**
   * The down type sfermion present when the vertex was evaluated. 
   */
  long _dlast;

  /**
   * The gauge boson present when the vertex was last evaluated. 
   */
  long _gblast;
  
  /**
   * The value of the mixing matrix dependent part when the vertex was
   * last evaluated
   */
  Complex _factlast;
  
};
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSWSSVertex. */
template <>
struct BaseClassTrait<Herwig::SSWSSVertex,1> {
  /** Typedef of the first base class of SSWSSVertex. */
  typedef ThePEG::Helicity::VSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSWSSVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSWSSVertex>
  : public ClassTraitsBase<Herwig::SSWSSVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSWSSVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSWSSVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSWSSVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSWSSVertex_H */
