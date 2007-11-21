// -*- C++ -*-
//
// SSGSSVertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSGSSVertex_H
#define HERWIG_SSGSSVertex_H
//
// This is the declaration of the SSGSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "SSGSSVertex.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * The SSGSSVertex implements the coupling of a gluon to 2 sfermions. 
 * It inherits from VSSVertex and implements the setCoupling method.
 * 
 * @see VSSVertex
 */
class SSGSSVertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  SSGSSVertex();

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
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSGSSVertex> initSSGSSVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSGSSVertex & operator=(const SSGSSVertex &);

 /**
   * Pointer to the standard model
   */
  tSusyBasePtr _theSS;
  
  /**
   * Store the value of the coupling when last evaluated
   */
  Complex _couplast;

  /**
   * Store the scale at which coupling was last evaluated
   */
  Energy2 _q2last;
};
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSGSSVertex. */
template <>
struct BaseClassTrait<Herwig::SSGSSVertex,1> {
  /** Typedef of the first base class of SSGSSVertex. */
  typedef ThePEG::Helicity::VSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSGSSVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSGSSVertex>
  : public ClassTraitsBase<Herwig::SSGSSVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSGSSVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSGSSVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSGSSVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}


#include "SSGSSVertex.icc"

#endif /* HERWIG_SSGSSVertex_H */
