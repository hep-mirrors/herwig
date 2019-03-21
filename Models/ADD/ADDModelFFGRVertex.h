// -*- C++ -*-
//
// ADDModelFFGRVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ADDModelFFGRVertex_H
#define HERWIG_ADDModelFFGRVertex_H
//
// This is the declaration of the ADDModelFFGRVertex class.

#include "ThePEG/Helicity/Vertex/Tensor/FFTVertex.h"
#include "ADDModel.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the Randell-Sundrum model fermion-antifermion
 *  tensor vertex for helicity amplitude calculations 
 *
 *  @see FFTVertex
 *  @see VertexBase
 */
class ADDModelFFGRVertex: public FFTVertex {
  
public:
  
  /**
   * Default constructor.
   */
  ADDModelFFGRVertex();
  
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
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
  
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Calculate the propagator for a diagram.
   * @param iopt The option for the Breit-Wigner shape
   * @param q2 The scale
   * @param part The ParticleData pointer for the off-shell particle.
   * @param mass The mass if not to be taken from the ParticleData object
   * @param width The width if not to be taken from the ParticleData object
   */
  virtual Complex propagator(int iopt, Energy2 q2,tcPDPtr part,
			     Energy mass=-GeV, Energy width=-GeV);
    
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ADDModelFFGRVertex> initADDModelFFGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  ADDModelFFGRVertex & operator=(const ADDModelFFGRVertex &) = delete;

  /**
   * The coupling.
   */
  InvEnergy kappa_;

  /**
   *  Mass ratio for the propagator
   */
  Energy r_;

};
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ADDModelFFGRVertex.
 */
template <>
struct BaseClassTrait<Herwig::ADDModelFFGRVertex,1> {
    /** Typedef of the base class of ADDModelFFGRVertex. */
  typedef ThePEG::Helicity::FFTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ADDModelFFGRVertex>
  : public ClassTraitsBase<Herwig::ADDModelFFGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::ADDModelFFGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwADDModel.so"; }

};

/** @endcond */

}
#endif /* HERWIG_ADDModelFFGRVertex_H */
