// -*- C++ -*-
//
// ADDModelWWWGRVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ADDModelWWWGRVertex_H
#define HERWIG_ADDModelWWWGRVertex_H
//
// This is the declaration of the ADDModelWWWGRVertex class.

#include "ThePEG/Helicity/Vertex/Tensor/VVVTVertex.h"
#include "ADDModel.h"

namespace Herwig {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The ADDModelWWWGRVertex class is the implementation of the 
 *  triple vector graviton couling in the ADD model. 
 *  It inherits from VVVTVertex and implements the setCoupling member.
 *
 *  @see VVVTVertex
 *  @see VertexBase
 */
class ADDModelWWWGRVertex: public VVVTVertex {
  
public:
  
  /**
   * Default constructor.
   */
  ADDModelWWWGRVertex();
  
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
   * @param part4 The ParticleData pointer for the foruth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4);

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
  static ClassDescription<ADDModelWWWGRVertex> initADDModelWWWGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  ADDModelWWWGRVertex & operator=(const ADDModelWWWGRVertex &) = delete;

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   * The graviton coupling.
   */
  InvEnergy kappa_;

  /**
   *  Mass ratio for the propagator
   */
  Energy r_;

  /**
   *  The last value of the coupling/
   */
  Complex couplast_;

  /**
   *  The last value of the scale, \f$q^2\f$.
   */
  Energy2 q2last_;

  /**
   *  The prefactor for the \f$Z\f$ vertex.
   */
  double zfact_;
  //@}
};
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ADDModelWWWGRVertex. 
 */
template <>
struct BaseClassTrait<Herwig::ADDModelWWWGRVertex,1> {
    /** Typedef of the base class of ADDModelWWWGRVertex. */
  typedef ThePEG::Helicity::VVVTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ADDModelWWWGRVertex>
  : public ClassTraitsBase<Herwig::ADDModelWWWGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::ADDModelWWWGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwADDModel.so"; }

};

/** @endcond */

}


#endif /* HERWIG_ADDModelWWWGRVertex_H */
