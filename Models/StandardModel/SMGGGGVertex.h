// -*- C++ -*-
//
// SMGGGGVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMGGGGVertex_H
#define HERWIG_SMGGGGVertex_H
//
// This is the declaration of the SMGGGGVertex class.
//
#include "ThePEG/Helicity/Vertex/Vector/VVVVVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The SMGGGGVertex class is the implementation of the 
 *  Standard Model quartic gluon vertex. It inherits from 
 *  VVVVVertex and implements the setCoupling member.
 *
 *  @see VVVVVertex
 *  @see VertexBase
 */
class SMGGGGVertex: public VVVVVertex {
  
public:
  
  /**
   * Default constructor.
   */
  SMGGGGVertex();
  
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
   * @param part4 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3,tcPDPtr part4);
  
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
   * Describe a concrete class with persistent data.
   */
  static NoPIOClassDescription<SMGGGGVertex> initSMGGGGVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMGGGGVertex & operator=(const SMGGGGVertex &) = delete;

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the strong coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;
  //@}
};
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SMGGGGVertex.
 */
template <>
struct BaseClassTrait<Herwig::SMGGGGVertex,1> {
    /** Typedef of the base class of SMGGGGVertex. */
  typedef ThePEG::Helicity::VVVVVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SMGGGGVertex>
  : public ClassTraitsBase<Herwig::SMGGGGVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::SMGGGGVertex"; }

};

/** @endcond */

}


#endif /* HERWIG_SMGGGGVertex_H */
