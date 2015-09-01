// -*- C++ -*-
//
// SMFFGVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMFFGVertex_H
#define HERWIG_SMFFGVertex_H
//
// This is the declaration of the SMFFGVertex class.
//
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Helicity
 * 
 *  This is the implementation of the Standard Model quark-antiquark gluon vertex.
 * 
 * @see FFVVertex
 * @see VertexBase
 */
class SMFFGVertex: public FFVVertex {
  
public:
  
  /**
   * Default constructor.
   */
  SMFFGVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

public:
  
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
  static NoPIOClassDescription<SMFFGVertex> initSMFFGVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMFFGVertex & operator=(const SMFFGVertex &);
  
private:

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
 * base class of SMFFGVertex.
 */
template <>
struct BaseClassTrait<Herwig::SMFFGVertex,1> {
  /** Typedef of the base class of SMFFGVertex. */
  typedef ThePEG::Helicity::FFVVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SMFFGVertex>
  : public ClassTraitsBase<Herwig::SMFFGVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::SMFFGVertex"; }
  
};
  
/** @endcond */
  
}

#endif /* HERWIG_SMFFGVertex_H */
