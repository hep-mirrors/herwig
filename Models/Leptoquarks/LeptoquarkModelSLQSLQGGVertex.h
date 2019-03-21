// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGGVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LeptoquarkModelSLQSLQGGVertex_H
#define HERWIG_LeptoquarkModelSLQSLQGGVertex_H
//
// This is the declaration of the LeptoquarkModelSLQSLQGGVertex class.

#include "ThePEG/Helicity/Vertex/Scalar/VVSSVertex.h"
#include "Herwig/Models/Leptoquarks/LeptoquarkModel.h"

namespace Herwig {
using namespace ThePEG;
    
/** \ingroup Helicity
 * 
 *  The LeptoquarkModelSLQSLQGGVertex class is the implementation of the gluon
 *  gluon coupling to pairs of scalar leptoquarks in the LeptoquarkModel. It inherits from the VVSSVertex 
 *  and implements the setCoupling member
 *
 *  @see VVSSVertex
 *  @see VertexBase
 */
class LeptoquarkModelSLQSLQGGVertex: public VVSSVertex {
  
public:
  
  /**
   * Default constructor.
   */
  LeptoquarkModelSLQSLQGGVertex();
    
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
   * @param part4 The ParticleData pointer for the third  particle.

   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3, tcPDPtr part4);

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
    
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<LeptoquarkModelSLQSLQGGVertex> initLeptoquarkModelSLQSLQGGVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  LeptoquarkModelSLQSLQGGVertex & operator=(const LeptoquarkModelSLQSLQGGVertex &) = delete;
  
  /**
   * Pointer to the model.
   */
  tcSMPtr _theModel;
  
  /**
   * The energy at which the coupling was last evaluated
   */
  Energy2 _q2last;

  /**
   * The coupling when it was last evaluated
   */
  Complex _couplast;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of LeptoquarkModelSLQSLQGGVertex.
 */
template <>
struct BaseClassTrait<Herwig::LeptoquarkModelSLQSLQGGVertex,1> {
  /** Typedef of the base class of LeptoquarkModelSLQSLQGGVertex. */
  typedef ThePEG::Helicity::VVSSVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::LeptoquarkModelSLQSLQGGVertex>
  : public ClassTraitsBase<Herwig::LeptoquarkModelSLQSLQGGVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::LeptoquarkModelSLQSLQGGVertex"; }

  /** 
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwLeptoquarkModel.so"; }

};

/** @endcond */
  
}


#endif /* HERWIG_LeptoquarkModelSLQSLQGGVertex_H */
