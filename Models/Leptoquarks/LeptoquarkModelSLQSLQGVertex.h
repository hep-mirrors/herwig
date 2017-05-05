// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LeptoquarkModelSLQSLQGVertex_H
#define HERWIG_LeptoquarkModelSLQSLQGVertex_H
//
// This is the declaration of the LeptoquarkModelSLQSLQGVertex class.

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig/Models/Leptoquarks/LeptoquarkModel.h"

namespace Herwig {
using namespace ThePEG;
    
/** \ingroup Helicity
 * 
 *  The LeptoquarkModelSLQSLQGVertex class is the implementation of the gluon
 *  coupling to pairs of scalar leptoquarks in the LeptoquarkModel. It inherits from the VSSVertex 
 *  and implements the setCoupling member
 *
 *  @see VSSVertex
 *  @see VertexBase
 */
class LeptoquarkModelSLQSLQGVertex: public VSSVertex {
  
public:
  
  /**
   * Default constructor.
   */
  LeptoquarkModelSLQSLQGVertex();
    
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
    
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<LeptoquarkModelSLQSLQGVertex> initLeptoquarkModelSLQSLQGVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  LeptoquarkModelSLQSLQGVertex & operator=(const LeptoquarkModelSLQSLQGVertex &);
  
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
 * base class of LeptoquarkModelSLQSLQGVertex.
 */
template <>
struct BaseClassTrait<Herwig::LeptoquarkModelSLQSLQGVertex,1> {
  /** Typedef of the base class of LeptoquarkModelSLQSLQGVertex. */
  typedef ThePEG::Helicity::VSSVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::LeptoquarkModelSLQSLQGVertex>
  : public ClassTraitsBase<Herwig::LeptoquarkModelSLQSLQGVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::LeptoquarkModelSLQSLQGVertex"; }

  /** 
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwLeptoquarkModel.so"; }

};

/** @endcond */
  
}


#endif /* HERWIG_LeptoquarkModelSLQSLQGVertex_H */
