// -*- C++ -*-
//
// SMWWWVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMWWWVertex_H
#define HERWIG_SMWWWVertex_H
//
// This is the declaration of the SMWWWVertex class.
//
#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::Direction;

/** \ingroup Helicity
 *
 * This is the implementation fo the vertex for the coupling of three
 * standard Model electroweak bosons.
 *
 * @see VVVVertex
 * @see VertexBase
 */
class SMWWWVertex: public Helicity::VVVVertex {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  SMWWWVertex();
  //@}  

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
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SMWWWVertex> initSMWWWVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMWWWVertex & operator=(const SMWWWVertex &);
  
private:
  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The factor for the \f$Z\f$ vertex.
   */
  double _zfact;

  /**
   *  The last value of the electroweak coupling calculated.
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
 * base class of SMWWWVertex.
 */
template <>
struct BaseClassTrait<Herwig::SMWWWVertex,1> {
    /** Typedef of the base class of SMWWWVertex. */
  typedef ThePEG::Helicity::VVVVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SMWWWVertex>
  : public ClassTraitsBase<Herwig::SMWWWVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::SMWWWVertex"; }

};

/** @endcond */

}


#endif /* HERWIG_SMWWWVertex_H */
