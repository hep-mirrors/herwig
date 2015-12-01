// -*- C++ -*-
//
// SMFFHVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMFFHVertex_H
#define HERWIG_SMFFHVertex_H
//
// This is the declaration of the SMFFHVertex class.

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the vertex coupling the Standard Model Higgs
 *  to the Standard Model fermions for helicity amplitude calculations
 *
 *  @see FFSVertex
 *  @see VertexBase
 */
class SMFFHVertex: public FFSVertex {
  
public:
  
  /**
   * Default constructor.
   */
  SMFFHVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
  */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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

protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
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
  static ClassDescription<SMFFHVertex> initSMFFHVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMFFHVertex & operator=(const SMFFHVertex &);

private:

  /**
   * Pointer to the SM object.
   */
  tcHwSMPtr _theSM;

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  Last evaluation of the coupling
   */
  complex<InvEnergy> _couplast;

  /**
   *  The PDG code of the last fermion the coupling was evaluated for.
   */
  int _idlast;

  /**
   *  The last \f$q^2\f$ the coupling was evaluated at.
   */
  Energy2 _q2last;

  /**
   * The mass of the last fermion for which the coupling was evaluated.
   */
  Energy _masslast;

  /**
   *  The mass of the \f$W\f$ boson.
   */
  Energy _mw;

  /**
   * A definite fermion flavour to couple to. All other flavours are ignored,
   * if this is different from zero.
   */
  int _fermion;
  //@}
};  

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SMFFHVertex.
 */
template <>
struct BaseClassTrait<Herwig::SMFFHVertex,1> {
  /** Typedef of the base class of SMFFHVertex. */
  typedef ThePEG::Helicity::FFSVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::SMFFHVertex>
  : public ClassTraitsBase<Herwig::SMFFHVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::SMFFHVertex"; }
  
};

/** @endcond */
  
}


#endif /* HERWIG_SMFFHVertex_H */
