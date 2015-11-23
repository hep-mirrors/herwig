// -*- C++ -*-
//
// TTbAModelZPQQVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TTbAModelZPQQVertex_H
#define HERWIG_TTbAModelZPQQVertex_H
//
// This is the declaration of the TTbAModelZPQQVertex class.

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/Models/TTbAsymm/TTbAModel.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the vertex coupling the Standard Model Higgs
 *  to the Standard Model fermions for helicity amplitude calculations
 *
 *  @see FFVVertex
 *  @see VertexBase
 */
class TTbAModelZPQQVertex: public FFVVertex {
  
public:
  
  /**
   * Default constructor.
   */
  TTbAModelZPQQVertex();
  
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
  static ClassDescription<TTbAModelZPQQVertex> initTTbAModelZPQQVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  TTbAModelZPQQVertex & operator=(const TTbAModelZPQQVertex &);

   /**
   * Pointer to the model object.
   */
  tcSMPtr _theModel;


private:

  /**
   * Storage of the couplings.
   */
  //@{

  /**
   *  Z prime coupling to top-up (left-handed)
   */
  double _cZPTU_L;
  

  /**
   *  Z prime coupling to top-up (right-handed)
   */
  double _cZPTU_R;

  /**
   *  Z prime coupling to up-upbar (left-handed)
   */
  double _cZPUU_L;
  

  /**
   *  Z prime coupling to up-upbar (right-handed)
   */
  double _cZPUU_R;

  /**
   *  Z prime coupling to charm-charmbar (left-handed)
   */
  double _cZPCC_L;
  

  /**
   *  Z prime coupling to charm-charmbar (right-handed)
   */
  double _cZPCC_R;


  /**
   * Model selector
   */
  int _models;


  //@}
};  

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of TTbAModelZPQQVertex.
 */
template <>
struct BaseClassTrait<Herwig::TTbAModelZPQQVertex,1> {
  /** Typedef of the base class of TTbAModelZPQQVertex. */
  typedef ThePEG::Helicity::FFVVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
  template <>
  
struct ClassTraits<Herwig::TTbAModelZPQQVertex>
  : public ClassTraitsBase<Herwig::TTbAModelZPQQVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::TTbAModelZPQQVertex"; }
  
};

/** @endcond */
  
}


#endif /* HERWIG_TTbAModelZPQQVertex_H */
