// -*- C++ -*-
//
// LeptoquarkModelSLQFFVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LeptoquarkModelSLQFFVertex_H
#define HERWIG_LeptoquarkModelSLQFFVertex_H
//
// This is the declaration of the LeptoquarkModelSLQFFVertex class.

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig/Models/Leptoquarks/LeptoquarkModel.h"
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
class LeptoquarkModelSLQFFVertex: public FFSVertex {
  
public:
  
  /**
   * Default constructor.
   */
  LeptoquarkModelSLQFFVertex();
  
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
  static ClassDescription<LeptoquarkModelSLQFFVertex> initLeptoquarkModelSLQFFVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  LeptoquarkModelSLQFFVertex & operator=(const LeptoquarkModelSLQFFVertex &);

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
   *  Overall coupling to fermions
   */
  double _CFF;
  
  
  /**
   *  Overall coupling to left-handed leptons 
   */
  //double _cL;
  Complex _cL;

  /**
   *  Overall coupling to right-handed leptons 
   */
  //  double _cR;
  Complex _cR;

  /**
   *  Overall coupling to left-handed leptons for S0
   */
  // double _cL0;
  Complex _cL0;
  

  /**
   *  Overall coupling to right-handed leptons for S0
   */
  //  double _cR0;
  Complex _cR0;

  /**
   *  Overall coupling to right-handed leptons for ~S0
   */
  // double _cR0t;
  Complex _cR0t;
  
  /**
   *  Overall coupling to left-handed leptons for ~S1 triplet
   */
  // double _cL1;
  Complex _cL1;

  /**
   *  Overall coupling to left-handed leptons for S1/2 triplet
   */
  // double _cL12;
  Complex _cL12;
  

  /**
   *  Overall coupling to right-handed leptons for S1/2 triplet
   */
  // double _cR12;
  Complex _cR12;

  /**
   *  Overall coupling to left-handed leptons for ~S1/2 triplet
   */
  //  double _cL12t;
  Complex _cL12t;
  
  /**
   *  Overall coupling to left-handed leptons 
   */
  //double _dcL;
  Complex _dcL;

  /**
   *  Overall coupling to right-handed leptons 
   */
  // double _dcR;
  Complex _dcR;


  /**
   *  Overall coupling to left-handed leptons for dS0
   */
  //  double _dcL0;
  Complex _dcL0;

  /**
   *  Overall coupling to right-handed leptons for dS0
   */
  //  double _dcR0;
  Complex _dcR0;

  /**
   *  Overall coupling to right-handed leptons for ~dS0
   */
  double _dcR0t;

  /**
   *  Overall coupling to left-handed leptons for ~dS1 triplet
   */
  double _dcL1;

  /**
   *  Overall coupling to left-handed leptons for dS1/2 triplet
   */
  //double _dcL12;
  Complex _dcL12;

  /**
   *  Overall coupling to right-handed leptons for dS1/2 triplet
   */
  //double _dcR12;
  Complex _dcR12;


  /**
   *  Overall coupling to left-handed leptons for ~dS1/2 triplet
   */
  //  double _dcL12t;
  Complex _dcL12t;
  
  /**
   *  Suppression scale for derivatively coupled scalar leptoquarks
   */
  Energy _derivscale;





  //@}
};  

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of LeptoquarkModelSLQFFVertex.
 */
template <>
struct BaseClassTrait<Herwig::LeptoquarkModelSLQFFVertex,1> {
  /** Typedef of the base class of LeptoquarkModelSLQFFVertex. */
  typedef ThePEG::Helicity::FFSVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::LeptoquarkModelSLQFFVertex>
  : public ClassTraitsBase<Herwig::LeptoquarkModelSLQFFVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::LeptoquarkModelSLQFFVertex"; }
  
};

/** @endcond */
  
}


#endif /* HERWIG_LeptoquarkModelSLQFFVertex_H */
