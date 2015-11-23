// -*- C++ -*-
//
// ZprimeModelZPQQVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ZprimeModelZPQQVertex_H
#define HERWIG_ZprimeModelZPQQVertex_H
//
// This is the declaration of the ZprimeModelZPQQVertex class.

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/Models/Zprime/ZprimeModel.h"
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
class ZprimeModelZPQQVertex: public FFVVertex {
  
public:
  
  /**
   * Default constructor.
   */
  ZprimeModelZPQQVertex();
  
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
  static ClassDescription<ZprimeModelZPQQVertex> initZprimeModelZPQQVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  ZprimeModelZPQQVertex & operator=(const ZprimeModelZPQQVertex &);

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
   *  Z prime coupling to top-anti-top (left-handed)
   */
  double _cZPTT_L;
  

  /**
   *  Z prime coupling to top-anti-top (right-handed)
   */
  double _cZPTT_R;

  /**
   *  Z prime coupling to u-ubar (left-handed)
   */
  double _cZPUU_L;
  

  /**
   *  Z prime coupling to u-ubar (right-handed)
   */
  double _cZPUU_R;

  /**
   *  Z prime coupling to c-cbar (left-handed)
   */
  double _cZPCC_L;
  

  /**
   *  Z prime coupling to c-cbar (right-handed)
   */
  double _cZPCC_R;

    /**
   *  Z prime coupling to s-sbar (left-handed)
   */
  double _cZPSS_L;
  

  /**
   *  Z prime coupling to s-sbar (right-handed)
   */
  double _cZPSS_R;
  
 /**
   *  Z prime coupling to d-dbar (left-handed)
   */
  double _cZPDD_L;

  /**
   *  Z prime coupling to d-dbar (right-handed)
   */
  double _cZPDD_R;

  /**
   *  Z prime coupling to d-dbar (left-handed)
   */
  double _cZPBB_L;
  
  /**
   *  Z prime coupling to d-dbar (right-handed)
   */
  double _cZPBB_R;

    /**
   *  Z prime coupling to e+e- (left-handed)
   */
  double _cZPee_L;
  
  /**
   *  Z prime coupling to e+e- (right-handed)
   */
  double _cZPee_R;
  
 /**
   *  Z prime coupling to mu+mu- (left-handed)
   */
  double _cZPmm_L;

  /**
   *  Z prime coupling to mu+mu- (right-handed)
   */
  double _cZPmm_R;

  /**
   *  Z prime coupling to tau+tau- (left-handed)
   */
  double _cZPtt_L;
  
  /**
   *  Z prime coupling to tau+tau- (right-handed)
   */
  double _cZPtt_R;


    /**
   *  Z prime coupling to nu_e nu_ebar (left-handed)
   */
  double _cZPnuenue_L;
  
  /**
   *  Z prime coupling to nu_e nu_ebar (right-handed)
   */
  double _cZPnuenue_R;
  
 /**
   *  Z prime coupling to nu_mu nu_mubar (left-handed)
   */
  double _cZPnumnum_L;

  /**
   *  Z prime coupling to nu_mu nu_mubar (right-handed)
   */
  double _cZPnumnum_R;

  /**
   *  Z prime coupling to nu_tau nu_taubar (left-handed)
   */
  double _cZPnutnut_L;
  
  /**
   *  Z prime coupling to nu_tau nu_taubar (right-handed)
   */
  double _cZPnutnut_R;
   
  /**
   *  Z prime overall coupling
   */
  double _cZP_o;


  //@}
};  

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ZprimeModelZPQQVertex.
 */
template <>
struct BaseClassTrait<Herwig::ZprimeModelZPQQVertex,1> {
  /** Typedef of the base class of ZprimeModelZPQQVertex. */
  typedef ThePEG::Helicity::FFVVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
  template <>
  
struct ClassTraits<Herwig::ZprimeModelZPQQVertex>
  : public ClassTraitsBase<Herwig::ZprimeModelZPQQVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::ZprimeModelZPQQVertex"; }
  
};

/** @endcond */
  
}


#endif /* HERWIG_ZprimeModelZPQQVertex_H */
