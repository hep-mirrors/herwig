// -*- C++ -*-
//
// TTbAModelWPTDVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TTbAModelWPTDVertex_H
#define HERWIG_TTbAModelWPTDVertex_H
//
// This is the declaration of the TTbAModelWPTDVertex class.

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
class TTbAModelWPTDVertex: public FFVVertex {
  
public:
  
  /**
   * Default constructor.
   */
  TTbAModelWPTDVertex();
  
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
   * Private and non-existent assignment operator.
   */
  TTbAModelWPTDVertex & operator=(const TTbAModelWPTDVertex &) = delete;

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
   *  W prime coupling to top-down (left-handed)
   */
  double _cWPTD_L;
  

  /**
   *  W prime coupling to top-down (right-handed)
   */
  double _cWPTD_R;

  /**
   * Model selector
   */
  int _models;


  //@}
};  

}

#endif /* HERWIG_TTbAModelWPTDVertex_H */
