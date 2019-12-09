// -*- C++ -*-
//
// SSNFSVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSNFSVertex_H
#define HERWIG_SSNFSVertex_H
//
// This is the declaration of the SSNFSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "MSSM.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This is the implementation of the coupling of the neutralinos to 
 * fermion-sfermions. It inherits from FFSVertex and implements the 
 * virtual setCoupling() method.
 * 
 * @see FFSVertex
 */
class SSNFSVertex: public FFSVertex {

public:

   /**
   * The default constructor.
   */
  SSNFSVertex();

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1,
                           tcPDPtr part2, tcPDPtr part3);

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSNFSVertex & operator=(const SSNFSVertex &) = delete;

  /**
   * Pointer to the stop mixing matrix
   */
  tMixingMatrixPtr _stop;

  /**
   * Pointer to the sbottom mixing matrix
   */
  tMixingMatrixPtr _sbot;

  /**
   * Pointer to the stau mixing matrix
   */
  tMixingMatrixPtr _stau;
  
  /**
   * Pointer to the neutralino mixing matrix
   */
  tMixingMatrixPtr _nmix;
 
  /**
   * Pointer to the Susy Model object
   */
  tMSSMPtr _theSS;

  /**
   * \f$\sin(\theta_w)\f$
   */
  double _sw;

  /**
   * \f$\cos(\theta_w)\f$
   */
  double _cw;
  
  /**
   * Mass of the W
   */
  Energy _mw;

  /**
   * \f$\sin(\beta)\f$
   */
  double _sb;
  
  /**
   * \f$\cos(\beta)\f$
   */
  double _cb;

  /**
   * The scale at which the coupling was last evaluated. 
   */
  Energy2 _q2last;

  /**
   * The value of the normalisation when it was evaluated at _q2last 
   */
  Complex _couplast;
  
  /**
   * Store the value of the left coupling when it was last evaluated
   */
  Complex _leftlast;

  /**
   * Store the value of the right coupling when it was last evaluated
   */
  Complex _rightlast;

  /**
   * Store the id of the last neutralino to be evaluate
   */
  long _id1last;
  
  /**
   * Store the id of the last SM fermion to be evaluate
   */
  long _id2last;

  /**
   * Store the id of the last scalar to be evaluate
   */
  long _id3last;

  /**
   *  Include Yukawa's ?
   */
  unsigned int yukawa_;
};
}

#endif /* HERWIG_SSNFSVertex_H */
