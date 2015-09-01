// -*- C++ -*-
//
// SSCFSVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSCFSVertex_H
#define HERWIG_SSCFSVertex_H
//
// This is the declaration of the SSCFSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "MSSM.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This is the implementation of the coupling of the charginos to fermion-
 * sfermions. It inherits from FFSVertex and implements the virtual 
 * setCoupling member
 *
 * @see FFSVertex
 */
class SSCFSVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  SSCFSVertex();

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
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
                           tcPDPtr part2,tcPDPtr part3);
  
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSCFSVertex & operator=(const SSCFSVertex &);

  /**
   * Pointer to the MSSM object
   */
  tMSSMPtr _theSS;

  /**
   * \f$\sin(\beta)\f$
   */
  double _sb;

  /**
   * \f$\cos(\beta)\f$
   */
  double _cb;

  /**
   * Mass of the W
   */
  Energy _mw;
  
  /**
   * Pointer to stop mixing matrix
   */
  tMixingMatrixPtr _stop;
  
  /**
   * Pointer to sbottom mixing matrix
   */
  tMixingMatrixPtr _sbot;

  /**
   * Pointer to stau mixing matrix 
   */
  tMixingMatrixPtr _stau;

  /**
   * Pointer to U chargino mixing matrix 
   */
  tMixingMatrixPtr _umix;

  /**
   * Pointer to V chargino mixing matrix
   */
  tMixingMatrixPtr _vmix;

  /**
   * The energy scale at which the coupling 
   * was last evaluated 
   */
  Energy2 _q2last;

  /**
   * The value of the coupling at the scale last evaluated
   */
  Complex _couplast;

  /**
   * The value of the left coupling when last evaluated
   */
  Complex _leftlast;

  /**
   * The value of the right coupling when last evaluated
   */
  Complex _rightlast;

  /**
   * Id of the last sm fermion that the coupling was evaluated for 
   */
  long _id1last; 

  /**
   * Id of the last scalar that the coupling was evaluated for 
   */
  long _id2last; 
  
  /**
   * Id of the last chargino that the coupling was evaluated for 
   */
  long _id3last; 

  /**
   *  Include Yukawa's ?
   */
  unsigned int yukawa_;
};
}

#endif /* HERWIG_SSCFSVertex_H */
