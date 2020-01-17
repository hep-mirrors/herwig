// -*- C++ -*-
//
// SSNNZVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSNNZVertex_H
#define HERWIG_SSNNZVertex_H
//
// This is the declaration of the SSNNZVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "SusyBase.h"
#include "MixingMatrix.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the coupling of a Z-boson to a pair of neutralinos.
 * It inherits from FFVertex and implements the setCoupling method.
 *
 * @see \ref SSNNZVertexInterfaces "The interfaces"
 * defined for SSNNZVertex.
 * @see FFVertex
 */
class SSNNZVertex: public FFVVertex {

public:

  /**
   * The default constructor.
   */
  SSNNZVertex();

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
  SSNNZVertex & operator=(const SSNNZVertex &) = delete;

private:

  /**
   * The value of \f$sin(\theta_W)\f$
   */
  double _sw;
  
  /**
   * The value of \f$cos(\theta_W)\f$
   */
  double _cw;

  /**
   * Store the neutralino mixing matrix
   */
  tMixingMatrixPtr _theN;

  /**
   * Store the id of the 1st neutralino when the coupling was last evaluated
   */
  long _id1last;

  /**
   * Store the id of the 1st neutralino when the coupling was last evaluated
   */
  long _id2last;

  /**
   * Store the value at which the coupling when it was last evalutated
   */
  Energy2 _q2last;

  /**
   * Store the value of the coupling when it was last evalutated
   */
  Complex _couplast;

  /**
   * Store the value of the left-coupling when it was last evalutated
   */
  Complex _leftlast;

  /**
   * Store the value of the right-coupling when it was last evalutated
   */
  Complex _rightlast;
};
}

#endif /* HERWIG_SSNNZVertex_H */
