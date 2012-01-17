// -*- C++ -*-
//
// SSCCZVertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSCCZVertex_H
#define HERWIG_SSCCZVertex_H
//
// This is the declaration of the SSCCZVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "SusyBase.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the coupling of a \f$\gamma/Z^0\f$ to a pair of 
 * charginos. It inherits from FFVVertex and implements the setCoupling method.
 *
 * @see \ref SSCCZVertexInterfaces "The interfaces"
 * defined for SSCCZVertex.
 * @see FFVVertex
 */
class SSCCZVertex: public FFVVertex {

public:

  /**
   * The default constructor.
   */
  SSCCZVertex();

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
  SSCCZVertex & operator=(const SSCCZVertex &);

private:

  /**
   * Value of \f$sin^2(\theta_W)\f$
   */
  double _sw2;

  /**
   * Value of \f$cos(\theta_W)\f$
   */
  double _cw;

  /**
   * The U mixing matrix
   */
  tMixingMatrixPtr _theU;

  /**
   * The V mixing matrix
   */
  tMixingMatrixPtr _theV;

  /**
   * The value of the coupling when it was last evaluated
   */
  Complex _couplast;
  
  /**
   * The scale at which the coupling was last evaluated
   */
  Energy2 _q2last;

  /**
   * The id of the first chargino the last time the vertex was evaluated
   */
  long _id1last;

  /**
   * The id of the second chargino the last time the vertex was evaluated
   */
  long _id2last;

  /**
   * The value of the left coupling when it was last evaluated
   */
  Complex _leftlast;

  /**
   * The value of the right coupling when it was last evaluated
   */
  Complex _rightlast;

  /**
   * The ID of the gauge boson when the vertex was last evaluated
   */
  long _gblast;
};
}

#endif /* HERWIG_SSCCZVertex_H */
