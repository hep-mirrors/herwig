// -*- C++ -*-
//
// SSNNPVertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSNNPVertex_H
#define HERWIG_SSNNPVertex_H
//
// This is the declaration of the SSNNPVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/GeneralFFVVertex.h"
#include "Herwig++/Models/Susy/MSSM.h"
#include "Herwig++/Models/Susy/MixingMatrix.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the coupling of a photon to a pair of neutralinos
 * via loop diagrams
 * It inherits from GeneralFFVertex and implements the setCoupling method.
 *
 * @see \ref SSNNPVertexInterfaces "The interfaces"
 * defined for SSNNPVertex.
 * @see FFVertex
 */
class SSNNPVertex: public GeneralFFVVertex {

public:

  /**
   * The default constructor.
   */
  SSNNPVertex();

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

  /**
   *  Evaluate the loop integrals
   */
  void loopIntegrals(Energy Mi, Energy Mj, Energy M, Energy m,
		     complex<InvEnergy2> & I, complex<InvEnergy2> & J,
		     complex<InvEnergy2> & K, complex<InvEnergy2> & I2);

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
  SSNNPVertex & operator=(const SSNNPVertex &);

private:

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
   * The value of \f$sin(\theta_W)\f$
   */
  double _sw;
  
  /**
   * The value of \f$cos(\theta_W)\f$
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
   * Store the neutralino mixing matrix
   */
  tMixingMatrixPtr _theN;

  /**
   * Store the U-type chargino mixing matrix
   */
  tMixingMatrixPtr _theU;

  /**
   * Store the V-type chargino mixing matrix
   */
  tMixingMatrixPtr _theV;

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
  complex<InvEnergy> _leftlast;

  /**
   * Store the value of the right-coupling when it was last evalutated
   */
  complex<InvEnergy> _rightlast;
};
}

#endif /* HERWIG_SSNNPVertex_H */
