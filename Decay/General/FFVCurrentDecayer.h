// -*- C++ -*-
//
// FFVCurrentDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FFVCurrentDecayer_H
#define HERWIG_FFVCurrentDecayer_H
//
// This is the declaration of the FFVCurrentDecayer class.
//

#include "GeneralCurrentDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::FFVVertexPtr;

/**
 * Here is the documentation of the FFVCurrentDecayer class.
 *
 * @see \ref FFVCurrentDecayerInterfaces "The interfaces"
 * defined for FFVCurrentDecayer.
 */
class FFVCurrentDecayer: public GeneralCurrentDecayer {

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const tPDVector & outgoing,
	     const vector<Lorentz5Momentum> & momenta,
	     MEOption meopt) const;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const;
  
  /**
   * Function to return partial Width
   * @param inpart Pointer to incoming particle data object
   * @param outa Pointer to first outgoing particle data object
   * @param currout The outgoing particles from the current
   */
  virtual Energy partialWidth(tPDPtr inpart, tPDPtr outa,
			      vector<tPDPtr> currout);
  //@}

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FFVCurrentDecayer & operator=(const FFVCurrentDecayer &) = delete;

private:
  
  /**
   * Pointer to FFVVertex
   */
  FFVVertexPtr FFVPtr_;

  /**
   *  Spinr density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Spinor wavefunction
   */
  mutable vector<SpinorWaveFunction>    wave_   ;

  /**
   *  Barred spinor wavefunction
   */
  mutable vector<SpinorBarWaveFunction> wavebar_;
};

}

#endif /* HERWIG_FFVCurrentDecayer_H */
