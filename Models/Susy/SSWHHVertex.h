// -*- C++ -*-
//
// SSWHHVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSWHHVertex_H
#define HERWIG_SSWHHVertex_H
//
// This is the declaration of the SSWHHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "MSSM.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The coupling of a pair of higgs to the SM gauge bosons
 * in the MSSM.
 *
 * @see \ref SSWHHVertexInterfaces "The interfaces"
 * defined for SSWHHVertex.
 */
class SSWHHVertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  SSWHHVertex();

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
   * Calculate the coupling for the vertex
   * @param q2 The scale to at which evaluate the coupling.
   * @param particle1 The first particle in the vertex.
   * @param particle2 The second particle in the vertex.
   * @param particle3 The third particle in the vertex.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3);

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
  SSWHHVertex & operator=(const SSWHHVertex &);

private:
  
  /**
   * The value of \f$\sin\theta_W\f$ 
   */
  double theSw;

  /**
   * The value of \f$\sin2\theta_W\f$ 
   */
  double theS2w;

  /**
   * The value of \f$\cos2\theta_W\f$ 
   */
  double theC2w;

  /**
   * The value of \f$\sin(\beta - \alpha)\f$ 
   */
  double thesbma;

  /**
   * The value of \f$\cos(\beta - \alpha)\f$ 
   */
  double thecbma;
  
  /**
   * The scale at which the coupling  was last evaluated.
   */
  Energy2 theq2last;
  
  /**
   * The value of the \f$\sqrt{4\pi\alpha}\f$  when last evaluated.
   */
  double theElast;

};

}

#endif /* HERWIG_SSWHHVertex_H */
