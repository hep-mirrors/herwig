// -*- C++ -*-
//
// SSFFHVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSFFHVertex_H
#define HERWIG_SSFFHVertex_H
//
// This is the declaration of the SSFFHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "MSSM.h"

namespace Herwig {

/**
 * The is the coupling of higgs bosons in the MSSM to a pair
 * of SM fermions.
 *
 * @see \ref SSFFHVertexInterfaces "The interfaces"
 * defined for SSFFHVertex.
 */
class SSFFHVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  SSFFHVertex();

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
  SSFFHVertex & operator=(const SSFFHVertex &);

private:

  /**
   * Pointer to the SM object.
   */
  tcMSSMPtr theMSSM;

  /**
   * The value of \f$\tan\beta\f$.
   */
  double thetanb;

  /**
   * The mass of the \f$W\f$.
   */
  Energy theMw;

  /**
   * The value of \f$\sin\theta_W\f$
   */
  double theSw;

  /**
   * The value of \f$\sin\alpha\f$ 
   */
  double theSa;

  /**
   * The value of \f$\sin\beta\f$ 
   */
  double theSb;

  /**
   * The value of \f$\cos\alpha\f$ 
   */
  double theCa;

  /**
   * The value of \f$\cos\beta\f$ 
   */
  double theCb; 

  /**
   * The ID of the last fermion for which the vertex was evaluated
   */
  pair<long,long> theFLast;
    
  /**
   * The value of \f$ \frac{e}{\sin\theta_W} \f$ when it was last evaluated.
   */
  double theGlast;

    /**
   * The scale at which then coupling was last evaluated. 
   */
  Energy2 theq2last;

  /**
   *  Values of the masses
   */
  pair<Energy,Energy> theMassLast;
};
}

#endif /* HERWIG_SSFFHVertex_H */
