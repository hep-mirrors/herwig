// -*- C++ -*-
//
// SSGOGOHVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSGOGOHVertex_H
#define HERWIG_SSGOGOHVertex_H
//
// This is the declaration of the SSGOGOHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "MSSM.h"

namespace Herwig {

/**
 * The is the coupling of higgs bosons in the MSSM to a pair
 * of SM fermions.
 *
 * @see \ref SSGOGOHVertexInterfaces "The interfaces"
 * defined for SSGOGOHVertex.
 */
class SSGOGOHVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  SSGOGOHVertex();

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
  SSGOGOHVertex & operator=(const SSGOGOHVertex &) = delete;

private:

  /**
   * The mass of the \f$W\f$.
   */
  Energy theMw;

  /**
   * The matrix \f$S_{ij}\f$ 
   */
  vector<vector<Complex> > theSij;

  /**
   * The matrix \f$Q_{ij}\f$ 
   */
  vector<vector<Complex> > theQij;

  /**
   * The matrix \f$Q_{ij}^{L'}\f$ 
   */
  vector<vector<Complex> > theQijLp;

  /**
   * The matrix \f$Q_{ij}^{R'}\f$ 
   */
  vector<vector<Complex> > theQijRp;

  /**
   * The matrix \f$S_{ij}^{''}\f$ 
   */
  vector<vector<Complex> > theSijdp;

  /**
   * The matrix \f$Q_{ij}^{''}\f$ 
   */
  vector<vector<Complex> > theQijdp; 

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
   * The value of the coupling when it was last evaluated.
   */
  Complex theCoupLast;

  /**
   * The value of the left-coupling when it was last evaluated.
   */
  Complex theLLast;
  
  /**
   * The value of the right-coupling when it was last evaluated.
   */
  Complex theRLast;

  /**
   * The ID of the last higgs for which the vertex was evaluated
   */
  long theHLast;

  /**
   * The ID of the first gaugino when the coupling was las evaluated
   */
  long theID1Last;

  /**
   * The ID of the first gaugino when the coupling was las evaluated
   */
  long theID2Last;

  /**
   * The scale at which the coupling was last evaluated 
   */
  Energy2 theq2last;
};
}

#endif /* HERWIG_SSGOGOHVertex_H */
