// -*- C++ -*-
//
// SSWWHVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSWWHVertex_H
#define HERWIG_SSWWHVertex_H
//
// This is the declaration of the SSWWHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "MSSM.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the coupling of a higgs in the MSSM
 * to a pair of SM gauge bosons.
 *
 * @see \ref SSWWHVertexInterfaces "The interfaces"
 * defined for SSWWHVertex.
 */
class SSWWHVertex: public VVSVertex {

public:

  /**
   * The default constructor.
   */
  SSWWHVertex();

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSWWHVertex> initSSWWHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSWWHVertex & operator=(const SSWWHVertex &) = delete;

private:

  
  /**
   * The value of the factor \f$\frac{m_W \sin(\beta-\alpha)}{\sin\theta_W}\f$
   */
  Energy theh0Wfact;

  /**
   * The value of the factor \f$\frac{m_W \cos(\beta-\alpha)}{\sin\theta_W}\f$
   */
  Energy theH0Wfact;

  /**
   * The value of the factor 
   * \f$\frac{m_Z\sin(\beta-\alpha)}{\sin\theta_W\cos\theta_W}\f$
   */
  Energy theh0Zfact;

  /**
   * The value of the factor 
   * \f$\frac{m_Z\cos(\beta-\alpha)}{\sin\theta_W\cos\theta_W}\f$
   */
  Energy theH0Zfact;
    
  /**
   * The value of the coupling when it was last evaluated
   */
  complex<Energy> theCoupLast;

  /**
   * The value of \f$\sqrt{4\pi \alpha(q2)}\f$ when it was last evaluated.
   */
  Complex theElast;


  /**
   * The scale at which the coupling was last evaluated
   */
  Energy2 theq2last;

  /**
   * The ID of the last higgs for which the vertex was evaluated
   */
  long theHlast;

  /**
   * The ID of the last gauge boson for which the vertex was evaluated
   */
  long theGBlast;
};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSWWHVertex. */
template <>
struct BaseClassTrait<Herwig::SSWWHVertex,1> {
  /** Typedef of the first base class of SSWWHVertex. */
  typedef ThePEG::Helicity::VVSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSWWHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSWWHVertex>
  : public ClassTraitsBase<Herwig::SSWWHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSWWHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSWWHVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSWWHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSWWHVertex_H */
