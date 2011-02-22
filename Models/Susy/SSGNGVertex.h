// -*- C++ -*-
//
// SSGNGVertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSGNGVertex_H
#define HERWIG_SSGNGVertex_H
//
// This is the declaration of the SSGNGVertex class.
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
 * @see \ref SSGNGVertexInterfaces "The interfaces"
 * defined for SSGNGVertex.
 * @see FFVertex
 */
class SSGNGVertex: public GeneralFFVVertex {

public:

  /**
   * The default constructor.
   */
  SSGNGVertex();

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSGNGVertex> initSSGNGVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSGNGVertex & operator=(const SSGNGVertex &);

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
   * Store the id of the neutralino when the coupling was last evaluated
   */
  long _idlast;

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


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSGNGVertex. */
template <>
struct BaseClassTrait<Herwig::SSGNGVertex,1> {
  /** Typedef of the first base class of SSGNGVertex. */
  typedef ThePEG::Helicity::GeneralFFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSGNGVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSGNGVertex>
  : public ClassTraitsBase<Herwig::SSGNGVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SSGNGVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSGNGVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSGNGVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SSGNGVertex_H */
