// -*- C++ -*-
#ifndef HERWIG_SSCNWVertex_H
#define HERWIG_SSCNWVertex_H
//
// This is the declaration of the SSCNWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "SSCNWVertex.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is implements the coupling of a W boson to a chargino and a 
 * neutralino. It inherits from FFVVertex and implements the setCoupling method.
 *
 * @see \ref SSCNWVertexInterfaces "The interfaces"
 * defined for SSCNWVertex.
 * @see FFVVertex 
 */
class SSCNWVertex: public FFVVertex {

public:

  /**
   * The default constructor.
   */
  SSCNWVertex();

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SSCNWVertex> initSSCNWVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSCNWVertex & operator=(const SSCNWVertex &);

private:

  /**
   * A pointer to the Susy Model object
   */
  tSusyBasePtr _theSS;

  /**
   * Store \f$sin(\theta_W)\f$
   */
  double _sw;

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
};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SSCNWVertex. */
template <>
struct BaseClassTrait<Herwig::SSCNWVertex,1> {
  /** Typedef of the first base class of SSCNWVertex. */
  typedef ThePEG::Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSCNWVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SSCNWVertex>
  : public ClassTraitsBase<Herwig::SSCNWVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SSCNWVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SSCNWVertex is implemented. It may also include several, space-separated,
   * libraries if the class SSCNWVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#include "SSCNWVertex.icc"

#endif /* HERWIG_SSCNWVertex_H */
