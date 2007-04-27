// -*- C++ -*-
#ifndef HERWIG_SSZSSVertex_H
#define HERWIG_SSZSSVertex_H
//
// This is the declaration of the SSZSSVertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/VSSVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "SSZSSVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
/**
 * This is the implementation of the \f$Z^0\f$ coupling to 2 sfermions.
 * It inherits from VSSVertex and implements the setCoupling() method.
 *
 * @see VSSVertex
 */
class SSZSSVertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  SSZSSVertex();

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
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
                           tcPDPtr part2,tcPDPtr part3);

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
   * Initialize this object after the setup phase before saving and
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
  static ClassDescription<SSZSSVertex> initSSZSSVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SSZSSVertex & operator=(const SSZSSVertex &);

  /**
   * Pointer to the Standard Model object
   */
  tSusyBasePtr _theSS;

  /**
   * Value of \f$sin(\theta_w)\f$
   */
  double _sw;
  
  /**
   * Value of \f$cos(\theta_w)\f$
   */
  double _cw;

  /**
   * Value of factor including mixing matrices when last evaluated
   */
  Complex _zfact;
  
  /**
   * Stop mixing matrix
   */
  tMixingMatrixPtr _stop;

  /**
   * Sbottom mixing matrix
   */
  tMixingMatrixPtr _sbottom;
  
  /**
   * Stau mixing matrix
   */
  tMixingMatrixPtr _stau;

  /**
   * Value of coupling when last evaluated
   */
  Complex _couplast;

  /**
   * Scale at which the coupling was last evaluated
   */
  Energy2 _q2last;

  /**
   * id of sfermion that coupling was last evaluated for 
   */
  long _idlast;
};
}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/// \if TRAITSPECIALIZATIONS

/** This template specialization informs ThePEG about the
 *  base classes of SSZSSVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::SSZSSVertex,1> {
  /** Typedef of the first base class of SSZSSVertex. */
  typedef Herwig::Helicity::VSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SSZSSVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::SSZSSVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSZSSVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SSZSSVertex"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SSZSSVertex class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSusy.so HwSusyVertex.so"; }
};

/// \endif

}

#include "SSZSSVertex.icc"

#endif /* HERWIG_SSZSSVertex_H */
