// -*- C++ -*-
#ifndef HERWIG_SMFFHVertex_H
#define HERWIG_SMFFHVertex_H
//
// This is the declaration of the SMFFHVertex class.

#include "Herwig++/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the vertex coupling the Standard Model Higgs
 *  to the Standard Model fermions for helicity amplitude calculations
 *
 *  @see FFSVertex
 *  @see VertexBase
 */
class SMFFHVertex: public FFSVertex {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline SMFFHVertex();

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param ioff Which particle is off-shell
  */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3, int ioff);

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
  
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SMFFHVertex> initSMFFHVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMFFHVertex & operator=(const SMFFHVertex &);

  /**
   * Pointer to the SM object.
   */
  tcHwSMPtr _theSM;

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  Last evaluation of the coupling
   */
  Complex _couplast;

  /**
   *  The value of \f$\sin\theta_w\f$.
   */
  double _sw;

  /**
   *  The PDG code of the last fermion the coupling was evaluated for.
   */
  int _idlast;

  /**
   *  The last \f$q^2\f$ the coupling was evaluated at.
   */
  Energy2 _q2last;

  /**
   * The mass of the last fermion for which the coupling was evaluated.
   */
  Energy _masslast;

  /**
   *  The mass of the \f$W\f$ boson.
   */
  Energy _mw;
  //@}
};  

}
} 
#include "SMFFHVertex.icc"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SMFFHVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::SMFFHVertex,1> {
  /** Typedef of the base class of SMFFHVertex. */
  typedef Herwig::Helicity::FFSVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::SMFFHVertex>
  : public ClassTraitsBase<Herwig::Helicity::SMFFHVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::SMFFHVertex"; }
  
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwSMVertex.so"; }
  
};

/** @endcond */
  
}


#endif /* HERWIG_SMFFHVertex_H */
