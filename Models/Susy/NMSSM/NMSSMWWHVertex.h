// -*- C++ -*-
#ifndef HERWIG_NMSSMWWHVertex_H
#define HERWIG_NMSSMWWHVertex_H
//
// This is the declaration of the NMSSMWWHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "Herwig++/Models/Susy/MixingMatrix.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Helicity
 * 
 * The NMSSMWWHVertex class is the implementation of the coupling of two electroweak
 * gauge bosons to the Higgs bosons of the NMSSM. It inherits from VVSVertex and
 * implements the setCoupling member.
 *
 * @see \ref NMSSMWWHVertexInterfaces "The interfaces"
 * @see VVSVertex
 * @see SMWWHVertex
 * defined for NMSSMWWHVertex.
 */
class NMSSMWWHVertex: public VVSVertex {

public:

  /**
   * The default constructor.
   */
  NMSSMWWHVertex();

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
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<NMSSMWWHVertex> initNMSSMWWHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMWWHVertex & operator=(const NMSSMWWHVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;

  /**
   *  The mass of the \f$W\f$ boson.
   */
  Energy _mw;

  /**
   *  The factor for the \f$Z\f$ vertex.
   */
  double _zfact;

  /**
   *  \f$\sin\beta\f$
   */
  double _sinb;

  /**
   *  \f$\cos\beta\f$
   */
  double _cosb;

  /**
   *  Mixing matrix for the CP-even Higgs bosons
   */
  MixingMatrixPtr _mixS;
  //@}

};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NMSSMWWHVertex. */
template <>
struct BaseClassTrait<Herwig::NMSSMWWHVertex,1> {
  /** Typedef of the first base class of NMSSMWWHVertex. */
  typedef ThePEG::Helicity::VVSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NMSSMWWHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NMSSMWWHVertex>
  : public ClassTraitsBase<Herwig::NMSSMWWHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NMSSMWWHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NMSSMWWHVertex is implemented. It may also include several, space-separated,
   * libraries if the class NMSSMWWHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwNMSSM.so"; }
};

/** @endcond */

}

#endif /* HERWIG_NMSSMWWHVertex_H */
