// -*- C++ -*-
#ifndef HERWIG_NMSSMFFHVertex_H
#define HERWIG_NMSSMFFHVertex_H
//
// This is the declaration of the NMSSMFFHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Helicity
 *
 * The NMSSMFFHVertex class implements the interactions of the NMSSM Higgs bosons
 * with the  Standard Model fermions.
 *
 * @see \ref NMSSMFFHVertexInterfaces "The interfaces"
 * defined for NMSSMFFHVertex.
 */
class NMSSMFFHVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  NMSSMFFHVertex();

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
  static ClassDescription<NMSSMFFHVertex> initNMSSMFFHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMFFHVertex & operator=(const NMSSMFFHVertex &) = delete;

private:

  /**
   *  Mixing matrix for the CP-even Higgs bosons
   */
  MixingMatrixPtr _mixS;

  /**
   *  Mixing matrix for the CP-odd  Higgs bosons
   */
  MixingMatrixPtr _mixP;

  /**
   * Pointer to the SM object.
   */
  tcHwSMPtr _theSM;

  /**
   *  Mass of the \f$W\f$ boson
   */
  Energy _mw;

  /**
   *  \f$\sin\beta\f$
   */
  double _sinb;

  /**
   *  \f$\cos\beta\f$
   */
  double _cosb;

  /**
   *  \f$\tan\beta\f$
   */
  double _tanb;

  /**
   *  \f$\sin\theta_W\f$
   */
  double _sw;

  /**
   *  The PDG code of the last fermion the coupling was evaluated for.
   */
  pair<int,int> _idlast;

  /**
   *  The last \f$q^2\f$ the coupling was evaluated at.
   */
  Energy2 _q2last;

  /**
   * The mass of the last fermion for which the coupling was evaluated.
   */
  pair<Energy,Energy> _masslast;

  /**
   *  The last value of the coupling
   */
  double _couplast;

};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NMSSMFFHVertex. */
template <>
struct BaseClassTrait<Herwig::NMSSMFFHVertex,1> {
  /** Typedef of the first base class of NMSSMFFHVertex. */
  typedef ThePEG::Helicity::FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NMSSMFFHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NMSSMFFHVertex>
  : public ClassTraitsBase<Herwig::NMSSMFFHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NMSSMFFHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NMSSMFFHVertex is implemented. It may also include several, space-separated,
   * libraries if the class NMSSMFFHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwNMSSM.so"; }
};

/** @endcond */

}

#endif /* HERWIG_NMSSMFFHVertex_H */
