// -*- C++ -*-
#ifndef HERWIG_NMSSMGOGOHVertex_H
#define HERWIG_NMSSMGOGOHVertex_H
//
// This is the declaration of the NMSSMGOGOHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * Here is the documentation of the NMSSMGOGOHVertex class.
 *
 * @see \ref NMSSMGOGOHVertexInterfaces "The interfaces"
 * defined for NMSSMGOGOHVertex.
 */
class NMSSMGOGOHVertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  inline NMSSMGOGOHVertex();

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
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3);

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
  static ClassDescription<NMSSMGOGOHVertex> initNMSSMGOGOHVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMGOGOHVertex & operator=(const NMSSMGOGOHVertex &);

private:

  /**
   *  The various mixing matrices and couplings
   */
  //@{
  /**
   *  The V chargino mixing matrix
   */
  MixingMatrixPtr _mixV;

  /**
   *  The U chargino mixing matrix
   */
  MixingMatrixPtr _mixU;

  /**
   *  The neutralino mixing matrix 
   */
  MixingMatrixPtr _mixN;

  /**
   *  The CP-even neutral Higgs mixing matrix
   */
  MixingMatrixPtr _mixS;

  /**
   *  The CP-odd neutral Higgs mixing matrix
   */
  MixingMatrixPtr _mixP;

  /**
   *  The tri-linear \f$\lambda\f$ coupling
   */
  double _lambda;

  /**
   *  The tri-linear \f$\kappa\f$ coupling
   */
  double _kappa;

  /**
   *  \f$\sin\beta\f$
   */
  double _sinb;

  /**
   *  \f$\cos\beta\f$
   */
  double _cosb;

  /**
   *  \f$\sin\theta_W\f$
   */
  double _sw;

  /**
   *  \f$\cos\theta_W\f$
   */
  double _cw;

  /**
   *  The last \f$q^2\f$ the coupling was evaluated at.
   */
  Energy2 _q2last;

  /**
   *  The last value of the coupling
   */
  double _couplast;
  //@}

};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NMSSMGOGOHVertex. */
template <>
struct BaseClassTrait<Herwig::NMSSMGOGOHVertex,1> {
  /** Typedef of the first base class of NMSSMGOGOHVertex. */
  typedef ThePEG::Helicity::FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NMSSMGOGOHVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NMSSMGOGOHVertex>
  : public ClassTraitsBase<Herwig::NMSSMGOGOHVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NMSSMGOGOHVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NMSSMGOGOHVertex is implemented. It may also include several, space-separated,
   * libraries if the class NMSSMGOGOHVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwNMSSM.so"; }
};

/** @endcond */

}

#endif /* HERWIG_NMSSMGOGOHVertex_H */
