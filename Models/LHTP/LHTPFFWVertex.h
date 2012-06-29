// -*- C++ -*-
#ifndef HERWIG_LHTPFFWVertex_H
#define HERWIG_LHTPFFWVertex_H
//
// This is the declaration of the LHTPFFWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "LHTPFFWVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHTPFFWVertex class implements the coupling of the \f$W^\pm\f$
 * and \f$W^\pm_H\f$ bosons of the Little Higgs model with T-parity to fermions.
 * For simplicity the coupling are assumed to be flavour diagonal.
 *
 * @see \ref LHTPFFWVertexInterfaces "The interfaces"
 * defined for LHTPFFWVertex.
 */
class LHTPFFWVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LHTPFFWVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
  virtual void doinit();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LHTPFFWVertex> initLHTPFFWVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHTPFFWVertex & operator=(const LHTPFFWVertex &);

private:

  /**
   * @name Storage of the couplings.
   */
  //@{
  /**
   *  Top mixing angles
   */
  //@{
  /**
   *  \f$\sin\alpha\f$, taken as an input
   */
  double _salpha;

  /**
   *  \f$\cos\alpha\f$
   */
  double _calpha;

  /**
   *  \f$\sin\alpha\f$
   */
  double _sbeta;

  /**
   *  \f$\cos\beta\f$
   */
  double _cbeta;

  /**
   *  The elements of the CKM matrix.
   */
  vector<vector<Complex> > _ckm;

  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LHTPFFWVertex. */
template <>
struct BaseClassTrait<Herwig::LHTPFFWVertex,1> {
  /** Typedef of the first base class of LHTPFFWVertex. */
  typedef Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LHTPFFWVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LHTPFFWVertex>
  : public ClassTraitsBase<Herwig::LHTPFFWVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LHTPFFWVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LHTPFFWVertex is implemented. It may also include several, space-separated,
   * libraries if the class LHTPFFWVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLHTPModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LHTPFFWVertex_H */
