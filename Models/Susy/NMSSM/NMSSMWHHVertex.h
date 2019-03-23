// -*- C++ -*-
#ifndef HERWIG_NMSSMWHHVertex_H
#define HERWIG_NMSSMWHHVertex_H
//
// This is the declaration of the NMSSMWHHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The NMSSMWHHVertex class implements the coupling of an electroweak"
 * gauge boson with two Higgs bosons in the NMSSM.
 *
 * @see \ref NMSSMWHHVertexInterfaces "The interfaces"
 * defined for NMSSMWHHVertex.
 */
class NMSSMWHHVertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  NMSSMWHHVertex();

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMWHHVertex & operator=(const NMSSMWHHVertex &) = delete;

private:

  /**
   * \f$\sin\beta\f$
   */
  double _sinb;

  /**
   * \f$\cos\beta\f$
   */
  double _cosb;

  /**
   * \f$\sin\theta_W\f$
   */
  double _sw;

  /**
   * \f$\cos\theta_W\f$
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

  /**
   *  Mixing matrix for the CP-even Higgs bosons
   */
  MixingMatrixPtr _mixS;

  /**
   *  Mixing matrix for the CP-odd  Higgs bosons
   */
  MixingMatrixPtr _mixP;

};
}


#endif /* HERWIG_NMSSMWHHVertex_H */
