// -*- C++ -*-
#ifndef HERWIG_NMSSMWWHHVertex_H
#define HERWIG_NMSSMWWHHVertex_H
//
// This is the declaration of the NMSSMWWHHVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VVSSVertex.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the NMSSMWWHHVertex class.
 *
 * @see \ref NMSSMWWHHVertexInterfaces "The interfaces"
 * defined for NMSSMWWHHVertex.
 */
class NMSSMWWHHVertex: public Helicity::VVSSVertex {

public:

  /**
   * The default constructor.
   */
  NMSSMWWHHVertex();
  
  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3, tcPDPtr part4);

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
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSMWWHHVertex & operator=(const NMSSMWWHHVertex &) = delete;

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex couplast_;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 q2last_;
  //@}

  /**
   * \f$\sin\theta_W\f$
   */
  double sw_;

  /**
   * \f$\cos\theta_W\f$
   */
  double cw_;

  /**
   * \f$\sin\beta\f$
   */
  double sb_;

  /**
   * \f$\cos\beta\f$
   */
  double cb_;

  /**
   * The CP-even Higgs mixing matrix 
   */
  MixingMatrixPtr mixS_;

  /**
   * The CP-odd Higgs mixing matrix 
   */
  MixingMatrixPtr mixP_;
};

}

#endif /* HERWIG_NMSSMWWHHVertex_H */
