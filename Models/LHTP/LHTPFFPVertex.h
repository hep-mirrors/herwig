// -*- C++ -*-
#ifndef HERWIG_LHTPFFPVertex_H
#define HERWIG_LHTPFFPVertex_H
//
// This is the declaration of the LHTPFFPVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHTPFFPVertex class implements the couplings of the photon
 * and its heavy partner \f$A_H\f$ in the Little Higgs model with T-parity.
 */
class LHTPFFPVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LHTPFFPVertex();
  
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHTPFFPVertex & operator=(const LHTPFFPVertex &) = delete;

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The charge of the Standard Model fermions.
   */
  vector<double> charge_;

  /**
   *  The last value of the coupling calculated.
   */
  Complex coupLast_;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 q2Last_;
  //@}

  /**
   *  Couplings of the fermion and T-fermions to the \f$A_H\f$
   */
  //@{
  /**
   *  Coupling of \f$dd_-A_H\f$
   */
  double coupd_;

  /**
   *  Coupling of \f$uu_-A_H\f$
   */
  double coupu_;

  /**
   *  Coupling of \f$ee_-A_H\f$
   */
  double coupe_;

  /**
   *  Coupling of \f$\nu\nu_-A_H\f$
   */
  double coupnu_;

  /**
   *  Coupling of heavy top
   */
  double TPreFactor_;

  /**
   *  Left mixings
   */
  double sL_,cL_;

  /**
   *  Right mixings
   */
  double sR_,cR_;
  //@}
};

}

#endif /* HERWIG_LHTPFFPVertex_H */
