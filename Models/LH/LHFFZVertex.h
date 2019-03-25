// -*- C++ -*-
#ifndef HERWIG_LHFFZVertex_H
#define HERWIG_LHFFZVertex_H
//
// This is the declaration of the LHFFZVertex class.
//

#include "LHModel.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

namespace Herwig {

/**
 * The LHFFZVertex class implements the couplings of the fermion,
 *  both of the Standard Model and the additional heavy top to the 
 *  \f$Z^0\f$ and \f$Z_H\f$ bosons in the Little Higgs model.
 */
class LHFFZVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LHFFZVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LHFFZVertex & operator=(const LHFFZVertex &) = delete;

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The left couplings of the Standard Model fermions to the Z.
   */
  vector<double> _gl;

  /**
   *  The right couplings of the Standard Model fermions to the Z.
   */
  vector<double> _gr;

  /**
   *  The left couplings of the Standard Model fermions to the heavy Z.
   */
  vector<double> _glH;

  /**
   *  The right couplings of the Standard Model fermions to the heavy Z.
   */
  vector<double> _grH;

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

#endif /* HERWIG_LHFFZVertex_H */
