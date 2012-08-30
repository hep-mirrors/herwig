// -*- C++ -*-
#ifndef HERWIG_LHTPFFWVertex_H
#define HERWIG_LHTPFFWVertex_H
//
// This is the declaration of the LHTPFFWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHTPFFWVertex class implements the coupling of the \f$W^\pm\f$
 * and \f$W^\pm_H\f$ bosons of the Little Higgs model with T-parity to fermions.
 * For simplicity the coupling are assumed to be flavour diagonal.
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
   *  Left mixings
   */
  double sL_,cL_;

  /**
   *  The elements of the CKM matrix.
   */
  vector<vector<Complex> > ckm_;

  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex coupLast_;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 q2Last_;
  //@}
};

}

#endif /* HERWIG_LHTPFFWVertex_H */
