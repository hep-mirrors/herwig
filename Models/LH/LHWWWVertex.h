// -*- C++ -*-
#ifndef HERWIG_LHWWWVertex_H
#define HERWIG_LHWWWVertex_H
//
// This is the declaration of the LHWWWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "LHModel.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::Direction;

/**
 * The LHWWWVertex class implements the triple boson coupling in the Little
 * Higgs model for the electroweak bosons of the Standard Model and the
 * additional \f$A_H\f$, \f$Z_H\f$ and \f$W_H^\pm\f$ bosons.
 */
class LHWWWVertex: public Helicity::VVVVertex {

public:

  /**
   * The default constructor.
   */
  inline LHWWWVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param d1 The direction for the first  particle.
   * @param d2 The direction for the second particle.
   * @param d3 The direction for the third  particle.
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
  LHWWWVertex & operator=(const LHWWWVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The correction factors for the different interacting particles
   */
  vector<double> _corr;

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

#endif /* HERWIG_LHWWWVertex_H */
