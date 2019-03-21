// -*- C++ -*-
#ifndef HERWIG_LHFFHVertex_H
#define HERWIG_LHFFHVertex_H
//
// This is the declaration of the LHFFHVertex class.
//

#include "LHModel.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHFFHVertex class implements the couplings of the fermions to 
 * the Higgs bosons of the Little Higgs model.
 */
class LHFFHVertex: public Helicity::FFSVertex {

public:

  /**
   * The default constructor.
   */
  LHFFHVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param ioff Which particle is off-shell
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
  LHFFHVertex & operator=(const LHFFHVertex &) = delete;

private:

  /**
   * Pointer to the SM object.
   */
  tcLHModelPtr _model;

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The PDG code of the last fermion the coupling was evaluated for.
   */
  int _idlast[2];

  /**
   *  The last \f$q^2\f$ the coupling was evaluated at.
   */
  Energy2 _q2last;

  /**
   * The mass of the last fermion for which the coupling was evaluated.
   */
  Energy _masslast[2];

  /**
   *  The factors for the individual interactions
   */
  vector<complex<InvEnergy> > _coup;
  //@}
};

}

#endif /* HERWIG_LHFFHVertex_H */
