// -*- C++ -*-
#ifndef HERWIG_LHTPWWWVertex_H
#define HERWIG_LHTPWWWVertex_H
//
// This is the declaration of the LHTPWWWVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::Direction;

/**
 * This is the coupling of vector bosons to each other in the
 * Littlest Higgs model with T-parity. It inherits from the 
 * Standard Model WWW vertex to use its setCoupling member
 * for the SM gauge boson self-couplings.
 */
class LHTPWWWVertex: public VVVVertex {

public:

  /**
   * The default constructor.
   */
  LHTPWWWVertex();

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
   * @param a The ParticleData pointer for the first  particle.
   * @param b The ParticleData pointer for the second particle.
   * @param c The ParticleData pointer for the third  particle.
   * @param d1 The direction for the first  particle.
   * @param d2 The direction for the second particle.
   * @param d3 The direction for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
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
  LHTPWWWVertex & operator=(const LHTPWWWVertex &) = delete;

private:
  
  /**
   * The value of the coupling when it was last evaluated.
   */
  Complex coupLast_;

  /**
   * The scale where the coulpling was last evaluated 
   */
  Energy2 q2Last_;

  /**
   * The couplings for the various possible interactions
   */
  vector<double> couplings_;
};

}

#endif /* HERWIG_LHTPWWWVertex_H */
