// -*- C++ -*-
#ifndef RADIATIVEZPRIME_HiddenValleyFFZPrimeVertex_H
#define RADIATIVEZPRIME_HiddenValleyFFZPrimeVertex_H
//
// This is the declaration of the HiddenValleyFFZPrimeVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiddenValleyFFZPrimeVertex class.
 *
 * @see \ref HiddenValleyFFZPrimeVertexInterfaces "The interfaces"
 * defined for HiddenValleyFFZPrimeVertex.
 */
class HiddenValleyFFZPrimeVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  inline HiddenValleyFFZPrimeVertex();
  
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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
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
  HiddenValleyFFZPrimeVertex & operator=(const HiddenValleyFFZPrimeVertex &);


  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The left couplings of the Standard Model fermions.
   */
  vector<double> _gl;

  /**
   *  The right couplings of the Standard Model fermions.
   */
  vector<double> _gr;
  /**
   *  The left couplings of the new fermions.
   */
  vector<double> _gql;

  /**
   *  The right couplings of the new fermions.
   */
  vector<double> _gqr;

  /**
   *  Coupling
   */
  double _gPrime;
  //@}

};

}

#endif /* RADIATIVEZPRIME_HiddenValleyFFZPrimeVertex_H */
