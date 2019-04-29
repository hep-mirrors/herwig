// -*- C++ -*-
#ifndef HERWIG_SextetFFSVertex_H
#define HERWIG_SextetFFSVertex_H
//
// This is the declaration of the SextetFFSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SextetFFSVertex class.
 *
 * @see \ref SextetFFSVertexInterfaces "The interfaces"
 * defined for SextetFFSVertex.
 */
  class SextetFFSVertex: public Helicity::FFSVertex {

public:

  /**
   * The default constructor.
   */
  SextetFFSVertex() {
    colourStructure(ColourStructure::SU3K6);
  }

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1,
                           tcPDPtr part2, tcPDPtr part3);

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
  SextetFFSVertex & operator=(const SextetFFSVertex &) = delete;

private:

  /**
   *   The various couplings
   */
  //@{
  /**
   * The \f$SU(2)\f$ quark-doublet coupling to \f$\Phi_{6,1,1/3}\f$
   */
    vector<double> g1L_;
    
    /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,1/3}\f$
   */
    vector<double> g1R_;
    
    /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,-2/3}\f$
   */
    vector<double> g1pR_;
    
  /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,4/3}\f$
   */
    vector<double> g1ppR_;
    
  /**
   * Coupling to \f$\Phi_{6,3,1/3}\f$
   */
    vector<double> g3L_;
  //@}

  };
  
}

#endif /* HERWIG_SextetFFSVertex_H */
