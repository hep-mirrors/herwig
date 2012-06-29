// -*- C++ -*-
#ifndef HERWIG_FtoFVVDecayer_H
#define HERWIG_FtoFVVDecayer_H
//
// This is the declaration of the FtoFVVDecayer class.
//

#include "GeneralThreeBodyDecayer.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The FtoFVVDecayer class provides the general matrix elements for the
 * decay of a fermion to a fermion and two vector bosons.
 *
 * @see \ref FtoFVVDecayerInterfaces "The interfaces"
 * defined for FtoFVVDecayer.
 */
class FtoFVVDecayer: public GeneralThreeBodyDecayer {

public:

  /**
   * Return the matrix element squared for a given mode and phase-space channel
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;
  
  /**
   * Method to return an object to calculate the 3 (or higher body) partial width
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FtoFVVDecayer> initFtoFVVDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FtoFVVDecayer & operator=(const FtoFVVDecayer &);

private:
  
  /**
   * Store the vector of scalar intermediates
   */
  vector<pair<AbstractFFSVertexPtr, AbstractVVSVertexPtr> > _sca;

  /**
   * Store the vector for fermion intermediates
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFVVertexPtr> > _fer;

  /**
   * Store the vector for gauge boson intermediates
   */
  vector<pair<AbstractFFVVertexPtr, AbstractVVVVertexPtr> > _vec;

  /**
   * Store the vector of tensor intermediates
   */
  vector<pair<AbstractFFTVertexPtr, AbstractVVTVertexPtr> > _ten;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Spinor wavefunctions
   */
  mutable vector<SpinorWaveFunction> _fwave;

  /**
   *  Barred spinor wavefunctions
   */
  mutable vector<SpinorBarWaveFunction> _fbwave;

  /**
   *  Vector wavefunctions
   */
  mutable pair<vector<VectorWaveFunction>, vector<VectorWaveFunction> > _vwave;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FtoFVVDecayer. */
template <>
struct BaseClassTrait<Herwig::FtoFVVDecayer,1> {
  /** Typedef of the first base class of FtoFVVDecayer. */
  typedef Herwig::GeneralThreeBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FtoFVVDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FtoFVVDecayer>
  : public ClassTraitsBase<Herwig::FtoFVVDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::FtoFVVDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_FtoFVVDecayer_H */
