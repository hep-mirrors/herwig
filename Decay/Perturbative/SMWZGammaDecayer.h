// -*- C++ -*-
#ifndef HERWIG_SMWZGammaDecayer_H
#define HERWIG_SMWZGammaDecayer_H
//
// This is the declaration of the SMWZGammaDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "SMWZGammaDecayer.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The SMWZGammaDecayer class generates the decayers of the W and Z including
 * an extra hard photon. It is solely intended for comparision with SOPTHY
 * approach.
 *
 * @see \ref SMWZGammaDecayerInterfaces "The interfaces"
 * defined for SMWZGammaDecayer.
 */
class SMWZGammaDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  inline SMWZGammaDecayer();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const PDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMWZGammaDecayer> initSMWZGammaDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMWZGammaDecayer & operator=(const SMWZGammaDecayer &);

private:

  /**
   * Pointer to the vertex objects
   */
  //@{
  /**
   * Pointer to the Z vertex
   */
  Helicity::FFVVertexPtr _zvertex;

  /**
   * Pointer to the W vertex
   */
  Helicity::FFVVertexPtr _wvertex;

  /**
   * Pointer to the \f$\gamma\f$ vertex
   */
  Helicity::FFVVertexPtr _pvertex;

  /**
   *  Pointer to the triple electroweak gauge boson vertex
   */
  Helicity::VVVVertexPtr _wwwvertex;
  //@}

  /**
   * maximum weights for the different integrations
   */
  //@{
  /**
   *  Weights for the Z to quarks decays.
   */
  vector<double> _Zquarkwgt;

  /**
   *  Weights for the Z to leptons decays.
   */
  vector<double> _Zleptonwgt;

  /**
   *  Weights for the W to quarks decays.
   */
  vector<double> _Wquarkwgt;

  /**
   *  Weights for the W to leptons decays.
   */
  vector<double> _Wleptonwgt;
  //@}

  /**
   * Weights for phase-space channels for the different integrations
   */
  //@{
  /**
   *  Weights for the Z to quarks decays.
   */
  vector<double> _Zquarkchannels;

  /**
   *  Weights for the Z to leptons decays.
   */
  vector<double> _Zleptonchannels;

  /**
   *  Weights for the W to quarks decays.
   */
  vector<double> _Wquarkchannels;

  /**
   *  Weights for the W to leptons decays.
   */
  vector<double> _Wleptonchannels;
  //@}

  /**
   *  Minimum energy to avoid soft singularity
   */
  Energy _emin;

  /**
   *  Option for the calculation of the matrix element
   */
  unsigned int _iopt;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMWZGammaDecayer. */
template <>
struct BaseClassTrait<Herwig::SMWZGammaDecayer,1> {
  /** Typedef of the first base class of SMWZGammaDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMWZGammaDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMWZGammaDecayer>
  : public ClassTraitsBase<Herwig::SMWZGammaDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMWZGammaDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SMWZGammaDecayer is implemented. It may also include several, space-separated,
   * libraries if the class SMWZGammaDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPerturbativeDecay.so"; }
};

/** @endcond */

}

#include "SMWZGammaDecayer.icc"

#endif /* HERWIG_SMWZGammaDecayer_H */
