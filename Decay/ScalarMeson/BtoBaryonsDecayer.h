// -*- C++ -*-
#ifndef HERWIG_BtoBaryonsDecayer_H
#define HERWIG_BtoBaryonsDecayer_H
//
// This is the declaration of the BtoBaryonsDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "BtoBaryonsDecayer.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The BtoBaryonsDecayer class implements the results of PRD65 054028, hep-ph/0210275
 * and hep-ph/0208185 for the decay of B mesons to baryons.
 * 
 * @see DecayIntegrator
 * @see \ref BtoBaryonsDecayerInterfaces "The interfaces"
 * defined for BtoBaryonsDecayer.
 */
class BtoBaryonsDecayer: public DecayIntegrator {

public:

  //@{
  /**
   * The default constructor.
   */
  inline BtoBaryonsDecayer();

public:
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param dm The decay mode
   */
  virtual int modeNumber(bool & cc,const DecayMode & dm) const;
  
  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(bool vertex, const int ichan,const Particle & part,
	     const ParticleVector & decay) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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
  static ClassDescription<BtoBaryonsDecayer> initBtoBaryonsDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BtoBaryonsDecayer & operator=(const BtoBaryonsDecayer &);

private:

  /**
   * The two body decay matrix element.
   * @param vertex Generate the information for spin correlations.
   * @param ichan The phase-space channel.
   * @param inpart The decaying particle.
   * @param decay The decay products.
   * @return The matrix element squared.
   */
  double twoBodyME(bool vertex, const int ichan,const Particle & inpart,
		   const ParticleVector & decay) const;

private:

  /**
   *  The Fermi constant, \f$G_F\f$
   */
  InvEnergy2 _gf;

  /** @name Perturbative coeffficients */
  //@{
  /**
   * The perturbative coefficient \f$c_1^{\rm eff}\f$
   */
  double _c1eff;

  /**
   * The perturbative coefficient \f$c_2^{\rm eff}\f$
   */
  double _c2eff;
  //@}

  /** @name Wave function parameters */
  //@{
  /**
   *  Four-quark overlap bag intergral, \f$X\f$.
   */
  Energy3 _x;

  /**
   *  Four-quark overlap bag intergral, \f$X_1\f$.
   */
  Energy3 _x1;

  /**
   *  Four-quark overlap bag intergral, \f$X_2\f$.
   */
  Energy3 _x2;
  //@}

  /** @name Strong couplings */
  //@{
  /**
   * Strong coupling for \f$\Sigma_b^+\to \bar{B}^0p\f$
   */
  double _gsigmabBbar0p;

  /**
   * Strong coupling for \f$\Lambda^0_b\to B^-p\f$.
   */
  double _glambdabBminusp;

  /**
   * Strong coupling for \f$\Sigma_b^0\to B^-p\f$.
   */
  double _gsigmabBminusp;

  /**
   * Strong coupling for \f$\Lambda^0_b\to \bar{B}^0n\f$.
   */
  double _glambdabBbar0n;

  /**
   * Strong coupling for \f$\Sigma_b^0\to \bar{B}^0n\f$.
   */
  double _gsigmabBbar0n;

  /**
   * Strong coupling for \f$\Sigma_b^+\to B^-\Delta^{++}\f$.
   */
  double _gsigmabBminusDelta;
  //@}

  /**
   *  The \f$A\f$ coefficient for the decays
   */ 
  vector<Complex> _a;

  /**
   *  The \f$B\f$ coefficient for the decays.
   */
  vector<Complex> _b;

  /**
   *  The PDG code for the incoming baryon
   */
  vector<int> _incoming;

  /**
   *  The PDG code for the outgoing baryon
   */
  vector<int> _outgoingB;

  /**
   *  The PDG code for the outgoing antibaryon
   */
  vector<int> _outgoingA;

  /**
   *  The PDG code for the outgoing meson
   */
  vector<int> _outgoingM;

  /**
   * location of the weights
   */
  vector<int> _wgtloc;

  /**
   * the maximum weights
   */
  vector<double> _wgtmax;

  /**
   * the weight for different channels
   */
  vector<double> _weights;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BtoBaryonsDecayer. */
template <>
struct BaseClassTrait<Herwig::BtoBaryonsDecayer,1> {
  /** Typedef of the first base class of BtoBaryonsDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BtoBaryonsDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BtoBaryonsDecayer>
  : public ClassTraitsBase<Herwig::BtoBaryonsDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BtoBaryonsDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * BtoBaryonsDecayer is implemented. It may also include several, space-separated,
   * libraries if the class BtoBaryonsDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwWeakCurrents.so HwSMDecay.so"; }
};

/** @endcond */

}

#include "BtoBaryonsDecayer.icc"

#endif /* HERWIG_BtoBaryonsDecayer_H */
