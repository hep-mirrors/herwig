// -*- C++ -*-
#ifndef HERWIG_TensorMeson2PScalarDecayer_H
#define HERWIG_TensorMeson2PScalarDecayer_H
//
// This is the declaration of the TensorMeson2PScalarDecayer class.
//
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
// #include "TensorMeson2PScalarDecayer.fh"
// #include "TensorMeson2PScalarDecayer.xh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>TensorMeson2PScalarDecayer</code> class is designed for the decay
 *  of a tensor meson to two pseudoscalars via matrix element which takes the form
 *  \f[ \mathcal{M} =  g \epsilon^{\mu\nu} p_{1,\mu}p_{2,\nu} \f]
 *  where \f$\epsilon^{\mu\nu}\f$ is the polarization tensor of the decaying tensor
 *  meson, $p_{1,2}$ are the momenta of the decay products and \f$g\f$ is the coupling.
 *
 *  It can also be used for the decay of a tensor to two scalar mesons although
 *  this rarely happens in practice.
 *
 *  The incoming tensor mesons together with their decay products and the coupling 
 *  \f$g\f$ can be specified using the interfaces for the class. The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 *  The incoming and outgoing particles, couplings and maximum weights for
 *  many of the common \f$T\to PP\f$ decays are specified in the default
 *  constructor.
 *
 * @see DecayIntegrator
 *
 * @see \ref TensorMeson2PScalarDecayerInterfaces "The interfaces"
 * defined for TensorMeson2PScalarDecayer.
 * 
 */
class TensorMeson2PScalarDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  TensorMeson2PScalarDecayer();

  /**
   * Copy-constructor.
   */
  inline TensorMeson2PScalarDecayer(const TensorMeson2PScalarDecayer &);

  /**
   * Destructor.
   */
  virtual ~TensorMeson2PScalarDecayer();
  //@}

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
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 7.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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
   * Standard Init function used to initialize the interfaces.
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
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<TensorMeson2PScalarDecayer> initTensorMeson2PScalarDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  TensorMeson2PScalarDecayer & operator=(const TensorMeson2PScalarDecayer &);

private:

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> _incoming;

  /**
   * the PDG codes for first outgoing particle
   */
  vector<int> _outgoing1;

  /**
   * the PDG codes for second outgoing particle
   */
  vector<int> _outgoing2;

  /**
   * the coupling for the decay
   */
  vector<InvEnergy> _coupling;

  /**
   * the maximum weight for the decay
   */
  vector<double> _maxweight;

  /**
   *  Initial size of the vectors
   */
  unsigned int _initsize;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of TensorMeson2PScalarDecayer.
 */
template <>
struct BaseClassTrait<Herwig::TensorMeson2PScalarDecayer,1> {
    /** Typedef of the base class of TensorMeson2PScalarDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::TensorMeson2PScalarDecayer>
  : public ClassTraitsBase<Herwig::TensorMeson2PScalarDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::TensorMeson2PScalarDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwTMDecay.so"; }

};

}

#include "TensorMeson2PScalarDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMeson2PScalarDecayer.tcc"
#endif

#endif /* HERWIG_TensorMeson2PScalarDecayer_H */
