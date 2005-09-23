// -*- C++ -*-
#ifndef HERWIG_TensorMesonVectorVectorDecayer_H
#define HERWIG_TensorMesonVectorVectorDecayer_H
//
// This is the declaration of the TensorMesonVectorVectorDecayer class.
//
#include "TensorMesonDecayerBase.h"
// #include "TensorMesonVectorVectorDecayer.fh"
// #include "TensorMesonVectorVectorDecayer.xh"

namespace Herwig {
using namespace ThePEG; 

/** \ingroup Decay
 *
 *  The <code>TensorMesonVectorVectorDecayer</code> class is designed to simulate
 *  the decay of a tensor meson to two spin-1 particles. In practice, at least
 *  for the light tensor mesons, this is only the decay of a tensor meson to two
 *  photons. In principle for bottom and charm tensors this may be the decay to
 *  two vector mesons. 
 *
 *  The form of the matrix element is based on the perturbative matrix element
 *  for the decay of a graviton to two vector bosons with the neglect of a mass term
 *
 *  \f[T^{\mu\nu} = g\left[
 *   \left(\epsilon_{1\alpha} p_1^\mu - \epsilon_1^\mu p_{1\alpha}\right)
 *   \left(\epsilon_2^\alpha  p_2^\nu - \epsilon_2^\nu p_2^\alpha\right)
 *  +\left(\epsilon_{1\alpha} p_1^\nu - \epsilon_1^\nu p_{1\alpha}\right)
 *   \left(\epsilon_2^\alpha  p_2^\mu - \epsilon_2^\nu p_2^\alpha\right)
 * -\frac12g^{\mu\nu}
 *    \left(\epsilon_{1\alpha} p_{1\beta}- \epsilon_{1\beta} p_{1\alpha}\right)
 *    \left(\epsilon_2^\alpha  p_2^\beta - \epsilon_2^\beta p_2^\alpha\right)\right] \f]
 *
 *  in such a way that it vanishes if the polarizations of the outgoing vectors are
 *  replaced with their momenta.
 *
 *  The incoming tensor mesons together with their decay products and the coupling 
 *  \f$g\f$ can be specified using the interfaces for the class. The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 *  The incoming and outgoing particles, couplings and maximum weights for
 *  many of the common \f$T\to VV\f$ decays are specified in the default
 *  constructor.
 *
 * @see TensorMesonDecayerBase
 * 
 */
class TensorMesonVectorVectorDecayer: public TensorMesonDecayerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  TensorMesonVectorVectorDecayer();

  /**
   * Copy-constructor.
   */
  inline TensorMesonVectorVectorDecayer(const TensorMesonVectorVectorDecayer &);

  /**
   * Destructor.
   */
  virtual ~TensorMesonVectorVectorDecayer();
  //@}

public:

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This version checks the particles against the 
   * list of allowed incoming  and outgoing mesons.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. This version uses PDG codes to
   * work which mode is being simulated and the generate member of the 
   * DecayIntegrator class for the phase-space  generation.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;

  /**
   * The hadronic tensor. This returns the tensor described above.
   * @param vertex Construct the information for spin correlations.
   * @param ichan The phase-space channel to calculate the current for.
   * @param inpart The decaying particle
   * @param outpart The decay products
   * @return The hadronic currents for the decay.
   */
  virtual vector<LorentzTensor> 
  decayTensor(const bool vertex, const int ichan,const Particle & inpart, 
	       const ParticleVector & outpart) const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 9.
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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<TensorMesonVectorVectorDecayer> initTensorMesonVectorVectorDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  TensorMesonVectorVectorDecayer & operator=(const TensorMesonVectorVectorDecayer &);

private:

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> _incoming;

  /**
   * the PDG codes for the first outgoing particle
   */
  vector<int> _outgoing1;

  /**
   * the PDG codes for the second outgoing particle
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
 * base class of Herwig::TensorMesonVectorVectorDecayer.
 */
template <>
struct BaseClassTrait<Herwig::TensorMesonVectorVectorDecayer,1> {
  /** Typedef of the base class of TensorMesonVectorVectorDecayer. */
  typedef Herwig::TensorMesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::TensorMesonVectorVectorDecayer>
  : public ClassTraitsBase<Herwig::TensorMesonVectorVectorDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::TensorMesonVectorVectorDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwTMDecay.so"; }
  
 };
  
}

#include "TensorMesonVectorVectorDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonVectorVectorDecayer.tcc"
#endif

#endif /* HERWIG_TensorMesonVectorVectorDecayer_H */
