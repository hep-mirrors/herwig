// -*- C++ -*-
#ifndef HERWIG_VectorMesonPScalarFermionsDecayer_H
#define HERWIG_VectorMesonPScalarFermionsDecayer_H
//
// This is the declaration of the VectorMesonPScalarFermionsDecayer class.
//
#include "VectorMesonDecayerBase.h"
// #include "VectorMesonPScalarFermionsDecayer.fh"
// #include "VectorMesonPScalarFermionsDecayer.xh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>VectorMesonPScalarFermionsDecayer</code> class is designed to perform the
 * decay of a vector meson to a pesudo scalar and a fermion-antifermion pair according
 * to a current which is the \f$V\to VP\f$ vertex combined with the branching of the
 * vector into a fermion-antifermion pair.
 *
 *  The current is
 *  \f[J^\mu = \frac{g}{(p_f+p_{\bar f})^2}\epsilon^{\mu\nu\alpha\beta} p_{0\nu}
 *             (p_f+p_{\bar f})_\alpha
 *             \bar{u}(p_f)\gamma_\beta v(p_{\bar f})
 *  \f]
 *
 *  It includes the option of a vector meson dominance (VMD) type form factor  
 *  \f$\frac{-M^2+i\Gamma M}{(m^2_{ff}-M^2+i\Gamma M)}\f$.
 *
 *  The incoming and outgoing meson together with the types of fermions can be
 *  specified using the interfaces.
 *
 * @see VectorMesonDecayerBase
 * @see VectorMesonVectorPScalarDecayer
 * 
 *  \author Peter Richardson
 *
 */
class VectorMesonPScalarFermionsDecayer: public VectorMesonDecayerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  VectorMesonPScalarFermionsDecayer();

  /**
   * Copy-constructor.
   */
  inline VectorMesonPScalarFermionsDecayer(const VectorMesonPScalarFermionsDecayer &);

  /**
   * Destructor.
   */
  virtual ~VectorMesonPScalarFermionsDecayer();
  //@}

public:

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This version checks the particles against the 
   * list of allowed incoming and outgoing mesons and fermions.
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
   * The hadronic current. This returns the current 
   *  described above.
   * @param vertex Construct the information for spin correlations.
   * @param ichan The phase-space channel to calculate the current for.
   * @param inpart The decaying particle
   * @param outpart The decay products
   * @return The hadronic currents for the decay.
   */
  virtual vector<LorentzPolarizationVector> 
  decayCurrent(const bool vertex, const int ichan,const Particle & inpart, 
	       const ParticleVector & outpart) const;

  /**
   * Method to return an object to calculate the 3 body partial width.
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;

  /**
   * The differential three body decay rate with one integral performed.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s  The invariant mass which still needs to be integrate over.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The differential rate \f$\frac{d\Gamma}{ds}\f$
   */
  virtual double threeBodydGammads(int imode,Energy q2, Energy2 s,Energy m1,Energy m2,
				   Energy m3);

  /**
   * Output the setup information for the particle database
   */
  void dataBaseOutput(ofstream &);

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
  virtual void doinit() throw(InitException);

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
  static ClassDescription<VectorMesonPScalarFermionsDecayer> initVectorMesonPScalarFermionsDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  VectorMesonPScalarFermionsDecayer & operator=(const VectorMesonPScalarFermionsDecayer &);

private:
  
  /**
   * coupling for a decay
   */
  vector<InvEnergy> _coupling;

  /**
   * PDG codes for the incoming particle
   */
  vector<int> _incoming;

  /**
   * PDG codes for the outgoing meson.
   */
  vector<int> _outgoingP;

  /**
   * PDG codes for the outgoing fermion.
   */
  vector<int> _outgoingf;

  /**
   * PDG codes for the outgoing antifermion.
   */
  vector<int> _outgoinga;

  /**
   * Maximum weight for a decay
   */
  vector<double> _maxweight;

  /**
   * Relative weights for the two channels
   */
  vector<double> _weight;

  /**
   * Include the VMD form factor.
   */
  vector<int> _includeVMD;

  /**
   * PDG code for the particle mass and width to use for the VMD form factor.
   */
  vector<int> _VMDid;

  /**
   * Mass for the VMD form factor.
   */
  vector<Energy> _VMDmass;

  /**
   * Width for the VMD form factor.
   */
  vector<Energy> _VMDwidth;

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
 * base class of Herwig::VectorMesonPScalarFermionsDecayer.
 */
template <>
struct BaseClassTrait<Herwig::VectorMesonPScalarFermionsDecayer,1> {
    /** Typedef of the base class of VectorMesonPScalarFermionsDecayer. */
  typedef Herwig::VectorMesonDecayerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::VectorMesonPScalarFermionsDecayer>
  : public ClassTraitsBase<Herwig::VectorMesonPScalarFermionsDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::VectorMesonPScalarFermionsDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVNDecay.so"; }

};

}

#include "VectorMesonPScalarFermionsDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonPScalarFermionsDecayer.tcc"
#endif

#endif /* HERWIG_VectorMesonPScalarFermionsDecayer_H */
