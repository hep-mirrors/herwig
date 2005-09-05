// -*- C++ -*-
#ifndef THEPEG_PScalarVectorFermionsDecayer_H
#define THEPEG_PScalarVectorFermionsDecayer_H
//
// This is the declaration of the PScalarVectorFermionsDecayer class.
//
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
// #include "PScalarVectorFermionsDecayer.fh"
// #include "PScalarVectorFermionsDecayer.xh"

namespace Herwig {
using namespace ThePEG;

/**  \ingroup Decay
 *
 * The <code>PScalarVectorFermionsDecayer</code> class is designed for the decay of a 
 * pseudoscalar meson to a spin-1 particle and a fermion-antifermion pair. In practice
 * these decays are of the form \f$\gamma\ell^+\ell^-\f$ and the propagator of
 * the off-shell boson is taken to be \f$\frac1{m^2_{f\bar{f}}}\f$.
 * There is also the option of including a vector meson dominance
 * form-factor.
 *
 *  In this case the matrix element is
 *  \f[\mathcal{M} = \frac{g}{m^2_{f\bar{f}}}
 *                   \epsilon^{\mu\nu\alpha\beta}p_{V\mu}\epsilon_{V\nu}
 *                   \bar{u}(p_f)\gamma_\alpha v(p_{\bar{f}}) p_{f\bar{f}\beta}
 *  \f]
 *  It includes the option of a vector meson dominance (VMD) type form factor  
 *  \f$\frac{-M^2+i\Gamma M}{(m^2_{f\bar{f}}-M^2+i\Gamma M)}\f$.
 *    
 *  The incoming pseudoscalar meson, the outgoing vector, the fermion and antifermion
 *  and the coupling can be specified using the relevant interfaces.
 *
 * @see DecayIntegrator
 * @see PScalarVectorVectorDecayer
 * @see PScalar4FermionsDecayer
 * 
 *  \author Peter Richardson
 *
 */
class PScalarVectorFermionsDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PScalarVectorFermionsDecayer();

  /**
   * Copy-constructor.
   */
  inline PScalarVectorFermionsDecayer(const PScalarVectorFermionsDecayer &);

  /**
   * Destructor.
   */
  virtual ~PScalarVectorFermionsDecayer();
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
  void dataBaseOutput(ofstream &) const;

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
  static ClassDescription<PScalarVectorFermionsDecayer> initPScalarVectorFermionsDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  PScalarVectorFermionsDecayer & operator=(const PScalarVectorFermionsDecayer &);

private:

  /**
   * coupling for a decay
   */
  vector<InvEnergy> _coupling;

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> _incoming;

  /**
   * the PDG codes for the outgoing vector
   */
  vector<int> _outgoingV;

  /**
   * the PDG codes for the outgoing fermion
   */
  vector<int> _outgoingf;

  /**
   * the PDG codes for the outgoing antifermion
   */
  vector<int> _outgoinga;

  /**
   * maximum weight for a decay
   */
  vector<double> _maxweight;

  /**
   * Include the VMD factor
   */
  vector<int> _includeVMD;

  /**
   * PDG code for thte particle to use in the VMD factor.
   */
  vector<int> _VMDid;

  /**
   * Mass to use in the VMD factor.
   */
  vector<Energy> _VMDmass;

  /**
   * Width to use in the VMD factor.
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
 * base class of PScalarVectorFermionsDecayer.
 */
template <>
struct BaseClassTrait<Herwig::PScalarVectorFermionsDecayer,1> {
    /** Typedef of the base class of PScalarVectorFermionsDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::PScalarVectorFermionsDecayer>
  : public ClassTraitsBase<Herwig::PScalarVectorFermionsDecayer> {
   /** Return the class name.*/
   static string className() { return "Herwig++::PScalarVectorFermionsDecayer"; }
   /**
    * Return the name of the shared library to be loaded to get
    * access to this class and every other class it uses
    * (except the base class).
    */
   static string library() { return "libHwSMDecay.so"; }

};

}

#include "PScalarVectorFermionsDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarVectorFermionsDecayer.tcc"
#endif

#endif /* THEPEG_PScalarVectorFermionsDecayer_H */
