// -*- C++ -*-
#ifndef HERWIG_EtaPiPiPiDecayer_H
#define HERWIG_EtaPiPiPiDecayer_H
// This is the declaration of the EtaPiPiPiDecayer class.

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "EtaPiPiPiDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>EtaPiPiPiDecayer</code> class is designed for the simulation of
 * the decay of the \f$\eta\f$ or \f$\eta'\f$ to either \f$\pi^+\pi^-\pi^0\f$ 
 * or \f$\pi^0\pi^0\pi^0\f$ and the decay of the \f$\eta'\f$ to 
 * \f$\pi^+\pi^-\eta\f$ or \f$\pi^0\pi^0\eta\f$
 *
 *  The matrix element takes the form 
 * \f[ |\mathcal{M}|^2 = N\left[1+ay+by^2+cx^2\right],\f]
 *  where 
 * \f[x = \frac{\sqrt{3}(u-t)}{2M_0(M_0-m_1-m_2-m_3)},\f]
 * \f[y = \frac{(m_1+m_2+m_3)((M_0-m_3)^2-s)}{2M_0(m_1+m_2)(M_0-m_1-m_2-m_3)}-1,\f]
 *  where 
 * - \f$m_{1,2,3}\f$ are the masses of the outgoing mesons
 * - \f$u = (p_0-p_1)^2\f$,
 * - \f$t = (p_0-p_2)^2\f$,
 * - \f$s = (p_0-p_3)^2\f$.
 *
 *  This form is taken from hep-ph/0301058 as are the experimental results for 
 * the constants which are used where available and the theory results which are
 * used when there is no experimental data. 
 *
 * @see DecayIntegrator
 * 
 */
class EtaPiPiPiDecayer: public DecayIntegrator {

public:


  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  EtaPiPiPiDecayer();

  /**
   * Copy-constructor.
   */
  inline EtaPiPiPiDecayer(const EtaPiPiPiDecayer &);

  /**
   * Destructor.
   */
  virtual ~EtaPiPiPiDecayer();
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
  static ClassDescription<EtaPiPiPiDecayer> initEtaPiPiPiDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  EtaPiPiPiDecayer & operator=(const EtaPiPiPiDecayer &);

private:

  /**
   * the id of the incoming particle
   */
  vector<int> _incoming;

  /**
   * the id of the last neutral meson
   */
  vector<int> _outgoing;

  /**
   * whether the pions are charged or neutral
   */
  vector<bool> _charged;

  /**
   * the prefactor for the decay
   */
  vector<double> _prefactor;

  /**
   * the constants for the matrix elements
   */
  vector<double> _a,_b,_c;

  /**
   * maximum weights
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
 * base class of EtaPiPiPiDecayer.
 */
template <>
struct BaseClassTrait<Herwig::EtaPiPiPiDecayer,1> {
    /** Typedef of the base class of  EtaPiPiPiDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::EtaPiPiPiDecayer>
  : public ClassTraitsBase<Herwig::EtaPiPiPiDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::EtaPiPiPiDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSMDecay.so"; }

};

}

#include "EtaPiPiPiDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EtaPiPiPiDecayer.tcc"
#endif

#endif /* HERWIG_EtaPiPiPiDecayer_H */
