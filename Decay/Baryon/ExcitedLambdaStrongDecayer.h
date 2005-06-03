// -*- C++ -*-
#ifndef HERWIG_ExcitedLambdaStrongDecayer_H
#define HERWIG_ExcitedLambdaStrongDecayer_H
//
// This is the declaration of the ExcitedLambdaStrongDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ExcitedLambdaStrongDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * The ExcitedLambdaStrongDecayer class uses results from heavy quark chiral
 * pertubration theory for the decays of the excited \f$\Lambda^{(*)}_{c,b}\f$ 
 * baryons to two pions and the \f$\Lambda_{b,c}\f$ or decay of \f$\Xi^{(*)}_{c,b1}\f$
 * to two pions and the \f$Xi_{c,b}\f$. For the  \f$\Lambda^{(*)}_{c,b}\f$ both \f$s\f$
 * and \f$d\f$-wave couplings are included while for \f$\Xi^{(*)}_{c,b1}\f$ only the
 * \f$s\f$-wave is included. The matrix elements are based on PRD56, 5483 but we get
 * slightly different results as our use of helicity amplitudes gives some additional
 * contributions which are subleading in the heavy baryon limit.
 *
 * @see DecayIntegrator
 */
class ExcitedLambdaStrongDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ExcitedLambdaStrongDecayer();

  /**
   * The copy constructor.
   */
  inline ExcitedLambdaStrongDecayer(const ExcitedLambdaStrongDecayer &);

  /**
   * The destructor.
   */
  virtual ~ExcitedLambdaStrongDecayer();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * This function is purely virtual and must be implemented in classes inheriting
   * from DecayIntegrator.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		      const ParticleVector & decay) const;

  /**
   * Method to return an object to calculate the 3 (or higher body) partial width
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;
  
  /**
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(int imode,Energy2 q2, Energy2 s3,Energy2 s2,
					Energy2 s1,Energy m1,Energy m2,Energy m3);  
  //@}

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
   * Initialize this object. Called in the run phase just before
   * a run begins.
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
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ExcitedLambdaStrongDecayer> initExcitedLambdaStrongDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ExcitedLambdaStrongDecayer & operator=(const ExcitedLambdaStrongDecayer &);

private:

  /**
   *  The pion decay constant \f$f_\pi\f$.
   */
  Energy _fpi;
  
  /**
   *  The \f$g_2\f$ coupling.
   */
  double _g2;

  /**
   * The\f$h_2\f$ coupling.
   */
  double _h2;

  /**
   * The \f$h_8\f$ coupling
   */
  InvEnergy _h8;
  
  /**
   * The PDG code for the incoming baryon
   */
  vector<int> _incoming;

  /**
   *  The PDG code for the outgoing baryon
   */
  vector<int> _outgoing;

  /**
   *  Whether the pions are neutral or charged
   */
  vector<int> _charged;

  /**
   * location of the weights
   */
  vector<int> _wgtloc;

  /**
   * the maximum weights and the weights
   */
  vector<double> _wgtmax,_weights;
};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of ExcitedLambdaStrongDecayer. */
template <>
struct BaseClassTrait<Herwig::ExcitedLambdaStrongDecayer,1> {
  /** Typedef of the first base class of ExcitedLambdaStrongDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ExcitedLambdaStrongDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ExcitedLambdaStrongDecayer>
  : public ClassTraitsBase<Herwig::ExcitedLambdaStrongDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::ExcitedLambdaStrongDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the ExcitedLambdaStrongDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "libHwBaryonDecay.so"; }
};

}

#include "ExcitedLambdaStrongDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ExcitedLambdaStrongDecayer.tcc"
#endif

#endif /* HERWIG_ExcitedLambdaStrongDecayer_H */
