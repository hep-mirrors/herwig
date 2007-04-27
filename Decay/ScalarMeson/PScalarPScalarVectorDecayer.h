// -*- C++ -*-
#ifndef HERWIG_PScalarPScalarVectorDecayer_H
#define HERWIG_PScalarPScalarVectorDecayer_H
//
// This is the declaration of the PScalarPScalarVectorDecayer class.
//
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "PScalarPScalarVectorDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>PScalarPScalarVectorDecayer</code> class is designed to perform the decay
 * of a pseudoscalar meson to another pseudoscalar meson and a vector meson.
 * There are only a few of these decays. In this case the matrix element has the 
 * form
 *  \f[\mathcal{M} = g\epsilon_2^\mu(p_0+p_1)_\mu,\f]
 *  where
 * - \f$p_0\f$ is the momentum of the incoming pseudoscalar meson.
 * - \f$p_1\f$ is the momentum of the outgoing pseudoscalar meson.
 * - \f$\epsilon_2\f$ is the polarization vector of the vector meson.
 * - \f$g\f$ is the coupling for the decay.
 *
 * @see DecayIntegrator
 *
 * \author Peter Richardson
 * 
 */
class PScalarPScalarVectorDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  PScalarPScalarVectorDecayer();

  /**
   * Copy-constructor.
   */
  inline PScalarPScalarVectorDecayer(const PScalarPScalarVectorDecayer &);

  /**
   * Destructor.
   */
  virtual ~PScalarPScalarVectorDecayer();

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
   *               in the GenericWidthGenerator class, in this case 10.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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
  static ClassDescription<PScalarPScalarVectorDecayer> initPScalarPScalarVectorDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  PScalarPScalarVectorDecayer & operator=(const PScalarPScalarVectorDecayer &);

private:

  /**
   * the PDG code for the incoming particle
   */
  vector<int> _incoming;

  /**
   * the PDG code for the outgoing pseudoscalar
   */
  vector<int> _outgoingP;

  /**
   * the PDG code for the outgoing vector
   */
  vector<int> _outgoingV;

  /**
   * the coupling for the decay
   */
  vector<double> _coupling;

  /**
   * the maximum weight for the decay
   */
  vector<double> _maxweight;

  /**
   *  initial number of modes
   */
  unsigned int _initsize;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of PScalarPScalarVectorDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::PScalarPScalarVectorDecayer,1> {
    /** Typedef of the base class of PScalarPScalarVectorDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::PScalarPScalarVectorDecayer>
  : public ClassTraitsBase<Herwig::PScalarPScalarVectorDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::PScalarPScalarVectorDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwWeakCurrents.so HwSMDecay.so"; }

};

/** @endcond */

}

#include "PScalarPScalarVectorDecayer.icc"

#endif /* HERWIG_PScalarPScalarVectorDecayer_H */
