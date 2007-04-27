// -*- C++ -*-
#ifndef HERWIG_ExcitedLambdaStrongDecayer_H
#define HERWIG_ExcitedLambdaStrongDecayer_H
//
// This is the declaration of the ExcitedLambdaStrongDecayer class.
//

#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Decay/DecayIntegrator.h"
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
 * @see \ref ExcitedLambdaStrongDecayerInterfaces "The interfaces"
 * defined for ExcitedLambdaStrongDecayer.
 */
class ExcitedLambdaStrongDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  ExcitedLambdaStrongDecayer();

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
   * Method to return an object to calculate the 3 body partial width.
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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();
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

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

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
  /**
   * The name of a file containing the dynamic library where the class
   * ExcitedLambdaStrongDecayer is implemented. It may also include several, space-separated,
   * libraries if the class ExcitedLambdaStrongDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwBaryonDecay.so"; }
};

/** @endcond */

}

#include "ExcitedLambdaStrongDecayer.icc"

#endif /* HERWIG_ExcitedLambdaStrongDecayer_H */
