// -*- C++ -*-
#ifndef HERWIG_Hw64Decayer_H
#define HERWIG_Hw64Decayer_H
//
// This is the declaration of the Hw64Decayer class.
//
#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/PDT/Decayer.h>
#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/Strategy.fh>
#include <fstream>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Decay
 *
 * <code>Hw64Decayer</code> is a class that defines all the general routines 
 * used in HERWIG++ to imitate the HERWIG 6.4 decays. The goal is to have an exact
 * copy of HERWIG 6.4 decay routines. This will allow for easy 'callibration'
 * of the new C++ code with the old Fortran code.
 *
 *  This class handles the non-partonic decays. In general it is used for
 *  exclusive meson and baryon decays. Three different matrix elements are supported
 *
 *  - MECode=0 flat-phase space.
 *  - MECode=100 free V-A matrix element
 *  - MECode=101 bound V-A matrix element
 *
 * @see HeavyDecayer
 * @see QuarkoniumDecayer
 * @see Decayer
 * 
 */
class Hw64Decayer: public Decayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline Hw64Decayer();

  /**
   * Copy constructor
   */
  inline Hw64Decayer(const Hw64Decayer &);

  /**
   * Destructor
   */
  virtual ~Hw64Decayer();
  //@}

public:

  /**
   * return true if this decayer can perfom the decay specified by the
   * given decay mode.
   */
  virtual bool accept(const DecayMode &) const;

  /**
   * for a given decay mode and a given particle instance, perform the
   * decay and return the decay products.
   */
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

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

private:

  /**
   * Perform a one body decay, used for \f$K^0,\bar{K}^0\to K_{L,S}\f$.
   * Two body decay is handled in static class Kinematics
   */
  static void oneBodyDecay(Lorentz5Momentum, Lorentz5Momentum &);

  /**
   * Weighting of phase space for V-A matrix elements
   */
  static double VAWt(double*);

  /**
   * Take an array of momenta and set the momentum member of the particles.
   * @param moms The input momenta to be assigned to the particles.
   * @param particles The particles whose momenta is to be set.
   * @param out The particles outputted with their momenta set.
   */
  void setParticleMomentum(ParticleVector & out, cPDVector particles, 
			   vector<Lorentz5Momentum> moms) const;

  /**
   *  Describe a concrete class with persistant data.
   */
  static ClassDescription<Hw64Decayer> initHw64Decayer;

  /**
   *  Private and non-existent assignment operator.
   */
  const Hw64Decayer & operator=(const Hw64Decayer &);

private:

  /**
   *  The code for the matrix element being used.
   */
  int MECode;

  /**
   *  Maximum number of attempts to generate the off-shell masses
   */
  unsigned int _masstry;
};

}

namespace ThePEG {

/**
 * This template specialization informs ThePEG about the base class of
 * Hw64Decayer.
 */
template <>
struct BaseClassTrait<Herwig::Hw64Decayer,1> {
  /** Typedef of the base class of Hw64Decayer. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Hw64Decayer class.
 */
template <>
struct ClassTraits<Herwig::Hw64Decayer>: public ClassTraitsBase<Herwig::Hw64Decayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::Hw64Decayer"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwDecay.so"; }
};

}

#include "Hw64Decayer.icc"

#endif /* HERWIG_Hw64Decayer_H */
