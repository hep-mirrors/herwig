// -*- C++ -*-
#ifndef HERWIG_SMWZDecayer_H
#define HERWIG_SMWZDecayer_H
//
// This is the declaration of the SMWZDecayer class.
//
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "SMWZDecayer.fh"
#include "Herwig++/Models/StandardModel/StandardModel.h"

namespace Herwig {
using namespace ThePEG;
using namespace Herwig::Helicity;
typedef Ptr<Herwig::Helicity::FFVVertex>::pointer FFVPtr;

/** \ingroup Decay
 *
 *  The <code>SMWZDecayer</code> is designed to perform the decay of the 
 *  W and Z bosons to the Standard Model fermions. In principle it can also
 *  be used for these decays in any model.
 *
 * @see DecayIntegrator
 * 
 */
class SMWZDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline SMWZDecayer();

  /**
   * Copy-constructor.
   */
  inline SMWZDecayer(const SMWZDecayer &);

  /**
   * Destructor.
   */
  virtual ~SMWZDecayer();
  //@}

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const;

public:

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. 
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   *
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. 
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;

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
  static ClassDescription<SMWZDecayer> initSMWZDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SMWZDecayer & operator=(const SMWZDecayer &);

 private:

  /**
   * Pointer to the Z vertex
   */
  mutable FFVPtr _Zvertex;

  /**
   * Pointer to the W vertex
   */
  mutable FFVPtr _Wvertex;

  /**
   * maximum weights for the different integrations
   */
  //@{
  /**
   *  Weights for the Z to quarks decays.
   */
  vector<double> _Zquarkwgt;

  /**
   *  Weights for the Z to leptons decays.
   */
  vector<double> _Zleptonwgt;

  /**
   *  Weights for the W to quarks decays.
   */
  vector<double> _Wquarkwgt;

  /**
   *  Weights for the W to leptons decays.
   */
  vector<double> _Wleptonwgt;
  //@}
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of SMWZDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SMWZDecayer,1> {
    /** Typedef of the base class of SMWZDecayer. */
   typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SMWZDecayer>
  : public ClassTraitsBase<Herwig::SMWZDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig++::SMWZDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwPerturbativeDecayer.so"; }

};

}

#include "SMWZDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SMWZDecayer.tcc"
#endif

#endif /* HERWIG_SMWZDecayer_H */
