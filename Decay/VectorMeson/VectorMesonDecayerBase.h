// -*- C++ -*-
#ifndef HERWIG_VectorMesonDecayerBase_H
#define HERWIG_VectorMesonDecayerBase_H
//
// This is the declaration of the VectorMesonDecayerBase class.
//
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

/** \ingroup Decay
 *
 *  This class is designed to be the base class for the HERWIG++ decays of
 *  vector mesons. It handles the generation of the phase space and the calculation
 *  of the matrix element. All the implementations should inherit from this class
 *  and implement the decayCurrent method to return the current for a given
 *  phase space point. This current is then contracted with the polarization
 *  vectors for the decaying meson.
 *
 *  It uses the DecayIntegrator class for the phase space and implements the 
 *  me2 member of the DecayIntegrator to give the matrix element as
 *
 *  \f[ \sum_{\rm spins}|\mathcal{M}|^2 = \sum_{\rm spins} \epsilon^\mu J_\mu \f]
 *
 *  where \f$\epsilon^\mu\f$ is the polarization vector of the decaying vector
 *  meson and \f$J_\mu\f$ is the hadronic current supplied by the decayCurrent member.
 *
 * @see DecayIntegrator
 * 
 *  \author Peter Richardson
 */
class VectorMesonDecayerBase: public DecayIntegrator {

public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline VectorMesonDecayerBase();

  /**
   * Copy-constructor.
   */
  inline VectorMesonDecayerBase(const VectorMesonDecayerBase &);

  /**
   * Destructor.
   */
  virtual ~VectorMesonDecayerBase();
  //@}
  
public:
  
  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. As this is the base class it returns false and
   * should be overridden in class implementing the decays.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;
  
public:
  
  /**
   * The hadronic current. This must be implemented in all classes inheriting from
   * this one and is therefore purely virtual. The currents should be specified so 
   * that the spin of the first outgoing particle is interated over first, then 
   * then second particle and so on.
   * @param vertex Construct the information for spin correlations.
   * @param ichan The phase-space channel to calculate the current for.
   * @param inpart The decaying particle
   * @param outpart The decay products
   * @return The hadronic currents for the decay.
   */
  virtual vector<LorentzPolarizationVector> 
  decayCurrent(const bool vertex, const int ichan,const Particle & inpart, 
	       const ParticleVector & outpart) const =0;
  
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
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<VectorMesonDecayerBase> initVectorMesonDecayerBase;
  
  /**
   * Private and non-existent assignment operator.
   */
  VectorMesonDecayerBase & operator=(const VectorMesonDecayerBase &);
  
};
  
}


namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of VectorMesonDecayerBase.
   */
  template <>
  struct BaseClassTrait<Herwig::VectorMesonDecayerBase,1> {
    /** Typedef of the base class of VectorMesonDecayerBase. */
    typedef Herwig::DecayIntegrator NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::VectorMesonDecayerBase>
    : public ClassTraitsBase<Herwig::VectorMesonDecayerBase> {
    /** Return the class name. */
    static string className() { return "Herwig++::VectorMesonDecayerBase"; }
    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwVMDecay.so"; }

  };
  
}

#include "VectorMesonDecayerBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonDecayerBase.tcc"
#endif

#endif /* HERWIG_VectorMesonDecayerBase_H */
