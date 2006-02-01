// -*- C++ -*-
#ifndef HERWIG_SMHiggsGGHiggsPPDecayer_H
#define HERWIG_SMHiggsGGHiggsPPDecayer_H
//
// This is the declaration of the SMHiggsGGHiggsPPDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Helicity/Vertex/StandardModel/SMHGGVertex.h"
#include "Herwig++/Helicity/Vertex/StandardModel/SMHPPVertex.h"
#include "SMHiggsGGHiggsPPDecayer.fh"

namespace Herwig {
  using namespace ThePEG;
  using namespace Herwig::Helicity;

  typedef Ptr<Herwig::Helicity::SMHGGVertex>::pointer HGGPtr;
  typedef Ptr<Herwig::Helicity::SMHPPVertex>::pointer HPPPtr;
  
  /**
   * The <code>SMHiggsGGHiggsPPDecayer</code> class performs the
   * of a Standard Model Higgs boson to either a pair
   * of photons or a pair of gluons.
   *
   * @see DecayIntegrator
   */
  
  class SMHiggsGGHiggsPPDecayer: public DecayIntegrator {
    
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SMHiggsGGHiggsPPDecayer();

  /**
   * The copy constructor.
   */
  inline SMHiggsGGHiggsPPDecayer(const SMHiggsGGHiggsPPDecayer &);

  /**
   * The destructor.
   */
  virtual ~SMHiggsGGHiggsPPDecayer();
 
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
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param dm The decay mode
   */
  virtual int modeNumber(bool & cc,const DecayMode & dm) const {return -1;}

    /**
     * Perform a decay for a given DecayMode and a given Particle instance.
     * @param dm the DecayMode describing the decay.
     * @param p the Particle instance to be decayed.
     * @return a ParticleVector containing the decay products.
     */
    virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
    //@}

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
    inline virtual void doinit() throw(InitException);

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
    static ClassDescription<SMHiggsGGHiggsPPDecayer> 
    initSMHiggsGGHiggsPPDecayer;

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    SMHiggsGGHiggsPPDecayer & operator=(const SMHiggsGGHiggsPPDecayer &);

    /**
     * Pointer to h->gluon,gluon vertex
     */
    HGGPtr _hggvertex;
    
    /**
     * Pointer to h->gamma,gamma vertex
     */
    HPPPtr _hppvertex;
    
    /**
     * Maximum weight for integration
     */
    vector<double> _h0wgt;
  
};
  
}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of SMHiggsGGHiggsPPDecayer. */
template <>
struct BaseClassTrait<Herwig::SMHiggsGGHiggsPPDecayer,1> {
  /** Typedef of the first base class of SMHiggsGGHiggsPPDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHiggsGGHiggsPPDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHiggsGGHiggsPPDecayer>
  : public ClassTraitsBase<Herwig::SMHiggsGGHiggsPPDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SMHiggsGGHiggsPPDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SMHiggsGGHiggsPPDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSMVertex.so HwPerturbativeDecay.so"; }
};

}

#include "SMHiggsGGHiggsPPDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#endif

#endif /* HERWIG_SMHiggsGGHiggsPPDecayer_H */
