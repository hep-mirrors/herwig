// -*- C++ -*-
#ifndef HERWIG_VV_ME_Analysis_H
#define HERWIG_VV_ME_Analysis_H
//
// This is the declaration of the VV_ME_Analysis class.
//
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "VV_ME_Analysis.fh"
#include "Herwig++/Utilities/Histogram.h"


namespace Herwig {
  
  using namespace ThePEG;
  /**
   * Here is the documentation of the VV_ME_Analysis class.
   *
   * @see \ref VV_ME_AnalysisInterfaces "The interfaces"
   * defined for VV_ME_Analysis.
   */
  class VV_ME_Analysis: public AnalysisHandler {
    
  public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
    inline VV_ME_Analysis();
    
    /**
     * The copy constructor.
     */
    inline VV_ME_Analysis(const VV_ME_Analysis &);
    
    /**
     * The destructor.
     */
    virtual ~VV_ME_Analysis();
    //@}
    
  public:
    
    /** @name Virtual functions required by the AnalysisHandler class. */
    //@{
    /**
     * Analyze a given Event. Note that a fully generated event
     * may be presented several times, if it has been manipulated in
     * between. The default version of this function will call transform
     * to make a lorentz transformation of the whole event, then extract
     * all final state particles and call analyze(tPVector) of this
     * analysis object and those of all associated analysis objects. The
     * default version will not, however, do anything on events which
     * have not been fully generated, or have been manipulated in any
     * way.
     * @param event pointer to the Event to be analyzed.
     * @param ieve the event number.
     * @param loop the number of times this event has been presented.
     * If negative the event is now fully generated.
     * @param state a number different from zero if the event has been
     * manipulated in some way since it was last presented.
     */
    virtual void analyze(tEventPtr event, long ieve, int loop, int state);
    
    /**
     * Analyze the given vector of particles. The default version calls
     * analyze(tPPtr) for each of the particles.
     * @param particles the vector of pointers to particles to be analyzed
     */
    virtual void analyze(const tPVector & particles);
    
    
    /**
     * Analyze the given particle.
     * @param particle pointer to the particle to be analyzed.
     */
    virtual void analyze(tPPtr particle);
    //@}
    
  public:
    
    /**
     * A pair of the incoming hadrons.
     */
    PPair inbound;

    /**
     * A vector of the final state leptons. 
     */
    tPVector leptons;

    /**
     * A vector of the final state leptons. 
     */
    tPPtr emitted;

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
    
    /**
     * Finalize this object. Called in the run phase just after a
     * run has ended. Used eg. to write out statistics.
     */
    virtual void dofinish();
    
  protected:
    
    /**
     * Function to return a specific boost to the rest frame of the 
     * leptons to histogram the polar angles in this frame. This is 
     * verbatim c++ translation of mcfm/src/Singletop/boostx.f.
     */
    Lorentz5Momentum boostx(Lorentz5Momentum p_in, Lorentz5Momentum pt,
			    Lorentz5Momentum ptt);
    
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
    
    
    // If needed, insert declarations of virtual function defined in the
    // InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).
    
    
  private:
    
    /**
     * The static object used to initialize the description of this class.
     * Indicates that this is a concrete class with persistent data.
     */
    static ClassDescription<VV_ME_Analysis> initVV_ME_Analysis;
    
    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    VV_ME_Analysis & operator=(const VV_ME_Analysis &);
    
  };
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {
  
  /** @cond TRAITSPECIALIZATIONS */
  
  /** This template specialization informs ThePEG about the
   *  base classes of VV_ME_Analysis. */
  template <>
 struct BaseClassTrait<Herwig::VV_ME_Analysis,1> {
    /** Typedef of the first base class of VV_ME_Analysis. */
    typedef AnalysisHandler NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the VV_ME_Analysis class and the shared object where it is defined. */
  template <>
  struct ClassTraits<Herwig::VV_ME_Analysis>
    : public ClassTraitsBase<Herwig::VV_ME_Analysis> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::VV_ME_Analysis"; }
    /**
     * The name of a file containing the dynamic library where the class
     * VV_ME_Analysis is implemented. It may also include several, space-separated,
     * libraries if the class VV_ME_Analysis depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwVV_ME_Analysis.so"; }
  };
  
  /** @endcond */

}

#include "VV_ME_Analysis.icc"

#endif /* HERWIG_VV_ME_Analysis_H */

