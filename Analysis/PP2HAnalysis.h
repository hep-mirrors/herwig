// -*- C++ -*-
#ifndef HERWIG_PP2HAnalysis_H
#define HERWIG_PP2HAnalysis_H
//
// This is the declaration of the PP2HAnalysis class.
//
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/EventRecord/Event.h"
#include "Herwig++/Interfaces/KtJetInterface.h"
#include "KtJet/KtJet.h"
#include "KtJet/KtLorentzVector.h"

namespace Herwig {
  
  using namespace ThePEG;
  /**
   * Here is the documentation of the PP2HAnalysis class.
   *
   * @see \ref PP2HAnalysisInterfaces "The interfaces"
   * defined for PP2HAnalysis.
   */
  class PP2HAnalysis: public AnalysisHandler {
    
  public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
    inline PP2HAnalysis() {}

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
    
    /**
     * Finalize this object. Called in the run phase just after a
     * run has ended. Used eg. to write out statistics.
     */
    virtual void dofinish();
    
  protected:
    
     /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const{ return new_ptr(*this); }
  //@}

  
  private:
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PP2HAnalysis> initPP2HAnalysis;
    
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PP2HAnalysis & operator=(const PP2HAnalysis &);
    
  public:
    
    /**
     * A vector of the final state particles. 
     */
    //    tPVector particles;
    tPVector particles;

    /**
     * The Higgs. 
     */
    tPPtr Higgs;

  private:

  /**
   *  The interface between Herwig++ and KtJet
   */
  Herwig::KtJetInterface _kint;



  };
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {
  
  /** @cond TRAITSPECIALIZATIONS */
  
  /** This template specialization informs ThePEG about the
   *  base classes of PP2HAnalysis. */
  template <>
 struct BaseClassTrait<Herwig::PP2HAnalysis,1> {
    /** Typedef of the first base class of PP2HAnalysis. */
    typedef AnalysisHandler NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the PP2HAnalysis class and the shared object where it is defined. */
  template <>
  struct ClassTraits<Herwig::PP2HAnalysis>
    : public ClassTraitsBase<Herwig::PP2HAnalysis> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::PP2HAnalysis"; }
    /**
     * The name of a file containing the dynamic library where the class
     * PP2HAnalysis is implemented. It may also include several, space-separated,
     * libraries if the class PP2HAnalysis depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwKtJet.so HwPP2HAnalysis.so"; }
  };
  
  /** @endcond */

}

#endif /* HERWIG_PP2HAnalysis_H */


