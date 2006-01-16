#ifndef HERWIG_FORCED_SPLITTING_H_
#define HERWIG_FORCED_SPLITTING_H_

#include <ThePEG/Handlers/MultipleInteractionHandler.h>
#include <ThePEG/Handlers/EventHandler.h>
#include "Herwig++/Shower/BackwardEvolver.h" 

namespace Herwig {

using namespace ThePEG;

/** \ingroup UnderlyingEvent
 *
 *  This is the definition of a simple forced splitting algorithm.
 *  This takes the Remnant object produced from the PDF and backward
 *  evolution (hadron - parton) and produce partons with the remaining 
 *  flavours and with the correct colour connections.
 */

class ForcedSplitting : public MultipleInteractionHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
   ForcedSplitting();

  /**
   * The copy constructor.
   */
  ForcedSplitting(const ForcedSplitting &);

  /**
   * The destructor.
   */
   virtual ~ForcedSplitting();
   //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

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

public:

  /**
   * This is the routine that starts the algorithm.
   */
  virtual void handle(EventHandler &ch, const tPVector &tagged,
		      const Hint &hint) 
    throw(Veto,Stop,Exception);

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
   * Initialize this object after the setup phase before saving an
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
   double kinCutoff;

   static ClassDescription<ForcedSplitting> initForcedSplitting;

   /**
    * This is never defined and since it can never be called it isn't 
    * needed. The prototype is defined so the compiler doesn't use the 
    * default = operator.
    */
  ForcedSplitting& operator=(const ForcedSplitting &);

  // Using the remnant given as the first argument and the last parton in
  // the shower, split the remnant into the correct flavours and add them
  // to the step.
  void split(const tPPtr ,const tPPtr, const tStepPtr,const double x );
  
  // This takes the particle and find a splitting for np -> p + child and 
  // creates the correct kinematics and connects for such a split. This
  // Splitting has an upper bound on qtilde given by the energy argument
  // The momentum pf is the momentum of the last reconstructed parton in the
  // backward chain and the momentum p is the momentum which is extracted from
  // the remnant in each step.
  PPtr forceSplit(const tPPtr rem, long child, Energy &oldQ, double &oldx, 
		  Lorentz5Momentum &pf, Lorentz5Momentum &p,const tStepPtr);


  // This computes the momentum of the emitted parton. par is the momentum of
  // the beam particle, lastQ and lastx are the values of qtilde and x after
  // the last emission, emittedm2 is the mass of the parton being emitted
  // and pf is the momentum of the previously reconstructed momentum.
  Lorentz5Momentum emit(const Lorentz5Momentum &par, Energy &lastQ, 
			double &lastx, double emittedm2, Lorentz5Momentum &pf);

  // This creates a parton from the remaining flavours of the hadron. The
  // last parton used was a valance parton, so only 2 (or 1, if meson) flavours
  // remain to be used.
  PPtr finalSplit(const tPPtr rem, int maxIdx, long q[3], int, Lorentz5Momentum,
		  const tStepPtr );


};

}

namespace ThePEG {

/** This template specialization informs ThePEG about the
 *  base classes of ForcedSplitting. */
template<>
struct BaseClassTrait<Herwig::ForcedSplitting,1> { 
  /** Typedef of the first base class of ForcedSplitting. */
  typedef MultipleInteractionHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ForcedSplitting class and the shared object where it is defined. */
template<>
struct ClassTraits<Herwig::ForcedSplitting> :
  public ClassTraitsBase<Herwig::ForcedSplitting> {
  /** Return a platform-independent class name */
    static string className() { return "/Herwig++/ForcedSplitting"; }
  /** Return the name of the shared library be loaded to get
   *  access to the ForcedSplitting class and every other class it uses
   *  (except the base class). */
    static string library() { return "HwForcedSplitting.so"; }
};

}

#include "ForcedSplitting.icc"

#endif
