// -*- C++ -*-
#ifndef HERWIG_MPIHandler_H
#define HERWIG_MPIHandler_H
//
// This is the declaration of the MPIHandler class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/MultipleInteractionHandler.h"
#include "ThePEG/Handlers/SamplerBase.h"

#include <cassert>

#include "MPIHandler.fh"
#include "stat.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup UnderlyingEvent
   * \class MPIHandler
   * This class is responsible for generating additional 
   * Parton interactions.
   * 
   * \author Manuel Bahr
   *
   * @see \ref MPIHandlerInterfaces "The interfaces"
   * defined for MPIHandler.
   * @see MPISampler
   */

class MPIHandler: public Interfaced, public LastXCombInfo<> {

  /**
   * Class for the integration is a friend to access private members
   */
  friend class Eikonalization;


public:

  /** A vector of <code>SubProcessHandler</code>s. */
  typedef vector<SubHdlPtr> SubHandlerList;

  /** A weighted list of pointers to StandardXComb objects. */
  typedef Selector<StdXCombPtr> XSelector;

  /** A vector of pointers to StandardXComb objects. */
  typedef vector<StdXCombPtr> XVector;

  /** A vector of cross sections. */
  typedef vector<CrossSection> XSVector;

  /** Map of pointers to StandardXComb objects indexed by pointers to
   *  the corresponding MEBase object. */
  typedef map<tMEPtr,XVector> MEXMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MPIHandler();

  /**
   * The copy constructor.
   */
  MPIHandler(const MPIHandler &);

  /**
   * The destructor.
   */
  virtual ~MPIHandler();
  //@}

public:

  /** @name Methods for the MPI generation. */
  //@{
  /**
   * Sample from the pretabulated multiplicity distribution.
   * @return the number of extra events in this collision
   */
  inline unsigned int multiplicity() const;

  /**
   * Select a StandardXComb according to it's weight
   * @return that StandardXComb Object
   */
  inline tStdXCombPtr generate();
  //@}


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
   * Initialize this Multiple Interaction handler and all related objects needed to
   * generate additional events.
   */
  void initialize();

  /**
   * Write out accumulated statistics about intergrated cross sections
   * and stuff.
   */
  void statistics(ostream &, Stat &) const;

  /** @name Functions used for the actual generation */
  //@{
  /**
   * Return the cross section for the chosen phase space point.
   * @param r a vector of random numbers to be used in the generation
   * of a phase space point.
   */
  virtual CrossSection dSigDR(const vector<double> & r);


  /** @name Simple access functions. */
  //@{

  /**
   * Return a reference to the Cuts of this
   * MultipleInteractionHandler. Note that these cuts may be overridden by the
   * SubProcess chosen.
   */
  inline tCutsPtr cuts() const;

  /**
   * The level of statistics. Controlls the amount of statistics
   * written out after each run to the <code>EventGenerator</code>s
   * <code>.out</code> file. Simply the EventHandler method is called here.
   */
  inline int statLevel() const;

  /**
   * Return the ThePEG::EventHandler assigned to this handler.
   * This methods shadows ThePEG::StepHandler::eventHandler(), because
   * it is not virtual in ThePEG::StepHandler. This is ok, because this
   * method would give a null-pointer at some stages, whereas this method
   * gives access to the explicitely copied pointer (in doinitrun()) 
   * to the ThePEG::EventHandler.
   */
  inline tEHPtr eventHandler() const;

  /**
   * Return the ThePEG::StandardEventHandler assigned to this handler.
   */
  inline tStdEHPtr stdeventHandler() const;

  /**
   * Return the sampler assigned to this handler.
   */
  inline tSamplerPtr sampler();

  /**
   * Return the sampler assigned to this handler.
   */
  inline tcSamplerPtr sampler() const;

  /**
   * The pair of incoming particle types obtained via the EventHandler
   */
  inline const cPDPair & incoming() const;

  /**
   * Access the luminosity function via the EventHandler.
   */
  inline const LuminosityFunction & lumiFn() const;

  /**
   * The number of phase space dimensions used by the luminosity
   * function. Calls the corresponding StandardEventHandler method.
   */
  inline int lumiDim() const;

  /**
   * Return the number of separate bins of StandardXComb objects to
   * sample.
   */
  int nBins() const;

  /**
   * Return the number of phase space dimensions needed for the
   * sampling of indicated bin of StandardXComb objects.
   */
  inline int maxDim(int bin) const;

  /**
   * The number of dimensions of the basic phase space to generate
   * sub-processes in for a given bin of StandardXComb objects.
   */
  inline int nDim(int bin) const;

  /**
   * Return the maximum number attemts allowed to select a sub-process
   * for each event. Calls the corresponding StandardEventHandler method.
   */
  inline long maxLoop() const;


protected:

  /**
   * Generate a phase space point and return the corresponding cross
   * section. Is called from sSigDR(const vector<double> &).
   * @param ll a pair of doubles giving the logarithms of the (inverse
   * energy fractions of the maximum CMS energy of the incoming
   * particles.
   * @param maxS the maximum squared CMS energy of the incoming particles.
   * @param ibin the preselected bin of StandardXComb objects to choose
   * sub-process from
   * @param nr the number of random numbers availiable in \a r.
   * @param r an array of random numbers to be used to generate a
   * phase-space point.
   */
  virtual CrossSection dSigDR(const pair<double,double> ll, Energy2 maxS,
			      int ibin, int nr, const double * r);


  /**
   * Select an StandardXComb. Given a preselected bin, \a ibin of
   * StandardXComb objects pick one to generate the corresponding
   * sub-process with the given \a weight.
   */
  tStdXCombPtr select(int bin, double weight);

  /**
   * Create and add <code>StandardXComb</code> objects.
   *
   * @param maxEnergy the maximum CMS energy of the incoming particles.
   * @param sub a pointer to the SubProcessHandler object.
   * @param extractor a pointer to the PartonExtractor object.
   * @param cuts a pointer to the Cuts object.
   * @param ckkw a currently empty pointer to a CascadeHandler to be used for CKKW reweighting.
   * @param me a pointer to the MEBase object.
   * @param pBins a pair of <code>PartonBin</code>s describing the
   * partons extracted from the particles
   */
  void addME(Energy maxEnergy, tSubHdlPtr sub, tPExtrPtr extractor, tCutsPtr cuts, 
	     tCascHdlPtr ckkw, tMEPtr me, const PBPair & pBins);

  /**
   * Return the vector of StandardXComb objects.
   */
  inline const XVector & xCombs() const;

  /**
   * Return the vector of StandardXComb objects.
   */
  inline XVector & xCombs();

  /**
   * Return the vector of cross sections.
   */
  inline const XSVector & xSecs() const;

  /**
   * Return the vector of cross sections.
   */
  inline XSVector & xSecs();

  /**
   * Return the strategy to be used when sampling different StandardXComb
   * objects.
   * @return 0 if all StandardXComb objects are sampled together. 1 if
   * all StandardXComb objects which have the same matrix element object are
   * sampled together. 2 if all StandardXComb objects are sampled separately.
   */
  inline int binStrategy() const;


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

private:

  /**
   * Access the list of sub-process handlers.
   */
  inline const SubHandlerList & subProcesses() const;

  /**
   * Access the list of sub-process handlers.
   */
  inline SubHandlerList & subProcesses();

  /**
   *  Method to calculate the individual probabilities for N scatters in the event.
   *  @param UEXSecs is(are) the inclusiv cross section(s) for the UE process(es).
   */
  void Probs(XSVector UEXSecs);
  
  /**
   * Return the value of the Overlap function A(b) for a given impact 
   * parameter \a b.
   *  @param b impact parameter
   *  @return inverse area.
   */
  InvArea OverlapFunction(Length b);

  /**
   *  Return n!
   */
  double factorial (unsigned int n);

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MPIHandler> initMPIHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MPIHandler & operator=(const MPIHandler &);

  /**
   * The phase space sampler responsible for generating phase space
   * points according to the cross section given by this handler.
   */
  SamplerPtr theSampler;

  /**
   * A pointer to the EventHandler that calls us. Has to be saved, because the
   * method eventHandler() inherited from ThePEG::StepHandler returns a null-pointer
   * sometimes. Leif changed that in r1053 so that a valid pointer is present, when
   * calling doinitrun().
   */
  tEHPtr theHandler;

  /**
   * The kinematical cuts used for this collision handler.
   */
  CutsPtr theCuts;


  /**
   * The list of <code>SubProcessHandler</code>s.
   */
  SubHandlerList theSubProcesses;

  /**
   * The StandardXComb objects.
   */
  XVector theXCombs;

  /**
   * The (incrementally summed) cross sections associated with the
   * StandardXComb objects for the last selected phase space point.
   */
  XSVector theXSecs;

  /**
   * The strategy to be used when sampling different StandardXComb
   * objects. 0 means all StandardXComb objects are sampled
   * together. 1 means all StandardXComb objects which have the same
   * matrix element object are sampled together. 2 means all
   * StandardXComb objects are sampled separately.
   */
  int theBinStrategy;

  /**
   * The map used to store all XBins with the same matrix element for
   * option 1 in theBinStrategy.
   */
  MEXMap theMEXMap;


  /**
   * The number of degrees of freedom needed to generate the phase
   * space for the different bins.
   */
  vector<int> theMaxDims;

  /**
   * A ThePEG::Selector where the individual Probabilities P_N are stored
   * and the actual Multiplicities can be selected.
   */
  Selector<unsigned int> theMultiplicities;

  /**
   * Switch to be set from outside to determine the algorithm used for 
   * UE activity.
   */
  int theJmueo;

  /**
   * Inverse Radius squared \f$ (\mu^2) \f$. Used inside the overlap function.  
   */
  Energy2 theRadius;

protected:

  /** @cond EXCEPTIONCLASSES */

  /**
   * Exception class used by the MultipleInteractionHandler, when something
   * during initialization went wrong.
   * \todo understand!!!
   */
  class InitError: public Exception {};

  /** @endcond */

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MPIHandler. */
template <>
struct BaseClassTrait<Herwig::MPIHandler,1> {
  /** Typedef of the first base class of MPIHandler. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MPIHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MPIHandler>
  : public ClassTraitsBase<Herwig::MPIHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MPIHandler"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MPIHandler class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwMPI.so"; }
};

/** @endcond */

}

namespace Herwig {

  /**
   * A struct for the eikonalization of the inclusive cross section. 
   */
  struct Eikonalization {

    /**
     *  The constructor
     *  @param handler is the pointer to the MPIHandler to get access to 
     *  MPIHandler::OverlapFunction and member variables of the MPIHandler.
     *  @param xsec is the cross section to be eikonalized.
     *  @param option is a flag, whether the inelastic or the total 
     *  cross section should be returned (-2 or -1). For option = N > 0 the integrand
     *  is N*(A(b)*sigma)^N/N! exp(-A(b)*sigma) this is the P_N*sigma where
     *  P_N is the Probability of having exactly N interaction (including the hard one)
     *  This is equation 14 from "Jimmy4: Multiparton Interactions in HERWIG for the LHC"
     */
    inline Eikonalization(tMPIHPtr handler, CrossSection xsec, int option);

    /**
     * Get the function value
     */
    Length operator ()(Length argument) const;
    typedef Length ValType;
    typedef Length ArgType;
    /**
     * Pointer to the Handler that calls this integrand
     */
    tMPIHPtr theHandler;

    /**
     * The cross section that is eikonalized
     */
    CrossSection theUneikXSec;

    /**
     * A flag to switch between the calculation of total and inelastig cross section
     * or calculations for the individual probabilities. See the constructor
     */
    int theoption;
  };
}


#include "MPIHandler.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MPIHandler.tcc"
#endif

#endif /* HERWIG_MPIHandler_H */
