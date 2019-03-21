// -*- C++ -*-
//
// MPIHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MPIHandler_H
#define HERWIG_MPIHandler_H
//
// This is the declaration of the MPIHandler class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "Herwig/Utilities/GSLBisection.h"
//#include "Herwig/Utilities/GSLMultiRoot.h"
#include "Herwig/Utilities/GSLIntegrator.h"
#include "Herwig/Shower/UEBase.h"

#include <cassert>
#include "ProcessHandler.h"
#include "MPIHandler.fh"


namespace Herwig {
using namespace ThePEG;

  /** \ingroup UnderlyingEvent
   * \class MPIHandler
   * This class is responsible for generating additional 
   * semi hard partonic interactions.
   * 
   * \author Manuel B\"ahr
   *
   * @see \ref MPIHandlerInterfaces "The interfaces"
   * defined for MPIHandler.
   * @see ProcessHandler
   * @see ShowerHandler
   * @see HwRemDecayer
   */

class MPIHandler: public UEBase {

  /**
   *  Maximum number of scatters
   */
  static const unsigned int maxScatters_ = 99;

  /**
   * Class for the integration is a friend to access private members
   */
  friend struct Eikonalization;
  friend struct TotalXSecBisection;
  friend struct slopeAndTotalXSec;
  friend struct slopeInt;
  friend struct slopeBisection;

public:

  /** A vector of <code>SubProcessHandler</code>s. */
  typedef vector<SubHdlPtr> SubHandlerList;

  /** A vector of <code>Cut</code>s. */
  typedef vector<CutsPtr> CutsList;

  /** A vector of <code>ProcessHandler</code>s. */
  typedef vector<ProHdlPtr> ProcessHandlerList;

  /** A vector of cross sections. */
  typedef vector<CrossSection> XSVector;

  /** A pair of multiplicities: hard, soft. */
  typedef pair<unsigned int, unsigned int> MPair;

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MPIHandler(): softMult_(0), identicalToUE_(-1), 
		PtOfQCDProc_(-1.0*GeV), Ptmin_(-1.0*GeV), 
		hardXSec_(0*millibarn), softXSec_(0*millibarn), 
		totalXSecExp_(0*millibarn),
		softMu2_(ZERO), beta_(100.0/GeV2), 
		algorithm_(2), numSubProcs_(0), 
		colourDisrupt_(0.0), softInt_(true), twoComp_(true),
		DLmode_(2), avgNhard_(0.0), avgNsoft_(0.0),
                energyExtrapolation_(2), EEparamA_(0.6*GeV),
                EEparamB_(37.5*GeV), refScale_(7000.*GeV),
		pT0_(3.11*GeV), b_(0.21) {}

  /**
   * The destructor.
   */
  virtual ~MPIHandler(){}
  //@}

public:

  /** @name Methods for the MPI generation. */
  //@{

  /*
   * @return true if for this beam setup MPI can be generated
   */
  virtual bool beamOK() const;

  /**
   * Return true or false depending on whether soft interactions are enabled.
   */
  virtual bool softInt() const {return softInt_;}

  /**
   * Get the soft multiplicity from the pretabulated multiplicity
   * distribution. Generated in multiplicity in the first place.
   * @return the number of extra soft events in this collision
   */
  virtual unsigned int softMultiplicity() const {return softMult_;} 

  /**
   * Sample from the pretabulated multiplicity distribution.
   * @return the number of extra events in this collision
   */
  virtual unsigned int multiplicity(unsigned int sel=0); 

  /**
   * Select a StandardXComb according to it's weight
   * @return that StandardXComb Object
   * @param sel is the subprocess that should be returned,
   * if more than one is specified.
   */
  virtual tStdXCombPtr generate(unsigned int sel=0);
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
  virtual void initialize();

  /**
   * Finalize this Multiple Interaction handler and all related objects needed to
   * generate additional events.
   */
  virtual void finalize();

  /**
   * Clean up the XCombs from our subprocesses after each event.
   * ThePEG cannot see them, so the usual cleaning misses these.
   */
  virtual void clean();

  /**
   * Write out accumulated statistics about integrated cross sections.
   */
  void statistics() const;

  /**
   * The level of statistics. Controlls the amount of statistics
   * written out after each run to the <code>EventGenerator</code>s
   * <code>.out</code> file. Simply the EventHandler method is called here.
   */
  int statLevel() const {return eventHandler()->statLevel();}

  /**
   * Return the hard cross section above ptmin
   */
  CrossSection hardXSec() const { return hardXSec_; }

  /**
   * Return the soft cross section below ptmin
   */
  CrossSection softXSec() const { return softXSec_; }

  /**
   * Return the inelastic cross section
   */
  CrossSection inelasticXSec() const { return inelXSec_; }

  /** @name Simple access functions. */
  //@{

  /**
   * Return the ThePEG::EventHandler assigned to this handler.
   * This methods shadows ThePEG::StepHandler::eventHandler(), because
   * it is not virtual in ThePEG::StepHandler. This is ok, because this
   * method would give a null-pointer at some stages, whereas this method
   * gives access to the explicitely copied pointer (in initialize()) 
   * to the ThePEG::EventHandler.
   */
  tEHPtr eventHandler() const {return theHandler;}

  /**
   * Return the current handler
   */
  static const MPIHandler * currentHandler() {
    return currentHandler_;
  }

  /**
   * Return theAlgorithm.
   */
  virtual int Algorithm() const {return algorithm_;}

  /**
   * Return the ptmin parameter of the model
   */
  virtual Energy Ptmin() const {
    if(Ptmin_ > ZERO)
      return Ptmin_;
    else
      throw Exception() << "MPIHandler::Ptmin called without initialize before"
			<< Exception::runerror;
  }

  /**
   * Return the slope of the soft pt spectrum as calculated.
   */
  virtual InvEnergy2 beta() const {
    if(beta_ != 100.0/GeV2)
      return beta_;
    else
      throw Exception() << "MPIHandler::beta called without initialization"
			<< Exception::runerror;
  }

  /**
   * Return the pt Cutoff of the Interaction that is identical to the UE
   * one.
   */
  virtual Energy PtForVeto() const {return PtOfQCDProc_;}
  
  /**
   * Return the number of additional "hard" processes ( = multiple
   * parton scattering)
   */
  virtual unsigned int additionalHardProcs() const {return numSubProcs_-1;}

  /**
   * Return the fraction of colour disrupted connections to the
   * suprocesses.
   */
  virtual double colourDisrupt() const {return colourDisrupt_;}

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
   * Access the list of sub-process handlers.
   */
  const SubHandlerList & subProcesses() 
    const {return theSubProcesses;}

  /**
   * Access the list of sub-process handlers.
   */
  SubHandlerList & subProcesses() {return theSubProcesses;}

  /**
   * Access the list of cuts.
   */
  const CutsList & cuts() const {return theCuts;}

  /**
   * Access the list of cuts.
   */
  CutsList & cuts() {return theCuts;}

  /**
   * Access the list of sub-process handlers.
   */
  const ProcessHandlerList & processHandlers() 
    const {return theProcessHandlers;}

  /**
   * Access the list of sub-process handlers.
   */
  ProcessHandlerList & processHandlers() {return theProcessHandlers;}


  /**
   *  Method to calculate the individual probabilities for N scatters in the event.
   *  @param UEXSecs is(are) the inclusiv cross section(s) for the UE process(es).
   */
  void Probs(XSVector UEXSecs);

  /**
   * Debug method to check the individual probabilities.
   * @param filename is the file the output gets written to
   */
  void MultDistribution(string filename) const;
  
  /**
   * Return the value of the Overlap function A(b) for a given impact 
   * parameter \a b.
   *  @param b impact parameter
   *  @param mu2 = inv hadron radius squared. 0 will use the value of
   *  invRadius_
   *  @return inverse area.
   */
  InvArea OverlapFunction(Length b, Energy2 mu2=ZERO) const;

  /**
   * Method to calculate the poisson probability for expectation value
   * \f$<n> = A(b)\sigma\f$, and multiplicity N.
   */
  double poisson(Length b, CrossSection sigma, 
		 unsigned int N, Energy2 mu2=ZERO) const;

  /**
   *  Return n!
   */
  double factorial (unsigned int n) const;

  /**
   * Returns the total cross section for the current CMenergy.  The
   * decision which parametrization will be used is steered by a
   * external parameter of this class.
   */
  CrossSection totalXSecExp() const;

  /**
   * Difference of the calculated total cross section and the
   * experimental one from totalXSecExp.
   * @param softXSec = the soft cross section that is used
   * @param softMu2 = the soft radius, if 0 the hard radius will be used
   */
  CrossSection totalXSecDiff(CrossSection softXSec, 
			     Energy2 softMu2=ZERO) const;

  /**
   * Difference of the calculated elastic slope and the
   * experimental one from slopeExp.
   * @param softXSec = the soft cross section that is used
   * @param softMu2 = the soft radius, if 0 the hard radius will be used
   */
  InvEnergy2 slopeDiff(CrossSection softXSec, 
			 Energy2 softMu2=ZERO) const;

  /**
   * Returns the value of the elastic slope for the current CMenergy.
   * The decision which parametrization will be used is steered by a
   * external parameter of this class.
   */
  InvEnergy2 slopeExp() const;


  /**
   * Calculate the minimal transverse momentum from the extrapolation
   */
  void overrideUECuts();


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
  MPIHandler & operator=(const MPIHandler &) = delete;

  /**
   * A pointer to the EventHandler that calls us. Has to be saved, because the
   * method eventHandler() inherited from ThePEG::StepHandler returns a null-pointer
   * sometimes. Leif changed that in r1053 so that a valid pointer is present, when
   * calling doinitrun().
   */
  tEHPtr theHandler;

  /**
   * The list of <code>SubProcessHandler</code>s.
   */
  SubHandlerList theSubProcesses;

  /**
   * The kinematical cuts used for this collision handler.
   */
  CutsList theCuts;

  /**
   * List of ProcessHandler used to sample different processes independently
   */
  ProcessHandlerList theProcessHandlers;

  /**
   * A ThePEG::Selector where the individual Probabilities P_N are stored
   * and the actual Multiplicities can be selected.
   */
  Selector<MPair> theMultiplicities;

  /**
   * Variable to store the soft multiplicity generated for a event. This
   * has to be stored as it is generated at the time of the hard
   * additional interactions but used later on.
   */
  unsigned int softMult_;

  /**
   * Variable to store the multiplicity of the second hard process
   */
  vector<int> additionalMultiplicities_;

  /**
   * Variable to store the information, which process is identical to
   * the UE one (QCD dijets).
   * 0 means "real" hard one
   * n>0 means the nth additional hard scatter
   * -1 means no one!
   */
  int identicalToUE_;

  /**
   * Variable to store the minimal pt of the process that is identical
   * to the UE one. This only has to be set, if it can't be determined
   * automatically (i.e. when reading QCD LesHouches files in).
   */
  Energy PtOfQCDProc_;

  /**
   * Variable to store the parameter ptmin
   */
  Energy Ptmin_;

  /**
   * Variable to store the hard cross section above ptmin
   */
  CrossSection hardXSec_;

  /**
   * Variable to store the final soft cross section below ptmin
   */
  CrossSection softXSec_;

  /**
   * Variable to store the inelastic cross section
   */
  CrossSection inelXSec_;

  /**
   * Variable to store the total pp cross section (assuming rho=0!) as
   * measured at LHC. If this variable is set, this value is used in the
   * subsequent run instead of any of the Donnachie-Landshoff
   * parametrizations.
   */
  CrossSection totalXSecExp_;

  /**
   * Variable to store the soft radius, that is calculated during
   * initialization for the two-component model.
   */
  Energy2 softMu2_;

  /**
   * slope to the non-perturbative pt spectrum: \f$d\sigma/dp_T^2 = A \exp
   * (- beta p_T^2)\f$. Its value is determined durint initialization.
   */
  InvEnergy2 beta_;
  /**
   * Switch to be set from outside to determine the algorithm used for 
   * UE activity.
   */
  int algorithm_;

  /**
   * Inverse hadron Radius squared \f$ (\mu^2) \f$. Used inside the overlap function.  
   */
  Energy2 invRadius_;

  /**
   * Member variable to store the actual number of separate SubProcesses
   */ 
  unsigned int numSubProcs_;

  /**
   * Variable to store the relative number of colour disrupted
   * connections to additional subprocesses. This variable is used in
   * Herwig::HwRemDecayer but store here, to have access to all
   * parameters through one Object.
   */
  double colourDisrupt_;

  /** 
   * Flag to store whether soft interactions, i.e. pt < ptmin should be
   * simulated.
   */
  bool softInt_;

  /** 
   * Flag to steer wheather the soft part has a different radius, that
   * will be dynamically fixed.
   */
  bool twoComp_;
  
  /**
   * Switch to determine which Donnachie & Landshoff parametrization
   * should be used.
   */
  unsigned int DLmode_;

  /**
   * Variable to store the average hard multiplicity.
   */
  double avgNhard_;

  /**
   * Variable to store the average soft multiplicity.
   */
  double avgNsoft_;

  /**
   * The current handler
   */
  static MPIHandler * currentHandler_;

  /**
   * Flag to store whether to calculate the minimal UE pt according to an
   * extrapolation formula or whether to use MPIHandler:Cuts[0]:OneCuts[0]:MinKT
   */
  unsigned int energyExtrapolation_;

  /**
   * Parameters for the energy extrapolation formula
   */
  Energy EEparamA_;
  Energy EEparamB_;
  Energy refScale_;
  Energy pT0_;
  double b_;

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
  static string library() { return "JetCuts.so SimpleKTCut.so HwMPI.so"; }
};

/** @endcond */

}

namespace Herwig {

  /**
   * A struct for the 2D root finding that is necessary to determine the
   * soft cross section and the soft radius that is needed to describe
   * the total cross section correctly.
   * NOT IN USE CURRENTLY
   */
  struct slopeAndTotalXSec : public GSLHelper<CrossSection, CrossSection> {

  public:

    /**
     *  Constructor
     */
    slopeAndTotalXSec(tcMPIHPtr handler): handler_(handler) {}

    /** second argument type */
    typedef Energy2 ArgType2;

    /** second value type */
    typedef InvEnergy2 ValType2;

    /** first element of the vector like function to find root for 
     * @param softXSec soft cross-section
     * @param softMu2 \f$\mu^2\f$ 
     */
    CrossSection f1(ArgType softXSec, ArgType2 softMu2) const {
      return handler_->totalXSecDiff(softXSec, softMu2);
    }

    /** second element of the vector like function to find root for 
     * @param softXSec soft cross-section
     * @param softMu2 \f$\mu^2\f$ 
     */
    InvEnergy2 f2(ArgType softXSec, ArgType2 softMu2) const {
      return handler_->slopeDiff(softXSec, softMu2);
    }

    /** provide the actual units of use */
    virtual ValType vUnit() const {return 1.0*millibarn;}
    
    /** otherwise rounding errors may get significant */
    virtual ArgType aUnit() const {return 1.0*millibarn;}

    /** provide the actual units of use */
    ValType2 vUnit2() const {return 1.0/GeV2;}
    
    /** otherwise rounding errors may get significant */
    ArgType2 aUnit2() const {return GeV2;}

  private: 

    /**
     *  Pointer to the handler
     */
    tcMPIHPtr handler_;

  };
  
  /**
   * A struct for the root finding that is necessary to determine the
   * slope of the soft pt spectrum to match the soft cross section
   */
  struct betaBisection : public GSLHelper<Energy2, InvEnergy2>{
  public:
    /**
     * Constructor.
     * @param soft = soft cross section, i.e. the integral of the soft
     * pt spectrum f(u=p_T^2) = dsig exp(-beta*u/u_min)
     * @param dsig = dsigma_hard/dp_T^2 at the p_T cutoff
     * @param ptmin = p_T cutoff
     */
    betaBisection(CrossSection soft, DiffXSec dsig, Energy ptmin) 
      : softXSec_(soft), dsig_(dsig), ptmin_(ptmin) {}
   
    /**
     * Operator that is used inside the GSLBisection class
     */
    virtual Energy2 operator ()(InvEnergy2 beta) const
    {
      if( fabs(beta*GeV2) < 1.E-4 )
	beta = (beta > ZERO) ? 1.E-4/GeV2 : -1.E-4/GeV2;

      return (exp(beta*sqr(ptmin_)) - 1.0)/beta - softXSec_/dsig_;
    }

    /** provide the actual units of use */
    virtual ValType vUnit() const {return 1.0*GeV2;}

    /** provide the actual units of use */
    virtual ArgType aUnit() const {return 1.0/GeV2;}

  private: 

    /** soft cross section */
    CrossSection softXSec_;

    /** dsigma/dp_T^2 at ptmin */
    DiffXSec dsig_;

    /** pt cutoff */
    Energy ptmin_;
  };

  /**
   * A struct for the root finding that is necessary to determine the
   * soft cross section and soft mu2 that are needed to describe the
   * total cross section AND elastic slope correctly.
   */
  struct slopeBisection : public GSLHelper<InvEnergy2, Energy2> {
  public:
    /** Constructor */
    slopeBisection(tcMPIHPtr handler) : handler_(handler) {}

    /** 
     * Return the difference of the calculated elastic slope to the
     * experimental one for a given value of the soft mu2. During that,
     * the soft cross section get fixed.
     */
    InvEnergy2 operator ()(Energy2 arg) const;
    
    /** Return the soft cross section that has been calculated */
    CrossSection softXSec() const {return softXSec_;}

  private:
    /** const pointer to the MPIHandler to give access to member functions.*/
    tcMPIHPtr handler_;
    /** soft cross section that is determined on the fly.*/
    mutable CrossSection softXSec_;
  };

  /**
   * A struct for the root finding that is necessary to determine the
   * soft cross section that is needed to describe the total cross
   * section correctly.
   */
  struct TotalXSecBisection : public GSLHelper<CrossSection, CrossSection> {
  public:

    /**
     *  Constructor
     * @param handler The handler
     * @param softMu2 \f$\mu^2\f$
     */
    TotalXSecBisection(tcMPIHPtr handler, Energy2 softMu2=ZERO): 
      handler_(handler), softMu2_(softMu2) {}

    /**
     *  operator to return the cross section
     * @param argument input cross section
     */
    CrossSection operator ()(CrossSection argument) const {
      return handler_->totalXSecDiff(argument, softMu2_);
    }

    /** provide the actual units of use */
    virtual ValType vUnit() const {return 1.0*millibarn;}
    
    /** otherwise rounding errors may get significant */
    virtual ArgType aUnit() const {return 1.0*millibarn;}

  private: 

    /**
     *  The handler
     */
    tcMPIHPtr handler_;

    /**
     *  \f$\mu^2\f$
     */
    Energy2 softMu2_;

  };

  /**
   *  Typedef for derivative of the length
   */
  typedef Qty<1,-2,0> LengthDiff;

  /**
   *  A struct for the integrand for the slope
   */
  struct slopeInt : public GSLHelper<LengthDiff, Length>{

  public:
    /** Constructor 
     * @param handler The handler
     * @param hard The hard cross section
     * @param soft The soft cross section
     * @param softMu2 \f$\mu^2\f$
     */
    slopeInt(tcMPIHPtr handler, CrossSection hard, 
	      CrossSection soft=0*millibarn, Energy2 softMu2=ZERO)
      : handler_(handler), hardXSec_(hard), 
	softXSec_(soft), softMu2_(softMu2) {}

    /**
     *  Operator to return the answer
     * @param arg The argument
     */
    ValType operator ()(ArgType arg) const;

  private:

    /**
     * Pointer to the Handler that calls this integrand
     */
    tcMPIHPtr handler_;

    /**
     * The hard cross section to be eikonalized
     */
    CrossSection hardXSec_;

    /**
     * The soft cross section to be eikonalized. Default is zero
     */
    CrossSection softXSec_;

    /**
     * The inv radius^2 of the soft interactions.
     */
    Energy2 softMu2_;

  };

  /**
   * A struct for the eikonalization of the inclusive cross section. 
   */
  struct Eikonalization : public GSLHelper<Length, Length>{

    /**
     *  The constructor
     *  @param handler is the pointer to the MPIHandler to get access to 
     *  MPIHandler::OverlapFunction and member variables of the MPIHandler.
     *  @param option is a flag, whether the inelastic or the total 
     *  @param handler The handler
     *  @param hard The hard cross section
     *  @param soft The soft cross section
     *  @param softMu2 \f$\mu^2\f$
     *  cross section should be returned (-2 or -1). For option = N > 0 the integrand
     *  is N*(A(b)*sigma)^N/N! exp(-A(b)*sigma) this is the P_N*sigma where
     *  P_N is the Probability of having exactly N interaction (including the hard one)
     *  This is equation 14 from "Jimmy4: Multiparton Interactions in HERWIG for the LHC"
     */
    Eikonalization(tcMPIHPtr handler, int option, CrossSection hard, 
		   CrossSection soft=0*millibarn, Energy2 softMu2=ZERO) 
      : theHandler(handler), theoption(option), hardXSec_(hard), 
	softXSec_(soft), softMu2_(softMu2) {}

    /**
     * Get the function value
     */
    Length operator ()(Length argument) const;

  private:
    /**
     * Pointer to the Handler that calls this integrand
     */
    tcMPIHPtr theHandler;

    /**
     * A flag to switch between the calculation of total and inelastic cross section
     * or calculations for the individual probabilities. See the constructor
     */
    int theoption;

    /**
     * The hard cross section to be eikonalized
     */
    CrossSection hardXSec_;

    /**
     * The soft cross section to be eikonalized. Default is zero
     */
    CrossSection softXSec_;

    /**
     * The inv radius^2 of the soft interactions.
     */
    Energy2 softMu2_;

  };
}


#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MPIHandler.tcc"
#endif

#endif /* HERWIG_MPIHandler_H */
