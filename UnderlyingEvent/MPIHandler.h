// -*- C++ -*-
//
// MPIHandler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/PDT/StandardMatchers.h"
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

class MPIHandler: public Interfaced {

  /**
   * Class for the integration is a friend to access private members
   */
  friend class Eikonalization;


public:

  /** A vector of <code>SubProcessHandler</code>s. */
  typedef vector<SubHdlPtr> SubHandlerList;

  /** A vector of <code>Cut</code>s. */
  typedef vector<CutsPtr> CutsList;

  /** A vector of <code>ProcessHandler</code>s. */
  typedef vector<ProHdlPtr> ProcessHandlerList;

  /** A vector of cross sections. */
  typedef vector<CrossSection> XSVector;

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

  /*
   * @return true if for this beam setup MPI can be generated
   */
  inline bool beamOK() const;

  /**
   * Sample from the pretabulated multiplicity distribution.
   * @return the number of extra events in this collision
   */
  inline unsigned int multiplicity() const;

  /**
   * Select a StandardXComb according to it's weight
   * @return that StandardXComb Object
   * @param sel is the subprocess that should be returned,
   * if more than one is specified.
   */
  inline tStdXCombPtr generate(unsigned int sel=0);
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
   * Finalize this Multiple Interaction handler and all related objects needed to
   * generate additional events.
   */
  void finalize();

  /**
   * Write out accumulated statistics about intergrated cross sections
   * and stuff.
   */
  void statistics(string file) const;

  /**
   * The level of statistics. Controlls the amount of statistics
   * written out after each run to the <code>EventGenerator</code>s
   * <code>.out</code> file. Simply the EventHandler method is called here.
   */
  inline int statLevel() const;

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
  inline tEHPtr eventHandler() const;

  /**
   * Return theAlgorithm.
   */
  inline int Algorithm() const;

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
   * Access the list of cuts.
   */
  inline const CutsList & cuts() const;

  /**
   * Access the list of cuts.
   */
  inline CutsList & cuts();

  /**
   * Access the list of sub-process handlers.
   */
  inline const ProcessHandlerList & processHandlers() const;

  /**
   * Access the list of sub-process handlers.
   */
  inline ProcessHandlerList & processHandlers();


  /**
   *  Method to calculate the individual probabilities for N scatters in the event.
   *  @param UEXSecs is(are) the inclusiv cross section(s) for the UE process(es).
   */
  void Probs(XSVector UEXSecs);
  
  /**
   * Method to calculate the poisson probability for expectation value
   * <n> = A(b)*sigma, and multiplicity mult.
   */
  double poisson(Length b, CrossSection sigma, unsigned int mult) const;

  /**
   * Return the value of the Overlap function A(b) for a given impact 
   * parameter \a b.
   *  @param b impact parameter
   *  @return inverse area.
   */
  InvArea OverlapFunction(Length b) const;

  /**
   *  Return n!
   */
  double factorial (unsigned int n) const;

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
  Selector<unsigned int> theMultiplicities;

  /**
   * Switch to be set from outside to determine the algorithm used for 
   * UE activity.
   */
  int theAlgorithm;

  /**
   * Inverse hadron Radius squared \f$ (\mu^2) \f$. Used inside the overlap function.  
   */
  Energy2 theInvRadius;

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

    /** Resulting integral is of type Length, because the integrand has no dimension */
    typedef Length ValType;

    /** Integration variable is of type Length */
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
