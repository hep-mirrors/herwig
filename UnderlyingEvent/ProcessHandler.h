// -*- C++ -*-
//
// ProcessHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ProcessHandler_H
#define HERWIG_ProcessHandler_H
//
// This is the declaration of the ProcessHandler class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/CascadeHandler.h"

#include <cassert>

#include "ProcessHandler.fh"
#include "stat.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup UnderlyingEvent
   * \class ProcessHandler
   * This class is for handling the sampling of 
   * semi hard partonic interactions. If several types of partonic interactions
   * are needed. Each of them has it's own ProcessHandler. A reference to them is
   * stored in MPIHandler, which administers everything.
   * 
   * \author Manuel B\"ahr
   *
   * @see \ref ProcessHandlerInterfaces "The interfaces"
   * defined for ProcessHandler.
   * @see MPISampler
   * @see MPIHandler
   */

class ProcessHandler: public Interfaced, public LastXCombInfo<> {

public:

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
  ProcessHandler();

  /**
   * The copy constructor.
   */
  ProcessHandler(const ProcessHandler &);

  /**
   * The destructor.
   */
  virtual ~ProcessHandler();
  //@}

public:

  /** @name Methods for the MPI generation. */
  //@{

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
  void initialize(tSubHdlPtr sub, tCutsPtr cuts, tEHPtr eh);

  /**
   * Return the integrated cross section.
   */
  CrossSection integratedXSec() const;

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
   * The level of statistics. Controlls the amount of statistics
   * written out after each run to the <code>EventGenerator</code>s
   * <code>.out</code> file. Simply the EventHandler method is called here.
   */
  inline int statLevel() const;

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
   * Return the ThePEG::EventHandler assigned to this handler.
   * This methods shadows ThePEG::StepHandler::eventHandler(), because
   * it is not virtual in ThePEG::StepHandler. This is ok, because this
   * method would give a null-pointer at some stages, whereas this method
   * gives access to the explicitely copied pointer (in doinitrun()) 
   * to the ThePEG::EventHandler.
   */
  inline tEHPtr eventHandler() const;

  /**
   * Return the sampler assigned to this handler.
   */
  inline tSamplerPtr sampler();

  /**
   * Return the sampler assigned to this handler.
   */
  inline tcSamplerPtr sampler() const;

  /**
   * Return a reference to the Cuts of this
   * MultipleInteractionHandler. Note that these cuts may be overridden by the
   * SubProcess chosen.
   */
  inline tCutsPtr cuts() const;

  /**
   * Access the sub-process handler.
   */
  inline tSubHdlPtr subProcess();


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
  static ClassDescription<ProcessHandler> initProcessHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ProcessHandler & operator=(const ProcessHandler &);

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
  tCutsPtr theCuts;//used a transient pointer, 
  //because regular pointer is already in MPIHandler


  /**
   * The SubProcessHandler that is connected to this ProcessHandler.
   */
  tSubHdlPtr theSubProcess;//used a transient pointer, 
  //because regular pointer is already in MPIHandler

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
 *  base classes of ProcessHandler. */
template <>
struct BaseClassTrait<Herwig::ProcessHandler,1> {
  /** Typedef of the first base class of ProcessHandler. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ProcessHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ProcessHandler>
  : public ClassTraitsBase<Herwig::ProcessHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ProcessHandler"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the ProcessHandler class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "JetCuts.so SimpleKTCut.so HwMPI.so"; }
};

/** @endcond */

}

#include "ProcessHandler.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ProcessHandler.tcc"
#endif

#endif /* HERWIG_ProcessHandler_H */
