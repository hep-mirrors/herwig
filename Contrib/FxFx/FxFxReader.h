// -*- C++ -*-
//
// FxFxReader.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_FxFxReader_H
#define THEPEG_FxFxReader_H
// This is the declaration of the FxFxReader class.

#include "FxFx.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Utilities/ObjectIndexer.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Utilities/XSecStat.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "ThePEG/PDF/PartonBin.fh"
#include "ThePEG/MatrixElement/ReweightBase.h"
#include "FxFxEventHandler.fh"
#include "FxFxReader.fh"
#include "ThePEG/Utilities/CFile.h"
#include <cstdio>
#include <cstring>

namespace ThePEG {

/**
 * FxFxReader is an abstract base class to be used for objects
 * which reads event files or streams from matrix element
 * generators. Derived classes must at least implement the open() and
 * doReadEvent() methods to read in information about the whole run into
 * the HEPRUP variable and next event into the HEPEUP variable
 * respectively. Also the close() function to close the file or stream
 * read must be implemented. Although these functions are named as if
 * we are reading from event files, they could just as well implement
 * the actual generation of events.
 *
 * After filling the HEPRUP and HEPEUP variables, which are protected
 * and easily accesible from the sub-class, this base class will then
 * be responsible for transforming this data to the ThePEG Event
 * record in the getEvent() method. <code>FxFxReader</code>s can
 * only be used inside FxFxEventHandler objects.
 *
 * In the initialization the virtual open() and scan() functions are
 * called. Here the derived class must provide the information about
 * the processes in the variables corresponding to the HEPRUP common
 * block. Note that the IDWTUP is required to be +/- 1, and sub
 * classes are required to change the information accordingly to
 * ensure the correct corss section sampling. Note also that the
 * controlling FxFxEventHandler may choose to generate weighted
 * events even if IDWTUP is 1.
 *
 * Note that the information given per process in e.g. the XSECUP and
 * XMAXUP vectors is not used by the FxFxEventHandler and by
 * default the FxFxReader is not assumed to be able to actively
 * choose between the sub-processes. Instead, the
 * FxFxEventHandler can handle several FxFxReader objects
 * and choose between them. However, a sub-class of FxFxReader
 * may set the flag isActive, in which case it is assumed to be able
 * to select between its sub-processes itself.
 *
 * The FxFxReader may be assigned a number ReweightBase objects
 * which either completely reweights the events produced (in the
 * reweights vector), or only biases the selection without influencing
 * the cross section (in the preweights vector). Note that it is the
 * responsibility of a sub-class to call the reweight() function and
 * multiply the weight according to its return value (typically done
 * in the readEvent() function).
 *
 * @see \ref FxFxReaderInterfaces "The interfaces"
 * defined for FxFxReader.
 * @see Event
 * @see FxFxEventHandler
 */
class FxFxReader: public HandlerBase, public LastXCombInfo<> {

  /**
   * FxFxEventHandler should have access to our private parts.
   */
  friend class FxFxEventHandler;

  /**
   * Map for accumulating statistics of cross sections per process
   * number.
   */
  typedef map<int,XSecStat> StatMap;

  /**
   * Map of XComb objects describing the incoming partons indexed by
   * the corresponding PartonBin pair.
   */
  typedef map<tcPBPair,XCombPtr> XCombMap;

  /**
   * A vector of pointers to ReweightBase objects.
   */
  typedef vector<ReweightPtr> ReweightVector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor. If the optional argument is true, the reader
   * is assumed to be able to produce events on demand for a given
   * process.
   */
  FxFxReader(bool active = false);

  /**
   * Copy-constructor.
   */
  FxFxReader(const FxFxReader &);

  /**
   * Destructor.
   */
  virtual ~FxFxReader();
  //@}

public:

  /** @name Main virtual fuctions to be overridden in
   *  sub-classes. They are named as if we are reading from event
   *  files, but could equally well implement the actual generation of
   *  events. */
  //@{
  /**
   * Open a file or stream with events and read in the run information
   * into the heprup variable.
   */
  virtual void open() = 0;

  /**
   * Read the next event from the file or stream into the
   * corresponding protected variables. Return false if there is no
   * more events.
   */
  virtual bool doReadEvent() = 0;

  /**
   * Close the file or stream from which events have been read.
   */
  virtual void close() = 0;

  /**
   * return the weight names 
   */
  //  virtual vector<string> optWeightsNamesFunc();
  virtual vector<string> optWeightsNamesFunc() = 0;
  //virtual vector<string*> optWeightNamesFunc() = 0;
  vector<string> optionalWeightsNames;
  
  /** 
   * The ID (e.g. 100x, 2001) for the weight
   */ 

  // vector<string> optionalWeightsNames;

  
  //@}

  /** @name Other important function which may be overridden in
   *  sub-classes which wants to bypass the basic HEPRUP or HEPEUP
   *  variables or otherwise facilitate the conversion to ThePEG
   *  objects. */
  //@{
  /**
   * Initialize. This function is called by the FxFxEventHandler
   * to which this object is assigned.
   */
  virtual void initialize(FxFxEventHandler & eh);

  /**
   * Calls readEvent() or uncacheEvent() to read information into the
   * FxFx common block variables. This function is called by the
   * FxFxEventHandler if this reader has been selectod to
   * produce an event.
   *
   * @return the weight asociated with this event. If negative weights
   * are allowed it should be between -1 and 1, otherwise between 0
   * and 1. If outside these limits the previously estimated maximum
   * is violated. Note that the estimated maximum then should be
   * updated from the outside.
   */
  virtual double getEvent();

  /**
   * Calls doReadEvent() and performs pre-defined reweightings. A
   * sub-class overrides this function it must make sure that the
   * corresponding reweightings are done.
   */
  virtual bool readEvent();

  /**
   * Skip \a n events. Used by FxFxEventHandler to make sure
   * that a file is scanned an even number of times in case the events
   * are not ramdomly distributed in the file.
   */
  virtual void skip(long n);

  /**
   * Get an XComb object. Converts the information in the Les Houches
   * common block variables to an XComb object describing the sub
   * process. This is the way information is conveyed from the reader
   * to the controlling FxFxEventHandler.
   */
  tXCombPtr getXComb();

  /**
   * Get a SubProcess object corresponding to the information in the
   * Les Houches common block variables.
   */
  tSubProPtr getSubProcess();

  /**
   * Scan the file or stream to obtain information about cross section
   * weights and particles etc. This function should fill the
   * variables corresponding to the /HEPRUP/ common block. The
   * function returns the number of events scanned.
   */
  virtual long scan();

  /**
   * Take the information corresponding to the HEPRUP common block and
   * initialize the statistics for this reader.
   */
  virtual void initStat();

  /**
   * Reweights the current event using the reweights and preweights
   * vectors. It is the responsibility of the sub-class to call this
   * function after the HEPEUP information has been retrieved.
   */
  double reweight();

  /**
   * Converts the information in the Les Houches common block
   * variables into a Particle objects.
   */
  virtual void fillEvent();

  /**
   * Removes the particles created in the last generated event,
   * preparing to produce a new one.
   */
  void reset();

  /**
   * Possibility for subclasses to recover from non-conformant
   * settings of XMAXUP when an event file has been scanned with \a
   * neve events. Should set weightScale so that the average XMAXUP
   * times weightScale gives the cross section for a process. (This is
   * needed for MadEvent).
   */
  virtual void setWeightScale(long neve);

  //@}

  /** @name Access information about the current event. */
  //@{

  /**
   * Return the size of this event in bytes. To be used for the cache
   * file. \a npart is the number of particles. If \a npart is 0, the
   * number is taken from NUP.
   */
  static size_t eventSize(int N) {
    return (N + 1)*sizeof(int) +       // IDPRUP, ISTUP
      (7*N + 4)*sizeof(double) +       // XWGTUP, SCALUP, AQEDUP, AQCDUP, PUP,
                                       // VTIMUP, SPINUP
      N*sizeof(long) +                 // IDUP
      2*N*sizeof(pair<int,int>) +      // MOTHUP, ICOLUP
      sizeof(pair<double,double>) +    // XPDWUP.
      2*sizeof(double);                // lastweight and preweight
  }

  /**
   * The current event weight given by XWGTUP times possible
   * reweighting. Note that this is not necessarily the same as what
   * is returned by getEvent(), which is scaled with the maximum
   * weight.
   */
  double eventWeight() const { return hepeup.XWGTUP*lastweight; }

  /**
   * Return the optional named weights associated to the current event.
   */
  const map<string,double>& optionalEventWeights() const { return optionalWeights; }

  /**
   * Return the optional npLO and npNLO
   */
  const int& optionalEventnpLO() const { return optionalnpLO; }
  const int& optionalEventnpNLO() const { return optionalnpNLO; }
  
  /**
   * The pair of PartonBinInstance objects describing the current
   * incoming partons in the event.
   */
  const PBIPair & partonBinInstances() const { return thePartonBinInstances; }
  /**
   * Return the instances of the beam particles for the current event.
   */
  const PPair & beams() const { return theBeams; }
  /**
   * Return the instances of the incoming particles to the sub process
   * for the current event.
   */
  const PPair & incoming() const { return theIncoming; }
  /**
   * Return the instances of the outgoing particles from the sub process
   * for the current event.
   */
  const PVector & outgoing() const { return theOutgoing; }
  /**
   * Return the instances of the intermediate particles in the sub
   * process for the current event.
   */
  const PVector & intermediates() const { return theIntermediates; }
  /**
   * If this reader is to be used (possibly together with others) for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the highest multiplicity matrix element in
   * the group.
   */
  int maxMultCKKW() const { return theMaxMultCKKW; }
  /**
   * If this reader is to be used (possibly together with others) for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the lowest multiplicity matrix element in
   * the group.
   */
  int minMultCKKW() const { return theMinMultCKKW; }  //@}

  /** @name Other inlined access functions. */
  //@{
  /**
   * The number of events found in this reader. If less than zero the
   * number of events are unlimited.
   */
  long NEvents() const { return theNEvents; }

  /**
   * The number of events produced so far. Is reset to zero if an
   * event file is reopened.
   */
  long currentPosition() const { return position; }

  /**
   * The maximum number of events to scan to collect information about
   * processes and cross sections. If less than 0, all events will be
   * scanned.
   */
  long maxScan() const { return theMaxScan; }

  /**
   * Return true if this reader is active.
   */
  bool active() const { return isActive; }

  /**
   * True if negative weights may be produced.
   */
  bool negativeWeights() const { return heprup.IDWTUP < 0; }

  /**
   * The collected cross section statistics for this reader.
   */
  const XSecStat & xSecStats() const { return stats; }

  /**
   * Collected statistics about the individual processes.
   */
  const StatMap & processStats() const { return statmap; }

  /**
   * Select the current event. It will later be rejected with a
   * probability given by \a weight.
   */
  void select(double weight) {
    stats.select(weight);
    statmap[hepeup.IDPRUP].select(weight);
  }

  /**
   * Accept the current event assuming it was previously selcted.
   */
  void accept() {
    stats.accept();
    statmap[hepeup.IDPRUP].accept();
  }

  /**
   * Reject the current event assuming it was previously accepted.
   */
  void reject(double w) {
    stats.reject(w);
    statmap[hepeup.IDPRUP].reject(w);
  }

  /**
   * Increase the overestimated cross section for this reader.
   */
  virtual void increaseMaxXSec(CrossSection maxxsec);

  /**
   * The PartonExtractor object used to construct remnants.
   */
  tPExtrPtr partonExtractor() const { return thePartonExtractor; }

  /**
   * Return a possibly null pointer to a CascadeHandler to be used for
   * CKKW-reweighting.
   */
  tCascHdlPtr CKKWHandler() const { return theCKKW; }

  /**
   * The pairs of PartonBin objects describing the partons which can
   * be extracted by the PartonExtractor object.
   */
  const PartonPairVec & partonBins() const { return thePartonBins; }

  /**
   * The map of XComb objects indexed by the corresponding PartonBin
   * pair.
   */
  const XCombMap & xCombs() const { return theXCombs; }

  /**
   * The Cuts object to be used for this reader.
   */
  const Cuts & cuts() const { return *theCuts; }

  //@}

protected:

  /** @name Functions for manipulating cache files. */
  //@{

  /**
   * Name of file used to cache the events form the reader in a
   * fast-readable form. If empty, no cache file will be generated.
   */
  string cacheFileName() const { return theCacheFileName; }

  /**
   * Determines whether to apply cuts to events converting them to
   * ThePEG format.
   */
  bool cutEarly() const { return doCutEarly; }

  /**
   * File stream for the cache.
   */
  CFile cacheFile() const { return theCacheFile;}

  /**
   * Open the cache file for reading.
   */
  void openReadCacheFile();

  /**
   * Open the cache file for writing.
   */
  void openWriteCacheFile();

  /**
   * Close the cache file;
   */
  void closeCacheFile();

  /**
   * Write the current event to the cache file.
   */
  void cacheEvent() const;

  /**
   * Read an event from the cache file. Return false if something went wrong.
   */
  bool uncacheEvent();

  /**
   * Reopen a reader. If we have reached the end of an event file,
   * reopen it and issue a warning if we have used up a large fraction
   * of it.
   */
  void reopen();

  /**
   * Helper function to write a variable to a memory location
   */
  template <typename T>
  static char * mwrite(char * pos, const T & t, size_t n = 1) {
    std::memcpy(pos, &t, n*sizeof(T));
    return pos + n*sizeof(T);
  }

  /**
   * Helper function to read a variable from a memory location
   */
  template <typename T>
  static const char * mread(const char * pos, T & t, size_t n = 1) {
    std::memcpy(&t, pos, n*sizeof(T));
    return pos + n*sizeof(T);
  }

  //@}

  /** @name Auxilliary virtual methods which may be verridden by sub-classes. */
  //@{
  /**
   * Check the existence of a pair of PartonBin objects corresponding
   * to the current event.
   *
   * @return false if no pair of suitable PartonBin objects was found.
   */
  virtual bool checkPartonBin();

  /**
   * Create instances of all particles in the event and store them
   * in particleIndex.
   */
  virtual void createParticles();

  /**
   * Using the already created particles create a pair of
   * PartonBinInstance objects corresponding to the incoming
   * partons. Return the corresponding PartonBin objects.
   */
  virtual tcPBPair createPartonBinInstances();

  /**
   * Create instances of the incoming beams in the event and store
   * them in particleIndex. If no beam particles are included in the
   * event they are created from the run info.
   */
  virtual void createBeams();

  /**
   * Go through the mother indices and connect up the Particles.
   */
  virtual void connectMothers();
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /** @name Set functions for some variables not in the Les Houches accord. */
  //@{
  /**
   * The number of events in this reader. If less than zero the number
   * of events is unlimited.
   */
  void NEvents(long x) { theNEvents = x; }

  /**
   * The map of XComb objects indexed by the corresponding PartonBin
   * pair.
   */
  XCombMap & xCombs() { return theXCombs; }  
  //@}

  /** @name Standard (and non-standard) Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish() {
    close();
    HandlerBase::dofinish();
  }

  /**
   * Return true if this object needs to be initialized before all
   * other objects because it needs to extract PDFs from the event file.
   */
  virtual bool preInitialize() const;

  /**
   * Called from doinit() to extract PDFs from the event file and add
   * the corresponding objects to the current EventGenerator.
   */
  virtual void initPDFs();
  //@}

protected:

  /**
   * The HEPRUP common block.
   */
  HEPRUP heprup;

  /**
   * The HEPEUP common block.
   */
  HEPEUP hepeup;

  /**
   * The ParticleData objects corresponding to the incoming particles.
   */
  tcPDPair inData;

  /**
   * The PDFBase objects which has been used for the beam particle
   * when generating the events being read. Specified in the interface
   * or derived from PDFGUP and PDFSUP.
   */
  pair<PDFPtr,PDFPtr> inPDF;

  /**
   * The PDFBase object to be used in the subsequent generation.
   */
  pair<cPDFPtr,cPDFPtr> outPDF;

  /**
   * The PartonExtractor object used to construct remnants.
   */
  PExtrPtr thePartonExtractor;

  /**
   * A pointer to a CascadeHandler to be used for CKKW-reweighting.
   */
  tCascHdlPtr theCKKW;

  /**
   * The pairs of PartonBin objects describing the partons which can
   * be extracted by the PartonExtractor object.
   */
  PartonPairVec thePartonBins;

  /**
   * The map of XComb objects indexed by the corresponding PartonBin
   * pair.
   */
  XCombMap theXCombs;

  /**
   * The Cuts object to be used for this reader.
   */
  CutsPtr theCuts;

  /**
   * The number of events in this reader. If less than zero the number
   * of events is unlimited.
   */
  long theNEvents;

  /**
   * The number of events produced by this reader so far. Is reset
   * every time an event file is reopened.
   */
  long position;

  /**
   * The number of times this reader has been reopened.
   */
  int reopened;

  /**
   * The maximum number of events to scan to collect information about
   * processes and cross sections. If less than 0, all events will be
   * scanned.
   */
  long theMaxScan;

  /**
   * Flag to tell whether we are in the process of scanning.
   */
  bool scanning;

  /**
   * True if this is an active reader.
   */
  bool isActive;

  /**
   * Name of file used to cache the events form the reader in a
   * fast-readable form. If empty, no cache file will be generated.
   */
  string theCacheFileName;

  /**
   * Determines whether to apply cuts to events before converting them
   * to ThePEG format.
   */
  bool doCutEarly;

  /**
   * Collect statistics for this reader.
   */
  XSecStat stats;

  /**
   * Collect statistics for each individual process.
   */
  StatMap statmap;

  /**
   * The pair of PartonBinInstance objects describing the current
   * incoming partons in the event.
   */
  PBIPair thePartonBinInstances;

  /**
   * Association between ColourLines and colour indices in the current
   * translation.
   */
  ObjectIndexer<long,ColourLine> colourIndex;

  /**
   * Association between Particles and indices in the current
   * translation.
   */
  ObjectIndexer<long,Particle> particleIndex;

  /**
   * The instances of the beam particles for the current event.
   */
  PPair theBeams;

  /**
   * The instances of the incoming particles to the sub process for
   * the current event.
   */
  PPair theIncoming;

  /**
   * The instances of the outgoing particles from the sub process for
   * the current event.
   */
  PVector theOutgoing;

  /**
   * The instances of the intermediate particles in the sub process for
   * the current event.
   */
  PVector theIntermediates;

  /**
   * File stream for the cache.
   */
  CFile theCacheFile;

  /**
   * The reweight objects modifying the weights of this reader.
   */
  ReweightVector reweights;

  /**
   * The preweight objects modifying the weights of this reader.
   */
  ReweightVector preweights;

  /**
   * The factor with which this reader was last pre-weighted.
   */
  double preweight;

  /**
   * Should the event be reweighted by PDFs used by the PartonExtractor?
   */
  bool reweightPDF;

  /**
   * Should PDFBase objects be constructed from the information in the
   * event file in the initialization?
   */
  bool doInitPDFs;

  /**
   * If this reader is to be used (possibly together with others) for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the highest multiplicity matrix element in
   * the group.
   */
  int theMaxMultCKKW;

  /**
   * If this reader is to be used (possibly together with others) for
   * CKKW reweighting and veto, this should give the multiplicity of
   * outgoing particles in the lowest multiplicity matrix element in
   * the group.
   */
  int theMinMultCKKW;

  /**
   * The weight multiplying the last read event due to PDF
   * reweighting, CKKW reweighting or assigned reweight and preweight
   * objects.
   */
  double lastweight;

  /**
   * The optional weights associated to the last read events.
   */
  map<string,double> optionalWeights;

  /**
   * npLO for FxFx merging
   */
  int optionalnpLO;

 /**
   * npNLO for FxFx merging
   */
  int optionalnpNLO;

  /**
   * If the maximum cross section of this reader has been increased
   * with increaseMaxXSec(), this is the total factor with which it
   * has been increased.
   */
  double maxFactor;

  /**
   * The (reweighted) XWGTUP value should be scaled with this cross
   * section when compared to the overestimated cross section.
   */
  CrossSection weightScale;

  /**
   * Individual scales for different sub-processes if reweighted.
   */
  vector<double> xSecWeights;

  /**
   * Individual maximum weights for individual (possibly reweighted)
   * processes.
   */
  map<int,double> maxWeights;

  /**
   * Is set to true when getEvent() is called from skip(int).
   */
  bool skipping;

  /**
   *  Option for the treatment of the momenta supplied
   */
  unsigned int theMomentumTreatment;

  /**
   * Set to true if warnings about possible weight incompatibilities
   * should be issued.
   */
  bool useWeightWarnings;

  /**
   *  Option to allow reopening of the file
   */
  bool theReOpenAllowed;

  /**
   *  Use the spin information
   */
  bool theIncludeSpin;

private:

  /** Access function for the interface. */
  void setBeamA(long id);
  /** Access function for the interface. */
  long getBeamA() const;
  /** Access function for the interface. */
  void setBeamB(long id);
  /** Access function for the interface. */
  long getBeamB() const;
  /** Access function for the interface. */
  void setEBeamA(Energy e);
  /** Access function for the interface. */
  Energy getEBeamA() const;
  /** Access function for the interface. */
  void setEBeamB(Energy e);
  /** Access function for the interface. */
  Energy getEBeamB() const;
  /** Access function for the interface. */
  void setPDFA(PDFPtr);
  /** Access function for the interface. */
  PDFPtr getPDFA() const;
  /** Access function for the interface. */
  void setPDFB(PDFPtr);
  /** Access function for the interface. */
  PDFPtr getPDFB() const;

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<FxFxReader> initFxFxReader;

  /**
   * Private and non-existent assignment operator.
   */
  FxFxReader & operator=(const FxFxReader &);

public:

  /** @cond EXCEPTIONCLASSES */
  /** Exception class used by FxFxReader in case inconsistencies
   *  are encountered. */
  class FxFxInconsistencyError: public Exception {};
  
  /** Exception class used by FxFxReader in case more events
      than available are requested. */
  class FxFxReopenWarning: public Exception {};

  /** Exception class used by FxFxReader in case reopening an
      event file fails. */
  class FxFxReopenError: public Exception {};

  /** Exception class used by FxFxReader in case there is
      information missing in the initialization phase. */
  class FxFxInitError: public InitException {};
  /** @endcond */

};

/// Stream output for HEPEUP
ostream & operator<<(ostream & os, const HEPEUP & h);

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of FxFxReader.
 */
template <>
struct BaseClassTrait<FxFxReader,1>: public ClassTraitsType {
  /** Typedef of the base class of FxFxReader. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * FxFxReader class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<FxFxReader>
  : public ClassTraitsBase<FxFxReader> {
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::FxFxReader"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the FxFxReader class and every other class it uses
   * (except the base class).
   */
  static string library() { return "FxFx.so"; }

};

/** @endcond */

}

#endif /* THEPEG_FxFxReader_H */
