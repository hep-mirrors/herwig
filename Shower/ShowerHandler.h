// -*- C++ -*-
//
// ShowerHandler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerHandler_H
#define HERWIG_ShowerHandler_H
//
// This is the declaration of the ShowerHandler class.
//

#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/PDT/RemnantDecayer.h"
#include "ThePEG/EventRecord/RemnantParticle.h"

#include "Herwig++/UnderlyingEvent/MPIHandler.h" 
#include "Herwig++/Shower/Base/Evolver.fh"
#include "Herwig++/Shower/Base/ShowerParticle.fh"
#include "Herwig++/Shower/Base/ShowerTree.fh"
#include "Herwig++/PDF/HwRemDecayer.h"
#include "Herwig++/Shower/CKKW/Clustering/CascadeReconstructor.fh"
#include "Herwig++/Shower/CKKW/Reweighting/Reweighter.fh"
#include "ShowerHandler.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is the main driver of the shower: it is responsible for 
 *  the proper handling of all other specific collaborating classes
 *  and for the storing of the produced particles in the event record.
 * 
 *  @see \ref ShowerHandlerInterfaces "The interfaces"
 *
 *  @see ThePEG::CascadeHandler
 *  @see MPIHandler
 *  @see HwRemDecayer
 */
class ShowerHandler: public CascadeHandler {

public:
  
  /** Typedef for a pair of ThePEG::RemnantParticle pointers. */
  typedef pair<tRemPPtr, tRemPPtr> RemPair;

  /**
   * The default constructor.
   */
  ShowerHandler();

  /**
   *  Destructor
   */
  virtual ~ShowerHandler();

public:

  /**
   * The main method which manages the multiple interactions and starts the shower by calling
   * cascade(sub, lastXC).
   */
  virtual void cascade();

  /**
   * It returns true if the particle with the specified id
   * is in the list of those that should be decayed during the showering
   * showering.
   */
  inline bool decayInShower(const long id) const;

public:

  /**@name Methods related to ME/PS merging */
  //@{

  /**
   * Perform CKKW reweighting
   */
  virtual double reweightCKKW(int minMult, int maxMult);

  /**
   * Return the cascade reconstructor
   */
  inline tCascadeReconstructorPtr cascadeReconstructor () const;

  /**
   * Return the reweighter
   */
  inline tReweighterPtr reweighter () const;

  //@}

public:

  /**@name Methods related to PDF freezing */
  //@{

  /**
   * Set the PDF freezing scale
   */
  inline void pdfFreezingScale (Energy scale) {
    _pdfFreezingScale = scale;
  }

  /**
   * Get the PDF freezing scale
   */
  inline Energy pdfFreezingScale () const {
    return _pdfFreezingScale;
  }

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


  /** @name Functions to access information. */
  //@{

  /**
   * Return true if currently the primary subprocess is showered.
   */
  inline bool FirstInt() const;

  /**
   * Return the currently used SubProcess.
   */
  inline tSubProPtr currentSubProcess() const;

  /**
   * Return true if hard multiple parton interactions are ordered 
   * according to their scale.
   */
  inline bool IsOrdered() const;

  /**
   * Return true if multiple parton interactions are switched on 
   * and can be used for this beam setup.
   */
  inline bool IsMPIOn() const;
  //@}


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

protected:

  /**
   * The main method which manages the showering of a subprocess.
   */
  tPPair cascade(tSubProPtr sub);

  /**
   * At the end of the Showering, transform ShowerParticle objects
   * into ThePEG particles and fill the event record with them.
   * Notice that the parent/child relationships and the 
   * transformation from ShowerColourLine objects into ThePEG
   * ColourLine ones must be properly handled.
   */
  void fillEventRecord();

  /**
   * Identify the particles in the hard process and decayed particles
   * which need to be showered
   */
  void findShoweringParticles();

  /**
   * Find the final unstable time-like parent of a particle
   * @param parent The ultimate parent for the decaying particle
   * @param isHard Whether nay particles in chain are from the hard process
   * @param outgoing The outgoing particles from the hard process
   */
  PPtr findParent(PPtr parent, bool & isHard, set<PPtr> outgoing) const;

  /**
   * Find the parton extracted from the incoming particle after ISR
   */
  PPtr findFirstParton(tPPtr seed, tPPair incoming) const;

  /**
   * Fix Remnant connections after ISR
   */
  tPPair remakeRemnant(tPPair oldp); 

  /**
   * Get the remnants from the ThePEG::PartonBinInstance es and 
   * do some checks.
   */
  RemPair getRemnants(PBIPair incbins);

  /**
   *  Make the remnant after the shower
   */
  void makeRemnants();

  /**
   *  Test for decay products
   */
  bool decayProduct(tPPtr) const;

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Called at the end of the run phase.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ShowerHandler> initShowerHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerHandler & operator=(const ShowerHandler &);

private:
  /**
   * Access function for the MPIHandler. If it is a zero pointer
   * an exception is thrown, because it should only be called after
   * checking with IsMPIOn.
   */
  inline tMPIHPtr getMPIHandler() const;

  /**
   * Switch for Multi Parton Interactions to be ordered
   */
  bool theOrderSecondaries;

  /**
   * Switch for Multi Parton Interactions. Not used any more.
   */
  bool theMPIOnOff;

  /**
   * a MPIHandler to administer the creation of several (semihard) 
   * partonic interactions.
   */
  MPIHPtr theMPIHandler;

private:

  /**
   *  Pointer to the evolver
   */
  EvolverPtr _evolver;

  /**
   *  Pointer to the HwRemDecayer
   */
  HwRemDecPtr theRemDec;

  /**
   * The PDF freezing scale
   */
  Energy _pdfFreezingScale;

  /**
   *  Maximum number of attempts for the
   *   main showering loop
   */
  unsigned int _maxtry;

  /**
   *  Maximum number of attempts for the regeneration of an additional
   *  scattering, before the number of scatters is reduced.
   */
  unsigned int _maxtryMPI;

  /**
   *  PDG codes of the particles which decay during showering
   *  this is fast storage for use during running
   */
  set<long> _particlesDecayInShower;

  /**
   *  PDG codes of the particles which decay during showering
   *  this is a vector that is interfaced so they can be changed
   */
  vector<long> _inputparticlesDecayInShower;

  /**
   *  The ShowerTree for the hard process
   */
  ShowerTreePtr _hard;

  /**
   *  The ShowerTree for the decays
   */
  multimap<Energy,ShowerTreePtr> _decay;

  /**
   *  The ShowerTrees for which the initial shower 
   */
  vector<ShowerTreePtr> _done;

  /**
   *  Const pointer to the current step
   */
  tcStepPtr _current;

  /**
   *  Const pointer to the currently handeled ThePEG::SubProcess
   */
  tSubProPtr theSubProcess;

  /**
   *  pointer to "this", the current ShowerHandler.
   */
  static ShowerHandler * theHandler;

public:

  /** 
   * struct that is used to catch exceptions which are thrown
   * due to energy conservation issues of additional scatters
   */
  struct ExtraScatterVeto {};

  /** 
   * struct that is used to catch exceptions which are thrown
   * due to fact that the Shower has been invoked more than
   * a defined threshold on a certain configuration
   */
  struct ShowerTriesVeto {
    /** variable to store the number of attempts */
    int theTries;

    /** constructor */
    ShowerTriesVeto(int tries){theTries = tries;}
  };

  /**
   *  pointer to "this", the current ShowerHandler.
   */
  static inline const ShowerHandler * currentHandler();

private:

  /**
   * Wether or not to use CKKW
   */
  bool _useCKKW;

  /**
   * The cascade reconstructor used for ME/PS merging
   */
  CascadeReconstructorPtr _reconstructor;

  /**
   * The reweighter used for ME/PS merging
   */
  ReweighterPtr _reweighter;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerHandler. */
template <>
struct BaseClassTrait<Herwig::ShowerHandler,1> {
  /** Typedef of the first base class of ShowerHandler. */
  typedef CascadeHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerHandler>
  : public ClassTraitsBase<Herwig::ShowerHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ShowerHandler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerHandler is implemented. It may also include several, space-separated,
   * libraries if the class ShowerHandler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "ShowerHandler.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerHandler.tcc"
#endif

#endif /* HERWIG_ShowerHandler_H */
