// -*- C++ -*-
//
// ShowerHandler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
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
#include "Herwig++/Shower/UEBase.h" 
#include "Herwig++/Shower/Base/Evolver.fh"
#include "Herwig++/Shower/Base/ShowerParticle.fh"
#include "Herwig++/Shower/Base/ShowerTree.fh"
#include "Herwig++/Shower/Base/HardTree.fh"
#include "Herwig++/PDF/HwRemDecayer.fh"
#include "ThePEG/EventRecord/RemnantParticle.fh"
#include "ShowerHandler.fh"

namespace Herwig {

/**
 *  Typedef for the ShowerTree for the decays
 */
typedef multimap<Energy,ShowerTreePtr,std::greater<Energy> > ShowerDecayMap;

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
   * The main method which manages the multiple interactions and starts
   * the shower by calling cascade(sub, lastXC).
   */
  virtual void cascade();

  /**
   * Hook to allow vetoing of event after showering hard sub-process
   * as in e.g. MLM merging.
   */
  virtual bool showerHardProcessVeto() { return false; };

  /**
   * It returns true if the particle with the specified id
   * is in the list of those that should be decayed during the showering
   * showering.
   */
  bool decaysInShower(const long id) const {
    return ( particlesDecayInShower_.find( abs(id) ) != 
	     particlesDecayInShower_.end() ); 
  }

public:

  /**@name Methods related to PDF freezing */
  //@{
  /**
   * Get the PDF freezing scale
   */
  Energy pdfFreezingScale() const {return pdfFreezingScale_;}
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

public:

  /** @name Functions to access information. */
  //@{

  /**
   * Return true if currently the primary subprocess is showered.
   */
  bool firstInteraction() const {
    return ( subProcess_ == 
	     eventHandler()->currentCollision()->primarySubProcess() );
  }

  /**
   * Return the currently used SubProcess.
   */
  tSubProPtr currentSubProcess() const {
    assert(subProcess_);
    return subProcess_;
  }

  /**
   * Return true if multiple parton interactions are switched on 
   * and can be used for this beam setup.
   */
  bool isMPIOn() const {
    return MPIHandler_ && MPIHandler_->beamOK();
  }

  /**
   * Return the remnant decayer.
   */
  tHwRemDecPtr remnantDecayer() const { return remDec_; }
  //@}

  /**
   *  Access to the Evolver
   */
  tEvolverPtr evolver() const {return evolver_;}

  /**
   *  Generate hard emissions for CKKW etc
   */
  virtual HardTreePtr generateCKKW(ShowerTreePtr tree) const;

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
   * Prepare to shower the given subprocess
   */
  void prepareCascade(tSubProPtr sub);

  /**
   * The main method which manages the showering of a subprocess.
   */
  virtual tPPair cascade(tSubProPtr sub, XCPtr xcomb);

  /**
   * Return the maximum number of attempts for showering
   * a given subprocess.
   */
  unsigned int maxtry() const { return maxtry_; }

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
  PPtr findParent(PPtr parent, bool & isHard, set<PPtr> outgoing,
		  bool checkStep, tSubProPtr process=SubProPtr()) const;

  /**
   * Find the parton extracted from the incoming particle after ISR
   */
  PPtr findFirstParton(tPPtr seed) const;

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
   *  Reset the PDF's after the hard collision has been showered
   */
  void setMPIPDFs();

  /**
   *  Test for decay products
   */
  bool decayProduct(tPPtr,tSubProPtr process=SubProPtr()) const;

  /**
   *  Boost all the particles in the collision so that the collision always occurs
   * in the rest frame with the incoming particles along the z axis
   */
  void boostCollision(bool boost);

  /**
   *  Is a beam particle where hadronic structure is resolved
   */
  bool isResolvedHadron(tPPtr);

protected:

  /** @name Standard Interfaced functions. */
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
   * Access function for the MPIHandler, it should only be called after
   * checking with isMPIOn.
   */
  tUEBasePtr getMPIHandler() const {  
    assert(MPIHandler_);
    return MPIHandler_;    
  }

private:

  /**
   * a MPIHandler to administer the creation of several (semihard) 
   * partonic interactions.
   */
  UEBasePtr MPIHandler_;

  /**
   *  Pointer to the evolver
   */
  EvolverPtr evolver_;

  /**
   *  Pointer to the HwRemDecayer
   */
  HwRemDecPtr remDec_;

  /**
   * The PDF for beam particle A. Overrides the particle's own PDF setting.
   */
  PDFPtr PDFA_;

  /**
   * The PDF for beam particle B. Overrides the particle's own PDF setting.
   */
  PDFPtr PDFB_;

  /**
   * The PDF freezing scale
   */
  Energy pdfFreezingScale_;

  /**
   *  Maximum number of attempts for the
   *   main showering loop
   */
  unsigned int maxtry_;

  /**
   *  Maximum number of attempts for the regeneration of an additional
   *  scattering, before the number of scatters is reduced.
   */
  unsigned int maxtryMPI_;

  /**
   *  Maximum number of attempts for the regeneration of an additional
   *  hard scattering, before this event is vetoed.
   */
  unsigned int maxtryDP_;

  /**
   *  PDG codes of the particles which decay during showering
   *  this is fast storage for use during running
   */
  set<long> particlesDecayInShower_;

  /**
   *  PDG codes of the particles which decay during showering
   *  this is a vector that is interfaced so they can be changed
   */
  vector<long> inputparticlesDecayInShower_;

  /**
   *  The ShowerTree for the hard process
   */
  ShowerTreePtr hard_;

  /**
   *  The incoming beam particles for the current collision
   */
  tPPair incoming_;

  /**
   *  The ShowerTree for the decays
   */
  ShowerDecayMap decay_;

  /**
   *  The ShowerTrees for which the initial shower 
   */
  vector<ShowerTreePtr> done_;

  /**
   *  Const pointer to the current step
   */
  tcStepPtr current_;

  /**
   *  Const pointer to the currently handeled ThePEG::SubProcess
   */
  tSubProPtr subProcess_;

  /**
   *  pointer to "this", the current ShowerHandler.
   */
  static ShowerHandler * currentHandler_;

  /**
   *  Boost to get back to the lab
   */
  LorentzRotation boost_;

  /**
   * The MPI PDF's to be used for secondary scatters.
   */
  pair <PDFPtr, PDFPtr> mpipdfs_;

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
    const int tries;

    /** constructor */
    ShowerTriesVeto(int t) : tries(t) {}
  };

  /**
   *  pointer to "this", the current ShowerHandler.
   */
  static const ShowerHandler * currentHandler() {
    assert(currentHandler_);
    return currentHandler_;
  }

protected:

  /**
   *  Set the current handler
   */
  void setCurrentHandler() {
    currentHandler_ = this;
  }

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
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ShowerHandler_H */
