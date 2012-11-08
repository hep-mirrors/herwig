// -*- C++ -*-
//
// DipoleShowerHandler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleShowerHandler_H
#define HERWIG_DipoleShowerHandler_H
//
// This is the declaration of the DipoleShowerHandler class.
//

#include "Herwig++/Shower/ShowerHandler.h"

#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingReweight.h"
#include "Herwig++/DipoleShower/Kernels/DipoleSplittingKernel.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingGenerator.h"
#include "Herwig++/DipoleShower/Base/DipoleEventRecord.h"
#include "Herwig++/DipoleShower/Base/DipoleEvolutionOrdering.h"
#include "Herwig++/DipoleShower/Utility/ConstituentReshuffler.h"
#include "Herwig++/DipoleShower/Utility/IntrinsicPtGenerator.h"

namespace Herwig {

using namespace ThePEG;

/** 
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief The DipoleShowerHandler class manages the showering using
 * the dipole shower algorithm.
 *
 * @see \ref DipoleShowerHandlerInterfaces "The interfaces"
 * defined for DipoleShowerHandler.
 */
class DipoleShowerHandler: public ShowerHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleShowerHandler();

  /**
   * The destructor.
   */
  virtual ~DipoleShowerHandler();
  //@}

public:

  /**
   * Indicate a problem in the shower.
   */
  struct RedoShower {};

  /**
   * Insert an additional splitting kernel.
   */
  void addSplitting(Ptr<DipoleSplittingKernel>::ptr sp) {
    kernels.push_back(sp);
  }

  /**
   * Reset the alpha_s for all splitting kernels.
   */
  void resetAlphaS(Ptr<AlphaSBase>::tptr);

  /**
   * Reset the splitting reweight for all splitting kernels.
   */
  void resetReweight(Ptr<DipoleSplittingReweight>::tptr);

protected:

  typedef multimap<DipoleIndex,Ptr<DipoleSplittingGenerator>::ptr> GeneratorMap;

  /**
   * The main method which manages the showering of a subprocess.
   */
  virtual tPPair cascade(tSubProPtr sub, XCPtr xcomb);

  /**
   * Build splitting generators for the given
   * dipole index.
   */
  void getGenerators(const DipoleIndex&,
		     Ptr<DipoleSplittingReweight>::tptr rw =
		     Ptr<DipoleSplittingReweight>::tptr());

  /**
   * Setup the hard scales.
   */
  void hardScales();

  /**
   * Return the evolution ordering
   */
  Ptr<DipoleEvolutionOrdering>::tptr evolutionOrdering() const { return theEvolutionOrdering; }

  /**
   * Reshuffle to constituent mass shells
   */
  void constituentReshuffle();

  /**
   * Access the generator map
   */
  GeneratorMap& generators() { return theGenerators; }

  /**
   * Access the event record
   */
  DipoleEventRecord& eventRecord() { return theEventRecord; }

  /**
   * Return the event record
   */
  const DipoleEventRecord& eventRecord() const { return theEventRecord; }

  /**
   * Return the splitting kernels.
   */
  const vector<Ptr<DipoleSplittingKernel>::ptr>& splittingKernels() const {
    return kernels;
  }
  
  /**
   * Realign the event such as to have the incoming partons along thre
   * beam axes.
   */
  bool realign();

private:

  /**
   * Perform the cascade.
   */
  void doCascade(unsigned int& emDone);

  /**
   * Get the winning splitting for the
   * given dipole and configuration.
   */
  Energy getWinner(DipoleSplittingInfo& winner,
		   const Dipole& dip,
		   pair<bool,bool> conf);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


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
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}


private:

  /**
   * The splitting kernels to be used.
   */
  vector<Ptr<DipoleSplittingKernel>::ptr> kernels;

  /**
   * The evolution ordering considered
   */
  Ptr<DipoleEvolutionOrdering>::ptr theEvolutionOrdering;

  /**
   * The ConstituentReshuffler to be used
   */
  Ptr<ConstituentReshuffler>::ptr constituentReshuffler;

  /**
   * The intrinsic pt generator to be used.
   */
  Ptr<IntrinsicPtGenerator>::ptr intrinsicPtGenerator;

  /**
   * A global alpha_s to be used for all splitting kernels.
   */
  Ptr<AlphaSBase>::ptr theGlobalAlphaS;

  /**
   * Apply chain ordering to events from matrix
   * element corrections.
   */
  bool chainOrderVetoScales;

  /**
   * Limit the number of emissions.
   * Limit applied if > 0.
   */
  unsigned int nEmissions;

  /**
   * Discard events which did not radiate.
   */
  bool discardNoEmissions;

  /**
   * Perform the first MC@NLO emission only.
   */
  bool firstMCatNLOEmission;

  /**
   * Switch on or off final state radiation.
   */
  bool doFSR;

  /**
   * Switch on or off initial state radiation.
   */
  bool doISR;

  /**
   * The realignment scheme
   */
  int realignmentScheme;

  /**
   * True, if first emission should use the available phase space
   */
  bool hardFirstEmission;

private:

  /**
   * The verbosity level.
   * 0 - print no info
   * 1 - print diagnostic information on setting up
   *     splitting generators etc.
   * 2 - print detailed event information for up to
   *     printEvent events.
   * 3 - print dipole chains after each splitting.
   */
  int verbosity;

  /**
   * See verbosity.
   */
  int printEvent;

private:

  /**
   * The splitting generators indexed by the dipole
   * indices they can work on.
   */
  GeneratorMap theGenerators;

  /**
   * The evnt record used.
   */
  DipoleEventRecord theEventRecord;

  /**
   * The number of shoer tries so far.
   */
  unsigned int nTries;

  /**
   * Whether or not we did radiate anything
   */
  bool didRadiate;

  /**
   * Whether or not we did realign the event
   */
  bool didRealign;

private:

  /**
   * The factorization scale factor.
   */
  double theFactorizationScaleFactor;

  /**
   * The renormalization scale factor.
   */
  double theRenormalizationScaleFactor;

  /**
   * The scale factor for the hard scale
   */
  double theHardScaleFactor;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipoleShowerHandler> initDipoleShowerHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleShowerHandler & operator=(const DipoleShowerHandler &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleShowerHandler. */
template <>
struct BaseClassTrait<Herwig::DipoleShowerHandler,1> {
  /** Typedef of the first base class of DipoleShowerHandler. */
  typedef Herwig::ShowerHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleShowerHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleShowerHandler>
  : public ClassTraitsBase<Herwig::DipoleShowerHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleShowerHandler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleShowerHandler is implemented. It may also include several, space-separated,
   * libraries if the class DipoleShowerHandler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleShowerHandler_H */
