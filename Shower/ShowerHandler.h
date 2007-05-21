// -*- C++ -*-
#ifndef HERWIG_ShowerHandler_H
#define HERWIG_ShowerHandler_H
//
// This is the declaration of the ShowerHandler class.
//

#include "ThePEG/Handlers/CascadeHandler.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig++/Shower/Base/Evolver.fh"
#include "Herwig++/Shower/Base/ShowerParticle.fh"
#include "Herwig++/Shower/Base/ShowerTree.fh"
#include "ShowerHandler.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is the main driver of the shower: it is responsible for 
 *  the proper handling of all other specific collaborating classes
 *  and for the storing of the produced particles in the event record.
 * 
 *  @see CascadeHandler
 */
class ShowerHandler: public CascadeHandler {

public:

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
   * The main method which manages the all showering.
   */
  virtual void cascade();

  /**
   * It returns true if the particle with the specified id
   * is in the list of those that should be decayed during the showering
   * showering.
   */
  inline bool decayInShower(const long id) const;

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

protected:

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
  inline virtual void doinit() throw(InitException);
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
   *  Pointer to the evolver
   */
  EvolverPtr _evolver;

  /**
   *  Maximum number of attempts for the
   *   main showering loop
   */
  unsigned int _maxtry;

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
  static string className() { return "Herwig++::ShowerHandler"; }
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

#include "ShowerHandler.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerHandler.tcc"
#endif

#endif /* HERWIG_ShowerHandler_H */
