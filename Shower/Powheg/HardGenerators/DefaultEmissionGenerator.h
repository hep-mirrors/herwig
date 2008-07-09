// -*- C++ -*-
#ifndef HERWIG_DefaultEmissionGenerator_H
#define HERWIG_DefaultEmissionGenerator_H
//
// This is the declaration of the DefaultEmissionGenerator class.
//

#include "Herwig++/Shower/Powheg/HardestEmissionGenerator.h"
#include "pTSudakov.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Utilities/Histogram.h"
#include "DefaultEmissionGenerator.fh"

namespace Herwig {

using namespace ThePEG;

/**
 *  typedef to pair the SudakovFormFactor and the particles in a branching
 */
typedef pair<pTSudakovPtr,IdList> pTBranchingElement;

/**
 *  typedef to pair the PDG code of the particle and the BranchingElement
 */
typedef multimap<long,pTBranchingElement> pTBranchingList;

/**
 *  typedef to create a structure which can be inserted into a BranchingList
 */
typedef pair<long, pTBranchingElement> pTBranchingInsert; 

/**
 * Here is the documentation of the DefaultEmissionGenerator class.
 *
 * @see \ref DefaultEmissionGeneratorInterfaces "The interfaces"
 * defined for DefaultEmissionGenerator.
 */
class DefaultEmissionGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  inline DefaultEmissionGenerator();

  /**
   *  Implementation of virtual functions from the base class
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual NasonTreePtr generateHardest(ShowerTreePtr);

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr);
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

protected:

  /**
   *  Member to generate the hardest emission for hard processes and decays
   *  respectively.
   */
  //@{
  /**
   *  Member to generate the hardest emission for a hard process
   */
  NasonTreePtr generateHard(ShowerTreePtr);

  /**
   *  Member to generate the hardest emission for a decay process
   */
  NasonTreePtr generateDecay(ShowerTreePtr);
  //@}

  /**
   *  Members to generate the emission
   */
  //@{
  /**
   * Choose a new forward branching for a time-like particle
   * The method returns:
   * - a pointer to a ShowerKinematics object, which 
   *     contains the information about the new scale and all other
   *     kinematics variables that need to be generated simultaneously;
   * - a pointer to the SudakovFormFactor object associated 
   *     with the chosen emission.
   * - The PDG codes of the particles in the branching,
   * as a Branching struct.
   *
   * In the case no branching has been generated, both the returned 
   * pointers are null ( ShoKinPtr() , tSudakovFFPtr() ).
   *
   * @param particle The particle to be evolved
   * @return The Branching struct for the branching
   */
  Branching chooseForwardBranching(ShowerParticle & particle) const; 
  //@}

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

  /**
   *  Construct a new \f$p_T\f$ ordered Sudakov
   */
  pTSudakovPtr constructSudakov(tSudakovPtr);

  /**
   *  Reconstruct a final emission
   */
  bool reconstructFinal(tShowerParticlePtr emitter,tShowerParticlePtr spectator,
			ShoKinPtr, IdList) const;

  /**
   * Compute the boost to get from the the old momentum to the new 
   */
  inline LorentzRotation solveBoost(const Lorentz5Momentum & newq, 
				    const Lorentz5Momentum & oldq) const;

protected:

  /** @name Standard Interfaced functions. */
  //@{
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DefaultEmissionGenerator> initDefaultEmissionGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DefaultEmissionGenerator & operator=(const DefaultEmissionGenerator &);

private:

  /**
   *  List of the branchings and the appropriate Sudakovs for forward branchings
   */
  pTBranchingList _fbranchings;

  /**  
   * Lists of the branchings and the appropriate Sudakovs for backward branchings.
   */
  pTBranchingList _bbranchings;

  /**
   *  Histogram for the thrust
   */
  HistogramPtr _thrust[2];

  /**
   *  Histograms for the pt
   */
  HistogramPtr _pthist[2];

  /**
   * storage for the scatter plot
   */
  vector<double> _xq,_xqbar;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DefaultEmissionGenerator. */
template <>
struct BaseClassTrait<Herwig::DefaultEmissionGenerator,1> {
  /** Typedef of the first base class of DefaultEmissionGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DefaultEmissionGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DefaultEmissionGenerator>
  : public ClassTraitsBase<Herwig::DefaultEmissionGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DefaultEmissionGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DefaultEmissionGenerator is implemented. It may also include several, space-separated,
   * libraries if the class DefaultEmissionGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#include "DefaultEmissionGenerator.icc"

#endif /* HERWIG_DefaultEmissionGenerator_H */
