// -*- C++ -*-
#ifndef HERWIG_Evolver_H
#define HERWIG_Evolver_H
//
// This is the declaration of the Evolver class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/MECorrections/MECorrectionBase.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Shower/Kinematics/KinematicsReconstructor.h"
#include "PartnerFinder.h"
#include "ShowerTree.fh"
#include "ShowerProgenitor.fh"
#include "ShowerHandler.fh"
#include "Evolver.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * Here is the documentation of the Evolver class.
 *
 * @see \ref EvolverInterfaces "The interfaces"
 * defined for Evolver.
 */
class Evolver: public Interfaced {

/**
 *  The ShowerHandler is a friend to set some parameters at initialisation
 */
friend class ShowerHandler;

public:

  /**
   *  Default Constructor
   */
  inline Evolver();
  /**
   *  Is there any showering switched on
   */
  inline bool showeringON() const;

  /**
   * It returns true/false if the initial-state radiation is on/off.
   */
  inline bool isISRadiationON() const;  

  /**
   * It returns true/false if the final-state radiation is on/off.
   */
  inline bool isFSRadiationON() const;  

  /**
   * Perform the shower of the hard process
   */
  virtual void showerHardProcess(ShowerTreePtr);


  /**
   * Perform the shower of a decay
   */
  virtual void showerDecay(ShowerTreePtr);

  /**
   *
   */
  void makeRemnants(ShowerTreePtr);

  /**
   *  Get the ShowerVariables
   */
  inline ShowerVarsPtr showerVariables() const;

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
   *  Generate the hard matrix element correction
   */
  virtual void hardMatrixElementCorrection();

  /**
   *  Extract the particles to be showered, set the evolution scales
   *  and apply the hard matrix element correction
   * @param Whether this is a hard process or decay
   * @return The particles to be showered
   */
  vector<ShowerProgenitorPtr> setupShower(bool hard);

  /**
   *  set the colour partners
   */
  void setColourPartners(bool hard);

  /**
   *  Methods to perform the evolution of an individual particle, including
   *  recursive calling on the products
   */
  //@{
  /**
   * It does the forward evolution of the time-like input particle
   * (and recursively for all its radiation products).
   * accepting only emissions which conforms to the showerVariables
   * and soft matrix element correction pointed by meCorrectionPtr.
   * If at least one emission has occurred then the method returns true.
   * @param particle The particle to be showered
   */
  virtual bool timeLikeShower(tShowerParticlePtr particle); 

  /**
   * It does the backward evolution of the space-like input particle 
   * (and recursively for all its time-like radiation products).
   * accepting only emissions which conforms to the showerVariables.
   * If at least one emission has occurred then the method returns true
   * @param particle The particle to be showered
   */
  virtual bool spaceLikeShower(tShowerParticlePtr particle); 

  /**
   * If does the forward evolution of the input on-shell particle
   * involved in a decay 
   * (and recursively for all its time-like radiation products).
   * accepting only emissions which conforms to the showerVariables.
   * @param particle    The particle to be showered
   * @param maxscale    The maximum scale for the shower.
   * @param minimumMass The minimum mass of the final-state system
   */
  virtual bool spaceLikeDecayShower(tShowerParticlePtr particle,
				    vector<Energy> maxscale,
				    Energy minimumMass);
  //@}

  /**
   * This method sets all the properties of the new particles from the 
   * splitting: it fixes the hadron parent/children relations due to the 
   * splitting and the colour information for a backward branching
   * @param The particle being evolved
   * @param newparent The new initial-state particle
   * @param otherChild The new final-state particle
   * @param scale The scale of the branching
   * @param z     The energy fraction
   * @param inter The interaction
   */
  void createBackwardBranching(ShowerParticlePtr part, ShowerParticlePtr newparent, 
			       ShowerParticlePtr otherChild, Energy scale,
			       double z,
			       ShowerIndex::InteractionType inter);

  /**
   * Set up the colour for the backward evolution
   * @param newParent The new initial-state particle
   * @param oldParent The old initial-state particle
   * @param otherChild The outgoing final-state particle
   */
  void setBackwardColour(ShowerParticlePtr &newParent,
			 ShowerParticlePtr &oldParent,
			 ShowerParticlePtr &otherChild);

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
  inline virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Evolver> initEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Evolver & operator=(const Evolver &);
  
private:

  /**
   * Pointer to the splitting generator
   */
  SplittingGeneratorPtr _splittingGenerator;

  /**
   * Pointer to the PartnerFinfer
   */
  PartnerFinderPtr _partnerFinder;

  /**
   * Vector of reference to the possible matrix element corrections
   */
  vector<MECorrectionPtr> _mecorrections;

  /**
   *  The progenitor of the current shower
   */
  ShowerProgenitorPtr _progenitor;

  /**
   *  Matrix element correction used for the current shower
   */
  MECorrectionPtr _currentme;
  
  /**
   *  Pointer to the ShowerVariables object
   */
  ShowerVarsPtr _showerVariables;

  /**
   *  Pointer to the object responsible for the kinematical reconstruction
   */
  KinematicsReconstructorPtr _kinematicsrecon;

  /**
   * The ShowerTree currently being showered
   */
  ShowerTreePtr _currenttree;

  /**
   *  Maximum number of tries to generate the shower of a particular tree
   */
  unsigned int _maxtry;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Evolver. */
template <>
struct BaseClassTrait<Herwig::Evolver,1> {
  /** Typedef of the first base class of Evolver. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Evolver class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Evolver>
  : public ClassTraitsBase<Herwig::Evolver> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::Evolver"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Evolver is implemented. It may also include several, space-separated,
   * libraries if the class Evolver depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "Evolver.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Evolver.tcc"
#endif

#endif /* HERWIG_Evolver_H */
