// -*- C++ -*-
#ifndef HERWIG_Evolver_H
#define HERWIG_Evolver_H
//
// This is the declaration of the Evolver class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "ShowerModel.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ShowerTree.fh"
#include "MECorrectionBase.fh"
#include "ShowerProgenitor.fh"
#include "Herwig++/Shower/ShowerHandler.fh"
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

/**
 *  The MECorrectionBase class is a friend to access some protected
 *  member functions to set the radiation enhancement factors
 */
friend class MECorrectionBase;

public:

  /**
   *  Default Constructor
   */
  inline Evolver();

  /**
   *  Member to perform the shower
   */
  //@{
  /**
   * Perform the shower of the hard process
   */
  virtual void showerHardProcess(ShowerTreePtr);

  /**
   * Perform the shower of a decay
   */
  virtual void showerDecay(ShowerTreePtr);
  //@}

  /**
   *  Access to the flags and shower variables
   */
  //@{
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
   *  Get the ShowerModel
   */ 
  inline ShowerModelPtr showerModel() const;
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
   *  Generate the hard matrix element correction
   */
  virtual void hardMatrixElementCorrection();

  /**
   * Extract the particles to be showered, set the evolution scales
   * and apply the hard matrix element correction
   * @param hard Whether this is a hard process or decay
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
   * @param beam The beam particle
   */
  virtual bool spaceLikeShower(tShowerParticlePtr particle,PPtr beam); 

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
   *  Switches for matrix element corrections
   */
  //@{
  /**
   * Any ME correction?   
   */
  inline bool MECOn() const;

  /**
   * Any hard ME correction? 
   */
  inline bool hardMEC() const;

  /**
   * Any soft ME correction? 
   */
  inline bool softMEC() const;
  //@}

  /**
   *  Switches for vetoing hard emissions
   */
  //@{
  /**
   * Vetos on? 
   */
  inline bool hardVetoOn() const;

  /**
   * veto hard emissions in IS shower?
   */
  inline bool hardVetoIS() const;

  /**
   * veto hard emissions in FS shower?
   */
  inline bool hardVetoFS() const;
  //@}

  /**
   *  Enhancement factors for radiation needed to generate the soft matrix
   *  element correction.
   */
  //@{
  /**
   *  Access the enhancement factor for initial-state radiation
   */
  inline double initialStateRadiationEnhancementFactor() const;

  /**
   *  Access the enhancement factor for final-state radiation
   */
  inline double finalStateRadiationEnhancementFactor() const;

  /**
   *  Set the enhancement factor for initial-state radiation
   */
  inline void initialStateRadiationEnhancementFactor(double);

  /**
   *  Set the enhancement factor for final-state radiation
   */
  inline void finalStateRadiationEnhancementFactor(double);
  //@}

  /**
   * Access/set the beam particle for the current initial-state shower
   */
  //@{
  /**
   *  Get the beam particle data
   */
  inline Ptr<BeamParticleData>::const_pointer beamParticle() const;

  /**
   *  Set the beam particle data
   */
  inline void setBeamParticle(Ptr<BeamParticleData>::const_pointer);
  //@}

  /**
   *  Calculate the intrinsic \f$p_T\f$.
   */
  virtual void generateIntrinsicpT(vector<ShowerProgenitorPtr>);

  /**
   *  find the maximally allowed pt acc to the hard process. 
   */
  void setupMaximumScales(ShowerTreePtr, vector<ShowerProgenitorPtr>);

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
  static ClassDescription<Evolver> initEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Evolver & operator=(const Evolver &);
  
private:

  /**
   *  Pointer to the model for the shower evolution model
   */
  ShowerModelPtr _model;

  /**
   * Pointer to the splitting generator
   */
  SplittingGeneratorPtr _splittingGenerator;

  /**
   *  Maximum number of tries to generate the shower of a particular tree
   */
  unsigned int _maxtry;

  /**
   * Matrix element correction switch
   */
  unsigned int _meCorrMode; 

  /**
   * Hard emission veto switch
   */
  unsigned int _hardVetoMode; 

  /**
   *  The progenitor of the current shower
   */
  ShowerProgenitorPtr _progenitor;

  /**
   *  Matrix element correction used for the current shower
   */
  MECorrectionPtr _currentme;

  /**
   * The ShowerTree currently being showered
   */
  ShowerTreePtr _currenttree;

  /**
   *  Radiation enhancement factors for use with the veto algorithm
   *  if needed by the soft matrix element correction 
   */
  //@{
  /**
   *  Enhancement factor for initial-state radiation
   */
  double _initialenhance;

  /**
   *  Enhancement factor for final-state radiation
   */
  double _finalenhance;
  //@}

  /**
   *  The beam particle data for the current initial-state shower
   */
  Ptr<BeamParticleData>::const_pointer _beam;

  /**
   *  Storage of the intrinsic \f$p_t\f$ of the particles
   */
  map<tShowerProgenitorPtr,pair<Energy,double> > _intrinsic;
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
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "Evolver.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Evolver.tcc"
#endif

#endif /* HERWIG_Evolver_H */
