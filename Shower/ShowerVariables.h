// -*- C++ -*-
#ifndef HERWIG_ShowerVariables_H
#define HERWIG_ShowerVariables_H
//
// This is the declaration of the ShowerVariables class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerIndex.h"
#include "ShowerVariables.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is responsible for keeping all the constraint information 
 *  on the shower evolution. In particular, it has the scale value at 
 *  which to stop the shower. Here "scale" can be either the mass scale or the 
 *  \f$\tilde{q}\f$ (ordering variable) scale: 
 *  this class is also responsible for the
 *  conversion between these two different scale definitions.
 *
 *  Furthermore, this class can also have a veto for emission above a certain 
 *  \f$p_T\f$ scale, or a veto for emission below a certain \f$p_T\f$ scale, 
 *  where \f$p_T\f$ is the "resolution" variable. 
 *  This class has also two important switches: one for switching on/off 
 *  the multi-scale showering; and one for switching on/off the decay 
 *  of particles (mainly Susy ones) before showering. These two switches 
 *  are set by default, in Herwig++, to:
 *  -  multi-scale shower  1  (ON);
 *  -  decay before shower  0  (OFF).
 *  However, if you want the same behaviour as in Fortran Herwig, then set: 
 *  -  multi-scale shower   0  (OFF);
 *  -  decay before shower  1  (ON).
 * 
 *  In the case decay before shower is ON, the set of particles to be 
 *  decayed before showering are contained in the initialize() method: 
 *  you need to change this method if you want add/remove a particle. 
 *  Finally, this class has also three parameters to set the low energy 
 *  cutoff mass scales for respectively QCD, QED, EWK radiation. 
 *  The class provides also set/access to the upper scale for all 
 *  interaction types and events: it is supposed to be set, at 
 *  initialization time, by some other class, to the center of mass 
 *  energy of the beam-beam interaction, and used as upper scale value 
 *  for the numerically evaluation of Sudakov form factors. 
 *
 *  Notice that:
 *  -      to be more general, one should define an abstract class 
 *         AbsShowerVariables, which has, exactly like the present 
 *         ShowerVariables class, the definition of all methods, 
 *         but only the following two, which are declared as pure virtual methods 
 *         without implementation: 
 *       -- virtual ... convertMassScaleToQScale() = 0; 
 *       -- virtual ... convertQScaleToMassScale() = 0; 
 *         because it is only here that specific choice must be made 
 *         on which ordering variable we want to use. Then, the concrete class 
 *         ShowerVariables inherits from AbsShowerVariables 
 *         and provides a definition for those virtual methods. 
 *         Therefore, if we wanted a different choice, we could define another class,
 *         AlternativeShowerVariables, which also inherits from 
 *         AbsShowerVariables, but provides a different definition of those methods.
 *
 *
 * @see ShowerIndex
 *
 * @see \ref ShowerVariablesInterfaces "The interfaces"
 * defined for ShowerVariables.
 */
class ShowerVariables: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline ShowerVariables();

  /**
   * The copy constructor.
   */
  inline ShowerVariables(const ShowerVariables &);

  /**
   * The destructor.
   */
  virtual ~ShowerVariables();
  //@}

public:

  /**
   *  Access to the various switches
   */
  //@{
  /**
   * Access the multi-scale showering mode switch: <em>0 (OFF), 1 (ON).</em>
   * By choosing <em>0 (OFF)</em>, one gets a similar behaviour to  
   * Fortran Herwig, in which the showering is done in one go, 
   * from the starting scale to the cutoff.
   * The default for Herwig++ is <em>1 (ON)</em>: multi-scale showering.
   */
  inline int isMultiScaleShowerON() const;

  /**
   * Access the decay before shower mode switch: <em>0 (OFF), 1 (ON).</em>
   * By choosing <em>1 (ON)</em>, one gets a similar behaviour like in
   * the Fortran Herwig, in which some particles (mainly Susy
   * particles like gluinos, squarks,...) decay before showering.
   * The default for Herwig++ is <em>0 (OFF)</em>: decay and shower intermixed.
   */
  inline int isDecayBeforeShowerON() const;

  /**
   * It returns true if the particle with the specified id
   * is in the list of those that should be decayed before
   * showering. This method should be invoked only when the
   * decay before shower mode is <em>1 (ON)</em>. 
   */
  inline bool hasToDecayBeforeShower(const long id) const;

  /**
   * It returns the low energy cutoff <em>mass </em> scale for the 
   * interaction type specified in input.
   */
  Energy cutoffMassScale(const ShowerIndex::InteractionType interaction) const;

  /**
   * It returns the low energy cutoff \f$\tilde{q}\f$ scale for the 
   * interaction type specified in input.
   */
  Energy cutoffQScale(const ShowerIndex::InteractionType interaction) const;

  /**
   * Specifies a kinematic cutoff used in the parton shower phase space. 
   */
  inline Energy kinScale() const;
  //@}

  /**
   * It resets all the scales, and vetos.
   */
  void reset();

  /**
   *  Conversion begtween the scales.
   */
  //@{
  /**
   * Conversion between <em>mass</em> and \f$\tilde{q}\f$ scale.
   */
  inline Energy convertMassScaleToQScale(const Energy inputMassScale) const;

  /**
   * Conversion between \f$\tilde{q}\f$ and <em>mass</em> scale.
   */
  inline Energy convertQScaleToMassScale(const Energy inputQScale) const;
  //@}

  /**
   *   Access/set the <em>mass</em> and \f$\tilde{q}\f$ scales
   *   at which to shower the showering.
   */
  //@{
  /**
   *  Get the <em>Mass</em> scale at which to stop the showering
   */
  inline Energy stopShowerAtMassScale() const;

  /**
   *  Set the <em>Mass</em> scale at which to stop the showering
   * @param inputStopShowerAtMassScale scale at which to stop
   */
  inline void stopShowerAtMassScale(const Energy inputStopShowerAtMassScale);

  /**
   *  Get the \f$\tilde{q}\f$ scale at which to stop the showering
   */
  inline Energy stopShowerAtQScale() const;

  /**
   *  Set the \f$\tilde{q}\f$ scale at which to stop the showering
   * @param inputStopShowerAtQScale scale at which to stop
   */
  inline void stopShowerAtQScale(const Energy inputStopShowerAtQScale);
  //@}

  /**
   * Access/set the Veto in \f$p_T\f$ (resolution) scale.
   */
  //@{
  /**
   * Get \f$p_T\f$ scale above which to veto emissions.
   */
  inline Energy vetoAbovePtScale() const;

  /**
   * Set \f$p_T\f$ scale above which to veto emissions.
   * @param inputVetoAbovePtScale scale for the veto
   */
  inline void vetoAbovePtScale(const Energy inputVetoAbovePtScale);

  /**
   * Get \f$p_T\f$ scale below which to veto emissions.
   */
  inline Energy vetoBelowPtScale() const;

  /**
   * Set \f$p_T\f$ scale below which to veto emissions.
   * @param inputVetoBelowPtScale scale for the veto
   */
  inline void vetoBelowPtScale(const Energy inputVetoBelowPtScale);
  //@}

  /**
   * Use to initialize some scales.
   */
  static Energy HUGEMASS; 

  /**
   * Access/set \f$p_T\f$ of hardest emission so far.
   */
  //@{
  /**
   *  Get the scale of hardest emission from the quark.
   */
  inline Energy largestPtQ() const;

  /**
   *  Set the scale of hardest emission from the quark.
   * @param pt The scale
   */
  inline void setLargestPtQ(const Energy pt);

  /**
   *  Get the scale of hardest emission from the antiquark.
   */
  inline Energy largestPtQbar() const;

  /**
   *  Set the scale of hardest emission from the antiquark.
   * @param pt The scale
   */
  inline void setLargestPtQbar(const Energy pt);
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

  /**
   * Assign asymmetric initial condition to parton shower, random or
   * not? If not random, then quark gets larger initial scale.
   */
  inline bool asyPS() const;

  /**
   * Asymmetric parton shower phase space, random choice for jet with
   * large initial scale?
   */
  inline bool rndPS() const;
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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**  
   * Build the set of particles which should decay before shower
   * (only when the decay before shower mode is 1 (ON) ) 
   */
  void initialize();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ShowerVariables> initShowerVariables;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerVariables & operator=(const ShowerVariables &);

private:

  /**
   * The switch for on/off multi-scale shower
   */
  int _multiScaleShowerMode; 

  /**
   * The switch for on/off decay before shower
   */
  int _decayBeforeShowerMode;

  /** 
   * Low-energy cutoff mass scale for QCD radiation
   */
  Energy _cutoffQCDMassScale;
 
  /**
   * Low-energy cutoff mass scale for QED radiation
   */
  Energy _cutoffQEDMassScale;

  /**
   * Low-energy cutoff mass scale for EWK radiation
   */
  Energy _cutoffEWKMassScale; 

  /**
   * Kinematic cutoff used in the parton shower phase space. 
   */
  Energy _kinCutoffScale; 

  /**
   * Matrix element correction switch
   */
  int _meCorrMode; 

  /**
   *  Initial conditions for the shower
   */
  int _qqgPSMode; 

  /**
   *  Mass cut-off for the shower
   */
  Energy _stopShowerAtMassScale;

  /**
   *  Veto emissions above this \f$p_T\f$ scale.
   */
  Energy _vetoAbovePtScale;

  /**
   *  Veto emissions below this \f$p_T\f$ scale.
   */
  Energy _vetoBelowPtScale;

  /**
   *  \f$p_T\f$ scale of the hardest emission from the quark
   */
  Energy _largestPtQ;

  /**
   * \f$p_T\f$ scale of the hardest emission from the antiquark
   */
  Energy _largestPtQbar;

  /**
   *  PDG codes of the particles which decay before showering
   */
  set<long> _particlesDecayBeforeShower;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerVariables. */
template <>
struct BaseClassTrait<Herwig::ShowerVariables,1> {
  /** Typedef of the first base class of ShowerVariables. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerVariables class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerVariables>
  : public ClassTraitsBase<Herwig::ShowerVariables> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::ShowerVariables"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the ShowerVariables class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ShowerVariables.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerVariables.tcc"
#endif

#endif /* HERWIG_ShowerVariables_H */
