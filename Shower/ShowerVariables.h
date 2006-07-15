// -*- C++ -*-
#ifndef HERWIG_ShowerVariables_H
#define HERWIG_ShowerVariables_H
//
// This is the declaration of the ShowerVariables class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerConfig.h"
#include "Couplings/ShowerIndex.h"
#include "ShowerVariables.fh"
#include "ThePEG/PDF/BeamParticleData.h"

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
 *
 *  The ShowerHandler by default will decay all the unstable particles
 *  specified in the relevant interfacing during the shower however the 
 *  treatment of the shower from these particles can be either using
 *  a multi-scale approach of not depending on the switch setting
 * 
 *  Finally, this class has also three parameters to set the low energy 
 *  cutoff mass scales for respectively QCD, QED, EWK radiation. 
 *  The class provides also set/access to the upper scale for all 
 *  interaction types and events: it is supposed to be set, at 
 *  initialization time, by some other class, to the center of mass 
 *  energy of the beam-beam interaction, and used as upper scale value 
 *  for the numerically evaluation of Sudakov form factors. 
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
   * It returns true if the particle with the specified id
   * is in the list of those that should be decayed during the showering
   * showering.
   */
  inline bool decayInShower(const long id) const;

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
   * Specifies the kinematic cutoff used in the parton shower phase space. 
   */
  inline Energy kinScale() const;
  //@}

  /**
   *  Conversion between the scales.
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
   * Access/set the PDF and beam particle for the current initial-state shower
   */
  //@{
  /**
   *  Get the PDF
   */
  inline tcPDFPtr currentPDF() const;

  /**
   *  Set the PDF
   */
  inline void setCurrentPDF(tcPDFPtr);

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

  /**
   * The virtuality cut-off on the gluon \f$Q_g=\frac{\delta-am_q}{b}\f$
   * @param scale The scale \f$\delta\f$
   * @param mq The quark mass \f$m_q\f$.
   */
  inline Energy kinematicCutOff(Energy scale, Energy mq) const;

  /**
   *  Set/Get the gluon mass which should be used in the reconstruction
   */
  //@{
  /**
   *  Get the mass
   */
  inline Energy gluonMass() const;

  /**
   *  Set the mass
   */
  inline void gluonMass(Energy);
  //@}

  /**
   *  Get the GlobalParameters
   */
  inline GlobalParametersPtr globalParameters() const;

  /**
   * Set the gluon mass 
   * @param final If final is true gluon mass will be set to 
   * 0 or the effective mass depending on the choice of hadronisation model.
   * otherwise if the multi-scale shower is on the mass will be set to zero
   */
  inline void setGluonMass(bool final);


  /**
   *  Enhancement factors for radiation needed to generate the soft matrix element
   *  correction.
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
   *  Access functions for the type of shower phase space partition.
   *  These set/return the whether the so-called 'symmetric'/'maximal'
   *  /'smooth' choice was used (see _decay_shower_partition below).
   *  Also we have a similar function which returns whether the T2
   *  region is to be populated by the ME correction or the shower from
   *  the decaying particle.
   */
  //@{
  /**
   *  Access the option which determines the type of phase space partitioning 
   *  for the decay_shower.
   */
  inline unsigned int decay_shower_partition() const;

  /**
   *  Access the option denoting whether the T2 region of the decay phase
   *  space is populated by the shower (default, false) or the ME correction. 
   */
  inline bool use_me_for_t2();

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

  /**
   * Use to initialize some scales.
   */
  static const Energy HUGEMASS; 

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

  /** @name Standard Interfaced functions. */  //@{
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
   *  PDG codes of the particles which decay during showering
   *  this is fast storage for use during running
   */
  set<long> _particlesDecayInShower;

  /**
   *  PDG codes of the particles which decay during showering
   *  this is a vector that is interfaced so they can be changed
   */
  vector<long> _inputparticlesDecayInShower;

  //@{
  /**
   *  The PDF being used for the current initial-state shower
   */
  tcPDFPtr _pdf;

  /**
   *  The beam particle data for the current initial-state shower
   */
  Ptr<BeamParticleData>::const_pointer _beam;
  //@}

  /**
   *  Parameters for the \f$Q_g=\max(\frac{\delta-am_q}{b},c)\f$ kinematic cut-off
   */
  //@{
  /**
   *  The \f$a\f$ parameter
   */
  double _a;

  /**
   *  The \f$b\f$ parameter
   */
  double _b;

  /**
   *  The \f$c\f$ parameter
   */
  Energy _c;
  //@}

  /**
   *  The gluon mass
   */
  Energy _gluonMass;

  /**
   *  The global variables
   */
  GlobalParametersPtr _globalParameters;

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
   *  The following variables relate to the decay shower and
   *  its associated ME corrections.
   */

  /**
   *  Here we hold the option which determines the type of phase space 
   *  partitioning to do for the decay shower. Depending on the value
   *  of the option PartnerFinder will set different bounds on the starting
   *  $\tilde{q}$ values for the showers of the decaying particle and its
   *  charged child. This is done according to the top decay colour 
   *  connection calculation in JHEP12(2003)_045. The options act as follows:
   *  0: This is the default 'symmetric' choice which more or less divides
   *     the phase space evenly between the parent and its charged child.
   *  1: This 'maximal' choice maximises the phase space available for 
   *     gluons emitted from the charged child.
   *  2: This (experimental) 'smooth' choice does not suffer from
   *     a discontinuity at the boundary between the region populated by
   *     emissions from the charged child and the region populated by emissions
   *     from the parent. This does, however, mean that the phase space 
   *     available for emissions from the charged child is fairly minimal.
   */
  unsigned int _decay_shower_partition;

  /**
   *  This flag determines whether the T2 region in the decay shower
   *  (JHEP12(2003)_045) is populated by the ME correction (true) or
   *  the shower from the decaying particle.
   */
  bool _use_me_for_t2;

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
  /**
   * The name of a file containing the dynamic library where the class
   * ShowerVariables is implemented. It may also include several, space-separated,
   * libraries if the class ShowerVariables depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "ShowerVariables.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerVariables.tcc"
#endif

#endif /* HERWIG_ShowerVariables_H */
