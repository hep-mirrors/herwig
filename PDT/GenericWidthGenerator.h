// -*- C++ -*-
#ifndef HERWIG_GenericWidthGenerator_H
#define HERWIG_GenericWidthGenerator_H
//
// This is the declaration of the GenericWidthGenerator class.
//
#include "ThePEG/PDT/WidthGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "GenericWidthGenerator.fh"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "GenericMassGenerator.h"
#include <iostream>

namespace Herwig {
using namespace ThePEG;

  /** \ingroup PDT
   *
   * The <code>GenericWidthGenerator</code> class is designed to automatically
   * calculate the running width for a given particle using information from
   * the decayModes and the Decayers to construct the running width.
   *
   *  It also gives us the option of selecting the decay modes for a particle
   *  based on the mass.
   *
   * @see WidthGenerator
   * @see DecayIntegrator
   * @see GenericMassGenerator
   */
class GenericWidthGenerator: public WidthGenerator {

public:

  /**
   * Typedef to define a DecayMoap
   */
  typedef Selector<tDMPtr> DecayMap;

  /**
   * A friend class so the off-shell matrix elements can be integrated.
   */
  friend class TwoBodyAllOnCalculator;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline GenericWidthGenerator();

  /**
   * Copy constructor
   */
  inline GenericWidthGenerator(const GenericWidthGenerator &);

  /**
   * Destructor
   */
  virtual ~GenericWidthGenerator();
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

public:

  /**
   * Return true if this width generator can handle the given particle type.
   * @param part The particle data pointer of the particle.
   * @return True if this class can handle the particle and false otherwise
   */
  virtual bool accept(const ParticleData & part) const;

  /** @name Members to calculate the width and decay modes. */
  //@{
  /**
   * Calculate the width.
   * @param part The particle data pointer of the particle.
   * @param m The scale for the width calculation
   * @return The width at the mass given.
   */
  virtual Energy width(const ParticleData & part, Energy m) const;

  /**
   * Initialize the given decay map for the given particle type.
   * @param part The particle data pointer of the particle.
   * @return The decay map
   */
  virtual DecayMap rate(const ParticleData & part) const;

  /**
   * Return a decay map for a given particle instance. This allows us to
   * vary the branching ratios as a function of the particles mass.
   * @param  part The particle instance
   * @return The decay map
   */
  virtual DecayMap rate (const Particle & part);

  /**
   * The partial width for a given mode
   * @param m The mass, or scale, for the calculation
   * @param iloc The location of the mode in the list
   * @return The partial width for the mode.
   */
  Energy partialWidth(int iloc,Energy m) const;
  //@}

  /**
   * Output the initialisation info for the database
   * @param output The stream to output the information to
   * @param header output the header.
   **/
  virtual void dataBaseOutput(ofstream & output, bool header=true);

protected:

  /**
   * The \f$1\to2\f$ width for on-shell particles
   * @param m The mass, or scale, for the calculation
   * @param iloc The location of the mode in the list.
   * @return The partial width.
   */
  inline Energy partial2BodyWidth(int iloc,Energy m) const;

  /**
   * The \f$1\to2\f$ width for outgoing particles which can be off-shell.
   * @param iloc The location of the mode in the list.
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the first outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @return The partial width.
   */
  inline virtual Energy partial2BodyWidth(int iloc,Energy m0,Energy m1,Energy m2) const;

  /**
   * Perform the set up for a mode in classes inheriting from this one
   * @param mode The decay mode
   * @param decayer The decayer for the mode.
   * @param imode The number of the mode.
   */
  virtual void setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer, unsigned int imode);

  /**
   *  Access to the particle dat for inheriting classes
   */
  inline tPDPtr particle() const;

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
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
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

  /**
   * set up the interpolators
   */
  void setInterpolators();

  /**
   *  Matrix element code for a given mode 
   * @param imode The mode.
   */
  inline int MEcode(int imode) const;

  /**
   *  Coupling for a given mode
   * @param imode The mode.
   */
  inline double MEcoupling(int imode) const;


  /**
   *  The on-shell mass of the particle
   */
  inline Energy mass() const;

  /**
   * Initialization option for use by the inheriting classes
   */
  inline bool initialize() const;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<GenericWidthGenerator> initGenericWidthGenerator;

  /**
   * Private and non-existent assignment operator.
   */
  GenericWidthGenerator & operator=(const GenericWidthGenerator &);

private:

  /**
   * The pointer to the ParticleData object for the particle for this width generator.
   */
  PDPtr _theParticle;

  /**
   * The decaymodes
   */
  vector<DMPtr> _decaymodes;

  /**
   *  The minimum mass of the decaying particle for which this decay mode is possible
   */
  vector<Energy> _minmass;

  /**
   * The on-shell mass of the particle
   */
  Energy _mass;

  /**
   * Prefactor to get the on-shell width
   */
  double _prefactor;

  /**
   * The type of ME, whether it is fixed, calculated by this class or interpolation
   */
  vector<int> _MEtype;

  /**
   * The code for the matrix element
   */
  vector<int> _MEcode;

  /**
   *  Mass of the first outgoing particle for the simple \f$1\to2\f$ ME's
   */
  vector<Energy> _MEmass1;
  /**
   *  Mass of the second outgoing particle for the simple \f$1\to2\f$ ME's
   */
  vector<Energy> _MEmass2;

  /**
   * the coupling for a given me
   */
  vector<double> _MEcoupling; 

  /**
   * is this mode used for the running width
   */
  vector<bool> _modeon;

  /**
   * storage of the massesto set up the interpolation tables
   */
  vector<Energy> _intermasses;

  /**
   * storage of the widths to set up the interpolation tables
   */
  vector<Energy> _interwidths;

  /**
   * the number of entries in the decay table for a particular mode
   */
  vector<int> _noofentries;

  /**
   * initialize the generator using the particle data object
   */
  bool _initialize;

  /**
   * normalise the terms so that the partial widths for an on-shell particle are correct
   */
  bool _BRnorm;

  /**
   * number of points to use for interpolation tables
   */
  int _npoints;

  /**
   * intepolators for the running width
   */
  vector<Interpolator*> _interpolators;

  /**
   * minimum branching ratio for the inclusion in the total running width
   */
  double _BRminimum;

  /**
   *
   */
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of GenericWidthGenerator.
 */
template <>
struct BaseClassTrait<Herwig::GenericWidthGenerator,1> {
  /** Typedef of the base class of GenericWidthGenerator. */
  typedef WidthGenerator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::GenericWidthGenerator>
  : public ClassTraitsBase<Herwig::GenericWidthGenerator> {
   /** Return the class name. */
   static string className() { return "Herwig++::GenericWidthGenerator"; }
   /**
    * Return the name of the shared library to be loaded to get
    * access to this class and every other class it uses
    * (except the base class).
    */
   static string library() { return "libHwPDT.so"; }
};

}

#include "GenericWidthGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GenericWidthGenerator.tcc"
#endif

#endif /* HERWIG_GenericWidthGenerator_H */
