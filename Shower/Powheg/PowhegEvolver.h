// -*- C++ -*-
#ifndef HERWIG_PowhegEvolver_H
#define HERWIG_PowhegEvolver_H
//
// This is the declaration of the PowhegEvolver class.
//

#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/HardBranching.fh"
#include "PowhegEvolver.h"
#include "HardestEmissionGenerator.h"
#include "Herwig++/Utilities/Histogram.h"


namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the PowhegEvolver class.
 *
 * @see \ref PowhegEvolverInterfaces "The interfaces"
 * defined for PowhegEvolver.
 */
class PowhegEvolver: public Evolver {

public:

  /**
   * The default constructor.
   */
  PowhegEvolver() : _hardonly(false), _trunc_Mode(true) {}

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
   *  Start the shower of a timelike particle
   */
  virtual bool startTimeLikeShower(ShowerInteraction::Type);

  /**
   *  Start the shower of a spacelike particle
   */
  virtual bool startSpaceLikeShower(PPtr, ShowerInteraction::Type);

  /**
   *  Start the shower of a spacelike decaying aparticle
   */
  virtual bool startSpaceLikeDecayShower(Energy maxscale,Energy minimumMass,
					 ShowerInteraction::Type);

protected:
  
  /**
   * Is the truncated shower on?
   */
  bool isTruncatedShowerON() const {return _trunc_Mode;}

  /**
   * Extract the particles to be showered, set the evolution scales
   * and apply the hard matrix element correction
   * @param hard Whether this is a hard process or decay
   * @return The particles to be showered
   */
  virtual vector<ShowerProgenitorPtr> setupShower(bool hard);

  /**
   * Implementation of checks on momentum reconstruction
   */  
  virtual bool checkShowerMomentum( vector<ShowerProgenitorPtr> particlesToShower );

  /**
   *  set the colour partners
   */
  virtual void setEvolutionPartners(bool hard,ShowerInteraction::Type);

  /**
   *  Generate the hardest emission
   */
  virtual void hardestEmission();

  /**
   * Truncated shower from a time-like particle
   */
  virtual bool truncatedTimeLikeShower(tShowerParticlePtr particle,
				       HardBranchingPtr branch,
				       ShowerInteraction::Type type);

  /**
   * Truncated shower from a time-like particle
   */
  virtual bool truncatedSpaceLikeDecayShower(tShowerParticlePtr particle,
					     Energy maxscale, Energy minimumMass,
					     HardBranchingPtr branch,
					     ShowerInteraction::Type type);
 
  /**
   * Truncated shower from a space-like particle
   */
  virtual bool truncatedSpaceLikeShower(tShowerParticlePtr particle,PPtr beam,
					HardBranchingPtr branch,
					ShowerInteraction::Type type);

  /**
   *  Access to set/get the HardTree currently beinging showered
   */
  //@{
  /**
   *  The HardTree currently being showered
   */
  inline tHardTreePtr hardTree() {return _nasontree;}

  /**
   *  The HardTree currently being showered
   */
  inline void hardTree(tHardTreePtr in) {_nasontree = in;}
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}
  
  virtual void doinitrun();

  virtual void dofinish();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PowhegEvolver> initPowhegEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegEvolver & operator=(const PowhegEvolver &);

private:

  /**
   *  Vector of objects responisble for generating the hardest emission
   */
  vector<HardestEmissionGeneratorPtr> _hardgenerator;

  /**
   *  The HardTree currently being showered
   */
  HardTreePtr _nasontree;

  /**
   *  Only generate the emission from the hardest emission
   *  generate for testing only
   */
  bool _hardonly;

 /**
   *  Truncated shower switch
   */
  bool _trunc_Mode;
  
  /**
   *  Count of the number of truncated emissions
   */
  unsigned int _truncEmissions;
  
  /**
   *  Histograms of momentum differences in reconstructed momenta
   */
  HistogramPtr _h_Xdiff;
  HistogramPtr _h_Ydiff;
  HistogramPtr _h_Zdiff;
  HistogramPtr _h_Ediff;

  /**
   * Count of events passing momenta reconstruction acceptance
   */
  int _no_events;
  int _mom_fails;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PowhegEvolver. */
template <>
struct BaseClassTrait<Herwig::PowhegEvolver,1> {
  /** Typedef of the first base class of PowhegEvolver. */
  typedef Herwig::Evolver NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PowhegEvolver class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PowhegEvolver>
  : public ClassTraitsBase<Herwig::PowhegEvolver> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PowhegEvolver"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PowhegEvolver is implemented. It may also include several, space-separated,
   * libraries if the class PowhegEvolver depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_PowhegEvolver_H */

