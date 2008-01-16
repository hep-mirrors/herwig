// -*- C++ -*-
#ifndef HERWIG_NasonEvolver_H
#define HERWIG_NasonEvolver_H
//
// This is the declaration of the NasonEvolver class.
//

#include "Herwig++/Shower/Base/Evolver.h"
#include "NasonEvolver.h"
#include "HardestEmissionGenerator.h"
#include "Herwig++/Utilities/Histogram.h"


namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the NasonEvolver class.
 *
 * @see \ref NasonEvolverInterfaces "The interfaces"
 * defined for NasonEvolver.
 */
class NasonEvolver: public Evolver {

public:

  /**
   * The default constructor.
   */
  inline NasonEvolver();

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
 
  /**
   * Is the truncated shower on?
   */
  inline bool isTruncatedShowerON() const;

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
   * Extract the particles to be showered, set the evolution scales
   * and apply the hard matrix element correction
   * @param hard Whether this is a hard process or decay
   * @return The particles to be showered
   */
  vector<ShowerProgenitorPtr> setupShower(bool hard);

  /**
   *  Generate the hardest emission
   */
  virtual void hardestEmission();

  /**
   * It does the forward evolution of the time-like input particle
   * (and recursively for all its radiation products).
   * accepting only emissions which conforms to the showerVariables
   * and soft matrix element correction pointed by meCorrectionPtr.
   * If at least one emission has occurred then the method returns true.
   * @param particle The particle to be showered
   */
  virtual bool truncatedTimeLikeShower(tShowerParticlePtr particle,
				       NasonBranchingPtr branch);
 
  virtual bool truncatedSpaceLikeShower(tShowerParticlePtr particle,PPtr beam,
					NasonBranchingPtr branch); 

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
  inline virtual void doinit() throw(InitException);
  //@}
    /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<NasonEvolver> initNasonEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NasonEvolver & operator=(const NasonEvolver &);

private:

  /**
   *  Vector of objects responisble for generating the hardest emission
   */
  vector<HardestEmissionGeneratorPtr> _hardgenerator;

  /**
   *  The NasonTree currently being showered
   */
  NasonTreePtr _nasontree;

 /**
   *  Truncated shower switch
   */
  bool _trunc_Mode;

   /**
   *  Histogram object to record the number of truncated emissions
   */
  HistogramPtr _hTrunc;
  
  
  /**
   *  Count of the number of truncated emissions
   */
  unsigned int _truncEmissions;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NasonEvolver. */
template <>
struct BaseClassTrait<Herwig::NasonEvolver,1> {
  /** Typedef of the first base class of NasonEvolver. */
  typedef Herwig::Evolver NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NasonEvolver class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NasonEvolver>
  : public ClassTraitsBase<Herwig::NasonEvolver> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::NasonEvolver"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NasonEvolver is implemented. It may also include several, space-separated,
   * libraries if the class NasonEvolver depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#include "NasonEvolver.icc"

#endif /* HERWIG_NasonEvolver_H */
