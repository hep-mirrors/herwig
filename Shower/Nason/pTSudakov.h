// -*- C++ -*-
#ifndef HERWIG_pTSudakov_H
#define HERWIG_pTSudakov_H
//
// This is the declaration of the pTSudakov class.
//

#include "Herwig++/Shower/Base/SudakovFormFactor.h"
#include "pTSudakov.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the pTSudakov class.
 *
 * @see \ref pTSudakovInterfaces "The interfaces"
 * defined for pTSudakov.
 */
class pTSudakov: public SudakovFormFactor {

public:

  /**
   * The default constructor.
   */
  inline pTSudakov();

  /**
   * Constructor from anyother SudakovFormFactor
   */
  inline pTSudakov(tSudakovPtr);

  /**
   *  Members to generate the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param cc Whether this is the charge conjugate of the branching
   * @param enhance The radiation enhancement factor
   * defined.
   */
  virtual ShoKinPtr generateNextTimeBranching(const Energy startingScale,
					      const IdList &ids,const bool cc,
					      double enhance);

  /**
   * Return the scale of the next space-like decay branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param stoppingScale stopping scale for the evolution
   * @param minmass The minimum mass allowed for the spake-like particle.
   * @param ids The PDG codes of the particles in the splitting
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   * @param enhance The radiation enhancement factor
   */
  virtual ShoKinPtr generateNextDecayBranching(const Energy startingScale,
					       const Energy stoppingScale,
					       const Energy minmass,
					       const IdList &ids,
					       const bool cc,
					       double enhance);

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   * @param beam The beam particle
   * @param enhance THe radiation enhancement factor
   */
  virtual ShoKinPtr generateNextSpaceBranching(const Energy startingScale,
					       const IdList &ids,double x,
					       const bool cc,double enhance,
					       Ptr<BeamParticleData>::transient_const_pointer beam);

  /**
   *  set the maximum mass of the branching
   */
  inline void setQ2Max(Energy2);

  /**
   *  Method to return the evolution scale given the
   *  transverse momentum, \f$p_T\f$ and \f$z\f$.
   */
  virtual Energy calculateScale(double z, Energy pt, IdList ids,unsigned int iopt);

  /**
   *  Method to create the ShowerKinematics object for a final-state branching
   */
  virtual ShoKinPtr createFinalStateBranching(Energy scale,double z,
					      double phi, Energy pt);

  /**
   *  Method to create the ShowerKinematics object for an initial-state branching
   */
  virtual ShoKinPtr createInitialStateBranching(Energy scale,double z,
						double phi, Energy pt);

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
   *  Initialize the values of the cut-offs and scales
   * @param tmin The minimum scale
   * @param ids  The ids of the partics in the branching
   * @param cc Whether this is the charge conjugate of the branching
   */
  void initialize(const IdList & ids,Energy2 &tmin, Energy2 tmax,const bool cc);
  
  bool guessTimeLike(Energy2 &t,Energy2 tmin);

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<pTSudakov> initpTSudakov;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  pTSudakov & operator=(const pTSudakov &);

private:

  /**
   *  The evolution scale, \f$\tilde{q}\f$.
   */
  Energy _q;

  /**
   *  The Ids of the particles in the current branching
   */
  IdList _ids;

  /**
   *  The masses of the particles in the current branching
   */
  vector<Energy> _masses;

  /**
   *  The mass squared of the particles in the current branching
   */
  vector<Energy2> _masssquared;

  /**
   *  Kinematic cut-off
   */
  Energy _kinCutoff;

  /**
   *  The maximum off-shell mass
   */
  Energy2 _q2max;

  /**
   *  Transverse momentum squared
   */
  Energy2 _pt2;

  /**
   *  Mass squared
   */
  Energy2 _q2;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of pTSudakov. */
template <>
struct BaseClassTrait<Herwig::pTSudakov,1> {
  /** Typedef of the first base class of pTSudakov. */
  typedef Herwig::SudakovFormFactor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the pTSudakov class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::pTSudakov>
  : public ClassTraitsBase<Herwig::pTSudakov> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::pTSudakov"; }
  /**
   * The name of a file containing the dynamic library where the class
   * pTSudakov is implemented. It may also include several, space-separated,
   * libraries if the class pTSudakov depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#include "pTSudakov.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "pTSudakov.tcc"
#endif

#endif /* HERWIG_pTSudakov_H */
