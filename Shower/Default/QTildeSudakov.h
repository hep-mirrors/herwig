// -*- C++ -*-
//
// QTildeSudakov.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_QTildeSudakov_H
#define HERWIG_QTildeSudakov_H
//
// This is the declaration of the QTildeSudakov class.
//

#include "Herwig++/Shower/Base/SudakovFormFactor.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 * The QTildeSudakov class implements the Sudakov form factor for evolution in
 * \f$\tilde{q}^2\f$ using the veto algorithm.
 *
 * @see \ref QTildeSudakovInterfaces "The interfaces"
 * defined for QTildeSudakov.
 */
class QTildeSudakov: public SudakovFormFactor {

public:

  /**
   * The default constructor.
   */
  inline QTildeSudakov() {}
  
  /**
   *  Members to generate the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   * @param enhance The radiation enhancement factor
   */
  virtual ShoKinPtr generateNextTimeBranching(const Energy startingScale,
					      const IdList &ids,const bool cc,
					      double enhance);

  /**
   * Return the scale of the next space-like decay branching. If there is no 
   * branching then it returns ZERO.
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
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   * @param enhance The radiation enhancement factor
   * @param beam The beam particle
   */
  virtual ShoKinPtr generateNextSpaceBranching(const Energy startingScale,
					       const IdList &ids,double x,
					       const bool cc, double enhance,
					       tcBeamPtr beam);
  //@}

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

  /**
   *  Method to create the ShowerKinematics object for a decay branching
   */
  virtual ShoKinPtr createDecayBranching(Energy scale,double z,
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
   *  Methods to provide the next value of the scale before the vetos
   *  are applied.
   */
  //@{
  /**
   *  Value of the energy fraction and scale for time-like branching
   * @param t  The scale
   * @param tmin The minimum scale
   * @param enhance The radiation enhancement factor
   * @return False if scale less than minimum, true otherwise
   */
  bool guessTimeLike(Energy2 &t, Energy2 tmin, double enhance);

  /**
   * Value of the energy fraction and scale for time-like branching
   * @param t  The scale
   * @param tmax The maximum scale
   * @param minmass The minimum mass of the particle after the branching
   * @param enhance The radiation enhancement factor
   */
  bool guessDecay(Energy2 &t, Energy2 tmax,Energy minmass,
		  double enhance);

  /**
   * Value of the energy fraction and scale for space-like branching
   * @param t  The scale
   * @param tmin The minimum scale
   * @param x Fraction of the beam momentum.
   * @param enhance The radiation enhancement factor
   */
  bool guessSpaceLike(Energy2 &t, Energy2 tmin, const double x,
		      double enhance);
  //@}

  /**
   *  Initialize the values of the cut-offs and scales
   * @param tmin The minimum scale
   * @param ids  The ids of the partics in the branching
   * @param cc Whether this is the charge conjugate of the branching
   */
  void initialize(const IdList & ids,Energy2 &tmin, const bool cc);

  /**
   *  Phase Space veto member to implement the \f$\Theta\f$ function as a veto
   *  so that the emission is within the allowed phase space.
   * @param t  The scale
   * @return true if vetoed
   */
  bool PSVeto(const Energy2 t);

  /**
   * Compute the limits on \f$z\f$ for time-like branching
   * @param scale The scale of the particle
   * @return True if lower limit less than upper, otherwise false
   */
  bool computeTimeLikeLimits(Energy2 & scale);

  /**
   * Compute the limits on \f$z\f$ for space-like branching
   * @param scale The scale of the particle
   * @param x The energy fraction of the parton
   * @return True if lower limit less than upper, otherwise false
   */
  bool computeSpaceLikeLimits(Energy2 & scale, double x);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class with persistent data.
   */
  static NoPIOClassDescription<QTildeSudakov> initQTildeSudakov;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeSudakov & operator=(const QTildeSudakov &);

private:
  
  /**
   *  The evolution scale, \f$\tilde{q}\f$.
   */
  Energy q_;

  /**
   *  The Ids of the particles in the current branching
   */
  IdList ids_;

  /**
   *  The masses of the particles in the current branching
   */
  vector<Energy> masses_;

  /**
   *  The mass squared of the particles in the current branching
   */
  vector<Energy2> masssquared_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QTildeSudakov. */
template <>
struct BaseClassTrait<Herwig::QTildeSudakov,1> {
  /** Typedef of the first base class of QTildeSudakov. */
  typedef Herwig::SudakovFormFactor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QTildeSudakov class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QTildeSudakov>
  : public ClassTraitsBase<Herwig::QTildeSudakov> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QTildeSudakov"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QTildeSudakov is implemented. It may also include several, space-separated,
   * libraries if the class QTildeSudakov depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_QTildeSudakov_H */
