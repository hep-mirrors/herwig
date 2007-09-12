// -*- C++ -*-
#ifndef HERWIG_ForcedSplitting_H
#define HERWIG_ForcedSplitting_H
//
// This is the declaration of the ForcedSplitting class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/EventRecord/RemnantParticle.h" 
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/PDT/RemnantDecayer.h"

#include "ForcedSplitting.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Hadronization
 *
 *  This is the definition of a simple forced splitting algorithm.
 *  This takes the Remnant object produced from the PDF and backward
 *  evolution (hadron - parton) and produce partons with the remaining 
 *  flavours and with the correct colour connections.
 *
 *  The algorithim operates by starting with the parton which enters the hard process.
 *  If this is from the sea there is a forced branching to produce the antiparticle
 *  from a gluon branching. If the parton entering the hard process was a gluon, or
 *  a gluon was produced from the first step of the algorithm, there is then a further
 *  branching back to a valence parton. After these partons have been produced a quark or
 *  diquark is produced to give the remaining valence content of the incoming hadron.
 *
 *  The forced branching are generated using a scale between QSpac and EmissionRange times
 *  the minimum scale. The energy fractions are then distributed using
 *  \f[\frac{\alpha_S}{2\pi}\frac{P(z)}{z}f(x/z,\tilde{q})\f]
 *  with the massless splitting functions.
 *
 * @see \ref ForcedSplittingInterfaces "The interfaces"
 * defined for ForcedSplitting.
 */
class ForcedSplitting: public Interfaced {

public:

  /**
   * The default constructor.
   */
  inline ForcedSplitting();

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
   * This takes the particle and find a splitting for np -> p + child and 
   * creates the correct kinematics and connects for such a split. This
   * Splitting has an upper bound on qtilde given by the energy argument
   * @param rem The Remnant
   * @param child The PDG code for the outgoing particle
   * @param oldQ  The maximum scale for the evolution
   * @param oldx  The fraction of the hadron's momentum carried by the last parton
   * @param pf    The momentum of the last parton at input and after branching at output
   * @param p     The total emitted momentum
   * @param iopt  Whether to use the \f$q\to gq\f$ or \f$g\to q\bar{q}\f$ splitting function.
   * @param step The step into which the new particles are inserted
   * @param first true if split belongs to extraction of partons from the hard process.
   */
  PPtr forceSplit(const tRemPPtr rem, long child, Energy &oldQ, double &oldx, 
		  Lorentz5Momentum &pf, Lorentz5Momentum &p,const unsigned int iopt,
		  const tStepPtr step, const bool first) const;

  /**
   * This computes the momentum of the emitted parton. 
   * @param par Momentum of the beam particle 
   * @param lastQ The maximum scale for the branching
   * @param lastx \f$x\f$ after the last emission
   * @param parton The parton in the last emission
   * @param pf at input is the momentum of the last parton and the emitted momentum on output
   * @param iopt  Whether to use the \f$q\to gq\f$ or \f$g\to q\bar{q}\f$ splitting function.
   * @param first true if split belongs to extraction of partons from the hard process.
   */
  Lorentz5Momentum emit(const Lorentz5Momentum &par, Energy &lastQ, 
			double &lastx, PPtr parton, Lorentz5Momentum &pf,
			const unsigned int iopt, const bool first) const;

  /**
   * This creates a parton from the remaining flavours of the hadron. The
   * last parton used was a valance parton, so only 2 (or 1, if meson) flavours
   * remain to be used.
   */
  PPtr finalSplit(const tRemPPtr rem, long remID, Lorentz5Momentum,
		  const tStepPtr ) const;

  /**
   * Set the beam particle.
   */
  inline void setBeam(tcPPtr beam) const;

  /**
   * Set the ThePEG::BeamParticleData object of the current hadron.
   */
  inline Energy getQspac() const;

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ForcedSplitting> initForcedSplitting;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ForcedSplitting & operator=(const ForcedSplitting &);
  
private:

  /**
   *  The kinematic cut-off
   */
   Energy _kinCutoff;

  /**
   *  Range for emission
   */
  double _range;

  /**
   *  Start of the evolution
   */
  Energy _qspac;

  /**
   *  Size of the bins in z for the interpolation
   */
  double _zbin;

  /**
   *  Size of the bins in y for the interpolation
   */
  double _ybin;

  /**
   *  Maximum number of bins for the z interpolation
   */
  int _nbinmax;

  /**
   *  Pointer to the object calculating the QCD coupling
   */
  ShowerAlphaPtr _alpha;

  /**
   *  The beam particle data for the current initial-state shower
   */
  mutable tcPPtr _beam;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ForcedSplitting. */
template <>
struct BaseClassTrait<Herwig::ForcedSplitting,1> {
  /** Typedef of the first base class of ForcedSplitting. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ForcedSplitting class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ForcedSplitting>
  : public ClassTraitsBase<Herwig::ForcedSplitting> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ForcedSplitting"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ForcedSplitting is implemented. It may also include several, space-separated,
   * libraries if the class ForcedSplitting depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwHadronization.so"; }
};

/** @endcond */

}

#include "ForcedSplitting.icc"

#endif /* HERWIG_ForcedSplitting_H */
