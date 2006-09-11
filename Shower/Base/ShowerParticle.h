// -*- C++ -*-
#ifndef HERWIG_ShowerParticle_H
#define HERWIG_ShowerParticle_H
//
// This is the declaration of the ShowerParticle class.
//

#include "ThePEG/EventRecord/Particle.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.fh"
#include "Herwig++/Shower/ShowerConfig.h"
#include "ShowerKinematics.h"
#include "Herwig++/Shower/Couplings/ShowerIndex.h"
#include "ShowerParticle.fh"

namespace Herwig {

using namespace ThePEG;

  /** \ingroup Shower
   *  This class represents a particle in the showering process.
   *  It inherits from the Particle class of ThePEG and has some  
   *  specifics information useful only during the showering process.
   * 
   *  Notice that: 
   *    - for forward evolution, it is clear what is meant by parent/child; 
   *      for backward evolution, however, it depends whether we want 
   *      to keep a physical picture or a Monte-Carlo effective one. 
   *      In the former case, an incoming particle (emitting particle)  
   *      splits into an emitted particle and the emitting particle after 
   *      the emission: the latter two are then children of the 
   *      emitting particle, the parent. In the Monte-Carlo effective  
   *      picture, we have that the particle close to the hard subprocess, 
   *      with higher (space-like) virtuality, splits into an emitted particle 
   *      and the emitting particle at lower virtuality: the latter two are, 
   *      in this case, the children of the first one, the parent. However we
   *      choose a more physical picture where the new emitting particle is the
   *      parented of the emitted final-state particle and the original emitting
   *      particle.
   *    - the pointer to a SplitFun object is set only in the case 
   *      that the particle has undergone a shower emission. This is similar to
   *      the case of the decay of a normal Particle where 
   *      the pointer to a Decayer object is set only in the case 
   *      that the particle has undergone to a decay. 
   *      In the case of particle connected directly to the hard subprocess, 
   *      there is no pointer to the hard subprocess, but there is a method 
   *      isFromHardSubprocess() which returns true only in this case.
   *
   *  @see Particle
   *  @see ShowerConfig
   *  @see ShowerIndex
   *  @see ShowerColourLine
   *  @see ShowerKinematics
   */
class ShowerParticle: public Particle {

public:

  /** @name Construction and descruction functions. */
  //@{
  /**
   * Particle uses the FixedSizeAllocator for (de)allocation.
   */
  inline void * operator new(size_t);
 
  /**
   * Particle uses the FixedSizeAllocator for (de)allocation.
   */
  inline void operator delete(void *, size_t);

  /**
   * Standard Constructor. Note that the default constructor is
   * private - there is no particle without a pointer to a
   * ParticleData object.
   */
  inline ShowerParticle(tcEventPDPtr);

  /**
   * Copy constructor from a ThePEG Particle
   * @param part ThePEG particle
   * @param pert Where the particle came from
   */
  inline ShowerParticle(const Particle & part,unsigned int pert);

  //@}

public:

  /**
   *   Access/Set various flags about the state of the particle
   */
  //@{
  /**
   * Access the flag that tells if the particle is final state
   * or initial state.
   */
  inline bool isFinalState() const;

  /**
   * Set the flag that tells if the particle is final state
   * or initial state.
   */
  inline void setFinalState( const bool );

  /**
   * Access the flag that tells if the particle is initiating a
   * time like shower when it has been emitted in an initial state shower.
   */
  inline bool initiatesTLS() const;

  /**
   * Set the flag that tells if the particle is initiating a
   * time like shower when it has been emitted in an initial state shower.
   */
  inline void setInitiatesTLS( const bool );

  /**
   * Access the flag which tells us where the particle came from
   * This is 0 for a particle produced in the shower, 1 if the particle came
   * from the hard sub-process and 2 is it came from a decay.
   */
  inline unsigned int perturbative() const;
  //@}

  /**
   * Access/Set Sudakov variables.
   * Notice that the ShowerKinematics object is logically
   * associated more with the branching vertex than with the radiating 
   * particle itself, although it is stored in the ShowerParticle 
   * object associated with the branching (radiating) particle. 
   * Furthermore, the branching products do not have a ShowerKinematics
   * object (they eventually will have one only later if they branch).
   * Therefore these Sudakov variables can be considered as a useful
   * representation of the temporarily, preliminary, momentum of the 
   * ShowerParticle object during the showering evolution, 
   * whereas the momentum member describes the "final", "real" one.
   */
  //@{
  /**
   *  Get the \f$\alpha\f$ parameter of the Sudakov decomposition
   */
  inline double sudAlpha() const;

  /**
   *  Set the \f$\alpha\f$ parameter of the Sudakov decomposition
   */
  inline void sudAlpha(const double);

  /**
   *  Get the \f$\beta\f$ parameter of the Sudakov decompostion
   */
  inline double sudBeta() const;

  /**
   *  Set the \f$\beta\f$ parameter of the Sudakov decompostion
   */
  inline void sudBeta(const double);

  /**
   *  Get the x-component of the transverse momentum in the Sudakov decompostion
   */
  inline Energy sudPx() const;

  /**
   *  Set the x-component of the transverse momentum in the Sudakov decompostion
   */
  inline void sudPx(const Energy);

  /**
   *  Get the y-component of the transverse momentum in the Sudakov decompostion
   */
  inline Energy sudPy() const;

  /**
   *  Get the y-component of the transverse momentum in the Sudakov decompostion
   */
  inline void sudPy(const Energy);

  /**
   *  Get the transverse momentum
   */
  inline Energy sudPperp() const;

  /**
   *  get the square of the transverse momentum
   */
  inline Energy2 sudPperp2() const;

  /**
   *  For an initial state particle get the fraction of the beam momentum
   */
  inline void x(double x);

  /**
   *  For an initial state particle set the fraction of the beam momentum
   */
  inline double x() const;
  //@}


  /**
   * Set/Get methods for the ShowerKinematics objects
   */
  //@{
  /**
   * Access/ the ShowerKinematics object.
   */
  inline const ShoKinPtr & showerKinematics() const;

  /**
   * Set the ShowerKinematics object.
   */
  void setShowerKinematics(const ShoKinPtr);
  //@}

  /**
   *  Members relating to the initial evolution scale and partner for the particle
   */
  //@{
  /**
   * Return (a const reference to) the vector of evolution scales
   * (\f$\tilde{q}\f$ scales)
   */
  inline const vector<Energy> & evolutionScales() const;

  /**
   *  Set the evolution \f$\tilde{q}\f$ scale for a given interaction type
   */
  inline void setEvolutionScale(const ShowerIndex::InteractionType, 
				const Energy);

  /** 
   * Return a vector of (pointers to) the partners corresponding to each 
   * considered interaction types (QCD, QED, EWK,...) defined in ShowerIndex. 
   * The vector of (pointers to) the partners is needed only as the
   * most general way to decide in which frame the shower is described.
   */
  inline const tShowerParticleVector & partners() const;

  /**
   * Set the partner for the specified interaction.
   */
  inline void setPartner(const ShowerIndex::InteractionType, 
			 const tShowerParticlePtr);
  //@}

  /**
   * Access/Set the flag that tells if the particle should be
   * treated in a special way during the kinematics reconstruction
   * (see KinematicsReconstructor class). 
   * In practice, it returns true when either the particle is childless, 
   * or is a on-shell decaying particle (in which case we have to set the flag to
   * true before the showering of this particle: it is not enough to check 
   * if decayer() is not null, because if it emits radiation
   * the decays products will be "transferred" to the particle
   * instance after the showering).
   */
  //@{
  /**
   *  Get the flag
   */
  inline bool isReconstructionFixedPoint() const;

  /**
   *  Set the flag
   */
  inline void setReconstructionFixedPoint(const bool);
  //@}

  /**
   *  If this particle came from the hard process get a pointer to ThePEG particle
   *  it came from
   */
  inline const tcPPtr getThePEGBase() const;

protected:

  /**
   * Standard clone function.
   */
  inline virtual PPtr clone() const;

  /**
   * Standard clone function.
   */
  inline virtual PPtr fullclone() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ShowerParticle> initShowerParticle;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerParticle & operator=(const ShowerParticle &);

private:

  /**
   *  Whether the particle is in the final or initial state
   */
  bool _isFinalState;

  /**
   *  Whether the particle is a reconstruction fixed point
   */
  bool _reconstructionFixedPoint;

  /**
   *  Whether the particle came from 
   */
  unsigned int _perturbative;

  /**
   *  Does a particle produced in the backward shower initiate a time-like shower 
   */
  bool _initiatesTLS;

  /**
   *  The \f$\alpha\f$ parameter of the Sudakov decomposition
   */
  double _sudAlpha;

  /**
   *  The \f$\beta\f$ parameter of the Sudakov decompostion
   */
  double _sudBeta;

  /**
   *  The x-component of the transverse momentum in the Sudakov decompostion
   */
  Energy _sudPx;

  /**
   *  The y-component of the transverse momentum in the Sudakov decompostion
   */
  Energy _sudPy;

  /**
   *  The beam energy fraction for particle's in the initial state
   */
  double _x;

  /**
   *  The shower kinematics for the particle
   */
  ShoKinPtr _showerKinematics;

  /**
   *  Evolution scales
   */
  vector<Energy> _scales;

  /**
   *  Partners
   */
  tShowerParticleVector _partners;

  /**
   *  Pointer to ThePEG Particle this ShowerParticle was created from
   */
  const tcPPtr _thePEGBase;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ShowerParticle. */
template <>
struct BaseClassTrait<Herwig::ShowerParticle,1> {
  /** Typedef of the first base class of ShowerParticle. */
  typedef Particle NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ShowerParticle class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ShowerParticle>
  : public ClassTraitsBase<Herwig::ShowerParticle> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::ShowerParticle"; }
  /** Create a Event object. */
  static TPtr create() { return TPtr::Create(Herwig::ShowerParticle(tcEventPDPtr())); }
};

/** @endcond */

}

#include "ShowerParticle.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerParticle.tcc"
#endif

#endif /* HERWIG_ShowerParticle_H */
