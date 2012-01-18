// -*- C++ -*-
//
// ShowerParticle.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerParticle_H
#define HERWIG_ShowerParticle_H
//
// This is the declaration of the ShowerParticle class.
//

#include "ThePEG/EventRecord/Particle.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.fh"
#include "Herwig++/Shower/ShowerConfig.h"
#include "ShowerKinematics.h"
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
   *  @see ShowerKinematics
   */
class ShowerParticle: public Particle {

public:

  /** @name Construction and descruction functions. */
  //@{

  /**
   * Standard Constructor. Note that the default constructor is
   * private - there is no particle without a pointer to a
   * ParticleData object.
   * @param x the ParticleData object
   * @param fs  Whether or not the particle is an inital or final-state particle
   * @param tls Whether or not the particle initiates a time-like shower
   */
  ShowerParticle(tcEventPDPtr x, bool fs, bool tls=false) 
    : Particle(x), _isFinalState(fs), _reconstructionFixedPoint( false ),
      _perturbative(0), _initiatesTLS(tls), _x(1.0), _showerKinematics(),
      _scale(ZERO), _vMass(ZERO), _thePEGBase(), _evolutionScale2(), _radiationLine() {}

  /**
   * Copy constructor from a ThePEG Particle
   * @param x ThePEG particle
   * @param pert Where the particle came from
   * @param fs Whether or not the particle is an inital or final-state particle
   * @param tls Whether or not the particle initiates a time-like shower
   */
  ShowerParticle(const Particle & x, unsigned int pert, bool fs, bool tls=false)
    : Particle(x), _isFinalState(fs), _reconstructionFixedPoint( false ),
    _perturbative(pert), _initiatesTLS(tls), _x(1.0), _showerKinematics(),
    _scale(ZERO), _vMass(ZERO), _thePEGBase(&x), _evolutionScale2(), _radiationLine() {}
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
  bool isFinalState() const { return _isFinalState; }

  /**
   * Access the flag that tells if the particle is initiating a
   * time like shower when it has been emitted in an initial state shower.
   */
  bool initiatesTLS() const { return _initiatesTLS; }

  /**
   * Access the flag which tells us where the particle came from
   * This is 0 for a particle produced in the shower, 1 if the particle came
   * from the hard sub-process and 2 is it came from a decay.
   */
  unsigned int perturbative() const { return _perturbative; }
  //@}

  /**
   * Set/Get the momentum fraction for initial-state particles
   */
  //@{
  /**
   *  For an initial state particle get the fraction of the beam momentum
   */
  void x(double x) { _x = x; }

  /**
   *  For an initial state particle set the fraction of the beam momentum
   */
  double x() const { return _x; }
  //@}

  /**
   * Set/Get methods for the ShowerKinematics objects
   */
  //@{
  /**
   * Access/ the ShowerKinematics object.
   */
  const ShoKinPtr & showerKinematics() const { return _showerKinematics; }


  /**
   * Set the ShowerKinematics object.
   */
  void setShowerKinematics(const ShoKinPtr in) { _showerKinematics = in; }
  //@}

  /**
   *  Members relating to the initial evolution scale and partner for the particle
   */
  //@{
  /**
   * Return the evolution scale \f$\tilde{q}\f$
   */
  Energy evolutionScale() const { return _scale; }



  /**
   *  Set the evolution \f$\tilde{q}\f$ scale
   */
  void setEvolutionScale(Energy scale) { _scale = scale; }

  /**
   * Return the virtual mass\f$
   */
  Energy virtualMass() const { return _vMass; }

  /**
   *  Set the virtual mass
   */
  void setVirtualMass(Energy mass) { _vMass = mass; }

  /** 
   * Return the partner
   */
  tShowerParticlePtr partner() const { return _partner; }


  /**
   * Set the partner
   */
  void setPartner(const tShowerParticlePtr partner) { _partner = partner; } 


  /**
   * Return the evolution scale \f$\tilde{q}\f$ belonging to the second partner
   */
  Energy evolutionScale2() const { return _evolutionScale2; }

  /**
   *  Set the evolution \f$\tilde{q}\f$ scale of the second partner for gluon
   */
  void setEvolutionScale2(Energy evolutionScale2) { _evolutionScale2 = evolutionScale2; }
  
  /**
   *  Return the radiation line of a gluon
   *  This is 0 for a particle with random radiation choice, 1 for the colour
   *  line and 2 for the anti-colour line.
   */
  int radiationLine() { return _radiationLine; }

  /**
   *  Set the radiation line of a gluon
   */   
  void setRadiationLine(int radiationLine) { _radiationLine = radiationLine; }
  
  
  /** 
   * Return the progenitor of the shower
   */
  tShowerParticlePtr progenitor() const { return _progenitor; }


  /**
   * Set the progenitor of the shower
   */
  void setProgenitor(const tShowerParticlePtr progenitor) { _progenitor = progenitor; } 
    

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
  bool isReconstructionFixedPoint() const { return _reconstructionFixedPoint || children().empty(); }

  /**
   *  Set the flag
   */
  void setReconstructionFixedPoint(const bool in) { _reconstructionFixedPoint = in; }
  //@}

  /**
   *  Members to store and provide access to variables for a specific
   *  shower evolution scheme
   */
  //@{
  /**
   *  Set the vector containing dimensionless variables
   */
  vector<double> & showerParameters() { return _parameters; }

  /**
   *  Set the vector containing dimensionful variables
   */
  vector<Energy> & showerVariables() { return _variables; }
  //@}

  /**
   *  If this particle came from the hard process get a pointer to ThePEG particle
   *  it came from
   */
  const tcPPtr getThePEGBase() const { return _thePEGBase; }

protected:

  /**
   * Standard clone function.
   */
  virtual PPtr clone() const;

  /**
   * Standard clone function.
   */
  virtual PPtr fullclone() const;

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
   * Dimensionless parameters
   */
  vector<double> _parameters;

  /**
   *  Dimensionful parameters
   */
  vector<Energy> _variables;

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
  Energy _scale;

  /**
   *  Virtual mass
   */
  Energy _vMass;

  /**
   *  Partners
   */
  tShowerParticlePtr _partner;

  /**
   *  Pointer to ThePEG Particle this ShowerParticle was created from
   */
  const tcPPtr _thePEGBase;
  
  /**
   *  Second evolution scale
   */  
  Energy _evolutionScale2; 
  
  /**
   *  Radiation Line
   */
  int _radiationLine;
  
  /**
   *  Progenitor
   */   
  tShowerParticlePtr _progenitor;
    
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
  static string className() { return "Herwig::ShowerParticle"; }
  /** Create a Event object. */
  static TPtr create() { return TPtr::Create(Herwig::ShowerParticle(tcEventPDPtr(),true)); }
};

/** @endcond */

}

#endif /* HERWIG_ShowerParticle_H */
