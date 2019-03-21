// -*- C++ -*-
//
// ShowerParticle.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerParticle_H
#define HERWIG_ShowerParticle_H
//
// This is the declaration of the ShowerParticle class.
//

#include "ThePEG/EventRecord/Particle.h"
#include "Herwig/Shower/Core/SplittingFunctions/SplittingFunction.fh"
#include "Herwig/Shower/Core/ShowerConfig.h"
#include "ShowerBasis.h"
#include "ShowerKinematics.h"
#include "ShowerParticle.fh"
#include <iosfwd>

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

  /**
   *  Struct for all the info on an evolution partner
   */
  struct EvolutionPartner {

    /**
     *  Constructor
     */
    EvolutionPartner(tShowerParticlePtr p,double w, ShowerPartnerType t,
		     Energy s) : partner(p), weight(w), type(t), scale(s)
    {}

    /**
     * The partner
     */
    tShowerParticlePtr partner;
    
    /**
     *  Weight
     */
    double weight;

    /**
     *  Type
     */
    ShowerPartnerType type;

    /**
     *  The assoicated evolution scale
     */
    Energy scale;
  };

  /**
   *  Struct to store the evolution scales
   */
  struct EvolutionScales {

    /**
     *  Constructor
     */
    EvolutionScales() : QED(),QCD_c(),QCD_ac(),
			QED_noAO(),QCD_c_noAO(),QCD_ac_noAO(),
			Max_Q2(Constants::MaxEnergy2)
    {}

    /**
     *  QED scale
     */
    Energy QED;

    /**
     * QCD colour scale
     */
    Energy QCD_c;

    /**
     *  QCD anticolour scale
     */
    Energy QCD_ac;

    /**
     *  QED scale
     */
    Energy QED_noAO;

    /**
     * QCD colour scale
     */
    Energy QCD_c_noAO;

    /**
     *  QCD anticolour scale
     */
    Energy QCD_ac_noAO;

    /**
     *  Maximum allowed virtuality of the particle
     */
    Energy2 Max_Q2;
  };


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
    : Particle(x), _isFinalState(fs),
      _perturbative(0), _initiatesTLS(tls), _x(1.0), _showerKinematics(),
      _vMass(ZERO), _thePEGBase() {}

  /**
   * Copy constructor from a ThePEG Particle
   * @param x ThePEG particle
   * @param pert Where the particle came from
   * @param fs Whether or not the particle is an inital or final-state particle
   * @param tls Whether or not the particle initiates a time-like shower
   */
  ShowerParticle(const Particle & x, unsigned int pert, bool fs, bool tls=false)
    : Particle(x), _isFinalState(fs),
    _perturbative(pert), _initiatesTLS(tls), _x(1.0), _showerKinematics(),
    _vMass(ZERO), _thePEGBase(&x) {}
  //@}

public:

  /**
   *  Set a preliminary momentum for the particle
   */
  void setShowerMomentum(bool timelike);

  /**
   *  Construct the spin info object for a shower particle
   */
  void constructSpinInfo(bool timelike);

  /**
   * Perform any initial calculations needed after the branching has been selected
   */
  void initializeDecay();

  /**
   * Perform any initial calculations needed after the branching has been selected
   * @param parent The beam particle
   */
  void initializeInitialState(PPtr parent);

  /**
   * Perform any initial calculations needed after the branching has been selected
   */
  void initializeFinalState();

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
  void showerKinematics(const ShoKinPtr in) { _showerKinematics = in; }
  //@}

  /**
   * Set/Get methods for the ShowerBasis objects
   */
  //@{
  /**
   * Access/ the ShowerBasis object.
   */
  const ShowerBasisPtr & showerBasis() const { return _showerBasis; }


  /**
   * Set the ShowerBasis object.
   */
  void showerBasis(const ShowerBasisPtr in, bool copy) {
    if(!copy) 
      _showerBasis = in;
    else {
      _showerBasis = new_ptr(ShowerBasis());
      _showerBasis->setBasis(in->pVector(),in->nVector(),in->frame());
    } 
  }
  //@}

  /**    
   *  Members relating to the initial evolution scale and partner for the particle
   */
  //@{
  /**
   *  Veto emission at a given scale 
   */
  void vetoEmission(ShowerPartnerType type, Energy scale);

  /**
   *  Access to the evolution scales
   */
  const EvolutionScales & scales() const {return scales_;} 

  /**
   *  Access to the evolution scales
   */
  EvolutionScales & scales() {return scales_;} 

  /**
   * Return the virtual mass\f$
   */
  Energy virtualMass() const { return _vMass; }

  /**
   *  Set the virtual mass
   */
  void virtualMass(Energy mass) { _vMass = mass; }

  /** 
   * Return the partner
   */
  tShowerParticlePtr partner() const { return _partner; }

  /**
   * Set the partner
   */
  void partner(const tShowerParticlePtr partner) { _partner = partner; } 

  /**
   *  Get the possible partners 
   */
  vector<EvolutionPartner> & partners() { return partners_; }

  /**
   *  Add a possible partners 
   */
  void addPartner(EvolutionPartner in );

  /**
   *  Clear the evolution partners
   */
  void clearPartners() { partners_.clear(); }
    
  /** 
   * Return the progenitor of the shower
   */
  tShowerParticlePtr progenitor() const { return _progenitor; }

  /**
   * Set the progenitor of the shower
   */
  void progenitor(const tShowerParticlePtr progenitor) { _progenitor = progenitor; } 
  //@}


  /**
   *  Members to store and provide access to variables for a specific
   *  shower evolution scheme
   */
  //@{
  struct Parameters {
    Parameters() : alpha(1.), beta(), ptx(), pty(), pt() {}
    double alpha;
    double beta;
    Energy ptx;
    Energy pty;
    Energy pt;
  };


  /**
   *  Set the vector containing dimensionless variables
   */
  Parameters & showerParameters() { return _parameters; }
  //@}

  /**
   *  If this particle came from the hard process get a pointer to ThePEG particle
   *  it came from
   */
  const tcPPtr thePEGBase() const { return _thePEGBase; }
 
public:

  /**
   *  Extract the rho matrix including mapping needed in the shower
   */
  RhoDMatrix extractRhoMatrix(bool forward);

protected:

  /**
   * For a particle which came from the hard process get the spin density and
   * the mapping required to the basis used in the Shower
   * @param rho The \f$\rho\f$ matrix
   * @param mapping The mapping
   * @param showerkin The ShowerKinematics object
   */
  bool getMapping(SpinPtr &, RhoDMatrix & map);

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
  ShowerParticle & operator=(const ShowerParticle &) = delete;

private:

  /**
   *  Whether the particle is in the final or initial state
   */
  bool _isFinalState;

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
  Parameters _parameters;

  /**
   *  The beam energy fraction for particle's in the initial state
   */
  double _x;

  /**
   *  The shower kinematics for the particle
   */
  ShoKinPtr _showerKinematics;

  /**
   *  The shower basis for the particle
   */
  ShowerBasisPtr _showerBasis;

  /**
   *  Storage of the evolution scales
   */
  EvolutionScales scales_;

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
   *  Progenitor
   */   
  tShowerParticlePtr _progenitor;

  /**
   *  Partners
   */
  vector<EvolutionPartner> partners_;
    
};

inline ostream & operator<<(ostream & os, const ShowerParticle::EvolutionScales & es) {
  os << "Scales: QED=" << es.QED / GeV
     << " QCD_c=" << es.QCD_c / GeV
     << " QCD_ac=" << es.QCD_ac / GeV
     << " QED_noAO=" << es.QED_noAO / GeV
     << " QCD_c_noAO=" << es.QCD_c_noAO / GeV
     << " QCD_ac_noAO=" << es.QCD_ac_noAO / GeV
     << '\n';
  return os;
}

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
