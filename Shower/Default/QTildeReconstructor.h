// -*- C++ -*-
//
// QTildeReconstructor.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_QTildeReconstructor_H
#define HERWIG_QTildeReconstructor_H
//
// This is the declaration of the QTildeReconstructor class.
//

#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "ThePEG/Repository/UseRandom.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  A simple struct to store the information we need on the 
 *  showering
 */
struct JetKinStruct {

  /**
   *  Parent particle of the jet
   */
  tShowerParticlePtr parent;

  /**
   *  Momentum of the particle before reconstruction
   */
  Lorentz5Momentum p;

  /**
   *  Momentum of the particle after reconstruction
   */  
  Lorentz5Momentum q;
};

/**
 * typedef for a vector of JetKinStruct
 */  
typedef vector<JetKinStruct> JetKinVect;

/** \ingroup Shower
 *
 * This class is responsible for the kinematical reconstruction 
 * after each showering step, and also for the necessary Lorentz boosts 
 * in order to preserve energy-momentum conservation in the overall collision,
 * and also the invariant mass and the rapidity of the hard subprocess system.
 * In the case of multi-step showering, there will be not unnecessary
 * kinematical reconstructions. 
 *
 * There is also the option of taking a set of momenta for the particles
 * and inverting the reconstruction to give the evolution variables for the
 * shower.
 *
 * Notice:
 * - although we often use the term "jet" in either methods or variables names,
 *   or in comments, which could appear applicable only for QCD showering,
 *   there is indeed no "dynamics" represented in this class: only kinematics 
 *   is involved, as the name of this class remainds. Therefore it can be used
 *   for any kind of showers (QCD-,QED-,EWK-,... bremsstrahlung).
 * 
 * @see ShowerParticle
 * @see ShowerKinematics
 * @see \ref QTildeReconstructorInterfaces "The interfaces"
 * defined for QTildeReconstructor.
 */
class QTildeReconstructor: public KinematicsReconstructor {

public:

  /**
   *  Default constructor
   */
  QTildeReconstructor() : _reconopt(0), _minQ(0.001*GeV) {};

  /**
   *  Methods to reconstruct the kinematics of a scattering or decay process
   */
  //@{
  /**
   * Given in input a vector of the particles which initiated the showers
   * the method does the reconstruction of such jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructHardJets(ShowerTreePtr hard,
				   const map<tShowerProgenitorPtr,
				   pair<Energy,double> > & pt) const;

  /**
   * Given in input a vector of the particles which initiated the showers
   * the method does the reconstruction of such jets,
   * including the appropriate boosts (kinematics reshufflings)  
   * needed to conserve the total energy-momentum of the collision
   * and preserving the invariant mass and the rapidity of the 
   * hard subprocess system.
   */
  virtual bool reconstructDecayJets(ShowerTreePtr decay) const;
  //@}

  /**
   *  Methods to invert the reconstruction of the shower for
   *  a scattering or decay process and calculate
   *  the variables used to generate the
   *  shower given the particles produced.
   *  This is needed for the CKKW and POWHEG approaches
   */
  //@{
  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the 
   * shower
   */
  virtual bool deconstructDecayJets(HardTreePtr, EvolverPtr,
				    ShowerInteraction::Type) const;

  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the shower
   *  for a hard process
   */
  virtual bool deconstructHardJets(HardTreePtr, EvolverPtr,
				   ShowerInteraction::Type) const;
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
   *  Methods to reconstruct the kinematics of individual jets
   */
  //@{
  /**
   * Given the particle (ShowerParticle object) that 
   * originates a forward (time-like) jet, this method reconstructs the kinematics 
   * of the jet. That is, by starting from the final grand-children (which 
   * originates directly or indirectly from particleJetParent, 
   * and which don't have children), and moving "backwards" (in a physical
   * time picture), towards the particleJetParent, the 
   * ShowerKinematics objects associated with the various particles, 
   * which have been created during the showering, are now completed. 
   * In particular, at the end, we get the mass of the jet, which is the 
   * main information we want.
   * This methods returns false if there was no radiation or rescaling required
   */
  virtual bool reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent,
				      unsigned int iopt) const;

  /**
   * Exactly similar to the previous one, but for a space-like jet.
   * Also in this case we start from the final grand-children (which
   * are childless) of the particle which originates the jet, but in
   * this case we proceed "forward" (in the physical time picture)
   * towards the particleJetParent.
   * This methods returns false if there was no radiation or rescaling required
   */
  bool reconstructSpaceLikeJet(const tShowerParticlePtr particleJetParent) const;

  /**
   * Exactly similar to the previous one, but for a decay jet
   * This methods returns false if there was no radiation or rescaling required
   */
  bool reconstructDecayJet(const tShowerParticlePtr particleJetParent) const;
  //@}

  /**
   *  Methods to perform the reconstruction of various types of colour
   *  singlet systems
   */
  //@{
  /**
   *  Perform the reconstruction of a system with one incoming and at least one
   *  outgoing particle
   */
  void reconstructInitialFinalSystem(vector<ShowerProgenitorPtr>) const;

  /**
   *  Perform the reconstruction of a system with only final-state
   *  particles
   */
  void reconstructFinalStateSystem(bool applyBoost, Boost toRest, Boost fromRest, 
				   vector<ShowerProgenitorPtr>) const;

  /**
   *  Reconstruction of a general coloured system
   */
  void reconstructGeneralSystem(vector<ShowerProgenitorPtr> & ShowerHardJets) const;

  /**
   *  Perform the reconstruction of a system with only final-state
   *  particles
   */
  void reconstructInitialInitialSystem(bool & applyBoost,
				       Boost & toRest, Boost & fromRest,
				       vector<ShowerProgenitorPtr>) const;
  //@}

  /**
   *  Methods to perform the inverse reconstruction of various types of
   *  colour singlet systems
   */
  //@{
  /**
   *  Perform the inverse reconstruction of a system with only final-state
   *  particles
   */
  void deconstructFinalStateSystem(Boost & toRest, Boost & fromRest,
				   HardTreePtr,
				   vector<HardBranchingPtr>,
				   EvolverPtr,
				   ShowerInteraction::Type) const;
  
  /**
   *  Perform the inverse reconstruction of a system with only initial-state
   *  particles
   */
  void deconstructInitialInitialSystem(bool & applyBoost,
				       Boost & toRest, Boost & fromRest,
				       HardTreePtr,
				       vector<HardBranchingPtr>,
				       ShowerInteraction::Type ) const;

  /**
   *  Perform the inverse reconstruction of a system with only initial-state
   *  particles
   */
  void deconstructInitialFinalSystem(HardTreePtr,
				     vector<HardBranchingPtr>,
				     EvolverPtr,
				     ShowerInteraction::Type ) const;

  bool deconstructGeneralSystem(HardTreePtr, EvolverPtr,
				ShowerInteraction::Type) const;
  //@}

  /**
   *  Various methods for the Lorentz transforms needed to do the 
   *  rescalings
   */
  //@{
  /**
   * Compute the boost to get from the the old momentum to the new 
   */
  LorentzRotation solveBoost(const double k, const Lorentz5Momentum & newq, 
			     const Lorentz5Momentum & oldp) const;
  
  /**
   * Compute the boost to get from the the old momentum to the new 
   */
  LorentzRotation solveBoost(const Lorentz5Momentum & newq, 
			     const Lorentz5Momentum & oldq) const;
  
  /**
   * Compute the boost to get from the the old momentum to the new 
   */
  LorentzRotation solveBoostZ(const Lorentz5Momentum & newq, 
			      const Lorentz5Momentum & oldq) const;
  
  /**
   *  Recursively boost the initial-state shower
   * @param p The particle
   * @param bv The boost
   * @param parent The parent of the chain
   */
  void boostChain(tPPtr p, const LorentzRotation & bv, tPPtr & parent) const;

  /**
   * Given a 5-momentum and a scale factor, the method returns the
   * Lorentz boost that transforms the 3-vector vec{momentum} --->
   * k*vec{momentum}. The method returns the null boost in the case no
   * solution exists. This will only work in the case where the
   * outgoing jet-momenta are parallel to the momenta of the particles
   * leaving the hard subprocess. 
   */
  Boost solveBoostBeta( const double k, const Lorentz5Momentum & newq, 
			  const Lorentz5Momentum & oldp);

  /**
   * Compute boost parameter along z axis to get (Ep, any perp, qp)
   * from (E, same perp, q).
   */
  double getBeta(const double E, const double q, 
		 const double Ep, const double qp) const
  {return (q*E-qp*Ep)/(sqr(qp)+sqr(E));}
  //@}

  /**
   *  Methods to calculate the various scaling factors
   */
  //@{
  /**
   * Given a vector of 5-momenta of jets, where the 3-momenta are the initial
   * ones before showering and the masses are reconstructed after the showering,
   * this method returns the overall scaling factor for the 3-momenta of the
   * vector of particles, vec{P}_i -> k * vec{P}_i, such to preserve energy-
   * momentum conservation, i.e. after the rescaling the center of mass 5-momentum 
   * is equal to the one specified in input, cmMomentum. 
   * The method returns 0 if such factor cannot be found.
   * @param root_s Centre-of-mass energy
   * @param jets The jets
   */
  double solveKfactor( const Energy & root_s, const JetKinVect & jets ) const;

  /**
   *  Calculate the rescaling factors for the jets in a particle decay where
   *  there was initial-state radiation
   * @param mb The mass of the decaying particle
   * @param n  The reference vector for the initial state radiation
   * @param pjet The momentum of the initial-state jet
   * @param jetKinematics The JetKinStruct objects for the jets
   * @param partner The colour partner
   * @param ppartner The momentum of the colour partner of the decaying particle
   * before and after radiation
   * @param k1 The rescaling parameter for the partner
   * @param k2 The rescaling parameter for the outgoing singlet
   * @param qt The transverse momentum vector
   */
  bool solveDecayKFactor(Energy mb, 
			 const Lorentz5Momentum & n, 
			 const Lorentz5Momentum & pjet, 
			 const JetKinVect & jetKinematics, 
			 ShowerParticlePtr partner,
			 Lorentz5Momentum ppartner[2],
			 double & k1, 
			 double & k2,
			 Lorentz5Momentum & qt) const;

  /**
   * Compute the momentum rescaling factor needed to invert the shower
   * @param pout The momenta of the outgoing particles
   * @param mon  The on-shell masses
   * @param roots The mass of the decaying particle
   */
  double inverseRescalingFactor(vector<Lorentz5Momentum> pout,
				vector<Energy> mon,Energy roots) const;

  /**
   * Compute the momentum rescaling factor needed to invert the shower
   * @param pout The momenta of the outgoing particles
   * @param mon  The on-shell masses
   * @param roots The mass of the decaying particle
   */
  bool inverseDecayRescalingFactor(vector<Lorentz5Momentum> pout,
				   vector<Energy> mon,Energy roots,
				   Lorentz5Momentum ppartner, Energy mbar,
				   double & k1, double & k2) const;

  /**
   * Check the rescaling conserves momentum
   * @param k The rescaling
   * @param root_s The centre-of-mass energy
   * @param jets The jets
   */
  Energy momConsEq(const double & k, const Energy & root_s,
			  const JetKinVect & jets) const;
  //@}

  /**
   *  Find the colour partners of a particle to identify the colour singlet
   *  systems for the reconstruction.
   */
  vector<unsigned int> findPartners(unsigned int ,vector<ShowerProgenitorPtr>) const;

  /**
   *  Find the colour partners for as branching  to identify the colour singlet
   *  systems for the inverse reconstruction.
   */
  void findPartners(HardBranchingPtr branch,set<HardBranchingPtr> & done,
		    const set<HardBranchingPtr> & branchings,
		    vector<HardBranchingPtr> & jets) const;

  /**
   *  Add the intrinsic \f$p_T\f$ to the system if needed
   */
  bool addIntrinsicPt(vector<ShowerProgenitorPtr>) const;

  void deepTransform(PPtr particle,const LorentzRotation & r,
		     bool match=true,PPtr original=PPtr()) const;

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static ClassDescription<QTildeReconstructor> initQTildeReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeReconstructor & operator=(const QTildeReconstructor &);

private:

  /**
   *  Option for handling the reconstruction
   */
  unsigned int _reconopt;

  /**
   * Minimum invariant mass for initial-final dipoles to allow the
   * reconstruction
   */
  Energy _minQ;

  /**
   *  The progenitor of the jet currently being reconstructed
   */
  mutable tShowerParticlePtr _progenitor;

  /**
   * Storage of the intrinsic \f$p_T\f$
   */
  mutable map<tShowerProgenitorPtr,pair<Energy,double> > _intrinsic;

  /**
   *  Current ShowerTree
   */
  mutable tShowerTreePtr _currentTree;

  /**
   * Particles which shouldn't have their masses rescaled as
   * vector for the interface
   */
  PDVector _noRescaleVector;

  /**
   * Particles which shouldn't have their masses rescaled as
   * set for quick access
   */
  set<cPDPtr> _noRescale;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QTildeReconstructor. */
template <>
struct BaseClassTrait<Herwig::QTildeReconstructor,1> {
  /** Typedef of the first base class of QTildeReconstructor. */
  typedef Herwig::KinematicsReconstructor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QTildeReconstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QTildeReconstructor>
  : public ClassTraitsBase<Herwig::QTildeReconstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QTildeReconstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QTildeReconstructor is implemented. It may also include several, space-separated,
   * libraries if the class QTildeReconstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_QTildeReconstructor_H */
