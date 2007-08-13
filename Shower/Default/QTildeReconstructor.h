// -*- C++ -*-
#ifndef HERWIG_QTildeReconstructor_H
#define HERWIG_QTildeReconstructor_H
//
// This is the declaration of the QTildeReconstructor class.
//

#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "ThePEG/Repository/UseRandom.h"
#include "QTildeReconstructor.fh"
#include <cassert>

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
				   map<tShowerProgenitorPtr,
				   pair<Energy,double> > pt) const;

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
   *  Methods to reconstruct the variables used to generate the
   *  shower given the particles produced.
   *  This is needed for the CKKW and Nason approaches
   */
  //@{
  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the 
   * shower
   */
  virtual bool reconstructDecayShower(NasonTreePtr decay,EvolverPtr) const;

  /**
   *  Given the particles, with a history which we wish to interpret
   *  as a shower reconstruct the variables used to generate the shower
   *  for a hard process
   */
  virtual bool reconstructHardShower(NasonTreePtr hard,EvolverPtr) const;
  //@}

protected:

  /**
   *  Reconstruct the initial state jets
   */
  virtual bool reconstructISJets(Lorentz5Momentum pcm,
				 const vector<ShowerProgenitorPtr> & ShowerHardJets,
				 map<tShowerProgenitorPtr,pair<Energy,double> > intrinsic,
				 Boost & boostRest, Boost & boostNewF) const;

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

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

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
  const double solveKfactor( const Energy & root_s, const JetKinVect & jets ) const;

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
  bool solveDecayKFactor(Energy mb, Lorentz5Momentum n, Lorentz5Momentum pjet, 
			 const JetKinVect & jetKinematics, ShowerParticlePtr partner,
			 Lorentz5Momentum ppartner[2],
			 double & k1, double & k2,Lorentz5Momentum & qt) const;


  /**
   * Compute the momentum rescaling factor needed to invert the shower
   * @param pout The momenta of the outgoing particles
   * @param mon  The on-shell masses
   * @param roots The mass of the decaying particle
   */
  double inverseRescaleingFactor(vector<Lorentz5Momentum> pout,
				 vector<Energy> mon,Energy roots) const;

  /**
   * Check the rescaling conserves momentum
   * @param k The rescaling
   * @param root_s The centre-of-mass energy
   * @param jets The jets
   */
  inline Energy momConsEq(const double & k, const Energy & root_s,
			  const JetKinVect & jets) const;

  /**
   * Compute the boost to get from the the old momentum to the new 
   */
  inline LorentzRotation solveBoost(const double k, const Lorentz5Momentum & newq, 
				    const Lorentz5Momentum & oldp) const;

  /**
   * Compute the boost to get from the the old momentum to the new 
   */
  inline LorentzRotation solveBoost(const Lorentz5Momentum & newq, 
				    const Lorentz5Momentum & oldq) const;

  /**
   *  Recursively boost the initial-state shower
   * @param p The particle
   * @param bv The boost
   * @param parent The parent of the chain
   */
  inline void boostChain(tPPtr p, const Boost &bv, tPPtr & parent) const;

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
  inline double getBeta(const double E, const double q, 
			const double Ep, const double qp) const {
    return (q*E-qp*Ep)/(sqr(qp)+sqr(E));
  }

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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<QTildeReconstructor> initQTildeReconstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeReconstructor & operator=(const QTildeReconstructor &);
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

#include "QTildeReconstructor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QTildeReconstructor.tcc"
#endif

#endif /* HERWIG_QTildeReconstructor_H */
