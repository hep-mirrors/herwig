// -*- C++ -*-
#ifndef HERWIG_KinematicsReconstructor_H
#define HERWIG_KinematicsReconstructor_H
//
// This is the declaration of the <!id>KinematicsReconstructor<!!id> class.
//
// CLASSDOC SUBSECTION Description:
// 
// This class is responsible for the kinematical reconstruction <BR>
// after each showering step, and also for the necessary Lorentz boosts <BR> 
// in order to preserve energy-momentum conservation in the overall collision, <BR>
// and also the invariant mass and the rapidity of the hard subprocess system. <BR>
// In the case of multi-step showering, there will be not unnecessary <BR>
// kinematical reconstructions. <BR> 
//
// Notice: <BR>
// <UL>
//   <LI> although we often use the term "jet" in either methods or variables names, <BR>
//        or in comments, which could appear applicable only for QCD showering, <BR>
//        there is indeed no "dynamics" represented in this class: only kinematics <BR>
//        is involved, as the name of this class remainds. Therefore it can be used <BR>
//        for any kind of showers (QCD-,QED-,EWK-,... bremsstrahlung).
// </UL> 
// 
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerParticle.html">ShowerParticle.h</a>, <BR>
// <a href="http:ShowerKinematics.html">ShowerKinematics.h</a>.
// 

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/MatrixElement/MEBase.h"


namespace Herwig {

using namespace ThePEG;

class Lorentz5Vector;

class KinematicsReconstructor: public ThePEG::HandlerBase {

public:

  typedef vector<Lorentz5Momentum*> VecMomentaPtr;
  typedef vector<Lorentz5Momentum> VecMomenta;
  typedef vector<const Lorentz5Momentum*> CVecMomentaPtr;

  typedef map<tShowerParticlePtr, bool> MapShower;
  // For a given (pointer to) shower particle, which is the parent particle 
  // of a (forward, or backward, or decaying) jet, the flag tells whether
  // or not such jet needs to be reconstructed, that is whether some 
  // radiation has been emitted or if some of the "leaves", childless
  // gluons have been forced on their effective mass shell. 
  // It is useful only for multi-scale showering, to avoid to repeat
  // unnecessary kinematics reconstructions.
  // You can distinguish between time-like, space-like, or decaying jets
  // by using the properties of the pointed <!class>ShowerParticle<!!class> object.

  inline KinematicsReconstructor();
  inline KinematicsReconstructor(const KinematicsReconstructor &);
  virtual ~KinematicsReconstructor();
  // Standard ctors and dtor.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

  bool reconstructHardJets( const MapShower & mapShowerHardJets, 
			    const Lorentz5Momentum & pBeamHadron1,
			    const Lorentz5Momentum & pBeamHadron2,
                            const tcMEPtr specialHardProcess = tcMEPtr() ) 
    throw (Veto, Stop, Exception);
  // Given in input a map of pointers to hard particles (<!class>ShowerParticle<!!class>
  // objects), that is of the particles which enters (as incoming or
  // outcoming) the hard subprocess, each associated with a flag (bool)
  // which tells whether or not the associated jet needs a kinematical 
  // reconstruction, and the momenta (although we need only the directions)
  // of the beam hadrons, the method does the reconstruction of such jets,
  // including the appropriate boosts (kinematics reshufflings)  
  // needed to conserve the total energy-momentum of the collision
  // and preserving the invariant mass and the rapidity of the 
  // hard subprocess system.
  // Notice that these flags are useful only for multi-scale showering
  // in order to avoid unnecessary kinematical reconstruction.
  // The last argument, if not equal to the default null value,
  // specifies a special hard subprocess (indeed it is a pointer
  // to the corresponding matrix element, but this is used only
  // as "id" of that process, whereas the matrix element itself 
  // is never used here), like Deep Inelastic Scattering, that needs 
  // a special treatment (at the moment, only D.I.S. needs this 
  // special treatment; but, to be more general, we allow the possibility 
  // to have also other special processes, which could be treated 
  // similarly as D.I.S. or even differently).
  // If something goes wrong, the method throws an Exception.

  void reconstructDecayJets( const MapShower & mapShowerDecayJets )
    throw (Veto, Stop, Exception);
  // Given in input a map associated with the showering of a decay, 
  // that is made of the decaying particle and its decay products 
  // (indeed pointers to <!class>ShowerParticle<!!class> objects 
  //  associated with such particles),   
  // each associated with a flag (bool) which tells whether or not 
  // the corresponding jet needs a kinematical reconstruction, 
  // the method does the reconstruction of such jets, including
  // the appropriate boosts (kinematics reshuffling) needed
  // to conserve the energy-momentum of the showering decay, 
  // If something goes wrong, the method throw an Exception.
 
protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<KinematicsReconstructor> initKinematicsReconstructor;
  // Describe a concrete class with persistent data.

  KinematicsReconstructor & operator=(const KinematicsReconstructor &);
  // Private and non-existent assignment operator.

  typedef struct {
    Lorentz5Momentum p, q; 
    tShowerParticlePtr parent; 
  } JetKinStruct;

  typedef vector<JetKinStruct> JetKinVect;

  bool reconstructTimeLikeJet( const tShowerParticlePtr particleJetParent );
  // Given the particle (<!class>ShowerParticle<!!class> object) that 
  // originates a forward (time-like) jet, this method reconstructs the kinematics 
  // of the jet. That is, by starting from the final grand-children (which 
  // originates directly or indirectly from <!id>particleJetParent<!!id>, 
  // and which don't have children), and moving "backwards" (in a physical
  // time picture), towards the <!id>particleJetParent<!!id>, the 
  // <!class>ShowerKinematics<!!class> objects associated with the various particles, 
  // which have been created during the showering, are now completed. 
  // In particular, at the end, we get the mass of the jet, which is the 
  // main information we want.
  // The method should always returns true, because this reconstructor
  // is physically always possible; however, because of possible 
  // programming bugs, we let the method return false if something 
  // goes wrong.

  bool reconstructSpaceLikeJet( const tShowerParticlePtr particleJetParent );
  // Exactly similar to the previous one, but for a space-like jet.
  // Also in this case we start from the final grand-children (which
  // are childless) of the particle which originates the jet, but in
  // this case we proceed "forward" (in the physical time picture)
  // towards the <!id>particleJetParent<!!id>.

  bool reconstructSpecialTimeLikeDecayingJet( const tShowerParticlePtr particleJetParent );
  // This is a special case of reconstruction, for a time-like jet
  // whose originating particle is a decaying particle. It is a 
  // special case, because the showering evolution is forward but
  // with reverse angular ordering, whereas the kinematic reconstruction
  // is done going forward, starting from <!id>particleJetParent<!!id>
  // and ending on the particle attached to the decaying vertex.
  // In this case, the result of the reconstruction, which is 
  // mainly the mass of the jet, is stored in the <!class>ShowerKinematics<!!class>
  // object associated with the particle attached to the decaying
  // vertex, although, by definition, it does not split.

  double momConsEq(const double & k, const Energy & root_s, const JetKinVect & jets);
  // the the term that is made to be zero for a <!id>k<!!id> that is
  // found by the next method.

  const double solveKfactor( const Energy & root_s, const JetKinVect & jets );
  // Given a vector of 5-momenta of jets, where the 3-momenta are the initial
  // ones before showering and the masses are reconstructed after the showering,
  // this method returns the overall scaling factor for the 3-momenta of the
  // vector of particles, vec{P}_i -> k * vec{P}_i, such to preserve energy-
  // momentum conservation, i.e. after the rescaling the center of mass 5-momentum 
  // is equal to the one specified in input, <!id>cmMomentum<!!id>. 
  // The method returns 0 if such factor cannot be found.
  
  Vector3 solveBoostBeta( const double k, const Lorentz5Momentum & newq, 
			  const Lorentz5Momentum & oldp);
  // Given a 5-momentum and a scale factor, the method returns the
  // Lorentz boost that transforms the 3-vector vec{momentum} --->
  // k*vec{momentum}. The method returns the null boost in the case no
  // solution exists. This will only work in the case where the
  // outgoing jet-momenta are parallel to the momenta of the particles
  // leaving the hard subprocess. 

  LorentzRotation solveBoost( const double k, const Lorentz5Momentum & newq, 
		      const Lorentz5Momentum & oldp);
  // More general solution, includes the case of non-parallel parent-
  // and jet momenta, involves a bit more numerical work, however.

  bool solveOverallCMframeBoost( const Lorentz5Momentum & pBeamHadron1,
                                 const Lorentz5Momentum & pBeamHadron2,
				 const Lorentz5Momentum & pBeamParton1,
				 const Lorentz5Momentum & pBeamParton2,
				 const Lorentz5Momentum & pHard1Initial,
				 const Lorentz5Momentum & pHard2Initial,
				 const Lorentz5Momentum & pHard1Intermediate,
				 const Lorentz5Momentum & pHard2Intermediate,
				 Lorentz5Momentum & pHard1Final,
				 Lorentz5Momentum & pHard2Final );
  // Given in input the following 5-momenta, all expressed in the Lab frame:
  // --- of the two beam hadrons: <!id>pBeamHadron1, pBeamHadron2<!!id>;
  // --- of the two beam partons: <!id>pBeamParton1, pBeamParton2<!!id>;
  // --- of the two incoming partons entering the hard subprocess at the 
  //     beginning (before showering): <!id>pHard1Initial, pHard2Initial<!!id>;
  // --- of the two incoming partons entering the hard subprocess immediately 
  //     after the reconstruction: <!id>pHard1Intermediate, pHard2Initermediate<!!id>;
  // the method calculates: 
  //  --- the corresponding final momenta of the two incoming partons 
  //      entering the hard subprocess: <!id>pHard1Final, pHard2Final<!!id>.
  // This is obtained by boosting <!id>pHard1Intermediate<!!id> and <!id>pHard2Intermediate<!!id>
  // longitudinally, i.e. along the respective beam directions, <!id>pBeamHadron1<!!id> 
  // and <!id>pBeamHadron2<!!id>, such to satisfy the two following conditions:
  //   1) the s_hat of the hard subprocess remains unchanged:
  //      <I>  ( pHard1Initial + pHard2Initial )^2 = 
  //           ( pHard1Final + pHard2Final)^2 </I>
  //   2) the rapidity of the hard subprocess remains unchanged:
  //      <I>  y_initial = y_final </I>
  //      where  y = 1/2 * log ( (e1+e2+pz1+pz2) / (e1+e2-pz1-pz2) )
  //      (built with the same <!id> pHard1Initial, pHard2Initial,
  //       pHard1Final, pHard2Final <!!id>).
  // Notice that, it is not enough to provide the beam partons, but
  // also the beam hadrons, because of the possible intrinsic transverse
  // momenta of partons inside the hadrons.
  // All the vectors are expressed in the Lab frame,
  // The returns false if something goes wrong, true otherwise.
 
  bool solveSpecialDIS_CMframeBoost( const Lorentz5Momentum & pLepton,
				     const Lorentz5Momentum & pBeamHadron,
				     const Lorentz5Momentum & pBeamParton,
				     const Lorentz5Momentum & pHardInitial,
				     const Lorentz5Momentum & pHardIntermediate,
				     Lorentz5Momentum & pHardFinal );
  // Given in input the following 5-momenta, all expressed in the Lab frame:
  // --- of the incoming lepton: <!id>pLepton<!!id> (assumed not radiating);
  // --- of the beam hadron: <!id>pBeamHadron<!!id>; 
  // --- of the beam parton: <!id>pBeamParton<!!id>;
  // --- of the incoming parton entering the hard subprocess at the 
  //     beginning (before showering): <!id>pHardInitial<!!id>;
  // --- of the incoming parton entering the hard subprocess immediately 
  //     after the reconstruction: <!id>pHardIntermediate<!!id>;
  // the method calculates: 
  //  --- the corresponding final momentum of the incoming parton 
  //      entering the hard subprocess: <!id>pHardFinal<!!id>.
  // This is obtained by boosting <!id>pHardIntermediate<!!id> longitudinally, 
  // i.e. along the beam direction <!id>pBeamHadron<!!id>, such to satisfy the 
  // two following condition:
  //   1) the s_hat of the hard subprocess remains unchanged:
  //      <I> ( pLepton + pHardInitial )^2 = ( pLepton + pHardFinal )^2 </I>
  // Notice that we are assuming that the lepton is not radiating
  // (QED or EWK radiation), therefore by not touching it at all
  // the <I> (x,Q^2) </I> is automatically preserved.
  // Notice that, it is not enough to provide the beam parton, but
  // also the beam hadron, because of the possible intrinsic transverse
  // momenta of partons inside the hadron.
  // All the vectors are expressed in the Lab frame,
  // The returns false if something goes wrong, true otherwise.

  vector<MEPtr> _specialProcesses;
 
};


}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of KinematicsReconstructor.
template <>
struct BaseClassTrait<Herwig::KinematicsReconstructor,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::KinematicsReconstructor>: public ClassTraitsBase<Herwig::KinematicsReconstructor> {
  static string className() { return "/Herwig++/KinematicsReconstructor"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "KinematicsReconstructor.icc"

#endif /* HERWIG_KinematicsReconstructor_H */
