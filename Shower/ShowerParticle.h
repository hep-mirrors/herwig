// -*- C++ -*-
#ifndef HERWIG_ShowerParticle_H
#define HERWIG_ShowerParticle_H
//
// This is the declaration of the <!id>ShowerParticle<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class represents a particle in the showering process. <BR> 
// It has much less information than the ThePEG Particle, and it has some <BR> 
// specifics information useful only during the showering process. <BR>
//
// Notice that:
// <UL>
//   <LI> <!id>ShowerParticle<!!id> needs to inherit from <!id>HandlerBase<!!id> <BR>
//        (rather than, more simply, from <!id>ReferenceCounted<!!id>) becuse it needs <BR>
//        to create a <!id>ThePEG::Particle<!!id> object, and therefore the method <BR>
//        <!id>HandlerBase::getParticle( id )<!!id> is used. 
//   <LI> it has been necessary to define a new class <!class>ShowerColourLine<!!class> <BR>
//        in order to represent colour lines between <!id>ShowerParticle<!!id> objects <BR>
//        because the similar ThePEG class <!id>ColourLine<!!id> can be used only <BR>
//        with ThePEG Particle. 
//   <LI> for forward evolution, it is clear what does mean parent/child; <BR>
//        for backward evolution, however, it depends whether we want <BR>
//        to keep a physical picture or a Monte-Carlo effective one. <BR>
//        In the former case, an incoming particle (emitting particle) <BR> 
//        splits into an emitted particle and the emitting particle after <BR>
//        the emission: the latter two are then children of the <BR>
//        emitting particle, the parent. In the Monte-Carlo effective <BR> 
//        picture, we have that the particle close to the hard subprocess, <BR>
//        with higher (space-like) virtuality, splits into an emitted particle <BR>
//        and the emitting particle at lower virtuality: the latter two are, <BR>
//        in this case, the children of the first one, the parent. For obvious, <BR>
//        practical programming reasons, we choose the Monte-Carlo effective picture. 
//   <LI> the pointer to a <!class>SplitFun<!!class> object is set only in the case <BR>
//        that the particle has undergone to a shower emission; similarly, <BR> 
//        the pointer to a <!id>Decayer<!!id> object is set only in the case <BR>
//        that the particle has undergone to a decay. <BR>
//        In the case of particle connected directly to the hard subprocess, <BR>
//        there is no pointer to the hard subprocess, but there is a method <BR>
//        (<!id>isFromHardSubprocess()<!!id>) which returns true only in this case.
//   <LI> the spin density matrix (rho), and the decay matrix (D) are <BR>
//        implemented as complex bivectors, that we <I>typedef</I> (in 
//        <!class>ShowerConfig.h<!!class>) <BR>
//        as <!id>ComplexMatrix<!!id>, rather than as pair (one for the real part <BR>
//        and one for the imaginary part) of real true Matrix object as defined in <BR>
//        CLHEP (<!id>CLHEP/Matrix/Matrix.h<!!id> : notice that there is no <BR>
//        complex matrices defined in CLHEP), becaue the rho-D propagation is not <BR>
//        indeed a true matrix multiplication: nested for loops are used instead. <BR>
//        This is exactly the same approach as in Fortran Herwig. Furthermore, <BR> 
//        again similarly to what is done in Fortran Herwig, we keep one <BR>
//        single "matrix" for both rho and D, because we never need both rho and D <BR>
//        matrices for the same particle at the same time. In other words, when we <BR>
//        start evolving, even in a multi-scale showering, a time-like particle <BR>
//        (forward evolution), the "matrix" represents rho, whereas at the end of <BR>
//        the forward showering the "matrix" represents D; vice versa, for a space-like <BR>
//        particle (backward evolution), at the beginning the "matrix" represents D, <BR>
//        whereas at the end of backward showering the "matrix" represents rho. 
// </UL>
//
// ***LOOKHERE*** <BR>
//   --- The decayer has been defined as a transient pointer,
//       but maybe it shouldn't be transient at all, if for example
//       such decayer is not referenced by any other object... <BR> 
// ***endLOOKHERE*** 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:Particle.html">Particle.h</a>, <BR>
// <a href="http:ShowerConfig.html">ShowerConfig.h</a>, <BR>
// <a href="http:ShowerIndex.html">ShowerIndex.h</a>, <BR>
// <a href="http:ShowerColourLine.html">ShowerColourLine.h</a>, <BR>
// <a href="http:ShowerKinematics.html">ShowerKinematics.h</a>.
// 

#include "ShowerConfig.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ShowerKinematics.h"
#include "ShowerIndex.h"
//#include "ShowerColourLine.h"


namespace Herwig {

using namespace ThePEG;

  //class ThePEG::Particle;   // forward declaration

class ShowerParticle: public Particle {
private:
  //inline ShowerParticle();
  PPtr temp;
public:
  inline void * operator new(size_t);
  inline void operator delete(void *, size_t);

  inline ShowerParticle(tcEventPDPtr);
  inline ShowerParticle(const ShowerParticle &);
  inline ShowerParticle(const Particle &);
  virtual ~ShowerParticle();
  // Standard ctors and dtor.

  //explicit ShowerParticle(const ThePEG::Particle & inputP7Particle);
  // Create a <!id>ShowerParticle<!!id> object from a <!id>ThePEG::Particle<!!id> object.
  // This is useful at the beginning of the Showering, when we receive
  // the ThePEG particles (in the event record) from the hard subprocess.

  //PPtr createThePEGParticle() const;
  // Create a <!id>ThePEG::Particle<!!id> object from this <!id>ShowerParticle<!!id> object.
  // This is useful at the end of the Showering, when we have to update
  // the event record with ThePEG particles, using some of the 
  // <!id>ShowerParticle<!!id> objects produced by the showering 
  // (which ones depends on the degree of detail of the showering we want to keep
  //  in the event record). 

  //inline const ParticleData & data() const;
  //inline tcPDPtr dataPtr() const { return _pdptr; }
  //inline void dataPtr(const tcPDPtr & inputParticleDataPtr);
  // Access/Set the particle data.

  //inline const Lorentz5Momentum & momentum() const;
  //inline void momentum(const Lorentz5Momentum & inputMomentum);
  // Access/Set the 5-momentum vector of the component (parton or diquark).

  //inline const LorentzPoint & position() const;
  //inline void position(const LorentzPoint & inputPosition);
  // Access/Set the particle 4-position. It is eventually used by the
  // Hadronization (for the colour reconnection).

  //inline void transform(const LorentzRotation & r);
  //void deepTransform(const LorentzRotation & r);
  //inline void deepBoost(double bx, double by, double bz);
  // Do Lorentz transformations on this particle and its decay products.
  // ***LOOKHERE*** TO BE DEFINED : MAY BE IS ENOUGH ONE OF THE TWO

  //inline tShoColinePtr antiColourLine() const;
  //inline tShoColinePtr colourLine(bool anti = false) const;
  // Return the colour lines to which this particle is connected.

  //inline void setAntiColourLine(ShoColinePtr);
  //inline void setColourLine(ShoColinePtr, bool anti = false);
  // Set the colour lines to which this particle is connected.

  //inline tShoParPtr parent() const;
  //inline void parent(const tShoParPtr inputParent);
  // Access/Set the parent particle. 
  // Notice that in the case of backward evolution splitting,
  // the parent is the one with higher (space-like) virtuality,
  // that is "closer" to the hard subprocess

  inline bool isFinalState() const;
  inline void setFinalState( const bool );
  // Access/Set the flag that tells if the particle is final state
  // or initial state.

  inline bool isFromHardSubprocess() const;
  inline void setFromHardSubprocess(const bool);
  // It returns true is the particle has no parent, which means
  // that it is connected directly to the hard subprocess
  // (incoming if it is space-like, outgoing it is time-like);
  // in all other cases, it returns false.  

  //inline const ShowerParticleVector & children() const;
  // Return a const reference to the collection of decay products.

  //inline void addChild(const tShowerParticlePtr);
  // Add a child, setting child's parent pointer accordingly.

  //inline void removeChildren();
  // remove collection of children.

  bool addChildren(const tShowerParticleVector &);
  // Add a collection of children, setting children's parent pointer
  // accordingly. The input children can be either the decay products
  // of this particle, or the radiated particle and emitting particle
  // after the emission. In the latter case, when the emitting particle
  // is a decayed particle, therefore having already children, the
  // method automatically transfers the decay products to the
  // decayed particle after the emission (which should be one of
  // particle in <!id>inputChildren<!!id>). If some inconsistency is found,
  // then the method returns false; otherwise true.

  inline double sudAlpha() const;
  inline void sudAlpha(const double);
  inline double sudBeta() const;
  inline void sudBeta(const double);
  inline Energy sudPx() const;
  inline void sudPx(const Energy);
  inline Energy sudPy() const;
  inline void sudPy(const Energy);
  inline Energy sudPperp() const;
  inline Energy2 sudPperp2() const;
  // Access/Set Sudakov variables.
  // Notice that the <!id>ShowerKinematics<!!id> object is logically
  // associated more with the branching vertex than with the radiating 
  // particle itself, although it is stored in the <!id>ShowerParticle<!!id> 
  // object associated with the branching (radiating) particle. 
  // Furthermore, the branching products have not <!id>ShowerKinematics<!!id>
  // object (they eventually will have one only later if they branch).
  // Therefore this Sudakov variables can be considered as a useful
  // representation of the temporarily, preliminary, momentum of the 
  // <!id>ShowerParticle<!!id> object during the showering evolution, 
  // whereas the <!id>momentum<!id> member describes the "final", "real" one. 

  inline tSplittingFnPtr splitFun() const;
  inline void setSplittingFn(const tSplittingFnPtr);
  // Access/Set the <!class>SplitFun<!!class> object responsible of the 
  // eventual branching of this particle.

  //  inline tDecayerPtr decayer() const;
  //inline void decayer(const tDecayerPtr inputDecayer);
  // Access/Set the <!id>Decayer<!!id> object responsible of the 
  // eventual decay of this particle.

  inline ShoKinPtr & showerKinematics();
  void setShowerKinematics(const ShoKinPtr);
  // Access/Set the <!class>ShowerKinematics<!!class> object.

  inline vector<Energy> evolutionScales() const;
  inline const tShowerParticleVector & partners() const;
  // Return (a const reference to) the vector of evolution scales
  // (q-tilda scales) and of (pointers to) the partners corresponding to each 
  // considered interaction types (QCD, QED, EWK,...) defined in <!class>ShowerIndex<!!class>. 
  // The vector of (pointers to) the partners is needed only as the
  // most general way to decide in which frame the shower is described.

  inline void setEvolutionScale(const ShowerIndex::InteractionType, 
				const Energy);
  inline void setPartner(const ShowerIndex::InteractionType, 
			 const tShowerParticlePtr);
  // Set the scale/partner for the specified interaction.

  //inline bool isRhoDUpdate() const;
  //inline void setRhoDUpdate(const bool);
  // Access/Set the flag that tells if the rho-D matrix has been updated.

  //inline ComplexMatrix & rhoD();
  // Access the rho-D (spin density matrix or decay matrix) of the particle.

  inline bool isReconstructionFixedPoint() const;
  inline void setReconstructionFixedPoint(const bool);
  // Access/Set the flag that tells if the particle should be
  // treated in a special way during the kinematics reconstruction
  // (see <!class>KinematicsReconstructor<!!class> class). 
  // In practice, it returns true when either the particle is childless, 
  // or is a on-shell decaying particle (in which case we have to set the flag to
  // true before the showering of this particle: it is not enough to check 
  // if <!id>decayer()<!!id> is not null, because if it emits radiation
  // the decays products will be "transferred" to the particle
  // instance after the showering).

  void printInfo(); 
  // prints info of a single particle in some predefined way, mainly
  // for debug purposes. Used by <!id>deepPrintInfo()</!id>.

  void deepPrintInfo();
  // print info of all children

  void addChildrenEvtRec(const tStepPtr);

  Lorentz5Momentum sumParentsMomenta();
  // get the sum of all parents momenta (only for plotting)

  tShowerParticleVector getFSChildren();
  // get a list of all children of _particle that are in the final
  // state

  inline tPPtr getThePEGBase();
  inline void setThePEGBase(const tPPtr& );

protected:

  inline virtual PPtr clone() const;
  inline virtual PPtr fullclone() const;
  // Standard clone methods.

private:

  ShowerParticle & operator=(const ShowerParticle &);
  //  Private and non-existent assignment operator.

  bool _isFinalState;
  //bool _rhoDUpdate;
  bool _reconstructionFixedPoint;
  bool _isFromHardSubprocess;

  double _sudAlpha;
  double _sudBeta;

  Energy _sudPx;
  Energy _sudPy;

  tSplittingFnPtr _splitFun;        
  ShoKinPtr _showerKinematics;

  vector<Energy> _scales;
  tShowerParticleVector _partners;
  //ComplexMatrix _rhoD;

  tPPtr _thePEGBase;

  static ClassDescription<ShowerParticle> initShowerParticle;
 
};

}

namespace ThePEG {

template <>
struct BaseClassTrait<Herwig::ShowerParticle,1> {
  typedef EventRecordBase NthBase;
};
 
template <>
struct ClassTraits<Herwig::ShowerParticle>: 
    public ClassTraitsBase<Herwig::ShowerParticle> {
  static string className() { return "/Herwig++/ShowerParticle"; }
  static TPtr create() { return TPtr::Create(Herwig::ShowerParticle(tcEventPDPtr())); }
};

}



#include "ShowerParticle.icc"

#endif /* HERWIG_ShowerParticle_H */
