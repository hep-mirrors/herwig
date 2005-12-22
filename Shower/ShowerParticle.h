// -*- C++ -*-
#ifndef HERWIG_ShowerParticle_H
#define HERWIG_ShowerParticle_H
//
// This is the declaration of the ShowerParticle class.

#include "ShowerConfig.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.fh"
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


  /** \ingroup Shower
   *  This class represents a particle in the showering process. <BR>
   *  It has much less information than the ThePEG Particle, and it has some <BR> 
   *  specifics information useful only during the showering process. <BR>
   * 
   *  Notice that:
   *  <UL>
   *    <LI> ShowerParticle needs to inherit from HandlerBase <BR>
   *         (rather than, more simply, from ReferenceCounted becuse it needs <BR>
   *         to create a ThePEG::Particle object, and therefore the method <BR>
   *         HandlerBase::getParticle( id ) is used. 
   *    <LI> it has been necessary to define a new class ShowerColourLine <BR>
   *         in order to represent colour lines between ShowerParticle objects <BR>
   *         because the similar ThePEG class ColourLine can be used only <BR>
   *         with ThePEG Particle. 
   *    <LI> for forward evolution, it is clear what does mean parent/child; <BR>
   *         for backward evolution, however, it depends whether we want <BR>
   *         to keep a physical picture or a Monte-Carlo effective one. <BR>
   *         In the former case, an incoming particle (emitting particle) <BR> 
   *         splits into an emitted particle and the emitting particle after <BR>
   *         the emission: the latter two are then children of the <BR>
   *         emitting particle, the parent. In the Monte-Carlo effective <BR> 
   *         picture, we have that the particle close to the hard subprocess, <BR>
   *         with higher (space-like) virtuality, splits into an emitted particle <BR>
   *         and the emitting particle at lower virtuality: the latter two are, <BR>
   *         in this case, the children of the first one, the parent. For obvious, <BR>
   *         practical programming reasons, we choose the Monte-Carlo effective picture. 
   *    <LI> the pointer to a SplitFun object is set only in the case <BR>
   *         that the particle has undergone to a shower emission; similarly, <BR> 
   *         the pointer to a Decayer object is set only in the case <BR>
   *         that the particle has undergone to a decay. <BR>
   *         In the case of particle connected directly to the hard subprocess, <BR>
   *         there is no pointer to the hard subprocess, but there is a method <BR>
   *         isFromHardSubprocess() which returns true only in this case.
   *    <LI> the spin density matrix (rho), and the decay matrix (D) are <BR>
   *         implemented as complex bivectors, that we <I>typedef</I> (in 
   *         ShowerConfig) <BR>
   *         as ComplexMatrix, rather than as pair (one for the real part <BR>
   *         and one for the imaginary part) of real true Matrix object 
   *         as defined in <BR>
   *         CLHEP (CLHEP/Matrix/Matrix.h : notice that there is no <BR>  
   *         complex matrices defined in CLHEP), becaue the rho-D propagation is not <BR>
   *         indeed a true matrix multiplication: nested for loops are used instead. <BR>
   *         This is exactly the same approach as in Fortran Herwig. Furthermore, <BR> 
   *         again similarly to what is done in Fortran Herwig, we keep one <BR>
   *         single "matrix" for both rho and D, because we never need both 
   *         rho and D <BR>
   *         matrices for the same particle at the same time. 
   *         In other words, when we <BR>
   *         start evolving, even in a multi-scale showering, a time-like particle <BR>
   *         (forward evolution), the "matrix" represents rho, whereas at the end of <BR>
   *         the forward showering the "matrix" represents D; 
   *         vice versa, for a space-like <BR>
   *         particle (backward evolution), at the beginning the "matrix" 
   *         represents D, <BR>
   *         whereas at the end of backward showering the "matrix" represents rho. 
   *  </UL>
   *
   *  ***LOOKHERE*** <BR>
   *    --- The decayer has been defined as a transient pointer,
   *        but maybe it shouldn't be transient at all, if for example
   *        such decayer is not referenced by any other object... <BR> 
   *
   *  ***endLOOKHERE*** 
   *
   *  @see Particle
   *  @see ShowerConfig
   *  @see ShowerIndex
   *  @see ShowerColourLine
   *  @see ShowerKinematics
   */
class ShowerParticle: public Particle {

private:

  //inline ShowerParticle();
  PPtr temp;

public:

  inline void * operator new(size_t);
  inline void operator delete(void *, size_t);

  /**
   * Standard ctors and dtor.
   */
  inline ShowerParticle(tcEventPDPtr);
  inline ShowerParticle(const ShowerParticle &);
  inline ShowerParticle(const Particle &);
  virtual ~ShowerParticle();

  // Create a ShowerParticle object from a ThePEG::Particle object.
  // This is useful at the beginning of the Showering, when we receive
  // the ThePEG particles (in the event record) from the hard subprocess.
  //explicit ShowerParticle(const ThePEG::Particle & inputP7Particle);

  // Create a ThePEG::Particle object from this ShowerParticle object.
  // This is useful at the end of the Showering, when we have to update
  // the event record with ThePEG particles, using some of the 
  // ShowerParticle objects produced by the showering 
  // (which ones depends on the degree of detail of the showering 
  //  we want to keep in the event record). 
  //PPtr createThePEGParticle() const;

  // Access/Set the particle data.
  //inline const ParticleData & data() const;
  //inline tcPDPtr dataPtr() const { return _pdptr; }
  //inline void dataPtr(const tcPDPtr & inputParticleDataPtr);

  // Access/Set the 5-momentum vector of the component (parton or diquark).
  //inline const Lorentz5Momentum & momentum() const;
  //inline void momentum(const Lorentz5Momentum & inputMomentum);

  // Access/Set the particle 4-position. It is eventually used by the
  // Hadronization (for the colour reconnection).
  //inline const LorentzPoint & position() const;
  //inline void position(const LorentzPoint & inputPosition);

  // Do Lorentz transformations on this particle and its decay products.
  // ***LOOKHERE*** TO BE DEFINED : MAY BE IS ENOUGH ONE OF THE TWO
  //inline void transform(const LorentzRotation & r);
  //void deepTransform(const LorentzRotation & r);
  //inline void deepBoost(double bx, double by, double bz);

  // Return the colour lines to which this particle is connected.
  //inline tShoColinePtr antiColourLine() const;
  //inline tShoColinePtr colourLine(bool anti = false) const;

  // Set the colour lines to which this particle is connected.
  //inline void setAntiColourLine(ShoColinePtr);
  //inline void setColourLine(ShoColinePtr, bool anti = false);

  // Access/Set the parent particle. 
  // Notice that in the case of backward evolution splitting,
  // the parent is the one with higher (space-like) virtuality,
  // that is "closer" to the hard subprocess
  //inline tShoParPtr parent() const;
  //inline void parent(const tShoParPtr inputParent);

  /**
   * Access/Set the flag that tells if the particle is final state
   * or initial state.
   */
  inline bool isFinalState() const;
  inline void setFinalState( const bool );

  /**
   * Access/Set the flag that tells if the particle is initiating a
   * time like shower when it has been emitted in an initial state shower.
   */
  inline bool initiatesTLS() const;
  inline void setInitiatesTLS( const bool );

  /**
   * It returns true is the particle has no parent, which means
   * that it is connected directly to the hard subprocess
   * (incoming if it is space-like, outgoing it is time-like);
   * in all other cases, it returns false.  
   */
  inline bool isFromHardSubprocess() const;
  inline void setFromHardSubprocess(const bool);

  //Return a const reference to the collection of decay products.
  //inline const ShowerParticleVector & children() const;

  // Add a child, setting child's parent pointer accordingly.
  //inline void addChild(const tShowerParticlePtr);

  // remove collection of children.
  //inline void removeChildren();

  /**
   * Add a collection of children, setting children's parent pointer
   * accordingly. The input children can be either the decay products
   * of this particle, or the radiated particle and emitting particle
   * after the emission. In the latter case, when the emitting particle
   * is a decayed particle, therefore having already children, the
   * method automatically transfers the decay products to the
   * decayed particle after the emission (which should be one of
   * particle in inputChildren). If some inconsistency is found,
   * then the method returns false; otherwise true.
   */
  bool addChildren(const tShowerParticleVector &);

  /**
   * Access/Set Sudakov variables.
   * Notice that the ShowerKinematics object is logically
   * associated more with the branching vertex than with the radiating 
   * particle itself, although it is stored in the ShowerParticle 
   * object associated with the branching (radiating) particle. 
   * Furthermore, the branching products have not ShowerKinematics
   * object (they eventually will have one only later if they branch).
   * Therefore this Sudakov variables can be considered as a useful
   * representation of the temporarily, preliminary, momentum of the 
   * ShowerParticle object during the showering evolution, 
   * whereas the momentum member describes the "final", "real" one.
   */
 
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

  /**
   * Access/Set the SplitFun object responsible of the 
   * eventual branching of this particle.
   */
  inline tSplittingFnPtr splitFun() const;
  inline void setSplittingFn(const tSplittingFnPtr);

  // Access/Set the Decayer object responsible of the 
  // eventual decay of this particle.
  //inline tDecayerPtr decayer() const;
  //inline void decayer(const tDecayerPtr inputDecayer);

  /**
   * Access/Set the ShowerKinematics object.
   */
  inline ShoKinPtr & showerKinematics();
  void setShowerKinematics(const ShoKinPtr);

  /**
   * Return (a const reference to) the vector of evolution scales
   * (q-tilda scales) and of (pointers to) the partners corresponding to each 
   * considered interaction types (QCD, QED, EWK,...) defined in ShowerIndex. 
   * The vector of (pointers to) the partners is needed only as the
   * most general way to decide in which frame the shower is described.
   */
  inline vector<Energy> evolutionScales() const;
  inline const tShowerParticleVector & partners() const;

  inline void setEvolutionScale(const ShowerIndex::InteractionType, 
				const Energy);

  /**
   * Set the scale/partner for the specified interaction.
   */
  inline void setPartner(const ShowerIndex::InteractionType, 
			 const tShowerParticlePtr);

  //Access/Set the flag that tells if the rho-D matrix has been updated.
  //inline bool isRhoDUpdate() const;
  //inline void setRhoDUpdate(const bool);

  // Access the rho-D (spin density matrix or decay matrix) of the particle.
  //inline ComplexMatrix & rhoD();

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
  inline bool isReconstructionFixedPoint() const;
  inline void setReconstructionFixedPoint(const bool);

  /**
   * Prints info of a single particle in some predefined way, mainly
   * for debug purposes. Used by deepPrintInfo().
   */
  void printInfo(); 

  /**
   * Print info of all children.
   */
  void deepPrintInfo();

  void addChildrenEvtRec(const tStepPtr);

  /**
   * Get the sum of all parents momenta (only for plotting).
   */
  Lorentz5Momentum sumParentsMomenta();

  /**
   * Get a list of all children of _particle that are in the final state.
   */
  tShowerParticleVector getFSChildren();

  inline tcPPtr getThePEGBase();
  //inline void setThePEGBase(const tcPPtr& );
  inline void x(double x) { _x = x; }
  inline double x() const { return _x; }

protected:

  /**
   * Standard clone methods.
   */
  inline virtual PPtr clone() const;
  inline virtual PPtr fullclone() const;

private:

  /**
   * Private and non-existent assignment operator.
   */
  ShowerParticle & operator=(const ShowerParticle &);

  bool _isFinalState;
  //bool _rhoDUpdate;
  bool _reconstructionFixedPoint;
  bool _isFromHardSubprocess;
  bool _initiatesTLS;

  double _sudAlpha;
  double _sudBeta;

  Energy _sudPx;
  Energy _sudPy;
  double _x;

  tSplittingFnPtr _splitFun;        
  ShoKinPtr _showerKinematics;

  vector<Energy> _scales;
  tShowerParticleVector _partners;
  //ComplexMatrix _rhoD;

  tcPPtr _thePEGBase;

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
