// -*- C++ -*-
#ifndef HERWIG_Evolver_H
#define HERWIG_Evolver_H
//
// This is the declaration of the InsideRangeShowerEvolver class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerConfig.h"
#include "PartnerFinder.h"
#include "ForwardEvolver.h"
#include "BackwardEvolver.h"
#include "KinematicsReconstructor.h"
#include "SplittingGenerator.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * 
 *  This class is responsible for the showering in a given scale range,
 *  including the kinematic reshuffling necessary for energy-momentum
 *  conservation to balance the recoil of the emissions.
 *  Furthermore, just for convenience (it has all the necessary information 
 *  to do that), it provides also the possibility of setting the final state 
 *  gluons on the effective gluon mass shell (rather than on the physical 
 *  massless shell), which is necessary at the end of the showering if the 
 *  cluster hadronization model has to be used.
 *
 *  @see GlobalParameters
 *  @see PartnerFinder
 *  @see ForwardShowerEvolver
 *  @see BackwardShowerEvolver
 *  @see KinematicsReconstructor
 *  @see SplittingGenerator
 */ 
class Evolver: public ThePEG::HandlerBase {

public:

  /**
   * See the comment on class KinematicsReconstructor about this typedef.
   */
  typedef map<tShowerParticlePtr, bool> MapShower;

  /**
   * Each element of this vector, described above, is associated with
   * a decaying jets, which includes the single on-shell decaying particle,
   * plus all its decay products (each of these is the parent of a
   * time-like jet). 
   */
  typedef vector<MapShower> MapShowerVector;

  /**
   * Standard ctors and dtor.
   */
  inline Evolver();
  inline Evolver(const Evolver &);
  virtual ~Evolver();

  /**
   * It should be called at the beginning of each collision, not
   * for the showering in each scale range. It does some cleaning
   * of internal, private, data.  
   */
  void clear();

  /**
   * It does the normal showering of the particles entering the hard 
   * subprocess. The ParticleCollisionHandler object is needed to 
   * access the PDF. If skipKinReco is true, then the kinematics 
   * reconstruction is skipped.
   */
  bool showerNormally(tEHPtr ch, 
		      const tShowerVarsPtr showerVariables, 
		      //const tMECorrectionPtr meCorrection,
		      ShowerParticleVector & particles,
		      bool skipKinReco=false) throw (Veto, Stop, Exception);

  /**
   * It does the (special) showering of a decay.
   */
  void showerDecay(tEHPtr ch, 
		   const tShowerVarsPtr showerVariables, 
		   //const tMECorrectionPtr meCorrection,
		   ShowerParticleVector &) throw (Veto, Stop, Exception);

  /**
   * It does the overall showering, between two width scales 
   * (or from a width scale and to the end). It is used only
   * in multi-scale showering when some radiating particle has
   * lifetime shorter than the typical hadronization time scale. 
   * The ParticleCollisionHandler object is needed to access the 
   * PDF. If skipKinReco is true, then the kinematics 
   * reconstruction is skipped.
   */
  void showerGlobally(tEHPtr & ch,  
		      const tShowerVarsPtr showerVariables, 
		      //const tMECorrectionPtr meCorrection,
		      ShowerParticleVector & particles,
		      bool skipKinReco = false) throw (Veto, Stop, Exception);

  /**
   * It forces the final state gluons on the effective gluon mass shell
   * (rather than on the physical massless shell). It also set properly
   * the various flags that are needed for the kinematics reconstruction
   * (but the latter is not done by this method).
   */
  void setEffectiveGluonMass(const Energy, const ShowerParticleVector &) 
    throw (Veto, Stop, Exception);

  /**
   * It does the kinematics reconstruction and the necessary 
   * reshuffling in order to conserve energy-momentum.
   * It is used mainly by the above methods, but is also used
   * once by ShowerHandler, at the end of the all showering procedure, 
   * in order to avoid the possibility of double reconstruction:
   * first with massless gluons (which is the standard procedure) 
   * and then with gluons put on the effective mass shell (in the
   * case we want to use Herwig++ cluster hadronization for the
   * hadronization rather than ThePEG string fragmentation).
   */
  bool reconstructISKinematics(tEHPtr & ch) 
    throw (Veto, Stop, Exception);
  bool reconstructKinematics(tEHPtr & ch) 
    throw (Veto, Stop, Exception);

  /**
   * Public access to the splitting generator
   */
  inline Ptr<SplittingGenerator>::transient_const_pointer splittingGenerator()
  { return _splittingGenerator; }

  /**
   * (PR remnant related?)
   */
  void makeRemnants(ShowerParticleVector &);

public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

protected:

  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();

  /** 
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<Evolver> initEvolver;

  /**
   * Private and non-existent assignment operator.
   */
  Evolver & operator=(const Evolver &);

  /**
   * Set false all the data of the input map. This is used for
   * resetting _mapShowerHardJets, and the elements of
   * the collection _collecMapShowerDecayJets, after that the
   * kinematical reconstruction has been performed.
   * This avoid, in multi-scale showering, to repeat unnecessary
   * kinematical reconstructions.
   */
  void setDoneMapShower(MapShower & mapShower);

  Ptr<PartnerFinder>::pointer _partnerFinder;
  Ptr<ForwardEvolver>::pointer _forwardEvolver;
  Ptr<BackwardEvolver>::pointer _backwardEvolver;
  Ptr<KinematicsReconstructor>::pointer _kinematicsReconstructor;
  Ptr<SplittingGenerator>::pointer _splittingGenerator;

  /**
   * This map is a collection of elements ( key = pointer, value = boolean flag), 
   * where the pointer points to a ShowerParticle object that enters
   * (as incoming or outcoming) the hard subprocess, and the boolean
   * flag tells whether or not the jet originated by such particle
   * needs the kinematical reconstruction. This is necessary when
   * some radiation has been emitted or if some of the "leaves", 
   * childless gluons have been forced on their effective mass shell. 
   */
  MapShower _mapShowerHardJets;

  /**
   * Each element of this collection is a map associated with a
   * decaying showering particle. More precisely, this map is 
   * made of elements ( key = pointer, value = boolean flag), 
   * where the pointer points to the decaying ShowerParticle object 
   * or one of the decay products, and the boolean flag tells whether 
   * or not the jet originated by such particle needs the kinematical 
   * reconstruction. This is necessary when some radiation has been 
   * emitted or if some of the "leaves", childless gluons have been 
   * forced on their effective mass shell. 
   */
  MapShowerVector _mapShowerDecayJets;

};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of Evolver.
 */
template <>
struct BaseClassTrait<Herwig::Evolver,1> {
  typedef ThePEG::HandlerBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Evolver>: public ClassTraitsBase<Herwig::Evolver> {

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/Evolver"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }

};

}

#include "Evolver.icc"

#endif /* HERWIG_Evolver_H */
