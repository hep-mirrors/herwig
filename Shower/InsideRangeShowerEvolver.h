// -*- C++ -*-
#ifndef HERWIG_InsideRangeShowerEvolver_H
#define HERWIG_InsideRangeShowerEvolver_H
//
// This is the declaration of the <!id>InsideRangeShowerEvolver<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is responsible for the showering in a given scale range, <BR>
// including the kinematic reshuffling necessary for energy-momentum <BR>
// conservation to balance the recoil of the emissions. <BR>
// Furthermore, just for convenience (it has all the necessary information to do that), <BR>
// it provides also the possibility of setting the final state gluons on the <BR>
// effective gluon mass shell (rather than on the physical massless shell), <BR>
// which is necessary at the end of the showering if the cluster hadronization <BR>
// model has to be used.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:GlobalParameters.html">GlobalParameters.h</a>, <BR>
// <a href="http:PartnerFinder.html">PartnerFinder.h</a>, <BR>
// <a href="http:ForwardShowerEvolver.html">ForwardShowerEvolver.h</a>, <BR>
// <a href="http:BackwardShowerEvolver.html">BackwardShowerEvolver.h</a>, <BR>
// <a href="http:KinematicsReconstructor.html">KinematicsReconstructor.h</a>, <BR>
// <a href="http:SplittingGenerator.html">SplittingGenerator.h</a>.
// 

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/PartialCollisionHandler.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerConfig.h"
#include "PartnerFinder.h"
#include "ForwardShowerEvolver.h"
#include "BackwardShowerEvolver.h"
#include "KinematicsReconstructor.h"
#include "SplittingGenerator.h"


namespace Herwig {

using namespace ThePEG;

class InsideRangeShowerEvolver: public ThePEG::HandlerBase {

public:

  typedef map<tShowerParticlePtr, bool> MapShower;
  // See the comment on class <!class>KinematicsReconstructor<!!class> about this typedef.

  typedef vector<MapShower> CollecMapShower;
  // Each element of this vector, described above, is associated with
  // a decaying jets, which includes the single on-shell decaying particle,
  // plus all its decay products (each of these is the parent of a
  // time-like jet). 

  inline InsideRangeShowerEvolver();
  inline InsideRangeShowerEvolver(const InsideRangeShowerEvolver &);
  virtual ~InsideRangeShowerEvolver();
  // Standard ctors and dtor.

  void clear();
  // It should be called at the beginning of each collision, not
  // for the showering in each scale range. It does some cleaning
  // of internal, private, data.  

  void showerNormally( tPartCollHdlPtr ch, 
		       const tShoConstrPtr showerConstrainer, 
		       const tMECorrectionPtr meCorrection,
		       ShowerParticleVector & particles,
		       bool skipKinReco = false ) throw (Veto, Stop, Exception);
  // It does the normal showering of the particles entering the hard subprocess.
  // The <!id>ParticleCollisionHandler<!!id> object is needed to access the PDF.
  // If <!id>skipKinReco<!!id> is true, then the kinematics reconstruction is skipped.

  void showerDecay( tPartCollHdlPtr ch, 
		    const tShoConstrPtr showerConstrainer, 
		    const tMECorrectionPtr meCorrection,
		    ShowerParticleVector & particles ) throw (Veto, Stop, Exception);
  // It does the (special) showering of a decay.

  void showerGlobally( tPartCollHdlPtr & ch,  
		       const tShoConstrPtr showerConstrainer, 
		       const tMECorrectionPtr meCorrection,
		       ShowerParticleVector & particles,
		       bool skipKinReco = false ) throw (Veto, Stop, Exception);
  // It does the overall showering, between two width scales 
  // (or from a width scale and to the end). It is used only
  // in multi-scale showering when some radiating particle has
  // lifetime shorter than the typical hadronization time scale. 
  // The <!id>ParticleCollisionHandler<!!id> object is needed to access the PDF.
  // If <!id>skipKinReco<!!id> is true, then the kinematics reconstruction is skipped.

  void setEffectiveGluonMass( const Energy effectiveGluonMass,
			      const ShowerParticleVector & particles ) throw (Veto, Stop, Exception);
  // It forces the final state gluons on the effective gluon mass shell
  // (rather than on the physical massless shell). It also set properly
  // the various flags that are needed for the kinematics reconstruction
  // (but the latter is not done by this method).

  bool reconstructKinematics( tPartCollHdlPtr & ch ) throw (Veto, Stop, Exception);
  // It does the kinematics reconstruction and the necessary 
  // reshuffling in order to conserve energy-momentum.
  // It is used mainly by the above methods, but is also used
  // once by <!class>ShowerHandler<!!class>, at the end of the all showering procedure, 
  // in order to avoid the possibility of double reconstruction:
  // first with massless gluons (which is the standard procedure) 
  // and then with gluons put on the effective mass shell (in the
  // case we want to use Herwig++ cluster hadronization for the
  // hadronization rather than ThePEG string fragmentation).

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

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

  static ClassDescription<InsideRangeShowerEvolver> initInsideRangeShowerEvolver;
  // Describe a concrete class with persistent data.

  InsideRangeShowerEvolver & operator=(const InsideRangeShowerEvolver &);
  // Private and non-existent assignment operator.

  void setDoneMapShower(MapShower & mapShower);
  // Set false all the data of the input map. This is used for
  // resetting _mapShowerHardJets, and the elements of
  // the collection _collecMapShowerDecayJets, after that the
  // kinematical reconstruction has been performed.
  // This avoid, in multi-scale showering, to repeat unnecessary
  // kinematical reconstructions.

  Ptr<PartnerFinder>::pointer _pointerPartnerFinder;
  Ptr<ForwardShowerEvolver>::pointer _pointerForwardShowerEvolver;
  Ptr<BackwardShowerEvolver>::pointer _pointerBackwardShowerEvolver;
  Ptr<KinematicsReconstructor>::pointer _pointerKinematicsReconstructor;
  Ptr<SplittingGenerator>::pointer _pointerSplittingGenerator;

  MapShower _mapShowerHardJets;
  // This map is a collection of elements ( key = pointer, value = boolean flag), 
  // where the pointer points to a ShowerParticle object that enters
  // (as incoming or outcoming) the hard subprocess, and the boolean
  // flag tells whether or not the jet originated by such particle
  // needs the kinematical reconstruction. This is necessary when
  // some radiation has been emitted or if some of the "leaves", 
  // childless gluons have been forced on their effective mass shell. 

  CollecMapShower _collecMapShowerDecayJets;
  // Each element of this collection is a map associated with a
  // decaying showering particle. More precisely, this map is 
  // made of elements ( key = pointer, value = boolean flag), 
  // where the pointer points to the decaying ShowerParticle object 
  // or one of the decay products, and the boolean flag tells whether 
  // or not the jet originated by such particle needs the kinematical 
  // reconstruction. This is necessary when some radiation has been 
  // emitted or if some of the "leaves", childless gluons have been 
  // forced on their effective mass shell. 

};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of InsideRangeShowerEvolver.
template <>
struct BaseClassTrait<Herwig::InsideRangeShowerEvolver,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::InsideRangeShowerEvolver>: public ClassTraitsBase<Herwig::InsideRangeShowerEvolver> {
  static string className() { return "/Herwig++/InsideRangeShowerEvolver"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "InsideRangeShowerEvolver.icc"

#endif /* HERWIG_InsideRangeShowerEvolver_H */
