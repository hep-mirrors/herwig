// -*- C++ -*-
#ifndef HERWIG_VectorMeson3PionDecayer_H
#define HERWIG_VectorMeson3PionDecayer_H
//
// This is the declaration of the <!id>VectorMeson3PionDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>VectorMeson3PionDecayer<!!id> class is designed to perform the
//  decay of an I=0 meson to three pions via rho mesons including the
//  option of higher rho resonaces and a constant term. It is mainly intended for
//  the decays
//
//  omega -> pi+pi-pi0
//  phi   -> pi+pi-pi0
//
//  The default for the omega is to only include the contributions of the rho(770)
//  without a constant term, whereas for the phi the parameters from hep-ex/0303016 
//  (KLOE) which includes the rho(770) and a constant term to represent the effects
//  of the higher rho resonances. (The KLOE paper also included a omega contribution
//  but this is assumed to be non-resonant.)
//
//  To allow the easy addition of further modes the parameters for additional modes
//  can be set. The following must be specified for each mode
//
//  Incoming       - the PDG code for the incoming particle
//  Coupling       - the overall coupling for the decay
//  DirectCoupling - the relative coupling for the direct term
//  DirectPhase    - the phase of the coupling for the direct term
//  Rho2Coupling   - the relative coupling for the second rho multiplet
//  Rho2Phase      - the phase of the coupling for the second rho multiplet
//  Rho3Coupling   - the relative coupling for the third rho multiplet
//  Rho3Phase      - the phase of the coupling for the third rho multiplet
//  MaximumWeight  - the maximum weight for the integration of the channel
//  Rho1Weight     - the weight for the first rho multiplet in the 
//                   multichannel integration
//  Rho2Weight     - the weight for the second rho multiplet in the
//                   multichannel integration
//  Rho3Weight     - the weight for the third rho multiplet in the
//                   multichannel integration
//  Rho1Mass       - mass  of the first  rho multiplet
//  Rho2Mass       - mass  of the second rho multiplet
//  Rho3Mass       - mass  of the third  rho multiplet
//  Rho1Width      - width of the first  rho multiplet
//  Rho2Width      - width of the second rho multiplet
//  Rho3Width      - width of the third  rho multiplet
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:VectorMesonDecayerBase.html">VectorMesonDecayerBase.h</a>.
// 

#include "VectorMesonDecayerBase.h"
#include "Herwig++/Utilities/Kinematics.h"
// #include "VectorMeson3PionDecayer.fh"
// #include "VectorMeson3PionDecayer.xh"

namespace Herwig {
using namespace ThePEG;

class VectorMeson3PionDecayer: public VectorMesonDecayerBase {

public:

  inline VectorMeson3PionDecayer();
  inline VectorMeson3PionDecayer(const VectorMeson3PionDecayer &);
  virtual ~VectorMeson3PionDecayer();
  // Standard ctors and dtor.

public:


    virtual vector<LorentzPolarizationVector> 
    decayCurrent(const bool, const int, const int, 
		 const Particle &, const ParticleVector &) const;
    // the hadronic currents    

public:

  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.

  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.


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
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<VectorMeson3PionDecayer> initVectorMeson3PionDecayer;
  // Describe a concrete class with persistent data.

  VectorMeson3PionDecayer & operator=(const VectorMeson3PionDecayer &);
  // Private and non-existent assignment operator.

private:

  vector<double> _incoming;
  // vector storing the decaying particles for the different modes
  vector<double> _coupling;
  // the overall couplings for the decay
  vector<double> _directcoupling,_directphase;
  // relative coupling and phase for the direct term 
  vector<double> _rho2coupling,_rho2phase;
  // relative coupling and phase for the second rho multiplet
  vector<double> _rho3coupling,_rho3phase;
  // relative coupling and phase for the third rho multiplet
  vector<double> _maxwgt;
  // maximum weight for the integration of the channel
  vector<double> _rho1wgt;
  // weight for the first  rho multiplet in the integration
  vector<double> _rho2wgt;
  // weight for the second rho multiplet in the integration
  vector<double> _rho3wgt;
  // weight for the third  rho multiplet in the integration
  vector<double> _rho1mass;
  // mass of the first rho multiplet
  vector<double> _rho2mass;
  // mass of the  second rho multiplet
  vector<double> _rho3mass;
  // mass of the third rho multiplet
  vector<double> _rho1width;
  // width of the first rho multiplet
  vector<double> _rho2width;
  // width of the  second rho multiplet
  vector<double> _rho3width;
  // width of the third rho multiplet
  vector<double> _defaultmass;
  // use the default parameters for the rho masses and widths
  vector<vector <double> > _rho0const,_rhocconst;
  // constants for the running widths
  vector<vector<Energy> > _rhomass;
  vector<vector<Energy2> > _rhomass2;
  // rho mass parameters
  vector<vector <Complex> > _ccoupling;
  // couplings as complex numbers
  vector<int> _firstchannel;
  // store the first integration channel for a mode
};

}

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of VectorMeson3PionDecayer.
template <>
struct BaseClassTrait<Herwig::VectorMeson3PionDecayer,1> {
  typedef Herwig::VectorMesonDecayerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::VectorMeson3PionDecayer>
  : public ClassTraitsBase<Herwig::VectorMeson3PionDecayer> {
  static string className() { return "/Herwig++/VectorMeson3PionDecayer"; }
  // Return the class name.
  static string library() { return "libHwVMDecay.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "VectorMeson3PionDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson3PionDecayer.tcc"
#endif

#endif /* HERWIG_VectorMeson3PionDecayer_H */
