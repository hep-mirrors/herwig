// -*- C++ -*-
#ifndef HERWIG_a1ThreePionDecayer_H
#define HERWIG_a1ThreePionDecayer_H
//
// This is the declaration of the <!id>a1ThreePionDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The  <!id>a1ThreePionDecayer<!!id> class is designed to implement the decay
//  of the a_1 to three pions. The model used is one based on the Novosibirsk
//  4 pion current used in TAUOLA. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TauDecayerBase.html">TauDecayerBase.h</a>,
// <a href="Tau4PionNovosibirskDecayer.html">Tau4PionNovosibirskDecayer.h</a>.
// 

#include "VectorMesonDecayerBase.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig++/Utilities/Kinematics.h"
// #include "a1ThreePionDecayer.fh"
// #include "a1ThreePionDecayer.xh"

namespace Herwig {

using namespace ThePEG;

class a1ThreePionDecayer: public VectorMesonDecayerBase {
  
public:
  
  inline a1ThreePionDecayer();
  inline a1ThreePionDecayer(const a1ThreePionDecayer &);
  virtual ~a1ThreePionDecayer();
  // Standard ctors and dtor.
  
public:
  
  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.
  
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.
  
  virtual vector<LorentzPolarizationVector> decayCurrent(const  bool, const int,
							 const int, const Particle &,
							 const ParticleVector &) const;
  // the hadronic currents    
  
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
  
  static ClassDescription<a1ThreePionDecayer> inita1ThreePionDecayer;
  // Describe a concrete class with persistent data.
  
  a1ThreePionDecayer & operator=(const a1ThreePionDecayer &);
  // Private and non-existent assignment operator.
  
private:
  
  inline complex<InvEnergy2> sigmaBreitWigner(Energy2) const;
  // breit-wigner for the sigma
  
  inline double a1FormFactor(Energy2) const;
  // the a_1 form factor

  inline Complex rhoBreitWigner(Energy2,int) const;
  // the rho breit wigner 

  inline Energy2 DParameter(int) const ;
  // the d parameter in rho the propagator

  inline double dhdq2Parameter(int) const ;
  // the dh/dq2 function in the rho propagator evaluated at q2=m2
  
  inline Energy2 hFunction(const Energy) const ;
  // the h function in the rho propagator
  
private:

  vector<Energy> _rhomass,_rhowidth,_prho;
  vector<Energy2>_hm2,_rhoD;
  vector<double> _dhdq2m2;
  // mass and width of the rho resonace and parameters for breit-wigner
  Energy _sigmamass,_sigmawidth,_psigma;
  // mass and width of the sigma meson and parameters for breit-wigner
  Energy _mpi;Energy2 _mpi2;
  // masses of the pions
  Energy2 _lambda2,_a1mass2;
  // parameter for the form-factor
  Complex _zsigma;
  // relative sigma rho coupling
  vector<Complex> _rhocoupling;
  // the couplings for the different rho resonances
  double _coupling;
  // overall coupling for the decay
  bool _localparameters;
  // use local values of the mass parameters

  vector<bool> _zerochan,_onechan,_twochan,_threechan;
  vector<double> _zerowgts,_onewgts,_twowgts,_threewgts;
  double _zeromax,_onemax,_twomax,_threemax;
  // parameters for the multi-channel integration

};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of a1ThreePionDecayer.
  template <>
  struct BaseClassTrait<Herwig::a1ThreePionDecayer,1> {
    typedef Herwig::VectorMesonDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::a1ThreePionDecayer>
    : public ClassTraitsBase<Herwig::a1ThreePionDecayer> {
    static string className() { return "/Herwig++/a1ThreePionDecayer"; }
    // Return the class name.
    static string library() { return "libHwVMDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "a1ThreePionDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "a1ThreePionDecayer.tcc"
#endif

#endif /* HERWIG_a1ThreePionDecayer_H */
