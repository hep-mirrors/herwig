// -*- C++ -*-
#ifndef HERWIG_a1ThreePionCLEODecayer_H
#define HERWIG_a1ThreePionCLEODecayer_H
//
// This is the declaration of the <!id>a1ThreePionCLEODecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The  <!id>a1ThreePionCLEODecayer<!!id> class is designed to implement the decay
//  of the a_1 to three pions using the model of Phys.Rev.D61:012002,2000,
//  (hep-ex/9902022) (CLEO) which was fitted to the one charged and two neutral pion
//  channel for the charged a_1 decay in tau -> a_1 -> three pi. The other modes
//  are then infered from this using isospin. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="VectorMesonDecayerBase.html">VectorMesonDecayerBase.h</a>.
// 

#include "VectorMesonDecayerBase.h"
#include "Herwig++/Utilities/Kinematics.h"
// #include "a1ThreePionCLEODecayer.fh"
// #include "a1ThreePionCLEODecayer.xh"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

class a1ThreePionCLEODecayer: public VectorMesonDecayerBase {
  
public:
  
  inline a1ThreePionCLEODecayer();
  inline a1ThreePionCLEODecayer(const a1ThreePionCLEODecayer &);
  virtual ~a1ThreePionCLEODecayer();
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
  
  static ClassDescription<a1ThreePionCLEODecayer> inita1ThreePionCLEODecayer;
  // Describe a concrete class with persistent data.
  
  a1ThreePionCLEODecayer & operator=(const a1ThreePionCLEODecayer &);
  // Private and non-existent assignment operator.
  
private:

  // breit wigner for the rho
  inline Complex rhoBreitWigner(int, Energy2,int) const;

  // breit wigner for the sigma
  inline Complex sigmaBreitWigner(Energy2,int) const;
  
  // breit wigner for the f_0
  inline Complex f0BreitWigner(Energy2,int) const;

  // breit wigner for the f_2
  inline Complex f2BreitWigner(Energy2,int) const;
  
private:
  
  vector<Energy> _rhomass,_rhowidth,_prhocc,_prhoc0;
  // masses and widths of the rho resonaces and parameters for breit-wigner
  Energy _f2mass,_f2width,_pf2cc,_pf200,_f0mass,_f0width,_pf0cc,_pf000;
  // masses and widths of the f_2(1270) and f_0(1370)
  Energy _sigmamass,_sigmawidth,_psigmacc,_psigma00;
  // mass and width of the sigma meson and parameters for breit-wigner
  Energy _mpi0,_mpic;
  // masses of the pions

  InvEnergy _coupling;
  // overall coupling for the decay

  vector<double> _rhomagP,_rhophaseP;
  vector<Complex> _rhocoupP;
  vector<InvEnergy2> _rhomagD;vector<double>_rhophaseD;
  vector<complex<InvEnergy2> > _rhocoupD;
  // couplings of the rho resonances (beta1-4 in CLEO paper)

  InvEnergy2 _f2mag;double _f2phase;
  complex<InvEnergy2> _f2coup;
  // couplings of the f_2 resonance (beta5)

  double _f0mag,_f0phase;
  Complex _f0coup;
  // couplings of the f_0 resonance (beta6)

  double _sigmamag,_sigmaphase;
  Complex _sigmacoup;
  // couplings of the sigma resonance (beta7)

  bool _localparameters;
  // use local values of the mass parameters
  
  mutable vector<bool> _zerochan,_onechan,_twochan,_threechan;
  mutable vector<double> _zerowgts,_onewgts,_twowgts,_threewgts;
  mutable double _zeromax,_onemax,_twomax,_threemax;
  // parameters for the multi-channel integration
};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of a1ThreePionCLEODecayer.
  template <>
  struct BaseClassTrait<Herwig::a1ThreePionCLEODecayer,1> {
    typedef Herwig::VectorMesonDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::a1ThreePionCLEODecayer>
    : public ClassTraitsBase<Herwig::a1ThreePionCLEODecayer> {
    static string className() { return "/Herwig++/a1ThreePionCLEODecayer"; }
    // Return the class name.
    static string library() { return "libHwVMDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "a1ThreePionCLEODecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "a1ThreePionCLEODecayer.tcc"
#endif

#endif /* HERWIG_a1ThreePionCLEODecayer_H */
