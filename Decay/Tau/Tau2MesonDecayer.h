// -*- C++ -*-
#ifndef HERWIG_Tau2MesonDecayer_H
#define HERWIG_Tau2MesonDecayer_H
//
// This is the declaration of the <!id>Tau2MesonDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is designed to implement the decay of the tau to two scalar
// mesons using the models of either Kuhn and Santamaria (Z. Phys. C48, 445 (1990))
// or Gounaris and Sakurai Phys. Rev. Lett. 21, 244 (1968). The mixing parameters
// are taken from Phys.Rev.D61:112002,2000 (CLEO) for the decay pi+/- pi0 and
// the CLEO version of TAUOLA for the K pi decays.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TauDecayerBase.html">TauDecayerBase.h</a>.
// 
//  Author: Peter Richardson
//

#include "TauDecayerBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Utilities/Kinematics.h"

namespace Herwig {

using namespace ThePEG;
  
class Tau2MesonDecayer: public TauDecayerBase {
  
public:
  
  vector<LorentzPolarizationVector>  hadronCurrent(bool vertex, const int, const int,
						   const Particle &,
						   const ParticleVector &) const;
    // the hadronic current for this decay mode

public:
  
  inline Tau2MesonDecayer();
  inline Tau2MesonDecayer(const Tau2MesonDecayer &);
  virtual ~Tau2MesonDecayer();
  // Standard ctors and dtor.
  
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
  
    static ClassDescription<Tau2MesonDecayer> initTau2MesonDecayer;
  // Describe a concrete class with persistent data.
  
  Tau2MesonDecayer & operator=(const Tau2MesonDecayer &);
  // Private and non-existent assignment operator.

private:

  inline Complex BreitWigner(Energy2, unsigned int, unsigned int, unsigned int) const;
  // p-wave breit wigner for form-factors
  
  inline double GSModelDParameter(const unsigned int)const ;
  // the d parameter in the GS model for the propagator
  
  inline double GSModeldhdq2Parameter(const unsigned int)const ;
  // the dh/dq2 functino in the GS model for the propagator evaluated at q2=m2
  
  inline double GSModelhFunction(const unsigned int,const Energy)const ;
  // the h function in the GS model
  
  inline Energy pcm(const unsigned int, const Energy) const;
  // the momentum of the decay products for the running width 

private:
  
  vector<double> _piwgt,_kwgt;
  // weights for the different resonances
  int _pimodel,_Kmodel;
  // model to use for the propagator
  
  mutable vector<bool> _pichannels,_K0channels,_Kpluschannels,_KKchannels;
  mutable vector<double> _piwgts,_K0wgts,_Kpluswgts,_KKwgts;
  mutable double _pimax,_K0max,_Kplusmax,_KKmax;
  // weights for the integration
  
  // options not to use the physical masses and widths for the rho
  bool _rhoparameters;
  vector<Energy> _rhomasses,_rhowidths;
  
  // options not to use the physical masses and widths for the K*
  bool _Kstarparameters;
  vector<Energy> _Kstarmasses,_Kstarwidths;

private:

  // parameters for the Breit-Wigners

  // existance of the resonance
  vector<Energy> _mass;
  // the masses of the resonances
  vector<Energy> _width;
  // the widths of the resonances
  vector<Energy2> _mass2;
  // squares of the masses
  vector<Energy2> _massw;
  // product of the mass and the width
  vector<Energy> _massa,_massb;
  // masses for momentum
  vector<Energy>  _mom;
  // momentum
  vector<InvEnergy2> _dhdq2;
  vector<InvEnergy2> _hm2;
  vector<double> _dparam;
  //parameters for the GS form of the Breit-Wigner
  
};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of Tau2MesonDecayer.
  template <>
  struct BaseClassTrait<Herwig::Tau2MesonDecayer,1> {
    typedef Herwig::TauDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Tau2MesonDecayer>
    : public ClassTraitsBase<Herwig::Tau2MesonDecayer> {
    static string className() { return "/Herwig++/Tau2MesonDecayer"; }
    // Return the class name.
    static string library() { return "libHwTauDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "Tau2MesonDecayer.icc"

#endif /* HERWIG_Tau2MesonDecayer_H */
