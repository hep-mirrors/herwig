// -*- C++ -*-
#ifndef HERWIG_Tau2MesonPhotonDecayer_H
#define HERWIG_Tau2MesonPhotonDecayer_H
//
// This is the declaration of the <!id>Tau2MesonPhotonDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class implements the decay tau+/- -> pi+/- pi0 gamma via
//  an intermediate omega. It inherits from the <!id>TauDecayerBase<!!id> class
//  and specifies the phase-space integration and implements the hadronic current.
//
//  The model is based on the one used in TAUOLA, Comput.Phys.Commun.76:361-380,1993.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TauDecayerBase.html">TauDecayerBase.h</a>.
// 
//  Author: Peter Richardson
//

#include "TauDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

class Tau2MesonPhotonDecayer: public TauDecayerBase {
    
public:
  
  vector<LorentzPolarizationVector> hadronCurrent(bool vertex, const int, const int,
						  const Particle &,
						  const ParticleVector &) const;
  // the hadronic current for this decay mode
  
public:
  
  inline Tau2MesonPhotonDecayer();
  inline Tau2MesonPhotonDecayer(const Tau2MesonPhotonDecayer &);
  virtual ~Tau2MesonPhotonDecayer();
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
  
  static ClassDescription<Tau2MesonPhotonDecayer> initTau2MesonPhotonDecayer;
  // Describe a concrete class with persistent data.
  
  Tau2MesonPhotonDecayer & operator=(const Tau2MesonPhotonDecayer &);
  // Private and non-existent assignment operator.
  
private:
  
  inline Complex FFunction(Energy2,int) const;
  // calculate the F function at a given scale 

  inline Complex BreitWigner(Energy2,unsigned int) const;
  // fixed width Breit wigner
  
private:
  
  Energy2 _grho;
  // coupling of the rho
  
  InvEnergy _grhoomegapi;
  // coupling of the rho to the omega and a pion
  
  vector<double> _resweights;
  // weights of the different rho resonances in the current
  
  mutable vector<bool> _on;
  mutable vector<double> _weights;
  mutable double _maxwgt;
  //weights for the integration    
  
  bool _rhoparameters;
  vector<Energy> _rhomasses,_rhowidths; 
  // use local parameters for the rho resonances rather than from the particle data
  // objects
  
  bool _omegaparameters;
  Energy _omegamass,_omegawidth;
  // use local parameters for the omega rather than from the particle data objects
  
};
  
}

// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of Tau2MesonPhotonDecayer.
  template <>
  struct BaseClassTrait<Herwig::Tau2MesonPhotonDecayer,1> {
    typedef Herwig::TauDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Tau2MesonPhotonDecayer>
    : public ClassTraitsBase<Herwig::Tau2MesonPhotonDecayer> {
    static string className() { return "/Herwig++/Tau2MesonPhotonDecayer"; }
    // Return the class name.
    static string library() { return "libHwTauDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "Tau2MesonPhotonDecayer.icc"

#endif /* HERWIG_Tau2MesonPhotonDecayer_H */
