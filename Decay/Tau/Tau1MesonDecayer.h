// -*- C++ -*-
#ifndef HERWIG_Tau1MesonDecayer_H
#define HERWIG_Tau1MesonDecayer_H
//
// This is the declaration of the <!id>Tau1MesonDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class implements the hadronic currents for the decay of the tau
//  to a single meson either a pseudoscalar, i.e. pi+/- and K+/-,
//  or a vector, rho or K*.
//
//  In practice it is only used for the decay into pseudoscalar mesons.
//  The <!id>Tau2MesonDecayer<!!id> which includes the decay products of the vector
//  mesons and the correct correlations should be used for rho and K*
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TauDecayerBase.html">TauDecayerBase.h</a>,
// <a href="Tau2MesonDecayer.html">Tau2MesonDecayer.h</a>.
// 
//  Author: Peter Richardson
//

#include "TauDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

class Tau1MesonDecayer: public TauDecayerBase {
  
public:
  
  inline Tau1MesonDecayer();
  inline Tau1MesonDecayer(const Tau1MesonDecayer &);
  virtual ~Tau1MesonDecayer();
  // Standard ctors and dtor.
  
public:
  
  vector<LorentzPolarizationVector> 
  hadronCurrent(bool vertex, const int imode, const int ichan,
		const Particle &, const ParticleVector &) const;
  // the hadronic current for this decay mode
  
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
  
  static ClassDescription<Tau1MesonDecayer> initTau1MesonDecayer;
  // Describe a concrete class with persistent data.
  
  Tau1MesonDecayer & operator=(const Tau1MesonDecayer &);
  // Private and non-existent assignment operator.
  
private:
  
  Energy _fpi,_fk;
  // pseudoscalar decay constants
  
  Energy2 _grho,_gkstar;
  // vector decay constants
  
  mutable double _piwgt,_kwgt;
  // maximum weight for scalar decays
  
  mutable double _rhowgt,_kstarwgt;
  // maximum weight for vector decays

};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of Tau1MesonDecayer.
  template <>
  struct BaseClassTrait<Herwig::Tau1MesonDecayer,1> {
    typedef Herwig::TauDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Tau1MesonDecayer>
    : public ClassTraitsBase<Herwig::Tau1MesonDecayer> {
    static string className() { return "/Herwig++/Tau1MesonDecayer"; }
    // Return the class name.
    static string library() { return "libHwTauDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "Tau1MesonDecayer.icc"

#endif /* HERWIG_Tau1MesonDecayer_H */
