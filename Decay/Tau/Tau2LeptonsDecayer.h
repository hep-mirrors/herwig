// -*- C++ -*-
#ifndef HERWIG_Tau2LeptonsDecayer_H
#define HERWIG_Tau2LeptonsDecayer_H
//
// This is the declaration of the <!id>Tau2LeptonsDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class implements the hadronic currents for the decay of the tau
//  to leptons, obviously in this case the hadronic current is actually
//  a leptonic current. In inherits from the <!id>TauDecayerBase<!!id> class,
//  implements the hadronic current and specifies the phase space integration.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:TauDecayerBase.html">.h</a>.
// 
//  Author: Peter Richardson 
//

#include "TauDecayerBase.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

class Tau2LeptonsDecayer: public TauDecayerBase {
    
public:
  
  inline Tau2LeptonsDecayer();
  inline Tau2LeptonsDecayer(const Tau2LeptonsDecayer &);
  virtual ~Tau2LeptonsDecayer();
  // Standard ctors and dtor.
  
public:
  
  vector<LorentzPolarizationVector>  hadronCurrent(bool vertex, const int,
						   const int, const Particle &,
						   const ParticleVector &) const;
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
  
  static ClassDescription<Tau2LeptonsDecayer> initTau2LeptonsDecayer;
  // Describe a concrete class with persistent data.
  
  Tau2LeptonsDecayer & operator=(const Tau2LeptonsDecayer &);
  // Private and non-existent assignment operator.
  
private:
  
  double _electronwgt,_muonwgt;
  // maximum weights for the different channels
  vector<bool> _elon,_muon;
  // which channels to use
  vector<double> _elwgt,_muwgt;
  // weights for the channels
};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of Tau2LeptonsDecayer.
  template <>
  struct BaseClassTrait<Herwig::Tau2LeptonsDecayer,1> {
    typedef Herwig::TauDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Tau2LeptonsDecayer>
    : public ClassTraitsBase<Herwig::Tau2LeptonsDecayer> {
    static string className() { return "/Herwig++/Tau2LeptonsDecayer"; }
    // Return the class name.
    static string library() { return "libHwnewDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "Tau2LeptonsDecayer.icc"

#endif /* HERWIG_Tau2LeptonsDecayer_H */
