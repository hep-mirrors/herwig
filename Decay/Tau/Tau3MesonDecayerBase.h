// -*- C++ -*-
#ifndef HERWIG_Tau3MesonDecayerBase_H
#define HERWIG_Tau3MesonDecayerBase_H
//
// This is the declaration of the <!id>Tau3MesonDecayerBase<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This is the base class for the three meson decay modes of the tau.
//  It implements the following decay modes of the tau.
//
//    pi-  pi-    pi+
//    pi0  pi0    pi-
//    K-   pi-    K+
//    K0   pi-    Kbar0
//    K-   pi0    K0
//    pi0  pi0    K-
//    K-   pi-    pi+
//    pi-  Kbar0  pi0.
//    pi-  pi0    eta
// 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:TauDecayerBase.html">.h</a>.
// 
// Author: Peter Richardson

#include "TauDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

class Tau3MesonDecayerBase: public TauDecayerBase {
  
public:
  
  inline Tau3MesonDecayerBase();
  inline Tau3MesonDecayerBase(const Tau3MesonDecayerBase &);
  virtual ~Tau3MesonDecayerBase();
  // Standard ctors and dtor.
  
public:
  
  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.
  
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.
  
  vector<LorentzPolarizationVector>  hadronCurrent(bool vertex, const int,
						   const int, const Particle &,
						   const ParticleVector &) const;
  // the hadronic current for this decay mode

  virtual void calculateFormFactors(const int, const int,
				    Energy2,Energy2,Energy2,Energy2,
				    Complex&,Complex&,Complex&,
				    Complex&,Complex&) const;
  // calculate the form factor for the current
  
  virtual bool acceptMode(int) const;
  // can a particular decayer handle this type of mode
  
  virtual int phaseSpaceMode(int) const;
  // mapping of the mode to the phase space 
  
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
  
  static ClassDescription<Tau3MesonDecayerBase> initTau3MesonDecayerBase;
  // Describe a concrete class with persistent data.
  
  Tau3MesonDecayerBase & operator=(const Tau3MesonDecayerBase &);
  // Private and non-existent assignment operator.
  
};
  
}

// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of Tau3MesonDecayerBase.
  template <>
  struct BaseClassTrait<Herwig::Tau3MesonDecayerBase,1> {
    typedef Herwig::TauDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Tau3MesonDecayerBase>
    : public ClassTraitsBase<Herwig::Tau3MesonDecayerBase> {
    static string className() { return "/Herwig/Tau3MesonDecayerBase"; }
    // Return the class name.
    static string library() { return "libHwTauDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "Tau3MesonDecayerBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Tau3MesonDecayerBase.tcc"
#endif

#endif /* HERWIG_Tau3MesonDecayerBase_H */
