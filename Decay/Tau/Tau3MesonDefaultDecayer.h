// -*- C++ -*-
#ifndef THEPEG_Tau3MesonDefaultDecayer_H
#define THEPEG_Tau3MesonDefaultDecayer_H
//
// This is the declaration of the <!id>Tau3MesonDefaultDecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  This class implements the decays for the three meson decay modes of the tau.
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
//  using the currents in TAUOLA
//
// CLASSDOC SUBSECTION See also:
//
// <a href="TauDecayerBase.html">TauDecayerBase.h</a>,
// <a href="Tau3MesonDecayerBase.html">Tau3MesonDecayerBase.h</a>.
// 
//  Author: Peter Richardson
//

#include "Tau3MesonDecayerBase.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Decay/ThreeBodyIntegrator.h"
#include "Herwig++/Utilities/Kinematics.h"
// #include "Tau3MesonDefaultDecayer.fh"
// #include "Tau3MesonDefaultDecayer.xh"

namespace Herwig {

using namespace ThePEG;

class Tau3MesonDefaultDecayer: public Tau3MesonDecayerBase {
  
public:
    
  inline Tau3MesonDefaultDecayer();
  inline Tau3MesonDefaultDecayer(const Tau3MesonDefaultDecayer &);
  virtual ~Tau3MesonDefaultDecayer();
  // Standard ctors and dtor.
  
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
  
public:    
  
  inline double a1MatrixElement(Energy2,Energy2,Energy2,Energy2);
  // the matrix element for the a1 decay to calculate the running width
  
public:
  
  virtual void calculateFormFactors(const int,const int,
				    Energy2,Energy2,Energy2,Energy2,
				    Complex&,Complex&,Complex&,
				    Complex&,Complex&) const;
  // calculate the form factor for the current
  
  virtual bool acceptMode(int) const;
  // can a particular decayer handle this type of mode
  
  virtual int phaseSpaceMode(int) const;
  // mapping of the mode to the phase space 
  
private:
  
  inline Complex BrhoF123(Energy2,int) const;
  inline Complex BrhoF5(Energy2,int) const;
  inline Complex BKstarF123(Energy2,int) const;
  inline Complex BKstarF5(Energy2,int) const;
  // the Breit wigner factors
  
  inline Complex FKrho(Energy2,Energy2,int) const;
  // mixed Breit Wigner
  
  inline Complex a1BreitWigner(Energy2) const;
  // a_1 breit wigner
  
  inline Complex K1BreitWigner(Energy2) const;
  // K_1 breit Wigner
  
  inline Energy a1width(Energy2) const ;
  // return the a1 running width
  
  inline void inita1width(int);
  // initialize the a_1 running width

  inline Complex rhoKBreitWigner(Energy2,unsigned int,unsigned int) const;
  // breit wigners for the rho and K*
  
private:
  
  static ClassDescription<Tau3MesonDefaultDecayer> initTau3MesonDefaultDecayer;
    // Describe a concrete class with persistent data.
  
  Tau3MesonDefaultDecayer & operator=(const Tau3MesonDefaultDecayer &);
  // Private and non-existent assignment operator.
  
private:
  
  vector<double> _rhoF123wgts,_KstarF123wgts;
  // parameters for the rho and K* Breit-Wigner in the F1,F2,F3 form factors
  
  vector<double> _rhoF5wgts,_KstarF5wgts;
  // parameters for the rho and K*Breit-Wigner in the F5 form factor
  
  double _rhoKstarwgt;
  // the relative weight of the rho and Kstar where needed
  
  vector<Energy> _a1runwidth;
  vector<Energy> _a1runq2;
  Interpolator *_a1runinter;
  bool _initializea1;
  // the a1 width and q2 for the running width calculation
  
  Energy _a1mass,_a1width;
  // the masses and width of the a_1 resonances
  Energy _K1mass,_K1width;
  // the masses and width of the K_1 resonances
  InvEnergy _sinfact,_cosfact;
  Energy _fpi;
  // constants
  Energy _mpi,_mK;
  // masses for the running widths in the Breit wigner

  
  bool _rhoparameters;
  vector<Energy> _rhoF123masses,_rhoF5masses,_rhoF123widths,_rhoF5widths;
  // use local values of the rho masses and widths
  
  bool _Kstarparameters;
  vector<Energy> _KstarF123masses,_KstarF5masses,_KstarF123widths,_KstarF5widths;
  // use local values of the K* resonances masses and widths
  
  bool _a1parameters;
  // use local values of the a_1 parameters
  
  bool _K1parameters;
  // use local values of the K_1 parameters
  
  mutable vector<bool> _pimpimpipchan,_pi0pi0pimchan,_KmpimKpchan,
    _K0pimK0chan,_Kmpi0K0chan,_pi0pi0Kmchan,_Kmpimpipchan,_pimK0pi0chan,
    _pimpi0etachan;
  mutable vector<double> _pimpimpipwgts,_pi0pi0pimwgts,_KmpimKpwgts,
    _K0pimK0wgts,_Kmpi0K0wgts,_pi0pi0Kmwgts,_Kmpimpipwgts,_pimK0pi0wgts,
    _pimpi0etawgts;
  mutable double  _pimpimpipmax,_pi0pi0pimmax,_KmpimKpmax,
    _K0pimK0max,_Kmpi0K0max,_pi0pi0Kmmax,_Kmpimpipmax,_pimK0pi0max,
    _pimpi0etamax;
  // weights for the integration
};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of Tau3MesonDefaultDecayer.
  template <>
  struct BaseClassTrait<Herwig::Tau3MesonDefaultDecayer,1> {
    typedef Herwig::Tau3MesonDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Tau3MesonDefaultDecayer>
    : public ClassTraitsBase<Herwig::Tau3MesonDefaultDecayer> {
    static string className() { return "/Herwig++/Tau3MesonDefaultDecayer"; }
    // Return the class name.
    static string library() { return "libHwTauDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}




// definitions of the functions to be integrated to give the running
#include "CLHEP/GenericFunctions/AbsFunction.hh"
// function to return the matrix element for the a1 decay to be
// integrated to give the a1 running width
namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

class tau3MesonDefaulta1MatrixElement : public Genfun::AbsFunction {
  
FUNCTION_OBJECT_DEF(tau3MesonDefaulta1MatrixElement)
    
  public:

// Constructor
  tau3MesonDefaulta1MatrixElement(Ptr<Herwig::Tau3MesonDefaultDecayer>::pointer);
  
  // Destructor
  virtual ~tau3MesonDefaulta1MatrixElement();
  
  virtual unsigned int dimensionality() const ;     
  
  // Copy constructor
  tau3MesonDefaulta1MatrixElement(const tau3MesonDefaulta1MatrixElement &right);
  
  // Retreive function value
  virtual double operator ()(double) const {return 0.;}
  virtual double operator ()(const Argument & a) const ;
  
  inline void setQ2(Energy2);
  
  
private:
  
  // It is illegal to assign a function
  const tau3MesonDefaulta1MatrixElement & 
  operator=(const tau3MesonDefaulta1MatrixElement &right);
  
private:
  // the decayer
  Ptr<Herwig::Tau3MesonDefaultDecayer>::pointer _decayer;
};
}












#include "Tau3MesonDefaultDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Tau3MesonDefaultDecayer.tcc"
#endif

#endif /* THEPEG_Tau3MesonDefaultDecayer_H */
