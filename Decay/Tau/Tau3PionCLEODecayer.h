// -*- C++ -*-
#ifndef HERWIG_Tau3PionCLEODecayer_H
#define HERWIG_Tau3PionCLEODecayer_H
//
// This is the declaration of the <!id>Tau3PionCLEODecayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  The <!id>Tau3PionCLEODecayer<!!id> class implements
//  the decay of the tau to three pions
//  and a neutrino using the currents from CLEO Phys. Rev. D 61,012002. This is
//  a model including two rho mesons in both s and p wave, a sigma the F_2 and f_0.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="Tau3MesonDecayerBase.html">Tau3MesonDecayerBase.h</a>,
// <a href="TauDecayerBase.html">TauDecayerBase.h</a>.
// 

#include "Tau3MesonDecayerBase.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Decay/ThreeBodyIntegrator.h"
#include "Herwig++/Utilities/Kinematics.h"
// #include "Tau3PionCLEODecayer.fh"
// #include "Tau3PionCLEODecayer.xh"

namespace Herwig {
using namespace Herwig;

class Tau3PionCLEODecayer: public Tau3MesonDecayerBase {
    
public:
  
  inline Tau3PionCLEODecayer();
  inline Tau3PionCLEODecayer(const Tau3PionCLEODecayer &);
  virtual ~Tau3PionCLEODecayer();
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
  
  inline double a1MatrixElement(int,Energy2,Energy2,Energy2,Energy2);
  // the matrix element for the running a_1 width
  
public:
  
  virtual void calculateFormFactors(const int,const int,
				      Energy2,Energy2,Energy2,Energy2,
				    complex<double>&,complex<double>&,complex<double>&,
				    complex<double>&,complex<double>&) const;
  // calculate the form factor for the current

  void CLEOFormFactor(int,int,Energy2 q2,Energy2 s1, Energy2 s2, Energy2 s3,
		      Complex & F1, Complex & F2, Complex & F3) const;
  // the CLEO form factors
  
  virtual bool acceptMode(int) const;
  // can a particular decayer handle this type of mode
  
  virtual int phaseSpaceMode(int) const;
  // mapping of the mode to the phase space 
  
private:
  
  inline Energy a1width(Energy2) const ;
  // return the a1 running width
  
  inline void inita1width(int);
  // initialize the a_1 running width

  inline Complex a1BreitWigner(Energy2 q2) const;
  
private:
  
  static ClassDescription<Tau3PionCLEODecayer> initTau3PionCLEODecayer;
  // Describe a concrete class with persistent data.
  
  Tau3PionCLEODecayer & operator=(const Tau3PionCLEODecayer &);
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
  Energy _a1mass,_a1width;
  // mass and width of the a_1
  Energy _mKstar,_mK;
  double _gammk;
  // parameters for the K K* contribution to the a_1 running width

  Energy _fpi; InvEnergy _fact;
  // pion decya constant

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
  
  mutable vector<bool> _onechan,_threechan;
  mutable vector<double> _onewgts,_threewgts;
  mutable double _onemax,_threemax;
  // parameters for the multi-channel integration
  
  vector<Energy> _a1runwidth;
  vector<Energy> _a1runq2;
  Interpolator *_a1runinter;
  bool _initializea1;
  // the a1 width and q2 for the running width calculation
};
  
}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of Tau3PionCLEODecayer.
  template <>
  struct BaseClassTrait<Herwig::Tau3PionCLEODecayer,1> {
    typedef Herwig::Tau3MesonDecayerBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::Tau3PionCLEODecayer>
    : public ClassTraitsBase<Herwig::Tau3PionCLEODecayer> {
    static string className() { return "/Herwig++/Tau3PionCLEODecayer"; }
    // Return the class name.
    static string library() { return "libHwTauDecay.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

// definitions of the functions to be integrated to give the running
#include "CLHEP/GenericFunctions/AbsFunction.hh"
// functions to return the matrix element for the a1 decay to be
// integrated to give the a1 running width
namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

  // one chraged 2 neutral channel
class Tau3PionCLEOa1MatrixElement: public Genfun::AbsFunction {

  FUNCTION_OBJECT_DEF(Tau3PionCLEOa1MatrixElement);

public:
  
  // Constructor
  Tau3PionCLEOa1MatrixElement(int,Ptr<Herwig::Tau3PionCLEODecayer>::pointer);
  
  // Destructor
  virtual ~Tau3PionCLEOa1MatrixElement();
  
  virtual unsigned int dimensionality() const ;     
  
  // Copy constructor
  Tau3PionCLEOa1MatrixElement(const Tau3PionCLEOa1MatrixElement &right);
  
  // Retreive function value
  virtual double operator ()(double) const {return 0.;}
  virtual double operator ()(const Argument & a) const ;
  
  inline void setQ2(Energy2);
  
  
private:
  
  // It is illegal to assign a function
  const Tau3PionCLEOa1MatrixElement & 
  operator=(const Tau3PionCLEOa1MatrixElement &right);
  
private:
  // the decayer
  Ptr<Herwig::Tau3PionCLEODecayer>::pointer _decayer;
  // the integer for the mode
  int _mode;
};

}
#include "Tau3PionCLEODecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Tau3PionCLEODecayer.tcc"
#endif


#endif /* HERWIG_Tau3PionCLEODecayer_H */
