// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LightBaryonQuarkModelFormFactor class.
//

#include "LightBaryonQuarkModelFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LightBaryonQuarkModelFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig{
using namespace ThePEG;

LightBaryonQuarkModelFormFactor::~LightBaryonQuarkModelFormFactor() {}

void LightBaryonQuarkModelFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _f1 << _f2 << _g1 << _g2 << _Lambdaf1 << _Lambdaf2 << _Lambdag1 << _Lambdag2;}

void LightBaryonQuarkModelFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _f1 >> _f2 >> _g1 >> _g2 >> _Lambdaf1 >> _Lambdaf2 >> _Lambdag1 >> _Lambdag2;}

ClassDescription<LightBaryonQuarkModelFormFactor> LightBaryonQuarkModelFormFactor::initLightBaryonQuarkModelFormFactor;
// Definition of the static class description member.

void LightBaryonQuarkModelFormFactor::Init() {

  static ClassDocumentation<LightBaryonQuarkModelFormFactor> documentation
    ("The \\classname{LightBaryonQuarkModelFormFactor} class implements"
     " the quark model calculation of hep-ph/9409272 for the form-factors"
     " for the light quarks");

  static ParVector<LightBaryonQuarkModelFormFactor,double> interfacef1
    ("f1",
     "The form-factor f1 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_f1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,double> interfaceg1
    ("g1",
     "The form-factor g1 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_g1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,InvEnergy> interfacef2
    ("f2",
     "The form-factor f2 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_f2,
     1./GeV, 0, 0, -10./GeV, 10./GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,InvEnergy> interfaceg2
    ("g2",
     "The form-factor g2 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_g2,
     1./GeV, 0, 0, -10./GeV, 10./GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdaf1
    ("Lambdaf1",
     "The first mass for the energy dependence of the f1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdaf1,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdaf2
    ("Lambdaf2",
     "The second mass for the energy dependence of the f1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdaf2,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdag1
    ("Lambdag1",
     "The first mass for the energy dependence of the g1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdag1,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdag2
    ("Lambdag2",
     "The second mass for the energy dependence of the g1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdag2,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);
}

// form factor for spin-1/2 to spin-1/2
void LightBaryonQuarkModelFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int mode,int id0, int id1, Energy m0, Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  // f_3 is zero
  f3v=0.; f3a=0.;
  // energy dependence of the f1 and g1 factors
  InvEnergy2 lam1(1./(_Lambdaf1[mode]*_Lambdaf1[mode]));
  InvEnergy2 lam2(1./(_Lambdaf2[mode]*_Lambdaf2[mode]));
  f1v= _f1[mode]/(1.-q2*lam1+q2*q2*lam2*lam2);
  lam1 = 1./(_Lambdag1[mode]*_Lambdag1[mode]);
  lam2 = 1./(_Lambdag2[mode]*_Lambdag2[mode]);
  f1a=-_g1[mode]/(1.-q2*lam1+q2*q2*lam2*lam2);
  // the f2 and g2 factors
  f2v =-(m0+m1)*_f2[mode];
  f2a = (m0+m1)*_g2[mode]; 
}

}
