// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WSBFormFactor class.
//

#include "WSBFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WSBFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

WSBFormFactor::~WSBFormFactor() {}

void WSBFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _F0 << _V << _A0 << _A1 << _A2 << _mS0 << _mS1 << _mV0 << _mV1;
}

void WSBFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _F0 >> _V >> _A0 >> _A1 >> _A2 >> _mS0 >> _mS1 >> _mV0 >> _mV1;
}

ClassDescription<WSBFormFactor> WSBFormFactor::initWSBFormFactor;
// Definition of the static class description member.

void WSBFormFactor::Init() {

  static ClassDocumentation<WSBFormFactor> documentation
    ("The\\classname{WSBFormFactor} class is the implementation of the form-factors of "
     "Z.Phys.C29,637.");

  static ParVector<WSBFormFactor,double> interfaceV
    ("V",
     "The form-factor V at zero q^2",
     &WSBFormFactor::_V,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA0
    ("A0",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_A0,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA1
    ("A1",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_A1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA2
    ("A2",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_A2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,Energy> interfaceScalarMass
    ("ScalarMass",
     "The scalar mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mS0,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoScalarMass
    ("PseudoScalarMass",
     "The pseudoscalar mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mS1,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfaceVectorMass
    ("VectorMass",
     "The vector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV0,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoVectorMass
    ("PseudoVectorMass",
     "The pseudovector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV1,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);
}

// form-factor for scalar to scalar
void WSBFormFactor::ScalarScalarFormFactor(Energy2 q2, int id0, int mode,int id1,
					   Energy m0, Energy m1,Complex & f0,
					   Complex & fp) const
{
  f0 = _F0[mode]/(1.-q2/_mS1[mode]/_mS1[mode]);
  fp = _F0[mode]/(1.-q2/_mV0[mode]/_mV0[mode]);
}
void WSBFormFactor::ScalarVectorFormFactor(Energy2 q2, int mode, int id0, int id1, 
					   Energy mY, Energy mX,Complex & A0,
					   Complex & A1,Complex & A2,Complex & V) const
{
  Complex ii(0.,1.);
  A0 = -_A0[mode]/(1.-q2/_mS0[mode]/_mS0[mode]);
  A1 = -_A1[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  A2 = -_A2[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  V  =   _V[mode]/(1.-q2/_mV0[mode]/_mV0[mode]);
}

}
