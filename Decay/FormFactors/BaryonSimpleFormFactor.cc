// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonSimpleFormFactor class.
//

#include "BaryonSimpleFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonSimpleFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

BaryonSimpleFormFactor::~BaryonSimpleFormFactor() {}

void BaryonSimpleFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _gA << _alphaD << _etaV << _etaA << _rhoE << _rhoM 
     << _f1 << _f2 << _g1 << _g2;}

void BaryonSimpleFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _gA >> _alphaD >> _etaV >> _etaA >> _rhoE >> _rhoM 
     >> _f1 >> _f2 >> _g1 >> _g2;}

ClassDescription<BaryonSimpleFormFactor> BaryonSimpleFormFactor::initBaryonSimpleFormFactor;
// Definition of the static class description member.

void BaryonSimpleFormFactor::Init() {

  static ClassDocumentation<BaryonSimpleFormFactor> documentation
    ("The \\classname{BaryonSimpleFormFactor} class implements the"
     " quark model calculation of the form-factors from PRD25, 206");

  static Parameter<BaryonSimpleFormFactor,double> interfaceg_A
    ("g_A",
     "The axial-vector coupling in neutron beta decay.",
     &BaryonSimpleFormFactor::_gA, 1.25, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacealpha_D
    ("alpha_D",
     "SU(3) breaking parameter which is the ratio D/(D+F). ",
     &BaryonSimpleFormFactor::_alphaD, 0.6, 0.0, 1.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfaceeta_V
    ("eta_V",
     "The eta_V SU(3) breaking parameter",
     &BaryonSimpleFormFactor::_etaV, .97, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfaceeta_A
    ("eta_A",
     "The eta_A SU(3) breaking parameter",
     &BaryonSimpleFormFactor::_etaA, 1.08, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacerho_E
    ("rho_E",
     "The SU(3) breaking parameter for the electric dipole moment.",
     &BaryonSimpleFormFactor::_rhoE, 0.86, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacerho_M
    ("rho_M",
     "The SU(3) breaking parameter for the magentic dipole moment.",
     &BaryonSimpleFormFactor::_rhoM, 0.094, 0.0, 10.0,
     false, false, true);
}

void BaryonSimpleFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc, int id0,int id1,Energy,Energy,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  f1v =  _f1[iloc];
  f1a = -_g1[iloc];
  f2v = -_f2[iloc];
  f3v =  _g2[iloc];
  f3v = 0.;
  f3a = 0.;
}

void BaryonSimpleFormFactor::dataBaseOutput(ofstream&) {;}
}
