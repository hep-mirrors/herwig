// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyScalarFormFactor class.
//

#include "BallZwickyScalarFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BallZwickyScalarFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

BallZwickyScalarFormFactor::~BallZwickyScalarFormFactor() {}

void BallZwickyScalarFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _r10 << _r20 << _r1plus << _r2plus << _r1T << _r2T << _m120 << _mfit20 
     << _m12plus << _mfit2plus << _m12T << _mfit2T;
}

void BallZwickyScalarFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _r10 >> _r20 >> _r1plus >> _r2plus >> _r1T >> _r2T >> _m120 >> _mfit20 
     >> _m12plus >> _mfit2plus >> _m12T >> _mfit2T;
}

ClassDescription<BallZwickyScalarFormFactor> BallZwickyScalarFormFactor::initBallZwickyScalarFormFactor;
// Definition of the static class description member.

void BallZwickyScalarFormFactor::Init() {

  static ClassDocumentation<BallZwickyScalarFormFactor> documentation
    ("The \\classname{BallZwickyScalarFormFactor} class implements the form-factors"
     " of hep-ph/0406232 for the form-factor for the decay of a B-meson to a"
     " light pseudoscalar meson");

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_10
    ("r_10",
     "The r_1 coefficient for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_r10,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_20
    ("r_20",
     "The r_2 coefficient for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_r20,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_1plus
    ("r_1plus",
     "The r_1 coefficient for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_r1plus,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_2plus
    ("r_2plus",
     "The r_2 coefficient for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_r2plus,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_1T
    ("r_1T",
     "The r_1 coefficient for the f_T form-factor",
     &BallZwickyScalarFormFactor::_r1T,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_2T
    ("r_2T",
     "The r_2 coefficient for the f_T form-factor",
     &BallZwickyScalarFormFactor::_r2T,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem120
    ("m_120",
     "The value of m_1^2 for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_m120,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit20
    ("m_120",
     "The value of m_fit^2 for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_mfit20,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem12plus
    ("m_12plus",
     "The value of m_1^2 for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_m12plus,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit2plus
    ("m_12plus",
     "The value of m_fit^2 for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_mfit2plus,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem12T
    ("m_12T",
     "The value of m_1^2 for the f_T form-factor",
     &BallZwickyScalarFormFactor::_m12T,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit2T
    ("m_12T",
     "The value of m_fit^2 for the f_T form-factor",
     &BallZwickyScalarFormFactor::_mfit2T,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);
}


// form-factor for scalar to scalar
void BallZwickyScalarFormFactor::ScalarScalarFormFactor(Energy2 q2, int mode,int id0,
							int id1,Energy m0, Energy m1,
						 	Complex & f0, Complex & fp) const
{
  Complex ii(0.,1.);
  // the F_0 form-factor
  if(_m120[mode]<0)
    {f0=_r20[mode]/(1.-q2/_mfit20[mode]);}
  else if(_mfit20[mode]<0)
    {f0=(_r10[mode]+_r20[mode]/(1.-q2/_m120[mode]))/(1.-q2/_m120[mode]);}
  else
    {f0=_r10[mode]/(1.-q2/_m120[mode])+_r20[mode]/(1.-_mfit20[mode]);}
  // the F_1 form-factor
  if(_m12plus[mode]<0)
    {fp = _r2plus[mode]/(1.-q2/_mfit2plus[mode]);}
  else if(_mfit2plus[mode]<0)
    {fp = (_r1plus[mode]+_r2plus[mode]/(1.-q2/_m12plus[mode]))/(1.-q2/_m12plus[mode]);}
  else
    {fp =_r1plus[mode]/(1.-q2/_m12plus[mode])+_r2plus[mode]/(1.-q2/_mfit2plus[mode]);}
}

void BallZwickyScalarFormFactor::ScalarScalarSigmaFormFactorEnergy2(Energy2 q2,int mode,
								    int id0,int id1,
								    Energy m0, Energy m1,
								    Complex & fT) const
{
  // the F_T form-factor
  if(_m12T[mode]<0)
    {fT = _r2T[mode]/(1.-q2/_mfit2T[mode]);}
  else if(_mfit2T[mode]<0)
    {fT = (_r1T[mode]+_r2T[mode]/(1.-q2/_m12T[mode]))/(1.-q2/_m12T[mode]);}
  else
    {fT =_r1T[mode]/(1.-q2/_m12T[mode])+_r2T[mode]/(1.-q2/_mfit2T[mode]);}
}
}
