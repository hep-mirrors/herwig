// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyVectorFormFactor class.
//

#include "BallZwickyVectorFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BallZwickyVectorFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

BallZwickyVectorFormFactor::~BallZwickyVectorFormFactor() {}

void BallZwickyVectorFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _Vr1 << _Vr2 << _A0r1 << _A0r2 << _A1r1 << _A1r2 << _A2r1 << _A2r2 << _T1r1
     << _T1r2 << _T2r1 << _T2r2 << _T3r1 << _T3r2 << _VmR2 << _Vmfit2 << _A0mR2 
     << _A0mfit2 << _A1mR2 << _A1mfit2 << _A2mR2 << _A2mfit2 << _T1mR2 << _T1mfit2 
     << _T2mR2 << _T2mfit2 << _T3mR2 << _T3mfit2;
}

void BallZwickyVectorFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _Vr1 >> _Vr2 >> _A0r1 >> _A0r2 >> _A1r1 >> _A1r2 >> _A2r1 >> _A2r2 >> _T1r1
     >> _T1r2 >> _T2r1 >> _T2r2 >> _T3r1 >> _T3r2 >> _VmR2 >> _Vmfit2 >> _A0mR2 
     >> _A0mfit2 >> _A1mR2 >> _A1mfit2 >> _A2mR2 >> _A2mfit2 >> _T1mR2 >> _T1mfit2 
     >> _T2mR2 >> _T2mfit2 >> _T3mR2 >> _T3mfit2;
}

ClassDescription<BallZwickyVectorFormFactor> BallZwickyVectorFormFactor::initBallZwickyVectorFormFactor;
// Definition of the static class description member.

void BallZwickyVectorFormFactor::Init() {

  static ClassDocumentation<BallZwickyVectorFormFactor> documentation
    ("The \\classname{BallZwickyVectorFormFactor} class implements the vector form"
     " factors of hep-ph/0412079 for the form-factor for the decay of a B-meson to a"
     " light pseudoscalar meson");

  static ParVector<BallZwickyVectorFormFactor,double> interfaceVr_1
    ("Vr_1",
     "The r_1 coefficient for the V form-factor",
     &BallZwickyVectorFormFactor::_Vr1,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceVr_2
    ("Vr_2",
     "The r_2 coefficient for the V form-factor",
     &BallZwickyVectorFormFactor::_Vr2,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA0r_1
    ("A0r_1",
     "The r_1 coefficient for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0r1,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA0r_2
    ("A0r_2",
     "The r_2 coefficient for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0r2,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA1r_1
    ("A1r_1",
     "The r_1 coefficient for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1r1,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA1r_2
    ("A1r_2",
     "The r_2 coefficient for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1r2,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA2r_1
    ("A2r_1",
     "The r_1 coefficient for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2r1,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA2r_2
    ("A2r_2",
     "The r_2 coefficient for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2r2,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT1r_1
    ("T1r_1",
     "The r_1 coefficient for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1r1,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT1r_2
    ("T1r_2",
     "The r_2 coefficient for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1r2,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT2r_1
    ("T2r_1",
     "The r_1 coefficient for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2r1,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT2r_2
    ("T2r_2",
     "The r_2 coefficient for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2r2,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT3r_1
    ("T3r_1",
     "The r_1 coefficient for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3r1,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT3r_2
    ("T3r_2",
     "The r_2 coefficient for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3r2,
     0, 0, 0, -100000., 1000000., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceVmR2
    ("VmR2",
     "The value of m_R^2 for the V form-factor",
     &BallZwickyVectorFormFactor::_VmR2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceVmfit2 
    ("Vmfit2 ",
     "The value of m_fit^2 for the V form-factor",
     &BallZwickyVectorFormFactor::_Vmfit2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA0mR2
    ("A0mR2",
     "The value of m_R^2 for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0mR2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA0mfit2 
    ("A0mfit2 ",
     "The value of m_fit^2 for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0mfit2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA1mR2
    ("A1mR2",
     "The value of m_R^2 for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1mR2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA1mfit2 
    ("A1mfit2 ",
     "The value of m_fit^2 for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1mfit2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA2mR2
    ("A2mR2",
     "The value of m_R^2 for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2mR2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA2mfit2 
    ("A2mfit2 ",
     "The value of m_fit^2 for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2mfit2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT1mR2
    ("T1mR2",
     "The value of m_R^2 for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1mR2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT1mfit2 
    ("T1mfit2 ",
     "The value of m_fit^2 for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1mfit2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT2mR2
    ("T2mR2",
     "The value of m_R^2 for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2mR2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT2mfit2 
    ("T2mfit2 ",
     "The value of m_fit^2 for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2mfit2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT3mR2
    ("T3mR2",
     "The value of m_R^2 for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3mR2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT3mfit2 
    ("T3mfit2 ",
     "The value of m_fit^2 for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3mfit2,
     0, 0, 0, -1000.*GeV2, 1000.*GeV2, false, false, true);

}

// form-factor for scalar to vector
void BallZwickyVectorFormFactor::ScalarVectorFormFactor(Energy2 q2,int mode,int id0,
							int id1,Energy m0, Energy m1,
							Complex & A0,Complex & A1,
							Complex & A2,Complex & V) const
{
  Complex ii(0.,1.);
  // the form-factors
  // A_0
  if(_A0mR2[mode]<0)
    {A0 = (_A0r1[mode]+_A0r2[mode]/(1.-q2/_A0mfit2[mode]))/(1.-q2/_A0mfit2[mode]);}
  else
    {A0 = _A0r1[mode]/(1.-q2/_A0mR2[mode])+_A0r2[mode]/(1.-q2/_A0mfit2[mode]);}
  // A_1
  if(_A1mR2[mode]<0)
    {A1 = (_A1r1[mode]+_A1r2[mode]/(1.-q2/_A1mfit2[mode]))/(1.-q2/_A1mfit2[mode]);}
  else
    {A1 = _A1r1[mode]/(1.-q2/_A1mR2[mode])+_A1r2[mode]/(1.-q2/_A1mfit2[mode]);}
  // A_2
  if(_A2mR2[mode]<0)
    {A2 = (_A1r2[mode]+_A2r2[mode]/(1.-q2/_A2mfit2[mode]))/(1.-q2/_A2mfit2[mode]);}
  else
    {A2 = _A1r2[mode]/(1.-q2/_A2mR2[mode])+_A2r2[mode]/(1.-q2/_A2mfit2[mode]);}
  // V
  if(_VmR2[mode]<0)
    {V = (_Vr2[mode]+_Vr2[mode]/(1.-q2/_Vmfit2[mode]))/(1.-q2/_Vmfit2[mode]);}
  else
    {V = _Vr2[mode]/(1.-q2/_VmR2[mode])+_Vr2[mode]/(1.-q2/_Vmfit2[mode]);}
}

void BallZwickyVectorFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,int mode,int id0,
							     int id1,Energy m0,Energy m1,
							     Complex & T1,Complex & T2,
							     Complex & T3) const
{
  // T_1
  if(_T1mR2[mode]<0)
    {T1 = (_T1r2[mode]+_T1r2[mode]/(1.-q2/_T1mfit2[mode]))/(1.-q2/_T1mfit2[mode]);}
  else
    {T1 = _T1r2[mode]/(1.-q2/_T1mR2[mode])+_T1r2[mode]/(1.-q2/_T1mfit2[mode]);}
  // T_2
  if(_T2mR2[mode]<0)
    {T2 = (_T2r2[mode]+_T2r2[mode]/(1.-q2/_T2mfit2[mode]))/(1.-q2/_T2mfit2[mode]);}
  else
    {T2 = _T2r2[mode]/(1.-q2/_T2mR2[mode])+_T2r2[mode]/(1.-q2/_T2mfit2[mode]);}
  // T_3
  if(_T3mR2[mode]<0)
    {T3 = (_T3r2[mode]+_T3r2[mode]/(1.-q2/_T3mfit2[mode]))/(1.-q2/_T3mfit2[mode]);}
  else
    {T3 = _T3r2[mode]/(1.-q2/_T3mR2[mode])+_T3r2[mode]/(1.-q2/_T3mfit2[mode]);}
  // convert for T_3tilde to T_3
  T3 = (m0*m0-m1*m1)/q2*(T3-T2);
}
}
