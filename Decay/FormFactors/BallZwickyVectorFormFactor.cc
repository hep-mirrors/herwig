// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyVectorFormFactor class.
//

#include "BallZwickyVectorFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BallZwickyVectorFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

BallZwickyVectorFormFactor::BallZwickyVectorFormFactor() 
{
  double ort(1./sqrt(2.));
  // parameters for the different form-factors
  // B to rho
  // B+ to rho0
  _Vr1.push_back(1.045*ort);_Vr2.push_back(-0.721*ort);
  _VmR2.push_back(5.32*5.32*GeV2);_Vmfit2.push_back(38.34*GeV2);
  _A0r1.push_back(1.527*ort);_A0r2.push_back(-1.220*ort);
  _A0mR2.push_back(5.28*5.28*GeV2);_A0mfit2.push_back(33.36*GeV2);
  _A1r1.push_back(0.240*ort);_A1r2.push_back(0.);
  _A1mR2.push_back(-1.0*GeV2);_A1mfit2.push_back(37.51*GeV2);
  _A2r1.push_back(0.009*ort);_A2r2.push_back(0.212*ort);
  _A2mR2.push_back(-1.0*GeV2);_A2mfit2.push_back(40.82*GeV2);
  _T1r1.push_back(0.897*ort);_T1r2.push_back(-0.629*ort);
  _T1mR2.push_back(5.32*5.32*GeV2);_T1mfit2.push_back(38.04*GeV2);
  _T2r1.push_back(0.267*ort);_T2r2.push_back(0.);
  _T2mR2.push_back(-1.0*GeV2);_T2mfit2.push_back(38.59*GeV2);
  _T3r1.push_back(0.022*ort);_T3r2.push_back(0.245*ort);
  _T3mR2.push_back(-1.0*GeV2);_T3mfit2.push_back(40.88*GeV2);
  addFormFactor(521,113,1,2,-5,-2);
  // B0 to rho+
  _Vr1.push_back(1.045);_Vr2.push_back(-0.721);
  _VmR2.push_back(5.32*5.32*GeV2);_Vmfit2.push_back(38.34*GeV2);
  _A0r1.push_back(1.527);_A0r2.push_back(-1.220);
  _A0mR2.push_back(5.28*5.28*GeV2);_A0mfit2.push_back(33.36*GeV2);
  _A1r1.push_back(0.240);_A1r2.push_back(0.);
  _A1mR2.push_back(-1.0*GeV2);_A1mfit2.push_back(37.51*GeV2);
  _A2r1.push_back(0.009);_A2r2.push_back(0.212);
  _A2mR2.push_back(-1.0*GeV2);_A2mfit2.push_back(40.82*GeV2);
  _T1r1.push_back(0.897);_T1r2.push_back(-0.629);
  _T1mR2.push_back(5.32*5.32*GeV2);_T1mfit2.push_back(38.04*GeV2);
  _T2r1.push_back(0.267);_T2r2.push_back(0.);
  _T2mR2.push_back(-1.0*GeV2);_T2mfit2.push_back(38.59*GeV2);
  _T3r1.push_back(0.022);_T3r2.push_back(0.245);
  _T3mR2.push_back(-1.0*GeV2);_T3mfit2.push_back(40.88*GeV2);
  addFormFactor(511,-213,1,1,-5,-2);
  // B to K*
  // B+ to K*0
  _Vr1.push_back(0.923);_Vr2.push_back(-0.511);
  _VmR2.push_back(5.32*5.32*GeV2);_Vmfit2.push_back(49.40*GeV2);
  _A0r1.push_back(1.364);_A0r2.push_back(-0.990);
  _A0mR2.push_back(5.28*5.28*GeV2);_A0mfit2.push_back(36.78*GeV2);
  _A1r1.push_back(0.290);_A1r2.push_back(0.0);
  _A1mR2.push_back(-1.0*GeV2);_A1mfit2.push_back(40.38*GeV2);
  _A2r1.push_back(-0.084);_A2r2.push_back(0.342);
  _A2mR2.push_back(-1.0*GeV2);_A2mfit2.push_back(52.00*GeV2);
  _T1r1.push_back(0.823);_T1r2.push_back(-0.491);
  _T1mR2.push_back(5.32*5.32*GeV2);_T1mfit2.push_back(46.31*GeV2);
  _T2r1.push_back(0.333);_T2r2.push_back(0.);
  _T2mR2.push_back(-1.0*GeV2);_T2mfit2.push_back(41.41*GeV2);
  _T3r1.push_back(-0.036);_T3r2.push_back(0.369);
  _T3mR2.push_back(-1.0*GeV2);_T3mfit2.push_back(48.10*GeV2);
  addFormFactor(521,323,1,2,-5,-3);
  // B0 to K*+
  _Vr1.push_back(0.923);_Vr2.push_back(-0.511);
  _VmR2.push_back(5.32*5.32*GeV2);_Vmfit2.push_back(49.40*GeV2);
  _A0r1.push_back(1.364);_A0r2.push_back(-0.990);
  _A0mR2.push_back(5.28*5.28*GeV2);_A0mfit2.push_back(36.78*GeV2);
  _A1r1.push_back(0.290);_A1r2.push_back(0.);
  _A1mR2.push_back(-1.0*GeV2);_A1mfit2.push_back(40.38*GeV2);
  _A2r1.push_back(-0.084);_A2r2.push_back(0.342);
  _A2mR2.push_back(-1.0*GeV2);_A2mfit2.push_back(52.00*GeV2);
  _T1r1.push_back(0.823);_T1r2.push_back(-0.491);
  _T1mR2.push_back(5.32*5.32*GeV2);_T1mfit2.push_back(46.31*GeV2);
  _T2r1.push_back(0.333);_T2r2.push_back(0.);
  _T2mR2.push_back(-1.0*GeV2);_T2mfit2.push_back(41.41*GeV2);
  _T3r1.push_back(-0.036);_T3r2.push_back(0.369);
  _T3mR2.push_back(-1.0*GeV2);_T3mfit2.push_back(48.10*GeV2);
  addFormFactor(511,313,1,2,-5,-3);
  // B to omega
  _Vr1.push_back(1.006*ort);_Vr2.push_back(-0.713*ort);
  _VmR2.push_back(5.32*5.32*GeV2);_Vmfit2.push_back(37.45*GeV2);
  _A0r1.push_back(1.321*ort);_A0r2.push_back(-1.040*ort);
  _A0mR2.push_back(5.28*5.28*GeV2);_A0mfit2.push_back(34.47*GeV2);
  _A1r1.push_back(0.217*ort);_A1r2.push_back(0.);
  _A1mR2.push_back(-1.0*GeV2);_A1mfit2.push_back(37.01*GeV2);
  _A2r1.push_back(0.006*ort);_A2r2.push_back(0.192*ort);
  _A2mR2.push_back(-1.0*GeV2);_A2mfit2.push_back(41.24*GeV2);
  _T1r1.push_back(0.865*ort);_T1r2.push_back(-0.622*ort);
  _T1mR2.push_back(5.32*5.32*GeV2);_T1mfit2.push_back(37.19*GeV2);
  _T2r1.push_back(0.242*ort);_T2r2.push_back(0.);
  _T2mR2.push_back(-1.0*GeV2);_T2mfit2.push_back(37.95*GeV2);
  _T3r1.push_back(0.023*ort);_T3r2.push_back(0.219*ort);
  _T3mR2.push_back(-1.0*GeV2);_T3mfit2.push_back(40.87*GeV2);
  addFormFactor(521,223,1,2,-5,2);
  // B_s to K*
  _Vr1.push_back(2.351);_Vr2.push_back(-2.039);
  _VmR2.push_back(5.42*5.42*GeV2);_Vmfit2.push_back(33.10*GeV2);
  _A0r1.push_back(2.813);_A0r2.push_back(-2.453);
  _A0mR2.push_back(5.37*5.37*GeV2);_A0mfit2.push_back(31.58*GeV2);
  _A1r1.push_back(0.231);_A1r2.push_back(0.);
  _A1mR2.push_back(-1.0*GeV2);_A1mfit2.push_back(32.94*GeV2);
  _A2r1.push_back(-0.011);_A2r2.push_back(0.192);
  _A2mR2.push_back(-1.0*GeV2);_A2mfit2.push_back(40.14*GeV2);
  _T1r1.push_back(2.047);_T1r2.push_back(-1.787);
  _T1mR2.push_back(5.42*5.42*GeV2);_T1mfit2.push_back(32.83*GeV2);
  _T2r1.push_back(0.260);_T2r2.push_back(0.);
  _T2mR2.push_back(-1.0*GeV2);_T2mfit2.push_back(33.01*GeV2);
  _T3r1.push_back(0.043);_T3r2.push_back(0.217);
  _T3mR2.push_back(-1.0*GeV2);_T3mfit2.push_back(39.38*GeV2);
  addFormFactor(531,-323,1,3,-5,-2); 
  // B_s to phi
  _Vr1.push_back(1.484);_Vr2.push_back(-1.049);
  _VmR2.push_back(5.42*5.42*GeV2);_Vmfit2.push_back(39.52*GeV2);
  _A0r1.push_back(3.310);_A0r2.push_back(-2.835);
  _A0mR2.push_back(5.37*5.37*GeV2);_A0mfit2.push_back(31.57*GeV2);
  _A1r1.push_back(0.308);_A1r2.push_back(0.);
  _A1mR2.push_back(-1.0*GeV2);_A1mfit2.push_back(36.54*GeV2);
  _A2r1.push_back(-0.054);_A2r2.push_back(0.288);
  _A2mR2.push_back(-1.0*GeV2);_A2mfit2.push_back(48.94*GeV2);
  _T1r1.push_back(1.303);_T1r2.push_back(-0.954);
  _T1mR2.push_back(5.42*5.42*GeV2);_T1mfit2.push_back(38.28*GeV2);
  _T2r1.push_back(0.349);_T2r2.push_back(0.);
  _T2mR2.push_back(-1.0*GeV2);_T2mfit2.push_back(37.21*GeV2);
  _T3r1.push_back(0.027);_T3r2.push_back(0.322);
  _T3mR2.push_back(-1.0*GeV2);_T3mfit2.push_back(45.56*GeV2);
  addFormFactor(531,333,1,3,-5,-3);
  initialModes(numberOfFactors());
  // cut-off parameter
  _cutoff=2.*GeV2;
}

BallZwickyVectorFormFactor::~BallZwickyVectorFormFactor() {}

void BallZwickyVectorFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _Vr1 << _Vr2 << _A0r1 << _A0r2 << _A1r1 << _A1r2 << _A2r1 << _A2r2 << _T1r1
     << _T1r2 << _T2r1 << _T2r2 << _T3r1 << _T3r2 << _VmR2 << _Vmfit2 << _A0mR2 
     << _A0mfit2 << _A1mR2 << _A1mfit2 << _A2mR2 << _A2mfit2 << _T1mR2 << _T1mfit2 
     << _T2mR2 << _T2mfit2 << _T3mR2 << _T3mfit2 << _cutoff;
}

void BallZwickyVectorFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _Vr1 >> _Vr2 >> _A0r1 >> _A0r2 >> _A1r1 >> _A1r2 >> _A2r1 >> _A2r2 >> _T1r1
     >> _T1r2 >> _T2r1 >> _T2r2 >> _T3r1 >> _T3r2 >> _VmR2 >> _Vmfit2 >> _A0mR2 
     >> _A0mfit2 >> _A1mR2 >> _A1mfit2 >> _A2mR2 >> _A2mfit2 >> _T1mR2 >> _T1mfit2 
     >> _T2mR2 >> _T2mfit2 >> _T3mR2 >> _T3mfit2 >> _cutoff;
}

ClassDescription<BallZwickyVectorFormFactor> BallZwickyVectorFormFactor::initBallZwickyVectorFormFactor;
// Definition of the static class description member.

void BallZwickyVectorFormFactor::Init() {

  static ClassDocumentation<BallZwickyVectorFormFactor> documentation
    ("The BallZwickyVectorFormFactor class implements the vector form"
     " factors of hep-ph/0412079 for the form-factor for the decay of a B-meson to a"
     " light pseudoscalar meson");

  static ParVector<BallZwickyVectorFormFactor,double> interfaceVr_1
    ("Vr_1",
     "The r_1 coefficient for the V form-factor",
     &BallZwickyVectorFormFactor::_Vr1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceVr_2
    ("Vr_2",
     "The r_2 coefficient for the V form-factor",
     &BallZwickyVectorFormFactor::_Vr2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA0r_1
    ("A0r_1",
     "The r_1 coefficient for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA0r_2
    ("A0r_2",
     "The r_2 coefficient for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA1r_1
    ("A1r_1",
     "The r_1 coefficient for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA1r_2
    ("A1r_2",
     "The r_2 coefficient for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA2r_1
    ("A2r_1",
     "The r_1 coefficient for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA2r_2
    ("A2r_2",
     "The r_2 coefficient for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT1r_1
    ("T1r_1",
     "The r_1 coefficient for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT1r_2
    ("T1r_2",
     "The r_2 coefficient for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT2r_1
    ("T2r_1",
     "The r_1 coefficient for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT2r_2
    ("T2r_2",
     "The r_2 coefficient for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT3r_1
    ("T3r_1",
     "The r_1 coefficient for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT3r_2
    ("T3r_2",
     "The r_2 coefficient for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceVmR2
    ("VmR2",
     "The value of m_R^2 for the V form-factor",
     &BallZwickyVectorFormFactor::_VmR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceVmfit2 
    ("Vmfit2",
     "The value of m_fit^2 for the V form-factor",
     &BallZwickyVectorFormFactor::_Vmfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA0mR2
    ("A0mR2",
     "The value of m_R^2 for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA0mfit2 
    ("A0mfit2",
     "The value of m_fit^2 for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA1mR2
    ("A1mR2",
     "The value of m_R^2 for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA1mfit2 
    ("A1mfit2",
     "The value of m_fit^2 for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA2mR2
    ("A2mR2",
     "The value of m_R^2 for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA2mfit2 
    ("A2mfit2",
     "The value of m_fit^2 for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT1mR2
    ("T1mR2",
     "The value of m_R^2 for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT1mfit2 
    ("T1mfit2",
     "The value of m_fit^2 for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT2mR2
    ("T2mR2",
     "The value of m_R^2 for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT2mfit2 
    ("T2mfit2",
     "The value of m_fit^2 for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT3mR2
    ("T3mR2",
     "The value of m_R^2 for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT3mfit2 
    ("T3mfit2",
     "The value of m_fit^2 for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static Parameter<BallZwickyVectorFormFactor,Energy2> interfaceCutOff
    ("CutOff",
     "Parameter controlling the value of q^2 where we switch from the fit "
     "to a small q^2 expansion for numerical stability.",
     &BallZwickyVectorFormFactor::_cutoff, GeV2, 2.0*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, true);

}

// form-factor for scalar to vector
void BallZwickyVectorFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
							int id0,int id1,Energy m0, 
							Energy m1,
							Complex & A0,Complex & A1,
							Complex & A2,Complex & V) const
{
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
    {A2 = (_A2r1[mode]+_A2r2[mode]/(1.-q2/_A2mfit2[mode]))/(1.-q2/_A2mfit2[mode]);}
  else
    {A2 = _A2r1[mode]/(1.-q2/_A2mR2[mode])+_A2r2[mode]/(1.-q2/_A2mfit2[mode]);}
  // V
  if(_VmR2[mode]<0)
    {V = (_Vr1[mode]+_Vr2[mode]/(1.-q2/_Vmfit2[mode]))/(1.-q2/_Vmfit2[mode]);}
  else
    {V = _Vr1[mode]/(1.-q2/_VmR2[mode])+_Vr2[mode]/(1.-q2/_Vmfit2[mode]);}
}

void BallZwickyVectorFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,
							     unsigned int mode,int id0,
							     int id1,Energy m0,Energy m1,
							     Complex & T1,Complex & T2,
							     Complex & T3) const
{
  // T_1
  if(_T1mR2[mode]<0)
    {T1 = (_T1r1[mode]+_T1r2[mode]/(1.-q2/_T1mfit2[mode]))/(1.-q2/_T1mfit2[mode]);}
  else
    {T1 = _T1r1[mode]/(1.-q2/_T1mR2[mode])+_T1r2[mode]/(1.-q2/_T1mfit2[mode]);}
  // T_2
  if(_T2mR2[mode]<0)
    {T2 = (_T2r1[mode]+_T2r2[mode]/(1.-q2/_T2mfit2[mode]))/(1.-q2/_T2mfit2[mode]);}
  else
    {T2 = _T2r1[mode]/(1.-q2/_T2mR2[mode])+_T2r2[mode]/(1.-q2/_T2mfit2[mode]);}
  // T_3
  if(q2>_cutoff)
    {
      if(_T3mR2[mode]<0)
	{T3 = (_T3r1[mode]+_T3r2[mode]/(1.-q2/_T3mfit2[mode]))/(1.-q2/_T3mfit2[mode]);}
      else
	{T3 = _T3r1[mode]/(1.-q2/_T3mR2[mode])+_T3r2[mode]/(1.-q2/_T3mfit2[mode]);}
      // convert for T_3tilde to T_3
      T3 = (m0*m0-m1*m1)/q2*(T3-T2);
    }
  else
    {
      InvEnergy2 smallT2,smallT3;
      if(_T2mR2[mode]<0)
	{
	  double a(q2/_T2mfit2[mode]);
	  smallT2=1./_T2mfit2[mode]*
	    (_T2r1[mode]+2.*_T2r2[mode]+a*(_T2r1[mode]+3.*_T2r2[mode]+
					   a*(_T2r1[mode]+4.*_T2r2[mode]+
					      a*(_T2r1[mode]+5.*_T2r2[mode]))));
	}
      else
	{smallT2=(_T2r1[mode]/_T2mR2[mode]+_T2r2[mode]/_T2mfit2[mode])
	    +q2*(+_T2r1[mode]/_T2mR2[mode]/_T2mR2[mode]
		 +_T2r2[mode]/_T2mfit2[mode]/_T2mfit2[mode])
	    +q2*q2*(+_T2r1[mode]/_T2mR2[mode]/_T2mR2[mode]/_T2mR2[mode]
		   +_T2r2[mode]/_T2mfit2[mode]/_T2mfit2[mode]/_T2mfit2[mode]);}
      if(_T3mR2[mode]<0)
	{
	  double a(q2/_T3mfit2[mode]);
	  smallT3=1./_T3mfit2[mode]*
	    (_T3r1[mode]+2.*_T3r2[mode]+a*(_T3r1[mode]+3.*_T3r2[mode]+
					   a*(_T3r1[mode]+4.*_T3r2[mode]+
					      a*(_T3r1[mode]+5.*_T3r2[mode]))));
	}
      else
	{smallT3=(_T3r1[mode]/_T3mR2[mode]+_T3r2[mode]/_T3mfit2[mode])
	    +q2*(+_T3r1[mode]/_T3mR2[mode]/_T3mR2[mode]
		 +_T3r2[mode]/_T3mfit2[mode]/_T3mfit2[mode])
	    +q2*q2*(+_T3r1[mode]/_T3mR2[mode]/_T3mR2[mode]/_T3mR2[mode]
		   +_T3r2[mode]/_T3mfit2[mode]/_T3mfit2[mode]/_T3mfit2[mode]);}
      T3 = (m0*m0-m1*m1)*(smallT3-smallT2);
    }
}

void BallZwickyVectorFormFactor::dataBaseOutput(ofstream & output,bool header,
						bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::BallZwickyVectorFormFactor " << fullName() << " \n";}
  output << "set " << fullName() << ":CutOff " << _cutoff/GeV2 << "\n";
  for(unsigned int ix=0;ix<_Vr1.size();++ix)
    {
      if(ix<initialModes())
	{
	  output << "set " << fullName() << ":Vr_1 "  << ix << "  " 
		 << _Vr1[ix]  << "\n";
	  output << "set " << fullName() << ":Vr_2 "  << ix << "  " 
		 << _Vr2[ix]  << "\n";
	  output << "set " << fullName() << ":A0r_1 "  << ix << "  " 
		 << _A0r1[ix] << "\n";
	  output << "set " << fullName() << ":A0r_2 "  << ix << "  " 
		 << _A0r2[ix] << "\n";
	  output << "set " << fullName() << ":A1r_1 "  << ix << "  " 
		 << _A1r1[ix] << "\n";
	  output << "set " << fullName() << ":A1r_2 "  << ix << "  " 
		 << _A1r2[ix] << "\n";
	  output << "set " << fullName() << ":A2r_1 "  << ix << "  " 
		 << _A2r1[ix] << "\n";
	  output << "set " << fullName() << ":A2r_2 "  << ix << "  " 
		 << _A2r2[ix] << "\n";
	  output << "set " << fullName() << ":T1r_1 "  << ix << "  " 
		 << _T1r1[ix] << "\n";
	  output << "set " << fullName() << ":T1r_2 "  << ix << "  " 
		 << _T1r2[ix] << "\n";
	  output << "set " << fullName() << ":T2r_1 "  << ix << "  " 
		 << _T2r1[ix] << "\n";
	  output << "set " << fullName() << ":T2r_2 "  << ix << "  " 
		 << _T2r2[ix] << "\n";
	  output << "set " << fullName() << ":T3r_1 "  << ix << "  " 
		 << _T3r1[ix] << "\n";
	  output << "set " << fullName() << ":T3r_2 "  << ix << "  " 
		 << _T3r2[ix] << "\n";
	  output << "set " << fullName() << ":VmR2 "  << ix 
		 << "  " << _VmR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":Vmfit2 "  << ix 
		 << "  " << _Vmfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A0mR2 "  << ix 
		 << "  " << _A0mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A0mfit2 "  << ix 
		 << "  " << _A0mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A1mR2 "  << ix 
		 << "  " << _A1mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A1mfit2 "  << ix 
		 << "  " << _A1mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A2mR2 "  << ix 
		 << "  " << _A2mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A2mfit2 "  << ix 
		 << "  " << _A2mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T1mR2 "  << ix 
		 << "  " << _T1mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T1mfit2 "  << ix 
		 << "  " << _T1mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T2mR2 "  << ix 
		 << "  " << _T2mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T2mfit2 "  << ix 
		 << "  " << _T2mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T3mR2 "  << ix 
		 << "  " << _T3mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T3mfit2 "  << ix 
		 << "  " << _T3mfit2[ix]/GeV2 << "\n";
	}
      else
	{
	  output << "set " << fullName() << ":Vr_1 "  << ix << "  " 
		 << _Vr1[ix] << "\n";
	  output << "set " << fullName() << ":Vr_2 "  << ix << "  " 
		 << _Vr2[ix] << "\n";
	  output << "set " << fullName() << ":A0r_1 "  << ix << "  " 
		 << _A0r1[ix] << "\n";
	  output << "set " << fullName() << ":A0r_2 "  << ix << "  " 
		 << _A0r2[ix] << "\n";
	  output << "set " << fullName() << ":A1r_1 "  << ix << "  " 
		 << _A1r1[ix] << "\n";
	  output << "set " << fullName() << ":A1r_2 "  << ix << "  " 
		 << _A1r2[ix] << "\n";
	  output << "set " << fullName() << ":A2r_1 "  << ix << "  " 
		 << _A2r1[ix] << "\n";
	  output << "set " << fullName() << ":A2r_2 "  << ix << "  " 
		 << _A2r2[ix] << "\n";
	  output << "set " << fullName() << ":T1r_1 "  << ix << "  " 
		 << _T1r1[ix] << "\n";
	  output << "set " << fullName() << ":T1r_2 "  << ix << "  " 
		 << _T1r2[ix] << "\n";
	  output << "set " << fullName() << ":T2r_1 "  << ix << "  " 
		 << _T2r1[ix] << "\n";
	  output << "set " << fullName() << ":T2r_2 "  << ix << "  " 
		 << _T2r2[ix] << "\n";
	  output << "set " << fullName() << ":T3r_1 "  << ix << "  " 
		 << _T3r1[ix] << "\n";
	  output << "set " << fullName() << ":T3r_2 "  << ix << "  " 
		 << _T3r2[ix] << "\n";
	  output << "set " << fullName() << ":VmR2 "  << ix 
		 << "  " << _VmR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":Vmfit2 "  << ix 
		 << "  " << _Vmfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A0mR2 "  << ix 
		 << "  " << _A0mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A0mfit2 "  << ix 
		 << "  " << _A0mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A1mR2 "  << ix 
		 << "  " << _A1mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A1mfit2 "  << ix 
		 << "  " << _A1mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A2mR2 "  << ix 
		 << "  " << _A2mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":A2mfit2 "  << ix 
		 << "  " << _A2mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T1mR2 "  << ix 
		 << "  " << _T1mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T1mfit2 "  << ix 
		 << "  " << _T1mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T2mR2 "  << ix 
		 << "  " << _T2mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T2mfit2 "  << ix 
		 << "  " << _T2mfit2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T3mR2 "  << ix 
		 << "  " << _T3mR2[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":T3mfit2 "  << ix 
		 << "  " << _T3mfit2[ix]/GeV2 << "\n";
	}
    }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
