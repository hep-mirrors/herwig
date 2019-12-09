// -*- C++ -*-
//
// BallZwickyVectorFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyVectorFormFactor class.
//

#include "BallZwickyVectorFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

BallZwickyVectorFormFactor::BallZwickyVectorFormFactor() 
  : _Vr1(9), _Vr2(9), _A0r1(9), _A0r2(9), _A1r1(9), _A1r2(9), 
    _A2r1(9), _A2r2(9), _T1r1(9), _T1r2(9), _T2r1(9), _T2r2(9), 
    _T3r1(9), _T3r2(9),  _VmR2(9), _Vmfit2(9), _A0mR2(9), _A0mfit2(9), 
    _A1mR2(9), _A1mfit2(9), _A2mR2(9), _A2mfit2(9), _T1mR2(9), _T1mfit2(9), 
    _T2mR2(9), _T2mfit2(9), _T3mR2(9), _T3mfit2(9) {
  double ort(1./sqrt(2.));
  // parameters for the different form-factors
  // B to rho
  addFormFactor(-521, 113,1,-2,5,2);
  addFormFactor(-511, 113,1,-2,5,1);
  addFormFactor(-511, 213,1,-2,5,2);  
  addFormFactor(-521, 213,1,-2,5,1);
  for(unsigned int ix=0;ix<4;++ix) {
    _Vr1[ix]   = 1.045         ; _Vr2[ix]     = -0.721; 
    _VmR2[ix]  = sqr(5.32)*GeV2; _Vmfit2[ix]  = 38.34*GeV2; 
    _A0r1[ix]  = 1.527         ; _A0r2[ix]    = -1.220; 
    _A0mR2[ix] = sqr(5.28)*GeV2; _A0mfit2[ix] = 33.36*GeV2; 
    _A1r1[ix]  = 0.240         ; _A1r2[ix]    = 0.; 
    _A1mR2[ix] = -1.0*GeV2     ; _A1mfit2[ix] = 37.51*GeV2; 
    _A2r1[ix]  = 0.009         ; _A2r2[ix]    = 0.212; 
    _A2mR2[ix] = -1.0*GeV2     ; _A2mfit2[ix] = 40.82*GeV2; 
    _T1r1[ix]  = 0.897         ; _T1r2[ix]    = -0.629; 
    _T1mR2[ix] = sqr(5.32)*GeV2; _T1mfit2[ix] = 38.04*GeV2; 
    _T2r1[ix]  = 0.267         ; _T2r2[ix]    = 0.; 
    _T2mR2[ix] = -1.0*GeV2     ; _T2mfit2[ix] = 38.59*GeV2; 
    _T3r1[ix]  = 0.022         ; _T3r2[ix]    = 0.245; 
    _T3mR2[ix] = -1.0*GeV2     ; _T3mfit2[ix] = 40.88*GeV2;
  }
  for(unsigned int ix=0;ix<2;++ix) {
    double fact = ix==0 ? ort : -ort;
    _Vr1[ix]  *= fact; _Vr2[ix]  *= fact;
    _A0r1[ix] *= fact; _A0r2[ix] *= fact;
    _A1r1[ix] *= fact; _A1r2[ix] *= fact;
    _A2r1[ix] *= fact; _A2r2[ix] *= fact;
    _T1r1[ix] *= fact; _T1r2[ix] *= fact;
    _T2r1[ix] *= fact; _T2r2[ix] *= fact;
    _T3r1[ix] *= fact; _T3r2[ix] *= fact;
  }
  // parameters for the B to K   form-factors
  addFormFactor(-521,-323,1,-2,5,3);  
  addFormFactor(-511,-313,1,-2,5,3); 
  for(unsigned int ix=4;ix<6;++ix) {
    _Vr1[ix]   = 0.923         ; _Vr2[ix] = -0.511; 
    _VmR2[ix]  = sqr(5.32)*GeV2; _Vmfit2[ix] = 49.40*GeV2; 
    _A0r1[ix]  = 1.364         ; _A0r2[ix] = -0.990; 
    _A0mR2[ix] = sqr(5.28)*GeV2; _A0mfit2[ix] = 36.78*GeV2; 
    _A1r1[ix]  = 0.290         ; _A1r2[ix] = 0.0; 
    _A1mR2[ix] = -1.0*GeV2     ; _A1mfit2[ix] = 40.38*GeV2; 
    _A2r1[ix]  = -0.084        ; _A2r2[ix] = 0.342; 
    _A2mR2[ix] = -1.0*GeV2     ; _A2mfit2[ix] = 52.00*GeV2; 
    _T1r1[ix]  = 0.823         ; _T1r2[ix] = -0.491; 
    _T1mR2[ix] = sqr(5.32)*GeV2; _T1mfit2[ix] = 46.31*GeV2; 
    _T2r1[ix]  = 0.333         ; _T2r2[ix] = 0.; 
    _T2mR2[ix] = -1.0*GeV2     ; _T2mfit2[ix] = 41.41*GeV2; 
    _T3r1[ix]  = -0.036        ; _T3r2[ix] = 0.369; 
    _T3mR2[ix] = -1.0*GeV2     ; _T3mfit2[ix] = 48.10*GeV2;
  }
  // B to omega
  addFormFactor(-521,223,1,-2,5,2);
  _Vr1[6]   = 1.006*ort     ; _Vr2[6]     = -0.713*ort; 
  _VmR2[6]  = 5.32*5.32*GeV2; _Vmfit2[6]  = 37.45*GeV2; 
  _A0r1[6]  = 1.321*ort     ; _A0r2[6]    = -1.040*ort; 
  _A0mR2[6] = 5.28*5.28*GeV2; _A0mfit2[6] = 34.47*GeV2; 
  _A1r1[6]  = 0.217*ort     ; _A1r2[6]    = 0.; 
  _A1mR2[6] = -1.0*GeV2     ; _A1mfit2[6] = 37.01*GeV2; 
  _A2r1[6]  = 0.006*ort     ; _A2r2[6]    = 0.192*ort; 
  _A2mR2[6] = -1.0*GeV2     ; _A2mfit2[6] = 41.24*GeV2; 
  _T1r1[6]  = 0.865*ort     ; _T1r2[6]    = -0.622*ort; 
  _T1mR2[6] = 5.32*5.32*GeV2; _T1mfit2[6] = 37.19*GeV2; 
  _T2r1[6]  = 0.242*ort     ; _T2r2[6]    = 0.; 
  _T2mR2[6] = -1.0*GeV2     ; _T2mfit2[6] = 37.95*GeV2; 
  _T3r1[6]  = 0.023*ort     ; _T3r2[6]    = 0.219*ort; 
  _T3mR2[6] = -1.0*GeV2     ; _T3mfit2[6] = 40.87*GeV2; 
  // B_s to K*
  addFormFactor(-531,323,1,-3,5,2); 
  _Vr1[7]   = 2.351         ; _Vr2[7]     = -2.039; 
  _VmR2[7]  = sqr(5.42)*GeV2; _Vmfit2[7]  = 33.10*GeV2; 
  _A0r1[7]  = 2.813         ; _A0r2[7]    = -2.450; 
  _A0mR2[7] = sqr(5.37)*GeV2; _A0mfit2[7] = 31.58*GeV2; 
  _A1r1[7]  = 0.231         ; _A1r2[7]    = 0.; 
  _A1mR2[7] = -1.0*GeV2     ; _A1mfit2[7] = 32.94*GeV2; 
  _A2r1[7]  = -0.011        ; _A2r2[7]    = 0.192; 
  _A2mR2[7] = -1.0*GeV2     ; _A2mfit2[7] = 40.14*GeV2; 
  _T1r1[7]  = 2.047         ; _T1r2[7]    = -1.787; 
  _T1mR2[7] = sqr(5.42)*GeV2; _T1mfit2[7] = 32.83*GeV2; 
  _T2r1[7]  = 0.260         ; _T2r2[7]    = 0.; 
  _T2mR2[7] = -1.0*GeV2     ; _T2mfit2[7] = 33.01*GeV2; 
  _T3r1[7]  = 0.043         ; _T3r2[7]    = 0.217; 
  _T3mR2[7] = -1.0*GeV2     ; _T3mfit2[7] = 39.38*GeV2; 
  // B_s to phi
  addFormFactor(-531,333,1,-3,5,3);
  _Vr1[8]   = 1.484         ; _Vr2[8]     = -1.049; 
  _VmR2[8]  = sqr(5.42)*GeV2; _Vmfit2[8]  = 39.52*GeV2; 
  _A0r1[8]  = 3.310         ; _A0r2[8]    = -2.835; 
  _A0mR2[8] = sqr(5.37)*GeV2; _A0mfit2[8] = 31.57*GeV2; 
  _A1r1[8]  = 0.308         ; _A1r2[8]    = 0.; 
  _A1mR2[8] = -1.0*GeV2     ; _A1mfit2[8] = 36.54*GeV2; 
  _A2r1[8]  = -0.054        ; _A2r2[8]    = 0.288; 
  _A2mR2[8] = -1.0*GeV2     ; _A2mfit2[8] = 48.94*GeV2; 
  _T1r1[8]  = 1.303         ; _T1r2[8]    = -0.954; 
  _T1mR2[8] = sqr(5.42)*GeV2; _T1mfit2[8] = 38.28*GeV2; 
  _T2r1[8]  = 0.349         ; _T2r2[8]    = 0.; 
  _T2mR2[8] = -1.0*GeV2     ; _T2mfit2[8] = 37.21*GeV2; 
  _T3r1[8]  = 0.027         ; _T3r2[8]    = 0.322; 
  _T3mR2[8] = -1.0*GeV2     ; _T3mfit2[8] = 45.56*GeV2; 
  initialModes(numberOfFactors());
  // cut-off parameter
  _cutoff=0.01*GeV2;
}

void BallZwickyVectorFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_Vr1.size()||isize!=_Vr2.size()||isize!=_A0r1.size()||isize!=_A0r2.size()||
     isize!=_A1r1.size()||isize!=_A1r2.size()||isize!=_A2r1.size()||isize!=_A2r2.size()||
     isize!=_T1r1.size()||isize!=_T1r2.size()||isize!=_T2r1.size()||isize!=_T2r2.size()||
     isize!=_T3r1.size()||isize!=_T3r2.size()||isize!=_VmR2.size()||
     isize!=_Vmfit2.size()||isize!=_A0mR2.size()||isize!=_A0mfit2.size()||
     isize!=_A1mR2.size()||isize!=_A1mfit2.size()||isize!=_A2mR2.size()||
     isize!=_A2mfit2.size()||isize!=_T1mR2.size()||isize!=_T1mfit2.size()||
     isize!=_T2mR2.size()||isize!=_T2mfit2.size()||isize!=_T3mR2.size()||
     isize!=_T3mfit2.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "BallZwickyScalarFormFactor::doinit()" 
			  << Exception::abortnow;
  // output some graphs to check the answers
//   int id0,id1;
//   unsigned int iz;
//   Energy m0,m1;
//   Energy2 q2,step(14./100.*GeV2);
//   tcPDPtr in,out;
//   Complex A0,A1,A2,V;
//   ofstream output("Ball.top");
//   for(unsigned int ix=0;ix<numberOfFactors();++ix) {
//     particleID(ix,id0,id1);
//     in = getParticleData(id0);
//     m0=in->mass();
//     out= getParticleData(id1);
//     m1=out->mass();
//     output << "new frame " << endl;
//     output << "newdef font duplex" << endl;
//     output << "title top \"" << in->PDGName() << " to " << out->PDGName() 
// 	   << " vector form factors \"" << endl;
//     output << "newdef limits x 0 14. y 0 1" << endl;
//     double rt(sqrt(2.));
//     for(iz=0;iz<4;++iz) {
//       q2=ZERO;
//       for( ;q2<14.*GeV2+step;q2+=step) {
// 	ScalarVectorFormFactor(q2,ix,id0,id1,m0,m1,A0,A1,A2,V);
// 	if(id1==113||id1==223) {
//  	  if((abs(id0)%100)/10==1) {
// 	    A0*=-rt;
// 	    A1*=-rt;
// 	    A2*=-rt;
// 	    V *=-rt;
// 	  }
// 	  else {
// 	    A0*=rt;
// 	    A1*=rt;
// 	    A2*=rt;
// 	    V*=rt;
// 	  }
// 	}
// 	if(iz==0)      output << q2/GeV2 << "   " << A0.real() << endl;
// 	else if(iz==1) output << q2/GeV2 << "   " << A1.real() << endl;
// 	else if(iz==2) output << q2/GeV2 << "   " << A2.real() << endl;
// 	else if(iz==3) output << q2/GeV2 << "   " << V.real()  << endl;
//       }
//       if(iz==0)      output << "join red"    << endl;
//       else if(iz==1) output << "join blue"   << endl;
//       else if(iz==2) output << "join green"  << endl;
//       else if(iz==3) output << "join yellow" << endl;
//     }
//     output << "new frame " << endl;
//     output << "newdef font duplex" << endl;
//     output << "title top \"" << in->PDGName() << " to " << out->PDGName() 
// 	   << " penguin form factors\" " << endl;
//     output << "newdef limits x 0 14. y 0 1" << endl;
//     for(iz=0;iz<3;++iz) {
//       q2=ZERO;
//       for( ;q2<14.*GeV2+step;q2+=step) {
// 	ScalarVectorSigmaFormFactor(q2,ix,id0,id1,m0,m1,A0,A1,A2);
// 	if(id1==113||id1==223) {
//  	  if((abs(id0)%100)/10==1) {
// 	    A0*=-rt;
// 	    A1*=-rt;
// 	    A2*=-rt;
// 	  }
// 	  else {
// 	    A0*=rt;
// 	    A1*=rt;
// 	    A2*=rt;
// 	  }
// 	}
// 	if(iz==0)      output << q2/GeV2 << "   " << A0.real() << endl;
// 	else if(iz==1) output << q2/GeV2 << "   " << A1.real() << endl;
// 	else if(iz==2) output << q2/GeV2 << "   " << A2.real() << endl;
//       }
//       if(iz==0){output      << "join red"    << endl;}
//       else if(iz==1){output << "join blue"   << endl;}
//       else if(iz==2){output << "join green"  << endl;}
//     }
//   }
}

void BallZwickyVectorFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _Vr1 << _Vr2 << _A0r1 << _A0r2 << _A1r1 << _A1r2 << _A2r1 << _A2r2 << _T1r1
     << _T1r2 << _T2r1 << _T2r2 << _T3r1 << _T3r2 
     << ounit(_VmR2,GeV2) << ounit(_Vmfit2,GeV2) << ounit(_A0mR2,GeV2) 
     << ounit(_A0mfit2,GeV2) << ounit(_A1mR2,GeV2) << ounit(_A1mfit2,GeV2) 
     << ounit(_A2mR2,GeV2) << ounit(_A2mfit2,GeV2) << ounit(_T1mR2,GeV2) 
     << ounit(_T1mfit2,GeV2) << ounit(_T2mR2,GeV2) << ounit(_T2mfit2,GeV2) 
     << ounit(_T3mR2,GeV2) << ounit(_T3mfit2,GeV2) << ounit(_cutoff,GeV2);
}

void BallZwickyVectorFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _Vr1 >> _Vr2 >> _A0r1 >> _A0r2 >> _A1r1 >> _A1r2 >> _A2r1 >> _A2r2 >> _T1r1
     >> _T1r2 >> _T2r1 >> _T2r2 >> _T3r1 >> _T3r2 
     >> iunit(_VmR2,GeV2) >> iunit(_Vmfit2,GeV2) >> iunit(_A0mR2 ,GeV2)
     >> iunit(_A0mfit2,GeV2) >> iunit(_A1mR2,GeV2) >> iunit(_A1mfit2,GeV2) 
     >> iunit(_A2mR2,GeV2) >> iunit(_A2mfit2,GeV2) >> iunit(_T1mR2,GeV2) 
     >> iunit(_T1mfit2,GeV2) >> iunit(_T2mR2,GeV2) >> iunit(_T2mfit2,GeV2) 
     >> iunit(_T3mR2,GeV2) >> iunit(_T3mfit2,GeV2) >> iunit(_cutoff,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BallZwickyVectorFormFactor,ScalarFormFactor>
describeHerwigBallZwickyVectorFormFactor("Herwig::BallZwickyVectorFormFactor", "HwFormFactors.so");

void BallZwickyVectorFormFactor::Init() {

  static ClassDocumentation<BallZwickyVectorFormFactor> documentation
    ("The BallZwickyVectorFormFactor class implements the vector form"
     " factors of hep-ph/0412079 for the form-factor for the decay of a B-meson to a"
     " light pseudoscalar meson",
     "The form factors of \\cite{Ball:2004rg} for $B_{d,s}\\to\\rho,\\omega,K^*,\\phi$"
     " were used.",
     "\\bibitem{Ball:2004rg} P.~Ball and R.~Zwicky, \n"
     "Phys.\\ Rev.\\  D {\\bf 71} (2005) 014029 [arXiv:hep-ph/0412079].\n"
     "%%CITATION = PHRVA,D71,014029;%%\n");

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
     &BallZwickyVectorFormFactor::_cutoff, GeV2, 2.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);

}

// form-factor for scalar to vector
void BallZwickyVectorFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
							int,int,Energy, 
							Energy,
							Complex & A0,Complex & A1,
							Complex & A2,Complex & V) const {
  useMe();
  // the form-factors
  // A_0
  if(_A0mR2[mode]<ZERO) {
    A0 = (_A0r1[mode]+_A0r2[mode]/(1.-q2/_A0mfit2[mode]))/(1.-q2/_A0mfit2[mode]);
  }
  else {
    A0 = _A0r1[mode]/(1.-q2/_A0mR2[mode])+_A0r2[mode]/(1.-q2/_A0mfit2[mode]);
  }
  // A_1
  if(_A1mR2[mode]<ZERO) {
    A1 = (_A1r1[mode]+_A1r2[mode]/(1.-q2/_A1mfit2[mode]))/(1.-q2/_A1mfit2[mode]);
  }
  else {
    A1 = _A1r1[mode]/(1.-q2/_A1mR2[mode])+_A1r2[mode]/(1.-q2/_A1mfit2[mode]);
  }
  // A_2
  if(_A2mR2[mode]<ZERO) {
    A2 = (_A2r1[mode]+_A2r2[mode]/(1.-q2/_A2mfit2[mode]))/(1.-q2/_A2mfit2[mode]);
  }
  else {
    A2 = _A2r1[mode]/(1.-q2/_A2mR2[mode])+_A2r2[mode]/(1.-q2/_A2mfit2[mode]);
  }
  // V
  if(_VmR2[mode]<ZERO) {
    V = (_Vr1[mode]+_Vr2[mode]/(1.-q2/_Vmfit2[mode]))/(1.-q2/_Vmfit2[mode]);
  }
  else {
    V = _Vr1[mode]/(1.-q2/_VmR2[mode])+_Vr2[mode]/(1.-q2/_Vmfit2[mode]);
  }
}

void BallZwickyVectorFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,
							     unsigned int mode,int,
							     int,Energy m0,Energy m1,
							     Complex & T1,Complex & T2,
							     Complex & T3) const {
  useMe();
  // T_1
  if(_T1mR2[mode]<ZERO) {
    T1 = (_T1r1[mode]+_T1r2[mode]/(1.-q2/_T1mfit2[mode]))/(1.-q2/_T1mfit2[mode]);
  }
  else {
    T1 = _T1r1[mode]/(1.-q2/_T1mR2[mode])+_T1r2[mode]/(1.-q2/_T1mfit2[mode]);
  }
  // T_2
  if(_T2mR2[mode]<ZERO) {
    T2 = (_T2r1[mode]+_T2r2[mode]/(1.-q2/_T2mfit2[mode]))/(1.-q2/_T2mfit2[mode]);
  }
  else {
    T2 = _T2r1[mode]/(1.-q2/_T2mR2[mode])+_T2r2[mode]/(1.-q2/_T2mfit2[mode]);
  }
  // T_3
  if(q2>_cutoff) {
    if(_T3mR2[mode]<ZERO) {
      T3 = (_T3r1[mode]+_T3r2[mode]/(1.-q2/_T3mfit2[mode]))/(1.-q2/_T3mfit2[mode]);
    }
    else {
      T3 = _T3r1[mode]/(1.-q2/_T3mR2[mode])+_T3r2[mode]/(1.-q2/_T3mfit2[mode]);
    }
    // convert for T_3tilde to T_3
    T3 = Complex((m0*m0-m1*m1)/q2*(T3-T2));
  }
  else {
    InvEnergy2 smallT2,smallT3;
    if(_T2mR2[mode]<ZERO) {
      double a(q2/_T2mfit2[mode]);
      smallT2=1./_T2mfit2[mode]*
	(_T2r1[mode]+2.*_T2r2[mode]+a*(_T2r1[mode]+3.*_T2r2[mode]+
				       a*(_T2r1[mode]+4.*_T2r2[mode]+
					  a*(_T2r1[mode]+5.*_T2r2[mode]))));
    }
    else {
      smallT2=(_T2r1[mode]/_T2mR2[mode]+_T2r2[mode]/_T2mfit2[mode])
	+q2*(+_T2r1[mode]/_T2mR2[mode]/_T2mR2[mode]
	     +_T2r2[mode]/_T2mfit2[mode]/_T2mfit2[mode])
	+q2*q2*(+_T2r1[mode]/_T2mR2[mode]/_T2mR2[mode]/_T2mR2[mode]
		+_T2r2[mode]/_T2mfit2[mode]/_T2mfit2[mode]/_T2mfit2[mode]);
    }
    if(_T3mR2[mode]<ZERO) {
      double a(q2/_T3mfit2[mode]);
      smallT3=1./_T3mfit2[mode]*
	(_T3r1[mode]+2.*_T3r2[mode]+a*(_T3r1[mode]+3.*_T3r2[mode]+
				       a*(_T3r1[mode]+4.*_T3r2[mode]+
					  a*(_T3r1[mode]+5.*_T3r2[mode]))));
    }
    else {
      smallT3=(_T3r1[mode]/_T3mR2[mode]+_T3r2[mode]/_T3mfit2[mode])
	+q2*(+_T3r1[mode]/_T3mR2[mode]/_T3mR2[mode]
	     +_T3r2[mode]/_T3mfit2[mode]/_T3mfit2[mode])
	+q2*q2*(+_T3r1[mode]/_T3mR2[mode]/_T3mR2[mode]/_T3mR2[mode]
		+_T3r2[mode]/_T3mfit2[mode]/_T3mfit2[mode]/_T3mfit2[mode]);
    }
    T3 = (m0+m1)*(m0-m1)*(smallT3-smallT2);
  }
}

void BallZwickyVectorFormFactor::dataBaseOutput(ofstream & output,bool header,
						bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BallZwickyVectorFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":CutOff " << _cutoff/GeV2 << "\n";
  for(unsigned int ix=0;ix<_Vr1.size();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":Vr_1 "  << ix << "  " 
	     << _Vr1[ix]  << "\n";
      output << "newdef " << name() << ":Vr_2 "  << ix << "  " 
	     << _Vr2[ix]  << "\n";
      output << "newdef " << name() << ":A0r_1 "  << ix << "  " 
	     << _A0r1[ix] << "\n";
      output << "newdef " << name() << ":A0r_2 "  << ix << "  " 
	     << _A0r2[ix] << "\n";
      output << "newdef " << name() << ":A1r_1 "  << ix << "  " 
	     << _A1r1[ix] << "\n";
      output << "newdef " << name() << ":A1r_2 "  << ix << "  " 
	     << _A1r2[ix] << "\n";
      output << "newdef " << name() << ":A2r_1 "  << ix << "  " 
	     << _A2r1[ix] << "\n";
      output << "newdef " << name() << ":A2r_2 "  << ix << "  " 
	     << _A2r2[ix] << "\n";
      output << "newdef " << name() << ":T1r_1 "  << ix << "  " 
	     << _T1r1[ix] << "\n";
      output << "newdef " << name() << ":T1r_2 "  << ix << "  " 
	     << _T1r2[ix] << "\n";
      output << "newdef " << name() << ":T2r_1 "  << ix << "  " 
	     << _T2r1[ix] << "\n";
      output << "newdef " << name() << ":T2r_2 "  << ix << "  " 
	     << _T2r2[ix] << "\n";
      output << "newdef " << name() << ":T3r_1 "  << ix << "  " 
	     << _T3r1[ix] << "\n";
      output << "newdef " << name() << ":T3r_2 "  << ix << "  " 
	     << _T3r2[ix] << "\n";
      output << "newdef " << name() << ":VmR2 "  << ix 
	     << "  " << _VmR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":Vmfit2 "  << ix 
	     << "  " << _Vmfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mR2 "  << ix 
	     << "  " << _A0mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mfit2 "  << ix 
	     << "  " << _A0mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mR2 "  << ix 
	     << "  " << _A1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mfit2 "  << ix 
	     << "  " << _A1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mR2 "  << ix 
	     << "  " << _A2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mfit2 "  << ix 
	     << "  " << _A2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mR2 "  << ix 
	     << "  " << _T1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mfit2 "  << ix 
	     << "  " << _T1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mR2 "  << ix 
	     << "  " << _T2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mfit2 "  << ix 
	     << "  " << _T2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mR2 "  << ix 
	     << "  " << _T3mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mfit2 "  << ix 
	     << "  " << _T3mfit2[ix]/GeV2 << "\n";
    }
    else {
      output << "newdef " << name() << ":Vr_1 "  << ix << "  " 
	     << _Vr1[ix] << "\n";
      output << "newdef " << name() << ":Vr_2 "  << ix << "  " 
	     << _Vr2[ix] << "\n";
      output << "newdef " << name() << ":A0r_1 "  << ix << "  " 
	     << _A0r1[ix] << "\n";
      output << "newdef " << name() << ":A0r_2 "  << ix << "  " 
	     << _A0r2[ix] << "\n";
      output << "newdef " << name() << ":A1r_1 "  << ix << "  " 
	     << _A1r1[ix] << "\n";
      output << "newdef " << name() << ":A1r_2 "  << ix << "  " 
	     << _A1r2[ix] << "\n";
      output << "newdef " << name() << ":A2r_1 "  << ix << "  " 
	     << _A2r1[ix] << "\n";
      output << "newdef " << name() << ":A2r_2 "  << ix << "  " 
	     << _A2r2[ix] << "\n";
      output << "newdef " << name() << ":T1r_1 "  << ix << "  " 
	     << _T1r1[ix] << "\n";
      output << "newdef " << name() << ":T1r_2 "  << ix << "  " 
	     << _T1r2[ix] << "\n";
      output << "newdef " << name() << ":T2r_1 "  << ix << "  " 
	     << _T2r1[ix] << "\n";
      output << "newdef " << name() << ":T2r_2 "  << ix << "  " 
	     << _T2r2[ix] << "\n";
      output << "newdef " << name() << ":T3r_1 "  << ix << "  " 
	     << _T3r1[ix] << "\n";
      output << "newdef " << name() << ":T3r_2 "  << ix << "  " 
	     << _T3r2[ix] << "\n";
      output << "newdef " << name() << ":VmR2 "  << ix 
	     << "  " << _VmR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":Vmfit2 "  << ix 
	     << "  " << _Vmfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mR2 "  << ix 
	     << "  " << _A0mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mfit2 "  << ix 
	     << "  " << _A0mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mR2 "  << ix 
	     << "  " << _A1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mfit2 "  << ix 
	     << "  " << _A1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mR2 "  << ix 
	     << "  " << _A2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mfit2 "  << ix 
	     << "  " << _A2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mR2 "  << ix 
	     << "  " << _T1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mfit2 "  << ix 
	     << "  " << _T1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mR2 "  << ix 
	     << "  " << _T2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mfit2 "  << ix 
	     << "  " << _T2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mR2 "  << ix 
	     << "  " << _T3mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mfit2 "  << ix 
	     << "  " << _T3mfit2[ix]/GeV2 << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
