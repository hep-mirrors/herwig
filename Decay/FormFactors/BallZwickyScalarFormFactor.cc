// -*- C++ -*-
//
// BallZwickyScalarFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyScalarFormFactor class.
//

#include "BallZwickyScalarFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

BallZwickyScalarFormFactor::BallZwickyScalarFormFactor()  
  : _r10(7), _r20(7), _r1plus(7), _r2plus(7), _r1T(7), _r2T(7), _m120(7),
    _mfit20(7), _m12plus(7), _mfit2plus(7), _m12T(7), _mfit2T(7) {
  // parameters for the B to pi  form-factors
  for(unsigned int ix=0;ix<4;++ix) {
    _r10[ix]     = 0.            ; _r20[ix]       = 0.258; 
    _m120[ix]    = -1.*GeV2      ; _mfit20[ix]    = 33.81*GeV2; 
    _r1plus[ix]  = 0.744         ; _r2plus[ix]    = -0.486; 
    _m12plus[ix] = sqr(5.32)*GeV2; _mfit2plus[ix] = 40.73*GeV2; 
    _r1T[ix]     = 1.387         ; _r2T[ix]       = -1.134; 
    _m12T[ix]    = sqr(5.32)*GeV2; _mfit2T[ix]    = 32.22*GeV2; 
  }
  addFormFactor(-521, 111,0,-2,5,2);
  addFormFactor(-511, 111,0,-2,5,1);
  addFormFactor(-511, 211,0,-2,5,2);  
  addFormFactor(-521, 211,0,-2,5,1);
  for(unsigned int ix=0;ix<2;++ix) {
    double fact(sqrt(0.5));
    if(ix==1) fact *= -1.;
    _r20[ix] *= fact; _r1plus[ix] *= fact; _r2plus[ix] *= fact; 
    _r1T[ix] *= fact; _r2T[ix]    *= fact; 
  }
  // parameters for the B to K   form-factors
  addFormFactor(-521,-321,0,-2,5,3);
  addFormFactor(-511,-311,0,-2,5,3);
  for(unsigned int ix=4;ix<6;++ix) {
    _r10[ix]     = 0.            ; _r20[ix]       = 0.330; 
    _m120[ix]    = -1.*GeV2      ; _mfit20[ix]    = 37.46*GeV2; 
    _r1plus[ix]  = 0.162         ; _r2plus[ix]    = 0.173; 
    _m12plus[ix] = sqr(5.41)*GeV2; _mfit2plus[ix] = -1.*GeV2; 
    _r1T[ix]     = 0.161         ; _r2T[ix]       = 0.198; 
    _m12T[ix]    = sqr(5.41)*GeV2; _mfit2T[ix]    = -1.*GeV2; 
  }
  // parameters for the B to eta form-factors
  addFormFactor(521,221,0,2,-5,-2); 
  _r10[6]     = 0.            ; _r20[6]       = 0.273; 
  _m120[6]    = -1.*GeV2      ; _mfit20[6]    = 31.03*GeV2; 
  _r1plus[6]  = 0.122         ; _r2plus[6]    = 0.155; 
  _m12plus[6] = sqr(5.32)*GeV2; _mfit2plus[6] = -1.*GeV2; 
  _r1T[6]     = 0.111         ; _r2T[6]       = 0.175; 
  _m12T[6]    = sqr(5.32)*GeV2; _mfit2T[6]    = -1.*GeV2; 
  // initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta = -Constants::pi/9.;
}

void BallZwickyScalarFormFactor::doinit() {
  ScalarFormFactor::doinit();
  // check all the vectors have the same size
  unsigned int isize=numberOfFactors();
  if(isize!=_r10.size()||isize!=_r20.size()||isize!=_r1plus.size()||
     isize!=_r2plus.size()||isize!=_r1T.size()||
     isize!=_r2T.size()||isize!=_m120.size()||isize!=_mfit20.size()||
     isize!=_m12plus.size()||isize!=_mfit2plus.size()||
     isize!=_m12T.size()||isize!=_mfit2T.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "BallZwickyScalarFormFactor::doinit()" 
			  << Exception::abortnow;
  // output some graphs to check the answers
//   int id0,id1;
//   unsigned int iz;
//   Energy m0,m1; 
//   Energy2 q2,step(14./100.*GeV2);
//   tcPDPtr in,out;
//   Complex f0,fp,ft;
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
// 	   << " scalar form factors \"" << endl;
//     output << "newdef limits x 0 14. y 0 1" << endl;
//     double rt(sqrt(2.));
//     for(iz=0;iz<3;++iz) {
//       q2=ZERO;
//       for( ;q2<14.*GeV2+step;q2+=step) {
// 	ScalarScalarFormFactor(q2,ix,id0,id1,m0,m1,f0,fp);
// 	ScalarScalarSigmaFormFactor(q2,ix,id0,id1,m0,m1,ft);
// 	if(id1==111) {
// 	  if((abs(id0)%100)/10==1) {
// 	    f0*=-rt;
// 	    fp*=-rt;
// 	    ft*=-rt;
// 	  }
// 	  else {
// 	    f0*=rt;
// 	    fp*=rt;
// 	    ft*=rt;
// 	  }
// 	}
// 	else if(id1==221) {
// 	  double fact(cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.));
// 	  f0/=fact;
// 	  fp/=fact;
// 	  ft/=fact;
// 	}
// 	if(iz==0)      output << q2/GeV2 << "   " << f0.real() << "\n";
// 	else if(iz==1) output << q2/GeV2 << "   " << fp.real() << "\n";
// 	else if(iz==2) output << q2/GeV2 << "   " << ft.real() << "\n";
//       }
//       if(iz==0)      output << "join red  \n";
//       else if(iz==1) output << "join blue \n";
//       else if(iz==2) output << "join green\n";
//     }
//   }
}

void BallZwickyScalarFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _r10 << _r20 << _r1plus << _r2plus << _r1T << _r2T 
     << ounit(_m120,GeV2) << ounit(_mfit20,GeV2) 
     << ounit(_m12plus,GeV2) << ounit(_mfit2plus,GeV2) 
     << ounit(_m12T,GeV2) << ounit(_mfit2T,GeV2) << _thetaeta;
}
  
void BallZwickyScalarFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _r10 >> _r20 >> _r1plus >> _r2plus >> _r1T >> _r2T 
     >> iunit(_m120,GeV2) >> iunit(_mfit20,GeV2) 
     >> iunit(_m12plus,GeV2) >> iunit(_mfit2plus,GeV2) 
     >> iunit(_m12T,GeV2) >> iunit(_mfit2T,GeV2) >> _thetaeta;
}

ClassDescription<BallZwickyScalarFormFactor> BallZwickyScalarFormFactor::initBallZwickyScalarFormFactor;
// Definition of the static class description member.

void BallZwickyScalarFormFactor::Init() {

  static ClassDocumentation<BallZwickyScalarFormFactor> documentation
    ("The BallZwickyScalarFormFactor class implements the form-factors"
     " of PRD71 014015 (2005) for the form-factor for the decay of a B-meson to a"
     " light pseudoscalar meson",
     "The form factors of \\cite{Ball:2004ye} for $B\\to\\pi, K, \\eta$ were used.",
     "\\bibitem{Ball:2004ye} P.~Ball and R.~Zwicky,\n "
     "Phys.\\ Rev.\\  D {\\bf 71} (2005) 014015 [arXiv:hep-ph/0406232].\n"
     "%%CITATION = PHRVA,D71,014015;%%\n");

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_10
    ("r_10",
     "The r_1 coefficient for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_r10,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_20
    ("r_20",
     "The r_2 coefficient for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_r20,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_1plus
    ("r_1plus",
     "The r_1 coefficient for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_r1plus,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_2plus
    ("r_2plus",
     "The r_2 coefficient for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_r2plus,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_1T
    ("r_1T",
     "The r_1 coefficient for the f_T form-factor",
     &BallZwickyScalarFormFactor::_r1T,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_2T
    ("r_2T",
     "The r_2 coefficient for the f_T form-factor",
     &BallZwickyScalarFormFactor::_r2T,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem120
    ("m_120",
     "The value of m_1^2 for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_m120,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit20
    ("mfit20",
     "The value of m_fit^2 for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_mfit20,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem12plus
    ("m_12plus",
     "The value of m_1^2 for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_m12plus,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit2plus
    ("mfit2plus",
     "The value of m_fit^2 for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_mfit2plus,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem12T
    ("m_12T",
     "The value of m_1^2 for the f_T form-factor",
     &BallZwickyScalarFormFactor::_m12T,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit2T
    ("mfit2T",
     "The value of m_fit^2 for the f_T form-factor",
     &BallZwickyScalarFormFactor::_mfit2T,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static Parameter<BallZwickyScalarFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &BallZwickyScalarFormFactor::_thetaeta, -Constants::pi/9.,
     -Constants::pi, Constants::pi,
     false, false, true);
}

// form-factor for scalar to scalar
void BallZwickyScalarFormFactor::
ScalarScalarFormFactor(Energy2 q2,unsigned  int mode,
		       int, int id1, Energy, Energy,
		       Complex & f0, Complex & fp) const {
  useMe();
  // the F_0 form-factor
  if(_m120[mode]<ZERO) {
    f0=_r20[mode]/(1.-q2/_mfit20[mode]);
  }
  else if(_mfit20[mode]<ZERO) {
    f0=(_r10[mode]+_r20[mode]/(1.-q2/_m120[mode]))/(1.-q2/_m120[mode]);
  }
  else {
    f0=_r10[mode]/(1.-q2/_m120[mode])+_r20[mode]/(1.-q2/_mfit20[mode]);
  }
  // the F_1 form-factor
  if(_m12plus[mode]<ZERO) {
    fp = _r2plus[mode]/(1.-q2/_mfit2plus[mode]);
  }
  else if(_mfit2plus[mode]<ZERO) {
    fp = (_r1plus[mode]+_r2plus[mode]/(1.-q2/_m12plus[mode]))/(1.-q2/_m12plus[mode]);
  }
  else {
    fp =_r1plus[mode]/(1.-q2/_m12plus[mode])+_r2plus[mode]/(1.-q2/_mfit2plus[mode]);
  }
  if(id1==ParticleID::eta) {
    double fact(cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.));
    fp *= fact;
    f0 *= fact;
  }
}

void BallZwickyScalarFormFactor::ScalarScalarSigmaFormFactor(Energy2 q2,
							     unsigned int mode,int,
							     int id1,Energy,
							     Energy,
							     Complex & fT) const {
  useMe();
  // the F_T form-factor
  if(_m12T[mode]<ZERO) {
    fT = _r2T[mode]/(1.-q2/_mfit2T[mode]);
  }
  else if(_mfit2T[mode]<ZERO) {
    fT = (_r1T[mode]+_r2T[mode]/(1.-q2/_m12T[mode]))/(1.-q2/_m12T[mode]);
  }
  else {
    fT =_r1T[mode]/(1.-q2/_m12T[mode])+_r2T[mode]/(1.-q2/_mfit2T[mode]);
  }
  if(id1==ParticleID::eta) {
    fT *=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
  }
}

void BallZwickyScalarFormFactor::dataBaseOutput(ofstream & output,bool header,
						bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BallZwickyScalarFormFactor "
		    << name() << " \n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":r_10 " << ix << " " << _r10[ix] << "\n";
      output << "newdef " << name() << ":r_20 " << ix << " " << _r20[ix] << "\n";
      output << "newdef " << name() << ":r_1plus " << ix << " " << _r1plus[ix] << "\n";
      output << "newdef " << name() << ":r_2plus " << ix << " " << _r2plus[ix] << "\n";
      output << "newdef " << name() << ":r_1T " << ix << " " << _r1T[ix] << "\n";
      output << "newdef " << name() << ":r_2T " << ix << " " << _r2T[ix] << "\n";
      output << "newdef " << name() << ":m_120 " 
	     << ix << " " << _m120[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":mfit20 " 
	     << ix << " " << _mfit20[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":m_12plus " 
	     << ix << " " << _m12plus[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":mfit2plus " 
	     << ix << " " << _mfit2plus[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":m_12T " 
	     << ix << " " << _m12T[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":mfit2T " 
	     << ix << " " << _mfit2T[ix]/GeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":r_10 " 
	     << ix << " " << _r10[ix] << "\n";
      output << "insert " << name() << ":r_20 " 
	     << ix << " " << _r20[ix] << "\n";
      output << "insert " << name() << ":r_1plus " 
	     << ix << " " << _r1plus[ix] << "\n";
      output << "insert " << name() << ":r_2plus " 
	     << ix << " " << _r2plus[ix] << "\n";
      output << "insert " << name() << ":r_1T " << ix << " " << _r1T[ix] << "\n";
      output << "insert " << name() << ":r_2T " << ix << " " << _r2T[ix] << "\n";
      output << "insert " << name() << ":m_120 " 
	     << ix << " " << _m120[ix]/GeV2 << "\n";
      output << "insert " << name() << ":mfit20 " 
	     << ix << " " << _mfit20[ix]/GeV2 << "\n";
      output << "insert " << name() << ":m_12plus " 
	     << ix << " " << _m12plus[ix]/GeV2 << "\n";
      output << "insert " << name() << ":mfit2plus " 
	     << ix << " " << _mfit2plus[ix]/GeV2 << "\n";
      output << "insert " << name() << ":m_12T " 
	     << ix << " " << _m12T[ix]/GeV2 << "\n";
      output << "insert " << name() << ":mfit2T " 
	     << ix << " " << _mfit2T[ix]/GeV2 << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
