// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyScalarFormFactor class.
//

#include "BallZwickyScalarFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BallZwickyScalarFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {
using namespace ThePEG;

BallZwickyScalarFormFactor::BallZwickyScalarFormFactor() 
{
  double ort(1./sqrt(2.));
  // parameters for the B to pi  form-factors
  // B+ to pi0
  _r10.push_back(0.);_r20.push_back(ort*0.258);
  _m120.push_back(-1.*GeV2);_mfit20.push_back(33.81*GeV2);
  _r1plus.push_back(ort*0.744);_r2plus.push_back(-ort*0.486);
  _m12plus.push_back(5.32*5.32*GeV2);_mfit2plus.push_back(40.73*GeV2);
  _r1T.push_back(ort*1.387);_r2T.push_back(-ort*1.134);
  _m12T.push_back(5.32*5.32*GeV2);_mfit2T.push_back(32.22*GeV2);
  addFormFactor(521, 111,0,2,-5,-2);
  // B0 to pi+
  _r10.push_back(0.);_r20.push_back(0.258);
  _m120.push_back(-1.*GeV2);_mfit20.push_back(33.81*GeV2);
  _r1plus.push_back(0.744);_r2plus.push_back(-0.486);
  _m12plus.push_back(5.32*5.32*GeV2);_mfit2plus.push_back(40.73*GeV2);
  _r1T.push_back(1.387);_r2T.push_back(-1.134);
  _m12T.push_back(5.32*5.32*GeV2);_mfit2T.push_back(32.22*GeV2);
  addFormFactor(511,-211,0,2,-5,-2);
  // parameters for the B to K   form-factors
  // B+ to K+
  _r10.push_back(0.);_r20.push_back(0.330);
  _m120.push_back(-1.*GeV2);_mfit20.push_back(37.46*GeV2);
  _r1plus.push_back(0.162);_r2plus.push_back(0.173);
  _m12plus.push_back(5.41*5.41*GeV2);_mfit2plus.push_back(-1.*GeV2);
  _r1T.push_back(0.161);_r2T.push_back(0.198);
  _m12T.push_back(5.41*5.41*GeV2);_mfit2T.push_back(-1.*GeV2);
  addFormFactor(521,321,0,2,-5,-3);
  // B0 to K0
  _r10.push_back(0.);_r20.push_back(0.330);
  _m120.push_back(-1.*GeV2);_mfit20.push_back(37.46*GeV2);
  _r1plus.push_back(0.162);_r2plus.push_back(0.173);
  _m12plus.push_back(5.41*5.41*GeV2);_mfit2plus.push_back(-1.*GeV2);
  _r1T.push_back(0.161);_r2T.push_back(0.198);
  _m12T.push_back(5.41*5.41*GeV2);_mfit2T.push_back(-1.*GeV2);
  addFormFactor(511,311,0,2,-5,-3);
  // parameters for the B to eta form-factors
  _r10.push_back(0.);_r20.push_back(ort*0.273);
  _m120.push_back(-1.*GeV2);_mfit20.push_back(31.03*GeV2);
  _r1plus.push_back(ort*0.122);_r2plus.push_back(ort*0.155);
  _m12plus.push_back(5.32*5.32*GeV2);_mfit2plus.push_back(-1.*GeV2);
  _r1T.push_back(ort*0.111);_r2T.push_back(ort*0.175);
  _m12T.push_back(5.32*5.32*GeV2);_mfit2T.push_back(-1.*GeV2);
  addFormFactor(521,221,0,2,-5,-2);
  // initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta = -pi/9.;
}

BallZwickyScalarFormFactor::~BallZwickyScalarFormFactor() {}

void BallZwickyScalarFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _r10 << _r20 << _r1plus << _r2plus << _r1T << _r2T << _m120 << _mfit20 
     << _m12plus << _mfit2plus << _m12T << _mfit2T << _thetaeta;
}

void BallZwickyScalarFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _r10 >> _r20 >> _r1plus >> _r2plus >> _r1T >> _r2T >> _m120 >> _mfit20 
     >> _m12plus >> _mfit2plus >> _m12T >> _mfit2T >> _thetaeta;
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
     &BallZwickyScalarFormFactor::_thetaeta, -pi/9., -pi, pi,
     false, false, true);
}


// form-factor for scalar to scalar
void BallZwickyScalarFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned  int mode,
							int id0,
							int id1,Energy m0, Energy m1,
						 	Complex & f0, Complex & fp) const
{
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
  if(id0==ParticleID::eta)
    {
      double fact(sqrt(1./3.)*cos(_thetaeta)-sqrt(2./3.)*sin(_thetaeta));
      fp *=fact;f0*=fact;
    }
}

void BallZwickyScalarFormFactor::ScalarScalarSigmaFormFactor(Energy2 q2,
							     unsigned int mode,int id0,
							     int id1,Energy m0,
							     Energy m1,
							     Complex & fT) const
{
  // the F_T form-factor
  if(_m12T[mode]<0)
    {fT = _r2T[mode]/(1.-q2/_mfit2T[mode]);}
  else if(_mfit2T[mode]<0)
    {fT = (_r1T[mode]+_r2T[mode]/(1.-q2/_m12T[mode]))/(1.-q2/_m12T[mode]);}
  else
    {fT =_r1T[mode]/(1.-q2/_m12T[mode])+_r2T[mode]/(1.-q2/_mfit2T[mode]);}
  if(id0==ParticleID::eta)
    {fT *=sqrt(1./3.)*cos(_thetaeta)-sqrt(2./3.)*sin(_thetaeta);}
}

void BallZwickyScalarFormFactor::dataBaseOutput(ofstream & output) const
{
  output << "create Herwig++::BallZwickyScalarFormFactor " << fullName() << " \n";
  output << "set " << fullName() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      if(ix<initialModes())
	{
	  output << "set " << fullName() << ":r_10 " << ix << " " << _r10[ix] << "\n";
	  output << "set " << fullName() << ":r_20 " << ix << " " << _r20[ix] << "\n";
	  output << "set " << fullName() << ":r_1plus " << ix << " " << _r1plus[ix] << "\n";
	  output << "set " << fullName() << ":r_2plus " << ix << " " << _r2plus[ix] << "\n";
	  output << "set " << fullName() << ":r_1T " << ix << " " << _r1T[ix] << "\n";
	  output << "set " << fullName() << ":r_2T " << ix << " " << _r2T[ix] << "\n";
	  output << "set " << fullName() << ":m_120 " 
		 << ix << " " << _m120[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":mfit20 " 
		 << ix << " " << _mfit20[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":m_12plus " 
		 << ix << " " << _m12plus[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":mfit2plus " 
		 << ix << " " << _mfit2plus[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":m_12T " 
		 << ix << " " << _m12T[ix]/GeV2 << "\n";
	  output << "set " << fullName() << ":mfit2T " 
		 << ix << " " << _mfit2T[ix]/GeV2 << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":r_10 " 
		 << ix << " " << _r10[ix] << "\n";
	  output << "insert " << fullName() << ":r_20 " 
		 << ix << " " << _r20[ix] << "\n";
	  output << "insert " << fullName() << ":r_1plus " 
		 << ix << " " << _r1plus[ix] << "\n";
	  output << "insert " << fullName() << ":r_2plus " 
		 << ix << " " << _r2plus[ix] << "\n";
	  output << "insert " << fullName() << ":r_1T " << ix << " " << _r1T[ix] << "\n";
	  output << "insert " << fullName() << ":r_2T " << ix << " " << _r2T[ix] << "\n";
	  output << "insert " << fullName() << ":m_120 " 
		 << ix << " " << _m120[ix]/GeV2 << "\n";
	  output << "insert " << fullName() << ":mfit20 " 
		 << ix << " " << _mfit20[ix]/GeV2 << "\n";
	  output << "insert " << fullName() << ":m_12plus " 
		 << ix << " " << _m12plus[ix]/GeV2 << "\n";
	  output << "insert " << fullName() << ":mfit2plus " 
		 << ix << " " << _mfit2plus[ix]/GeV2 << "\n";
	  output << "insert " << fullName() << ":m_12T " 
		 << ix << " " << _m12T[ix]/GeV2 << "\n";
	  output << "insert " << fullName() << ":mfit2T " 
		 << ix << " " << _mfit2T[ix]/GeV2 << "\n";
	}
    }
  ScalarFormFactor::dataBaseOutput(output);
}

}
