// -*- C++ -*-
//
// WSBFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WSBFormFactor class.
//

#include "WSBFormFactor.h"
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

void WSBFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_F0.size() ||isize!=_V.size()  ||isize!=_A0.size() ||
     isize!=_A1.size() ||isize!=_A2.size() ||isize!=_mS0.size()||
     isize!=_mS1.size()||isize!=_mV0.size()||isize!=_mV1.size())
    throw InitException() << "Inconsistent parameters in WSBFormFactor::doinit()" 
			  << Exception::abortnow;
}


WSBFormFactor::WSBFormFactor() 
  : _F0(51), _V(51), _A0(51), _A1(51), _A2(51), 
    _mS0(51), _mS1(51), _mV0(51), _mV1(51) {
  // modes handled by this and the parameters model
  // K to pi 
  addFormFactor(-321, 111,0,-2,3,2);
  addFormFactor(-311, 211,0,-1,3,2);
  addFormFactor(-321,-211,0,-2,3,1);
  addFormFactor(-311, 111,0,-1,3,1);
  for(unsigned int ix=0;ix<4;++ix) {
    _F0[ix] = 0.992; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 0.494*GeV; _mV0[ix] = 0.892*GeV; 
    _mS1[ix] = 1.430*GeV; _mV1[ix] = 1.273*GeV; 
  }
  // D to K 
  addFormFactor(421,-321,0,-2,4,3);
  addFormFactor(411,-311,0,-1,4,3);
  for(unsigned int ix=4;ix<6;++ix) {
    _F0[ix] = 0.762; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.97*GeV; _mV0[ix] = 2.11*GeV; 
    _mS1[ix] = 2.60*GeV; _mV1[ix] = 2.53*GeV;
  }
  // D to pi
  addFormFactor(421,-211,0,-2,4,1);
  addFormFactor(421, 111,0,-2,4,2);
  addFormFactor(411, 111,0,-1,4,1);
  addFormFactor(411, 211,0,-1,4,2);
  for(unsigned int ix=6;ix<10;++ix) {
    _F0[ix] = 0.692; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to eta
  addFormFactor(421,221,0,-2,4,2);
  addFormFactor(411,221,0,-1,4,1);
  for(unsigned int ix=10;ix<12;++ix) {
    _F0[ix] = 0.681; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to eta'
  addFormFactor(421,331,0,-2,4,2);
  addFormFactor(411,331,0,-1,4,1);
  for(unsigned int ix=12;ix<14;++ix) {
    _F0[ix] = 0.655; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to K*
  addFormFactor(421,-323,1,-2,4,3);
  addFormFactor(411,-313,1,-1,4,3);
  for(unsigned int ix=14;ix<16;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.226; _A0[ix] = 0.733; 
    _A1[ix] = 0.880; _A2[ix] = 1.147; 
    _mS0[ix] = 1.97*GeV; _mV0[ix] = 2.11*GeV; 
    _mS1[ix] = 2.60*GeV; _mV1[ix] = 2.53*GeV;
  } 
  // D to rho
  addFormFactor(421,-213,1,-2,4,1);
  addFormFactor(421, 113,1,-2,4,2);
  addFormFactor(411, 113,1,-1,4,1);
  addFormFactor(411, 213,1,-1,4,2);
  for(unsigned int ix=16;ix<20;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.225; _A0[ix] = 0.669; 
    _A1[ix] = 0.775; _A2[ix] = 0.923; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to omega
  addFormFactor(411,223,1,-1,4,1); 
  addFormFactor(421,223,1,-2,4,2); 
  for(unsigned int ix=20;ix<22;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.236; _A0[ix] = 0.669; 
    _A1[ix] = 0.772; _A2[ix] = 0.920; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D_s to eta
  addFormFactor(431,221,0,-3,4,3);
  _F0[22] = 0.723; _V[22]  = 0.000; _A0[22] = 0.000; 
  _A1[22] = 0.000; _A2[22] = 0.000; 
  _mS0[22] = 1.97*GeV; _mV0[22] = 2.11*GeV; 
  _mS1[22] = 2.60*GeV; _mV1[22] = 2.53*GeV; 
  // D_s to eta'
  addFormFactor(431,331,0,-3,4,3);
  _F0[23] = 0.704; _V[23]  = 0.000; _A0[23] = 0.000; 
  _A1[23] = 0.000; _A2[23] = 0.000; 
  _mS0[23] = 1.97*GeV; _mV0[23] = 2.11*GeV; 
  _mS1[23] = 2.60*GeV; _mV1[23] = 2.53*GeV; 
  // D_s to K
  addFormFactor(431,311,0,-3,4,1);
  addFormFactor(431,321,0,-3,4,2);
  for(unsigned int ix=24;ix<26;++ix) {
    _F0[ix] = 0.643; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D_s to K*
  addFormFactor(431,313,1,-3,4,1);
  addFormFactor(431,323,1,-3,4,2);
  for(unsigned int ix=26;ix<28;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.250; _A0[ix] = 0.634; 
    _A1[ix] = 0.717; _A2[ix] = 0.853; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D_s to phi
  addFormFactor(431,333,1,-3,4,3);
  _F0[28] = 0.000; _V[28]  = 1.319; _A0[28] = 0.700; 
  _A1[28] = 0.820; _A2[28] = 1.076; 
  _mS0[28] = 1.97*GeV; _mV0[28] = 2.11*GeV; 
  _mS1[28] = 2.60*GeV; _mV1[28] = 2.53*GeV; 
  // B to D
  addFormFactor(-521,421,0,-2,5,4);
  addFormFactor(-511,411,0,-2,5,4);
  for(unsigned int ix=29;ix<31;++ix) {
    _F0[ix] = 0.690; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 6.30*GeV; _mV0[ix] = 6.34*GeV; 
    _mS1[ix] = 6.80*GeV; _mV1[ix] = 6.73*GeV;
  }
  // B to K 
  addFormFactor(-521,-321,0,-2,5,3);
  addFormFactor(-511,-311,0,-1,5,3);
  for(unsigned int ix=31;ix<33;++ix) {
    _F0[ix] = 0.379; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.38*GeV; _mV0[ix] = 5.43*GeV; 
    _mS1[ix] = 5.89*GeV; _mV1[ix] = 5.82*GeV;
  }
  // B to pi
  addFormFactor(-521, 111,0,-2,5,2);
  addFormFactor(-511, 211,0,-1,5,2);
  addFormFactor(-521,-211,0,-2,5,1);
  addFormFactor(-511, 111,0,-1,5,1);
  for(unsigned int ix=33;ix<37;++ix) {
    _F0[ix] = 0.333; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to eta
  addFormFactor(-521,221,0,-2,5,2);
  addFormFactor(-511,221,0,-1,5,1);
  for(unsigned int ix=37;ix<39;++ix) {
    _F0[ix] = 0.307; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to eta'
  addFormFactor(-521,331,0,-2,5,2);
  addFormFactor(-511,331,0,-1,5,1);
  for(unsigned int ix=39;ix<41;++ix) {
    _F0[ix] = 0.254; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to D*
  addFormFactor(-521,423,1,-2,5,4);
  addFormFactor(-511,413,1,-1,5,4);
  for(unsigned int ix=41;ix<43;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 0.705; _A0[ix] = 0.623; 
    _A1[ix] = 0.651; _A2[ix] = 0.686; 
    _mS0[ix] = 6.30*GeV; _mV0[ix] = 6.34*GeV; 
    _mS1[ix] = 6.80*GeV; _mV1[ix] = 6.73*GeV;
  }
  // B to K* 
  addFormFactor(-521,-323,1,-2,5,3);
  addFormFactor(-511,-313,1,-1,5,3);
  for(unsigned int ix=43;ix<45;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 0.369; _A0[ix] = 0.321; 
    _A1[ix] = 0.328; _A2[ix] = 0.331; 
    _mS0[ix] = 5.38*GeV; _mV0[ix] = 5.43*GeV; 
    _mS1[ix] = 5.89*GeV; _mV1[ix] = 5.82*GeV;
  }
  // B to rho
  addFormFactor(-521, 113,1,-2,5,2);
  addFormFactor(-511, 213,1,-1,5,2);
  addFormFactor(-521,-213,1,-2,5,1);
  addFormFactor(-511, 113,1,-1,5,1); 
  for(unsigned int ix=45;ix<49;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 0.329; _A0[ix] = 0.281; 
    _A1[ix] = 0.283; _A2[ix] = 0.283; 
 _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to omega
  addFormFactor(-521,223,1,-2,5,2); 
  addFormFactor(-511,223,1,-1,5,1);
  for(unsigned int ix=49;ix<51;++ix) {
    _F0[ix] = 0.000; _V[ix] = 0.328; _A0[ix] = 0.280; 
    _A1[ix] = 0.281; _A2[ix] = 0.281; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  } 
  // set the initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta=-0.194;
}

void WSBFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _F0 << _V << _A0 << _A1 << _A2 << ounit(_mS0,GeV) 
     << ounit(_mS1,GeV) << ounit(_mV0,GeV) << ounit(_mV1,GeV) << _thetaeta;
}

void WSBFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _F0 >> _V >> _A0 >> _A1 >> _A2 >> iunit(_mS0,GeV) 
     >> iunit(_mS1,GeV) >> iunit(_mV0,GeV) >> iunit(_mV1,GeV) >> _thetaeta;
}

ClassDescription<WSBFormFactor> WSBFormFactor::initWSBFormFactor;
// Definition of the static class description member.

void WSBFormFactor::Init() {

  static ClassDocumentation<WSBFormFactor> documentation
    ("The WSBFormFactor class is the implementation of the form-factors of "
     "Z.Phys.C29,637.",
     "The form factor model of \\cite{Bauer:1986bm,Wirbel:1985ji} was used"
     " for either semi-leptonic or hadronic weak decays",
     "\\bibitem{Bauer:1986bm} M.~Bauer, B.~Stech and M.~Wirbel,\n"
     "Z.\\ Phys.\\  C {\\bf 34} (1987) 103.\n"
     "%%CITATION = ZEPYA,C34,103;%%\n"
     "\\bibitem{Wirbel:1985ji} M.~Wirbel, B.~Stech and M.~Bauer,"
     "Z.\\ Phys.\\  C {\\bf 29} (1985) 637.\n"
     "%%CITATION = ZEPYA,C29,637;%%\n");

  static ParVector<WSBFormFactor,double> interfaceF0
    ("F0",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_F0,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceV
    ("V",
     "The form-factor V at zero q^2",
     &WSBFormFactor::_V,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA0
    ("A0",
     "The form-factor A0 at zero q^2",
     &WSBFormFactor::_A0,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA1
    ("A1",
     "The form-factor A1 at zero q^2",
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
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoScalarMass
    ("PseudoScalarMass",
     "The pseudoscalar mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mS1,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfaceVectorMass
    ("VectorMass",
     "The vector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV0,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoVectorMass
    ("PseudoVectorMass",
     "The pseudovector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV1,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static Parameter<WSBFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &WSBFormFactor::_thetaeta, -0.194, -Constants::pi, Constants::pi,
     false, false, true);
}

// form-factor for scalar to scalar
void WSBFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
					   int,int id1,
					   Energy, Energy,Complex & f0,
					   Complex & fp) const {
  useMe();
  f0 = _F0[mode]/(1.-q2/sqr(_mS1[mode]));
  fp = _F0[mode]/(1.-q2/sqr(_mV0[mode]));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)) {
    double fact;
    if(id1==ParticleID::eta) {
      if(abs(outquark)==3) fact = -2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
      else                 fact =     cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
    }
    else if(id1==ParticleID::etaprime) {
      if(abs(outquark)==3) fact = -2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
      else                 fact =     sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
    }
    else if(id1==ParticleID::pi0&&abs(outquark)==1) fact=-sqrt(0.5);
    else                                            fact= sqrt(0.5);
    f0*=fact;
    fp*=fact;
  }
}

void WSBFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
					   int, int id1, 
					   Energy, Energy,Complex & A0,
					   Complex & A1,Complex & A2,Complex & V) const {
  A0 = -_A0[mode]/(1.-q2/_mS0[mode]/_mS0[mode]);
  A1 = -_A1[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  A2 = -_A2[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  V  =   _V[mode]/(1.-q2/_mV0[mode]/_mV0[mode]);
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3) {
    double fact = id1==ParticleID::rho0&&abs(outquark)==1 ? -sqrt(0.5) : sqrt(0.5);
    A0 *= fact;
    A1 *= fact;
    A2 *= fact;
    V  *= fact;
  }
}

void WSBFormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const {
  useMe();
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::WSBFormFactor " << name() << " \n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":F0 " 
	     << ix << "  " << _F0[ix] << "\n";
      output << "newdef " << name() << ":V  " 
	     << ix << "  " << _V[ix]  << "\n";
      output << "newdef " << name() << ":A0 " 
	     << ix << "  " << _A0[ix] << "\n";
      output << "newdef " << name() << ":A1 " 
	     << ix << "  " << _A1[ix] << "\n";
      output << "newdef " << name() << ":A2 " 
	     << ix << "  " << _A2[ix] << "\n";
      output << "newdef " << name() << ":ScalarMass " 
	     << ix << "  " << _mS0[ix]/GeV << "\n";
      output << "newdef " << name() << ":PseudoScalarMass " 
	     << ix << "  " << _mS1[ix]/GeV << "\n";
      output << "newdef " << name() << ":VectorMass " 
	     << ix << "  " << _mV0[ix]/GeV << "\n";
      output << "newdef " << name() << ":PseudoVectorMass " 
	     << ix << "  " << _mV1[ix]/GeV << "\n";
    }
    else {
      output << "insert " << name() << ":F0 " 
	     << ix << "  " << _F0[ix] << "\n";
      output << "insert " << name() << ":V  " 
	     << ix << "  " << _V[ix]  << "\n";
      output << "insert " << name() << ":A0 " 
	     << ix << "  " << _A0[ix] << "\n";
      output << "insert " << name() << ":A1 " 
	     << ix << "  " << _A1[ix] << "\n";
      output << "insert " << name() << ":A2 " 
	     << ix << "  " << _A2[ix] << "\n";
      output << "insert " << name() << ":ScalarMass " 
	     << ix << "  " << _mS0[ix]/GeV << "\n";
      output << "insert " << name() << ":PseudoScalarMass " 
	     << ix << "  " << _mS1[ix]/GeV << "\n";
      output << "insert " << name() << ":VectorMass " 
	     << ix << "  " << _mV0[ix]/GeV << "\n";
      output << "insert " << name() << ":PseudoVectorMass " 
	     << ix << "  " << _mV1[ix]/GeV << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
