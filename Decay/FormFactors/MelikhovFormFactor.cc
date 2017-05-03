// -*- C++ -*-
//
// MelikhovFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MelikhovFormFactor class.
//

#include "MelikhovFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

MelikhovFormFactor::MelikhovFormFactor() : 
  _ifit(1), _Rplus0(0), _Mplus(ZERO), _nplus(0), _RV0(0), _MV(ZERO), _nV(0),
  _R10(0), _M1(ZERO), _n1(0), _R20(0), _M2(ZERO), _n2(0) {
  // the possible modes
  // B to rho
  addFormFactor(-521, 113,1,2,5,2);
  addFormFactor(-511, 213,1,1,5,2);
  addFormFactor(-521,-213,1,2,5,1);
  addFormFactor(-511, 113,1,1,5,1);
  // B to pi
  addFormFactor(-521, 111,0,2,5,2);
  addFormFactor(-511, 211,0,1,5,2);
  addFormFactor(-521,-211,0,2,5,1);
  addFormFactor(-511, 111,0,1,5,1);
  // set the initial number of modes
  initialModes(numberOfFactors());
}

void MelikhovFormFactor::doinit() {
  ScalarFormFactor::doinit();
  // the parameters for the different fits
  double Rplus0[4]={0.29    ,0.20    ,0.21    ,0.26    };
  Energy Mplus[4] ={6.29*GeV,6.22*GeV,5.90*GeV,5.44*GeV};
  double nplus[4] ={2.35    ,2.45    ,2.33    ,1.72    };
  double RV[4]    ={0.30    ,0.20    ,0.21    ,0.29    };
  Energy MV[4]    ={6.28*GeV,6.22*GeV,5.90*GeV,5.46*GeV};
  double nv[4]    ={2.36    ,2.46    ,2.35    ,1.73    };
  double R10[4]   ={0.27    ,0.20    ,0.21    ,0.29    };
  Energy M1[4]    ={7.07*GeV,6.78*GeV,6.50*GeV,5.68*GeV};
  double n1[4]    ={2.65    ,2.65    ,2.70    ,1.67    };
  double R20[4]   ={0.25    ,0.19    ,0.20    ,0.28    };
  Energy M2[4]    ={6.13*GeV,6.00*GeV,5.90*GeV,5.36*GeV};
  double n2[4]    ={2.17    ,2.34    ,2.45    ,1.67    };
  // set the values
  _Rplus0=Rplus0[_ifit-1];
  _Mplus=Mplus[_ifit-1];
  _nplus=nplus[_ifit-1];
  _RV0=RV[_ifit-1];
  _MV=MV[_ifit-1];
  _nV=nv[_ifit-1];
  _R10=R10[_ifit-1];
  _M1=M1[_ifit-1];
  _n1=n1[_ifit-1];
  _R20=R20[_ifit-1];
  _M2=M2[_ifit-1];
  _n2=n2[_ifit-1];
}

void MelikhovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _ifit << _Rplus0 << ounit(_Mplus,GeV) << _nplus << _RV0 << ounit(_MV,GeV) 
     << _nV << _R10 << ounit(_M1,GeV) << _n1 << _R20 << ounit(_M2,GeV) << _n2; 
}

void MelikhovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _ifit >> _Rplus0 >> iunit(_Mplus,GeV) >> _nplus >> _RV0 >> iunit(_MV,GeV) 
     >> _nV >> _R10 >> iunit(_M1,GeV) >> _n1 >> _R20 >> iunit(_M2,GeV) >> _n2; 
}

ClassDescription<MelikhovFormFactor> MelikhovFormFactor::initMelikhovFormFactor;
// Definition of the static class description member.

void MelikhovFormFactor::Init() {

  static ClassDocumentation<MelikhovFormFactor> documentation
    ("The MelikhovFormFactor class implements the form factors from"
     " Phys. Lett. B 380 (1996) 363 for B to pi,rho",
     "The form factors of \\cite{Melikhov:1996ge} were used for $B\\to \\pi,\\rho$.",
     "\\bibitem{Melikhov:1996ge} D.~Melikhov, Phys.\\ Lett.\\  B {\\bf 380} (1996) 363\n"
     "[arXiv:hep-ph/9603340]. %%CITATION = PHLTA,B380,363;%%\n");

  static Parameter<MelikhovFormFactor,unsigned int> interfaceFit
    ("Fit",
     "Which of the fits from hep-ph/9603340 to use",
     &MelikhovFormFactor::_ifit, 1, 1, 4,
     false, false, true);

}

// form-factor for scalar to scalar
void MelikhovFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
						int,int id1,
						Energy, Energy,Complex & f0,
						Complex & fp) const {
  useMe();
  fp = _Rplus0/pow((1.-q2/_Mplus/_Mplus),_nplus);
  f0 = fp;
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)) {
    double fact = (id1==ParticleID::pi0&&abs(outquark)==1) ? -sqrt(0.5) : sqrt(0.5);
    f0 *= fact;
    fp *= fact;
  }
}

void MelikhovFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
						int, int id1, Energy mY, Energy mX,
						Complex & A0,Complex & A1,Complex & A2,
						Complex & V) const {
  useMe();
  // constants
  double r(mX/mY),y(q2/mY/mY);
  // V form-factor
  V  =-double((1.+r)*_RV0/pow((1.-q2/_MV/_MV),_nV));
  // A_1 form-factor
  A1 = (1.+r*r-y)/(1.+r)*_R10/pow((1.-q2/_M1/_M1),_n1);
  // A_2 form-factor
  A2 = (1.+r)*(1.-r*r-y)/((1.+r)*(1.+r)-y)*_R20/pow((1.-q2/_M2/_M2),_n2);
  // set the a_- factor to zero
//   A0 = 0.5/mX*((mY+mX)*A1-(mY-mX)*A2);
  A0 = 0.;
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3) {
    double fact = (id1==ParticleID::rho0&&abs(outquark)==1) ? -sqrt(0.5) : sqrt(0.5);
    A0 *= fact;
    A1 *= fact;
    A2 *= fact;
    V  *= fact;
  }
}

void MelikhovFormFactor::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::MelikhovFormFactor " << name() << " \n";
  output << "newdef " << name() << ":Fit " << _ifit << " \n";
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
