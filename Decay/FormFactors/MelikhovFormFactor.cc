// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MelikhovFormFactor class.
//

#include "MelikhovFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MelikhovFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void MelikhovFormFactor::doinit() throw(InitException) {
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

MelikhovFormFactor::~MelikhovFormFactor() {}

void MelikhovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _ifit << _Rplus0 << _Mplus << _nplus << _RV0 << _MV << _nV << _R10 << _M1
     << _n1 << _R20 << _M2 << _n2; 
}

void MelikhovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _ifit >> _Rplus0 >> _Mplus >> _nplus >> _RV0 >> _MV >> _nV >> _R10 >> _M1
     >> _n1 >> _R20 >> _M2 >> _n2; 
}

ClassDescription<MelikhovFormFactor> MelikhovFormFactor::initMelikhovFormFactor;
// Definition of the static class description member.

void MelikhovFormFactor::Init() {

  static ClassDocumentation<MelikhovFormFactor> documentation
    ("The MelikhovFormFactor class implements the form factors from"
     " hep-ph/9603340 for B to pi,rho");

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
						Complex & fp) const
{
  fp = _Rplus0/pow((1.-q2/_Mplus/_Mplus),_nplus);
  f0 = fp;
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect))
    {
      double fact(sqrt(0.5));
      if(id1==ParticleID::pi0&&abs(outquark)==1){fact=-fact;}
      f0*=fact;fp*=fact;
    }
}

void MelikhovFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
						int, int id1, Energy mY, Energy mX,
						Complex & A0,Complex & A1,Complex & A2,
						Complex & V) const
{
  // constants
  double r(mX/mY),y(q2/mY/mY);
  // V form-factor
  V  = (1.+r)*_RV0/pow((1.-q2/_MV/_MV),_nV);
  // A_1 form-factor
  A1 = (1.+r*r-y)/(1.+r)*_R10/pow((1.-q2/_M1/_M1),_n1);
  // A_2 form-factor
  A2 = (1.+r)*(1.-r*r-y)/((1.+r)*(1.+r)-y)*_R20/pow((1.-q2/_M2/_M2),_n2);
  // set the a_- factor to zero
  A0 = 0.5/mX*((mY+mX)*A1-(mY-mX)*A2);
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3)
    {
      double fact(sqrt(0.5));
      if(id1==ParticleID::rho0&&abs(outquark)==1){fact=-fact;}
      A0*=fact;A1*=fact;A2*=fact;V*=fact;
    }
}

void MelikhovFormFactor::dataBaseOutput(ofstream & output,bool header,
					bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::MelikhovFormFactor " << fullName() << " \n";}
  output << "set " << fullName() << ":Fit " << _ifit << " \n";
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
