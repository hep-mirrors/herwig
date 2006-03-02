// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the O2AlphaS class.
//

#include "O2AlphaS.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "O2AlphaS.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

O2AlphaS::~O2AlphaS() {}

void O2AlphaS::persistentOutput(PersistentOStream & os) const {
  os << _lambdaQCD << _copt;
  for(unsigned int ix=0;ix<6;++ix)
    {os << _bcoeff[ix] << _ccoeff[ix] << _lambdas[ix] 
	<< _threshold[ix] << _match[ix];}
}

void O2AlphaS::persistentInput(PersistentIStream & is, int) {
  is >> _lambdaQCD >> _copt;
  for(unsigned int ix=0;ix<6;++ix)
    {is >> _bcoeff[ix] >> _ccoeff[ix] >> _lambdas[ix] 
	>> _threshold[ix] >> _match[ix];}
}

ClassDescription<O2AlphaS> O2AlphaS::initO2AlphaS;
// Definition of the static class description member.

void O2AlphaS::Init() {

  static ClassDocumentation<O2AlphaS> documentation
    ("The O2AlphaS class implements the next-to-leading order alphaS in the same"
     " way as in FORTRAN HERWIG");

  static Parameter<O2AlphaS,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "The value of Lambda QCD",
     &O2AlphaS::_lambdaQCD, MeV, 180.*MeV, 50.*MeV, 500.0*MeV,
     false, false, Interface::limited);


  static Switch<O2AlphaS,unsigned int> interfaceLambdaType
    ("LambdaType",
     "Which type of Lambda to use",
     &O2AlphaS::_copt, 0, false, false);
  static SwitchOption interfaceLambdaTypeMonteCarlo
    (interfaceLambdaType,
     "MonteCarlo",
     "Use a Monte Carlo scheme as in the FORTRAN program",
     0);
  static SwitchOption interfaceLambdaTypeMSBar
    (interfaceLambdaType,
     "MSBar",
     "Use the MSBar scheme",
     1);
}


vector<Energy2> O2AlphaS::flavourThresholds() const {
  vector<Energy2> thresholds;
  for(unsigned int ix=0;ix<6;++ix)
    {thresholds.push_back(sqr(_threshold[ix]));}
  return thresholds;
}

void O2AlphaS::doinit() throw(InitException) {
  AlphaSBase::doinit();
  // thresholds
  for(unsigned int ix=1;ix<7;++ix)
    {
      PDPtr p = getParticleData(ix);
      _threshold[ix-1]=p->mass();
    }
  // beta function coefficients
  double ca(generator()->standardModel()->Nc()),cf((sqr(ca)-1.)/2./ca);
  cerr << "testing ca and cf " << ca << " " << cf << endl;
  for(unsigned int ix=3;ix<7;++ix)
    {
      _bcoeff[ix-1]=(11.*ca-2.*ix)/(12.*pi);
      _ccoeff[ix-1]=(17.*sqr(ca)-ix*(5.*ca+3.*cf))/(24.*sqr(pi))/sqr(_bcoeff[ix-1]);
    }
  if(_copt==0)
    {
      double kfac(ca*(67./18.-sqr(pi)/6.)-25./9.);
      _lambdas[4]=_lambdaQCD*exp(kfac/(4.*pi*_bcoeff[4]))/sqrt(2.);
    }
  else{_lambdas[4]=_lambdaQCD;}
  cerr << "testing lambda " << _lambdas[4] << endl;
  cerr << "testing b " 
       << _bcoeff[2] << " " << _bcoeff[3] << " " 
       << _bcoeff[4] << " " << _bcoeff[5] << endl;
  cerr << "testing c " 
       << _ccoeff[2] << " " << _ccoeff[3] << " " 
       << _ccoeff[4] << " " << _ccoeff[5] << endl;
  // calculate the threshold matching
  double rho=2.*log(_threshold[5]/_lambdas[4]);
  double rat=log(rho)/rho;
  _match[5]=(_bcoeff[4]/(1.-_ccoeff[4]*rat)-_bcoeff[5]/(1.-_ccoeff[5]*rat))*rho;
  rho=2.*log(_threshold[4]/_lambdas[4]);
  rat=log(rho)/rho;
  _match[4]=(_bcoeff[4]/(1.-_ccoeff[4]*rat)-_bcoeff[3]/(1.-_ccoeff[3]*rat))*rho;
  rho=2.*log(_threshold[3]/_lambdas[4]);
  rat=log(rho)/rho;
  _match[3]=(_bcoeff[3]/(1.-_ccoeff[3]*rat)-_bcoeff[2]/(1.-_ccoeff[2]*rat))*rho
    +_match[4];
  cerr << "testing matching parameter " 
       << _match[3] << " " << _match[4] << " " 
       << _match[5] << endl;
  // calculate the 4-flavour lambda
  _lambdas[3]=_lambdas[4]*pow(_threshold[4]/_lambdas[4],2./25.)*
    pow(2.*log(_threshold[4]/_lambdas[4]),963./14375.);
  // calculate the 3-flavour lambda
  double eps(1.e-6),d35(-1./(_bcoeff[2]*_match[3])),rlf,drh;
  unsigned int ix=0;
  do
    {
      rat=log(d35)/d35;
      rlf=_bcoeff[2]*d35/(1.-_ccoeff[2]*rat);
      drh=_bcoeff[2]*(rlf+_match[3])*sqr(d35)/
	((1.-2.*_ccoeff[2]*rat+_ccoeff[2]/d35)*sqr(rlf));
      d35=d35-drh;
      ++ix;
    }
  while(ix<100&&abs(drh)>eps*d35);
  _lambdas[2]=_lambdas[4]*exp(0.5*d35);

  Energy2 scale(1.*GeV2);
  for(;scale<100000*GeV2;scale+=50.*GeV2)
    {
      cout << sqrt(scale) << " " << value(scale,*(generator()->standardModel())) 
  	   << endl;
    }
  //cerr << "testing " << 20*GeV << "  " 
  //     << value(400.*GeV2,*(generator()->standardModel())) << endl;
}

vector<Energy> O2AlphaS::LambdaQCDs() const
{
  vector<Energy> output(3,_lambdas[2]);
  output.push_back(_lambdas[3]);
  output.push_back(_lambdas[4]);
  output.push_back(_lambdas[4]);
  return output;
}



double O2AlphaS::value(Energy2 scale, const StandardModelBase & sm) const
{
  Energy rs=sqrt(scale);
  if(scale<_lambdas[4])
    {
      cerr << "O2AlphaS called with scale less than Lambda QCD "
			 << "scale = " << rs << " MeV and "
			 << "Lambda = " << _lambdas[4] << " MeV\n";
      generator()->log() << "O2AlphaS called with scale less than Lambda QCD "
			 << "scale = " << rs << " MeV and "
			 << "Lambda = " << _lambdas[4] << " MeV\n";
      return 0.;
    }
  double rho=2.*log(rs/_lambdas[4]),rat(log(rho)/rho);
  double rlf;
  if(rs>_threshold[5])
    {
      //      cerr << "testing did A " 
      //	   << _bcoeff[5] << " " << rho << " " 
      //	   << (1.-_ccoeff[5]*rat) << " " << _match[5] << endl;
      rlf=_bcoeff[5]*rho/(1.-_ccoeff[5]*rat)+_match[5];
    }
  else if(rs>_threshold[4])
    {
      //      cerr << "testing did B " << _bcoeff[4] << " " <<rho << " " 
      //	   <<(1.-_ccoeff[4]*rat) << endl;
      rlf=_bcoeff[4]*rho/(1.-_ccoeff[4]*rat);
    }
  else if(rs>_threshold[3])
    {
      //      cerr << "testing did C " << _bcoeff[3] << " " <<rho << " " <<(1.-_ccoeff[3]*rat) << " " <<_match[4] << endl;
      rlf=_bcoeff[3]*rho/(1.-_ccoeff[3]*rat)+_match[4];
    }
  else
    {
      //      cerr << "testing did D " << _bcoeff[2] << " " <<rho << " " <<(1.-_ccoeff[2]*rat) << " " << _match[3] << endl;
      rlf=_bcoeff[2]*rho/(1.-_ccoeff[2]*rat)+_match[3];
    }
  // must be possible
  if(rlf<=0.)
    {
      generator()->log() << "O2AlphaS coupling less than zero \n";
      return 0.;
    }
  return 1./rlf;
}


