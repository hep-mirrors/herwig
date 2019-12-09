// -*- C++ -*-
//
// O2AlphaS.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the O2AlphaS class.
//

#include "O2AlphaS.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void O2AlphaS::persistentOutput(PersistentOStream & os) const {
  os << ounit(_lambdaQCD,GeV) << _bcoeff << _ccoeff << ounit(_lambdas,GeV) 
     << ounit(_threshold,GeV) << _match << _copt;
}

void O2AlphaS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_lambdaQCD,GeV) >> _bcoeff >> _ccoeff >> iunit(_lambdas,GeV) 
     >> iunit(_threshold,GeV) >> _match >> _copt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<O2AlphaS,AlphaSBase>
describeHerwigO2AlphaS("Herwig::O2AlphaS", "Herwig.so");

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
  vector<Energy2> thresholds(_threshold.size());
  transform(_threshold.begin(), _threshold.end(),
	    thresholds.begin(),
	    sqr<Energy>);
  return thresholds;
}

void O2AlphaS::doinit() {
  // thresholds
  for ( int ix=1; ix<7; ++ix ) {
    if ( quarkMasses().empty() ) {
      tPDPtr p = getParticleData(ix);
      _threshold[ix-1] = p->mass();
    } else
      _threshold[ix-1] = quarkMasses()[ix-1];
  }
  // d is heavier than u, need to swap
  swap(_threshold[0],_threshold[1]);

  // beta function coefficients
  const double ca = generator()->standardModel()->Nc();
  const double cf = (sqr(ca)-1.)/2./ca;
  for(unsigned int ix=3;ix<7;++ix)
    {
      _bcoeff[ix-1]=(11.*ca-2.*ix)/(12.*Constants::pi);
      _ccoeff[ix-1]=(17.*sqr(ca)-ix*(5.*ca+3.*cf))/(24.*sqr(Constants::pi))/sqr(_bcoeff[ix-1]);
    }
  if(_copt==0)
    {
      double kfac(ca*(67./18.-sqr(Constants::pi)/6.)-25./9.);
      _lambdas[5]=_lambdaQCD*exp(kfac/(4.*Constants::pi*_bcoeff[4]))/sqrt(2.);
    }
  else{_lambdas[5]=_lambdaQCD;}
  // calculate the threshold matching
  double rho=2.*log(_threshold[5]/_lambdas[5]);
  double rat=log(rho)/rho;
  _match[5]=(_bcoeff[4]/(1.-_ccoeff[4]*rat)-_bcoeff[5]/(1.-_ccoeff[5]*rat))*rho;
  rho=2.*log(_threshold[4]/_lambdas[5]);
  rat=log(rho)/rho;
  _match[4]=(_bcoeff[4]/(1.-_ccoeff[4]*rat)-_bcoeff[3]/(1.-_ccoeff[3]*rat))*rho;
  rho=2.*log(_threshold[3]/_lambdas[5]);
  rat=log(rho)/rho;
  _match[3]=(_bcoeff[3]/(1.-_ccoeff[3]*rat)-_bcoeff[2]/(1.-_ccoeff[2]*rat))*rho
    +_match[4];
  // calculate the 4-flavour lambda
  _lambdas[4]=_lambdas[5]*pow(_threshold[4]/_lambdas[5],2./25.)*
    pow(2.*log(_threshold[4]/_lambdas[5]),963./14375.);
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
  _lambdas[3]=_lambdas[5]*exp(0.5*d35);
  AlphaSBase::doinit();
}

vector<Energy> O2AlphaS::LambdaQCDs() const
{
  vector<Energy> output(4,_lambdas[3]);
  output.push_back(_lambdas[4]);
  output.push_back(_lambdas[5]);
  output.push_back(_lambdas[5]);
  return output;
}



double O2AlphaS::value(Energy2 scale, const StandardModelBase &) const
{
  Energy rs=sqrt(scale);
  if(scale<sqr(_lambdas[5])) {
    generator()->logWarning(Exception()
			    << "O2AlphaS called with scale less than Lambda QCD "
			    << "scale = " << rs/MeV << " MeV and "
			    << "Lambda = " << _lambdas[5]/MeV << " MeV"
			    << Exception::warning);
    return 0.;
  }
  double rho=2.*log(rs/_lambdas[5]),rat(log(rho)/rho);
  double rlf;
  if(rs>_threshold[5])      rlf=_bcoeff[5]*rho/(1.-_ccoeff[5]*rat)+_match[5];
  else if(rs>_threshold[4]) rlf=_bcoeff[4]*rho/(1.-_ccoeff[4]*rat);
  else if(rs>_threshold[3]) rlf=_bcoeff[3]*rho/(1.-_ccoeff[3]*rat)+_match[4];
  else                      rlf=_bcoeff[2]*rho/(1.-_ccoeff[2]*rat)+_match[3];
  // must be possible
  if(rlf<=0.) {
    generator()->logWarning(Exception() << "O2AlphaS coupling less than zero"
			    << Exception::warning) ;
    return 0.;
  }
  return 1./rlf;
}


