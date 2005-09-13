// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LambdabExcitedLambdacSumRuleFormFactor class.
//

#include "LambdabExcitedLambdacSumRuleFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LambdabExcitedLambdacSumRuleFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

LambdabExcitedLambdacSumRuleFormFactor::~LambdabExcitedLambdacSumRuleFormFactor() {}

void LambdabExcitedLambdacSumRuleFormFactor::
persistentOutput(PersistentOStream & os) const {
  os << _xi1 << _rho2;
}

void LambdabExcitedLambdacSumRuleFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _xi1 >> _rho2;
}

ClassDescription<LambdabExcitedLambdacSumRuleFormFactor> LambdabExcitedLambdacSumRuleFormFactor::initLambdabExcitedLambdacSumRuleFormFactor;
// Definition of the static class description member.

void LambdabExcitedLambdacSumRuleFormFactor::Init() {

  static ClassDocumentation<LambdabExcitedLambdacSumRuleFormFactor> documentation
    ("The \\classname{LambdabExcitedLambdacSumRuleFormFactor} class implements the"
     " form-factors for Lambda_b to Lambda_c1(*) from hep-ph/0012114.");

  static Parameter<LambdabExcitedLambdacSumRuleFormFactor,double> interfaceXi
    ("Xi",
     "The intercept for the Isgur-Wise form-factor",
     &LambdabExcitedLambdacSumRuleFormFactor::_xi1, 0.29, 0.0, 10.0,
     false, false, true);

  static Parameter<LambdabExcitedLambdacSumRuleFormFactor,double> interfaceRho2
    ("Rho2",
     "The slope parameter for the form-factor.",
     &LambdabExcitedLambdacSumRuleFormFactor::_rho2, 2.01, -10.0, 10.0,
     false, false, true);

}

void LambdabExcitedLambdacSumRuleFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc,int id0,int id1,Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  double omega(.5/m0/m1*(m0*m0+m1*m1-q2)),orr(1./sqrt(3.));
  // the universal form-factor
  double xi=_xi1*(1.-_rho2*(omega-1.));
  // the couplings in the velocity form
  Complex g1v,g1a,g2v,g2a,g3a,g3v;
  g1v = orr*(omega-1.)*xi;
  g1a = orr*(omega+1.)*xi;
  g2v =-2.*orr*xi;
  g3v = 0.;
  g2a =-2.*orr*xi;
  g3a = 0.;
  // convert to our form
  f1a = g1v-0.5*(m0-m1)*(g2v/m0+g3v/m1);
  f1v =-g1a-0.5*(m0+m1)*(g2a/m0+g3a/m1);
  f2a = 0.5*(m0+m1)*( g2v/m0+g3v/m1);
  f2v =-0.5*(m0+m1)*( g2a/m0+g3a/m1);
  f3a = 0.5*(m0+m1)*( g2v/m0-g3v/m1);
  f3v = 0.5*(m0+m1)*(-g2a/m0+g3a/m1);
}

void  LambdabExcitedLambdacSumRuleFormFactor::
 SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int iloc,int id0,int id1,Energy m0,Energy m1,
				 Complex & f1v,Complex & f2v,
				 Complex & f3v,Complex & f4v,
				 Complex & f1a,Complex & f2a,
				 Complex & f3a,Complex & f4a )
{
  // the omega value
  double omega(.5/m0/m1*(m0*m0+m1*m1-q2));
  // the universal form-factor
  double xi(_xi1*(1.-_rho2*(omega-1.)));
  // calculate the form factor
  // the couplings in the velocity form
  Complex N1,N2,N3,N4,K1,K2,K3,K4; 
  Energy msum(m0+m1);Energy2 msum2(msum*msum);
  // in the form of the heavy quark papers
  N1 = xi;
  K1 = xi;
  N2 = 0.;
  K2 = 0.;
  N3 = 0.;
  K3 = 0.;
  N4 = 0.;
  K4 = 0.;
  // convert to our form
  f1v =-N4;
  f1a = K4;
  f2v =-N1*msum/m0;
  f2a = K1*msum/m0;
  f3v =-msum2/m0*(N2/m0+N3/m1);
  f3a = msum2/m0*(K2/m0+K3/m1);
  f4v =-msum2/m0/m0*N2;
  f4a = msum2/m0/m0*K2;
}

void LambdabExcitedLambdacSumRuleFormFactor::dataBaseOutput(ofstream & output,
							    bool header,
							    bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create /Herwig++/LambdabExcitedLambdacSumRuleFormFactor " 
	    << fullName() << " \n";}
  output << "set " << fullName() << ":Xi          " << _xi1          << " \n";
  output << "set " << fullName() << ":Rho2        " << _rho2         << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}


}

