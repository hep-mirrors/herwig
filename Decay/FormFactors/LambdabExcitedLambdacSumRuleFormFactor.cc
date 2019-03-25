// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LambdabExcitedLambdacSumRuleFormFactor class.
//

#include "LambdabExcitedLambdacSumRuleFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

LambdabExcitedLambdacSumRuleFormFactor::LambdabExcitedLambdacSumRuleFormFactor() {
  _xi1=0.29;
  _rho2=2.01;
  // modes handled by this form-factor
  // lambda_b to lambda_c1
  addFormFactor(5122,14122,2,2,1,2,5,4);
  // lambda_b to lambda_c1*
  addFormFactor(5122,4124 ,2,4,1,2,5,4);
}

void LambdabExcitedLambdacSumRuleFormFactor::
persistentOutput(PersistentOStream & os) const {
  os << _xi1 << _rho2;
}

void LambdabExcitedLambdacSumRuleFormFactor::
persistentInput(PersistentIStream & is, int) {
  is >> _xi1 >> _rho2;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<LambdabExcitedLambdacSumRuleFormFactor,BaryonFormFactor>
describeHerwigLambdabExcitedLambdacSumRuleFormFactor("Herwig::LambdabExcitedLambdacSumRuleFormFactor", "HwFormFactors.so");

void LambdabExcitedLambdacSumRuleFormFactor::Init() {

  static ClassDocumentation<LambdabExcitedLambdacSumRuleFormFactor> documentation
    ("The LambdabExcitedLambdacSumRuleFormFactor class implements the"
     " form-factors for Lambda_b to Lambda_c1(*) from hep-ph/0012114.",
     "Lambda_b to Lambda_c1(*) used the formfactors from \\cite{Huang:2000xw}.",
     "%\\cite{Huang:2000xw}\n"
     "\\bibitem{Huang:2000xw}\n"
     "  M.~Q.~Huang, J.~P.~Lee, C.~Liu and H.~S.~Song,\n"
     "  %``Leading Isgur-Wise form factor of Lambda/b to Lambda/c1 transition  using\n"
     "  %QCD sum rules,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 502}, 133 (2001)\n"
     "  [arXiv:hep-ph/0012114].\n"
     "  %%CITATION = PHLTA,B502,133;%%\n"
     );

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
SpinHalfSpinHalfFormFactor(Energy2 q2,int,int,int,Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a) {
  useMe();
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
  f2a = Complex(0.5*(m0+m1)*( g2v/m0+g3v/m1));
  f2v =-Complex(0.5*(m0+m1)*( g2a/m0+g3a/m1));
  f3a = Complex(0.5*(m0+m1)*( g2v/m0-g3v/m1));
  f3v = Complex(0.5*(m0+m1)*(-g2a/m0+g3a/m1));
}

void  LambdabExcitedLambdacSumRuleFormFactor::
 SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int,int,int,Energy m0,Energy m1,
				 Complex & f1v,Complex & f2v,
				 Complex & f3v,Complex & f4v,
				 Complex & f1a,Complex & f2a,
				 Complex & f3a,Complex & f4a )
{
  useMe();
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
  f3v =-Complex(msum2/m0*(N2/m0+N3/m1));
  f3a = Complex(msum2/m0*(K2/m0+K3/m1));
  f4v =-Complex(msum2/m0/m0*N2);
  f4a = Complex(msum2/m0/m0*K2);
}

void LambdabExcitedLambdacSumRuleFormFactor::dataBaseOutput(ofstream & output,
							    bool header,
							    bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::LambdabExcitedLambdacSumRuleFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":Xi          " << _xi1          << " \n";
  output << "newdef " << name() << ":Rho2        " << _rho2         << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
