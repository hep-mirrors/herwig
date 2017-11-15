// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LightBaryonQuarkModelFormFactor class.
//

#include "LightBaryonQuarkModelFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LightBaryonQuarkModelFormFactor::doinit() {
  BaryonFormFactor::doinit();
  // check that the parameters are consistent
  unsigned int isize=numberOfFactors();
  if(isize!=_f1.size()||isize!=_f2.size()||isize!=_g1.size()||isize!=_g2.size()||
     isize!=_Lambdaf1.size()||isize!=_Lambdaf2.size()||
     isize!=_Lambdag1.size()||isize!=_Lambdag2.size())
    throw InitException() << "Inconsistent parameters in "
			  << "LightBaryonQuarkModelFormFactor::doinit()" 
			  << Exception::abortnow;
}

LightBaryonQuarkModelFormFactor::LightBaryonQuarkModelFormFactor() {
  // the various decay modes handled by the model and the parameters
  // neutron to proton
  addFormFactor(2112,2212,2,2,2,1,1,2);
  _f1.push_back(1.00);_f2.push_back(1.81/GeV);
  _g1.push_back(1.25);_g2.push_back(ZERO);
  _Lambdaf1.push_back(0.69*GeV);_Lambdaf2.push_back(0.96*GeV);
  _Lambdag1.push_back(0.76*GeV);_Lambdag2.push_back(1.04*GeV);
  // sigma+  to lambda
  addFormFactor(3222,3122,2,2,3,2,2,1);
  _f1.push_back(0.00);_f2.push_back(1.04/GeV);
  _g1.push_back(0.60);_g2.push_back(ZERO);
  _Lambdaf1.push_back(-0.32*GeV);_Lambdaf2.push_back(-1.72*GeV);
  _Lambdag1.push_back( 0.77*GeV);_Lambdag2.push_back( 1.05*GeV);
  // sigma-  to lambda
  addFormFactor(3112,3122,2,2,3,1,1,2);
  _f1.push_back(0.00);_f2.push_back(1.04/GeV);
  _g1.push_back(0.60);_g2.push_back(ZERO);
  _Lambdaf1.push_back(-0.32*GeV);_Lambdaf2.push_back(-1.72*GeV);
  _Lambdag1.push_back( 0.77*GeV);_Lambdag2.push_back( 1.05*GeV);
  // sigma-  to sigma0
  addFormFactor(3112,3212,2,2,3,1,1,2);
  _f1.push_back(1.41);_f2.push_back(0.76/GeV);
  _g1.push_back(0.69);_g2.push_back(ZERO);
  _Lambdaf1.push_back(0.60*GeV);_Lambdaf2.push_back(0.81*GeV);
  _Lambdag1.push_back(0.77*GeV);_Lambdag2.push_back(1.04*GeV);
  // sigma0  to sigma+
  addFormFactor(3212,3222,2,2,3,2,1,2);
  _f1.push_back(-1.41);_f2.push_back(-0.76/GeV);
  _g1.push_back(-0.69);_g2.push_back( ZERO);
  _Lambdaf1.push_back(0.60*GeV);_Lambdaf2.push_back(0.81*GeV);
  _Lambdag1.push_back(0.77*GeV);_Lambdag2.push_back(1.04*GeV);
  // Xi- to Xi0
  addFormFactor(3312,3322,2,2,3,3,1,2);
  _f1.push_back(-1.00);_f2.push_back(0.73/GeV);
  _g1.push_back( 0.24);_g2.push_back(ZERO);
  _Lambdaf1.push_back(0.56*GeV);_Lambdaf2.push_back(0.71*GeV);
  _Lambdag1.push_back(0.76*GeV);_Lambdag2.push_back(1.04*GeV);
  // lambda to proton
  addFormFactor(3122,2212,2,2,1,2,3,2);
  _f1.push_back(-1.19);_f2.push_back(-0.850/GeV);
  _g1.push_back(-0.99);_g2.push_back(-0.025/GeV);
  _Lambdaf1.push_back(0.71*GeV);_Lambdaf2.push_back(0.98*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.12*GeV);
  // sigma0 to proton
  addFormFactor(3212,2212,2,2,2,1,3,2);
  _f1.push_back(-0.69);_f2.push_back(0.44/GeV);
  _g1.push_back( 0.19);_g2.push_back(0.0043/GeV);
  _Lambdaf1.push_back(0.64*GeV);_Lambdaf2.push_back(0.84*GeV);
  _Lambdag1.push_back(0.83*GeV);_Lambdag2.push_back(1.16*GeV);
  // sigma- to neutron
  addFormFactor(3112,2112,2,2,1,1,3,2);
  _f1.push_back(-0.97);_f2.push_back(0.62/GeV);
  _g1.push_back(0.27);_g2.push_back(0.0061/GeV);
  _Lambdaf1.push_back(0.64*GeV);_Lambdaf2.push_back(0.90*GeV);
  _Lambdag1.push_back(0.83*GeV);_Lambdag2.push_back(1.16*GeV);
  // xi- to lambda
  addFormFactor(3312,3122,2,2,3,1,3,2);
  _f1.push_back(1.19);_f2.push_back(0.07/GeV);
  _g1.push_back(0.33);_g2.push_back(0.0076/GeV);
  _Lambdaf1.push_back(0.68*GeV);_Lambdaf2.push_back(0.89*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.10*GeV);
  // xi- to sigma0
  addFormFactor(3312,3212,2,2,3,1,3,2);
  _f1.push_back(0.69);_f2.push_back(0.98/GeV);
  _g1.push_back(0.94);_g2.push_back(0.022/GeV);
  _Lambdaf1.push_back(0.75*GeV);_Lambdaf2.push_back(1.05*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.12*GeV);
  // xi0 to sigma+
  addFormFactor(3322,3222,2,2,3,2,3,2);
  _f1.push_back(0.98);_f2.push_back(1.38/GeV);
  _g1.push_back(1.33);_g2.push_back(0.031/GeV);
  _Lambdaf1.push_back(0.75*GeV);_Lambdaf2.push_back(1.05*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.12*GeV);
  // set the inital number of form factors
  initialModes(numberOfFactors());
}

void LightBaryonQuarkModelFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _f1 << ounit(_f2,1/GeV) << _g1 << ounit(_g2,1/GeV) 
     << ounit(_Lambdaf1,GeV) << ounit(_Lambdaf2,GeV) 
     << ounit(_Lambdag1,GeV) << ounit(_Lambdag2,GeV);
}

void LightBaryonQuarkModelFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _f1 >> iunit(_f2,1/GeV) >> _g1 >> iunit(_g2,1/GeV) 
     >> iunit(_Lambdaf1,GeV) >> iunit(_Lambdaf2,GeV) 
     >> iunit(_Lambdag1,GeV) >> iunit(_Lambdag2,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<LightBaryonQuarkModelFormFactor,BaryonFormFactor>
describeHerwigLightBaryonQuarkModelFormFactor("Herwig::LightBaryonQuarkModelFormFactor", "HwFormFactors.so");

void LightBaryonQuarkModelFormFactor::Init() {

  static ClassDocumentation<LightBaryonQuarkModelFormFactor> documentation
    ("The LightBaryonQuarkModelFormFactor class implements"
     " the quark model calculation of hep-ph/9409272 for the form-factors"
     " for the light quarks",
     "The quark model calculation of \\cite{Schlumpf:1994fb} was used"
     "for the weak decay of the light baryons",
     "\\bibitem{Schlumpf:1994fb}\n"
     "F.~Schlumpf,\n"
     "Phys.\\ Rev.\\  D {\\bf 51} (1995) 2262 [arXiv:hep-ph/9409272].\n"
     "%%CITATION = PHRVA,D51,2262;%%\n");

  static ParVector<LightBaryonQuarkModelFormFactor,double> interfacef1
    ("f1",
     "The form-factor f1 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_f1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,double> interfaceg1
    ("g1",
     "The form-factor g1 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_g1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,InvEnergy> interfacef2
    ("f2",
     "The form-factor f2 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_f2,
     1./GeV, 0, ZERO, -10./GeV, 10./GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,InvEnergy> interfaceg2
    ("g2",
     "The form-factor g2 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_g2,
     1./GeV, 0, ZERO, -10./GeV, 10./GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdaf1
    ("Lambdaf1",
     "The first mass for the energy dependence of the f1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdaf1,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdaf2
    ("Lambdaf2",
     "The second mass for the energy dependence of the f1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdaf2,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdag1
    ("Lambdag1",
     "The first mass for the energy dependence of the g1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdag1,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdag2
    ("Lambdag2",
     "The second mass for the energy dependence of the g1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdag2,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);
}

// form factor for spin-1/2 to spin-1/2
void LightBaryonQuarkModelFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int mode,int, int, Energy m0, Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a) {
  useMe();
  // f_3 is zero
  f3v = 0.;
  f3a = 0.;
  // energy dependence of the f1 and g1 factors
  InvEnergy2 lam1,lam2;
  if(_Lambdaf1[mode]>ZERO) {
    lam1 = 1./sqr(_Lambdaf1[mode]);
    lam2 = 1./sqr(_Lambdaf2[mode]);
    f1v= _f1[mode]/(1.-q2*(lam1-q2*sqr(lam2)));
  }
  else {
    f1v = _f1[mode]*(1.+_Lambdaf1[mode]/GeV*q2/GeV2+_Lambdaf2[mode]/GeV*sqr(q2/GeV2));
  }
  lam1 = 1./sqr(_Lambdag1[mode]);
  lam2 = 1./sqr(_Lambdag2[mode]);
  f1a=-_g1[mode]/(1.-q2*(lam1-q2*sqr(lam2)));
  // the f2 and g2 factors
  f2v =-(m0+m1)*_f2[mode];
  f2a = (m0+m1)*_g2[mode]; 
}

void LightBaryonQuarkModelFormFactor::dataBaseOutput(ofstream& output,bool header,
						     bool create) const  {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::LightBaryonQuarkModelFormFactor " 
		    << name() << " \n";
  for(unsigned int ix=0;ix<_f1.size();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":f1 " << ix << "  " << _f1[ix] << "\n";
      output << "newdef " << name() << ":g1 " << ix << "  " << _g1[ix] << "\n";
      output << "newdef " << name() << ":f2 " << ix << "  " << _f2[ix]*GeV << "\n";
      output << "newdef " << name() << ":g2 " << ix << "  " << _g2[ix]*GeV << "\n";
      output << "newdef " << name() << ":Lambdaf1 " << ix << "  " 
	     << _Lambdaf1[ix]/GeV << "\n";
      output << "newdef " << name() << ":Lambdaf2 " << ix << "  " 
	     << _Lambdaf2[ix]/GeV << "\n";
      output << "newdef " << name() << ":Lambdag1 " << ix << "  " 
	     << _Lambdag1[ix]/GeV << "\n";
      output << "newdef " << name() << ":Lambdag2 " << ix << "  " 
	     << _Lambdag2[ix]/GeV << "\n";
    }
    else {
      output << "insert " << name() << ":f1 " << ix << "  " << _f1[ix] << "\n";
      output << "insert " << name() << ":g1 " << ix << "  " << _g1[ix] << "\n";
      output << "insert " << name() << ":f2 " << ix << "  " << _f2[ix]*GeV << "\n";
      output << "insert " << name() << ":g2 " << ix << "  " << _g2[ix]*GeV << "\n";
      output << "insert " << name() << ":Lambdaf1 " << ix << "  " 
	     << _Lambdaf1[ix]/GeV << "\n";
      output << "insert " << name() << ":Lambdaf2 " << ix << "  " 
	     << _Lambdaf2[ix]/GeV << "\n";
      output << "insert " << name() << ":Lambdag1 " << ix << "  " 
	     << _Lambdag1[ix]/GeV << "\n";
      output << "insert " << name() << ":Lambdag2 " << ix << "  " 
	     << _Lambdag2[ix]/GeV << "\n";
    }
  }
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
