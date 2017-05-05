// -*- C++ -*-
//
// HQETFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HQETFormFactor class.
//

#include "HQETFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HQETFormFactor::HQETFormFactor() {
  // parameters from arXiv:0705.4008
  // F(1) from F(1)V_cb and our value of V_cb
  _f1scalar = 1.0269328;
  _f1vector = 0.84;
  // slope
  _rho2scalar = 1.17;
  _rho2vector = 1.179;
  // R_1(1)
  _r1 = 1.417;
  // R_2(1)
  _r2 = 0.836;
  // allowed form factors
  addFormFactor(-521, 421  ,0,-2, 5, 4);
  addFormFactor(-521, 423  ,1,-2, 5, 4);
  addFormFactor(-511, 411  ,0, 1, 5, 4);
  addFormFactor(-511, 413  ,1, 1, 5, 4);
  // set the initial number of modes
  initialModes(numberOfFactors());
}

void HQETFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _f1scalar << _f1vector << _r1 << _r2 << _rho2scalar << _rho2vector;
}

void HQETFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _f1scalar >> _f1vector >> _r1 >> _r2 >> _rho2scalar >> _rho2vector;
}

ClassDescription<HQETFormFactor> HQETFormFactor::initHQETFormFactor;
// Definition of the static class description member.

void HQETFormFactor::Init() {

  static ClassDocumentation<HQETFormFactor> documentation
    ("The HQETFormFactor class uses the parameterisation of hep-ph/9712417"
     " of the form factor in the heavy quark limit.",
     "The parameterisation of \\cite{Caprini:1997mu} was used for the "
     "$B\\to D^{(*)}$ form factors together with the parameters from"
     "\\cite{Snyder:2007qn} for the $D$ and \\cite{Aubert:2007rs}"
     " for the $D^*$",
     "\\bibitem{Caprini:1997mu} I.~Caprini, L.~Lellouch and M.~Neubert,"
     "Nucl.\\ Phys.\\  B {\\bf 530} (1998) 153 [arXiv:hep-ph/9712417].\n"
     "%%CITATION = NUPHA,B530,153;%%\n"
     "\\bibitem{Aubert:2007rs} B.~Aubert {\\it et al.}  [BABAR Collaboration],"
     "arXiv:0705.4008 [hep-ex]. %%CITATION = ARXIV:0705.4008;%%\n"
     "\\bibitem{Snyder:2007qn} A.~E.~Snyder, [arXiv:hep-ex/0703035].\n"
     "%%CITATION = ECONF,C0610161,015;%%\n");

  static Parameter<HQETFormFactor,double> interfaceF1Scalar
    ("F1Scalar",
     "The normalisation factor for the scalar form factor",
     &HQETFormFactor::_f1scalar, 1.0269328, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceF1Vector
    ("F1Vector",
     "The normalisation factor for the vector form factor",
     &HQETFormFactor::_f1vector, 0.84,  0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceRho2Scalar
    ("Rho2Scalar",
     "The slope parameter for the scalar form factor",
     &HQETFormFactor::_rho2scalar, 1.17, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceRho2Vector
    ("Rho2Vector",
     "The slope parameter for the vector form factor",
     &HQETFormFactor::_rho2vector, 1.179, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceR1
    ("R1",
     "The ratio R_1 at omega=1",
     &HQETFormFactor::_r1, 1.417, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceR2
    ("R2",
     "The ratio R_2 at omega=1",
     &HQETFormFactor::_r2, 0.836, 0.0, 10.0,
     false, false, Interface::limited);

}

void HQETFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int ,int,int,
					    Energy m0,Energy m1,Complex & f0,
					    Complex & fp) const {
  useMe();
  double omega = 0.5*(sqr(m0)+sqr(m1)-q2)/m0/m1;
  double root = sqrt(1.+omega),rt2=sqrt(2.);
  double z = (root-rt2)/(root+rt2);
  double Rs = 2.*sqrt(m0*m1)/(m0+m1);
  fp = 1.-8.*_rho2scalar*z+((51.*_rho2scalar-10.)-(252*_rho2scalar-84.)*z)*sqr(z);
  fp *=_f1scalar/Rs;
  f0 = fp*(1.-q2/sqr(m0+m1));
}

void HQETFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int,
					    int, int, Energy m0, Energy m1,
					    Complex & A0, Complex & A1,Complex & A2,
					    Complex & V) const {
  useMe();
  double omega = 0.5*(sqr(m0)+sqr(m1)-q2)/m0/m1;
  double root = sqrt(1.+omega),rt2=sqrt(2.);
  double z = (root-rt2)/(root+rt2);
  double hA1 = _f1vector*(1.-8.*_rho2vector*z+((53.*_rho2vector-15.)
					       -(231.*_rho2vector-91.)*z)*sqr(z));
  double wmo=omega-1.;
  double R1 = _r1-0.12*wmo+0.05*sqr(wmo);
  double R2 = _r2+0.11*wmo-0.06*sqr(wmo);
  double Rs = 2.*sqrt(m0*m1)/(m0+m1);
  A1 = 0.5*(omega+1.)*Rs*hA1;
  A2 = R2*hA1/Rs;
  V  =-R1*hA1/Rs;
  Complex A3 = 0.5/m1*((m0+m1)*A1-(m0-m1)*A2);
  A0 = A3+0.5*A2*q2/m1/(m0+m1);
}

void HQETFormFactor::dataBaseOutput(ofstream & os,bool header,bool create) const {
  if(header) os << "update decayers set parameters=\"";
  if(create) os << "create Herwig::HQETFormFactor " << name() << "\n";
  ScalarFormFactor::dataBaseOutput(os,false,false);
  os << "newdef " << name() << ":F1Scalar   " << _f1scalar   << "\n";
  os << "newdef " << name() << ":F1Vector   " << _f1vector   << "\n";
  os << "newdef " << name() << ":Rho2Scalar " << _rho2scalar << "\n";
  os << "newdef " << name() << ":Rho2Vector " << _rho2vector << "\n";
  os << "newdef " << name() << ":R1         " << _r1         << "\n";
  os << "newdef " << name() << ":R2         " << _r2         << "\n";
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
