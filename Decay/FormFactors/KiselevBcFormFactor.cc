// -*- C++ -*-
//
// KiselevBcFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KiselevBcFormFactor class.
//
#include "KiselevBcFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KiselevBcFormFactor::KiselevBcFormFactor() :
  _fp  (16,0.    ), _fm  (16,0.    ), _FV( 16,ZERO),
  _F0A (16,ZERO), _FpA (16,ZERO), _FmA (16,ZERO),
  _Mfp (16,ZERO), _Mfm (16,ZERO), _MFV (16,ZERO),
  _MF0A(16,ZERO), _MFpA(16,ZERO), _MFmA(16,ZERO) {
  // B_c to B_s
  addFormFactor(541,531,0,-5,4,3);
  _fp[0] =   1.30    ;_Mfp[0] =   1.8*GeV;
  _fm[0] =  -5.80    ;_Mfm[0] =   1.8*GeV;
  // B_c to B_s*
  addFormFactor(541,533,1,-5,4,3);
  _FV[1] =   1.10/GeV;_MFV[1] =   1.8*GeV;
  _F0A[1] =  8.10*GeV;_MF0A[1] =  1.8*GeV;
  _FpA[1] =  0.20/GeV;_MFpA[1] =  1.8*GeV;
  _FmA[1] =  1.80/GeV;_MFmA[1] =  1.8*GeV;
  // B_c to B
  addFormFactor(541,511,0,-5,4,1);
  addFormFactor(541,521,0,-5,4,2);
  for(unsigned int ix=2;ix<4;++ix) {
    _fp[ix] =   1.27    ;_Mfp[ix] =   1.7*GeV;
    _fm[ix] =  -7.30    ;_Mfm[ix] =   1.7*GeV;
  }
  // B_c to B*
  addFormFactor(541,513,1,-5,4,1);
  addFormFactor(541,523,1,-5,4,2);
  for(unsigned int ix=4;ix<6;++ix) {
    _FV[ix] =   1.35/GeV; _MFV[ix] =   2.2*GeV;
    _F0A[ix] =  9.80*GeV; _MF0A[ix] =  3.2*GeV;
    _FpA[ix] =  0.35/GeV; _MFpA[ix] =  2.2*GeV;
    _FmA[ix] =  2.50/GeV; _MFmA[ix] =  3.2*GeV;
  }
  // B_c to D
  addFormFactor(-541,-411,0,-4,5,1);
  addFormFactor(-541,-421,0,-4,5,2);
  for(unsigned int ix=6;ix<8;++ix) {
    _fp[ix] =   0.32    ; _Mfp[ix] =   5.0*GeV;
    _fm[ix] =  -0.34    ; _Mfm[ix] =   5.0*GeV;
  }
  // B_c to D*
  addFormFactor(-541,-413,1,-4,5,1);
  addFormFactor(-541,-423,1,-4,5,2);
  for(unsigned int ix=8;ix<10;++ix) {
    _FV[ix] =   0.20 /GeV; _MFV[ix] =   6.2*GeV;
    _F0A[ix] =  3.60 *GeV; _MF0A[ix] = -1.0*GeV;
    _FpA[ix] = -0.062/GeV; _MFpA[ix] =  6.2*GeV;
    _FmA[ix] =  0.10 /GeV; _MFmA[ix] =  6.2*GeV;
  }
  // B_c to D_s
  addFormFactor(-541,-431,0,-4,5,3);
  _fp[10] =   0.45    ; _Mfp[10] =   5.0*GeV;
  _fm[10] =  -0.43    ; _Mfm[10] =   5.0*GeV;
  // B_c to D_s
  addFormFactor(-541,-433,1,-4,5,3);
  _FV[11] =   0.24 /GeV; _MFV[11] =   6.2*GeV;
  _F0A[11] =  4.70 *GeV; _MF0A[11] = -1.0*GeV;
  _FpA[11] = -0.077/GeV; _MFpA[11] =  6.2*GeV;
  _FmA[11] =  0.13 /GeV; _MFmA[11] =  6.2*GeV;
  // B_c to eta_c
  addFormFactor(-541,441,0,-4,5,4);
  _fp[12] =   0.66    ; _Mfp[12] =   4.5*GeV;
  _fm[12] =  -0.36    ; _Mfm[12] =   4.5*GeV;
  // B_c to J/psi
  addFormFactor(-541,443,1,-4,5,4);
  _FV[13] =   0.11 /GeV; _MFV[13] =   5.5*GeV;
  _F0A[13] =  5.90 *GeV; _MF0A[13] =  5.5*GeV;
  _FpA[13] = -0.074/GeV; _MFpA[13] =  5.5*GeV;
  _FmA[13] =  0.12 /GeV; _MFmA[13] =  5.5*GeV;
  // B_c to eta_c(2S)
  addFormFactor(-541,100441,0,-4,5,4);
  _fp[14] =   0.17    ; _Mfp[14] =   4.5*GeV;
  _fm[14] =  -0.16    ; _Mfm[14] =   4.5*GeV;
  // B_c to J/psi(2S)
  addFormFactor(-541,100443,1,-4,5,4);
  _FV[15] =   0.035/GeV; _MFV[15] =   4.2*GeV;
  _F0A[15] =  1.686*GeV; _MF0A[15] =  4.2*GeV;
  _FpA[15] = -0.015/GeV; _MFpA[15] =  4.2*GeV;
  _FmA[15] =  0.052/GeV; _MFmA[15] =  4.2*GeV;
  // set the inital number of modes
  initialModes(numberOfFactors());
}

void KiselevBcFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_fp.size() ||isize!=_fm.size()  ||isize!=_FV.size()  ||isize!=_F0A.size()||
     isize!=_FpA.size()||isize!=_FmA.size() ||isize!=_Mfp.size() ||isize!=_Mfm.size()||
     isize!=_MFV.size()||isize!=_MF0A.size()||isize!=_MFpA.size()||isize!=_MFmA.size())
    throw InitException() << "Inconsistent parameters in KiselevBcFormFactor::doinit()" 
			  << Exception::abortnow;

}

void KiselevBcFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _fp << _fm << ounit(_FV,1/GeV) << ounit(_F0A,GeV) << ounit(_FpA,1/GeV) 
     << ounit(_FmA,1/GeV) << ounit(_Mfp,GeV) << ounit(_Mfm,GeV) << ounit(_MFV,GeV) 
     << ounit(_MF0A,GeV) << ounit(_MFpA,GeV) << ounit(_MFmA,GeV);
}

void KiselevBcFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _fp >> _fm >> iunit(_FV,1/GeV) >> iunit(_F0A,GeV) >> iunit(_FpA,1/GeV) 
     >> iunit(_FmA,1/GeV) >> iunit(_Mfp,GeV) >> iunit(_Mfm,GeV) >> iunit(_MFV,GeV) 
     >> iunit(_MF0A,GeV) >> iunit(_MFpA,GeV) >> iunit(_MFmA,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KiselevBcFormFactor,ScalarFormFactor>
describeHerwigKiselevBcFormFactor("Herwig::KiselevBcFormFactor", "HwFormFactors.so");

void KiselevBcFormFactor::Init() {

  static ClassDocumentation<KiselevBcFormFactor> documentation
    ("The KiselevBcFormFactor class implements the form factors from hep-ph/0211021"
     " for the decay of the B_c",
     "The form factors of \\cite{Kiselev:2002vz} for the decay of the $B_c$ meson"
     " were used.",
     "\\bibitem{Kiselev:2002vz} V.~V.~Kiselev, arXiv:hep-ph/0211021.\n"
     "%%CITATION = HEP-PH/0211021;%%");

  static ParVector<KiselevBcFormFactor,double> interfaceFplus
    ("Fplus",
     "The value of the f_+ form factor at q^2=0",
     &KiselevBcFormFactor::_fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<KiselevBcFormFactor,double> interfaceFminus
    ("Fminus",
     "The value of the f_- form factor at q^2=0",
     &KiselevBcFormFactor::_fm, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFV
    ("FV",
     "The value of the F_V form factor at q^2=0",
     &KiselevBcFormFactor::_FV, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceF0A
    ("F0A",
     "The value of the F_0^A form factor at q^2=0",
     &KiselevBcFormFactor::_F0A, 1.*GeV, -1, ZERO, -10.*GeV, 10.*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFpA
    ("FplusA",
     "The value of the F_+^A form factor at q^2=0",
     &KiselevBcFormFactor::_FpA, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFmA
    ("FminusA",
     "The value of the F_-^A form factor at q^2=0",
     &KiselevBcFormFactor::_FmA, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFplus
    ("MpoleFplus",
     "The pole mass for the f_+ form factor",
     &KiselevBcFormFactor::_Mfp, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFminus
    ("MpoleFminus",
     "The pole mass for the f_- form factor",
     &KiselevBcFormFactor::_Mfm, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFV
    ("MpoleFV",
     "The pole mass for the f_V form factor",
     &KiselevBcFormFactor::_MFV, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleF0A
    ("MpoleF0A",
     "The pole mass for the f_0^A form factor",
     &KiselevBcFormFactor::_MF0A, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFpA
    ("MpoleFplusA",
     "The pole mass for the f_+^A form factor",
     &KiselevBcFormFactor::_MFpA, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFmA
    ("MpoleFminusA",
     "The pole mass for the f_-^A form factor",
     &KiselevBcFormFactor::_MFmA, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);
}

void KiselevBcFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,
						 int,int,Energy m0,Energy m1,
						 Complex & f0,Complex & fp) const {
  useMe();
  fp = _fp[iloc]/(1.-q2/_Mfp[iloc]/_Mfp[iloc]);
  f0 = Complex(q2/(m0+m1)/(m0-m1)*_fm[iloc]/(1.-q2/_Mfm[iloc]/_Mfm[iloc]))+fp;
}

void KiselevBcFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int,
						 int,Energy m0, Energy m1,
						 Complex & A0,Complex & A1,Complex & A2,
						 Complex & V) const {
  useMe();
  InvEnergy fv,fp,fm;
  Energy f0;
  if(_MFV[iloc]>ZERO)  fv = _FV[iloc]/(1.-q2/_MFV[iloc]/_MFV[iloc]);
  else                  fv = _FV[iloc];
  if(_MFmA[iloc]>ZERO) fm = _FmA[iloc]/(1.-q2/_MFmA[iloc]/_MFmA[iloc]);
  else                  fm = _FmA[iloc];
  if(_MFpA[iloc]>ZERO) fp = _FpA[iloc]/(1.-q2/_MFpA[iloc]/_MFpA[iloc]);
  else                  fp = _FpA[iloc];
  if(_MF0A[iloc]>ZERO) f0 = _F0A[iloc]/(1.-q2/_MF0A[iloc]/_MF0A[iloc]);
  else                  f0 = _F0A[iloc];
  Energy msum(m0+m1);
  V  = fv*msum;
  A1 =-f0/msum;
  A2 = fp*msum;
  A0 =-0.5/m1*(f0+msum*(m0-m1)*fp+q2*fm);
}

void KiselevBcFormFactor::dataBaseOutput(ofstream & output,
					 bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KiselevBcFormFactor " << name() << " \n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":Fplus "  << ix << "  " 
	     << _fp[ix]  << "\n";
      output << "newdef " << name() << ":Fminus "  << ix << "  " 
	     << _fm[ix]  << "\n";
      output << "newdef " << name() << ":FV "  << ix << "  " 
	     << _FV[ix]*GeV  << "\n";
      output << "newdef " << name() << ":F0A "  << ix << "  " 
	     << _F0A[ix]/GeV  << "\n";
      output << "newdef " << name() << ":FplusA "  << ix << "  " 
	     << _FpA[ix]*GeV  << "\n";
      output << "newdef " << name() << ":FminusA "  << ix << "  " 
	     << _FmA[ix]*GeV  << "\n";
      output << "newdef " << name() << ":MpoleFplus "  << ix << "  " 
	     << _Mfp[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFminus "  << ix << "  " 
	     << _Mfm[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFV "  << ix << "  " 
	     << _MFV[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleF0A "  << ix << "  " 
	     << _MF0A[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFplusA "  << ix << "  " 
	     << _MFpA[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFminusA "  << ix << "  " 
	     << _MFmA[ix]/GeV  << "\n";
    }
    else {
      output << "insert " << name() << ":Fplus "  << ix << "  " 
	     << _fp[ix]  << "\n";
      output << "insert " << name() << ":Fminus "  << ix << "  " 
	     << _fm[ix]  << "\n";
      output << "insert " << name() << ":FV "  << ix << "  " 
	     << _FV[ix]*GeV  << "\n";
      output << "insert " << name() << ":F0A "  << ix << "  " 
	     << _F0A[ix]/GeV  << "\n";
      output << "insert " << name() << ":FplusA "  << ix << "  " 
	     << _FpA[ix]*GeV  << "\n";
      output << "insert " << name() << ":FminusA "  << ix << "  " 
	     << _FmA[ix]*GeV  << "\n";
      output << "insert " << name() << ":MpoleFplus "  << ix << "  " 
	     << _Mfp[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFminus "  << ix << "  " 
	     << _Mfm[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFV "  << ix << "  " 
	     << _MFV[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleF0A "  << ix << "  " 
	     << _MF0A[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFplusA "  << ix << "  " 
	     << _MFpA[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFminusA "  << ix << "  " 
	     << _MFmA[ix]/GeV  << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
