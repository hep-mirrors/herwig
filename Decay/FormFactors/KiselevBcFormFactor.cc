// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KiselevBcFormFactor class.
//

#include "KiselevBcFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KiselevBcFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {

using namespace ThePEG;

inline KiselevBcFormFactor::KiselevBcFormFactor() 
{
  // B_c to B_s
  addFormFactor(541,531,0,-5,4,3);
  _fp.push_back(  1.30    );_fm.push_back( -5.80    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  1.8*GeV);_Mfm.push_back(  1.8*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  // B_c to B_s*
  addFormFactor(541,533,1,-5,4,3);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  1.10/GeV);
  _F0A.push_back( 8.10*GeV);_FpA.push_back( 0.20/GeV);_FmA.push_back( 1.80/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  1.8*GeV);
  _MF0A.push_back( 1.8*GeV);_MFpA.push_back( 1.8*GeV);_MFmA.push_back( 1.8*GeV);
  // B_c to B
  addFormFactor(541,511,0,-5,4,1);
  addFormFactor(541,521,0,-5,4,2);
  _fp.push_back(  1.27    );_fm.push_back( -7.30    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  1.7*GeV);_Mfm.push_back(  1.7*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  _fp.push_back(  1.27    );_fm.push_back( -7.30    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  1.7*GeV);_Mfm.push_back(  1.7*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  // B_c to B*
  addFormFactor(541,513,1,-5,4,1);
  addFormFactor(541,523,1,-5,4,2);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  1.35/GeV);
  _F0A.push_back( 9.80*GeV);_FpA.push_back( 0.35/GeV);_FmA.push_back( 2.50/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  2.2*GeV);
  _MF0A.push_back( 3.2*GeV);_MFpA.push_back( 2.2*GeV);_MFmA.push_back( 3.2*GeV);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  1.35/GeV);
  _F0A.push_back( 9.80*GeV);_FpA.push_back( 0.35/GeV);_FmA.push_back( 2.50/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  2.2*GeV);
  _MF0A.push_back( 3.2*GeV);_MFpA.push_back( 2.2*GeV);_MFmA.push_back( 3.2*GeV);
  // B_c to D
  addFormFactor(541,411,0,4,-5,-1);
  addFormFactor(541,421,0,4,-5,-2);
  _fp.push_back(  0.32    );_fm.push_back( -0.34    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  5.0*GeV);_Mfm.push_back(  5.0*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  _fp.push_back(  0.32    );_fm.push_back( -0.34    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  5.0*GeV);_Mfm.push_back(  5.0*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  // B_c to D*
  addFormFactor(541,413,1,4,-5,-1);
  addFormFactor(541,423,1,4,-5,-2);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  0.20/GeV);
  _F0A.push_back( 3.60*GeV);_FpA.push_back(-0.062/GeV);_FmA.push_back( 0.10/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  6.2*GeV);
  _MF0A.push_back(-1.0*GeV);_MFpA.push_back( 6.2*GeV);_MFmA.push_back( 6.2*GeV);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  0.20/GeV);
  _F0A.push_back( 3.60*GeV);_FpA.push_back(-0.062/GeV);_FmA.push_back( 0.10/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  6.2*GeV);
  _MF0A.push_back(-1.0*GeV);_MFpA.push_back( 6.2*GeV);_MFmA.push_back( 6.2*GeV);
  // B_c to D_s
  addFormFactor(541,431,0,4,-5,-3);
  _fp.push_back(  0.45    );_fm.push_back( -0.43    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  5.0*GeV);_Mfm.push_back(  5.0*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  // B_c to D_s
  addFormFactor(541,433,1,4,-5,-3);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  0.24/GeV);
  _F0A.push_back( 4.70*GeV);_FpA.push_back(-0.077/GeV);_FmA.push_back( 0.13/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  6.2*GeV);
  _MF0A.push_back(-1.0*GeV);_MFpA.push_back( 6.2*GeV);_MFmA.push_back( 6.2*GeV);
  // B_c to eta_c
  addFormFactor(541,441,0,4,-5,-4);
  _fp.push_back(  0.66    );_fm.push_back( -0.36    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  4.5*GeV);_Mfm.push_back(  4.5*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  // B_c to J/psi
  addFormFactor(541,443,1,4,-5,-4);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  0.11/GeV);
  _F0A.push_back( 5.90*GeV);_FpA.push_back(-0.074/GeV);_FmA.push_back( 0.12/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  5.5*GeV);
  _MF0A.push_back( 5.5*GeV);_MFpA.push_back( 5.5*GeV);_MFmA.push_back( 5.5*GeV);
  // B_c to eta_c(2S)
  addFormFactor(541,100441,0,4,-5,-4);
  _fp.push_back(  0.17    );_fm.push_back( -0.16    );_FV.push_back(  0.00/GeV);
  _F0A.push_back( 0.00*GeV);_FpA.push_back( 0.00/GeV);_FmA.push_back( 0.00/GeV);
  _Mfp.push_back(  4.5*GeV);_Mfm.push_back(  4.5*GeV);_MFV.push_back(  0.0*GeV);
  _MF0A.push_back( 0.0*GeV);_MFpA.push_back( 0.0*GeV);_MFmA.push_back( 0.0*GeV);
  // B_c to J/psi(2S)
  addFormFactor(541,100443,1,4,-5,-4);
  _fp.push_back(  0.00    );_fm.push_back(  0.00    );_FV.push_back(  0.032/GeV);
  _F0A.push_back( 1.97*GeV);_FpA.push_back(-0.015/GeV);_FmA.push_back( 0.052/GeV);
  _Mfp.push_back(  0.0*GeV);_Mfm.push_back(  0.0*GeV);_MFV.push_back(  4.5*GeV);
  _MF0A.push_back( 4.5*GeV);_MFpA.push_back( 4.5*GeV);_MFmA.push_back( 4.5*GeV);
  // set the inital number of modes
  initialModes(numberOfFactors());
}

KiselevBcFormFactor::~KiselevBcFormFactor() {}

void KiselevBcFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _fp << _fm << _FV << _F0A << _FpA << _FmA << _Mfp << _Mfm << _MFV << _MF0A
     << _MFpA << _MFmA;
}

void KiselevBcFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _fp >> _fm >> _FV >> _F0A >> _FpA >> _FmA >> _Mfp >> _Mfm >> _MFV >> _MF0A
     >> _MFpA >> _MFmA;
}

ClassDescription<KiselevBcFormFactor> KiselevBcFormFactor::initKiselevBcFormFactor;
// Definition of the static class description member.

void KiselevBcFormFactor::Init() {

  static ClassDocumentation<KiselevBcFormFactor> documentation
    ("The KiselevBcFormFactor class implements the form factors from hep-ph/0211021"
     " for the decay of the B_c");

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
     &KiselevBcFormFactor::_FV, 1./GeV, -1, 0./GeV, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceF0A
    ("F0A",
     "The value of the F_0^A form factor at q^2=0",
     &KiselevBcFormFactor::_F0A, 1.*GeV, -1, 0.*GeV, -10.*GeV, 10.*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFpA
    ("FplusA",
     "The value of the F_+^A form factor at q^2=0",
     &KiselevBcFormFactor::_FpA, 1./GeV, -1, 0./GeV, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFmA
    ("FminusA",
     "The value of the F_-^A form factor at q^2=0",
     &KiselevBcFormFactor::_FmA, 1./GeV, -1, 0./GeV, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFplus
    ("MpoleFplus",
     "The pole mass for the f_+ form factor",
     &KiselevBcFormFactor::_Mfp, GeV, -1, 0.0*GeV, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFminus
    ("MpoleFminus",
     "The pole mass for the f_- form factor",
     &KiselevBcFormFactor::_Mfm, GeV, -1, 0.0*GeV, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFV
    ("MpoleFV",
     "The pole mass for the f_V form factor",
     &KiselevBcFormFactor::_MFV, GeV, -1, 0.0*GeV, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleF0A
    ("MpoleF0A",
     "The pole mass for the f_0^A form factor",
     &KiselevBcFormFactor::_MF0A, GeV, -1, 0.0*GeV, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFpA
    ("MpoleFplusA",
     "The pole mass for the f_+^A form factor",
     &KiselevBcFormFactor::_MFpA, GeV, -1, 0.0*GeV, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFmA
    ("MpoleFminusA",
     "The pole mass for the f_-^A form factor",
     &KiselevBcFormFactor::_MFmA, GeV, -1, 0.0*GeV, -2.0*GeV, 10.0*GeV,
     false, false, true);
}

void KiselevBcFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,
						 int id0,int id1,Energy m0,Energy m1,
						 Complex & f0,Complex & fp) const
{
  fp = _fp[iloc]/(1.-q2/_Mfp[iloc]/_Mfp[iloc]);
  f0 = q2/(m0+m1)/(m0-m1)*_fm[iloc]/(1.-q2/_Mfm[iloc]/_Mfm[iloc])+fp;
}

void KiselevBcFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0,
						 int id1,Energy m0, Energy m1,
						 Complex & A0,Complex & A1,Complex & A2,
						 Complex & V) const
{
  InvEnergy fv,fp,fm;
  Energy f0;
  if(_MFV[iloc]>0){fv=_FV[iloc]/(1.-q2/_MFV[iloc]/_MFV[iloc]);}
  else{fv=_FV[iloc];}
  if(_MFmA[iloc]>0){fm=_FmA[iloc]/(1.-q2/_MFmA[iloc]/_MFmA[iloc]);}
  else{fm=_FmA[iloc];}
  if(_MFpA[iloc]>0){fp=_FpA[iloc]/(1.-q2/_MFpA[iloc]/_MFpA[iloc]);}
  else{fp=_FpA[iloc];}
  if(_MF0A[iloc]>0){f0=_F0A[iloc]/(1.-q2/_MF0A[iloc]/_MF0A[iloc]);}
  else{f0=_F0A[iloc];}
  Energy msum(m0+m1);
  V  =-fv*msum;
  A1 = f0/msum;
  A2 =-fp*msum;
  A0 = 0.5/m1*(f0+msum*(m0-m1)*fp+q2*fm);
}

void KiselevBcFormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create){output << "create Herwig++::KiselevBcFormFactor " << fullName() << " \n";}
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      if(ix<initialModes())
	{
	  output << "set " << fullName() << ":Fplus "  << ix << "  " 
		 << _fp[ix]  << "\n";
	  output << "set " << fullName() << ":Fminus "  << ix << "  " 
		 << _fm[ix]  << "\n";
	  output << "set " << fullName() << ":FV "  << ix << "  " 
		 << _FV[ix]*GeV  << "\n";
	  output << "set " << fullName() << ":F0A "  << ix << "  " 
		 << _F0A[ix]/GeV  << "\n";
	  output << "set " << fullName() << ":FplusA "  << ix << "  " 
		 << _FpA[ix]*GeV  << "\n";
	  output << "set " << fullName() << ":FminusA "  << ix << "  " 
		 << _FmA[ix]*GeV  << "\n";
	  output << "set " << fullName() << ":MpoleFplus "  << ix << "  " 
		 << _Mfp[ix]/GeV  << "\n";
	  output << "set " << fullName() << ":MpoleFminus "  << ix << "  " 
		 << _Mfm[ix]/GeV  << "\n";
	  output << "set " << fullName() << ":MpoleFV "  << ix << "  " 
		 << _MFV[ix]/GeV  << "\n";
	  output << "set " << fullName() << ":MpoleF0A "  << ix << "  " 
		 << _MF0A[ix]/GeV  << "\n";
	  output << "set " << fullName() << ":MpoleFplusA "  << ix << "  " 
		 << _MFpA[ix]/GeV  << "\n";
	  output << "set " << fullName() << ":MpoleFminusA "  << ix << "  " 
		 << _MFmA[ix]/GeV  << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Fplus "  << ix << "  " 
		 << _fp[ix]  << "\n";
	  output << "insert " << fullName() << ":Fminus "  << ix << "  " 
		 << _fm[ix]  << "\n";
	  output << "insert " << fullName() << ":FV "  << ix << "  " 
		 << _FV[ix]*GeV  << "\n";
	  output << "insert " << fullName() << ":F0A "  << ix << "  " 
		 << _F0A[ix]/GeV  << "\n";
	  output << "insert " << fullName() << ":FplusA "  << ix << "  " 
		 << _FpA[ix]*GeV  << "\n";
	  output << "insert " << fullName() << ":FminusA "  << ix << "  " 
		 << _FmA[ix]*GeV  << "\n";
	  output << "insert " << fullName() << ":MpoleFplus "  << ix << "  " 
		 << _Mfp[ix]/GeV  << "\n";
	  output << "insert " << fullName() << ":MpoleFminus "  << ix << "  " 
		 << _Mfm[ix]/GeV  << "\n";
	  output << "insert " << fullName() << ":MpoleFV "  << ix << "  " 
		 << _MFV[ix]/GeV  << "\n";
	  output << "insert " << fullName() << ":MpoleF0A "  << ix << "  " 
		 << _MF0A[ix]/GeV  << "\n";
	  output << "insert " << fullName() << ":MpoleFplusA "  << ix << "  " 
		 << _MFpA[ix]/GeV  << "\n";
	  output << "insert " << fullName() << ":MpoleFminusA "  << ix << "  " 
		 << _MFmA[ix]/GeV  << "\n";
	}
    }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
