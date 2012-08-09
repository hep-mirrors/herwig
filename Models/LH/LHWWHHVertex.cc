// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWHHVertex class.
//

#include "LHWWHHVertex.h"
#include "LHModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHWWHHVertex::LHWWHHVertex() : 
  couplast_(0.), q2last_(ZERO), coup_(107) {
  orderInGs(0);
  orderInGem(2);
}

IBPtr LHWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void LHWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWHHVertex,VVSSVertex>
describeHerwigLHWWHHVertex("Herwig::LHWWHHVertex", "HwLHModel.so");

void LHWWHHVertex::Init() {

  static ClassDocumentation<LHWWHHVertex> documentation
    ("The LHWWHHVertex class implements the couplings of a pair"
     " of electroweak gauge bosons and a pair of Higgs bosons in"
     " the Little Higgs model.");

}

void LHWWHHVertex::doinit() {
  // VVHH
  addToList(  24,  -24,  25,  25);
  addToList(  23,   23,  25,  25);
  addToList(  24,  -34,  25,  25);
  addToList(  34,  -24,  25,  25);
  addToList(  23,   32,  25,  25);
  addToList(  34,  -34,  25,  25);
  addToList(  33,   33,  25,  25);
  addToList(  32,   32,  25,  25);
  addToList(  23,   33,  25,  25);
  addToList(  32,   33,  25,  25);
  // VVH Phi_0
  addToList(  24,  -24,  25,  35);
  addToList(  23,   23,  25,  35);
  addToList(  24,  -34,  25,  35);
  addToList(  34,  -24,  25,  35);
  addToList(  23,   32,  25,  35);
  addToList(  34,  -34,  25,  35);
  addToList(  33,   33,  25,  35);
  addToList(  32,   32,  25,  35);
  addToList(  23,   33,  25,  35);
  addToList(  32,   33,  25,  35);
  // VV Phi_0 Phi_0
  addToList(  24,  -24,  35,  35);
  addToList(  23,   23,  35,  35);
  addToList(  24,  -34,  35,  35);
  addToList(  34,  -24,  35,  35);
  addToList(  23,   32,  35,  35);
  addToList(  34,  -34,  35,  35);
  addToList(  33,   33,  35,  35);
  addToList(  32,   32,  35,  35);
  addToList(  23,   33,  35,  35);
  addToList(  32,   33,  35,  35);
  // VV Phi_P Phi_P
  addToList(  24,  -24,  36,  36);
  addToList(  23,   23,  36,  36);
  addToList(  24,  -34,  36,  36);
  addToList(  34,  -24,  36,  36);
  addToList(  23,   32,  36,  36);
  addToList(  34,  -34,  36,  36);
  addToList(  33,   33,  36,  36);
  addToList(  32,   32,  36,  36);
  addToList(  23,   33,  36,  36);
  addToList(  32,   33,  36,  36);
  // VV Phi+ Phi-
  addToList(  24,  -24,  37, -37);
  addToList(  23,   23,  37, -37);
  addToList(  22,   22,  37, -37);
  addToList(  22,   23,  37, -37);
  addToList(  24,  -34,  37, -37);
  addToList(  34,  -24,  37, -37);
  addToList(  34,  -34,  37, -37);
  addToList(  33,   33,  37, -37);
  addToList(  32,   32,  37, -37);
  addToList(  32,   33,  37, -37);
  addToList(  22,   32,  37, -37);
  addToList(  23,   33,  37, -37);
  addToList(  23,   32,  37, -37);
  // VV Phi++ Phi--
  addToList(  24,  -24,  38, -38);
  addToList(  23,   23,  38, -38);
  addToList(  22,   22,  38, -38);
  addToList(  22,   23,  38, -38);
  addToList(  24,  -34,  38, -38);
  addToList(  34,  -24,  38, -38);
  addToList(  34,  -34,  38, -38);
  addToList(  33,   33,  38, -38);
  addToList(  32,   32,  38, -38);
  addToList(  32,   33,  38, -38);
  addToList(  22,   32,  38, -38);
  addToList(  23,   33,  38, -38);
  addToList(  23,   32,  38, -38);
  // VV H phi-  + cc
  addToList(  24,   22,  25, -37);
  addToList(  24,   23,  25, -37);
  addToList(  24,   32,  25, -37);
  addToList(  24,   33,  25, -37);
  addToList(  34,   22,  25, -37);
  addToList(  34,   23,  25, -37);
  addToList(  34,   32,  25, -37);
  addToList(  34,   33,  25, -37);
  addToList( -24,   22,  25,  37);
  addToList( -24,   23,  25,  37);
  addToList( -24,   32,  25,  37);
  addToList( -24,   33,  25,  37);
  addToList( -34,   22,  25,  37);
  addToList( -34,   23,  25,  37);
  addToList( -34,   32,  25,  37);
  addToList( -34,   33,  25,  37);
  // VV phi0  phi-  + cc
  addToList(  24,   22,  35, -37);
  addToList(  24,   23,  35, -37);
  addToList(  24,   32,  35, -37);
  addToList(  24,   33,  35, -37);
  addToList(  34,   22,  35, -37);
  addToList(  34,   23,  35, -37);
  addToList(  34,   32,  35, -37);
  addToList(  34,   33,  35, -37);
  addToList( -24,   22,  35,  37);
  addToList( -24,   23,  35,  37);
  addToList( -24,   32,  35,  37);
  addToList( -24,   33,  35,  37);
  addToList( -34,   22,  35,  37);
  addToList( -34,   23,  35,  37);
  addToList( -34,   32,  35,  37);
  addToList( -34,   33,  35,  37);
  // VV phiP  phi-  + cc
  addToList(  24,   22,  36, -37);
  addToList(  24,   23,  36, -37);
  addToList(  24,   32,  36, -37);
  addToList(  24,   33,  36, -37);
  addToList(  34,   22,  36, -37);
  addToList(  34,   23,  36, -37);
  addToList(  34,   32,  36, -37);
  addToList(  34,   33,  36, -37);
  addToList( -24,   22,  36,  37);
  addToList( -24,   23,  36,  37);
  addToList( -24,   32,  36,  37);
  addToList( -24,   33,  36,  37);
  addToList( -34,   22,  36,  37);
  addToList( -34,   23,  36,  37);
  addToList( -34,   32,  36,  37);
  addToList( -34,   33,  36,  37);
  // VV phi+ phi -- + cc
  addToList(  24,   22,  37, -38);
  addToList(  24,   23,  37, -38);
  addToList(  24,   32,  37, -38);
  addToList(  24,   33,  37, -38);
  addToList(  34,   22,  37, -38);
  addToList(  34,   23,  37, -38);
  addToList(  34,   32,  37, -38);
  addToList(  34,   33,  37, -38);
  addToList( -24,   22, -37,  38);
  addToList( -24,   23, -37,  38);
  addToList( -24,   32, -37,  38);
  addToList( -24,   33, -37,  38);
  addToList( -34,   22, -37,  38);
  addToList( -34,   23, -37,  38);
  addToList( -34,   32, -37,  38);
  addToList( -34,   33, -37,  38);
  // VV H phi-- + cc
  addToList(  24,   24,  25, -38);
  addToList( -24,  -24,  25,  38);
  addToList(  24,   34,  25, -38);
  addToList( -24,  -34,  25,  38);
  addToList(  34,   34,  25, -38);
  addToList( -34,  -34,  25,  38);
  // VV phi0 phi-- + cc
  addToList(  24,   24,  35, -38);
  addToList( -24,  -24,  35,  38);
  addToList(  24,   34,  35, -38);
  addToList( -24,  -34,  35,  38);
  addToList(  34,   34,  35, -38);
  addToList( -34,  -34,  35,  38);
  // VV phiP  phi-- + cc
  addToList(  24,   24,  36, -38);
  addToList( -24,  -24,  36,  38);
  addToList(  24,   34,  36, -38);
  addToList( -24,  -34,  36,  38);
  addToList(  34,   34,  36, -38);
  addToList( -34,  -34,  36,  38);
  VVSSVertex::doinit();
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWWVertex::doinit()"
			  << Exception::runerror;
  double sw2(sin2ThetaW()),cw2(1.-sw2);
  double sw(sqrt(sw2)),cw(sqrt(cw2));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double s0   (model->sinTheta0());
  double sPlus(model->sinThetaPlus());
  // VV HH
  coup_[  0] = 0.5/sw2;
  coup_[  1] = 0.5/sw2/cw2;
  coup_[  2] = 0.;
  coup_[  3] =-0.25/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[  4] =-0.25/sw/cw2*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[  5] =-0.5/sw2;
  coup_[  6] =-0.5/sw2;
  coup_[  7] =-0.5/cw2;
  coup_[  8] =-0.25/sw2/cw*(sqr(c)-sqr(s))/c/s;
  coup_[  9] =-0.25/sw/cw*(sqr(c*sp)+sqr(s*cp))/c/s/sp/cp;
  // VV H Phi_0
  coup_[ 10] = 0.5*s0/sw2;
  coup_[ 11] = 1.5*s0/sw2/cw2;
  coup_[ 12] = 0.;
  coup_[ 13] =-0.25/sw2*(sqr(c)-sqr(s))/c/s*s0;
  coup_[ 14] =-0.75/sw/cw2*(sqr(cp)-sqr(sp))/cp/sp*s0;
  coup_[ 15] = -0.5/sw2*s0;
  coup_[ 16] = 0.5/sw2*(1.+sqr(sqr(c)-sqr(s))/sqr(s*c))*s0;
  coup_[ 17] = 0.5/cw2*(1.+sqr(sqr(cp)-sqr(sp))/sqr(sp*cp))*s0;
  coup_[ 18] =-0.75/cw/sw2*(sqr(c)-sqr(s))/c/s*s0;
  coup_[ 19] = 0.25/sw/cw/c/s/sp/cp*((sqr(c*sp)+sqr(s*cp))
				    +2.*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp)))*s0;
  // VV phi0 phi0
  coup_[ 20] = 1./sw2;
  coup_[ 21] = 2./cw2/sw2;
  coup_[ 22] = 0.;
  coup_[ 23] =-0.5/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 24] =-1./sw/cw2*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 25] =-1./sw2;
  coup_[ 26] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 26] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 28] =-1./cw/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 29] = 0.5/cw/sw*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))/s/c/sp/cp;
  // VV phi_P phi_P
  coup_[ 30] = 1./sw2;
  coup_[ 31] = 2./sw2/cw2;
  coup_[ 32] = 0.;
  coup_[ 33] =-0.5/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 34] =-1./sw/cw2*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 35] =-1./sw2;
  coup_[ 36] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 37] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 38] =-1./cw/sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 39] = 0.5/cw/sw*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))/s/c/sp/cp;
  // VV phi+ phi-
  coup_[ 40] = 2./sw2;
  coup_[ 41] = 2.*sw2/cw2;
  coup_[ 42] = 2.;
  coup_[ 43] =-2.*sw/cw;
  coup_[ 44] =-1./sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 45] = 0.;
  coup_[ 46] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 47] =-0.5/sw2/sqr(s*c);
  coup_[ 48] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 49] = 0.;
  coup_[ 50] =-1./cw*(sqr(cp)-sqr(sp))/sp/cp;
  coup_[ 51] = 0.;
  coup_[ 52] = sw/cw2*(sqr(cp)-sqr(sp))/sp/cp;
  // VV phi++ phi--
  coup_[ 53] = 1./sw2;
  coup_[ 54] = 2./cw2/sw2*sqr(1.-2.*sw2);
  coup_[ 55] = 8.;
  coup_[ 56] = 4./sw/cw*(1.-2.*sw2);
  coup_[ 57] =-0.5/sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 58] = 2./sw*(sqr(c)-sqr(s))/s/c;
  coup_[ 59] =-1./sw2;
  coup_[ 60] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 61] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 62] =-0.5/cw/sw*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))/s/c/sp/cp;
  coup_[ 63] =-2./cw*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 64] = 1./sw2/cw*(sqr(c )-sqr(s ))/s /c *(1.-2.*sw2);
  coup_[ 65] =-1./cw2/sw*(sqr(cp)-sqr(sp))/sp/cp*(1.-2.*sw2);
  // VV h phi-
  coup_[ 66] =-0.5/sw*(sPlus-sqrt(2.)*s0);
  coup_[ 67] = 0.5/cw/sw2*(sPlus*sw2-sqrt(2.)*s0*(1.+sw2));
  coup_[ 68] =-0.25/sw/cw*(sqr(cp)-sqr(sp))/cp/sp*(sPlus-2.*sqrt(2.)*s0);
  coup_[ 69] = 0.25/sw2*(sqr(c)-sqr(s))/s/c*s0;
  coup_[ 70] = 0.25/sw*(sqr(c)-sqr(s))/s/c*(sPlus-sqrt(2.)*s0);
  coup_[ 71] =-0.25/sw2/cw*(sqr(c)-sqr(s))/s/c*(sPlus*sw2-sqrt(2.)*s0*(1.+sw2));
  coup_[ 72] =-0.25/sw/cw/s/c/sp/cp*(sPlus*(sqr(c*sp)+sqr(s*cp))
				    +sqrt(2.)*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp)));
  coup_[ 73] =-0.25/sw2*s0*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phi0 phi-
  coup_[ 74] =-sqrt(0.5)/sw;
  coup_[ 75] =-sqrt(0.5)/sw2/cw*(1.+sw2);
  coup_[ 76] = sqrt(0.5)/sw/cw*(sqr(cp)-sqr(sp))/sp/cp;
  coup_[ 77] = 0.5*sqrt(0.5)/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 78] = 0.5*sqrt(0.5)/sw*(sqr(c)-sqr(s))/c/s;
  coup_[ 79] = 0.5*sqrt(0.5)/sw2/cw*(sqr(c)-sqr(s))/c/s*(1.+sw2);
  coup_[ 80] =-0.5*sqrt(0.5)/sw/cw*(sqr(cp)-sqr(sp))/cp/sp*(sqr(c)-sqr(s))/c/s;
  coup_[ 81] =-0.5*sqrt(0.5)/sw2*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phi+ phi--
  coup_[ 82] = 3./sw;
  coup_[ 83] = (1.-3.*sw2)/cw/sw2;
  coup_[ 84] = 1./sw/cw*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 85] = Complex(0.,1.)*0.5*sqrt(0.5)/sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 86] =-3./sw*(sqr(c)-sqr(s))/s/c;
  coup_[ 87] =-0.5/sw2/cw*(sqr(c)-sqr(s))/s/c*(1.-3.*sw2);
  coup_[ 88] =-0.5/sw/cw*(sqr(c)-sqr(s))/s/c*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 89] =-Complex(0.,1.)*0.5*sqrt(0.5)/sw2*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phip phi-
  coup_[ 90] =-Complex(0.,1.)/sw*sqrt(0.5);
  coup_[ 91] =-Complex(0.,1.)/sw2/cw*(1.+sw2);
  coup_[ 92] = Complex(0.,1.)/sw/cw*sqrt(0.5)*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 93] = Complex(0.,1.)/sw2*sqrt(0.5)*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[ 94] =-Complex(0.,1.)*sqrt(0.5)*0.5/sw*(sqr(c)-sqr(s))/s/c;
  coup_[ 95] = Complex(0.,1.)*sqrt(0.5)*0.5/sw2/cw*(sqr(c)-sqr(s))/s/c*(1.+sw2);
  coup_[ 96] =-Complex(0.,1.)*sqrt(0.5)*0.5/sw/cw*(sqr(c )-sqr(s ))/s/c*
                                                 (sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 97] =-Complex(0.,1.)/sw2*sqrt(0.5)*0.5*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV H phi--
  coup_[ 98] = sqrt(2.)/sw2*s0;
  coup_[ 99] =-sqrt(2.)/sw2*0.5*(sqr(c)-sqr(s))/s/c*s0;
  coup_[100] = sqrt(2.)/sw2*0.5*(pow(c,4)+pow(s,4))/sqr(s*c)*s0;
  // VV phi0 phi--
  coup_[101] = sqrt(2.)/sw2;
  coup_[102] =-sqrt(2.)/sw2*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[103] = sqrt(2.)/sw2*0.5*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phip phi--
  coup_[104] = Complex(0.,1.)*sqrt(2.)/sw2;
  coup_[105] =-Complex(0.,1.)*sqrt(2.)/sw2*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[106] = Complex(0.,1.)*sqrt(2.)/sw2*0.5*(pow(c,4)+pow(s,4))/sqr(s*c);
}

void LHWWHHVertex::setCoupling(Energy2 q2,
			       tcPDPtr part1,tcPDPtr part2,
			       tcPDPtr part3,tcPDPtr part4) {
  if( q2 != q2last_ || couplast_==0.) {
    q2last_ = q2;
    couplast_ = sqr(electroMagneticCoupling(q2));
  }
  int ibos1 = part1->id();
  int ibos2 = part2->id();
  int isca1 = part3->id();
  int isca2 = part4->id();
  if( isca1 == isca2 || 
      (isca1==25&&isca2==35) || (isca1==35&&isca2==25)) {
    unsigned int ioff = 0;
    if     (isca1!=isca2) ioff = 10;
    else if(isca1==35   ) ioff = 20;
    else if(isca1==36   ) ioff = 30;
    if(ibos1==23&&ibos2==23)
      norm(coup_[1+ioff]*couplast_);
    else if(ibos1==33&&ibos2==33)
      norm(coup_[6+ioff]*couplast_);
    else if(ibos1==33&&ibos2==33)
      norm(coup_[7+ioff]*couplast_);
    else if(abs(ibos1)==24&&abs(ibos2)==24)
      norm(coup_[0+ioff]*couplast_);
    else if(abs(ibos1)==34&&abs(ibos2)==34)
      norm(coup_[5+ioff]*couplast_);
    else if(( abs(ibos1) == 24 && abs(ibos2) == 34) ||
	    ( abs(ibos1) == 34 && abs(ibos2) == 24))
      norm(coup_[3+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 32) ||
	    ( ibos1 == 32 && ibos2 == 23))
      norm(coup_[4+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 23))
      norm(coup_[8+ioff]*couplast_);
    else if(( ibos1 == 32 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 32))
      norm(coup_[9+ioff]*couplast_);
    else
      assert(false);
  }
  else if(isca1==-isca2) {
    unsigned int ioff = abs(isca1) == 37 ? 40 : 53;
    if(abs(ibos1)==24&&abs(ibos2)==24)
      norm(coup_[0+ioff]*couplast_);
    else if(ibos1==23&&ibos2==23)
      norm(coup_[1+ioff]*couplast_);
    else if(ibos1==22&&ibos2==22)
      norm(coup_[2+ioff]*couplast_);
    else if(( ibos1 == 22 && ibos2 == 23) ||
	    ( ibos1 == 23 && ibos2 == 22))
      norm(coup_[3+ioff]*couplast_);
    else if(( abs(ibos1) == 24 && abs(ibos2) == 34) ||
	    ( abs(ibos1) == 34 && abs(ibos2) == 24))
      norm(coup_[4+ioff]*couplast_);
    else if(( ibos1 == 22 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 22))
      norm(coup_[5+ioff]*couplast_);
    else if(abs(ibos1)==34&&abs(ibos2)==34)
      norm(coup_[6+ioff]*couplast_);
    else if(ibos1==33&&ibos2==33)
      norm(coup_[7+ioff]*couplast_);
    else if(ibos1==32&&ibos2==32)
      norm(coup_[8+ioff]*couplast_);
    else if(( ibos1 == 32 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 32))
      norm(coup_[9+ioff]*couplast_);
    else if(( ibos1 == 22 && ibos2 == 32) ||
	    ( ibos1 == 32 && ibos2 == 22))
      norm(coup_[10+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 23))
      norm(coup_[11+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 32) ||
	    ( ibos1 == 32 && ibos2 == 23))
      norm(coup_[12+ioff]*couplast_);
    else
      assert(false);
  }
  else if(((abs(ibos1) == 24 || abs(ibos1) == 34) && 
	   (abs(ibos2) != 24 && abs(ibos2) != 34)) ||
	  ((abs(ibos2) == 24 || abs(ibos2) == 34) && 
	   (abs(ibos1) != 24 && abs(ibos1) != 34))) {
    int iw,ineut;
    if(abs(ibos1) == 24 || abs(ibos1) == 34) {
      iw    = abs(ibos1);
      ineut = ibos2;
    }
    else {
      iw    = abs(ibos2);
      ineut = ibos1;
    }
    unsigned int ioff = 66;
    if((isca1 == 35 && abs(isca2) == 37) ||
       (isca2 == 35 && abs(isca1) == 37)) {
      ioff += 8;
    }
    else if ((abs(isca1) == 37 && abs(isca2) == 38) ||
	     (abs(isca2) == 37 && abs(isca1) == 38)) {
      ioff += 16;
    } 
    else if ((isca1 == 35 && abs(isca2) == 37) ||
	     (isca2 == 35 && abs(isca1) == 37)) {
      ioff += 24;
    }
    else
      assert(false);
    if(iw==34) ioff += 4;
    if(ineut==22)
      norm(coup_[0+ioff]*couplast_);
    else if(ineut==23)
      norm(coup_[1+ioff]*couplast_);
    else if(ineut==32)
      norm(coup_[2+ioff]*couplast_);
    else if(ineut==33)
      norm(coup_[3+ioff]*couplast_);
    else
      assert(false);
  }
  else {
    unsigned int ioff = 98;
    if(isca1==25||isca2==25)
      ioff += 0;
    else if(isca1==35||isca2==35)
      ioff += 3;
    else if(isca1==36||isca2==36)
      ioff += 6;
    else
      assert(false);
    if(ibos1==ibos2) {
      if(abs(ibos1)==24) {
	norm(coup_[0+ioff]*couplast_);
      }
      else if(abs(ibos1)==34) {
	norm(coup_[1+ioff]*couplast_);
      }
      else
	assert(false);
    }
    else {
      norm(coup_[2+ioff]*couplast_);
    }
  }
}
