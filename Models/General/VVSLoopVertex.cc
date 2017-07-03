// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSLoopVertex class.
//

#include "VVSLoopVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;
using namespace ThePEG;
namespace LT = Looptools;

IBPtr VVSLoopVertex::clone() const {
  return new_ptr(*this);
}

IBPtr VVSLoopVertex::fullclone() const {
  return new_ptr(*this);
}

void VVSLoopVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(masses,GeV) << type << couplings << Npart_;
}

void VVSLoopVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(masses,GeV) >> type >> couplings >> Npart_;
}
void VVSLoopVertex::dofinish() {
  if(loopToolsInit_) Looptools::ltexi();
  GeneralVVSVertex::dofinish();
}

ClassDescription<VVSLoopVertex> VVSLoopVertex::initVVSLoopVertex;
// Definition of the static class description member.

void VVSLoopVertex::Init() {

  static ClassDocumentation<VVSLoopVertex> documentation
    ("The VVSLoopVertex class calculates the tensor integral"
     " coefficients using Looptools.");

}

void VVSLoopVertex::setCoupling(Energy2, tcPDPtr p1, tcPDPtr,tcPDPtr) {
  if(!loopToolsInit_) {
    Looptools::ltini();
    loopToolsInit_ = true;
  }
  //Kinematic invariants
  double ps2 = invariant(0,0) / MeV2;
  double pv1s = invariant(1,1) / MeV2;
  double pv2s = invariant(2,2) / MeV2;
  Complex a(0.),b(0.),c(0.),d(0.),e(0.),f(0.);
  for(unsigned int i = 0; i< Npart_;++i) {
    double lmass = masses[i] / MeV;
    double mls = sqr(lmass);
    // left coupling for fermions or the total coupling
    // for scalars or vectors
    Complex lc = couplings[i].first;
    // double p1p2 = 0.5*(ps2-pv1s-pv2s);
    Complex C0A = LT::C0i(LT::cc0,pv2s,ps2,pv1s,mls,mls,mls);
    long theC = LT::Cget(pv2s,ps2,pv1s,    mls,mls,mls);
    // Complex C1A  = LT::Cval(LT::cc1,theC);
    // Complex C2A  = LT::Cval(LT::cc2,theC);
    // Complex C11A = LT::Cval(LT::cc11,theC);
    Complex C12A = LT::Cval(LT::cc12,theC);
    // Complex C22A = LT::Cval(LT::cc22,theC);
    // Complex C00A = LT::Cval(LT::cc00,theC);
    Complex C0B = LT::C0i(LT::cc0,pv1s,ps2,pv2s,mls,mls,mls);
    theC = LT::Cget(pv1s,ps2,pv2s,    mls,mls,mls);
    // Complex C1B  = LT::Cval(LT::cc1,theC);
    // Complex C2B  = LT::Cval(LT::cc2,theC);
    // Complex C11B = LT::Cval(LT::cc11,theC);
    Complex C12B = LT::Cval(LT::cc12,theC);
    // Complex C22B = LT::Cval(LT::cc22,theC);
    // Complex C00B = LT::Cval(LT::cc00,theC);
    if(type[i] == PDT::Spin1Half) {
      Complex lpr = lc + couplings[i].second;
      Complex loop = 2.*lpr*lmass*( - 4.*(C12A + C12B) + C0A  + C0B );
      a -= loop;
      d += loop;
      f +=  2.*(lc - couplings[i].second)*lmass*(C0A+C0B);
      // a +=  2.*lpr*lmass*(2.*C12A-C0A+2.*C12B-C0B+
      // 			  (1. - pv1s*(C22A+C11B)- pv2s*(C11A+C22B) + mls*(C0A+C0B) )/p1p2);
      // b += 8.*lpr*lmass*(2.*C22A + C2A  +2.*C11B + C1B);
      // c += 2.*lpr*lmass*(- 4.*(C12A + C12B) - 2.*(C2A + C1A + C2B + C1B) - C0A - C0B);
      // d +=  2.*lpr*lmass*( - 4.*(C12A + C12B) + C0A  + C0B );
      // e +=  4.*lpr*lmass*( 2.*(C11A + C22B) + C1A  + C2B);
    }
    else if(type[i] == PDT::Spin1) {
      Complex B0W(LT::B0(ps2 ,mls,mls));
      double mr(sqr(p1->mass()/masses[i]));
      Complex loop = 2.*lc*( -(6.+mr)*(C12A+C12B) + 4.*(C0A +C0B));
      a -= loop;
      d += loop;
      // a += lc*(-8.*(C0A+C0B) +(6.+mr)*(2.*(C00A+C00B)-B0W)/p1p2);
      // b += lc*(6.+mr)*(2.*(C22A+C11B)+C2A+C1B);
      // c += 0.5*lc*(6.+mr)*(  - 4.*(C12A+C12B) - 2.*(C1A+C2A+C1B+C2B) - C0A - C0B );  
      // d += 2.*lc*( -(6.+mr)*(C12A+C12B) + 4.*(C0A +C0B)); 
      // e += lc*(6.+mr)*( 2.*(C11A+C22B) + C1A + C2B );
    }    
    else if(type[i] == PDT::Spin0) {
      Complex loop = 4.*lc*(C12A+C12B);
      a -= loop;
      d += loop;
      // a += -2.*lc*( 2.*(C00A+C00B) - LT::B0(ps2,mls,mls))/p1p2;
      // b += -2.*lc* ( 2.*(C22A + C11B) + C2A + C1B);
      // c +=    -lc*(  - 4.*(C12A + C12B) - 2.*(C2A + C1A + C2B + C1B) - C0A - C0B);
      // d +=  4.*lc*(C12A+C12B);
      // e +=    -lc*( 2.*(C11A + C22B) + C1A + C2B );
    }
    else {
      throw Helicity::HelicityConsistencyError() 
	<< "SVVLoopVertex::setCoupling - Incorrect particle in SVV loop. "
	<< "Spin: " << type[i]
	<< Exception::runerror;
    }
  }
  //Looptools defines integrals differently
  double fact = 1./16./sqr(Constants::pi);
  a00(fact*a);
  a11(fact*b);
  a12(fact*c);
  a21(fact*d);
  a22(fact*e);
  aEp(fact*f);
}
