// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleVVSLoopVertex class.
//

#include "SimpleVVSLoopVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

IBPtr SimpleVVSLoopVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleVVSLoopVertex::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<SimpleVVSLoopVertex> SimpleVVSLoopVertex::initSimpleVVSLoopVertex;
// Definition of the static class description member.

void SimpleVVSLoopVertex::Init() {

  static ClassDocumentation<SimpleVVSLoopVertex> documentation
    ("The SimpleVVSLoopVertex class calculates the tensor integral"
     " coefficients using Looptools.");

}
  
void SimpleVVSLoopVertex::setCoupling(Energy2 q2, tcPDPtr, tcPDPtr,tcPDPtr) {
  Complex loop(0.);
  for(unsigned int i = 0; i < masses.size(); ++i) {
    loop += A1(q2,sqr(masses[i]));
//    tmp+=W3(sqrt(q2),masses[i]);
  }
  a00(loop);
  a11(0.0);
  a12(0.0);
  a21(-loop);
  a22(0.0);
  aEp(0.0);
}

Complex SimpleVVSLoopVertex::A1(Energy2 s,Energy2 mf2) const {
  return mf2/s*(4.-W2(s,mf2)*(1.-4.*mf2/s));
}

Complex SimpleVVSLoopVertex::W2(Energy2 s, Energy2 mf2) const {
  double pi = Constants::pi;
  Complex ac(0.);
  double root=0.5*sqrt(abs(s)/mf2);

  if(s < ZERO)
    ac = sqr(asinh(root));
  else if(root<1.)
    ac = -sqr(asin(root));
  else
    ac = sqr(acosh(root))-0.25*sqr(pi)-pi*acosh(root)*Complex(0.,1.);
  return 4.*ac;
}

Complex SimpleVVSLoopVertex::W3(Energy mh, Energy mf) const {
  double pi = Constants::pi;
  Complex ac(0.);
  Complex i(0.,1.);
  double ratio = mf/mh;
  double ratio2 = sqr(ratio);

  if(2.*mf > mh) {
    ac = -2.*sqr(asin(0.5/ratio));
  } else if(2.*mf < mh) {
    double etalog = log((1.0+sqrt(1.-4.*ratio2))/(1.-sqrt(1.-4.*ratio2)));
    ac = 0.5 * (sqr(etalog) - sqr(pi)) + i*(pi*etalog);
  } else {
    ac = 0.5 * (- sqr(pi));
  }
  ac = 3.*ratio2*(2. + (4.*ratio2-1.)*ac);

  return ac;
}
