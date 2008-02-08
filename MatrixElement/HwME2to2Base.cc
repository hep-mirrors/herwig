// -*- C++ -*-
//
// HwME2to2Base.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwME2to2Base class.
//
#include "HwME2to2Base.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"

using namespace Herwig;
using namespace ThePEG;

int HwME2to2Base::nDim() const {
  return 1 + bool(_massopt1==2) + bool(_massopt2==2);
}

void HwME2to2Base::setKinematics() {
  ME2to2Base::setKinematics();
}

bool HwME2to2Base::generateKinematics(const double * r) {
  Energy mass[2];
  // generate the masses of the particles
  if(_massopt1==2&&_massopt2==2) {
    throw Exception() << "Both off-shell not yet handled in "
		      << "HwME2to2Base::generateKinematics "
		      << Exception::runerror;
  }
  else if(_massopt1==2) {
    mass[1] = _massopt1==0 ? 0.*GeV : mePartonData()[0]->mass();
    throw Exception() << "First particle off-shell not yet handled in "
		      << "HwME2to2Base::generateKinematics "
		      << Exception::runerror;
  }
  else if(_massopt2==2) {
    mass[0] = _massopt1==0 ? 0.*GeV : mePartonData()[0]->mass();
    throw Exception() << "Second particle off-shell not yet handled in "
		      << "HwME2to2Base::generateKinematics "
		      << Exception::runerror;
  }
  else {
    mass[0] = _massopt1==0 ? 0.*GeV : mePartonData()[2]->mass();
    mass[1] = _massopt1==0 ? 0.*GeV : mePartonData()[3]->mass();
  }
  for ( int i = 2, N = meMomenta().size(); i < N; ++i ) {
    meMomenta()[i] = Lorentz5Momentum(mass[i-2]);
  }

  double ctmin = -1.0;
  double ctmax = 1.0;
  Energy q = 0.0*GeV;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } catch ( ImpossibleKinematics ) {
    return false;
  }

  Energy e = sqrt(sHat())/2.0;
     	    
  Energy2 m22 = meMomenta()[2].mass2();
  Energy2 m32 = meMomenta()[3].mass2();
  Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e0e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 e1e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 pq = 2.0*e*q;

  Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
  if ( thmin > 0.0*GeV2 ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
  if ( thmin > 0.0*GeV2 ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if ( thmin > 0.0*GeV2 ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
  if ( thmin > 0.0*GeV2 ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);

  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
   		     lastCuts().minKT(mePartonData()[3]));
  if ( ptmin > 0.0*GeV ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax, sqrt(ctm));
  }

  double ymin2 = lastCuts().minYStar(mePartonData()[2]);
  double ymax2 = lastCuts().maxYStar(mePartonData()[2]);
  double ymin3 = lastCuts().minYStar(mePartonData()[3]);
  double ymax3 = lastCuts().maxYStar(mePartonData()[3]);
  double ytot = lastCuts().Y() + lastCuts().currentYHat();
  if ( ymin2 + ytot > -0.9*Constants::MaxRapidity )
    ctmin = max(ctmin, sqrt(sqr(q) +  m22)*tanh(ymin2)/q);
  if ( ymax2 + ytot < 0.9*Constants::MaxRapidity )
    ctmax = min(ctmax, sqrt(sqr(q) +  m22)*tanh(ymax2)/q);
  if ( ymin3 + ytot > -0.9*Constants::MaxRapidity )
    ctmax = min(ctmax, sqrt(sqr(q) +  m32)*tanh(-ymin3)/q);
  if ( ymax3 + ytot < 0.9*Constants::MaxRapidity )
    ctmin = max(ctmin, sqrt(sqr(q) +  m32)*tanh(-ymax3)/q);
  
  if ( ctmin >= ctmax ) return false;
    
  double cth = getCosTheta(ctmin, ctmax, r);
  Energy pt = q*sqrt(1.0-sqr(cth));
  phi(rnd(2.0*Constants::pi));
  meMomenta()[2].setVect(Momentum3( pt*sin(phi()),  pt*cos(phi()),  q*cth));
  meMomenta()[3].setVect(Momentum3(-pt*sin(phi()), -pt*cos(phi()), -q*cth));

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  vector<LorentzMomentum> out(2);
  out[0] = meMomenta()[2];
  out[1] = meMomenta()[3];
  tcPDVector tout(2);
  tout[0] = mePartonData()[2];
  tout[1] = mePartonData()[3];
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;

  tHat(pq*cth + m22 - e0e2);
  uHat(m22 + m32 - sHat() - tHat());
  jacobian((pq/sHat())*Constants::pi*jacobian());
  return true;
}

void HwME2to2Base::persistentOutput(PersistentOStream & os) const {
  os << _massopt1 << _massopt2;
}

void HwME2to2Base::persistentInput(PersistentIStream & is, int) {
  is >> _massopt1 >> _massopt2;
}

AbstractClassDescription<HwME2to2Base> HwME2to2Base::initHwME2to2Base;
// Definition of the static class description member.

void HwME2to2Base::Init() {

  static ClassDocumentation<HwME2to2Base> documentation
    ("The ME2to2Base class may be used as a base class "
     "for all \\f$2\\rightarrow 2\\f$ matrix elements.");

}
