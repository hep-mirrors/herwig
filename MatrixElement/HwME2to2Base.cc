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
#include "Herwig++/PDT/GenericMassGenerator.h"

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
  double mjac(1.);
  if(_massopt1==2&&_massopt2==2) {
    // minimum masses of the outgoing particles
    Energy ecm = sqrt(sHat());
    Energy mmin[2]  = {mePartonData()[2]->massMin(),
		       mePartonData()[3]->massMin()};
    // not kinematically possible return
    if(ecm<mmin[0]+mmin[1]) return false;
    // maximum masses of the outgoing particles, including kinematic limit
    Energy mmax[2]  = {min(mePartonData()[2]->massMax(),ecm-mmin[1]),
		       min(mePartonData()[3]->massMax(),ecm-mmin[0])};
    // generate the mass of the first particle
    if(mmax[0]<mmin[0]) return false;
    tGenericMassGeneratorPtr gen1= mePartonData()[2]->massGenerator() ?
      dynamic_ptr_cast<tGenericMassGeneratorPtr>(mePartonData()[2]->massGenerator()) :
      tGenericMassGeneratorPtr();
    if(gen1) {
      double jtemp(0.);
      mass[0] = gen1->mass(jtemp,*mePartonData()[2],mmin[0],mmax[0],r[1]);
      mjac *= jtemp;
    }
    else {
      Energy mon(mePartonData()[2]->mass()),width(mePartonData()[2]->width());
      double rhomin = atan((sqr(mmin[0])-sqr(mon))/mon/width);
      double rhomax = atan((sqr(mmax[0])-sqr(mon))/mon/width);
      mass[0] = sqrt(mon*width*tan(rhomin+r[1]*(rhomax-rhomin))+sqr(mon));
      mjac *= (rhomax-rhomin)/Constants::pi;
    }
    // generate the mass of the second particle
    mmax[1] = min(mmax[1],ecm-mass[0]);
    if(mmax[1]<mmin[1]) return false;
    tGenericMassGeneratorPtr gen2 = mePartonData()[3]->massGenerator() ?
      dynamic_ptr_cast<tGenericMassGeneratorPtr>(mePartonData()[3]->massGenerator()) :
      tGenericMassGeneratorPtr();
    if(gen2) {
      double jtemp(0.);
      mass[1] = gen2->mass(jtemp,*mePartonData()[3],mmin[1],mmax[1],r[2]);
      mjac *= jtemp;
    }
    else {
      Energy mon(mePartonData()[3]->mass()),width(mePartonData()[3]->width());
      double rhomin = atan((sqr(mmin[1])-sqr(mon))/mon/width);
      double rhomax = atan((sqr(mmax[1])-sqr(mon))/mon/width);
      mass[1] = sqrt(mon*width*tan(rhomin+r[2]*(rhomax-rhomin))+sqr(mon));
      mjac *= (rhomax-rhomin);
    }
  }
  else if(_massopt1==2) {
    mass[1] = _massopt2==0 ? 0.*GeV : mePartonData()[3]->mass();
    Energy ecm = sqrt(sHat());
    Energy mmin  = mePartonData()[2]->massMin();
    // not kinematically possible return
    if(ecm<mmin+mass[1]) return false;
    // maximum masses of the outgoing particle, including kinematic limit
    Energy mmax  = min(mePartonData()[2]->massMax(),ecm-mass[1]);
    // generate the mass of the particle
    if(mmax<mmin) return false;
    tGenericMassGeneratorPtr gen = mePartonData()[2]->massGenerator() ?
      dynamic_ptr_cast<tGenericMassGeneratorPtr>(mePartonData()[2]->massGenerator()) :
      tGenericMassGeneratorPtr();
    if(gen) {
      double jtemp(0.);
      mass[0] = gen->mass(jtemp,*mePartonData()[2],mmin,mmax,r[1]);
      mjac *= jtemp;
    }
    else {
      Energy mon(mePartonData()[2]->mass()),width(mePartonData()[2]->width());
      double rhomin = atan((sqr(mmin)-sqr(mon))/mon/width);
      double rhomax = atan((sqr(mmax)-sqr(mon))/mon/width);
      mass[0] = sqrt(mon*width*tan(rhomin+r[1]*(rhomax-rhomin))+sqr(mon));
      mjac *= (rhomax-rhomin)/Constants::pi;
    }
  }
  else if(_massopt2==2) {
    mass[0] = _massopt1==0 ? 0.*GeV : mePartonData()[2]->mass();
    Energy ecm = sqrt(sHat());
    Energy mmin  = mePartonData()[3]->massMin();
    // not kinematically possible return
    if(ecm<mmin+mass[0]) return false;
    // maximum masses of the outgoing particle, including kinematic limit
    Energy mmax  = min(mePartonData()[2]->massMax(),ecm-mass[0]);
    // generate the mass of the particle
    if(mmax<mmin) return false;
    tGenericMassGeneratorPtr gen = mePartonData()[3]->massGenerator() ?
      dynamic_ptr_cast<tGenericMassGeneratorPtr>(mePartonData()[3]->massGenerator()) :
      tGenericMassGeneratorPtr();
    if(gen) {
      double jtemp(0.);
      mass[1] = gen->mass(jtemp,*mePartonData()[3],mmin,mmax,r[1]);
      mjac *= jtemp;
    }
    else {
      Energy mon(mePartonData()[3]->mass()),width(mePartonData()[3]->width());
      double rhomin = atan((sqr(mmin)-sqr(mon))/mon/width);
      double rhomax = atan((sqr(mmax)-sqr(mon))/mon/width);
      mass[1] = sqrt(mon*width*tan(rhomin+r[1]*(rhomax-rhomin))+sqr(mon));
      mjac *= (rhomax-rhomin)/Constants::pi;
    }
  }
  else {
    mass[0] = _massopt1==0 ? 0.*GeV : mePartonData()[2]->mass();
    mass[1] = _massopt2==0 ? 0.*GeV : mePartonData()[3]->mass();
  }
  // set up the momenta
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
  jacobian((pq/sHat())*Constants::pi*jacobian()*mjac);
  // compute the rescaled momenta 
  return rescaleMomenta(meMomenta(),mePartonData());
}

void HwME2to2Base::persistentOutput(PersistentOStream & os) const {
  os << _massopt1 << _massopt2 << _rescaleOption;
}

void HwME2to2Base::persistentInput(PersistentIStream & is, int) {
  is >> _massopt1 >> _massopt2 >> _rescaleOption;
}

AbstractClassDescription<HwME2to2Base> HwME2to2Base::initHwME2to2Base;
// Definition of the static class description member.

void HwME2to2Base::Init() {

  static ClassDocumentation<HwME2to2Base> documentation
    ("The ME2to2Base class may be used as a base class "
     "for all \\f$2\\rightarrow 2\\f$ matrix elements.");

}

bool HwME2to2Base::rescaleMomenta(const vector<Lorentz5Momentum> & momenta,
				  const cPDVector & data) {
  assert(momenta.size()==4&&data.size()==4);
  // default just use the ones we generated
  _rescaledMomenta=momenta;
  if(_rescaleOption==1) return true;
  Energy mnew[2];
  if(_rescaleOption==0) {
    mnew[0] = 0.*GeV;
    mnew[1] = 0.*GeV;
  }
  else if(_rescaleOption==2) {
    mnew[0] = data[2]->mass();
    mnew[1] = data[3]->mass();
  }
  else if(_rescaleOption==3) {
    if(abs(data[2]->id())!=abs(data[3]->id())) return true;
    mnew[0] = 0.5*(momenta[2].mass()+momenta[3].mass());
    mnew[1] = mnew[0];
  }
  Lorentz5Momentum pcm(momenta[2]+momenta[3]);
  Energy m0=pcm.m();
  if(m0<mnew[0]+mnew[1]) return false;
  Boost bv = pcm.boostVector();
  _rescaledMomenta[2].boost(bv);
  _rescaledMomenta[2].setMass(mnew[0]);
  _rescaledMomenta[2].setE(0.5*(sqr(m0)+sqr(mnew[0])-sqr(mnew[1]))/m0);
  _rescaledMomenta[2].rescaleRho();
  _rescaledMomenta[2].boost(-bv);
  _rescaledMomenta[3].boost(bv);
  _rescaledMomenta[3].setMass(mnew[1]);
  _rescaledMomenta[3].setE(0.5*(sqr(m0)-sqr(mnew[0])+sqr(mnew[1]))/m0);
  _rescaledMomenta[3].rescaleRho();
  _rescaledMomenta[2].boost(-bv);
  return true;
}
