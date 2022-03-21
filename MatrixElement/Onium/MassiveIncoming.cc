// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MassiveIncoming class.
//

#include "MassiveIncoming.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"

using namespace Herwig;

bool MassiveIncoming::generateKinematics(const double * r) {
  jacobian(1.);
  Energy ecm = sqrt(sHat());
  Energy mBc(ZERO);
  Energy mout = mePartonData()[3]->hardProcessMass();
  if(mOpt_==0 || !massGen_) {
    mBc=mePartonData()[2]->mass();
  }
  else {
    Energy mmin = mePartonData()[2]->massMin(), mmax = min(mePartonData()[2]->massMax(),ecm-mout);
    double jtemp(0.);
    mBc = massGen_->mass(jtemp,*mePartonData()[2],mmin,mmax,r[1]);
    jacobian(jacobian()*jtemp);
  }
  // set the masses
  meMomenta()[2].setMass(mBc);
  meMomenta()[3].setMass(mout );

  double ctmin = -1.0, ctmax = 1.0;
  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } 
  catch ( ImpossibleKinematics & e) {
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
  if ( thmin > ZERO ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if ( thmin > ZERO ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
  if ( thmin > ZERO ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);

  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
   		     lastCuts().minKT(mePartonData()[3]));
  if ( ptmin > ZERO ) {
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

  double cth = getCosTheta(ctmin, ctmax, r[0]);
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
  // now compute the rescaled momenta we need for the ME
  vector<Lorentz5Momentum> rescaled(4);
  if(mePartonData()[0]->id()==ParticleID::g) {
    double rr = mePartonData()[1]->mass()/mePartonData()[3]->mass();
    // masses
    rescaled[0].setMass(          ZERO);
    rescaled[1].setMass(rr/(1.+rr)*mBc);
    rescaled[2].setMass(           mBc);
    rescaled[3].setMass(1./(1.+rr)*mBc);
  }
  else {
    double rr = mePartonData()[1]->mass()/mePartonData()[0]->mass();
    // masses
    rescaled[0].setMass(1./(1.+rr)*mBc);
    rescaled[1].setMass(rr/(1.+rr)*mBc);
    rescaled[2].setMass(           mBc);
    rescaled[3].setMass(mout          );
  }
  // incoming
  Energy pin = SimplePhaseSpace::getMagnitude(sHat(), rescaled[0].mass(), rescaled[1].mass());
  rescaled[0].setZ(pin); rescaled[1].setZ(-pin);
  rescaled[0].setT(0.5*(sHat()+sqr(rescaled[0].mass())-sqr(rescaled[1].mass()))/ecm);
  rescaled[1].setT(0.5*(sHat()-sqr(rescaled[0].mass())+sqr(rescaled[1].mass()))/ecm);
  // outgoing
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), rescaled[2].mass(), rescaled[3].mass());
  } 
  catch ( ImpossibleKinematics & e) {
    return false;
  }
  pt = q*sqrt(1.0-sqr(cth));
  rescaled[2].setVect(Momentum3( pt*sin(phi()),  pt*cos(phi()),  q*cth));
  rescaled[3].setVect(Momentum3(-pt*sin(phi()), -pt*cos(phi()), -q*cth));
  rescaled[2].rescaleEnergy();
  rescaled[3].rescaleEnergy();
  rescaledMomenta(rescaled);
  return true;
}

CrossSection MassiveIncoming::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

void MassiveIncoming::persistentOutput(PersistentOStream & os) const {
  os << mOpt_ << massGen_;
}

void MassiveIncoming::persistentInput(PersistentIStream & is, int) {
  is >> mOpt_ >> massGen_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MassiveIncoming,HwMEBase>
describeHerwigMassiveIncoming("Herwig::MassiveIncoming",
			      "HwOniumParameters.so HwMEHadronOnium.so");

void MassiveIncoming::Init() {

  static ClassDocumentation<MassiveIncoming> documentation
    ("The MassiveIncoming class handles the kinematics if thw incoming partons are massive");
  
  static Switch<MassiveIncoming,unsigned int> interfaceMassOption
    ("MassOption",
     "Mass of the treatment of mas of the B_c state",
     &MassiveIncoming::mOpt_, 1, false, false);
  static SwitchOption interfaceMassOptionOnShell
    (interfaceMassOption,
     "OnShell",
     "Use the on-shell mass",
     0);
  static SwitchOption interfaceMassOptionOffShell
    (interfaceMassOption,
     "OffShell",
     "Use an off-shell mass generated by the MassGenerator object for the B_c state.",
     1);
}

