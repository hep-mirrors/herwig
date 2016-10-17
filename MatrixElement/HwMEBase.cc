// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwMEBase class.
//

#include "HwMEBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Shower/Core/Base/Branching.h"
#include "Herwig/Shower/RealEmissionProcess.h"

using namespace Herwig;

void HwMEBase::persistentOutput(PersistentOStream & os) const {
  os << massOption_ << rescaleOption_;
}

void HwMEBase::persistentInput(PersistentIStream & is, int) {
  is >> massOption_ >> rescaleOption_;
}

AbstractClassDescription<HwMEBase> HwMEBase::initHwMEBase;
// Definition of the static class description member.

void HwMEBase::Init() {

  static ClassDocumentation<HwMEBase> documentation
    ("The HwMEBase class is the base class for matrix elements in Herwig"
     " and provides the virtual members for hard radiation corrections in the"
     " shower.");

}

int HwMEBase::nDim() const {
  unsigned ndim = 1;
  for(unsigned int ix=0;ix<massOption_.size();++ix)
    if(massOption_[ix]==2) ++ndim;
  return ndim;
}

CrossSection HwMEBase::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

void HwMEBase::setKinematics() {
  MEBase::setKinematics();
  lastTHat_ = (meMomenta()[0] - meMomenta()[2]).m2();
  lastUHat_ = (meMomenta()[1] - meMomenta()[2]).m2();
  lastPhi_ = meMomenta()[2].phi();
}

bool HwMEBase::generateMasses(vector<Energy> & masses, double & mjac,
			      const double *r) {
  assert(massOption_.size()+2==mePartonData().size());
  mjac = 1.;
  masses.clear();
  masses.resize(massOption_.size(),ZERO);
  Energy ecm = sqrt(sHat());
  Energy emin(ZERO);
  int noff(0);
  for(unsigned int ix=0;ix<massOption_.size();++ix) {
    if(massOption_[ix]==1) {
      masses[ix] = mePartonData()[ix+2]->hardProcessMass();
      emin += masses[ix];
    }
    else if (massOption_[ix]==2) {
      emin += mePartonData()[ix+2]->massMin();
      ++noff;
    }
  }
  // check allowed
  if(emin>ecm) return false;
  // if nothing off-shell return
  if(noff==0) return true;
  int iloc = nDim()-noff;
  emin = ecm - emin;
  // generate the masses
  for(unsigned int ix=0;ix<massOption_.size();++ix) {
    if(massOption_[ix]!=2) continue;
    Energy mmin = mePartonData()[ix+2]->massMin();
    emin += mmin;
    Energy mmax = min(mePartonData()[ix+2]->massMax(),emin);
    if(mmin>mmax) return false;
    tGenericMassGeneratorPtr gen = mePartonData()[ix+2]->massGenerator() ?
      dynamic_ptr_cast<tGenericMassGeneratorPtr>(mePartonData()[ix+2]->massGenerator()) :
      tGenericMassGeneratorPtr();
    if(gen) {
      double jtemp(0.);
      masses[ix] = gen->mass(jtemp,*mePartonData()[ix+2],mmin,mmax,r[iloc]);
      mjac *= jtemp;
    }
    else {
      Energy mon(mePartonData()[ix+2]->hardProcessMass());
      Energy width(mePartonData()[ix+2]->width());
      double rhomin = atan2((sqr(mmin)-sqr(mon)), mon*width);
      double rhomax = atan2((sqr(mmax)-sqr(mon)), mon*width);
      masses[ix] = sqrt(mon*width*tan(rhomin+r[iloc]*(rhomax-rhomin))+sqr(mon));
      mjac *= (rhomax-rhomin)/Constants::pi;
    }
    emin -= masses[ix];
    if(emin<ZERO) return false;
    ++iloc;
  }
  return true;
}

bool HwMEBase::generateKinematics(const double * r) {
  jacobian(1.);
  vector<Energy> masses;
  double mjac(0.);
  if(!generateMasses(masses,mjac,r)) return false;
  // set up the momenta
  for ( int i = 2, N = meMomenta().size(); i < N; ++i ) {
    meMomenta()[i] = Lorentz5Momentum(masses[i-2]);
  }
  double ctmin = -1.0, ctmax = 1.0;
  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } 
  catch ( ImpossibleKinematics ) {
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
  jacobian((pq/sHat())*Constants::pi*jacobian()*mjac);
  // compute the rescaled momenta
  return rescaleMomenta(meMomenta(),mePartonData());
}

bool HwMEBase::rescaleMomenta(const vector<Lorentz5Momentum> & momenta,
			      const cPDVector & data) {
  assert(momenta.size()==4&&data.size()==4);
  // default just use the ones we generated
  rescaledMomenta_=momenta;
  if(rescaleOption_==1) return true;
  Energy mnew[2] = {0*MeV, ZERO};
  if(rescaleOption_==0) {
    mnew[0] = ZERO;
    mnew[1] = ZERO;
  }
  else if(rescaleOption_==2) {
    mnew[0] = data[2]->hardProcessMass();
    mnew[1] = data[3]->hardProcessMass();
  }
  else if(rescaleOption_==3) {
    if(abs(data[2]->id())!=abs(data[3]->id())) return true;
    mnew[0] = 0.5*(momenta[2].mass()+momenta[3].mass());
    mnew[1] = mnew[0];
  } 
  else {
    assert(false);
  }
  Lorentz5Momentum pcm(momenta[2]+momenta[3]);
  Energy m0=pcm.m();
  if(m0<mnew[0]+mnew[1]) return false;
  Boost bv = pcm.boostVector();
  rescaledMomenta_[2].boost(bv);
  rescaledMomenta_[2].setMass(mnew[0]);
  rescaledMomenta_[2].setE(0.5*(sqr(m0)+sqr(mnew[0])-sqr(mnew[1]))/m0);
  if(rescaledMomenta_[2].t()-rescaledMomenta_[2].mass()>1e-10*(rescaledMomenta_[2].t()+rescaledMomenta_[2].mass()))
    rescaledMomenta_[2].rescaleRho();
  else {
    rescaledMomenta_[2].setX(ZERO);
    rescaledMomenta_[2].setY(ZERO);
    rescaledMomenta_[2].setZ(ZERO);
  }
  rescaledMomenta_[2].boost(-bv);
  rescaledMomenta_[3].boost(bv);
  rescaledMomenta_[3].setMass(mnew[1]);
  rescaledMomenta_[3].setE(0.5*(sqr(m0)-sqr(mnew[0])+sqr(mnew[1]))/m0);
  if(rescaledMomenta_[3].t()-rescaledMomenta_[3].mass()>1e-10*(rescaledMomenta_[3].t()+rescaledMomenta_[3].mass()))
    rescaledMomenta_[3].rescaleRho();
  else {
    rescaledMomenta_[3].setX(ZERO);
    rescaledMomenta_[3].setY(ZERO);
    rescaledMomenta_[3].setZ(ZERO);
  }
  rescaledMomenta_[3].boost(-bv);
  return true;
}

double HwMEBase::getCosTheta(double ctmin, double ctmax, const double r) {
  double cth = 0.0;
  static const double eps = 1.0e-6;
  if ( 1.0 + ctmin <= eps && 1.0 - ctmax <= eps ) {
    jacobian(jacobian()*(ctmax - ctmin));
    cth = ctmin + r*(ctmax - ctmin);
  } else if (  1.0 + ctmin <= eps ) {
    cth = 1.0 - (1.0 - ctmax)*pow((1.0 - ctmin)/(1.0 - ctmax), r);
    jacobian(jacobian()*log((1.0 - ctmin)/(1.0 - ctmax))*(1.0 - cth));
  } else if (  1.0 - ctmax <= eps ) {
    cth = -1.0 + (1.0 + ctmin)*pow((1.0 + ctmax)/(1.0 + ctmin), r);
    jacobian(jacobian()*log((1.0 + ctmax)/(1.0 + ctmin))*(1.0 + cth));
  } else {
    double zmin = 0.5*(1.0 - ctmax);
    double zmax = 0.5*(1.0 - ctmin);
    double A1 = -ctmin/(zmax*(1.0-zmax));
    double A0 = -ctmax/(zmin*(1.0-zmin));
    double A = r*(A1 - A0) + A0;
    double z = A < 2.0? 2.0/(sqrt(sqr(A) + 4.0) + 2 - A):
      0.5*(A - 2.0 + sqrt(sqr(A) + 4.0))/A;
    cth = 1.0 - 2.0*z;
    jacobian(jacobian()*2.0*(A1 - A0)*sqr(z)*sqr(1.0 - z)/(sqr(z) + sqr(1.0 - z)));
  }
  return cth;
}

bool HwMEBase::softMatrixElementVeto(ShowerProgenitorPtr,
				     ShowerParticlePtr,Branching) {
  assert(false);
  return false;
}

RealEmissionProcessPtr HwMEBase::generateHardest(RealEmissionProcessPtr,ShowerInteraction) {
  assert(false);
  return RealEmissionProcessPtr();
}

RealEmissionProcessPtr HwMEBase::applyHardMatrixElementCorrection(RealEmissionProcessPtr) {
  assert(false);
  return RealEmissionProcessPtr();
}

void HwMEBase::initializeMECorrection(RealEmissionProcessPtr , double & ,
				      double & ) {
  assert(false);
}
