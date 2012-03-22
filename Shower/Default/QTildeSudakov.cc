// -*- C++ -*-
//
// QTildeSudakov.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeSudakov class.
//

#include "QTildeSudakov.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/Default/FS_QTildeShowerKinematics1to2.h"
#include "Herwig++/Shower/Default/IS_QTildeShowerKinematics1to2.h"
#include "Herwig++/Shower/Default/Decay_QTildeShowerKinematics1to2.h"

using namespace Herwig;

NoPIOClassDescription<QTildeSudakov> QTildeSudakov::initQTildeSudakov;
// Definition of the static class description member.

void QTildeSudakov::Init() {

  static ClassDocumentation<QTildeSudakov> documentation
    ("The QTildeSudakov class implements the Sudakov form factor for ordering it"
     " qtilde");
}

bool QTildeSudakov::guessTimeLike(Energy2 &t,Energy2 tmin,double enhance) {
  Energy2 told = t;
  // calculate limits on z and if lower>upper return
  if(!computeTimeLikeLimits(t)) return false;
  // guess values of t and z
  t = guesst(told,0,ids_,enhance,ids_[1]==ids_[2]);
  z(guessz(0,ids_)); 
  // actual values for z-limits
  if(!computeTimeLikeLimits(t)) return false;
  if(t<tmin) {
    t=-1.0*GeV2;
    return false;
  }
  else
    return true; 
} 

bool QTildeSudakov::guessSpaceLike(Energy2 &t, Energy2 tmin, const double x,
				   double enhance) {
  Energy2 told = t;
  // calculate limits on z if lower>upper return
  if(!computeSpaceLikeLimits(t,x)) return false;
  // guess values of t and z
  t = guesst(told,1,ids_,enhance,ids_[1]==ids_[2]); 
  z(guessz(1,ids_)); 
  // actual values for z-limits
  if(!computeSpaceLikeLimits(t,x)) return false;
  if(t<tmin) {
    t=-1.0*GeV2;
    return false;
  }
  else
    return true; 
} 

bool QTildeSudakov::PSVeto(const Energy2 t) {
  // still inside PS, return true if outside
  // check vs overestimated limits
  if(z() < zLimits().first || z() > zLimits().second) return true;
  // compute the pts
  Energy2 pt2=sqr(z()*(1.-z()))*t-masssquared_[1]*(1.-z())-masssquared_[2]*z();
  if(ids_[0]!=ParticleID::g) pt2+=z()*(1.-z())*masssquared_[0];
  // if pt2<0 veto
  if(pt2<pT2min()) return true;
  // otherwise calculate pt and return
  pT(sqrt(pt2));
  return false;
}
 
ShoKinPtr QTildeSudakov::generateNextTimeBranching(const Energy startingScale,
						   const IdList &ids,const bool cc,
						   double enhance) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  q_ = ZERO;
  z(0.);
  phi(0.); 
  // perform initialization
  Energy2 tmax(sqr(startingScale)),tmin;
  initialize(ids,tmin,cc);
  // check max > min
  if(tmax<=tmin) return ShoKinPtr();
  // calculate next value of t using veto algorithm
  Energy2 t(tmax);
  do {
    if(!guessTimeLike(t,tmin,enhance)) break;
  }
  while(PSVeto(t) || SplittingFnVeto(z()*(1.-z())*t,ids,true) || 
	alphaSVeto(sqr(z()*(1.-z()))*t));
  if(t > ZERO) q_ = sqrt(t);
  else q_ = -1.*MeV;
  phi(Constants::twopi*UseRandom::rnd());
  if(q_ < ZERO) return ShoKinPtr();
  // return the ShowerKinematics object
  return createFinalStateBranching(q_,z(),phi(),pT()); 
}

ShoKinPtr QTildeSudakov::
generateNextSpaceBranching(const Energy startingQ,
			   const IdList &ids,
			   double x,bool cc,
			   double enhance,
			   Ptr<BeamParticleData>::transient_const_pointer beam) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  q_ = ZERO;
  z(0.);
  phi(0.);
  // perform the initialization
  Energy2 tmax(sqr(startingQ)),tmin;
  initialize(ids,tmin,cc);
  // check max > min
  if(tmax<=tmin) return ShoKinPtr();
  // extract the partons which are needed for the PDF veto
  // Different order, incoming parton is id =  1, outgoing are id=0,2
  tcPDPtr parton0 = getParticleData(ids[0]);
  tcPDPtr parton1 = getParticleData(ids[1]);
  if(cc) {
    if(parton0->CC()) parton0 = parton0->CC();
    if(parton1->CC()) parton1 = parton1->CC();
  }
  // calculate next value of t using veto algorithm
  Energy2 t(tmax),pt2(ZERO);
  do {
    if(!guessSpaceLike(t,tmin,x,enhance)) break;
    pt2=sqr(1.-z())*t-z()*masssquared_[2];
  }
  while(z() > zLimits().second || 
	SplittingFnVeto((1.-z())*t/z(),ids,true) || 
	alphaSVeto(sqr(1.-z())*t) || 
	PDFVeto(t,x,parton0,parton1,beam) || pt2 < pT2min() );
  if(t > ZERO && zLimits().first < zLimits().second)  q_ = sqrt(t);
  else return ShoKinPtr();
  phi(Constants::twopi*UseRandom::rnd());
  pT(sqrt(pt2));
  // create the ShowerKinematics and return it
  return createInitialStateBranching(q_,z(),phi(),pT());
}

void QTildeSudakov::initialize(const IdList & ids, Energy2 & tmin,const bool cc) {
  ids_=ids;
  if(cc) {
    for(unsigned int ix=0;ix<ids.size();++ix) {
      if(getParticleData(ids[ix])->CC()) ids_[ix]*=-1;
    }
  }
  tmin = cutOffOption() != 2 ? ZERO : 4.*pT2min();
  masses_ = virtualMasses(ids);
  masssquared_.clear();
  for(unsigned int ix=0;ix<masses_.size();++ix) {
    masssquared_.push_back(sqr(masses_[ix]));
    if(ix>0) tmin=max(masssquared_[ix],tmin);
  }
}

ShoKinPtr QTildeSudakov::generateNextDecayBranching(const Energy startingScale,
						 const Energy stoppingScale,
						 const Energy minmass,
						 const IdList &ids,
						 const bool cc, 
						 double enhance) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to this method.
  q_ = Constants::MaxEnergy;
  z(0.);
  phi(0.); 
  // perform initialisation
  Energy2 tmax(sqr(stoppingScale)),tmin;
  initialize(ids,tmin,cc);
  tmin=sqr(startingScale);
  // check some branching possible
  if(tmax<=tmin) return ShoKinPtr();
  // perform the evolution
  Energy2 t(tmin),pt2(-MeV2);
  do {
    if(!guessDecay(t,tmax,minmass,enhance)) break;
    pt2 = sqr(1.-z())*(t-masssquared_[0])-z()*masssquared_[2];
  }
  while(SplittingFnVeto((1.-z())*t/z(),ids,true)|| 
	alphaSVeto(sqr(1.-z())*t) ||
	pt2<pT2min() ||
	t*(1.-z())>masssquared_[0]-sqr(minmass));
  if(t > ZERO) {
    q_ = sqrt(t);
    pT(sqrt(pt2));
  }
  else return ShoKinPtr();
  phi(Constants::twopi*UseRandom::rnd());
  // create the ShowerKinematics object
  return createDecayBranching(q_,z(),phi(),pT());
}

bool QTildeSudakov::guessDecay(Energy2 &t,Energy2 tmax, Energy minmass,
			       double enhance) {
  // previous scale
  Energy2 told = t;
  // overestimated limits on z
  if(tmax<masssquared_[0]) {
    t=-1.0*GeV2;
    return false;
  }
  Energy2 tm2 = tmax-masssquared_[0];
  Energy tm  = sqrt(tm2); 
  pair<double,double> limits=make_pair(sqr(minmass/masses_[0]),
				       1.-sqrt(masssquared_[2]+pT2min()+
					       0.25*sqr(masssquared_[2])/tm2)/tm
				       +0.5*masssquared_[2]/tm2);
  zLimits(limits);
  if(zLimits().second<zLimits().first) {
    t=-1.0*GeV2;
    return false;
  }
  // guess values of t and z
  t = guesst(told,2,ids_,enhance,ids_[1]==ids_[2]);
  z(guessz(2,ids_)); 
  // actual values for z-limits
  if(t<masssquared_[0])  {
    t=-1.0*GeV2;
    return false;
  }
  tm2 = t-masssquared_[0];
  tm  = sqrt(tm2); 
  limits=make_pair(sqr(minmass/masses_[0]),
		   1.-sqrt(masssquared_[2]+pT2min()+
			   0.25*sqr(masssquared_[2])/tm2)/tm
		   +0.5*masssquared_[2]/tm2);
  zLimits(limits);
  if(t>tmax||zLimits().second<zLimits().first) {
    t=-1.0*GeV2;
    return false;
  }
  else
    return true; 
} 

bool QTildeSudakov::computeTimeLikeLimits(Energy2 & t) {
  if (t < 1e-20 * GeV2) {
    t=-1.*GeV2;
    return false;
  }
  // special case for gluon radiating
  pair<double,double> limits;
  if(ids_[0]==ParticleID::g) {
    // no emission possible
    if(t<16.*masssquared_[1]) {
      t=-1.*GeV2;
      return false;
    }
    // overestimate of the limits
    limits.first  = 0.5*(1.-sqrt(1.-4.*sqrt((masssquared_[1]+pT2min())/t)));
    limits.second = 1.-limits.first;
  }
  // special case for radiated particle is gluon 
  else if(ids_[2]==ParticleID::g) {
    limits.first  =    sqrt((masssquared_[1]+pT2min())/t);
    limits.second = 1.-sqrt((masssquared_[2]+pT2min())/t);
  }
  else if(ids_[1]==ParticleID::g) {
    limits.second  =    sqrt((masssquared_[2]+pT2min())/t);
    limits.first   = 1.-sqrt((masssquared_[1]+pT2min())/t);
  }
  else {
    limits.first  =    (masssquared_[1]+pT2min())/t;
    limits.second = 1.-(masssquared_[2]+pT2min())/t; 
  }
  if(limits.first>=limits.second) {
    t=-1.*GeV2;
    return false;
  }
  zLimits(limits);
  return true;
}

bool QTildeSudakov::computeSpaceLikeLimits(Energy2 & t, double x) {
  if (t < 1e-20 * GeV2) {
    t=-1.*GeV2;
    return false;
  }
  pair<double,double> limits;
  // compute the limits
  limits.first = x;
  double yy = 1.+0.5*masssquared_[2]/t;
  limits.second = yy - sqrt(sqr(yy)-1.+pT2min()/t); 
  // return false if lower>upper
  zLimits(limits);
  if(limits.second<limits.first) {
    t=-1.*GeV2;
    return false;
  }
  else
    return true;
}

Energy QTildeSudakov::calculateScale(double zin, Energy pt, IdList ids,
				     unsigned int iopt) {
  Energy2 tmin;
  initialize(ids,tmin,false);
  // final-state branching
  if(iopt==0) {
    Energy2 scale=(sqr(pt)+masssquared_[1]*(1.-zin)+masssquared_[2]*zin);
    if(ids[0]!=ParticleID::g) scale -= zin*(1.-zin)*masssquared_[0];
    scale /= sqr(zin*(1-zin));
    return scale<=ZERO ? sqrt(tmin) : sqrt(scale);
  }
  else if(iopt==1) {
    Energy2 scale=(sqr(pt)+zin*masssquared_[2])/sqr(1.-zin);
    return scale<=ZERO ? sqrt(tmin) : sqrt(scale);
  }
  else if(iopt==2) {
    Energy2 scale = (sqr(pt)+zin*masssquared_[2])/sqr(1.-zin)+masssquared_[0];
    return scale<=ZERO ? sqrt(tmin) : sqrt(scale);
  }
  else {
    throw Exception() << "Unknown option in QTildeSudakov::calculateScale() "
		      << "iopt = " << iopt << Exception::runerror;
  }
}

ShoKinPtr QTildeSudakov::createFinalStateBranching(Energy scale,double z,
						   double phi, Energy pt) {
  ShoKinPtr showerKin = new_ptr(FS_QTildeShowerKinematics1to2());
  showerKin->scale(scale);
  showerKin->z(z);
  showerKin->phi(phi);
  showerKin->pT(pt);
  showerKin->SudakovFormFactor(this);
  return showerKin;
}

ShoKinPtr QTildeSudakov::createInitialStateBranching(Energy scale,double z,
						     double phi, Energy pt) {
  ShoKinPtr showerKin = new_ptr(IS_QTildeShowerKinematics1to2());
  showerKin->scale(scale);
  showerKin->z(z);
  showerKin->phi(phi);
  showerKin->pT(pt);
  showerKin->SudakovFormFactor(this);
  return showerKin;
}

ShoKinPtr QTildeSudakov::createDecayBranching(Energy scale,double z,
						     double phi, Energy pt) {
  ShoKinPtr  showerKin = new_ptr(Decay_QTildeShowerKinematics1to2());
  showerKin->scale(scale);
  showerKin->z(z);
  showerKin->phi(phi);
  showerKin->pT(pt);
  showerKin->SudakovFormFactor(this);
  return showerKin;
}
