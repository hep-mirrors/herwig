// -*- C++ -*-
//
// QTildeSudakov.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
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
#include "Herwig++/Shower/Default/FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/Default/IS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/Default/Decay_QtildaShowerKinematics1to2.h"

using namespace Herwig;

ClassDescription<QTildeSudakov> QTildeSudakov::initQTildeSudakov;
// Definition of the static class description member.

void QTildeSudakov::persistentOutput(PersistentOStream & os) const {
  os << a_ << b_ << ounit(c_,GeV) << ounit(kinCutoffScale_,GeV) << cutOffOption_
     << ounit(vgcut_,GeV) << ounit(vqcut_,GeV) 
     << ounit(pTmin_,GeV) << ounit(pT2min_,GeV2);
}

void QTildeSudakov::persistentInput(PersistentIStream & is, int) {
  is >> a_ >> b_ >> iunit(c_,GeV) >> iunit(kinCutoffScale_,GeV) >> cutOffOption_
     >> iunit(vgcut_,GeV) >> iunit(vqcut_,GeV) 
     >> iunit(pTmin_,GeV) >> iunit(pT2min_,GeV2);
}

void QTildeSudakov::Init() {

  static ClassDocumentation<QTildeSudakov> documentation
    ("The QTildeSudakov class implements the Sudakov form factor for ordering it"
     " qtilde");

  static Parameter<QTildeSudakov,double> interfaceaParameter
    ("aParameter",
     "The a parameter for the kinematic cut-off",
     &QTildeSudakov::a_, 0.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<QTildeSudakov,double> interfacebParameter
    ("bParameter",
     "The b parameter for the kinematic cut-off",
     &QTildeSudakov::b_, 2.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<QTildeSudakov,Energy> interfacecParameter
    ("cParameter",
     "The c parameter for the kinematic cut-off",
     &QTildeSudakov::c_, GeV, 0.3*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<QTildeSudakov,Energy>
    interfaceKinScale ("cutoffKinScale",
		       "kinematic cutoff scale for the parton shower phase"
		       " space (unit [GeV])",
		       &QTildeSudakov::kinCutoffScale_, GeV, 
		       2.3*GeV, 0.001*GeV, 10.0*GeV,false,false,false);

  static Switch<QTildeSudakov,unsigned int> interfaceCutOffOption
    ("CutOffOption",
     "The type of cut-off to use to end the shower",
     &QTildeSudakov::cutOffOption_, 0, false, false);
  static SwitchOption interfaceCutOffOptionDefault
    (interfaceCutOffOption,
     "Default",
     "Use the standard Herwig++ cut-off on virtualities with the minimum"
     " virtuality depending on the mass of the branching particle",
     0);
  static SwitchOption interfaceCutOffOptionFORTRAN
    (interfaceCutOffOption,
     "FORTRAN",
     "Use a FORTRAN-like cut-off on virtualities",
     1);
  static SwitchOption interfaceCutOffOptionpT
    (interfaceCutOffOption,
     "pT",
     "Use a cut on the minimum allowed pT",
     2);
  
  static Parameter<QTildeSudakov,Energy> interfaceGluonVirtualityCut
    ("GluonVirtualityCut",
     "For the FORTRAN cut-off option the minimum virtuality of the gluon",
     &QTildeSudakov::vgcut_, GeV, 0.85*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<QTildeSudakov,Energy> interfaceQuarkVirtualityCut
    ("QuarkVirtualityCut",
     "For the FORTRAN cut-off option the minimum virtuality added to"
     " the mass for particles other than the gluon",
     &QTildeSudakov::vqcut_, GeV, 0.85*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<QTildeSudakov,Energy> interfacepTmin
    ("pTmin",
     "The minimum pT if using a cut-off on the pT",
     &QTildeSudakov::pTmin_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
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
  // compute the pt
  Energy2 pt2=sqr(z()*(1.-z()))*t-masssquared_[1]*(1.-z())-masssquared_[2]*z();
  if(ids_[0]!=ParticleID::g) pt2+=z()*(1.-z())*masssquared_[0];
  // if pt2<0 veto
  if(pt2<pT2min_) return true;
  // otherwise calculate pt and return
  pT(sqrt(pt2));
  return false;
}

 
ShoKinPtr QTildeSudakov::generateNextTimeBranching(const Energy startingScale,
						   const IdList &ids,const bool cc,
						   double enhance) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  q_ = Energy();
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
  if(t > Energy2()) q_ = sqrt(t);
  else q_ = -1.*MeV;
  phi(Constants::twopi*UseRandom::rnd());
  if(q_ < Energy()) return ShoKinPtr();
  // construct the ShowerKinematics object
  return constructKinematics(0);
}

ShoKinPtr QTildeSudakov::
generateNextSpaceBranching(const Energy startingQ,
			   const IdList &ids,
			   double x,bool cc,
			   double enhance,
			   Ptr<BeamParticleData>::transient_const_pointer beam) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  q_ = Energy();
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
  // calculate next value of t using veto algorithm
  Energy2 t(tmax),pt2(0.*GeV2);
  do {
    if(!guessSpaceLike(t,tmin,x,enhance)) break;
    pt2=sqr(1.-z())*t-z()*masssquared_[2];
  }
  while(z() > zLimits().second || 
	SplittingFnVeto((1.-z())*t/z(),ids,true) || 
	alphaSVeto(sqr(1.-z())*t) || 
	PDFVeto(t,x,parton0,parton1,beam) || pt2 < pT2min_ );
  if(t > Energy2() && zLimits().first < zLimits().second)  q_ = sqrt(t);
  else return ShoKinPtr();
  phi(Constants::twopi*UseRandom::rnd());
  pT(sqrt(pt2));
  // create the ShowerKinematics and return it
  return constructKinematics(1);
}

void QTildeSudakov::initialize(const IdList & ids, Energy2 & tmin,const bool cc) {
  ids_=ids;
  if(cc) {
    for(unsigned int ix=0;ix<ids.size();++ix) {
      if(getParticleData(ids[ix])->CC()) ids_[ix]*=-1;
    }
  }
  masses_.clear();
  masssquared_.clear();
  tmin=Energy2();
  if(cutOffOption_ == 0) {
    for(unsigned int ix=0;ix<ids_.size();++ix)
      masses_.push_back(getParticleData(ids_[ix])->mass());
    Energy kinCutoff=
      kinematicCutOff(kinScale(),*std::max_element(masses_.begin(),masses_.end()));
    for(unsigned int ix=0;ix<masses_.size();++ix)
      masses_[ix]=max(kinCutoff,masses_[ix]);
  }
  else if(cutOffOption_ == 1) {
    for(unsigned int ix=0;ix<ids_.size();++ix) {
      masses_.push_back(getParticleData(ids_[ix])->mass());
      masses_.back() += ids_[ix]==ParticleID::g ? vgcut_ : vqcut_;
    }
  }
  else if(cutOffOption_ == 2) {
    for(unsigned int ix=0;ix<ids_.size();++ix) 
      masses_.push_back(getParticleData(ids_[ix])->mass());
    tmin = 4.*pT2min_;
  }
  else {
    throw Exception() << "Unknown option for the cut-off"
		      << " in QTildeSudakov::initialize()"
		      << Exception::runerror;
  }
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
  Energy2 t(tmin),pt2;
  do {
    if(!guessDecay(t,tmax,minmass,enhance)) break;
    pt2 = sqr(1.-z())*(t-masssquared_[0])-z()*masssquared_[2];
  }
  while(SplittingFnVeto((1.-z())*t/z(),ids,true)|| 
	alphaSVeto(sqr(1.-z())*t) ||
	pt2<pT2min_ ||
	t*(1.-z())>masssquared_[0]-sqr(minmass));
  if(t > Energy2()) {
    q_ = sqrt(t);
    pT(sqrt(pt2));
  }
  else return ShoKinPtr();
  phi(Constants::twopi*UseRandom::rnd());
  // create the ShowerKinematics object
  return constructKinematics(2);
}

bool QTildeSudakov::guessDecay(Energy2 &t,Energy2 tmax, Energy minmass,
			       double enhance) {
  // previous scale
  Energy2 told = t;
  // overestimated limits on z
  pair<double,double> limits=make_pair(sqr(minmass/masses_[0]),
				       1.-masses_[2]/sqrt(tmax-masssquared_[0])
				       +0.5*masssquared_[2]/(tmax-masssquared_[0]));
  if(zLimits().second<zLimits().first) {
    t=-1.0*GeV2;
    return false;
  }
  zLimits(limits);
  // guess values of t and z
  t = guesst(told,2,ids_,enhance,ids_[1]==ids_[2]); 
  z(guessz(2,ids_)); 
  // actual values for z-limits
  limits=make_pair(sqr(minmass/masses_[0]),
		   1.-masses_[2]/sqrt(t-masssquared_[0])
		   +0.5*masssquared_[2]/(t-masssquared_[0]));
  zLimits(limits);
  if(t>tmax||zLimits().second<zLimits().first) {
    t=-1.0*GeV2;
    return false;
  }
  else
    return true; 
} 

bool QTildeSudakov::computeTimeLikeLimits(Energy2 & t) {
  // special case for gluon radiating
  pair<double,double> limits;
  if(ids_[0]==ParticleID::g) {
    // no emission possible
    if(t<16.*masssquared_[1]) {
      t=-1.*GeV2;
      return false;
    }
    // overestimate of the limits
    limits.first  = 0.5*(1.-sqrt(1.-4.*sqrt((masssquared_[1]+pT2min_)/t)));
    limits.second = 1.-limits.first;
  }
  // special case for radiated particle is gluon 
  else if(ids_[2]==ParticleID::g) {
    limits.first  =    sqrt((masssquared_[1]+pT2min_)/t);
    limits.second = 1.-sqrt((masssquared_[2]+pT2min_)/t);
  }
  else if(ids_[1]==ParticleID::g) {
    limits.second  =    sqrt((masssquared_[2]+pT2min_)/t);
    limits.first   = 1.-sqrt((masssquared_[1]+pT2min_)/t);
  }
  else {
    limits.first  =    (masssquared_[1]+pT2min_)/t;
    limits.second = 1.-(masssquared_[2]+pT2min_)/t; 
  }
  if(limits.first>=limits.second) {
    t=-1.*GeV2;
    return false;
  }
  zLimits(limits);
  return true;
}

bool QTildeSudakov::computeSpaceLikeLimits(Energy2 & t, double x) {
  pair<double,double> limits;
  // compute the limits
  limits.first = x;
  double yy = 1.+0.5*masssquared_[2]/t;
  limits.second = yy - sqrt(sqr(yy)-1.+pT2min_/t); 
  // return false if lower>upper
  zLimits(limits);
  if(limits.second<limits.first) {
    t=-1.*GeV2;
    return false;
  }
  else
    return true;
}

ShoKinPtr QTildeSudakov::constructKinematics(int iopt) {
  ShoKinPtr showerKin;
  switch(iopt) {
    // time-like
  case 0:
    showerKin = new_ptr(FS_QtildaShowerKinematics1to2());
    break;
    // space-like
  case 1:
    showerKin = new_ptr(IS_QtildaShowerKinematics1to2());
    break;
    // decay
  case 2:
    showerKin = new_ptr(Decay_QtildaShowerKinematics1to2());
    break;
  default:
    throw Exception() << "Unknown type of branching in "
		      << "QTildeSudakov::constructKinematics()"
		      << Exception::runerror;
  }
  showerKin->scale(q_);
  showerKin->z(z());
  showerKin->phi(phi());
  showerKin->pT(pT());
  showerKin->splittingFn(splittingFn());
  return showerKin;
}

void QTildeSudakov::doinit() throw(InitException) {
  SudakovFormFactor::doinit();
  pT2min_ = cutOffOption_==2 ? sqr(pTmin_) : 0.*GeV2; 
}
