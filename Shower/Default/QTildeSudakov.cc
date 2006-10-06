// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeSudakov class.
//

#include "QTildeSudakov.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/Default/FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/Default/IS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Shower/Default/Decay_QtildaShowerKinematics1to2.h"

using namespace Herwig;

ClassDescription<QTildeSudakov> QTildeSudakov::initQTildeSudakov;
// Definition of the static class description member.

void QTildeSudakov::persistentOutput(PersistentOStream & os) const {
  os << _a << _b << _c << _kinCutoffScale << _cutoffQCDMassScale
     << _cutoffQEDMassScale << _cutoffEWKMassScale;
}

void QTildeSudakov::persistentInput(PersistentIStream & is, int) {
  is >> _a >> _b >> _c >> _kinCutoffScale >> _cutoffQCDMassScale
     >> _cutoffQEDMassScale >> _cutoffEWKMassScale;
}

void QTildeSudakov::Init() {

  static ClassDocumentation<QTildeSudakov> documentation
    ("The QTildeSudakov class implements the Sudakov form factor for ordering it"
     " qtilde");

  static Parameter<QTildeSudakov,double> interfaceaParameter
    ("aParameter",
     "The a parameter for the kinematic cut-off",
     &QTildeSudakov::_a, 0.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<QTildeSudakov,double> interfacebParameter
    ("bParameter",
     "The b parameter for the kinematic cut-off",
     &QTildeSudakov::_b, 2.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<QTildeSudakov,Energy> interfacecParameter
    ("cParameter",
     "The c parameter for the kinematic cut-off",
     &QTildeSudakov::_c, GeV, 0.3*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<QTildeSudakov,Energy>
    interfaceKinScale ("cutoffKinScale",
		       "kinematic cutoff scale for the parton shower phase"
		       " space (unit [GeV])",
		       &QTildeSudakov::_kinCutoffScale, GeV, 
		       0.75*GeV, 0.001*GeV, 10.0*GeV,false,false,false);

  static Parameter<QTildeSudakov,Energy>
    interfaceCutoffQCD ("CutoffQCDMassScale",
			"low energy cutoff mass scale for QCD radiation  (unit [GeV])",
			&QTildeSudakov::_cutoffQCDMassScale, GeV, 
			0.0*GeV, 0.0*GeV, 10.0*GeV,false,false,false);
  static Parameter<QTildeSudakov,Energy>
    interfaceCutoffQED ("CutoffQEDMassScale",
			"low energy cutoff mass scale for QED radiation  (unit [GeV])",
			&QTildeSudakov::_cutoffQEDMassScale, GeV, 
			0.0005*GeV, 0.0*GeV, 10.0*GeV,false,false,false);
  static Parameter<QTildeSudakov,Energy>
    interfaceCutoffEWK ("CutoffEWKMassScale",
			"low energy cutoff mass scale for EWK radiation  (unit [GeV])",
			&QTildeSudakov::_cutoffEWKMassScale, GeV, 
			91.0*GeV, 0.0*GeV, 1000.0*GeV,false,false,false);

}

bool QTildeSudakov::guessTimeLike(Energy2 &t,Energy2 tmin,double enhance) {
  Energy2 told = t;
  // calculate limits on z and if lower>upper return
  if(!computeTimeLikeLimits(t)) return false;
  // guess values of t and z
  t = guesst(told,0,enhance,_ids[1]==_ids[2]); 
  z(guessz()); 
  // actual values for z-limits
  if(!computeTimeLikeLimits(t)) return false;
  if(t<tmin) {
    t=-1.0*GeV;
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
  t = guesst(told,1,enhance,_ids[1]==_ids[2]); 
  z(guessz()); 
  // actual values for z-limits
  if(!computeSpaceLikeLimits(t,x)) return false;
  if(t<tmin) {
    t=-1.0*GeV;
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
  Energy2 pt2=sqr(z()*(1.-z()))*t-_masssquared[1]*(1.-z())-_masssquared[2]*z();
  if(_ids[0]!=ParticleID::g) pt2+=z()*(1.-z())*_masssquared[0];
  // if pt2<0 veto
  if(pt2<0.) return true;
  // otherwise calculate pt and return
  pT(sqrt(pt2));
  return false;
}

 
ShoKinPtr QTildeSudakov::generateNextTimeBranching(const Energy startingScale,
						   const IdList &ids,const bool cc,
						   double enhance) {
  if(startingScale<=cutoffQScale(interactionType())) ShoKinPtr();
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  _q = Energy();
  z(0.);
  phi(0.); 
  // perform initialization
  Energy2 tmax(sqr(startingScale)),tmin;
  initialize(ids,tmin,cc);
  // check max > min
  if(tmax<=tmin) return ShoKinPtr();
  // calculate next value of t using veto algorithm
  Energy2 t(tmax);
  do  
    if(!guessTimeLike(t,tmin,enhance)) break;
  while(PSVeto(t) || SplittingFnVeto(z()*(1.-z())*t,ids,true) || 
	alphaSVeto(sqr(z()*(1.-z()))*t));
  if(t > 0) _q = sqrt(t);
  else _q = -1.;
  phi(2.*pi*UseRandom::rnd());
  if(_q<0) return ShoKinPtr();
  // construct the ShowerKinematics object
  FS_QtildaShowerKinematics1to2Ptr showerKin = 
    new_ptr(FS_QtildaShowerKinematics1to2());
  showerKin->scale(_q);
  showerKin->z(z());
  showerKin->phi(phi());
  showerKin->pT(pT());
  showerKin->splittingFn(splittingFn());
  return showerKin;
}

ShoKinPtr QTildeSudakov::
generateNextSpaceBranching(const Energy startingQ,
			   const IdList &ids,
			   double x,bool cc,
			   double enhance,
			   Ptr<BeamParticleData>::transient_const_pointer beam) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  _q = Energy();
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
  Energy2 t(tmax),pt2(0.);
  do {
    if(!guessSpaceLike(t,tmin,x,enhance)) break;
    pt2=sqr(1.-z())*t-z()*sqr(_kinCutoff);
  }
  while(z() > zLimits().second || 
	SplittingFnVeto((1.-z())*t/z(),ids,true) || 
	alphaSVeto(sqr(1.-z())*t) || 
	PDFVeto(t,x,parton0,parton1,beam)||pt2<0);
  if(t > 0 && zLimits().first < zLimits().second)  _q = sqrt(t);
  else return ShoKinPtr();
  phi(2.*pi*UseRandom::rnd());
  pT(sqrt(pt2));
  // create the ShowerKinematics and return it
  IS_QtildaShowerKinematics1to2Ptr 
    showerKin = new_ptr(IS_QtildaShowerKinematics1to2());
  showerKin->scale(_q);
  showerKin->z(z());
  showerKin->phi(phi());
  showerKin->pT(pT());
  showerKin->splittingFn(splittingFn());
  return showerKin;
}

void QTildeSudakov::initialize(const IdList & ids, Energy2 & tmin,const bool cc)
{
  _ids=ids;
  if(cc) {
    for(unsigned int ix=0;ix<ids.size();++ix) {
      if(getParticleData(ids[ix])->CC()) _ids[ix]*=-1;
    }
  }
  _masses.clear();
  _masssquared.clear();
  tmin=0.;
  unsigned int ix;
  for(ix=0;ix<_ids.size();++ix)
    _masses.push_back(getParticleData(_ids[ix])->mass());
  _kinCutoff=
    kinematicCutOff(kinScale(),*std::max_element(_masses.begin(),_masses.end()));
  for(ix=0;ix<_masses.size();++ix)
    {
      _masses[ix]=max(_kinCutoff,_masses[ix]);
      _masssquared.push_back(sqr(_masses[ix]));
      if(ix>0) tmin=max(_masssquared[ix],tmin);
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
  _q = Constants::MaxEnergy;
  z(0.);
  phi(0.); 
  // perform initialisation
  Energy2 tmax(sqr(stoppingScale)),tmin;
  initialize(ids,tmin,cc);
  tmin=sqr(startingScale);
  // check some branching possible
  if(tmax<=tmin) return ShoKinPtr();
  // perform the evolution
  Energy2 t(tmin);
  do 
    if(!guessDecay(t,tmax,minmass,enhance)) break;
  while(SplittingFnVeto((1.-z())*t/z(),ids,true)|| 
	alphaSVeto(sqr(1.-z())*t)||
	sqr(1.-z())*(t-_masssquared[0])<z()*sqr(_kinCutoff)||
	t*(1.-z())>_masssquared[0]-sqr(minmass));
  if(t > 0) {
    _q = sqrt(t);
    pT(sqrt(sqr(1.-z())*(t-_masssquared[0])-z()*sqr(_kinCutoff)));
  }
  else return ShoKinPtr();
  phi(2.*pi*UseRandom::rnd());
  // create the ShowerKinematics object
  Decay_QtildaShowerKinematics1to2Ptr showerKin = 
    new_ptr(Decay_QtildaShowerKinematics1to2());
  showerKin->scale(_q);
  showerKin->z(z());
  showerKin->phi(phi());
  showerKin->pT(pT());
  showerKin->splittingFn(splittingFn());
  return showerKin;
}

bool QTildeSudakov::guessDecay(Energy2 &t,Energy2 tmax, Energy minmass,
			       double enhance) {
  // previous scale
  Energy2 told = t;
  // overestimated limits on z
  pair<double,double> limits=make_pair(sqr(minmass/_masses[0]),
				       1.-_kinCutoff/sqrt(tmax-_masssquared[0])
				       +0.5*sqr(_kinCutoff)/(tmax-_masssquared[0]));
  if(zLimits().second<zLimits().first) {
    t=-1.0*GeV;
    return false;
  }
  zLimits(limits);
  // guess values of t and z
  t = guesst(told,2,enhance,_ids[1]==_ids[2]); 
  z(guessz()); 
  // actual values for z-limits
  limits=make_pair(sqr(minmass/_masses[0]),
		   1.-_kinCutoff/sqrt(t-_masssquared[0])
		   +0.5*sqr(_kinCutoff)/(t-_masssquared[0]));
  zLimits(limits);
  if(t>tmax||zLimits().second<zLimits().first) {
    t=-1.0*GeV;
    return false;
  }
  else
    return true; 
} 

bool QTildeSudakov::computeTimeLikeLimits(Energy2 & t)
{
  // special case for gluon radiating
  pair<double,double> limits;
  if(_ids[0]==ParticleID::g) {
    // no emission possible
    if(t<16.*_masssquared[1]) {
      t=-1.;
      return false;
    }
    // overestimate of the limits
    limits.first  = 0.5*(1.-sqrt(1.-4.*sqrt(_masssquared[1]/t)));
    limits.second = 1.-limits.first;
  }
  // special case for radiated particle is gluon 
  else if(_ids[2]==ParticleID::g) {
    limits.first  = 0.5*sqrt(_masssquared[1]/t);
    limits.second = 1.-0.5*sqrt(_masssquared[2]/t);
  }
  else if(_ids[1]==ParticleID::g) {
    limits.second  = 0.5*sqrt(_masssquared[2]/t);
    limits.first   = 1.-0.5*sqrt(_masssquared[1]/t);
  }
  else
    {throw Exception() << "QTildeSudakov::computeTimeLikeLimits() " 
			<< "general case not implemented " << Exception::runerror;}
  zLimits(limits);
  return true;
}

bool QTildeSudakov::computeSpaceLikeLimits(Energy2 & t, double x)
{
  pair<double,double> limits;
  // compute the limits
  limits.first = x;
  double yy = 1.+0.5*sqr(_kinCutoff)/t;
  limits.second = yy - sqrt(sqr(yy)-1.); 
  // return false if lower>upper
  zLimits(limits);
  if(limits.second<limits.first) {
    t=-1.*GeV;
    return false;
  }
  else
    return true;
}
