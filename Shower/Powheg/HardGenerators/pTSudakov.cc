// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the pTSudakov class.
//

#include "pTSudakov.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Default/FS_QtildaShowerKinematics1to2.h"

using namespace Herwig;

void pTSudakov::persistentOutput(PersistentOStream & ) const {
}

void pTSudakov::persistentInput(PersistentIStream & , int) {
}

ClassDescription<pTSudakov> pTSudakov::initpTSudakov;
// Definition of the static class description member.

void pTSudakov::Init() {

  static ClassDocumentation<pTSudakov> documentation
    ("There is no documentation for the pTSudakov class");

}

ShoKinPtr pTSudakov::generateNextTimeBranching(const Energy startingScale,
					       const IdList &ids,const bool cc,
					       double) {
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to the method.
  _q = Energy();
  z(0.);
  phi(0.); 
  // perform initialization
  Energy2 tmax(sqr(startingScale)),tmin;
  initialize(ids,tmin,tmax,cc);
  // check max > min
  if(tmax<=tmin) return ShoKinPtr();
  // calculate next value of t using veto algorithm
  Energy2 t(tmax);
  do {
    if(!guessTimeLike(t,tmin)) break;
  }
  while(_pt2<0.*MeV2
	||_q2>_q2max||SplittingFnVeto(t/(z()*(1.-z())),ids,true) || 
	alphaSVeto(t));
  if(t > 0*MeV2) _q = sqrt(t);
  else _q = -1.*MeV;
  phi(Constants::twopi*UseRandom::rnd());
  if(_q<0*MeV) return ShoKinPtr();
  // construct the ShowerKinematics object
  FS_QtildaShowerKinematics1to2Ptr showerKin = 
    new_ptr(FS_QtildaShowerKinematics1to2());
  pT(sqrt(_pt2));
  showerKin->scale(_q);
  showerKin->z(z());
  showerKin->phi(phi());
  showerKin->pT(pT());
  showerKin->splittingFn(splittingFn());
  return showerKin;
}

ShoKinPtr pTSudakov::generateNextDecayBranching(const Energy,
						const Energy,
						const Energy,
						const IdList &,
						const bool,
						double) {
  throw Exception() << "pTSudakov::generateNextDecayBranching not yet implemented"
		    << Exception::abortnow;
}

ShoKinPtr pTSudakov::generateNextSpaceBranching(const Energy,
						const IdList &, double,
						const bool, double,
						Ptr<BeamParticleData>::transient_const_pointer) {
  throw Exception() << "pTSudakov::generateNextSpaceBranching not yet implemented"
		    << Exception::abortnow;
}

void pTSudakov::initialize(const IdList & ids, Energy2 & tmin,Energy2 tmax,
			   const bool cc) {
  _ids=ids;
  _masses.clear();
  _masssquared.clear();
  tmin=0.*MeV2;
  unsigned int ix;
  for(ix=0;ix<_ids.size();++ix) {
    tcPDPtr part=getParticleData(_ids[ix]);
    _masses.push_back(part->mass());
    if(cc&&part->CC()) _ids[ix]*=-1;
  }
  _kinCutoff=
    kinematicCutOff(kinScale(),*std::max_element(_masses.begin(),_masses.end()));
  for(ix=0;ix<_masses.size();++ix) {
    _masses[ix]=max(_kinCutoff,_masses[ix]);
    _masssquared.push_back(sqr(_masses[ix]));
  }
  Energy2 mdiff(abs(_masssquared[1]+_masssquared[2]));
  double lambda = (1.+sqr(_masssquared[1]/_q2max)+sqr(_masssquared[2]/_q2max)
		   -2.*(_masssquared[1]+_masssquared[2])/_q2max
		   +2.*_masssquared[1]*_masssquared[2]/sqr(_q2max));
  tmin = (1.-_masssquared[0]/_q2max)*(_masssquared[1]+_masssquared[2]
 				      -sqr(mdiff)/_q2max
 				      +mdiff*lambda);
  //cerr << "testing " << sqrt(tmin) << "\n";
  tmin=0.2*GeV*MeV;
}

bool pTSudakov::guessTimeLike(Energy2 &t,Energy2 tmin) {
  Energy2 told = t;
  // calculate limits on z and if lower>upper return
  if(tmin>0.25*(_q2max-_masssquared[0])) return false;
  double root = sqrt(1.-4.*tmin/(_q2max-_masssquared[0]));
  pair<double,double> limits(0.5*(1.-root),0.5*(1.+root));
  if(limits.first>limits.second) return false;
  zLimits(limits);
  // guess values of t and z
  t = guesst(told,0,_ids,1.,_ids[1]==_ids[2]); 
  z(guessz(0,_ids));
  // actual values for z-limits
  _pt2 = t -_masssquared[1]*(1.-z())-_masssquared[2]*z();
  _q2  = _pt2/(z()*(1.-z()))+_masssquared[1]/z()+_masssquared[2]/(1.-z());
  if(t<tmin) {
    t=-1.0*MeV2;
    return false;
  }
  else
    return true; 
} 

Energy pTSudakov::calculateScale(double z, Energy pt, IdList ids,
					 unsigned int iopt) {
  throw Exception() << "Base class SudakovFormFactor::calculateScale() called "
		    << "this should be overidden in the inheriting class"
		    << Exception::runerror;
}

ShoKinPtr pTSudakov::createFinalStateBranching(Energy scale,double z,
				    double phi, Energy pt) {
  throw Exception() << " pTSudakov::createFinalStateBranching() not implemented"
		    << Exception::runerror;
}

ShoKinPtr pTSudakov::createInitialStateBranching(Energy scale,double z,
						 double phi, Energy pt) {
  throw Exception() << " pTSudakov::createInitialStateBranching() not implemented"
		    << Exception::runerror;
}
