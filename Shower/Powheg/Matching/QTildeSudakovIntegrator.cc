// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeSudakovIntegrator class.
//

#include "QTildeSudakovIntegrator.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

double QTildeSudakovIntegrator::innerIntegrand(double z) const {
  using Constants::pi;
  // compute the pt
  Energy2 pt2=sqr(z*(1.-z))*sqr(qtilde_)-masssquared_[1]*(1.-z)-masssquared_[2]*z;
  if(ids_[0]!=ParticleID::g) pt2+=z*(1.-z)*masssquared_[0];

  if( jetMeasureMode_ == 0 ){
    if( pt2 < sqr( max( z, 1. - z ) * mergeScale_ ) ) return 0.;
  }
  else {
    if( pt2 < sqr( mergeScale_ ) ) return 0.;
  }

  // return the answer
  Energy2 t = z*(1.-z)*sqr(qtilde_);
  return 0.5/pi*coupling_->value(sqr(z*(1.-z)*qtilde_))*
    splittingFunction_->P(z,t,ids_,true);
}

QTildeSudakovIntegrator::QTildeSudakovIntegrator(const BranchingElement & br,
						 Energy MergeScale,
						 unsigned int jetMeasureMode) {
  mergeScale_ = MergeScale; 
  jetMeasureMode_ = jetMeasureMode;
  // get the id list
  ids_ = br.second;
  // get the coupling
  coupling_ = br.first->alpha();
  // get the splitting function
  splittingFunction_ = br.first->splittingFn();
  // get the parameters we'll need for the integration limits etc
  // minimum pT
  pTmin_ = sqrt(br.first->pT2min());
  // masses etc
  Energy2 tmin=Energy2();
  Energy kinCutoff;
  switch (br.first->cutOffOption()) {
  case 0:
    for(unsigned int ix=0;ix<ids_.size();++ix)
      masses_.push_back(CurrentGenerator::current().getParticleData(ids_[ix])->mass());
    kinCutoff = br.first->kinematicCutOff(br.first->kinScale(),
					  *std::max_element(masses_.begin(),
							    masses_.end()));
    for(unsigned int ix=0;ix<masses_.size();++ix)
      masses_[ix]=max(kinCutoff,masses_[ix]);
    break;
  case 1:
    for(unsigned int ix=0;ix<ids_.size();++ix) {
      masses_.push_back(CurrentGenerator::current().getParticleData(ids_[ix])->mass());
      masses_.back() += ids_[ix]==ParticleID::g ? br.first->vgCut() : br.first->vqCut();
    }
    break;
  case 2:
    for(unsigned int ix=0;ix<ids_.size();++ix) 
      masses_.push_back(CurrentGenerator::current().getParticleData(ids_[ix])->mass());
    tmin = 4.*br.first->pT2min();
    break;
  default:
    throw Exception() << "Unknown option for the cut-off"
		      << " in QTildeSudakovIntegrator::QTildeSudakovIntegrator()"
		      << Exception::runerror;
  }
  // squares of masses and minimum scale
  for(unsigned int ix=0;ix<masses_.size();++ix) {
    masssquared_.push_back(sqr(masses_[ix]));
    if(ix>0) tmin=max(masssquared_[ix],tmin);
  }
  qtildemin_=sqrt(tmin);
}

double QTildeSudakovIntegrator::value(Energy qtildemax, Energy qtildemin) {
  // create the integrand for the inner integral
  inner_ = InnerSudakovIntegrand(this);
  qtildeh_ = qtildemax;
  double rhomax(0.),rhomin(2.*log(qtildemin/qtildemax));
  return outerIntegrator_.value(*this,rhomin,rhomax);
}

double QTildeSudakovIntegrator::operator() (double rho) const {
  Energy2 qtilde2 = exp(rho)*sqr(qtildeh_);
  qtilde_ = sqrt(qtilde2);
  double zmin(0.),zmax(1.);
  if(ids_[0]==ParticleID::g) {
    // no emission possible
    if(qtilde2<16.*masssquared_[1]) return 0.;
    // overestimate of the limits
    zmin = 0.5*(1.-sqrt(1.-4.*sqrt(masssquared_[1]/qtilde2)));
    zmax = 1.-zmin;
  }
  // special case for radiated particle is gluon 
  else if(ids_[2]==ParticleID::g) {
    zmin = 0.5*sqrt(masssquared_[1]/qtilde2);
    zmax = 1.-0.5*sqrt(masssquared_[2]/qtilde2);
  }
  else if(ids_[1]==ParticleID::g) {
    zmax = 0.5*sqrt(masssquared_[2]/qtilde2);
    zmin = 1.-0.5*sqrt(masssquared_[1]/qtilde2);
  }
  else {
    zmin = masssquared_[1]/qtilde2;
    zmax = 1.-masssquared_[2]/qtilde2; 
  }
  if(zmax<=zmin) return 0.;
  return innerIntegrator_.value(inner_,zmin,zmax);
}
