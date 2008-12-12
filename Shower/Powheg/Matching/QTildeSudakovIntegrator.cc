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

  
  //massless veto approximations fo speed
  if( jetMeasureMode_ == 0 ){
    if( pt2 < sqr( max( z, 1. - z ) * mergeScale_ ) ) return 0.;
  }
  else {
    if( pt2 < sqr( mergeScale_ ) ) return 0.;
  }
  
  /*
  Energy2 kt_measure;

  //calculate mass dependent kt scale here
 
  if( jetMeasureMode_ == 0 ||  jetMeasureMode_ == 2 ){
    Energy2 s = s_;
    Energy pt = sqrt(pt2);

    Energy2 m0 = masssquared_[0];
    Energy2 m1 = masssquared_[1];
    Energy2 m2 = masssquared_[2];

    double lambda = sqrt( 1. - 4.*m0/s );
    double beta1 = 2.*( m1 - sqr(z)*m0 + sqr(pt) )
      / z / lambda / ( lambda + 1. ) / s;
    double beta2 = 2.*( m2 - sqr( 1. - z )*m0 + sqr(pt) )
      / ( 1. - z ) / lambda / ( lambda + 1. ) / s;

    Energy E1 = sqrt(s)/2.*( z + lambda*beta1 );
    Energy E2 = sqrt(s)/2.*( (1.-z) + lambda*beta2 );
    Energy Z1 = sqrt(s)/2.*lambda*( z - beta1 );
    Energy Z2 = sqrt(s)/2.*lambda*( (1.-z) - beta2 );;

    double costheta = ( Z1*Z2 - sqr(pt) )
      / sqrt( sqr(Z1)+sqr(pt) ) / sqrt( sqr(Z2)+sqr(pt) );

    //durham
    if( jetMeasureMode_ == 0 )
      kt_measure = 2.*min( sqr(E1), sqr(E2) )*( 1. - costheta );
    //luclus
    else if( jetMeasureMode_ == 2 )
      kt_measure = 2.*sqr(E1)*sqr(E2)/sqr(E1+E2)*( 1. - costheta );
  }
  else kt_measure = pt2;

  //return 0 if kt measure is less than the mergeScale_
  if( kt_measure < sqr( mergeScale_ ) ) return 0.;
  */  

  Energy2 t = z*(1.-z)*sqr(qtilde_);

  return 0.5/pi*coupling_->value(sqr(z*(1.-z)*qtilde_))*
    splittingFunction_->P(z,t,ids_,true);
}

QTildeSudakovIntegrator::QTildeSudakovIntegrator(const BranchingElement & br,
						 unsigned int jetMeasureMode,
						 Energy2 s) {
  jetMeasureMode_ = jetMeasureMode;
  // get the id list
  ids_ = br.second;
  s_ = s;
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

double QTildeSudakovIntegrator::value(Energy qtildemax, Energy qtildemin, Energy pt_cut) {
  // create the integrand for the inner integral
  
  inner_ = InnerSudakovIntegrand(this);
  qtildeh_ = qtildemax;
  mergeScale_ = pt_cut;
  // cerr<<"mergescale = "<< mergeScale_ / GeV <<"\t";
  // cerr<<"between "<< qtildemin / GeV <<" and "<< qtildemax / GeV <<"\t";

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
