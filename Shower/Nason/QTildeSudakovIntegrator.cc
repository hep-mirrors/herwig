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
  Energy2 pt2=sqr(z*(1.-z))*sqr(_qtilde)-_masssquared[1]*(1.-z)-_masssquared[2]*z;
  if(_ids[0]!=ParticleID::g) pt2+=z*(1.-z)*_masssquared[0];
  // if pt2<0 veto
  if(pt2<Energy2()) return 0.;
  // return the answer
  Energy2 t = z*(1.-z)*sqr(_qtilde);
  return 0.5/pi*_coupling->value(sqr(z*(1.-z)*_qtilde))*
    _splittingFunction->P(z,t,_ids,true);
}

QTildeSudakovIntegrator::QTildeSudakovIntegrator(const BranchingElement & br) {
  // get the id list
  _ids = br.second;
  // get the coupling
  _coupling = br.first->alpha();
  // get the splitting function
  _splittingFunction = br.first->splittingFn();
  // get the parameters we'll need for the integration limits etc
  Energy2 tmin=Energy2();
  for(unsigned int ix=0;ix<_ids.size();++ix)
    _masses.push_back(CurrentGenerator::current().getParticleData(_ids[ix])->mass());
  _kinCutoff=
    br.first->kinematicCutOff(br.first->kinScale(),
			      *std::max_element(_masses.begin(),
						_masses.end()));
  for(unsigned int ix=0;ix<_masses.size();++ix) {
    _masses[ix]=max(_kinCutoff,_masses[ix]);
    _masssquared.push_back(sqr(_masses[ix]));
    if(ix>0) tmin=max(_masssquared[ix],tmin);
  }
  _qtildemin=sqrt(tmin);
}

double QTildeSudakovIntegrator::value(Energy qtildemax, Energy qtildemin) {
  // create the integrand for the inner integral
  _inner = InnerSudakovIntegrand(this);
  _qtildeh = qtildemax;
  double rhomax(0.),rhomin(2.*log(qtildemin/qtildemax));
  return _outerIntegrator.value(*this,rhomin,rhomax);
}

double QTildeSudakovIntegrator::operator() (double rho) const {
  Energy2 qtilde2 = exp(rho)*sqr(_qtildeh);
  _qtilde = sqrt(qtilde2);
  double zmin(0.),zmax(1.);
  if(_ids[0]==ParticleID::g) {
    // no emission possible
    if(qtilde2<16.*_masssquared[1]) return 0.;
    // overestimate of the limits
    zmin = 0.5*(1.-sqrt(1.-4.*sqrt(_masssquared[1]/qtilde2)));
    zmax = 1.-zmin;
  }
  // special case for radiated particle is gluon 
  else if(_ids[2]==ParticleID::g) {
    zmin = 0.5*sqrt(_masssquared[1]/qtilde2);
    zmax = 1.-0.5*sqrt(_masssquared[2]/qtilde2);
  }
  else if(_ids[1]==ParticleID::g) {
    zmax = 0.5*sqrt(_masssquared[2]/qtilde2);
    zmin = 1.-0.5*sqrt(_masssquared[1]/qtilde2);
  }
  else {
    zmin = _masssquared[1]/qtilde2;
    zmax = 1.-_masssquared[2]/qtilde2; 
  }
  if(zmax<=zmin) return 0.;
  return _innerIntegrator.value(_inner,zmin,zmax);
}
