// -*- C++ -*-
//
// ThreeBodyAllOnCalculator.tcc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the ThreeBodyAllOnCalculator class.
//
using namespace Herwig;

// shift the variables for the outer integrand and give limits for the inner one
template <class T>
void ThreeBodyAllOnCalculator<T>::outerVariables(const double & x, Energy2 & low,
						 Energy2 & upp) const { 
  // first convert the value of x into the value of souter
  if(_mapping[_thechannel]==0) {
    _souter = _channelmass[_thechannel]*(_channelmass[_thechannel]+
					 _channelwidth[_thechannel]*tan(x));
  }
  else if(_mapping[_thechannel]==1) {
    _souter = sqr(_channelmass[_thechannel])*(1.+1./x);
  }
  else {
    _souter = UnitRemoval::E2 * pow(x,1./(_channelpower[_thechannel]+1.));
  }
  // now the limits of the inner integral
  Energy ea(ZERO),eb(ZERO);
  Energy rs=sqrt(_souter);
  Energy2 eam2(ZERO),ebm2(ZERO);
  switch(_channeltype[_thechannel]) {
  case 1:
    ea   = 0.5*(_souter-_m2[1]+_m2[2])/rs; 
    eam2 = sqr(ea)-_m2[2];
    eb   = 0.5*(_m2[0]-_souter-_m2[3])/rs; 
    ebm2 = sqr(eb)-_m2[3];
    break;
  case 2:
    ea   = 0.5*(_souter-_m2[1]+_m2[3])/rs; 
    eam2 = sqr(ea)-_m2[3];
    eb   = 0.5*(_m2[0]-_souter-_m2[2])/rs; 
    ebm2 = sqr(eb)-_m2[2];
    break;
  case 3:
    ea   = 0.5*(_souter-_m2[2]+_m2[3])/rs; 
    eam2 = sqr(ea)-_m2[3];
    eb   = 0.5*(_m2[0]-_souter-_m2[1])/rs; 
    ebm2 = sqr(eb)-_m2[1];
    break;
  default:
    assert(false);
  }
  Energy eam = sqrt(max(ZERO,eam2));
  Energy ebm = sqrt(max(ZERO,ebm2));
  Energy2 sum = sqr(ea+eb);
  // calculate the limits
  low = sum - sqr(eam+ebm);
  upp = sum - sqr(eam-ebm);
}

template <class T>
Energy2 ThreeBodyAllOnCalculator<T>::operator ()(Energy2 y) const {
  assert(!isnan(y/MeV2)); 
  // set up the values of the s variables
  Energy2 s12(ZERO),s23(ZERO),s13(ZERO),
    m2sum(_m2[0]+_m2[1]+_m2[2]+_m2[3]);
  switch(_channeltype[_thechannel]) {
  case 1:
    s12 = _souter;
    s23 = y;
    s13 = m2sum-s12-s23;
    break;
  case 2:
    s23 = y;
    s13 = _souter;
    s12 = m2sum-s23-s13;
    break;
  case 3:
    s23 = _souter;
    s13 = y;
    s12 = m2sum-s23-s13;
    break;
  }
  // compute the jacobian
  // computer the denominator for the jacobian
  InvEnergy2 jacdem = ZERO;
  Energy2 sjac(ZERO); 
  Energy2 rm2,rw2;
  for(unsigned int ix=0,N=_channeltype.size(); ix<N; ++ix) {
    switch(_channeltype[ix]) {
    case 1:
      sjac = s12;
      break;
    case 2:
      sjac = s13;
      break;
    case 3:
      sjac = s23;
      break;
    }
    assert(!isnan(sjac/MeV2));
    InvEnergy2 term; 

    if(_mapping[ix]==0) {
      rm2 = sqr(_channelmass[ix]);
      rw2 = sqr(_channelwidth[ix]);
      Energy4 tmp = sqr(sjac-rm2) + rw2*rm2;
      term = _channelweights[ix]*_channelmass[ix]*_channelwidth[ix]/tmp;
    }
    else if(_mapping[ix]==1) {
      term = _channelweights[ix]*sqr(_channelmass[ix]/(sjac-sqr(_channelmass[ix])));
    }
    else if(_mapping[ix]==2) {
      term = UnitRemoval::InvE2 * _channelweights[ix]*(_channelpower[ix]+1.)*
     	pow(sjac*UnitRemoval::InvE2, _channelpower[ix]);
    }
    else
      assert(false);
    jacdem += term;
  }
  // now computer the matrix element
  return _theME.threeBodyMatrixElement(_mode,_m2[0],s12,s13,
				       s23,_m[1],_m[2],_m[3])/jacdem;
}

// calculate the width for a given mass
template <class T>
Energy ThreeBodyAllOnCalculator<T>::partialWidth(Energy2 q2) const {
  Outer outer(this,_relerr);
  _m[0] = sqrt(q2);
  _m2[0]=q2;
  // check the decay is kinematically allowed
  if(_m[0]<_m[1]+_m[2]+_m[3]) return ZERO;
  // set up for the different channels
  unsigned int N = _channeltype.size();
  vector<double> rupp(N,0.),rlow(N,0.);
  for(unsigned int ix=0; ix<N; ++ix) {
    Energy2 upp(ZERO),low(ZERO);
    // work out the kinematic limits
    switch(_channeltype[ix]) {
    case 1:
      upp = sqr(_m[0]-_m[3]);
      low = sqr(_m[1]+_m[2]);
      break;
    case 2:
      upp = sqr(_m[0]-_m[2]);
      low = sqr(_m[1]+_m[3]);
      break;
    case 3:
      upp = sqr(_m[0]-_m[1]);
      low = sqr(_m[2]+_m[3]);
      break;
    default:
      assert(false);
    }
    // transform them
    if(_channelmass[ix] > ZERO) {
      if(_channelwidth[ix] > 1e-8*MeV) {
	rupp[ix] = atan2((upp-_channelmass[ix]*_channelmass[ix]),
			 _channelmass[ix]*_channelwidth[ix]);
	rlow[ix] =  atan2((low-_channelmass[ix]*_channelmass[ix]),
			  _channelmass[ix]*_channelwidth[ix]);
	_mapping[ix] = 0;
	if(rupp[ix]/rlow[ix]>0.&&_channelwidth[ix]/_channelmass[ix]<1e-6) {
	  _mapping[ix] = 1;
	  Energy2 m2=sqr(_channelmass[ix]);
	  rupp[ix] = m2/(low-m2);
	  rlow[ix] = m2/(upp-m2);
	}
      }
      else {
	_mapping[ix] = 1;
	Energy2 m2=sqr(_channelmass[ix]);
	rupp[ix] = m2/(low-m2);
	rlow[ix] = m2/(upp-m2);
      }
    }
    else {
      _mapping[ix] = 2;
      rupp[ix] = pow(upp*UnitRemoval::InvE2, _channelpower[ix]+1.);
      rlow[ix] = pow(low*UnitRemoval::InvE2, _channelpower[ix]+1.);
    }
  }
  // perform the integrals for all the different channels
  Energy4 sum(ZERO);
  for(unsigned int ix=0,N=_channeltype.size(); ix<N; ++ix) {
    // perform the integral using GSLIntegrator class
    _thechannel=ix;
    GSLIntegrator intb(1e-35,_relerr,1000);
    sum +=  _channelweights[ix] * intb.value(outer,rlow[ix],rupp[ix]);
  }
  // final factors
  Energy3 fact = pow<3,1>(Constants::twopi * _m[0]);
  return sum/fact/32.;
}
