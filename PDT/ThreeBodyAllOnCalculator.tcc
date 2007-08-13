// -*- C++ -*-
//
// This is the implementation of the non-inlined templated member
// functions of the ThreeBodyAllOnCalculator class.
//
namespace Herwig {
using namespace ThePEG;

// shift the variables for the outer integrand and give limits for the inner one
template <class T>
void ThreeBodyAllOnCalculator<T>::outerVariables(const double & x, Energy2 & low,
						 Energy2 & upp) const { 
  // first convert the value of x into the value of souter
  if(_channelmass[_thechannel] > 0*MeV) {
    _souter = _channelmass[_thechannel]*(_channelmass[_thechannel]+
					 _channelwidth[_thechannel]*tan(x));
  }
  else {
    _souter = UnitRemoval::E2 * pow(x,1./(_channelpower[_thechannel]+1.));
  }
  // now the limits of the inner integral
  Energy ea(0.*MeV),eb(0.*MeV),eam(0.*MeV),ebm(0.*MeV);
  Energy rs=sqrt(_souter);
  switch(_channeltype[_thechannel]) {
  case 1:
    ea = 0.5*(_souter-_m2[1]+_m2[2])/rs; eam=sqrt(ea*ea-_m2[2]);
    eb = 0.5*(_m2[0]-_souter-_m2[3])/rs; ebm=sqrt(eb*eb-_m2[3]);
    break;
  case 2:
    ea = 0.5*(_souter-_m2[1]+_m2[3])/rs; eam=sqrt(ea*ea-_m2[3]);
    eb = 0.5*(_m2[0]-_souter-_m2[2])/rs; ebm=sqrt(eb*eb-_m2[2]);
    break;
  case 3:
    ea = 0.5*(_souter-_m2[2]+_m2[3])/rs; eam=sqrt(ea*ea-_m2[3]);
    eb = 0.5*(_m2[0]-_souter-_m2[1])/rs; ebm=sqrt(eb*eb-_m2[1]);
    break;
  }
  Energy2 sum = sqr(ea+eb);
  // calculate the limits
  low = sum - sqr(eam+ebm);
  upp = sum - sqr(eam-ebm);
}

template <class T>
Energy2 ThreeBodyAllOnCalculator<T>::operator ()(Energy2 y) const {
  // set up the values of the s variables
  Energy2 s12(0.*MeV2),s23(0.*MeV2),s13(0.*MeV2),m2sum=_m2[0]+_m2[1]+_m2[2]+_m2[3];
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
  InvEnergy2 jacdem = 0./MeV2;
  Energy2 sjac(0.*MeV2); 
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
    InvEnergy2 term; 
    if(_channelmass[ix] > 0*MeV) {
      rm2 = sqr(_channelmass[ix]);
      rw2 = sqr(_channelwidth[ix]);
      Energy4 tmp = sqr(sjac-rm2) + rw2*rm2;
      term = _channelweights[ix]*_channelmass[ix]*_channelwidth[ix]/tmp;
    }
    else {
      term = UnitRemoval::InvE2 * _channelweights[ix]*(_channelpower[ix]+1.)*
	pow(sjac*UnitRemoval::InvE2, _channelpower[ix]);
    }
    jacdem += term;
  }
  // now computer the matrix element
  return _theME.threeBodyMatrixElement(_mode,_m2[0],s12,s13,
				       s23,_m[1],_m[2],_m[3])/jacdem;
}

// calculate the width for a given mass
template <class T>
Energy ThreeBodyAllOnCalculator<T>::partialWidth(Energy2 q2) const {
  Outer outer(this);
  _m[0] = sqrt(q2);
  _m2[0]=q2;
  // check the decay is kinematically allowed
  if(_m[0]<_m[1]+_m[2]+_m[3]){return 0.*MeV;}
  // perform the integrals for all the different channels
  Energy4 sum(0.*MeV2*MeV2);
  for(unsigned int ix=0,N=_channeltype.size(); ix<N; ++ix) {
    Energy2 upp(0.*MeV2),low(0.*MeV2);
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
    double rupp, rlow;
    // transform them
    if(_channelmass[ix] > 0*MeV) {
      rupp = atan((upp-_channelmass[ix]*_channelmass[ix])/
		  _channelmass[ix]/_channelwidth[ix]);
      rlow =  atan((low-_channelmass[ix]*_channelmass[ix])/
		   _channelmass[ix]/_channelwidth[ix]);
    }
    else {
      rupp = pow(upp*UnitRemoval::InvE2, _channelpower[ix]+1.);
      rlow = pow(low*UnitRemoval::InvE2, _channelpower[ix]+1.);
    }
    // perform the integral using my Gaussian quadature class
    _thechannel=ix;
    GaussianIntegrator intb;
    sum +=  _channelweights[ix] * intb.value(outer,rlow,rupp);
  }
  // final factors
  Energy3 fact = pow<3,1>(Constants::twopi * _m[0]);
  return sum/fact/32.;
}

}

  
