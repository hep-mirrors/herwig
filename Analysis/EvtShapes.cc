#include "EvtShapes.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

using namespace CLHEP;


using namespace Herwig;
using namespace ThePEG;

// live and die

EvtShapes::EvtShapes() {}
EvtShapes::~EvtShapes() {}

EvtShapes::EvtShapes(const tPVector& part) {
  _pv = part; 
  _linTenDone = false; 
  _spherDone = false; 
  _thrustDone = false; 
  _hemDone = false; 
  _broadDone = false; 
  _useCmBoost = false; 
}

// some convenient single particle variables
double EvtShapes::getX(const Lorentz5Momentum & p, 
			      const Energy & Ebeam) {
  return(Ebeam > Energy() ? double(p.vect().mag()/Ebeam) : -1.); 
}

double EvtShapes::getXi(const Lorentz5Momentum & p, 
			       const Energy & Ebeam) {
  return((Ebeam > Energy() && p.vect().mag() > Energy()) ? 
	 log(Ebeam/p.vect().mag()) : -1.); 
}

Energy EvtShapes::getPt(const Lorentz5Momentum & p) {
  return p.perp(); 
}

double EvtShapes::getRapidity(const Lorentz5Momentum & p) {
  return (p.t() > p.z() ? p.rapidity() : 1e99); 
}

// single particle variables related to shape axis
Energy EvtShapes::ptInT(const Lorentz5Momentum & p) {
  checkThrust(); 
  return p.vect()*_thrustAxis[1]; 
}

Energy EvtShapes::ptOutT(const Lorentz5Momentum & p) {
  checkThrust(); 
  return p.vect()*_thrustAxis[2]; 
}

double EvtShapes::yT(const Lorentz5Momentum & p) {
  checkThrust(); 
  return (p.t() > p.vect()*_thrustAxis[0] ? 
	  p.rapidity(_thrustAxis[0]) : 1e99);
}

Energy EvtShapes::ptInS(const Lorentz5Momentum & p) { 
  checkSphericity(); 
  return p.vect()*_spherAxis[1]; 
}

Energy EvtShapes::ptOutS(const Lorentz5Momentum & p) {
  checkSphericity(); 
  return p.vect()*_spherAxis[2]; 
}

double EvtShapes::yS(const Lorentz5Momentum & p) {
  checkSphericity(); 
  return (p.t() > p.vect()*_spherAxis[0] ? 
	  p.rapidity(_spherAxis[0]) : 1e99);
}

// thrust related 
double EvtShapes::thrust() {
  checkThrust(); 
  return _thrust[0];
}

double EvtShapes::thrustMajor() {
  checkThrust(); 
  return _thrust[1];
}

double EvtShapes::thrustMinor() {
  checkThrust(); 
  return _thrust[2];
}

double EvtShapes::oblateness() {
  checkThrust(); 
  return _thrust[1]-_thrust[2];
}

Axis EvtShapes::thrustAxis() {
  checkThrust(); 
  return _thrustAxis[0];
}

Axis EvtShapes::majorAxis() {
  checkThrust(); 
  return _thrustAxis[1];
}

Axis EvtShapes::minorAxis() {
  checkThrust(); 
  return _thrustAxis[2];
}

// linear tensor related
vector<double> EvtShapes::linTenEigenValues() {
  checkLinTen(); 
  return _linTen; 
}

vector<Axis> EvtShapes::linTenEigenVectors() {
  checkLinTen(); 
  return _linTenAxis; 
}

double EvtShapes::CParameter() {
  checkLinTen(); 
  return 3.*(_linTen[0]*_linTen[1]+_linTen[1]*_linTen[2]
	     +_linTen[2]*_linTen[0]); 
}

double EvtShapes::DParameter() {
  checkLinTen(); 
  return 27.*(_linTen[0]*_linTen[1]*_linTen[2]); 
}

// quadratic tensor related
double EvtShapes::sphericity() {
  checkSphericity(); 
  return 3./2.*(_spher[1]+_spher[2]); 
}

double EvtShapes::aplanarity() {
  checkSphericity(); 
  return 3./2.*_spher[2];
}

double EvtShapes::planarity() {
  checkSphericity(); 
  return _spher[1]-_spher[2]; 
}

Axis EvtShapes::sphericityAxis() {
  checkSphericity(); 
  return _spherAxis[0]; 
}

vector<double> EvtShapes::sphericityEigenValues() {
  checkSphericity(); 
  return _spher; 
}

vector<Axis> EvtShapes::sphericityEigenVectors() {
  checkSphericity(); 
  return _spherAxis; 
}

// jet mass related
double EvtShapes::Mhigh2() {
  checkHemispheres();
  return _mPlus; 
} 

double EvtShapes::Mlow2() {
  checkHemispheres();
  return _mMinus; 
} 

double EvtShapes::Mdiff2() {
  checkHemispheres();
  return _mPlus-_mMinus; 
} 

// jet broadening
double EvtShapes::Bmax() {
  checkBroadening(); 
  return _bPlus;
}

double EvtShapes::Bmin() {
  checkBroadening(); 
  return _bMinus;
}

double EvtShapes::Bsum() {
  checkBroadening(); 
  return _bPlus+_bMinus;
}

double EvtShapes::Bdiff() {
  checkBroadening(); 
  return _bPlus-_bMinus;
}

void EvtShapes::normalizeEEC(vector<double> & hi, long evts) {
  for (unsigned int bin = 0; bin < hi.size(); bin++) bin /= (hi.size()*evts);
}

void EvtShapes::bookEEC(vector<double> & hi) {
  // hi is the histogram.  It is understood that hi.front() contains
  // the bin [-1 < cos(chi) < -1+delta] and hi.back() the bin [1-delta
  // < cos(chi) < 1].  Here, delta = 2/hi.size().
  Energy Evis(0.*MeV);
  Lorentz5Momentum p_i, p_j; 
  for (unsigned int bin = 0; bin < hi.size(); bin++) {
    double delta = 2./((double) hi.size());
    double coschi = -1+bin*delta;
    if (_pv.size() > 1) {
      for (unsigned int i = 0; i < _pv.size()-1; i++) {
	p_i = _pv[i]->momentum();
	Evis += p_i.e(); 
	for (unsigned int j = i+1; j < _pv.size(); j++) {
	  p_j = _pv[j]->momentum();
	  double diff = abs(coschi-cos( p_i.vect().angle(p_j.vect()) )); 
	  if (delta > diff) 
	    hi[bin] += p_i.e()*p_j.e()/sqr(MeV);
	}
      }
    }
    hi[bin] /= (Evis*Evis)/sqr(MeV);
  }
}

double EvtShapes::AEEC(vector<double> & hi, double& coschi) {
  if (coschi > 0. && coschi <= 1.) {
    int i = (int) floor((-coschi+1.)/2.*hi.size()); 
    int j = (int) floor((coschi+1.)/2.*hi.size()); 
    return hi[i]-hi[j];
  } else {
    return 1e99;
  }
}

// private methods
inline void EvtShapes::checkLinTen() {
  if (!_linTenDone) {
    _linTenDone = true;
    diagonalizeTensors(true, _useCmBoost); 
  }
}

inline void EvtShapes::checkSphericity() {
  if (!_spherDone) {
    _spherDone = true;
    diagonalizeTensors(false, _useCmBoost); 
  }
}

void EvtShapes::diagonalizeTensors(bool linear, bool cmboost) {  
  // initialize
  HepSymMatrix Theta = HepSymMatrix(3);
  for(int i=0; i<3; i++) 
    for(int j=0; j<3; j++) 
      Theta[i][j] = 0.0;
  double sum = 0.; 
  Momentum3 sumvec = Momentum3();
  vector<double> lam;
  vector<Axis> n; 
  // get cm-frame
  Lorentz5Momentum pcm = Lorentz5Momentum(); 
  tPVector::const_iterator cit;
  Boost beta; 
  if (cmboost) {
    for(cit=_pv.begin(); cit != _pv.end(); ++cit) 
      pcm += (*cit)->momentum();    
    beta = pcm.findBoostToCM(); 
  }
  // get Theta_ij
  for(cit=_pv.begin(); cit != _pv.end(); ++cit) {
    Lorentz5Momentum dum = (*cit)->momentum();
    if (cmboost) dum.boost( beta );
    Momentum3 pvec = dum.vect();
    double pvec_MeV[3] = {pvec.x()/MeV, pvec.y()/MeV, pvec.z()/MeV};
    if (pvec.mag() > 0*MeV) {
      sumvec += pvec;
      if (linear) 
	sum += pvec.mag()/MeV;
      else 
	sum += pvec.mag2()/sqr(MeV);
      for(int i=0; i<3; i++) 
	for(int j=i; j<3; j++) 
	  if (linear) 
	    Theta[i][j] += (pvec_MeV[i])*(pvec_MeV[j])*MeV/(pvec.mag());      
	  else 
	    Theta[i][j] += (pvec_MeV[i])*(pvec_MeV[j]);      
    }
  }
  Theta /= sum;      
  // diagonalize it
  HepMatrix U = diagonalize(&Theta);
  for(int i=0; i<3; i++) {
    lam.push_back( Theta[i][i] );
    Axis ndum(U[0][i], U[1][i], U[2][i]);
    n.push_back( ndum ); 
  }
  // sort according to size of eigenvalues
  // such that lam[0] > lam[1] > lam[2]
  if (lam[0] < lam[1]) {
    swap(lam[0], lam[1]); 
    swap(n[0], n[1]);     
  }
  if (lam[0] < lam[2]) {
    swap(lam[0], lam[2]); 
    swap(n[0], n[2]);         
  }
  if (lam[1] < lam[2]) {
    swap(lam[1], lam[2]); 
    swap(n[1], n[2]);     
  }
  if (linear) {
    _linTen = lam; 
    _linTenAxis = n; 
  } else {
    _spher = lam; 
    _spherAxis = n; 
  }
}

inline void EvtShapes::checkThrust() {
  if (!_thrustDone) {
    _thrustDone = true;
    calculateThrust(); 
  }
}

void EvtShapes::calcT(const vector<Momentum3> &p, 
		      Energy2 &t, Axis & taxis) {
  Energy2 tval;
  t = Energy2();
  Vector3<Energy2> tv;
  Momentum3 ptot;
  for (unsigned int k=1; k < p.size(); k++) {
    for (unsigned int j=0; j<k; j++) {
      tv = p[j].cross(p[k]);
      ptot = Momentum3();
      for (unsigned int l=0; l<p.size(); l++) {
	if (l!=j && l!=k) {
	  if (p[l]*tv > Energy3()) { 
	    ptot += p[l];
	  } else {
	    ptot -= p[l];
	  }
	}
      }
      vector<Momentum3> cpm;
      cpm.push_back(ptot-p[j]-p[k]);
      cpm.push_back(ptot-p[j]+p[k]);
      cpm.push_back(ptot+p[j]-p[k]);
      cpm.push_back(ptot+p[j]+p[k]);
      for (vector<Momentum3>::iterator it = cpm.begin();
	   it != cpm.end(); it++) {
	tval = it->mag2();
	if (tval > t) {
	  t = tval;
	  taxis = it->unit();
	}
      }
    }
  }
}

void EvtShapes::calcM(const vector<Momentum3> &p, 
		      Energy2 &m, Axis & maxis) {
  Energy2 mval;
  m = Energy2();
  Momentum3 tv, ptot;
  for (unsigned int j=0; j < p.size(); j++) {
    tv = p[j];
    ptot = Momentum3();
    for (unsigned int l=0; l<p.size(); l++) {
      if (l!=j) {
	if (p[l]*tv > Energy2()) { 
	  ptot += p[l];
	} else {
	  ptot -= p[l];
	}
      }
    }
    vector<Momentum3> cpm;
    cpm.push_back(ptot-p[j]);
    cpm.push_back(ptot+p[j]);
    for (vector<Momentum3>::iterator it = cpm.begin();
	 it != cpm.end(); it++) {
      mval = it->mag2();
      if (mval > m) {
	m = mval;
	maxis = it->unit();
      }
    }
  }
}


void EvtShapes::calculateThrust() { 
  // explicitly calculate in units of MeV
  // algorithm based on Brandt/Dahmen Z Phys C1 (1978)
  // and 'tasso' code from HERWIG
  // assumes all momenta in cm system, no explicit boost performed here!
  // unlike for C and D

  _thrust.clear();
  _thrustAxis.clear(); 

  if (_pv.size() < 2) {
    for (int i=0; i<3; i++) {
      _thrust.push_back(-1);
      _thrustAxis.push_back(Axis());
    }
    return;
  }

  // thrust
  vector<Momentum3 > p;
  p.clear();
  Energy psum = 0*MeV;
  for (unsigned int l=0; l<_pv.size(); l++) {
    p.push_back(_pv[l]->momentum().vect());
    psum += p.back().mag();
  }

  Energy2 val; 
  Axis axis;
  if (p.size() == 2) {
    _thrust.push_back(1.0);
    _thrust.push_back(0.0);
    _thrust.push_back(0.0);
    axis = p[0].unit();
    if (axis.z() < 0) axis = -axis;
    _thrustAxis.push_back(axis);
    _thrustAxis.push_back(axis.orthogonal());
    axis = _thrustAxis[0].cross(_thrustAxis[1]);
    return;
  }

  if (p.size() == 3) {
    if (p[0].mag2() < p[1].mag2()) swap(p[0], p[1]);
    if (p[0].mag2() < p[2].mag2()) swap(p[0], p[2]);
    if (p[1].mag2() < p[2].mag2()) swap(p[1], p[2]);
    // thrust
    axis = p[0].unit();
    if (axis.z() < 0) axis = -axis;
    _thrust.push_back(2.*p[0].mag()/psum);
    _thrustAxis.push_back(axis);
    // major
    axis = (p[1] - (axis*p[1])*axis).unit();
    if (axis.x() < 0) axis = -axis;
    _thrust.push_back((abs(p[1]*axis) + abs(p[2]*axis))/psum);
    _thrustAxis.push_back(axis);
    // minor
    _thrust.push_back(0.0);
    axis = _thrustAxis[0].cross(_thrustAxis[1]);
    _thrustAxis.push_back(axis);
    return;
  }

  // ACHTUNG special case with >= 4 coplanar particles will still fail. 
  // probably not too important... 

  calcT(p, val, axis);
  _thrust.push_back(sqrt(val)/psum);
  if (axis.z() < 0) axis = -axis;
  _thrustAxis.push_back(axis.unit()); 

  //major 
  p.clear();
  Momentum3 par;
  for (unsigned int l=0; l<_pv.size(); l++) {
    par = (_pv[l]->momentum().vect() * axis.unit())*axis.unit();
    p.push_back(_pv[l]->momentum().vect() - par);
  }
  calcM(p, val, axis);
  _thrust.push_back(sqrt(val)/psum);
  if (axis.x() < 0) axis = -axis;
  _thrustAxis.push_back(axis.unit()); 
  
  // minor
  if (_thrustAxis[0]*_thrustAxis[1] < 1e-10) {
    Energy value;
    axis = _thrustAxis[0].cross(_thrustAxis[1]);
    _thrustAxis.push_back(axis); 
    for (unsigned int l=0; l<_pv.size(); l++) {
      value += abs(axis*(_pv[l]->momentum().vect()));
    }
    _thrust.push_back(value/psum);
  } else {
    _thrust.push_back(-1.0);
    _thrustAxis.push_back(Axis()); 
  }
}

inline void EvtShapes::checkHemispheres() {
  if (!_hemDone) {
    _hemDone = true;
    calcHemisphereMasses(); 
  }
}

void EvtShapes::calcHemisphereMasses() {
  Lorentz5Momentum pos, neg;
  tPVector::const_iterator cit;
  for(cit = _pv.begin(); cit != _pv.end(); cit++) 
    if ((*cit)->momentum().vect()*thrustAxis() > Energy()) 
      pos += (*cit)->momentum();
    else neg += (*cit)->momentum();
  _mPlus = pos.m()/(pos+neg).e();
  _mPlus *= _mPlus;
  _mMinus = neg.m()/(pos+neg).e();
  _mMinus *= _mMinus; 
  if (_mPlus < _mMinus) swap(_mPlus, _mMinus);
}

inline void EvtShapes::checkBroadening() {
  if (!_broadDone) {
    _broadDone = true; 
    calcBroadening(); 
  }
}

void EvtShapes::calcBroadening() {
  Energy pos, neg, den; 
  pos = neg = den = 0.0*MeV;
  tPVector::const_iterator cit;
  for(cit = _pv.begin(); cit != _pv.end(); cit++) {
    if ((*cit)->momentum().vect()*thrustAxis() > Energy()) 
      pos += (*cit)->momentum().perp(thrustAxis());
    else 
      neg += (*cit)->momentum().perp(thrustAxis());
    den += (*cit)->momentum().vect().mag();
  }
  _bPlus = 0.5*pos/den;
  _bMinus = 0.5*neg/den;
  if (_bPlus < _bMinus) swap(_bPlus, _bMinus);
}

