#include "EvtShapes.h"

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
  return(Ebeam > 0 ? p.vect().mag()/Ebeam : -1.); 
}

double EvtShapes::getXi(const Lorentz5Momentum & p, 
			       const Energy & Ebeam) {
  return((Ebeam > 0 && p.vect().mag() > 0) ? 
	 log(Ebeam/p.vect().mag()) : -1.); 
}

Energy EvtShapes::getPt(const Lorentz5Momentum & p) {
  return p.perp(); 
}

Energy EvtShapes::getRapidity(const Lorentz5Momentum & p) {
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

Vector3 EvtShapes::thrustAxis() {
  checkThrust(); 
  return _thrustAxis[0];
}

Vector3 EvtShapes::majorAxis() {
  checkThrust(); 
  return _thrustAxis[1];
}

Vector3 EvtShapes::minorAxis() {
  checkThrust(); 
  return _thrustAxis[2];
}

// linear tensor related
vector<double> EvtShapes::linTenEigenValues() {
  checkLinTen(); 
  return _linTen; 
}

vector<Vector3> EvtShapes::linTenEigenVectors() {
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

Vector3 EvtShapes::sphericityAxis() {
  checkSphericity(); 
  return _spherAxis[0]; 
}

vector<double> EvtShapes::sphericityEigenValues() {
  checkSphericity(); 
  return _spher; 
}

vector<Vector3> EvtShapes::sphericityEigenVectors() {
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
  for (int bin = 0; bin < hi.size(); bin++) 
    bin /= (hi.size()*evts);
}

void EvtShapes::bookEEC(vector<double> & hi) {
  // hi is the histogram.  It is understood that hi.front() contains
  // the bin [-1 < cos(chi) < -1+delta] and hi.back() the bin [1-delta
  // < cos(chi) < 1].  Here, delta = 2/hi.size().
  Energy Evis;
  Lorentz5Momentum p_i, p_j; 
  for (int bin = 0; bin < hi.size(); bin++) {
    double delta = 2./((double) hi.size());
    double coschi = -1+bin*delta;
    if (_pv.size() > 1) {
      for (int i = 0; i < _pv.size()-1; i++) {
	p_i = _pv[i]->momentum();
	Evis += p_i.e(); 
	for (int j = i+1; j < _pv.size(); j++) {
	  p_j = _pv[j]->momentum();
	  double diff = abs(coschi-cos( p_i.vect().angle(p_j.vect()) )); 
	  if (delta > diff) 
	    hi[bin] += p_i.e()*p_j.e();
	}
      }
    }
    hi[bin] /= (Evis*Evis);
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
  Vector3 sumvec = Vector3();
  vector<double> lam;
  vector<Vector3> n; 
  // get cm-frame
  Lorentz5Momentum pcm = Lorentz5Momentum(); 
  tPVector::const_iterator cit;
  Vector3 beta; 
  if (cmboost) {
    for(cit=_pv.begin(); cit != _pv.end(); ++cit) 
      pcm += (*cit)->momentum();    
    beta = pcm.findBoostToCM(); 
  }
  // get Theta_ij
  for(cit=_pv.begin(); cit != _pv.end(); ++cit) {
    Lorentz5Momentum dum = (*cit)->momentum();
    if (cmboost) dum.boost( beta );
    Vector3 pvec = dum.vect();
    if (pvec.mag() > 0) {
      sumvec += pvec;
      if (linear) 
	sum += pvec.mag();
      else 
	sum += pvec.mag2();
      for(int i=0; i<3; i++) 
	for(int j=i; j<3; j++) 
	  if (linear) 
	    Theta[i][j] += (pvec[i])*(pvec[j])/(pvec.mag());      
	  else 
	    Theta[i][j] += (pvec[i])*(pvec[j]);      
    }
  }
  Theta /= sum;      
  // diagonalize it
  HepMatrix U = diagonalize(&Theta);
  for(int i=0; i<3; i++) {
    lam.push_back( Theta[i][i] );
    Vector3 ndum;
    for(int j=0; j<3; j++) 
      ndum[j] = U[j][i]; 
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

void EvtShapes::calcT(const vector<Vector3> &p, double &t, Vector3 &taxis) {
  double tval;
  t = 0.0;
  Vector3 tv, ptot;
  vector<Vector3> cpm;
  for (int k=1; k < p.size(); k++) {
    for (int j=0; j<k; j++) {
      tv = p[j].cross(p[k]);
      ptot = Vector3();
      for (int l=0; l<p.size(); l++) {
	if (l!=j && l!=k) {
	  if (p[l]*tv > 0.0) { 
	    ptot += p[l];
	  } else {
	    ptot -= p[l];
	  }
	}
      }
      cpm.clear();
      cpm.push_back(ptot-p[j]-p[k]);
      cpm.push_back(ptot-p[j]+p[k]);
      cpm.push_back(ptot+p[j]-p[k]);
      cpm.push_back(ptot+p[j]+p[k]);
      for (vector<Vector3>::iterator it = cpm.begin();
	   it != cpm.end(); it++) {
	tval = (*it).mag2();
	if (tval > t) {
	  t = tval;
	  taxis = *it;
	}
      }
    }
  }
}

void EvtShapes::calcM(const vector<Vector3> &p, double &m, Vector3 &maxis) {
  double mval;
  m = 0.0;
  Vector3 tv, ptot;
  vector<Vector3> cpm;
  for (int j=0; j < p.size(); j++) {
    tv = p[j];
    ptot = Vector3();
    for (int l=0; l<p.size(); l++) {
      if (l!=j) {
	if (p[l]*tv > 0.0) { 
	  ptot += p[l];
	} else {
	  ptot -= p[l];
	}
      }
    }
    cpm.clear();
    cpm.push_back(ptot-p[j]);
    cpm.push_back(ptot+p[j]);
    for (vector<Vector3>::iterator it = cpm.begin();
	 it != cpm.end(); it++) {
      mval = (*it).mag2();
      if (mval > m) {
	m = mval;
	maxis = *it;
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
      _thrustAxis.push_back(Vector3());
    }
    return;
  }

  // thrust
  vector<Vector3> p;
  p.clear();
  double psum = 0.0;
  for (int l=0; l<_pv.size(); l++) {
    p.push_back(_pv[l]->momentum().vect()/MeV);
    psum += p.back().mag();
  }

  double val; 
  Vector3 axis;
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
  Vector3 par;
  for (int l=0; l<_pv.size(); l++) {
    par = ((_pv[l]->momentum().vect()/MeV)*axis.unit())*axis.unit();
    p.push_back(_pv[l]->momentum().vect()/MeV - par);
  }
  calcM(p, val, axis);
  _thrust.push_back(sqrt(val)/psum);
  if (axis.x() < 0) axis = -axis;
  _thrustAxis.push_back(axis.unit()); 
  
  // minor
  if (_thrustAxis[0]*_thrustAxis[1] < 1e-10) {
    val = 0.;
    axis = _thrustAxis[0].cross(_thrustAxis[1]);
    _thrustAxis.push_back(axis); 
    for (int l=0; l<_pv.size(); l++) {
      val += abs(axis*(_pv[l]->momentum().vect()/MeV));
    }
    _thrust.push_back(val/psum);
  } else {
    _thrust.push_back(-1.0);
    _thrustAxis.push_back(Vector3()); 
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
    if ((*cit)->momentum().vect()*thrustAxis() > 0) 
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
    if ((*cit)->momentum().vect()*thrustAxis() > 0) 
      pos += (*cit)->momentum().perp(thrustAxis());
    else 
      neg += (*cit)->momentum().perp(thrustAxis());
    den += (*cit)->momentum().vect().mag();
  }
  _bPlus = pos/den/2.;
  _bMinus = neg/den/2;
  if (_bPlus < _bMinus) swap(_bPlus, _bMinus);
}

