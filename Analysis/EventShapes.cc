// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EventShapes class.
//

#include "EventShapes.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EventShapes.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EventShapes::~EventShapes() {}

void EventShapes::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void EventShapes::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<EventShapes> EventShapes::initEventShapes;
// Definition of the static class description member.

void EventShapes::Init() {

  static ClassDocumentation<EventShapes> documentation
    ("There is no documentation for the EventShapes class");

}

void EventShapes::diagonalizeTensors(bool linear, bool cmboost) {  
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
    for(unsigned int ix=0;ix<_pv.size();++ix) pcm += _pv[ix];    
    beta = pcm.findBoostToCM();
  }
  // get Theta_ij
  for(unsigned int ix=0;ix<_pv.size();++ix)
    {
      Lorentz5Momentum dum(_pv[ix]);
      if (cmboost) dum.boost( beta );
      Vector3 pvec = dum.vect();
      if (pvec.mag() > 0) 
	{
	  sumvec += pvec;
	  if (linear) sum += pvec.mag();
	  else        sum += pvec.mag2();
	  for(int i=0; i<3; i++) 
	    {
	      for(int j=i; j<3; j++) 
		{
		  if (linear) Theta[i][j] += (pvec[i])*(pvec[j])/(pvec.mag());      
		  else        Theta[i][j] += (pvec[i])*(pvec[j]);
		}
	    }
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

void EventShapes::calculateThrust() { 
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
  double psum = 0.0;
  for(unsigned int l=0; l<_pv.size(); l++) 
    {
      p.push_back(_pv[l].vect()/MeV);
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
  Vector3 par;
  for (unsigned int l=0; l<_pv.size(); l++) 
    {
      par   = (p[l]*axis.unit())*axis.unit();
      p[l]  = p[l] - par;
    }
  calcM(p, val, axis);
  _thrust.push_back(sqrt(val)/psum);
  if (axis.x() < 0) axis = -axis;
  _thrustAxis.push_back(axis.unit()); 

  // minor
  if (_thrustAxis[0]*_thrustAxis[1] < 1e-10) 
    {
      val = 0.;
      axis = _thrustAxis[0].cross(_thrustAxis[1]);
      _thrustAxis.push_back(axis); 
      for (unsigned int l=0; l<_pv.size(); l++) val += abs(axis*_pv[l]);
      _thrust.push_back(val/psum);
    } 
  else 
    {
      _thrust.push_back(-1.0);
      _thrustAxis.push_back(Vector3()); 
    }
}

void EventShapes::calcT(const vector<Vector3> &p, double &t, Vector3 &taxis) {
  double tval;
  t = 0.0;
  Vector3 tv, ptot;
  vector<Vector3> cpm;
  for (unsigned int k=1; k < p.size(); k++) {
    for (unsigned int j=0; j<k; j++) {
      tv = p[j].cross(p[k]);
      ptot = Vector3();
      for (unsigned int l=0; l<p.size(); l++) {
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

void EventShapes::calcM(const vector<Vector3> &p, double &m, Vector3 &maxis) {
  double mval;
  m = 0.0;
  Vector3 tv, ptot;
  vector<Vector3> cpm;
  for (unsigned int j=0; j < p.size(); j++) {
    tv = p[j];
    ptot = Vector3();
    for (unsigned int l=0; l<p.size(); l++) {
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

void EventShapes::bookEEC(vector<double> & hi) {
  // hi is the histogram.  It is understood that hi.front() contains
  // the bin [-1 < cos(chi) < -1+delta] and hi.back() the bin [1-delta
  // < cos(chi) < 1].  Here, delta = 2/hi.size().
  Energy Evis(0.);
  for (unsigned int bin = 0; bin < hi.size(); bin++) {
    double delta = 2./((double) hi.size());
    double coschi = -1+bin*delta;
    if (_pv.size() > 1) {
      for (unsigned int i = 0; i < _pv.size()-1; i++) {
	Evis += _pv[i].e(); 
	for (unsigned int j = i+1; j < _pv.size(); j++) {
	  double diff = abs(coschi-cos( _pv[i].vect().angle(_pv[j].vect()) )); 
	  if (delta > diff) 
	    hi[bin] += _pv[i].e()*_pv[j].e();
	}
      }
    }
    hi[bin] /= (Evis*Evis);
  }
}

