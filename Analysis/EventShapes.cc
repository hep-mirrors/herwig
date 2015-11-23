// -*- C++ -*-
//
// EventShapes.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EventShapes class.
//

#include "EventShapes.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

NoPIOClassDescription<EventShapes> EventShapes::initEventShapes;
// Definition of the static class description member.

void EventShapes::Init() {

  static ClassDocumentation<EventShapes> documentation
    ("There is no documentation for the EventShapes class");

}

void EventShapes::calcHemisphereMasses() {
  Lorentz5Momentum pos, neg;
  Energy pden(ZERO),epos(ZERO),eneg(ZERO);
  for(unsigned int ix=0;ix<_pv.size();++ix)
    {
      if(_pv[ix].vect() * thrustAxis() > ZERO)
	{
	  pos  += _pv[ix];
	  // can be replaced with, once perp() is giving non-nan results
	  //	  for nearly parallel vectors.  
	  // epos += _pv[ix].perp(thrustAxis());
	  epos += _pv[ix].vect().cross(thrustAxis()).mag();
	}
      else
	{
	  neg  += _pv[ix];
	  // see above
	  //	  eneg += _pv[ix].perp(thrustAxis()); 
	  eneg += _pv[ix].vect().cross(thrustAxis()).mag();
	}
      pden += _pv[ix].vect().mag();	 
    }
  // denominator and masses
  Energy2 den(sqr(pos.e()+neg.e()));
  _mPlus = pos.m2()/den;
  _mMinus = neg.m2()/den;
  if (_mPlus < _mMinus) swap(_mPlus, _mMinus);
  // jet broadening
  _bPlus  = 0.5*epos/pden;
  _bMinus = 0.5*eneg/pden;
  if (_bPlus < _bMinus) swap(_bPlus, _bMinus);
}

vector<double> EventShapes::eigenvalues(const double T[3][3]) {

  // b, c, d are the coefficients of the characteristic polynomial, 
  // a lambda^3 + b lambda^2 + c lambda + d
  // where a is chosen to be +1.
  double t11, t12, t13, t22, t23, t33;
  t11 = T[0][0]; t12 = T[0][1]; t13 = T[0][2]; 
  t22 = T[1][1]; t23 = T[1][2]; t33 = T[2][2]; 
  double b = -(t11 + t22 + t33);
  double c = t11*t22 + t11*t33 + t22*t33 - sqr(t12) - sqr(t13) - sqr(t23);
  double d = - t11*t22*t33 - 2.*t12*t23*t13 
    + t11*sqr(t23) + t22*sqr(t13) + t33*sqr(t12); 
  
  // use Cardano's formula to compute the zeros 
  double p = (3.*c - sqr(b))/3.;
  double q = (2.*sqr(b)*b - 9.*b*c + 27.*d)/27.;
  // check diskriminant to double precision
  vector<double> lambda;
  if (4.*p*sqr(p) + 27.*sqr(q) > 2.0e-16) {
    for (unsigned int i=0; i<3; ++i) {
      lambda.push_back(-1.);
    }
    cerr << "EventShapes::eigenvalues: found D = "
	 << 4.*p*sqr(p) + 27.*sqr(q) 
	 << " > 0! No real Eigenvalues!\n";
  } else {
    // get solutions
    double alpha = acos(-q/2.*sqrt(-27./(p*p*p)))/3.;
    double w = sqrt(-4.*p/3.);
    lambda.push_back(w*cos(alpha) - b/3.);
    lambda.push_back(-w*cos(alpha+M_PI/3.) - b/3.);
    lambda.push_back(-w*cos(alpha-M_PI/3.) - b/3.);
  }

  // sort according to size of eigenvalues
  // such that lambda[0] > lambda[1] > lambda[2]
  if (lambda[0] < lambda[1]) {
    swap(lambda[0], lambda[1]); 
  }
  if (lambda[0] < lambda[2]) {
    swap(lambda[0], lambda[2]); 
  }
  if (lambda[1] < lambda[2]) {
    swap(lambda[1], lambda[2]); 
  }

  return lambda;
}


Axis EventShapes::eigenvector(const double T[3][3], const double &lam) {
  // set up matrix of system to be solved
  double a11, a12, a13, a23, a33;
  a11 = T[0][0] - lam; 
  a12 = T[0][1]; 
  a13 = T[0][2]; 
  a23 = T[1][2]; 
  a33 = T[2][2] - lam;

  // intermediate steps from gauss type algorithm
  double b1, b2, b4;
  b1 = a11*a33 - sqr(a13); 
  b2 = a12*a33 - a13*a23; 
  b4 = a11*a23 - a12*a13;

  // eigenvector
  Axis u(b2, -b1, b4);

  return u.unit();
}


vector<Axis> EventShapes::
eigenvectors(const double T[3][3], const vector<double> &lam) {
  vector<Axis> n;
  for (unsigned int i=0; i<3; ++i) {
    n.push_back(eigenvector(T, lam[i]));
  }
  return n;
}

void EventShapes::diagonalizeTensors(bool linear, bool cmboost) {  
  // initialize
  double Theta[3][3];
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      Theta[i][j] = 0.0;
    }
  }
  double sum = 0.; 
  Momentum3 sumvec;
  vector<double> lam;
  vector<Axis> n; 
  // get cm-frame
  Lorentz5Momentum pcm = Lorentz5Momentum(); 
  Boost beta; 
  if (cmboost) {
    for(unsigned int ix=0;ix<_pv.size();++ix) {
      pcm += _pv[ix];    
    }
    beta = pcm.findBoostToCM();
  }
  // get Theta_ij
  for(unsigned int ix=0;ix<_pv.size();++ix) {
    Lorentz5Momentum dum(_pv[ix]);
    if (cmboost) {
      dum.boost( beta );
    }
    Momentum3 pvec = dum.vect();
    double pvec_MeV[3] = {pvec.x()/MeV, pvec.y()/MeV, pvec.z()/MeV};
    if (pvec.mag() > ZERO) {
      sumvec += pvec;
      if (linear) {
	sum += pvec.mag()*UnitRemoval::InvE;
      } else {
	sum += pvec.mag2()*UnitRemoval::InvE2;
      }
      for(int i=0; i<3; ++i) {
	for(int j=i; j<3; ++j) {
	  if (linear) {
	    Theta[i][j] += (pvec_MeV[i])*(pvec_MeV[j])*MeV/(pvec.mag());      
	  } else {
	    Theta[i][j] += (pvec_MeV[i])*(pvec_MeV[j]);
	  }
	}
      }
    }
  }
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      Theta[i][j] /= sum;
    }
  }
  
  // diagonalize it
  lam = eigenvalues(Theta);
  n = eigenvectors(Theta, lam);

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
    for (int i=0; i<3; ++i) {
      _thrust.push_back(-1);
      _thrustAxis.push_back(Axis());
    }
    return;
  }

  // thrust
  vector<Momentum3> p;
  Energy psum = ZERO;
  for(unsigned int l=0; l<_pv.size(); ++l) 
    {
      p.push_back(_pv[l].vect());
      psum += p.back().mag();
    }

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
  Energy2 val;
  calcT(p, val, axis);
  _thrust.push_back(sqrt(val)/psum);
  if (axis.z() < 0) axis = -axis;
  _thrustAxis.push_back(axis.unit()); 

  //major 
  Momentum3 par;
  for (unsigned int l=0; l<_pv.size(); ++l) 
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
      Energy eval = ZERO;
      axis = _thrustAxis[0].cross(_thrustAxis[1]);
      _thrustAxis.push_back(axis); 
      for (unsigned int l=0; l<_pv.size(); ++l) 
	eval += abs(axis*_pv[l].vect());
      _thrust.push_back(eval/psum);
    } 
  else 
    {
      _thrust.push_back(-1.0);
      _thrustAxis.push_back(Axis()); 
    }
}

void EventShapes::calcT(const vector<Momentum3> &p, Energy2 &t, Axis &taxis) {
  Energy2 tval;
  t = ZERO;
  ThreeVector<Energy2> tv;
  Momentum3 ptot;
  vector<Momentum3> cpm;
  for (unsigned int k=1; k < p.size(); ++k) {
    for (unsigned int j=0; j<k; ++j) {
      tv = p[j].cross(p[k]);
      ptot = Momentum3();
      for (unsigned int l=0; l<p.size(); ++l) {
	if (l!=j && l!=k) {
	  if (p[l]*tv > ZERO) { 
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
      for (vector<Momentum3>::iterator it = cpm.begin();
	   it != cpm.end(); ++it) {
	tval = it->mag2();
	if (tval > t) {
	  t = tval;
	  taxis = it->unit();
	}
      }
    }
  }
}

void EventShapes::calcM(const vector<Momentum3> &p, Energy2 &m, Axis &maxis) {
  Energy2 mval;
  m = ZERO;
  Momentum3 tv, ptot;
  vector<Momentum3> cpm;
  for (unsigned int j=0; j < p.size(); ++j) {
    tv = p[j];
    ptot = Momentum3();
    for (unsigned int l=0; l<p.size(); ++l) {
      if (l!=j) {
	if (p[l]*tv > ZERO) { 
	  ptot += p[l];
	} else {
	  ptot -= p[l];
	}
      }
    }
    cpm.clear();
    cpm.push_back(ptot-p[j]);
    cpm.push_back(ptot+p[j]);
    for (vector<Momentum3>::iterator it = cpm.begin();
	 it != cpm.end(); ++it) {
      mval = it->mag2();
      if (mval > m) {
	m = mval;
	maxis = it->unit();
      }
    }
  }
}

void EventShapes::bookEEC(vector<double> & hi) {
  // hi is the histogram.  It is understood that hi.front() contains
  // the bin [-1 < cos(chi) < -1+delta] and hi.back() the bin [1-delta
  // < cos(chi) < 1].  Here, delta = 2/hi.size().
  Energy Evis(ZERO);
  for (unsigned int bin = 0; bin < hi.size(); ++bin) {
    double delta = 2./hi.size();
    double coschi = -1+bin*delta;
    if (_pv.size() > 1) {
      for (unsigned int i = 0; i < _pv.size()-1; ++i) {
	Evis += _pv[i].e(); 
	for (unsigned int j = i+1; j < _pv.size(); ++j) {
	  double diff = abs(coschi-cos( _pv[i].vect().angle(_pv[j].vect()) )); 
	  if (delta > diff) 
	    hi[bin] += _pv[i].e()*_pv[j].e() / MeV2;
	}
      }
    }
    hi[bin] /= (Evis*Evis) / MeV2;
  }
}

