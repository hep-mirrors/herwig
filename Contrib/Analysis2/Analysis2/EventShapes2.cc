// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de
// Copyright (C) 2002-2007 The Herwig Collaboration

#include "EventShapes2.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EventShapes2.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

EventShapes2::~EventShapes2() {}

void EventShapes2::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _omthr << _maj << _min << _obl << _sph << _apl << _pla
     << _c << _d << _mhi << _mlo << _mdiff << _bmax << _bmin << _bsum << _bdiff;
}

void EventShapes2::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _omthr >> _maj >> _min >> _obl >> _sph >> _apl >> _pla
     >> _c >> _d >> _mhi >> _mlo >> _mdiff >> _bmax >> _bmin >> _bsum >> _bdiff;
}

ClassDescription<EventShapes2> EventShapes2::initEventShapes2;
// Definition of the static class description member.

void EventShapes2::Init() {

  static ClassDocumentation<EventShapes2> documentation
    ("Analyse event shapes");


  static Parameter<EventShapes2,string> interfaceOneMinusThrust
    ("OneMinusThrust",
     "Options for 1-T",
     &EventShapes2::_omthr, "",
     false, false);
  
  static Parameter<EventShapes2,string> interfaceThrustMajor
    ("ThrustMajor",
     "Options for thrust major",
     &EventShapes2::_maj, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceThrustMinor
    ("ThrustMinor",
     "Options for thrust minor",
     &EventShapes2::_min, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceOblateness
    ("Oblateness",
     "Options for Oblateness",
     &EventShapes2::_obl, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceSphericity
    ("Sphericity",
     "Options for Sphericity",
     &EventShapes2::_sph, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceAplanarity
    ("Aplanarity",
     "Options for Aplanarity",
     &EventShapes2::_apl, "",
     false, false);

  static Parameter<EventShapes2,string> interfacePlanarity
    ("Planarity",
     "Options for Planarity",
     &EventShapes2::_pla, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceCParameter
    ("CParameter",
     "Options for C parameter",
     &EventShapes2::_c, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceDParameter
    ("DParameter",
     "Options for D parameter",
     &EventShapes2::_d, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceMHigh
    ("MHigh",
     "Options for high hemisphere mass",
     &EventShapes2::_mhi, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceMLow
    ("MLow",
     "Options for low hemisphere mass",
     &EventShapes2::_mlo, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceMDiff
    ("MDiff",
     "Options for difference in hemisphere masses",
     &EventShapes2::_mdiff, "",
     false, false);


  static Parameter<EventShapes2,string> interfaceBMax
    ("BMax",
     "Options for wide jet broadening",
     &EventShapes2::_bmax, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceBMin
    ("BMin",
     "Options for narrow jet broadening",
     &EventShapes2::_bmin, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceBSum
    ("BSum",
     "Options for sum of jet broadening measures",
     &EventShapes2::_bsum, "",
     false, false);

  static Parameter<EventShapes2,string> interfaceBDiff
    ("BDiff",
     "Options for difference of jet broadening measures",
     &EventShapes2::_bdiff, "",
     false, false);

}

inline void EventShapes2::dofinish() {

  finish("OneMinusThrust");
  finish("ThrustMajor");
  finish("ThrustMinor");
  finish("Oblateness");
  finish("Sphericity");
  finish("Aplanarity");
  finish("Planarity");
  finish("CParameter");
  finish("DParameter");
  finish("MHigh");
  finish("MLow");
  finish("MDiff");
  finish("BMax");
  finish("BMin");
  finish("BSum");
  finish("BDiff");

  Analysis2Base::dofinish();

}

inline void EventShapes2::doinit() throw(InitException) {
  Analysis2Base::doinit();

  int plotFlags = HistogramOutput::Ylog | HistogramOutput::Frame | HistogramOutput::Errorbars;

  Histogram2Options options (plotFlags);

  insert("OneMinusThrust", _omthr,options);
  insert("ThrustMajor", _maj,options);
  insert("ThrustMinor", _min,options);
  insert("Oblateness", _obl,options);
  insert("Sphericity", _sph,options);
  insert("Aplanarity", _apl,options);
  insert("Planarity", _pla,options);
  insert("CParameter", _c,options);
  insert("DParameter", _d,options);
  insert("MHigh", _mhi,options);
  insert("MLow", _mlo,options);
  insert("MDiff", _mdiff,options);
  insert("BMax", _bmax,options);
  insert("BMin", _bmin,options);
  insert("BSum", _bsum,options);
  insert("BDiff", _bdiff,options);

}

void EventShapes2::analyze(const tPVector &) {

  pair<vector<Lorentz5Momentum>,double> ev;

  while (*eventExtractor() >> ev) {

    reset(ev.first);  

    book(1.-thrust(), "OneMinusThrust", ev.second);
    book(thrustMajor(), "ThrustMajor", ev.second);
    book(thrustMinor(), "ThrustMinor", ev.second);
    book(oblateness(), "Oblateness", ev.second); 
    book(CParameter(), "CParameter", ev.second); 
    book(DParameter(), "DParameter", ev.second); 
    book(sphericity(), "Sphericity", ev.second);
    book(aplanarity(), "Aplanarity", ev.second);
    book(planarity(), "Planarity", ev.second); 
    book(Mhigh2(), "MHigh", ev.second);
    book(Mlow2(), "MLow", ev.second); 
    book(Mdiff2(), "MDiff", ev.second); 
    book(Bmax(), "BMax", ev.second); 
    book(Bmin(), "BMin", ev.second); 
    book(Bsum(), "BSum", ev.second); 
    book(Bdiff(), "BDiff", ev.second); 

  }
  
}

vector<double> EventShapes2::eigenvalues(const double T[3][3]) {

  vector<double> lambda;

  if (_pv.size() > 2) {

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
    // check diskriminant
    if (4.*p*sqr(p) + 27.*sqr(q) > 0) {
      for (unsigned int i=0; i<3; i++) {
	lambda.push_back(-1.);
      }
      cout << flush
	   << "EventShapes2::eigenvalues: found D > 0! \n"
	   << "Matrix doesn't have real Eigenvalues in this case\n"
	   << "(event with " << _pv.size() << " final state particles)\n";

      cout << flush
	   << "| " << T[0][0] << " " << T[0][1] << " " << T[0][2] << " |\n"
	   << "| " << T[1][0] << " " << T[1][1] << " " << T[1][2] << " |\n"
	   << "| " << T[2][0] << " " << T[2][1] << " " << T[2][2] << " |\n"
	   << flush;

    } else {
      // get solutions
      double alpha = acos(-q/2.*sqrt(-27./(p*p*p)))/3.;
      double w = sqrt(-4.*p/3.);
      lambda.push_back(w*cos(alpha) - b/3.);
      lambda.push_back(-w*cos(alpha+M_PI/3.) - b/3.);
      lambda.push_back(-w*cos(alpha-M_PI/3.) - b/3.);
    }

  } else {

    // linear and sphericity coincide for ideal back-to-back
    // as q_i q_j/(2 q^2) where q = p_1 = -p_2
    // eigenvalues trivial to get

    double mag2 = _pv[0].vect().mag2() / MeV2;
    double q1 = _pv[0].vect().x() / MeV;
    double q2 = _pv[0].vect().y() / MeV;
    double q3 = _pv[0].vect().z() / MeV;

    lambda.resize(3);

    lambda[0] = sqr(q1)/(2.*mag2);
    lambda[1] = sqr(q2)/(2.*mag2);
    lambda[2] = sqr(q3)/(2.*mag2);

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


Axis EventShapes2::eigenvector(const double T[3][3], const double &lam) {
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


vector<Axis> EventShapes2::
eigenvectors(const double T[3][3], const vector<double> &lam) {
  vector<Axis> n;
  for (unsigned int i=0; i<3; i++) {
    n.push_back(eigenvector(T, lam[i]));
  }
  return n;
}

void EventShapes2::diagonalizeTensors(bool linear, bool cmboost) {  
  // initialize
  double Theta[3][3];
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
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
    if (pvec.mag() > 0*MeV) {
      sumvec += pvec;
      if (linear) {
	sum += pvec.mag()*UnitRemoval::InvE;
      } else {
	sum += pvec.mag2()*UnitRemoval::InvE2;
      }
      for(int i=0; i<3; i++) {
	for(int j=i; j<3; j++) {
	  if (linear) {
	    Theta[i][j] += (pvec_MeV[i])*(pvec_MeV[j])*MeV/(pvec.mag());      
	  } else {
	    Theta[i][j] += (pvec_MeV[i])*(pvec_MeV[j]);
	  }
	}
      }
    }
  }
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
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

void EventShapes2::calculateThrust() { 
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
  vector<Momentum3> p;
  Energy psum = 0.0*MeV;
  for(unsigned int l=0; l<_pv.size(); l++) 
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
      Energy eval = 0.*MeV;
      axis = _thrustAxis[0].cross(_thrustAxis[1]);
      _thrustAxis.push_back(axis); 
      for (unsigned int l=0; l<_pv.size(); l++) 
	eval += abs(axis*_pv[l].vect());
      _thrust.push_back(eval/psum);
    } 
  else 
    {
      _thrust.push_back(-1.0);
      _thrustAxis.push_back(Axis()); 
    }
}

void EventShapes2::calcT(const vector<Momentum3> &p, Energy2 &t, Axis &taxis) {
  Energy2 tval;
  t = 0.0*MeV2;
  ThreeVector<Energy2> tv;
  Momentum3 ptot;
  vector<Momentum3> cpm;
  for (unsigned int k=1; k < p.size(); k++) {
    for (unsigned int j=0; j<k; j++) {
      tv = p[j].cross(p[k]);
      ptot = Momentum3();
      for (unsigned int l=0; l<p.size(); l++) {
	if (l!=j && l!=k) {
	  if (p[l]*tv > 0.0*MeV*MeV2) { 
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

void EventShapes2::calcM(const vector<Momentum3> &p, Energy2 &m, Axis &maxis) {
  Energy2 mval;
  m = 0.0 * MeV2;
  Momentum3 tv, ptot;
  vector<Momentum3> cpm;
  for (unsigned int j=0; j < p.size(); j++) {
    tv = p[j];
    ptot = Momentum3();
    for (unsigned int l=0; l<p.size(); l++) {
      if (l!=j) {
	if (p[l]*tv > 0.0*MeV2) { 
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
	 it != cpm.end(); it++) {
      mval = it->mag2();
      if (mval > m) {
	m = mval;
	maxis = it->unit();
      }
    }
  }
}


