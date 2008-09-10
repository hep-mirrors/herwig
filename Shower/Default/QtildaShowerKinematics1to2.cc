// -*- C++ -*-
//
// QtildaShowerKinematics1to2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

vector<Lorentz5Momentum> QtildaShowerKinematics1to2::getBasis() const {
  vector<Lorentz5Momentum> dum;
  dum.push_back( _pVector );
  dum.push_back( _nVector );
  return dum; 
}

void QtildaShowerKinematics1to2::setBasis(const Lorentz5Momentum &p,
					  const Lorentz5Momentum & n) {
  _pVector=p;
  _nVector=n;
}

Lorentz5Momentum QtildaShowerKinematics1to2::
sudakov2Momentum(double alpha, double beta, Energy px, 
		 Energy py,unsigned int iopt) const {
  if(isnan(beta)) 
    throw Exception() << "beta infinite in "
		      << "QtildaShowerKinematics1to2::sudakov2Momentum()"
		      << Exception::eventerror;
  Lorentz5Momentum dq;
  if(iopt==0) {
    const Boost beta_bb = -(_pVector + _nVector).boostVector();
    Lorentz5Momentum p_bb = _pVector;
    Lorentz5Momentum n_bb = _nVector; 
    p_bb.boost( beta_bb );
    n_bb.boost( beta_bb );
    // set first in b2b frame along z-axis (assuming that p and n are
    // b2b as checked above)
    dq=Lorentz5Momentum(0.0*MeV, 0.0*MeV, (alpha - beta)*p_bb.vect().mag(), 
			alpha*p_bb.t() + beta*n_bb.t());
    // add transverse components
    dq.setX(px);
    dq.setY(py);
    // rotate to have z-axis parallel to p
    // this rotation changed by PR to a different rotation with the same effect
    // but different azimuthal angle to make implementing spin correlations easier
    //    dq.rotateUz( unitVector(p_bb.vect()) );
    Axis axis(p_bb.vect().unit());
    if(axis.perp2()>0.) {
      LorentzRotation rot;
      double sinth(sqrt(1.-sqr(axis.z())));
      rot.setRotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      dq.transform(rot);
    }
    else if(axis.z()<0.) {
      dq.setZ(-dq.z());
    }
    // boost back 
    dq.boost( -beta_bb ); 
    dq.rescaleMass(); 
    // return the momentum
  }
  else {
    const Boost beta_bb = -pVector().boostVector();
    Lorentz5Momentum p_bb = pVector();
    Lorentz5Momentum n_bb = nVector(); 
    p_bb.boost( beta_bb );
    n_bb.boost( beta_bb );
    // set first in b2b frame along z-axis (assuming that p and n are
    // b2b as checked above)
    dq=Lorentz5Momentum (0.0*MeV, 0.0*MeV, 0.5*beta*pVector().mass(), 
			 alpha*pVector().mass() + 0.5*beta*pVector().mass());
    // add transverse components
    dq.setX(px);
    dq.setY(py);
    // rotate to have z-axis parallel to n
    dq.rotateUz( unitVector(n_bb.vect()) );
    // boost back 
    dq.boost( -beta_bb ); 
    dq.rescaleMass();
  }
  return dq; 
}
