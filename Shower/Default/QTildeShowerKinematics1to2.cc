// -*- C++ -*-
//
// QTildeShowerKinematics1to2.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeShowerKinematics1to2 class.
//

#include "QTildeShowerKinematics1to2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

vector<Lorentz5Momentum> QTildeShowerKinematics1to2::getBasis() const {
  vector<Lorentz5Momentum> dum;
  dum.push_back( _pVector );
  dum.push_back( _nVector );
  return dum; 
}

void QTildeShowerKinematics1to2::setBasis(const Lorentz5Momentum &p,
					  const Lorentz5Momentum & n,
					  Frame inframe) {
  _pVector=p;
  _nVector=n;
  frame(inframe);
}

Lorentz5Momentum QTildeShowerKinematics1to2::
sudakov2Momentum(double alpha, double beta, Energy px, Energy py) const {
  if(isnan(beta)||isinf(beta)) 
    throw Exception() << "beta infinite in "
		      << "QTildeShowerKinematics1to2::sudakov2Momentum()"
		      << Exception::eventerror;
  Lorentz5Momentum dq;
  if(frame()==BackToBack) {
    const Boost beta_bb = -(_pVector + _nVector).boostVector();
    Lorentz5Momentum p_bb = _pVector;
    Lorentz5Momentum n_bb = _nVector; 
    p_bb.boost( beta_bb );
    n_bb.boost( beta_bb );
    // set first in b2b frame along z-axis (assuming that p and n are
    // b2b as checked above)
    dq=Lorentz5Momentum(ZERO, ZERO, (alpha - beta)*p_bb.vect().mag(), 
			alpha*p_bb.t() + beta*n_bb.t());
    // add transverse components
    dq.setX(px);
    dq.setY(py);
    // rotate to have z-axis parallel to p
    // this rotation changed by PR to a different rotation with the same effect
    // but different azimuthal angle to make implementing spin correlations easier
    //    dq.rotateUz( unitVector(p_bb.vect()) );
    Axis axis(p_bb.vect().unit());
    LorentzRotation rot;
    if(axis.perp2()>0.) {
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    else if(axis.z()<0.) {
      rot.rotate(Constants::pi,Axis(1.,0.,0.));
    }
    dq.transform(rot);
    // boost back 
    dq.boost( -beta_bb ); 
    dq.rescaleMass(); 
    // return the momentum
  }
  else if(frame()==Rest) {
    const Boost beta_bb = -pVector().boostVector();
    Lorentz5Momentum p_bb = pVector();
    Lorentz5Momentum n_bb = nVector(); 
    p_bb.boost( beta_bb );
    n_bb.boost( beta_bb );
    // set first in b2b frame along z-axis (assuming that p and n are
    // b2b as checked above)
    dq=Lorentz5Momentum (ZERO, ZERO, 0.5*beta*pVector().mass(), 
			 alpha*pVector().mass() + 0.5*beta*pVector().mass());
    // add transverse components
    dq.setX(px);
    dq.setY(py);
    // changed to be same as other case
//     dq.rotateUz( unitVector(n_bb.vect()) );
    Axis axis(n_bb.vect().unit());
    LorentzRotation rot;
    if(axis.perp2()>0.) {
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    else if(axis.z()<0.) {
      rot.rotate(Constants::pi,Axis(1.,0.,0.));
    }
    dq.transform(rot);
    // boost back 
    dq.boost( -beta_bb ); 
    dq.rescaleMass();
  }
  else
    assert(false);
  return dq; 
}

void QTildeShowerKinematics1to2::constructSpinInfo(tShowerParticlePtr particle,
						   bool timeLike) const {
  Energy mass = particle->data().mass(); 
  // calculate the momentum of the assuming on-shell
  Energy2 pt2 = sqr(particle->showerParameters().pt);
  double alpha = timeLike ? particle->showerParameters().alpha : particle->x();
  double beta = 0.5*(sqr(mass) + pt2 - sqr(alpha)*pVector().m2())/(alpha*p_dot_n());
  Lorentz5Momentum porig=sudakov2Momentum(alpha,beta,
					  particle->showerParameters().ptx,
					  particle->showerParameters().pty);
  porig.setMass(mass);
  // now construct the required spininfo and calculate the basis states
  PDT::Spin spin(particle->dataPtr()->iSpin());
  if(spin==PDT::Spin0) {
    assert(false);
  }
  // calculate the basis states and construct the SpinInfo for a spin-1/2 particle
  else if(spin==PDT::Spin1Half) {
    // outgoing particle
    if(particle->id()>0) {
      vector<LorentzSpinorBar<SqrtEnergy> > stemp;
      SpinorBarWaveFunction::calculateWaveFunctions(stemp,particle,outgoing);
      SpinorBarWaveFunction::constructSpinInfo(stemp,particle,outgoing,timeLike);
    }
    // outgoing antiparticle
    else {
      vector<LorentzSpinor<SqrtEnergy> > stemp;
      SpinorWaveFunction::calculateWaveFunctions(stemp,particle,outgoing);
      SpinorWaveFunction::constructSpinInfo(stemp,particle,outgoing,timeLike);
    }
  }
  // calculate the basis states and construct the SpinInfo for a spin-1 particle
  else if(spin==PDT::Spin1) {
    bool massless(particle->id()==ParticleID::g||particle->id()==ParticleID::gamma);
    vector<Helicity::LorentzPolarizationVector> vtemp;
    VectorWaveFunction::calculateWaveFunctions(vtemp,particle,outgoing,massless);
    VectorWaveFunction::constructSpinInfo(vtemp,particle,outgoing,timeLike,massless);
  }
  else {
    throw Exception() << "Spins higher than 1 are not yet implemented in " 
		      << "FS_QtildaShowerKinematics1to2::constructVertex() "
		      << Exception::runerror;
  }
  particle->set5Momentum(porig);
}
