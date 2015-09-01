// -*- C++ -*-
//
// QTildeShowerKinematics1to2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeShowerKinematics1to2 class.
//

#include "QTildeShowerKinematics1to2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
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
  Boost beta_bb;
  if(frame()==BackToBack) {
    beta_bb = -(_pVector + _nVector).boostVector();
  }
  else if(frame()==Rest) {
    beta_bb = -pVector().boostVector();
  }
  else
    assert(false);
  Lorentz5Momentum p_bb = pVector();
  Lorentz5Momentum n_bb = nVector(); 
  p_bb.boost( beta_bb );
  n_bb.boost( beta_bb );
  // rotate to have z-axis parallel to p/n
  Axis axis;
  if(frame()==BackToBack) {
    axis = p_bb.vect().unit();
  }
  else if(frame()==Rest) {
    axis = n_bb.vect().unit();
  }
  else
    assert(false);
  LorentzRotation rot;
  if(axis.perp2()>1e-10) {
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  }
  else if(axis.z()<0.) {
    rot.rotate(Constants::pi,Axis(1.,0.,0.));
  }
  _xPerp=LorentzVector<double>(1.,0.,0.,0.);
  _yPerp=LorentzVector<double>(0.,1.,0.,0.);
  _xPerp.transform(rot);
  _yPerp.transform(rot);
  // boost back 
  _xPerp.boost( -beta_bb );
  _yPerp.boost( -beta_bb );
}

void QTildeShowerKinematics1to2::setMomentum(tShowerParticlePtr particle,
					     bool timeLike) const {
  Energy mass = particle->mass() > ZERO ? particle->mass() : particle->data().mass();
  // calculate the momentum of the assuming on-shell
  Energy2 pt2 = sqr(particle->showerParameters().pt);
  double alpha = timeLike ? particle->showerParameters().alpha : particle->x();
  double beta = 0.5*(sqr(mass) + pt2 - sqr(alpha)*pVector().m2())/(alpha*p_dot_n());
  Lorentz5Momentum porig=sudakov2Momentum(alpha,beta,
					  particle->showerParameters().ptx,
					  particle->showerParameters().pty);
  porig.setMass(mass);
  particle->set5Momentum(porig);
}

void QTildeShowerKinematics1to2::constructSpinInfo(tShowerParticlePtr particle,
						   bool timeLike) const {
  // now construct the required spininfo and calculate the basis states
  PDT::Spin spin(particle->dataPtr()->iSpin());
  if(spin==PDT::Spin0) {
    ScalarWaveFunction::constructSpinInfo(particle,outgoing,timeLike);
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
}
void QTildeShowerKinematics1to2::transform(const LorentzRotation & r) {
  _pVector *= r;
  _nVector *= r;
  _xPerp   *= r;
  _yPerp   *= r;
}
