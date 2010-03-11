// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardTree class.
//

#include "HardBranching.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

HardBranching::HardBranching(ShowerParticlePtr particle, SudakovPtr sudakov,
			     tHardBranchingPtr parent,bool incoming) 
  : _particle(particle), _incoming(incoming), _parent(parent),
    _sudakov(sudakov)
{}

void HardBranching::setMomenta(LorentzRotation R,double aparent,
			       Lorentz5Momentum ptparent,
			       bool setMomentum) {
  if(setMomentum) _original=R*_particle->momentum();
  // compute the shower variables
  Energy2 dot = _n*_p;
  if(dot==ZERO) return;
  double alpha = (_original*_n)/dot;
  _z=alpha/aparent;
  double beta = ((_original*_p)-alpha*sqr(_p.mass()))/dot;
  _qt=_original-alpha*_p-beta*_n-ptparent;
  _pt=sqrt(max(-_qt*_qt,ZERO));
  // reconstruct children
  for(unsigned int ix=0;ix<_children.size();++ix) {
    _children[ix]->_p=_p;
    _children[ix]->_n=_n;
    _children[ix]->setMomenta(R,alpha,_qt,setMomentum);
  }
  // calculate the evolution scale and phi
  if(!_children.empty()) {
    double z = _children[0]->_z;
    Energy pt = _children[0]->_pt;
    IdList ids(3);
    ids[0]=_particle->id();
    ids[1]=_children[0]->_particle->id();
    ids[2]=_children[1]->_particle->id();
    _scale=_sudakov->calculateScale(z,pt,ids,_incoming ? 1 : 0);
    // get the pt vector
    Lorentz5Momentum vect=_children[0]->_qt;
    Boost beta_bb = -(_p+ _n).boostVector();
    Lorentz5Momentum p_bb = _p;
    vect.boost(beta_bb);
    p_bb.boost( beta_bb );
    Axis axis(p_bb.vect().unit());
    LorentzRotation rot;
    if(axis.perp2()>0.) {
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      rot.setRotate(-acos(axis.z()),
		    Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      vect.transform(rot);
    }
    else if(axis.z()<0.) {
      vect.setZ(vect.z());
    }
    _phi= atan2(vect.y(),vect.x());
    if(_phi<0.)                 _phi+=Constants::twopi;
    if(_children[1]->_incoming) _phi+=Constants::pi;
  }
}

void HardBranching::fixColours() {
  if(!_sudakov) return;
  if(!_incoming&&_children.empty()) return;
  if(_incoming && !_parent) return;
  if(_incoming)
    _sudakov->splittingFn()->
      colourConnection(_parent->_particle,_particle,
		       _parent->children()[1]->_particle,true);
  else
    _sudakov->splittingFn()->
      colourConnection(_particle,_children[0]->_particle,
		       _children[1]->_particle,false);
}
