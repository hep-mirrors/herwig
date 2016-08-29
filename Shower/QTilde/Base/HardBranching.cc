// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardTree class.
//

#include "HardBranching.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void HardBranching::setMomenta(LorentzRotation R,double aparent,
			       Lorentz5Momentum ptparent,
			       bool setMomentum) {
  if(setMomentum) _original=R*_particle->momentum();
  // compute the shower variables
  Energy2 dot = _n*_p;
  if(dot==ZERO) return;
  double alpha = (_original*_n)/dot;
  //set x for incoming partons
  if( status() == HardBranching::Incoming ) x_frac( alpha ); 
  _z=alpha/aparent;
  double beta = ((_original*_p)-alpha*sqr(_p.mass()))/dot;
  _qt = _original - alpha*_p - beta*_n - _z*ptparent;
  _pt=sqrt(max(-_qt*_qt,ZERO));
  // reconstruct children
  for(unsigned int ix=0;ix<_children.size();++ix) {
    _children[ix]->_p=_p;
    _children[ix]->_n=_n;
    _children[ix]->setMomenta( R, alpha, _qt + _z*ptparent, setMomentum);
  }
  // calculate the evolution scale and phi
  if(!_children.empty()) {
    IdList ids(3);
    ids[0]=_particle->dataPtr();
    ids[1]=_children[0]->_particle->dataPtr();
    ids[2]=_children[1]->_particle->dataPtr();
    double z;
    Energy pt;
    Lorentz5Momentum vect;
    if( _status==Outgoing ||
	( (_status==Incoming || _status==Decay ) && 
	  _children[1]->_status == Outgoing ) ) { 
      z  = _children[0]->_z ;
      pt = _children[0]->_pt;
      vect=_children[0]->_qt;
    }
    else {
      z  = _children[1]->_z ;
      pt = _children[1]->_pt;
      swap(ids[1],ids[2]);
      vect=_children[1]->_qt;
    }
    _scale=_sudakov->calculateScale(z,pt,ids,_status);
    // get the pt vector
    if(_status!=Decay) {
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
	vect.setZ( vect.z());
	vect.setY(-vect.y());
      }
      _phi= atan2(vect.y(),vect.x());
      if(_phi<0.) _phi+=Constants::twopi;
    }      
    else {
      const Boost beta_bb = -pVector().boostVector();
      Lorentz5Momentum p_bb = pVector();
      Lorentz5Momentum n_bb = nVector(); 
      p_bb.boost( beta_bb );
      n_bb.boost( beta_bb );
      vect.boost( beta_bb);
      Axis axis(n_bb.vect().unit());
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
      if(_phi<0) _phi+=Constants::twopi;
    }
  }
}

HardBranching::HardBranching(ShowerParticlePtr particle, SudakovPtr sudakov,
			     tHardBranchingPtr parent,Status status) 
  : _particle(particle), _original(), _p(particle->momentum()), _n(), _qt(),
    _shower(particle->momentum()), _pt(ZERO), _x_frac(0.),
    _status(status), _scale(ZERO), _z(0.),_phi(0.), _parent(parent),
    _sudakov(sudakov), type_(ShowerPartnerType::Undefined)
{}
