// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardTree class.
//

#include "HardTree.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

HardTree::HardTree(vector<HardBranchingPtr> branchings,
		     vector<HardBranchingPtr> spacelike) {
  _branchings.insert(branchings.begin(),branchings.end());
  _spacelike .insert(spacelike .begin(),spacelike .end());
}

void HardBranching::setMomenta(LorentzRotation R,double aparent,
				Lorentz5Momentum ptparent) {
  _original=R*_particle->momentum();
  // compute the shower variables
  Energy2 dot=_n*_p;
  double alpha = (_original*_n)/dot;
  double beta = ((_original*_p)-alpha*sqr(_p.mass()))/dot;
  _qt=_original-alpha*_p-beta*_n-ptparent;
  _pt=sqrt(max(-_qt*_qt,0.*MeV2));
  _z=alpha/aparent;
  // reconstruct children
  for(unsigned int ix=0;ix<_children.size();++ix) {
    _children[ix]->_p=_p;
    _children[ix]->_n=_n;
    _children[ix]->setMomenta(R,alpha,_qt);
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
    if(axis.perp2()>0.) {
      LorentzRotation rot;
      double sinth(sqrt(1.-sqr(axis.z())));
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

