// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonTree class.
//

#include "NasonTree.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

NasonTree::NasonTree(vector<NasonBranchingPtr> branchings,
		     vector<NasonBranchingPtr> spacelike) {
  for(unsigned int ix=0;ix<branchings.size();++ix) {
    _branchings.insert(branchings[ix]);
//     cerr << "testing hard " 
// 	 << *branchings[ix]->_particle << " " 
// 	 << branchings[ix]->_sudakov << "\n";
//     for(unsigned int iy=0;iy<branchings[ix]->_children.size();++iy) {
//       cerr << "testing children " << *branchings[ix]->_children[iy]->_particle
// 	   << "\n";
//     }
  }
  for(unsigned int ix=0;ix<spacelike.size();++ix) {
    _spacelike.insert(spacelike[ix]);
//     cerr << "testing spacelike " 
// 	 << *spacelike[ix]->_particle << " " 
// 	 << spacelike[ix]->_sudakov << "\n";
//     for(unsigned int iy=0;iy<spacelike[ix]->_children.size();++iy) {
//       cerr << "testing children " << *spacelike[ix]->_children[iy]->_particle
// 	   << "\n";
//     }
  }
}

void NasonBranching::setMomenta(LorentzRotation R,double aparent,
				Lorentz5Momentum ptparent) {
  _original=R*_particle->momentum();
  // compute the shower variables
  Energy2 dot=_n*_p;
  double alpha = (_original*_n)/dot;
  double beta = ((_original*_p)-alpha*sqr(_p.mass()))/dot;
  _z=alpha/aparent;
  _qt=_original-alpha*_p-beta*_n-ptparent;
  _pt=sqrt(max(-_qt*_qt,0.));
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
    if(_incoming) {
      _scale=_sudakov->calculateScale(z,pt,ids,1);
    }
    else {
      _scale=_sudakov->calculateScale(z,pt,ids,1);
    }
    // get the pt vector
    Lorentz5Momentum vect=_children[0]->_qt;
    Hep3Vector beta_bb = -(_p+ _n).boostVector();
    Lorentz5Momentum p_bb = _p;
    vect.boost(beta_bb);
    p_bb.boost( beta_bb );

    Hep3Vector axis(p_bb.vect().unit());
    if(axis.perp2()>0.) {
      LorentzRotation rot;
      double sinth(sqrt(1.-sqr(axis.z())));
      rot.setRotate(-acos(axis.z()),Hep3Vector(-axis.y()/sinth,axis.x()/sinth,0.));
      vect.transform(rot);
    }
    else if(axis.z()<0.) {
      vect.setPz(vect.pz());
    }
    _phi= atan2(vect.y(),vect.x());
    if(_phi<0.) _phi+=2.*pi;
    if(_children[1]->_incoming) _phi+=pi;
    
  }
}

