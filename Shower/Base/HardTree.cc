// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardTree class.
//

#include "HardTree.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

HardTree::HardTree(vector<HardBranchingPtr> branchings,
		   vector<HardBranchingPtr> spacelike) 
  : _branchings(branchings.begin(),branchings.end()),
    _spacelike (spacelike .begin(),spacelike .end())
{}

void HardBranching::setMomenta(LorentzRotation R,double aparent,
			       Lorentz5Momentum ptparent,
			       bool setMomentum) {
  if(setMomentum) _original=R*_particle->momentum();
  // compute the shower variables
  Energy2 dot = _n*_p;
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

bool HardTree::connect(ShowerTreePtr shower) {
  _particles.clear();
  // extract the progenitors from the ShowerTree
  vector<ShowerProgenitorPtr> progenitors = shower->extractProgenitors();
  // connect the trees up
  for(set<HardBranchingPtr>::iterator it=branchings().begin();
      it!=branchings().end();++it) {
    Energy2 dmin(1e30*GeV2);
    tShowerParticlePtr partner;   
    for(unsigned int ix=0;ix<progenitors.size();++ix) {
      if((**it).branchingParticle()->id()!=progenitors[ix]->progenitor()->id()) continue;
      if((**it).branchingParticle()->isFinalState()!=
	 progenitors[ix]->progenitor()->isFinalState()) continue;
      Energy2 dtest = 
	sqr(progenitors[ix]->progenitor()->momentum().x()-(**it).showerMomentum().x())+
	sqr(progenitors[ix]->progenitor()->momentum().y()-(**it).showerMomentum().y())+
	sqr(progenitors[ix]->progenitor()->momentum().z()-(**it).showerMomentum().z())+
	sqr(progenitors[ix]->progenitor()->momentum().t()-(**it).showerMomentum().t());
      if(dtest<dmin) {
	partner = progenitors[ix]->progenitor();
	dmin = dtest;
      }
    }
    if(!partner) return false;
    connect(partner,*it);
    if((**it).incoming()) {
      double z((**it).z());
      tHardBranchingPtr parent=(**it).parent();
      while (parent) {
	z *= parent->z();
	parent=parent->parent();
      }
      partner->x(z);
    }
  }
  // return false if not matched
  return particles().size()==progenitors.size();
}

HardBranching::HardBranching(ShowerParticlePtr particle, SudakovPtr sudakov,
			     tHardBranchingPtr parent,bool incoming) 
  : _particle(particle), _incoming(incoming), _parent(parent),
    _sudakov(sudakov)
{}

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

ostream & Herwig::operator<<(ostream & os, const HardTree & x) {
  os << "Output of HardTree " << &x << "\n";
  for(set<HardBranchingPtr>::const_iterator it=x._branchings.begin();
      it!=x._branchings.end();++it) {
    os << "Hard Particle: " << *(**it).branchingParticle() << " has Sudakov " 
       << (**it).sudakov() << "\n";
    os << "It's colour lines are " << (**it).branchingParticle()->colourLine() << "\t" 
       <<  (**it).branchingParticle()->antiColourLine() << "\n";
    for(unsigned int iy=0;iy<(**it).children().size();++iy) {
      os << "\t Children : " << *(**it).children()[iy]->branchingParticle()
	 << "\n";
      os << "It's colour lines are " 
	 << (**it).children()[iy]->branchingParticle()->colourLine() << "\t" 
	 <<  (**it).children()[iy]->branchingParticle()->antiColourLine() << "\n";
    }
  }
  for(set<HardBranchingPtr>::const_iterator it=x._spacelike.begin();
      it!=x._spacelike.end();++it) {
    os << "SpaceLike: " << *(**it).branchingParticle() << " has Sudakov" 
       << (**it).sudakov() << "\n";
    os << "It's colour lines are " 
       << (**it).branchingParticle()->colourLine() << "\t" 
       << (**it).branchingParticle()->antiColourLine() << "\n";
    for(unsigned int iy=0;iy<(**it).children().size();++iy) {
      os << "\t Children: " << *(**it).children()[iy]->branchingParticle()
	 << "\n";
      os << "It's colour lines are " 
	 << (**it).children()[iy]->branchingParticle()->colourLine() << "\t" 
	 << (**it).children()[iy]->branchingParticle()->antiColourLine() << "\n";
    }
  }
  return os;
}
