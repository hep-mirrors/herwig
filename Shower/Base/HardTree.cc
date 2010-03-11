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
