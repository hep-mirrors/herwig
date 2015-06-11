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
		   vector<HardBranchingPtr> spacelike,
		   ShowerInteraction::Type type) 
  : interaction_(type),
    branchings_(branchings.begin(),branchings.end()),
    spacelike_ (spacelike .begin(),spacelike .end()),
    partnersSet_(false)
{}

bool HardTree::connect(ShowerTreePtr shower) {
  particles_.clear();
  // extract the progenitors from the ShowerTree
  vector<ShowerProgenitorPtr> progenitors = shower->extractProgenitors();
  vector<bool> connectedProgenitors(progenitors.size(),false);
  // connect the trees up
  for( set<HardBranchingPtr>::iterator it = branchings().begin();
       it != branchings().end(); ++it) {
    Energy2 dmin( 1e30*GeV2 );
    tShowerParticlePtr partner;   
    unsigned int progenitorsIndex(999);
    for( unsigned int ix = 0; ix < progenitors.size(); ++ix ) {
      if( connectedProgenitors[ix] ) continue;
      if( (**it).branchingParticle()->id() != progenitors[ix]->progenitor()->id() ) continue;
      if( (**it).branchingParticle()->isFinalState() !=
	  progenitors[ix]->progenitor()->isFinalState() ) continue;
      Energy2 dtest = 
	sqr( progenitors[ix]->progenitor()->momentum().x() - (**it).showerMomentum().x() ) +
	sqr( progenitors[ix]->progenitor()->momentum().y() - (**it).showerMomentum().y() ) +
	sqr( progenitors[ix]->progenitor()->momentum().z() - (**it).showerMomentum().z() ) +
	sqr( progenitors[ix]->progenitor()->momentum().t() - (**it).showerMomentum().t() );
      if( dtest < dmin ) {
	partner = progenitors[ix]->progenitor();
	progenitorsIndex = ix;
	dmin = dtest;
      }
    }
    if( !partner ) return false;
    connectedProgenitors[progenitorsIndex] = true;
    connect( partner, *it );
    if( (**it).status() == HardBranching::Incoming ) {
      double z( (**it).z() );
      tHardBranchingPtr parent = (**it).parent();
      while (parent) {
	z *= parent->z();
	parent = parent->parent();
      }
      partner->x(z);
    }
  }
  if( particles().size() == progenitors.size() ) return  true;
  else                                           return false;
}

ostream & Herwig::operator<<(ostream & os, const HardTree & x) {
  os << "Output of HardTree " << &x << "\n";
  for(set<HardBranchingPtr>::const_iterator it=x.branchings_.begin();
      it!=x.branchings_.end();++it) {
    os << "Hard Particle: " << *(**it).branchingParticle() << " has Sudakov " 
       << (**it).sudakov() << " pT = " << (**it).pT()/GeV
       << " scale = " << (**it).scale()/GeV << "\n";
    os << "Its colour lines are " << (**it).branchingParticle()->colourLine() << "\t" 
       <<  (**it).branchingParticle()->antiColourLine() << "\n";
    os << "Its basis vectors are " << (**it).pVector()/GeV 
       << " " << (**it).nVector()/GeV << "\n";
    os << "Its shower momentum is " << (**it).showerMomentum()/GeV << "\n";
    for(unsigned int iy=0;iy<(**it).children().size();++iy) {
      os << "\t Children : " << *(**it).children()[iy]->branchingParticle()
	 << "\n";
      os << "It's colour lines are " 
	 << (**it).children()[iy]->branchingParticle()->colourLine() << "\t" 
	 <<  (**it).children()[iy]->branchingParticle()->antiColourLine() << "\n";
    }
  }
  for(set<HardBranchingPtr>::const_iterator it=x.spacelike_.begin();
      it!=x.spacelike_.end();++it) {
    os << "SpaceLike: " << *(**it).branchingParticle() << " has Sudakov" 
       << (**it).sudakov() << " pT = " << (**it).pT()/GeV
       << " scale = " << (**it).scale()/GeV << "\n";
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
