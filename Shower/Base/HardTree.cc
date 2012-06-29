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
    spacelike_ (spacelike .begin(),spacelike .end())
{}

bool HardTree::connect(ShowerTreePtr shower) {
  particles_.clear();
  // extract the progenitors from the ShowerTree
  vector<ShowerProgenitorPtr> progenitors = shower->extractProgenitors();
  // KMH - 120809 - Added boolean vector to hold on to which progenitors have
  // already been connected to a branching. If connectedProgenitors[ix] = true 
  // it means progenitors[ix] was already connected to something and so it is
  // skipped in the loop over progenitors. This guards against the possiblility
  // of using the same progenitor twice. This can, wrongly, cause events to fail.
  // This was noticed at a rate of around 1 event in 20 for Powheg ZZ production
  // - both Z's would sometimes get associated with the same progenitor. There
  // is still no doubt room for further improvement but using this bool vector
  // already seemed like a good start.
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
  assert( particles().size() == progenitors.size() );
  return true;
}

ostream & Herwig::operator<<(ostream & os, const HardTree & x) {
  os << "Output of HardTree " << &x << "\n";
  for(set<HardBranchingPtr>::const_iterator it=x.branchings_.begin();
      it!=x.branchings_.end();++it) {
    os << "Hard Particle: " << *(**it).branchingParticle() << " has Sudakov " 
       << (**it).sudakov() << " pT = " << (**it).pT()/GeV
       << " scale = " << (**it).scale()/GeV << "\n";
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
