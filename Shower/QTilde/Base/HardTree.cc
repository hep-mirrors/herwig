// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardTree class.
//

#include "HardTree.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "Herwig/Shower/RealEmissionProcess.h"

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

HardTree::HardTree(RealEmissionProcessPtr real) 
  : interaction_(real->interaction()),
    partnersSet_(false) {
  // check the type
  unsigned int pType = real->incoming().size() == 2 ? 1 : 2;
  // create the Branchings for the incoming particles
  vector<HardBranchingPtr> spaceBranchings;
  ShowerParticleVector incoming;
  for(unsigned int ix=0;ix<real->incoming().size();++ix) {
    unsigned int type = ix!=real->emitter() ? pType : 0;
    ShowerParticlePtr part(new_ptr(ShowerParticle(*real->incoming()[ix],type,false)));
    incoming.push_back(part);
    spaceBranchings.push_back(new_ptr(HardBranching(part,SudakovPtr(),HardBranchingPtr(),
						    HardBranching::Incoming)));
  }
  // create the outgoing particles
  ShowerParticleVector outgoing;
  vector<HardBranchingPtr> timeBranchings;
  for(unsigned int ix=0;ix<real->outgoing().size();++ix) {
    unsigned int type = 
      int(ix)!=int(real->emitted()-real->incoming().size()) &&
      int(ix)!=int(real->emitter()-real->incoming().size()) ? pType : 0;
    ShowerParticlePtr part(new_ptr(ShowerParticle(*real->outgoing()[ix],type,true)));
    outgoing.push_back(part);
    timeBranchings.push_back(new_ptr(HardBranching(part,SudakovPtr(),HardBranchingPtr(),
						   HardBranching::Outgoing)));
  }
  // create the branching for the radiating particle
  HardBranchingPtr emitterBranch;
  // II
  if(real->emitter()   < real->incoming().size() &&
     real->spectator() < real->incoming().size()) {
    assert(false);
  }
  // FF
  else if(real->emitter()   >= real->incoming().size() &&
	  real->spectator() >= real->incoming().size()) {
    unsigned int iemitter = real->emitter()  -real->incoming().size();
    unsigned int iemitted = real->emitted()  -real->incoming().size();
    unsigned int ispect   = real->spectator()-real->incoming().size();
    ShowerParticlePtr parent(new_ptr(ShowerParticle(real->bornOutgoing()[iemitter]->dataPtr(),true)));
    Lorentz5Momentum parentMomentum = 
      real->outgoing()[iemitter]->momentum() +
      real->outgoing()[iemitted]->momentum();
    parentMomentum.setMass(real->bornOutgoing()[iemitter]->dataPtr()->mass());
    parent->set5Momentum(parentMomentum);
    emitterBranch = new_ptr(HardBranching(parent,SudakovPtr(),HardBranchingPtr(),
					  HardBranching::Outgoing));
    // QED radiation
    if(real->interaction()==ShowerInteraction::QED) {
      emitterBranch->type(ShowerPartnerType::QED);
      if(parent->id()!=ParticleID::gamma) {
	if(timeBranchings[iemitter]->branchingParticle()->    colourLine())
	  timeBranchings[iemitter]->branchingParticle()->    colourLine()->addColoured(parent);
	if(timeBranchings[iemitter]->branchingParticle()->antiColourLine())
	  timeBranchings[iemitter]->branchingParticle()->antiColourLine()->addAntiColoured(parent);
      }
    }
    else {
      if(real->outgoing()[iemitted]->dataPtr()->iColour()==PDT::Colour8) {
	// g -> g g 
	if(real->outgoing()[iemitter]->dataPtr()->iColour()==PDT::Colour8) {
	  cerr << "g to gg\n";
	  assert(false);
	}
	// q -> q g 
	else if(real->outgoing()[iemitter]->dataPtr()->iColour()==PDT::Colour3) {
	  timeBranchings[iemitted]->branchingParticle()->colourLine()->addColoured(parent);
      emitterBranch->type(ShowerPartnerType::QCDColourLine);
	}
	// qbar -> qbar g
	else if(real->outgoing()[iemitter]->dataPtr()->iColour()==PDT::Colour3bar) {
	  timeBranchings[iemitted]->branchingParticle()->antiColourLine()->addAntiColoured(parent);
	  emitterBranch->type(ShowerPartnerType::QCDAntiColourLine);
	}
      }
      else
	assert(false);
    }
    emitterBranch->addChild(timeBranchings[iemitter]);
    emitterBranch->addChild(timeBranchings[iemitted]);
    // set partners
    emitterBranch         ->colourPartner(timeBranchings[ispect]);
    timeBranchings[ispect]->colourPartner(emitterBranch);
    emitterBranch         ->branchingParticle()->partner(timeBranchings[ispect]->branchingParticle());
    timeBranchings[ispect]->branchingParticle()->partner(emitterBranch         ->branchingParticle());
  }
  // IF
  else {
    assert(false);
  }
  assert(emitterBranch);
  spacelike_=set<HardBranchingPtr>(spaceBranchings.begin(),spaceBranchings.end());
  // insert the incoming branchings
  for(unsigned int ix=0;ix<spaceBranchings.size();++ix) {
    HardBranchingPtr inBranch;
    if(real->emitter()<real->incoming().size() && real->emitter()==ix)
      inBranch = emitterBranch;
    else
      inBranch = spaceBranchings[ix];
    inBranch->beam(real->hadrons()[ix]);
    branchings_.insert(inBranch);
  }
  // and the outgoing ones
  for(unsigned int ix=0;ix<timeBranchings.size();++ix) {
    if(real->emitter()>=real->incoming().size() && 
       real->emitter()==ix+real->incoming().size())
      branchings_.insert(emitterBranch);
    else if(real->emitted()==ix+real->incoming().size())
      continue;
    else
      branchings_.insert(timeBranchings[ix]);
  }
}

  // incoming progenitors
//   ShowerProgenitorPtr 
//     ePProgenitor = tree->incomingLines().begin() ->first,
//     eMProgenitor = tree->incomingLines().rbegin()->first;
//   if(eMProgenitor->id()<0) swap(eMProgenitor,ePProgenitor);
//   // outgoing progenitors
//   ShowerProgenitorPtr 
//     qkProgenitor = tree->outgoingLines().begin() ->first,
//     qbProgenitor = tree->outgoingLines().rbegin()->first;

//   if(force==ShowerInteraction::QED) {
//     spaceBranchings[0]->branchingParticle()->partner(spaceBranchings[1]->branchingParticle());
//     spaceBranchings[1]->branchingParticle()->partner(spaceBranchings[0]->branchingParticle());
//   }
//   // Make the HardTree from the HardBranching vectors.
//   HardTreePtr hardtree = new_ptr(HardTree(allBranchings,spaceBranchings,
// 					  force));
//   hardtree->partnersSet(true);
//   // Connect the particles with the branchings in the HardTree
//   hardtree->connect( eMProgenitor->progenitor(), allBranchings[0] );
//   hardtree->connect( ePProgenitor->progenitor(), allBranchings[1] );
//   hardtree->connect( qkProgenitor->progenitor(), allBranchings[2] );
//   hardtree->connect( qbProgenitor->progenitor(), allBranchings[3] );
