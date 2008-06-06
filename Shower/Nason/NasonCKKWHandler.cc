// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonCKKWHandler class.
//

#include "NasonCKKWHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "QTildeSudakovIntegrator.h"
#include "NasonTree.h"


using namespace Herwig;

NasonCKKWHandler::~NasonCKKWHandler() {}

IBPtr NasonCKKWHandler::clone() const {
  return new_ptr(*this);
}

IBPtr NasonCKKWHandler::fullclone() const {
  return new_ptr(*this);
}

void NasonCKKWHandler::persistentOutput(PersistentOStream & os) const {
}

void NasonCKKWHandler::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<NasonCKKWHandler> NasonCKKWHandler::initNasonCKKWHandler;
// Definition of the static class description member.

void NasonCKKWHandler::Init() {

  static ClassDocumentation<NasonCKKWHandler> documentation
    ("There is no documentation for the NasonCKKWHandler class");

}

double NasonCKKWHandler::reweightCKKW(int minMult, int maxMult) {
  PPair in = lastXCombPtr()->subProcess()->incoming();
  ParticleVector out  = lastXCombPtr()->subProcess()->outgoing();
  PPtr vBoson = lastXCombPtr()->subProcess()->intermediates()[0];
  _s = lastXCombPtr()->lastS();
  // cluster the event into a NasonTree
  NasonTreePtr nasonTree = doClustering( out, vBoson );
  if(!nasonTree) return 0.;
  // compute the CKKW weight
  //
  //   need to implement this here
  //
  // change the hard process
  generator()->log() << *lastXCombPtr()->subProcess() << "\n";
  generator()->log() << "testing got to the exit!!!!!!! " << flush;
  exit(0);
  return 1.;
}

void NasonCKKWHandler::doinitrun() {
  ShowerHandler::doinitrun();
//   // integrator for the outer integral
//   GaussianIntegrator outer;
//   // get the final-state branchings from the evolver
//   ofstream output("test.top");
//   output << "SET FONT DUPLEX\n";
//   for(BranchingList::const_iterator 
//         it = evolver()->splittingGenerator()->finalStateBranchings().begin();
//         it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
//     Ptr<QTildeSudakovIntegrator>::pointer integrator = 
//       new_ptr(QTildeSudakovIntegrator(it->second));
//     cerr << "testing sudakov " << it->second.first->fullName() << "\t"
// 	 << it->second.second[0] << "\t"
// 	 << it->second.second[1] << "\t"
// 	 << it->second.second[2] << "\n";
//     Energy qtildemax=generator()->maximumCMEnergy();
//     Energy qtildemin=integrator->minimumScale();
//     vector<double> sud;
//     vector<Energy> scale;
//     sud.push_back(0.); scale.push_back(qtildemin);
//     Energy currentScale=qtildemin;
//     double fact = pow(qtildemax/qtildemin,1./(_npoint-1));
//     for(unsigned int ix=1;ix<_npoint;++ix) {
//       currentScale *= fact;
//       double currentSud = integrator->value(currentScale,scale.back());
//       scale.push_back(currentScale);
//       sud.push_back(sud.back()+currentSud);
//       cerr << "testing values " << scale.back()/GeV << "\t" << sud.back() << " " << exp(-sud.back()) << "\n";
//     }
//     // convert to the Sudakov
//     for(unsigned int ix=0;ix<sud.size();++ix) {
//       sud[ix] = exp(-sud[ix]);
//     }
//     // construct the Interpolators
//     Interpolator<double,Energy>::Ptr intq = new_ptr(Interpolator<double,Energy>(sud,scale,3));
//     Interpolator<Energy,double>::Ptr ints = new_ptr(Interpolator<Energy,double>(scale,sud,3));
//     _fbranchings.insert(make_pair(it->first,make_pair(intq,ints)));




//     output << "NEWFRAME\n";
//     output << "TITLE TOP \"Sudakov for " << getParticleData(it->second.second[0])->PDGName() << " -> "
// 	   << getParticleData(it->second.second[1])->PDGName() << " "
// 	   << getParticleData(it->second.second[2])->PDGName() << "\"\n";
//     for(unsigned int ix=0;ix<sud.size();++ix)
//       output << scale[ix]/GeV << " " << sud[ix] << "\n";
//     output << "JOIN RED\n" << flush;
// //     HistogramPtr temp(new_ptr(Histogram(0.,100.,200)));

// //     double slst = (*intq)(91.2*GeV);
// //     for(unsigned int ix=0;ix<100000000;++ix) {
// //       double snow = slst/UseRandom::rnd();
// //       if(snow>=1.) continue;
// //       Energy qnow = (*ints)(snow);
// //       *temp +=qnow/GeV;
// //     }
// //     using namespace HistogramOptions;
// //     temp->topdrawOutput(output,Frame);
//   }
}

void NasonCKKWHandler::cascade() {
  ShowerHandler::cascade();
}

double NasonCKKWHandler::getJetMeasure(ShowerParticlePtr part_i,
				       ShowerParticlePtr part_j){
  double yij;
  double costheta = part_i->momentum().vect().dot( part_j->momentum().vect() ) 
    / part_i->momentum().vect().mag() / part_j->momentum().vect().mag();
  if( sqr( part_i->momentum().e() ) > sqr( part_j->momentum().e() ) )
    yij = 2. * sqr( part_j->momentum().e() ) * ( 1. - costheta ) / _s ;
  else
    yij = 2. * sqr( part_i->momentum().e() ) * ( 1. - costheta ) / _s ;
  return yij;
}

//given two particles returns value of durham jet algorithm
bool NasonCKKWHandler:: splittingAllowed( ShowerParticlePtr part_i,
					  ShowerParticlePtr part_j,
					  int  qq_pairs ) {
  // g 2 q qbar or an incorrect qq type
  if ( abs ( part_i->id() ) < 7 && abs ( part_j->id() ) < 7 ) { 
    if ( abs ( part_i->id() ) != abs ( part_j->id() )  ) return false;
    if ( ( part_i->id() < 0 &&  part_j->id() < 0 ) ||
	 ( part_i->id() > 0 &&  part_j->id() > 0 ) ) return false;
    if ( qq_pairs < 2 ) return false;
  }
  
  return true;
}

// finds the id of the emitting particle and sudakov for the desired clustering

SudakovPtr NasonCKKWHandler:: getSud( int & qq_pairs, long & emmitter_id,
				      ShowerParticlePtr part_i, 
				      ShowerParticlePtr part_j ) {
  // g 2 q qbar or an incorrect qq type
  if ( abs ( part_i->id() ) < 7 && abs ( part_j->id() ) < 7 ) { 
    if ( abs ( part_i->id() ) != abs ( part_j->id() )  ) return SudakovPtr();
    if ( ( part_i->id() < 0 &&  part_j->id() < 0 ) ||
	 ( part_i->id() > 0 &&  part_j->id() > 0 ) ) return SudakovPtr();
    if ( qq_pairs < 2 ) return SudakovPtr();
    qq_pairs -= 1;
    emmitter_id = 21;
  }
  // q/qbar 2 q/qbar g
  else if ( abs ( part_i->id() ) < 7 || abs ( part_j->id() ) < 7 ) {
    if( abs ( part_i->id() ) < 7 ){
      emmitter_id = part_i->id();
    }
    else {
      emmitter_id = part_j->id();
    }
  }
  // g 2 g g
  else {
    emmitter_id = 21;
  }

  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();

  //cycle through list of all branchings with the correct abs ( emmitter_id )
  for(BranchingList::const_iterator cit = branchings.lower_bound( abs(emmitter_id) );
      cit != branchings.upper_bound( abs(emmitter_id) ); ++cit ) {
    IdList ids = cit->second.second;
    if( abs( ids[0] ) == abs( emmitter_id ) 
	&& ( ( abs(ids[1]) == abs( part_i->id() ) 
	       && abs( ids[2] ) == abs( part_j->id() ) ) ||
	     ( abs( ids[1] ) == abs( part_j->id() ) 
	       && abs( ids[2] ) == abs( part_i->id() ) ) ) ) {
      return cit->second.first;
    }
  }
  return SudakovPtr();
}


NasonTreePtr NasonCKKWHandler::doClustering( ParticleVector theParts, 
					     PPtr vb ) {
  int qq_pairs = 0;
  map <ShowerParticlePtr,NasonBranchingPtr> theParticles;
  tcPDPtr particle_data;
  ShowerParticlePtr vBoson = new_ptr( ShowerParticle( *vb, 1, false, false ) );
  //loops through the FS particles and create naon branchings
  for( unsigned int i = 0; i < theParts.size(); i++){
    ShowerParticlePtr currentParticle = 
      new_ptr( ShowerParticle( *theParts[i], 1, true, false ) );
    //    currentParticle->rescaleMass();
    if( currentParticle->id() > 0 && currentParticle->id() < 7 ) qq_pairs++;
    theParticles.insert(make_pair(currentParticle, 
				  new_ptr( NasonBranching( currentParticle, SudakovPtr(),
							   NasonBranchingPtr(),false ) ) ) );  
    if(currentParticle->dataPtr()->iColour()==PDT::Colour3||
       currentParticle->dataPtr()->iColour()==PDT::Colour8) {
      ColinePtr newline = new_ptr(ColourLine());
      newline->addColoured(currentParticle);
    }
    if(currentParticle->dataPtr()->iColour()==PDT::Colour3bar||
       currentParticle->dataPtr()->iColour()==PDT::Colour8) {
      ColinePtr newline = new_ptr(ColourLine());
      newline->addAntiColoured(currentParticle);
    }
  }
  // loops clustering until we get down to qqbar
  while( theParticles.size() > 2 ){
    double yij_min = 1.;
    pair< ShowerParticlePtr, ShowerParticlePtr > clusterPair;
    //loops over all pairs of particles in theParticles
    for( map<ShowerParticlePtr, NasonBranchingPtr>::iterator ita = theParticles.begin();
	 ita != theParticles.end() ; ita++ ) {
      for( map<ShowerParticlePtr, NasonBranchingPtr>::iterator itb = theParticles.begin();
	   itb != ita; itb++) {
	double yij = getJetMeasure( ita->first, itb->first );
	if( yij < yij_min && splittingAllowed( ita->first, itb->first, qq_pairs )  ) {
	  clusterPair.first  = ita->first;
	  clusterPair.second = itb->first;
	  yij_min = yij;
	}	
      }
    }
    long thePartId;
    SudakovPtr theSudakov=getSud( qq_pairs, thePartId,
				  clusterPair.first, clusterPair.second ); 
    if( !theSudakov ){
      cerr << "can't find the sudakov in: \n";
      cerr << *clusterPair.first<<"\n"
	   << *clusterPair.second<<"\n";
      cerr << "with qq_pairs = " << qq_pairs <<"\n";
    }
    generator()->log() << "testing sudakov?? " << theSudakov << "\n";
    Lorentz5Momentum pairMomentum = clusterPair.first->momentum() + 
      clusterPair.second->momentum();
    pairMomentum.setMass(0.*MeV);
    particle_data = getParticleData( thePartId );
    
    //creates emitter particle
    ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, true ) );
    clustered->set5Momentum( pairMomentum );

    NasonBranchingPtr clusteredBranch(new_ptr(NasonBranching( clustered, theSudakov,
							      NasonBranchingPtr(),false)));
    fixColours(clustered,clusterPair.first,clusterPair.second);
    theParticles.insert( make_pair( clustered, clusteredBranch ) );
    clusteredBranch->addChild( theParticles.find( clusterPair.first )->second  );
    clusteredBranch->addChild( theParticles.find( clusterPair.second )->second );
    theParticles.erase( clusterPair.first );
    theParticles.erase( clusterPair.second );
  
  }
  vector<NasonBranchingPtr> theBranchings;
  for(  map<ShowerParticlePtr, NasonBranchingPtr>::iterator it = theParticles.begin(); 
	it != theParticles.end(); ++it){ 
    theBranchings.push_back( it->second );
  }
  // fix for e+e- to match up the colours of the q qbar pair
  if(theBranchings[0]->_particle->dataPtr()->iColour()==PDT::Colour3) {
    ColinePtr temp = theBranchings[1]->_particle->antiColourLine();
    theBranchings[0]->_particle->colourLine()->join(temp);
  }
  else {
    ColinePtr temp = theBranchings[0]->_particle->antiColourLine();
    theBranchings[1]->_particle->colourLine()->join(temp);
  }
  vector<NasonBranchingPtr> spaceBranchings;
  spaceBranchings.push_back( new_ptr( NasonBranching( vBoson, SudakovPtr(),
						      NasonBranchingPtr(), 
						      true ) ) );
  theBranchings.push_back( spaceBranchings.back() );
  NasonTreePtr nasontree = new_ptr( NasonTree( theBranchings,
					       spaceBranchings ) );

  //should I connect the branchings to the particles
  //doesn't do anything unless showerParticles come from showertree
  for(  map<ShowerParticlePtr, NasonBranchingPtr>::iterator it = theParticles.begin(); 
	it != theParticles.end(); ++it){ 
    nasontree->connect( it->first, it->second );
  }
  generator()->log() << *nasontree << "\n" << flush;

  // Calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->
    reconstructDecayShower(nasontree,evolver());
  generator()->log() << "testing hard momenta for the shower\n";
  for(unsigned int ix=0;ix<theBranchings.size();++ix) {
    generator()->log() << "testing " << theBranchings[ix]->_shower/GeV << "\t"
		       << theBranchings[ix]->_shower.m()/GeV << "\n";
  }
  return nasontree;
}

void NasonCKKWHandler::fixColours(tPPtr parent, tPPtr child1, tPPtr child2) {
  // the different possible cases
  if(parent->dataPtr()->iColour()==PDT::Colour3&&
     child1->dataPtr()->iColour()==PDT::Colour3&&
     child2->dataPtr()->iColour()==PDT::Colour8) {
    child2->colourLine()->addColoured(parent);
    ColinePtr temp = child2->antiColourLine();
    child1->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3&&
	  child2->dataPtr()->iColour()==PDT::Colour3&&
	  child1->dataPtr()->iColour()==PDT::Colour8) {
    child1->colourLine()->addColoured(parent);
    ColinePtr temp = child1->antiColourLine();
    child2->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3bar&&
	  child1->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour8) {
    child2->antiColourLine()->addAntiColoured(parent);
    ColinePtr temp = child1->antiColourLine();
    child2->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour3bar&&
	  child1->dataPtr()->iColour()==PDT::Colour8) {
    child1->antiColourLine()->addAntiColoured(parent);
    ColinePtr temp = child2->antiColourLine();
    child1->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour8&&
	  child1->dataPtr()->iColour()==PDT::Colour8&&
	  child2->dataPtr()->iColour()==PDT::Colour8) {
    if(UseRandom::rndbool(0.5)) {
      child1->colourLine()->addColoured(parent);
      child2->antiColourLine()->addAntiColoured(parent);
      ColinePtr temp = child1->antiColourLine();
      child2->colourLine()->join(temp);
    }
    else {
      child2->colourLine()->addColoured(parent);
      child1->antiColourLine()->addAntiColoured(parent);
      ColinePtr temp = child2->antiColourLine();
      child1->colourLine()->join(temp);
    }
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour8&&
	  child1->dataPtr()->iColour()==PDT::Colour3&&
	  child2->dataPtr()->iColour()==PDT::Colour3bar) {
    child1->colourLine()->addColoured(parent);
    child2->antiColourLine()->addAntiColoured(parent);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour8&&
	  child1->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour3) {
    child2->colourLine()->addColoured(parent);
    child1->antiColourLine()->addAntiColoured(parent);
  }  
  else {
    generator()->log() << "testing in fixcolurs " << *parent << "\n" << *child1 << "\n"
		       << *child2 << "\n" << flush;
    throw Exception() << "Unknown colour in NasonCKKWHandler::fixColours()"
		      << Exception::runerror;
  }
}
