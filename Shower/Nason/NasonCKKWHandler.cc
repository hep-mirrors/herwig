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
  
  if( lastXCombPtr()->subProcess()->intermediates().size() > 1 )
    cerr<< "more than one intermediate \n";

  PPtr vBoson = lastXCombPtr()->subProcess()->intermediates()[0];
  _s = lastXCombPtr()->lastS();

  for(unsigned int ix=0;ix<out.size();++ix) {
    generator()->log() << *out[ix] << "\n";
  }
  generator()->log() << *lastXCombPtr()->subProcess() << "\n";

  doClustering( out, vBoson );
  return 1.;
}

void NasonCKKWHandler::doinitrun() {
  ShowerHandler::doinitrun();
  // integrator for the outer integral
  GaussianIntegrator outer;
  // get the final-state branchings from the evolver
  ofstream output("test.top");
  output << "SET FONT DUPLEX\n";
  for(BranchingList::const_iterator 
        it = evolver()->splittingGenerator()->finalStateBranchings().begin();
        it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
    Ptr<QTildeSudakovIntegrator>::pointer integrator = 
      new_ptr(QTildeSudakovIntegrator(it->second));
    cerr << "testing sudakov " << it->second.first->fullName() << "\t"
	 << it->second.second[0] << "\t"
	 << it->second.second[1] << "\t"
	 << it->second.second[2] << "\n";
    Energy qtildemax=generator()->maximumCMEnergy();
    Energy qtildemin=integrator->minimumScale();
    vector<double> sud;
    vector<Energy> scale;
    sud.push_back(0.); scale.push_back(qtildemin);
    Energy currentScale=qtildemin;
    double fact = pow(qtildemax/qtildemin,1./(_npoint-1));
    for(unsigned int ix=1;ix<_npoint;++ix) {
      currentScale *= fact;
      double currentSud = integrator->value(currentScale,scale.back());
      scale.push_back(currentScale);
      sud.push_back(sud.back()+currentSud);
      cerr << "testing values " << scale.back()/GeV << "\t" << sud.back() << " " << exp(-sud.back()) << "\n";
    }
    // convert to the Sudakov
    for(unsigned int ix=0;ix<sud.size();++ix) {
      sud[ix] = exp(-sud[ix]);
    }
    // construct the Interpolators
    Interpolator<double,Energy>::Ptr intq = new_ptr(Interpolator<double,Energy>(sud,scale,3));
    Interpolator<Energy,double>::Ptr ints = new_ptr(Interpolator<Energy,double>(scale,sud,3));
    _fbranchings.insert(make_pair(it->first,make_pair(intq,ints)));




    output << "NEWFRAME\n";
    output << "TITLE TOP \"Sudakov for " << getParticleData(it->second.second[0])->PDGName() << " -> "
	   << getParticleData(it->second.second[1])->PDGName() << " "
	   << getParticleData(it->second.second[2])->PDGName() << "\"\n";
    for(unsigned int ix=0;ix<sud.size();++ix)
      output << scale[ix]/GeV << " " << sud[ix] << "\n";
    output << "JOIN RED\n" << flush;
    HistogramPtr temp(new_ptr(Histogram(0.,100.,200)));

    double slst = (*intq)(91.2*GeV);
    for(unsigned int ix=0;ix<100000000;++ix) {
      double snow = slst/UseRandom::rnd();
      if(snow>=1.) continue;
      Energy qnow = (*ints)(snow);
      *temp +=qnow/GeV;
    }
    using namespace HistogramOptions;
    temp->topdrawOutput(output,Frame);
  }
}

void NasonCKKWHandler::cascade() {
  reweightCKKW( 1, 5);
  ShowerHandler::cascade();
  return;
}

double NasonCKKWHandler::getJetMeasure(ShowerParticlePtr part_i, ShowerParticlePtr part_j){

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

bool NasonCKKWHandler:: getSud( int * qq_pairs, long * emmitter_id,
				  SudakovPtr clusterSudakov,
				  ShowerParticlePtr part_i, 
				  ShowerParticlePtr part_j ) {
  // g 2 q qbar or an incorrect qq type
  if ( abs ( part_i->id() ) < 7 && abs ( part_j->id() ) < 7 ) { 
    if ( abs ( part_i->id() ) != abs ( part_j->id() )  ) return false;
    if ( ( part_i->id() < 0 &&  part_j->id() < 0 ) ||
	 ( part_i->id() > 0 &&  part_j->id() > 0 ) ) return false;
    if ( *qq_pairs < 2 ) return false;
    *qq_pairs -= 1;
    *emmitter_id = 21;
  }
  // q/qbar 2 q/qbar g
  else if ( abs ( part_i->id() ) < 7 || abs ( part_j->id() ) < 7 ) {
    if( abs ( part_i->id() ) < 7 ){
      *emmitter_id = part_i->id();
    }
    else {
      *emmitter_id = part_j->id();
    }
  }
  // g 2 g g
  else {
    *emmitter_id = 21;
  }

  BranchingList branchings = 
    evolver()->splittingGenerator()->finalStateBranchings();

  //cycle through list of all branchings with the correct abs ( emmitter_id )
  for(BranchingList::const_iterator cit = branchings.lower_bound( abs(*emmitter_id) );
      cit != branchings.upper_bound( abs(*emmitter_id) ); ++cit ) {
    IdList ids = cit->second.second;
    if( abs( ids[0] ) == abs( *emmitter_id ) 
	&& ( ( abs(ids[1]) == abs( part_i->id() ) 
	       && abs( ids[2] ) == abs( part_j->id() ) ) ||
	     ( abs( ids[1] ) == abs( part_j->id() ) 
	       && abs( ids[2] ) == abs( part_i->id() ) ) ) ) {
      clusterSudakov = cit->second.first;
      return true; 	    
    }
  }
  return false;
}


bool NasonCKKWHandler::doClustering( ParticleVector theParts, 
				     PPtr vb ) {

  int qq_pairs = 0;
  map <ShowerParticlePtr,NasonBranchingPtr> theParticles;
  tcPDPtr particle_data;
 
  ShowerParticlePtr vBoson = new_ptr( ShowerParticle( *vb, 1, false, false ) );

  //loops through the FS particles
  for( int i = 0; i < theParts.size(); i++){
    ShowerParticlePtr currentParticle = new_ptr( ShowerParticle
						 ( *theParts[i], 1, true, false ) );
    //    currentParticle->rescaleMass();

    if( currentParticle->id() > 0 && currentParticle->id() < 7 ) qq_pairs++;

    theParticles.insert( make_pair( currentParticle, 
				    new_ptr( NasonBranching( currentParticle, SudakovPtr(),
							     NasonBranchingPtr(),false ) ) ) );
    
  }

  double yij;

  //loops clustering until we get down to qqbar
  while( theParticles.size() > 2 ){
    //cerr<<"the number of particles in map = " << theParticles.size()<< " \n";
    double yij_min = 1.;
    
    pair< ShowerParticlePtr, ShowerParticlePtr > clusterPair;

    //loops over all pairs of particles in theParticles
    for( map<ShowerParticlePtr, NasonBranchingPtr>::iterator ita = theParticles.begin();
	 ita != theParticles.end() ; ita++ ) {
      
      for( map<ShowerParticlePtr, NasonBranchingPtr>::iterator itb = theParticles.begin();
	   itb != ita; itb++) {
	yij = getJetMeasure( ita->first, itb->first );
	if( yij < yij_min && splittingAllowed( ita->first, itb->first, qq_pairs )  ) {
	  clusterPair.first = ita->first;
	  clusterPair.second = itb->first;
	  yij_min = yij;
	}
	//	cerr<< "\n \n Particle pair = \n" << *(ita->first)<<"and  \n"<< *(itb->first) <<" \n";
	//	cerr<< "yij for pair = " << yij <<" \n";
	
      }
      
    }
    
    long thePartId;
    SudakovPtr theSudakov; 
    if( !getSud( & qq_pairs, & thePartId, theSudakov,
		 clusterPair.first, clusterPair.second ) ){
      cerr<< "can't find the sudakov in: \n";
      cerr<< *clusterPair.first<<"\n"
	  << *clusterPair.second<<"\n";
      cerr<< "with qq_pairs = " << qq_pairs <<"\n";
    }
    Lorentz5Momentum pairMomentum = clusterPair.first->momentum() + 
      clusterPair.second->momentum();
    pairMomentum.setMass(0.*MeV);
    particle_data = getParticleData( thePartId );
    
    //creates emitter particle
    ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, true ) );
    clustered->set5Momentum( pairMomentum );

    NasonBranchingPtr clusteredBranch(new_ptr(NasonBranching( clustered, theSudakov,
							      NasonBranchingPtr(),false)));
    
    theParticles.insert( make_pair( clustered, clusteredBranch ) );
    clusteredBranch->addChild( theParticles.find( clusterPair.first )->second  );
    clusteredBranch->addChild( theParticles.find( clusterPair.second )->second );
    theParticles.erase( clusterPair.first );
    theParticles.erase( clusterPair.second );
  
  }
  //check the particles are q qbar
  //for(  map<ShowerParticlePtr, NasonBranchingPtr>::iterator it = theParticles.begin(); 
  //	it != theParticles.end(); ++it){ 
  // cerr<< "\n \n Particle = \n" << *(it->first)<<"\n \n";
  // }
  vector<NasonBranchingPtr> theBranchings;
  for(  map<ShowerParticlePtr, NasonBranchingPtr>::iterator it = theParticles.begin(); 
	it != theParticles.end(); ++it){ 
    theBranchings.push_back( it->second );
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
  
  return true;
}
