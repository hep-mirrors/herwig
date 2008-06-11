// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonCKKWHandler class.
//

#include "NasonCKKWHandler.h"
#include "ThePEG/Utilities//CFileLineReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "QTildeSudakovIntegrator.h"
#include "NasonTree.h"

using namespace Herwig;

IBPtr NasonCKKWHandler::clone() const {
  return new_ptr(*this);
}

IBPtr NasonCKKWHandler::fullclone() const {
  return new_ptr(*this);
}

void NasonCKKWHandler::persistentOutput(PersistentOStream & os) const {
  os  << _alphaS << _sudopt << _sudname << _jetMeasureMode;
}

void NasonCKKWHandler::persistentInput(PersistentIStream & is, int) {
  is  >> _alphaS >> _sudopt >> _sudname >> _jetMeasureMode;
}

ClassDescription<NasonCKKWHandler> NasonCKKWHandler::initNasonCKKWHandler;
// Definition of the static class description member.

void NasonCKKWHandler::Init() {

  static ClassDocumentation<NasonCKKWHandler> documentation
    ("The NasonCKKWHandler class manages the implementation of the CKKW approach using"
     "the truncated shower.");

  static Reference<NasonCKKWHandler,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &NasonCKKWHandler::_alphaS, false, false, true, false, false);

  static Switch<NasonCKKWHandler,unsigned int> interfaceSudakovOption
    ("SudakovOption",
     "Option for the initialisation of the Sudakov tables",
     &NasonCKKWHandler::_sudopt, 0, false, false);
  static SwitchOption interfaceSudakovOptionWrite
    (interfaceSudakovOption,
     "Write",
     "Calculate the Sudakov and write the table to a file",
     1);
  static SwitchOption interfaceSudakovOptionRead
    (interfaceSudakovOption,
     "Read",
     "Read the Sudakov table from a file",
     2);
  static SwitchOption interfaceSudakovOptionCompute
    (interfaceSudakovOption,
     "Compute",
     "Calculate the Sudakov but don't write the table",
     0);

  static Parameter<NasonCKKWHandler,string> interfaceSudakovName
    ("SudakovName",
     "Name for the file containing the Sudakov form factors",
     &NasonCKKWHandler::_sudname, "sudakov.data",
     false, false);

  static Switch<NasonCKKWHandler, unsigned int> ifaceJetMeasureMode
    ("JetMeasure",
     "Choice of the jet measure algorithm",
     &NasonCKKWHandler::_jetMeasureMode, 0, false, false);
  
  static SwitchOption Durham
    (ifaceJetMeasureMode,"Durham","Durham jet algorithm", 0);
  
  static SwitchOption LUCLUS
    (ifaceJetMeasureMode,"LUCLUS","LUCLUS jet algorithm", 1);
}

double NasonCKKWHandler::reweightCKKW(int minMult, int maxMult) {

  PPair in = lastXCombPtr()->subProcess()->incoming();
  ParticleVector out  = lastXCombPtr()->subProcess()->outgoing();
  PPtr vBoson = lastXCombPtr()->subProcess()->intermediates()[0];
  _s = lastXCombPtr()->lastS();

  NasonTreePtr _theNasonTree = doClustering( out, vBoson );
  if(!_theNasonTree) return 0.;

  // cerr<<"finished do clustering \n";

  //the jet resolution parameter used in MG generation
  //Q1 is the lower scale and should be used as a cut off 
  //  double y_ini = 0.001;
  Energy Q = sqrt( _s );
  //  Energy Q1 = sqrt( y_ini ) * Q;

  //the Sudakov weight factor
  double SudWgt = 1.;
  
  //include the sud factor for each external line
  
  for( map<ShowerParticlePtr,double>::const_iterator cit 
	 = _theExternals.begin();
       cit != _theExternals.end(); ++cit ) {
    //itererate over all matching sudakovs (only 1 except for gluon)
    multimap< long, pair < Interpolator<double,Energy>::Ptr,
      Interpolator<Energy,double>::Ptr >  >::const_iterator cjt;
    for( cjt =  _fbranchings.lower_bound( abs( cit->first->id() ) );
	 cjt != _fbranchings.upper_bound( abs( cit->first->id() ) );
	 ++cjt ) {
      SudWgt *= (* cjt->second.first )( sqrt( cit->second ) * Q );
    }
  }


  //include the intermediate line wgts
  for( map< long, pair< double, double > >::const_iterator cit 
	 = _theIntermediates.begin();
       cit != _theIntermediates.end(); ++cit ) {

    //itererate over all matching sudakovs (only 1 except for gluon)
    multimap< long, pair < Interpolator< double, Energy >::Ptr, Interpolator<Energy,double>::Ptr > >::const_iterator cjt;
    for( cjt =  _fbranchings.lower_bound( abs( cit->first ) );
	 cjt != _fbranchings.upper_bound( abs( cit->first ) );
	 ++cjt ) {
      SudWgt *= (* cjt->second.first )( sqrt( cit->second.first  ) * Q );
      SudWgt /= (* cjt->second.first )( sqrt( cit->second.second ) * Q );
    }
  }
 
  double alphaWgt = 1.;
  //need to add the alphaS weight
  //the alphaS ratio evaluated at all nodal values
  for( map<NasonBranchingPtr,double>::const_iterator cit 
	 = _theNodes.begin();
       cit != _theNodes.end(); ++cit ) {
    alphaWgt *= _alphaS->ratio( cit->second * sqr( Q ) );
  }

  //update the sub process 
  ParticleVector outgoing = lastXCombPtr()->subProcess()->outgoing();
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    lastXCombPtr()->subProcess()->removeEntry(outgoing[ix]);
    tParticleVector parents=outgoing[ix]->parents();
    for(unsigned int iy=0;iy<parents.size();++iy)
      parents[iy]->abandonChild(outgoing[ix]);
  }
  // add new ones based on the NasonTree
  map<ColinePtr,ColinePtr> colourMap;
  for(set<NasonBranchingPtr>::const_iterator it=_theNasonTree->branchings().begin();
      it!=_theNasonTree->branchings().end();++it) {
    if((**it)._incoming) continue;
    generator()->log() << *(**it)._particle << "\n";
    PPtr newParticle = new_ptr(Particle((**it)._particle->dataPtr()));
    newParticle->set5Momentum((**it)._shower);
    if((**it)._particle->colourLine()) {
      map<ColinePtr,ColinePtr>::iterator loc 
	= colourMap.find((**it)._particle->colourLine());
      if(loc!=colourMap.end()) loc->second->addColoured(newParticle);
      else {
	ColinePtr newLine=new_ptr(ColourLine());
	colourMap[(**it)._particle->colourLine()]=newLine;
	newLine->addColoured(newParticle);
      }
    }
    if((**it)._particle->antiColourLine()) {
      map<ColinePtr,ColinePtr>::iterator loc 
	= colourMap.find((**it)._particle->antiColourLine());
      if(loc!=colourMap.end()) loc->second->addAntiColoured(newParticle);
      else {
	ColinePtr newLine=new_ptr(ColourLine());
	colourMap[(**it)._particle->antiColourLine()]=newLine;
	newLine->addAntiColoured(newParticle);
      }
    }
    lastXCombPtr()->subProcess()->addOutgoing(newParticle);
  }
  return SudWgt * alphaWgt;
}

void NasonCKKWHandler::doinitrun() {
  ShowerHandler::doinitrun();
  // integrator for the outer integral
  GaussianIntegrator outer;
  // get the final-state branchings from the evolver
  if(_sudopt!=2) {
    ofstream sudFileOutput;
    if(_sudopt==1) sudFileOutput.open(_sudname.c_str());
    for(BranchingList::const_iterator 
	  it = evolver()->splittingGenerator()->finalStateBranchings().begin();
	it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
      // class to return the integrated factor in the exponent
      Ptr<QTildeSudakovIntegrator>::pointer integrator = 
	new_ptr(QTildeSudakovIntegrator(it->second));
      if(_sudopt==1) sudFileOutput << it->second.first->fullName() << "\t"
				   << it->second.second[0] << "\t"
				   << it->second.second[1] << "\t"
				   << it->second.second[2] << "\n";
      Energy qtildemax = generator()->maximumCMEnergy();
      Energy qtildemin = integrator->minimumScale();
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
      }
      // convert to the Sudakov
      for(unsigned int ix=0;ix<sud.size();++ix) sud[ix] = exp(-sud[ix]);
      // construct the Interpolators
      Interpolator<double,Energy>::Ptr 
	intq = new_ptr(Interpolator<double,Energy>(sud,scale,3));
      Interpolator<Energy,double>::Ptr 
	ints = new_ptr(Interpolator<Energy,double>(scale,sud,3));
      _fbranchings.insert( make_pair( it->first, make_pair( intq, ints ) ) );
      if(_sudopt==1) {
	sudFileOutput << scale.size() << "\n";
	for(unsigned int ix=0;ix<scale.size();++ix)
	  sudFileOutput << setprecision(18) << scale[ix]/GeV << "\t" << sud[ix] << "\n";
      }
    }
    sudFileOutput.close();
  }
  else {
    CFileLineReader file(_sudname);
    while(file.readline()) {
      string line = file.getline(), name;
      istringstream is;
      is.str(line);
      IdList ids(3);
      is >> name >> ids[0] >> ids[1] >> ids[2];
      file.readline();
      vector<Energy> scale;
      vector<double> sud;
      unsigned int isize;
      double temp[2];
      is.str(file.getline());
      is >> isize;
      for(unsigned int ix=0;ix<isize;++ix) {
	file.readline();
	is.str(file.getline());
	is >> temp[0] >> temp[1];
	scale.push_back(temp[0]*GeV);
	sud.push_back(temp[1]);
      }
      BranchingList::const_iterator it,
	start = evolver()->splittingGenerator()->finalStateBranchings().lower_bound(ids[0]),
	end   = evolver()->splittingGenerator()->finalStateBranchings().upper_bound(ids[0]);
      for(it=start;it!=end;++it) {
	if(it->second.first->fullName()==name&&
	   it->second.second[0]==ids[0]&&
	   it->second.second[1]==ids[1]&&
	   it->second.second[2]==ids[2]) {
	  // construct the Interpolators
	  Interpolator<double,Energy>::Ptr 
	    intq = new_ptr(Interpolator<double,Energy>(sud,scale,3));
	  Interpolator<Energy,double>::Ptr 
	    ints = new_ptr(Interpolator<Energy,double>(scale,sud,3));
	  _fbranchings.insert( make_pair( it->first, make_pair( intq, ints ) ) );
	  break;
	}
      }
      if(it==end) {
	cerr << "testing fails\n";
      }
    }
  }
//   ofstream output("sudakovs.top");
//   for(BranchingList::const_iterator 
// 	it = evolver()->splittingGenerator()->finalStateBranchings().begin();
//       it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
//      output << "NEWFRAME\n";
//      output << "TITLE TOP \"Sudakov for " << getParticleData(it->second.second[0])->PDGName() << " -> "
//  	   << getParticleData(it->second.second[1])->PDGName() << " "
//  	   << getParticleData(it->second.second[2])->PDGName() << "\"\n";
//      for(unsigned int ix=0;ix<sud.size();++ix)
//        output << scale[ix]/GeV << " " << sud[ix] << "\n";
//      output << "JOIN RED\n" << flush;
 //     HistogramPtr temp(new_ptr(Histogram(0.,100.,200)));

 //     double slst = (*intq)(91.2*GeV);
 //     for(unsigned int ix=0;ix<100000000;++ix) {
 //       double snow = slst/UseRandom::rnd();
 //       if(snow>=1.) continue;
 //       Energy qnow = (*ints)(snow);
 //       *temp +=qnow/GeV;
 //     }
 //     using namespace HistogramOptions;
 //     temp->topdrawOutput(output,Frame);
}

void NasonCKKWHandler::doinit() throw(InitException) {
  
  ShowerHandler::doinit();
 
}

void NasonCKKWHandler::cascade() {
  ShowerHandler::cascade();
}

double NasonCKKWHandler::getJetMeasure(ShowerParticlePtr part_i,
				       ShowerParticlePtr part_j){
  double yij;
  double costheta = part_i->momentum().vect().dot( part_j->momentum().vect() ) 
    / part_i->momentum().vect().mag() / part_j->momentum().vect().mag();
  switch( _jetMeasureMode ){
  case 0:
    if( sqr( part_i->momentum().e() ) > sqr( part_j->momentum().e() ) )
      yij = 2. * sqr( part_j->momentum().e() ) * ( 1. - costheta ) / _s ;
    else
      yij = 2. * sqr( part_i->momentum().e() ) * ( 1. - costheta ) / _s ;
    break;
  case 1:
    yij = 2. * sqr( part_i->momentum().e() * part_j->momentum().e() /
		    ( part_i->momentum().e() + part_j->momentum().e() ) )/_s
      * ( 1. - costheta );
    break;
  default:
    yij = 1.;
    break;
  }
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
  _theNodes.clear();
  _theExternals.clear();
  _theIntermediates.clear();

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
    //insert all particles into externals and initialise all
    //jet res parameters to 1
    cerr<< "external: \n"
	<< currentParticle << "\n";
    _theExternals.insert( make_pair( currentParticle, 1. ) );
    
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
    SudakovPtr theSudakov = getSud( qq_pairs, thePartId,
				    clusterPair.first, clusterPair.second ); 
    if( !theSudakov ){
      cerr << "can't find the sudakov in: \n";
      cerr << *clusterPair.first<<"\n"
	   << *clusterPair.second<<"\n";
      cerr << "with qq_pairs = " << qq_pairs <<"\n";
    }

    cerr << "clustering with yij = "<< yij_min <<":  "<< *clusterPair.first<<"\n"
	 << *clusterPair.second<<"\n";

    generator()->log() << "testing sudakov?? " << theSudakov << "\n";
    Lorentz5Momentum pairMomentum = clusterPair.first->momentum() + 
      clusterPair.second->momentum();
    pairMomentum.setMass(0.*MeV);
    particle_data = getParticleData( thePartId );
    
    //creates emitter particle
    ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, true ) );
    clustered->set5Momentum( pairMomentum );   
 
    //decide which particle was created (ie which was only resolvable 
    //at this scale) from ids (if q->qg) or magnitudes
    int created;
    //are we clustering quarks
    if( abs( clusterPair.first->id() )  < 7 || 
	abs( clusterPair.second->id() ) < 7 ){
      if( abs( thePartId ) ==
	  abs( clusterPair.first->id() ) )  created = 1;
      else if( abs( thePartId ) ==
	       abs( clusterPair.second->id() ) )  created = 2;
      //clustering q qbar
      else{
	if( clusterPair.first->momentum().mag() > 
	    clusterPair.second->momentum().mag() ) created = 1;
	else created = 2;
      }
    }
    else{
      	if( clusterPair.first->momentum().mag() > 
	    clusterPair.second->momentum().mag() ) created = 1;
	else created = 2;
    }
  
    NasonBranchingPtr clusteredBranch( new_ptr( NasonBranching( clustered, theSudakov,
								NasonBranchingPtr(), false ) ) );
    fixColours( clustered, clusterPair.first, clusterPair.second );
    theParticles.insert( make_pair( clustered, clusteredBranch ) );
    
    //add children in the correct order
    if( created == 2 ){
      clusteredBranch->addChild( theParticles.find( clusterPair.first )
				 ->second  );
      clusteredBranch->addChild( theParticles.find( clusterPair.second )
				 ->second );
    }
    else{
      clusteredBranch->addChild( theParticles.find( clusterPair.second )
				 ->second  );
      clusteredBranch->addChild( theParticles.find( clusterPair.first )
				 ->second );
    }

    _theNodes.insert( make_pair( clusteredBranch, yij_min ) );

    theParticles.erase( clusterPair.first  );
    theParticles.erase( clusterPair.second );
  
  }
  vector<NasonBranchingPtr> theBranchings;
  for(  map<ShowerParticlePtr, NasonBranchingPtr>::iterator it = 
	  theParticles.begin(); 
	it != theParticles.end(); ++it ) 
    theBranchings.push_back( it->second );
  
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
  for( unsigned int jx = 0; jx < theBranchings.size(); ++jx ) {
    generator()->log() << "testing " << theBranchings[jx]->_shower/GeV << "\t"
		       << theBranchings[jx]->_shower.m()/GeV << "\n";
  }

 

  //set the scales at which the externals are resolved
  for( map<NasonBranchingPtr,double>::const_iterator cit 
	 = _theNodes.begin();
       cit != _theNodes.end(); ++cit ) {
    NasonBranchingPtr currentExternal = cit->first;
    //loop until the  assosiated external particle is found
    do {
      if( ! currentExternal ) break;
      currentExternal = currentExternal->_children[0];
    } while( _theExternals.find( currentExternal->branchingParticle() ) ==
	     _theExternals.end() );

    if( _theExternals.find( currentExternal->branchingParticle() ) !=
	_theExternals.end() )
      _theExternals.find( currentExternal->branchingParticle() )->
	second = cit->second;
    else cerr<<"doClustering()::can't find an external \n"
	     << currentExternal->branchingParticle()
	     << "\n";    
  }
  //get the intermediates
  for( map<NasonBranchingPtr,double>::const_iterator cit 
	 = _theNodes.begin();
       cit != _theNodes.end(); ++cit ) {
    //end scale and intermediate are given by theNodes
    double endScale = cit->second;
    double startScale = 1.;
    long intID = cit->first->branchingParticle()->id();
    //get start scale from parton
    if( cit->first->_parent ) {
      NasonBranchingPtr intParent = cit->first;
      if( _theNodes.find( cit->first ) !=
	  _theNodes.end() )
	startScale = _theNodes.find( cit->first )->second;
      else cerr<<"doClustering()::can't find start of intermediate with parents \n";
    }
    _theIntermediates.insert( make_pair( intID, make_pair( startScale, endScale ) ) );
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
