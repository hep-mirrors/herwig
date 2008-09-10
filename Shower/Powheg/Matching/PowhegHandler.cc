// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegHandler class.
//

#include "PowhegHandler.h"
#include "ThePEG/Utilities//CFileLineReader.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "QTildeSudakovIntegrator.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/PDF/HwRemDecayer.h"
#include <queue>

using namespace Herwig;

IBPtr PowhegHandler::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegHandler::fullclone() const {
  return new_ptr(*this);
}

void PowhegHandler::persistentOutput(PersistentOStream & os) const {
  os  << _alphaS << _sudopt << _sudname << _jetMeasureMode << _allowedInitial
      << _allowedFinal << _matrixElement << _lepton <<  _yini << _alphaSMG
      << _npoint <<  ounit( _max_qtilde, GeV );
}

void PowhegHandler::persistentInput(PersistentIStream & is, int) {
  is  >> _alphaS >> _sudopt >> _sudname >> _jetMeasureMode >> _allowedInitial
      >> _allowedFinal >> _matrixElement >> _lepton >> _yini >> _alphaSMG
      >> _npoint >> iunit( _max_qtilde, GeV );
}

ClassDescription<PowhegHandler> PowhegHandler::initPowhegHandler;
// Definition of the static class description member.

void PowhegHandler::Init() {

  static ClassDocumentation<PowhegHandler> documentation
    ("The PowhegHandler class manages the implementation of the CKKW approach using"
     "the truncated shower.");

  static Reference<PowhegHandler,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &PowhegHandler::_alphaS, false, false, true, false, false);

  static Switch<PowhegHandler,unsigned int> interfaceSudakovOption
    ("SudakovOption",
     "Option for the initialisation of the Sudakov tables",
     &PowhegHandler::_sudopt, 0, false, false);
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

  static Parameter<PowhegHandler,string> interfaceSudakovName
    ("SudakovName",
     "Name for the file containing the Sudakov form factors",
     &PowhegHandler::_sudname, "sudakov.data",
     false, false);

  static Switch<PowhegHandler, unsigned int> ifaceJetMeasureMode
    ("JetMeasure",
     "Choice of the jet measure algorithm",
     &PowhegHandler::_jetMeasureMode, 0, false, false);
  
  static SwitchOption Durham
    (ifaceJetMeasureMode,"Durham","Durham jet algorithm", 0);
  
  static SwitchOption LUCLUS
    (ifaceJetMeasureMode,"LUCLUS","LUCLUS jet algorithm", 1);

  static Parameter<PowhegHandler,double> interfaceMergeScale
    ("MergeScale",
     "The CKKW merging scale, yini",
     &PowhegHandler::_yini, 0.001, 0.0, 1.0,
     false, false, Interface::limited );

   static Parameter<PowhegHandler,double> interfaceAlphaSMG
    ("alphaSMG",
     "The fixed alphas used in MG event generation",
     &PowhegHandler::_alphaSMG, 0.118, 0.0, 1.0,
     false, false, Interface::limited );

  static Switch<PowhegHandler,bool> interfaceLepton
    ("Lepton",
     "Whether is a hadron-hadron or lepton-lepton collision",
     &PowhegHandler::_lepton, true, false, false);
  static SwitchOption interfaceLeptonLeptonic
    (interfaceLepton,
     "Leptonic",
     "Leptonic collision",
     true);
  static SwitchOption interfaceLeptonHadronic
    (interfaceLepton,
     "Hadronic",
     "Hadronic collision",
     false);

  static Reference<PowhegHandler,MEBase> interfaceMatrixElement
    ("MatrixElement",
     "The matrix element class for the core 2->2 process",
     &PowhegHandler::_matrixElement, false, false, true, true, false);

  static Parameter<PowhegHandler,unsigned int> interfaceInterpPoints
    ("InterpolatorPoints",
     "The number of points used for sudakov interpolation tables",
     &PowhegHandler::_npoint, 100, 0, 1000000,
     false, false, Interface::limited );

  static Parameter<PowhegHandler, Energy> interfaceMaxQTilde
    ("maxQTilde",
     "The maximum QTilde scale for sudakov interpolation tables",
     &PowhegHandler::_max_qtilde, GeV, 91.2*GeV, 1.*GeV, 1000000.*GeV,
     false, false, Interface::limited);
}

double PowhegHandler::reweightCKKW(int minMult, int maxMult) {
  // cluster the event
  _theHardTree = doClustering();
  // return if fails
  if(!_theHardTree) {
    return 0.;
  }
  // compute the Sudakov weight
  double SudWgt = _lepton ? sudakovWeight() : 1.;
  //update the sub process
  if(_lepton) {
    ParticleVector outgoing = lastXCombPtr()->subProcess()->outgoing();
    for(unsigned int ix=0;ix<outgoing.size();++ix) {
      lastXCombPtr()->subProcess()->removeEntry(outgoing[ix]);
      tParticleVector parents=outgoing[ix]->parents();
      for(unsigned int iy=0;iy<parents.size();++iy)
	parents[iy]->abandonChild(outgoing[ix]);
    }
    // add new ones based on the HardTree
    map<ColinePtr,ColinePtr> colourMap;
    for(set<HardBranchingPtr>::const_iterator it=_theHardTree->branchings().begin();
	it!=_theHardTree->branchings().end();++it) {
      if((**it).incoming()) continue;
      PPtr newParticle = new_ptr(Particle((**it).branchingParticle()->dataPtr()));
      newParticle->set5Momentum((**it).showerMomentum());
      //do colour connections
      if((**it).branchingParticle()->colourLine()) {
	map<ColinePtr,ColinePtr>::iterator loc 
	  = colourMap.find((**it).branchingParticle()->colourLine());
	if(loc!=colourMap.end()) loc->second->addColoured(newParticle);
	else {
	  ColinePtr newLine=new_ptr(ColourLine());
	  colourMap[(**it).branchingParticle()->colourLine()]=newLine;
	  newLine->addColoured(newParticle);
	}
      }
      if((**it).branchingParticle()->antiColourLine()) {
	map<ColinePtr,ColinePtr>::iterator loc 
	  = colourMap.find((**it).branchingParticle()->antiColourLine());
	if(loc!=colourMap.end()) loc->second->addAntiColoured(newParticle);
	else {
	  ColinePtr newLine=new_ptr(ColourLine());
	  colourMap[(**it).branchingParticle()->antiColourLine()]=newLine;
	  newLine->addAntiColoured(newParticle);
	}
      }
      lastXCombPtr()->subProcess()->addOutgoing(newParticle);
    }
  }
  else {
    set<HardBranchingPtr>::const_iterator it; 
    map<ColinePtr,ColinePtr> colourMap;
    ParticleVector outgoing;
    PPair incoming;
    for(it=_theHardTree->branchings().begin();
	it!=_theHardTree->branchings().end();++it) {
      PPtr newParticle = new_ptr(Particle((**it).branchingParticle()->dataPtr()));
      newParticle->set5Momentum((**it).showerMomentum());
      if((**it).branchingParticle()->colourLine()) {
	map<ColinePtr,ColinePtr>::iterator loc 
	  = colourMap.find((**it).branchingParticle()->colourLine());
	if(loc!=colourMap.end()) loc->second->addColoured(newParticle);
	else {
	  ColinePtr newLine=new_ptr(ColourLine());
	  colourMap[(**it).branchingParticle()->colourLine()]=newLine;
	  newLine->addColoured(newParticle);
	}
      }
      if((**it).branchingParticle()->antiColourLine()) {
	map<ColinePtr,ColinePtr>::iterator loc 
	  = colourMap.find((**it).branchingParticle()->antiColourLine());
	if(loc!=colourMap.end()) loc->second->addAntiColoured(newParticle);
	else {
	  ColinePtr newLine=new_ptr(ColourLine());
	  colourMap[(**it).branchingParticle()->antiColourLine()]=newLine;
	  newLine->addAntiColoured(newParticle);
	}
      }
      if((**it).incoming()) {
	if(lastXCombPtr()->subProcess()->incoming().first->momentum().z()/
	   newParticle->momentum().z()>0.)
	  incoming.first = newParticle;
	else
	  incoming.second = newParticle;
      }
      else
	outgoing.push_back(newParticle);
    }
    SubProPtr newSubProcess=
      new_ptr(SubProcess(incoming,
			 lastXCombPtr()->subProcess()->collision(),
			 lastXCombPtr()->subProcess()->handler()));
    for(unsigned int ix=0;ix<outgoing.size();++ix)
      newSubProcess->addOutgoing(outgoing[ix]);
    lastXCombPtr()->subProcess(newSubProcess);
  }
  return SudWgt;
}

void PowhegHandler::dofinish() {
  ShowerHandler::dofinish();
  string fname = generator()->filename() + string("-") 
    + string("wgts.top");
  ofstream output(fname.c_str());

  using namespace HistogramOptions;

  _hSud->topdrawOutput(output,Frame,
		       "RED",
		       "Sudakov wgts",
		       "",
		       "freq",
		       "",
		       "wgt",
		       "");
 
     
  _halphaS->topdrawOutput(output,Frame,
		       "RED",
		       "AlphaS wgts",
		       "",
		       "freq",
		       "",
		       "wgt",
		       "");
}


void PowhegHandler::doinitrun() {
  ShowerHandler::doinitrun();  
  _s = sqr( generator()->maximumCMEnergy() );

  _hSud = new_ptr(Histogram(0.,2.,100));
  _halphaS = new_ptr(Histogram(0.,2.,100));
 
  // integrator for the outer integral
  GaussianIntegrator outer;
  // get the final-state branchings from the evolver
  if(_sudopt!=2) {
    ofstream sudFileOutput;
    if(_sudopt==1) sudFileOutput.open(_sudname.c_str());
    for(BranchingList::const_iterator 
	  it = evolver()->splittingGenerator()->finalStateBranchings().begin();
	it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
   
      Ptr<QTildeSudakovIntegrator>::pointer integrator = 
	new_ptr( QTildeSudakovIntegrator(it->second, sqrt( _yini * _s ), _jetMeasureMode ) );
      if(_sudopt==1) sudFileOutput << it->second.first->fullName() << "\t"
				   << it->second.second[0] << "\t"
				   << it->second.second[1] << "\t"
				   << it->second.second[2] << "\n";
      Energy qtildemax = _max_qtilde;
      Energy qtildemin = integrator->minimumScale();
    
      vector<double> sud;
    
      vector<Energy> scale;
      sud.push_back(0.); scale.push_back(qtildemin);

      Energy currentScale = qtildemin;
      double fact = pow(qtildemax/qtildemin,1./double(_npoint-1));
      for(unsigned int ix=1;ix<_npoint;++ix) {
	currentScale *= fact;
	double currentSud = integrator->value(currentScale,scale.back());
	scale.push_back(currentScale);
	sud.push_back(sud.back()+currentSud);
      }
      // convert to the Sudakov
      for(unsigned int ix=0;ix<sud.size();++ix) {
	sud[ix] = exp(-sud[ix]);
      }
   
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

}

void PowhegHandler::doinit() throw(InitException) {
  ShowerHandler::doinit();
  // extract the allowed branchings
  // final-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->finalStateBranchings().begin();
      it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
    pair<long,long> prod(make_pair(it->second.second[1],it->second.second[2]));
    _allowedFinal.insert(make_pair(prod,it->second));
    swap(prod.first,prod.second);
    _allowedFinal.insert(make_pair(prod,it->second));
  }
  // initial-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->initialStateBranchings().begin();
      it != evolver()->splittingGenerator()->initialStateBranchings().end(); ++it) {
    _allowedInitial.insert(make_pair(it->second.second[0],it->second));
  }

}

void PowhegHandler::cascade() {
  ShowerHandler::cascade();
}

double PowhegHandler::getJetMeasure(ShowerParticlePtr part_i,
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
bool PowhegHandler:: splittingAllowed( ShowerParticlePtr part_i,
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

SudakovPtr PowhegHandler:: getSud( int & qq_pairs, long & emmitter_id,
				      ShowerParticlePtr & part_i, 
				      ShowerParticlePtr & part_j ) {
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
    if( abs( ids[0] ) == abs( emmitter_id ) ) {
      if( abs(ids[1]) == abs(part_i->id()) && 
	  abs(ids[2]) == abs(part_j->id()) ) {
	return cit->second.first;
      }
      if( abs(ids[1]) == abs(part_j->id()) && 
	  abs(ids[2]) == abs(part_i->id())  ) {
	swap(part_i,part_j);
	return cit->second.first;
      }
    }
  }
  return SudakovPtr();
}


HardTreePtr PowhegHandler::doClustering() {
  if(!_lepton) {
    return generalClustering();
  }


  ParticleVector theParts  = lastXCombPtr()->subProcess()->outgoing();

  //make an intermediate and add to subprocess if not read in
  if(lastXCombPtr()->subProcess()->intermediates().empty()) {
    return HardTreePtr();
    PPair theIncomings =  lastXCombPtr()->subProcess()->incoming();
    //set intermediate to Z
    long intermediate_id = 23;
    PPtr theIntermediate = new_ptr( Particle( getParticleData( intermediate_id ) ) );
    theIntermediate->set5Momentum( theIncomings.first->momentum() +
				   theIncomings.second->momentum() );
    //add the intermediate - parent/child relations should be updated
    lastXCombPtr()->subProcess()->addIntermediate( theIntermediate );
    cerr<<"added intermediate\n"
	<< *theIntermediate<<"\n";

  }
 
    
  PPtr vb = lastXCombPtr()->subProcess()->intermediates()[0];
  _theNodes.clear();
  _theExternals.clear();
  _theIntermediates.clear();

  int qq_pairs = 0;
  map <ShowerParticlePtr,HardBranchingPtr> theParticles;
  tcPDPtr particle_data;
  ShowerParticlePtr vBoson = new_ptr( ShowerParticle( *vb, 1, false, false ) );
  //loops through the FS particles and create naon branchings
  for( unsigned int i = 0; i < theParts.size(); i++){
    ShowerParticlePtr currentParticle = 
      new_ptr( ShowerParticle( *theParts[i], 1, true, false ) );
    //    currentParticle->rescaleMass();
    if( currentParticle->id() > 0 && currentParticle->id() < 7 ) qq_pairs++;
    theParticles.insert(make_pair(currentParticle, 
				  new_ptr( HardBranching( currentParticle, SudakovPtr(),
							   HardBranchingPtr(),false ) ) ) );
    //insert all particles into externals and initialise all
    //jet res parameters to 1
    _theExternals.insert( make_pair( currentParticle, 
				     make_pair( 1., HardBranchingPtr() ) ) );
    
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
    for( map<ShowerParticlePtr, HardBranchingPtr>::iterator ita = theParticles.begin();
	 ita != theParticles.end() ; ita++ ) {
      for( map<ShowerParticlePtr, HardBranchingPtr>::iterator itb = theParticles.begin();
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
    Lorentz5Momentum pairMomentum = clusterPair.first->momentum() + 
      clusterPair.second->momentum();
    pairMomentum.setMass(0.*MeV);
    particle_data = getParticleData( thePartId );
    
    //creates emitter particle
    ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, true ) );
    clustered->set5Momentum( pairMomentum );
  
    HardBranchingPtr clusteredBranch( new_ptr( HardBranching( clustered, theSudakov,
								HardBranchingPtr(), false ) ) );
    fixColours( clustered, clusterPair.first, clusterPair.second );
    theParticles.insert( make_pair( clustered, clusteredBranch ) );
    
    //add children
    clusteredBranch->addChild(theParticles.find(clusterPair.first )->second);
    clusteredBranch->addChild(theParticles.find(clusterPair.second)->second);
    
    _theNodes.insert( make_pair( clusteredBranch, yij_min ) );
    
    theParticles.erase( clusterPair.first  );
    theParticles.erase( clusterPair.second );
    /*
      cerr<<"\n\n clustered particles: \n"
	  << * clusterPair.first <<" \n and\n"
	  << * clusterPair.second << "\n into \n"
	  << * clustered <<"\n "
	  <<"scale of clustering was "<<yij_min<<"\n";
    */
    //search the externals for each of the 2 particles being clustered
    //if found update the yij scale and nason branching of each external
    if( _theExternals.find( clusterPair.first ) != _theExternals.end() ){
      _theExternals.find( clusterPair.first )->second.first = yij_min;
      _theExternals.find( clusterPair.first )->second.second = clusteredBranch;
    }
    if( _theExternals.find( clusterPair.second ) != _theExternals.end() ){
      _theExternals.find( clusterPair.second )->second.first = yij_min;
      _theExternals.find( clusterPair.second )->second.second = clusteredBranch;
    }
  }
  vector<HardBranchingPtr> theBranchings;
  for(  map<ShowerParticlePtr, HardBranchingPtr>::iterator it = 
	  theParticles.begin(); 
	it != theParticles.end(); ++it ) 
    theBranchings.push_back( it->second );
  // fix for e+e- to match up the colours of the q qbar pair
  if(theBranchings[0]->branchingParticle()->dataPtr()->iColour()==PDT::Colour3) {
    ColinePtr temp = theBranchings[1]->branchingParticle()->antiColourLine();
    temp->addColoured(theBranchings[0]->branchingParticle());
    theBranchings[0]->branchingParticle()->colourLine()->join(temp);
  }
  else {
    ColinePtr temp = theBranchings[0]->branchingParticle()->antiColourLine();
    temp->addColoured(theBranchings[1]->branchingParticle());
    theBranchings[1]->branchingParticle()->colourLine()->join(temp);
  }
  theBranchings[0]->colourPartner(theBranchings[1]);
  theBranchings[1]->colourPartner(theBranchings[0]);
  vector<HardBranchingPtr> spaceBranchings;
  spaceBranchings.push_back( new_ptr( HardBranching( vBoson, SudakovPtr(),
						      HardBranchingPtr(), 
						      true ) ) );
  theBranchings.push_back( spaceBranchings.back() );
  HardTreePtr powhegtree = new_ptr( HardTree( theBranchings,
					       spaceBranchings ) );

  // Calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructDecayJets(powhegtree,evolver());

  generator()->log() << "testing hard momenta for the shower\n";
  for( unsigned int jx = 0; jx < theBranchings.size(); ++jx ) {
    generator()->log() << "testing " << theBranchings[jx]->showerMomentum()/GeV << "\t"
		       << theBranchings[jx]->showerMomentum().m()/GeV << "\n";
  }
  //get the intermediates - there is one intermediate for each node
  for( map<HardBranchingPtr,double>::const_iterator cit 
	 = _theNodes.begin();
       cit != _theNodes.end(); ++cit ) {

    //end scale and intermediate are given by theNodes
    double endY = cit->second;
    double startY = 1.;
    Energy endScale = cit->first->scale();
    Energy startScale = sqrt( _s );

    long intID = cit->first->branchingParticle()->id();
    //get start scale from parent of the hardBranching
    //is there a parent hardbranching
    if( cit->first->parent() ) {
      HardBranchingPtr intParent = cit->first->parent();
      //find the parent
      if( _theNodes.find( intParent ) != _theNodes.end() ){
	startScale = intParent->scale() * cit->first->z();
	startY = _theNodes.find( cit->first )->second;
      }
      else cerr<<"doClustering()::can't find the parent of the intermediate in nason branchings\n";
    }
    //the parent was a q or qbar from hard sub process
    //set scale to that of the showerparticle
    else {
      startScale = cit->first->branchingParticle()->evolutionScale();
    }
    if ( startScale < endScale ) endScale = startScale;
    _theIntermediates.insert( make_pair( intID, 
					 make_pair( make_pair( startY, startScale ), 
						    make_pair( endY, endScale  ) ) ) );
  }
  return powhegtree;
}

void PowhegHandler::fixColours(tPPtr parent, tPPtr child1, tPPtr child2) {
  // the different possible cases
  if(parent->dataPtr()->iColour()==PDT::Colour3&&
     child1->dataPtr()->iColour()==PDT::Colour3&&
     child2->dataPtr()->iColour()==PDT::Colour8) {
    child2->colourLine()->addColoured(parent);
    ColinePtr temp = child2->antiColourLine();
    temp->addColoured(child1);
    child1->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3&&
	  child2->dataPtr()->iColour()==PDT::Colour3&&
	  child1->dataPtr()->iColour()==PDT::Colour8) {
    child1->colourLine()->addColoured(parent);
    ColinePtr temp = child1->antiColourLine();
    temp->addColoured(child2);
    child2->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3bar&&
	  child1->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour8) {
    child2->antiColourLine()->addAntiColoured(parent);
    ColinePtr temp = child1->antiColourLine();
    temp->addColoured(child2);
    child2->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour3bar&&
	  child2->dataPtr()->iColour()==PDT::Colour3bar&&
	  child1->dataPtr()->iColour()==PDT::Colour8) {
    child1->antiColourLine()->addAntiColoured(parent);
    ColinePtr temp = child2->antiColourLine();
    temp->addColoured(child1);
    child1->colourLine()->join(temp);
  }
  else if(parent->dataPtr()->iColour()==PDT::Colour8&&
	  child1->dataPtr()->iColour()==PDT::Colour8&&
	  child2->dataPtr()->iColour()==PDT::Colour8) {
    if(UseRandom::rndbool(0.5)) {
      child1->colourLine()->addColoured(parent);
      child2->antiColourLine()->addAntiColoured(parent);
      ColinePtr temp = child1->antiColourLine();
      temp->addColoured(child2);
      child2->colourLine()->join(temp);
    }
    else {
      child2->colourLine()->addColoured(parent);
      child1->antiColourLine()->addAntiColoured(parent);
      ColinePtr temp = child2->antiColourLine();
      temp->addColoured(child1);
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
    throw Exception() << "Unknown colour in PowhegHandler::fixColours()"
		      << Exception::runerror;
  }
}

HardTreePtr PowhegHandler::generalClustering() {
  if(!_matrixElement) 
    throw Exception() << "PowhegHandler::generalClustering()"
		      << " must have a MatrixElement object for the core "
		      << "2->2 process" << Exception::runerror;
  PPair incoming = lastXCombPtr()->subProcess()->incoming();
  ParticleVector outgoing = lastXCombPtr()->subProcess()->outgoing();
  _s = lastXCombPtr()->lastS();
  // queue with the prototype trees
  std::queue<PrototypeTree> potentialTrees;
  // the base tree we'll make the others from
  PrototypeTree root;
  ShowerParticlePtr newParticle = 
    new_ptr(ShowerParticle(incoming.first->dataPtr(),false));
  newParticle->set5Momentum(incoming.first->momentum());
  root.incoming.insert(new_ptr(PrototypeBranching(newParticle)));
  newParticle = new_ptr(ShowerParticle(incoming.second->dataPtr(),false));
  newParticle->set5Momentum(incoming.second->momentum());
  root.incoming.insert(new_ptr(PrototypeBranching(newParticle)));
  for(set<PrototypeBranchingPtr>::const_iterator it=root.incoming.begin();
      it!=root.incoming.end();++it) {
    root.currentSpaceLike.insert(*it);
  }
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    newParticle = new_ptr(ShowerParticle(outgoing[ix]->dataPtr(),true));
    newParticle->set5Momentum(outgoing[ix]->momentum());
    root.outgoing.insert(new_ptr(PrototypeBranching(newParticle)));
  }
  potentialTrees.push(root);
  // store the final potential trees
  list<PrototypeTree> trees;
  while (!potentialTrees.empty()) {
    PrototypeTree current = potentialTrees.front();
    bool found(false);
    // potential final-final mergings
    set<PrototypeBranchingPtr>::iterator it,jt;
    for(it=current.outgoing.begin();it!=current.outgoing.end();++it) {
      jt = it;
      ++jt;
      for( ; jt!=current.outgoing.end();++jt) {
	pair<PrototypeBranchingPtr,PrototypeBranchingPtr> 
	  branch = make_pair(*it,*jt);
	BranchingElement allowed = allowedFinalStateBranching(branch);
	if(!allowed.first) continue;
	// copy the tree
	PrototypeTree newTree = current;
	map<PrototypeBranchingPtr,PrototypeBranchingPtr> pmap = newTree.reset();
	branch.first  = pmap[branch.first ];
	branch.second = pmap[branch.second];
	// make the new branching	
	// new particle first
	tcPDPtr newData = getParticleData(allowed.second[0]);
	Lorentz5Momentum newMomentum(branch.first ->particle->momentum()+
				     branch.second->particle->momentum());
	if(!newData->CC()||
	   (branch.first ->particle->id()==allowed.second[1]&&
	    branch.second->particle->id()==allowed.second[2])) {
	  newParticle = new_ptr(ShowerParticle(newData,true));
	}
	else {
	  newParticle = new_ptr(ShowerParticle(newData->CC(),true));
	}
	newParticle->set5Momentum(newMomentum);
	// then the branching
	PrototypeBranchingPtr newBranching(new_ptr(PrototypeBranching(newParticle)));
	branch.first ->parent =newBranching;
	branch.second->parent =newBranching;
	newBranching->children.push_back(branch.first );
	newBranching->children.push_back(branch.second);
	newBranching->sudakov = allowed.first;
	newTree.outgoing.erase(branch.first );
	newTree.outgoing.erase(branch.second);
	newTree.outgoing.insert(newBranching);
	// jet measure
	newTree.scales.push_back(hadronJetMeasure(branch.first ->particle->momentum(),
						  branch.second->particle->momentum()));
	// insert in the relevant list
	if(newTree.outgoing.size()==2) trees.push_back(newTree);
	else                           potentialTrees.push(newTree);
	found = true;
      }
    }
    // initial-final mergings
    for(it=current.outgoing.begin();it!=current.outgoing.end();++it) {
      for(jt=current.currentSpaceLike.begin();
	  jt!=current.currentSpaceLike.end();++jt) {
	pair<PrototypeBranchingPtr,PrototypeBranchingPtr> 
	  branch = make_pair(*jt,*it);
	BranchingElement allowed = allowedInitialStateBranching(branch);
	if(!allowed.first) continue;
	// copy the tree
	PrototypeTree newTree = current;
	map<PrototypeBranchingPtr,PrototypeBranchingPtr> pmap = newTree.reset();
	branch.first  = pmap[branch.first ];
	branch.second = pmap[branch.second];
	// make the new branching	
	// new particle first
	tcPDPtr newData = getParticleData(allowed.second[1]);
	Lorentz5Momentum newMomentum(branch.first ->particle->momentum()-
				     branch.second->particle->momentum());
	if(!newData->CC()||
	   (branch.first ->particle->id()==allowed.second[0]&&
	    branch.second->particle->id()==allowed.second[2])) {
	  newParticle = new_ptr(ShowerParticle(newData,false));
	}
	else {
	  newParticle = new_ptr(ShowerParticle(newData->CC(),false));
	}
	newParticle->set5Momentum(newMomentum);
	// then the branching
	PrototypeBranchingPtr newBranching(new_ptr(PrototypeBranching(newParticle)));
	newBranching->parent  = branch.first;
	branch.second->parent = branch.first;
	branch.first->children.push_back(newBranching);
	branch.first->children.push_back(branch.second);
	newBranching->parent->sudakov = allowed.first;
	newTree.currentSpaceLike.erase(branch.first );
	newTree.outgoing        .erase(branch.second);
	newTree.currentSpaceLike.insert(newBranching);
	// jet measure
	newTree.scales.push_back(hadronJetMeasure(branch.second->particle->momentum(),
						  branch.second->particle->momentum(),false));
	if(branch.first ->particle->momentum().z()/
	   branch.second->particle->momentum().z()>0.) newTree.scales.back()-=0.001*MeV2; 
	// insert in the relevant list
	if(newTree.outgoing.size()==2) trees.push_back(newTree);
	else                           potentialTrees.push(newTree);
	found = true;
      }
    }
    // treated one branching so pop from the queue
    if(!found) trees.push_back(current);
    potentialTrees.pop();
  }
  // check the core process is allowed using the matrix element
  // and reomve ones which aren't allowed
  list<PrototypeTree>::iterator it=trees.begin(),jt;
  while(it!=trees.end()) {
    DiagPtr diagram = getDiagram(*it);
    if(!diagram) it = trees.erase(it);
    else {
      it->diagram = diagram;
      ++it;
    }
  }
  // finally for the moment select the one with the smallest pt for the first branching
  // now find the one with the minimum pt
  HardTreePtr newTree;
  while(!trees.empty()) {
    jt=trees.end();
    Energy2 minkT =1e30*GeV2;
    for(it=trees.begin();it!=trees.end();++it) {
      if(it->scales.back()<minkT) {
	minkT = it->scales.back();
	jt=it;
      }
    }
    // construct the hard tree
    newTree = (*jt).convert();
    // assign the beam particles
    setBeams(newTree);
    // construct the colour flow
    createColourFlow(newTree,jt->diagram);
    // Calculate the shower variables
    evolver()->showerModel()->kinematicsReconstructor()->
      deconstructDecayJets(newTree,evolver());
    if(checkTree(newTree)) break; 
    trees.erase(jt);
  }
  // if no tree return an empty one
  if(trees.empty()) return HardTreePtr();
  // return the tree
  return newTree;
}

BranchingElement PowhegHandler::
allowedFinalStateBranching(pair<PrototypeBranchingPtr,PrototypeBranchingPtr> & br) {
  // check with normal ID's
  pair<long,long> ptest = make_pair(br.first->particle->id(),br.second->particle->id());
  map<pair<long,long>,pair<SudakovPtr,IdList> >::const_iterator 
    split = _allowedFinal.find(ptest);
  if(split!=_allowedFinal.end()) {
    if(split->second.second[1]!=ptest.first) swap(br.first,br.second);
    return split->second;
  }
  // check with CC
  if(br.first ->particle->dataPtr()->CC()) ptest.first  *= -1;
  if(br.second->particle->dataPtr()->CC()) ptest.second *= -1;
  _allowedFinal.find(ptest);
  if(split!=_allowedFinal.end()) {
    if(split->second.second[1]!=ptest.first) swap(br.first,br.second);
    return split->second;
  }
  // not found found null pointer
  return make_pair(SudakovPtr(),IdList());
}

BranchingElement PowhegHandler::
allowedInitialStateBranching(pair<PrototypeBranchingPtr,PrototypeBranchingPtr> & br) {
  // veto top
  if(abs(br.first ->particle->id())==ParticleID::t||
     abs(br.second->particle->id())==ParticleID::t)
    return make_pair(SudakovPtr(),IdList());
  bool cc = br.first->particle->id()<0;
  pair<multimap<long, pair<SudakovPtr,IdList> >::const_iterator,
    multimap<long, pair<SudakovPtr,IdList> >::const_iterator>
    location = _allowedInitial.equal_range(abs(br.first->particle->id()));
  for(multimap<long, pair<SudakovPtr,IdList> >::const_iterator it=location.first;
      it!=location.second;++it) {
    long idtest = it->second.second[2];
    if(cc&&getParticleData(idtest)->CC()) idtest *= -1;
    if(idtest==br.second->particle->id()) return it->second;
    if(idtest==-br.second->particle->id()&&
       !br.first->particle->dataPtr()->CC()) return it->second;
  }
  // not found found null pointer
  return make_pair(SudakovPtr(),IdList());
}

DiagPtr PowhegHandler::getDiagram(const PrototypeTree & tree) {
  // extract the incoming particles
  set<PrototypeBranchingPtr>::const_iterator it=tree.currentSpaceLike.begin();
  tcPDPair incoming;
  incoming.first  = (**it).particle->dataPtr();
  ++it;
  incoming.second = (**it).particle->dataPtr();
  // and the outgoing particles
  multiset<tcPDPtr> outgoing;
  for(it=tree.outgoing.begin();it!=tree.outgoing.end();++it)
    outgoing.insert((**it).particle->dataPtr());
  // see if the process is allowed
  for(MEBase::DiagramVector::const_iterator dt = _matrixElement->diagrams().begin();
      dt!=_matrixElement->diagrams().end();++dt) {
    const cPDVector partons=(**dt).partons();
    // check incoming particles
    if(!((incoming.first==partons[0]&&incoming.second==partons[1])||
	 (incoming.first==partons[1]&&incoming.second==partons[0]))) continue;
    // check the number of outgoing
    if(partons.size()!=tree.outgoing.size()+2) return DiagPtr();
    // check the outgoing
    multiset<tcPDPtr> otemp(outgoing);
    multiset<tcPDPtr>::iterator it;
    for(unsigned int ix=2;ix<partons.size();++ix) {
      it=otemp.find(partons[ix]);
      if(it!=otemp.end()) otemp.erase(it);
    }
    if(!otemp.empty()) continue;
    return *dt;
  }
  return DiagPtr();
}

Energy2 PowhegHandler::hadronJetMeasure(const Lorentz5Momentum & p1,
					const Lorentz5Momentum & p2,
					bool final) {
  Energy2 output;
  if(final) {
    double deltay   = p1.rapidity()-p2.rapidity();
    double deltaphi = p1.phi()-p2.phi();
    if(deltaphi<-Constants::pi) deltaphi += Constants::twopi;
    if(deltaphi> Constants::pi) deltaphi -= Constants::twopi;
    double deltaR = sqr(deltay)+sqr(deltaphi);
    output = min(p1.perp2(),p2.perp2())*deltaR;
  }
  else {
    output = p1.perp2();
  }
  return output;
}

HardBranchingPtr PrototypeBranching::convert() {
  if(!particle) {
    cerr << "testing don't have particle for the branching shit" << "\n";
    exit(0);
  }
  // create the new particle
  HardBranchingPtr hard=new_ptr(HardBranching(particle,sudakov,
					      tHardBranchingPtr(),
					      !particle->isFinalState()));
  // and the children
  for(unsigned int ix=0;ix<children.size();++ix) {
    hard->addChild(children[ix]->convert());
    hard->children().back()->parent(hard);
  }
  return hard;
}

HardTreePtr PrototypeTree::convert() {
  vector<HardBranchingPtr> branchings,spacelike;
  set<PrototypeBranchingPtr>::const_iterator it,jt;
  // incoming lines and spacelike inot the hard process
  for(it=incoming.begin();it!=incoming.end();++it) {
    spacelike.push_back((**it).convert());
    HardBranchingPtr br(spacelike.back());
    while (!br->children().empty()) {
      for(unsigned int ix=0;ix<br->children().size();++ix) {
	if(br->children()[ix]->incoming()) {
	  br = br->children()[ix];
	  break;
	}
      }
    }
    branchings.push_back(br);
  }
  // outgoing particles
  for(it=outgoing.begin();it!=outgoing.end();++it) {
    branchings.push_back((**it).convert());
  }
  HardTreePtr newTree = new_ptr(HardTree(branchings,spacelike));
  return newTree;
}

map<PrototypeBranchingPtr,PrototypeBranchingPtr> PrototypeTree::reset() {
  map<PrototypeBranchingPtr,PrototypeBranchingPtr> output;
  set<PrototypeBranchingPtr> newOutgoing;
  set<PrototypeBranchingPtr> newIncoming;
  set<PrototypeBranchingPtr> newSpaceLike;
  set<PrototypeBranchingPtr>::iterator it,jt;
  for(it=incoming.begin();it!=incoming.end();++it) {
    PrototypeBranchingPtr newBr = (**it).reset(PrototypeBranchingPtr(),output); 
    newIncoming.insert(newBr);
    PrototypeBranchingPtr br=newBr;
    while(!br->children.empty()) {
      for(unsigned int ix=0;ix<br->children.size();++ix) {
	if(!br->children[ix]->particle->isFinalState()) {
	  br = br->children[ix];
	  break;
	}
      }
    }
    newSpaceLike.insert(br);
  }
  for(it=outgoing.begin();it!=outgoing.end();++it) {
    newOutgoing.insert((**it).reset(PrototypeBranchingPtr(),output));
  }
  outgoing  = newOutgoing;
  incoming  = newIncoming;
  currentSpaceLike = newSpaceLike;
  return output;
}

PrototypeBranchingPtr PrototypeBranching::
reset(PrototypeBranchingPtr newParent,
	map<PrototypeBranchingPtr,PrototypeBranchingPtr> & pmap) {
  PrototypeBranchingPtr output(new_ptr(PrototypeBranching(particle)));
  pmap[this] = output;
  output->sudakov  = sudakov;
  output->parent   = newParent;
  for(unsigned int ix=0;ix<children.size();++ix) {
    output->children.push_back(children[ix]->reset(output,pmap));
  }
  return output;
}

void PowhegHandler::createColourFlow(HardTreePtr tree,
				     DiagPtr diagram) {
  // first construct a set of on-shell momenta for the hard collison
  vector<Lorentz5Momentum> meMomenta;
  vector<tcPDPtr> mePartonData;
  PVector particles;
  set<HardBranchingPtr>::const_iterator it; 
  for(it=tree->branchings().begin();it!=tree->branchings().end();++it) {
    if((**it).incoming()) {
      meMomenta.push_back((**it).branchingParticle()->momentum());
      mePartonData.push_back((**it).branchingParticle()->dataPtr());
      particles.push_back((**it).branchingParticle());
    }
  }
  for(it=tree->branchings().begin();it!=tree->branchings().end();++it) {
    if(!(**it).incoming()) {
      meMomenta.push_back((**it).branchingParticle()->momentum());
      mePartonData.push_back((**it).branchingParticle()->dataPtr());
      particles.push_back((**it).branchingParticle());
    }
  }
//   cerr << "testing number of partons\n";
//   for(unsigned int ix=0;ix<meMomenta.size();++ix) {
//     cerr << *particles[ix] << "\n";
//   }
  // boost the momenta to the CMF frame
  // compte boost to reset frame
  Lorentz5Momentum prest(meMomenta[0]+meMomenta[1]);
  LorentzRotation R(-prest.boostVector());
  // and then to put beams along the axis
  Lorentz5Momentum ptest = R*meMomenta[0];
  Axis axis(ptest.vect().unit());
  if(axis.perp2()>0.) {
    R.rotateZ(-axis.phi());
    R.rotateY(-acos(axis.z()));
  }
  const cPDVector partons=diagram->partons();
  // order of the incoming partons
  if(mePartonData[0]!=partons[0]) {
    swap(mePartonData[0],mePartonData[1]);
    swap(meMomenta[0],meMomenta[1]);
    swap(particles[0],particles[1]);
  }
  // order of the outgoing partons
  for(unsigned int ix=2;ix<partons.size();++ix) {
    for(unsigned int iy=ix;iy<meMomenta.size();++iy) {
      if(partons[ix]==mePartonData[iy]) {
	if(ix!=iy) {
	  swap(mePartonData[ix],mePartonData[iy]);
	  swap(meMomenta[ix],meMomenta[iy]);
	  swap(particles[ix],particles[iy]);
	}
	break;
      }
    }
  }
  for(unsigned int ix=0;ix<meMomenta.size();++ix)
    meMomenta[ix].transform(R);
  PPair in(mePartonData[0]->produceParticle(meMomenta[0]),
	   mePartonData[1]->produceParticle(meMomenta[1]));
  PVector out;
  for(unsigned int ix=2;ix<meMomenta.size();++ix) {
    out.push_back(mePartonData[ix]->produceParticle(meMomenta[ix]));
  }
  _matrixElement->setKinematics(in,out);
  _matrixElement->dSigHatDR();
  const ColourLines & cl = _matrixElement->selectColourGeometry(diagram);
  PVector slike;
  tPVector ret;
  slike.push_back(particles[0]);
  Ptr<Tree2toNDiagram>::pointer diagram2 = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::pointer>(diagram);
  for ( int i = 1; i < diagram2->nSpace() - 1; ++i )
    slike.push_back(diagram2->allPartons()[i]->produceParticle());
  slike.push_back(particles[1]);
  ret = tPVector(slike.begin(), slike.end());
  int io = particles.size();
  PVector tlike(diagram2->allPartons().size() - diagram2->nSpace());
  for ( int i = diagram2->allPartons().size() - 1; i >=  diagram2->nSpace(); --i ) {
    int it = i - diagram2->nSpace();
    pair<int,int> ch = diagram2->children(i);
    bool iso = ch.first < 0;
    if ( iso ) {
      tlike[it] = particles[--io];
    } 
    else {
      Lorentz5Momentum p = tlike[ch.first - diagram2->nSpace()]->momentum() +
 	tlike[ch.second - diagram2->nSpace()]->momentum();
      tlike[it] = diagram2->allPartons()[i]->produceParticle(p);
    }
  }
  ret.insert(ret.end(), tlike.begin(), tlike.end());
  cl.connect(ret);
  for(unsigned int ix=0;ix<ret.size();++ix) {
    PVector::iterator it = find(particles.begin(),particles.end(),ret[ix]);
    if(it==particles.end()) {
      ColinePtr line = ret[ix]->colourLine();
      if(line) line->removeColoured(ret[ix]);
      line = ret[ix]->antiColourLine();
      if(line) line->removeAntiColoured(ret[ix]);
    }
  }
  // now the colours of the rest of the particles
  for(set<HardBranchingPtr>::const_iterator it=tree->branchings().begin();
      it!=tree->branchings().end();++it) (**it).fixColours();
}

double PowhegHandler::Sud( Energy scale, long id ){
  //upper limit on scale 
  double sudwgt = 1.;
  Energy scale_cut = 1000*GeV;
  multimap< long, pair < Interpolator<double,Energy>::Ptr,
    Interpolator<Energy,double>::Ptr >  >::const_iterator cjt;
  for( cjt =  _fbranchings.lower_bound( abs( id ) );
       cjt != _fbranchings.upper_bound( abs( id ) );
       ++cjt ) {
    if( scale < scale_cut ) sudwgt *= (* cjt->second.first )( scale );
    else sudwgt *= (* cjt->second.first )( scale_cut );
  }
  return sudwgt;
}

double PowhegHandler::sudakovWeight() {
  double SudWgt = 1.;
  //include the sud factor for each external line
  for( map<ShowerParticlePtr,pair<double,HardBranchingPtr> >::const_iterator cit = _theExternals.begin();
       cit != _theExternals.end(); ++cit ) {
  
    Energy scale= sqrt( _s );
    if( cit->second.second ){
      if( cit->first == 
	  cit->second.second->children()[0]->branchingParticle() )
	scale =  cit->second.second->scale()* 
	  cit->second.second->children()[0]->z();
      else if(  cit->first == 
		cit->second.second->children()[1]->branchingParticle() )
	scale =  cit->second.second->scale()* 
	  cit->second.second->children()[1]->z();
      else cerr<<"could not find child in external HardBranching \n";
    }
    else{
      scale = cit->first->evolutionScale();
    }
    SudWgt *= Sud( scale, cit->first->id() );  
  }
  
  if(SudWgt > 1.1) cerr<<"\n\n\nsudakov from externals > 1!!\n\n\n";
  //include the intermediate line wgts
  for( map< long, pair< pair< double, Energy >, pair< double, Energy > > >::const_iterator cit 
	 = _theIntermediates.begin();
       cit != _theIntermediates.end(); ++cit ) {
    
    Energy scale =  cit->second.first.second;
    
    double internal_wgt = Sud( scale, cit->first );
    scale =  cit->second.second.second;
    if(scale > cit->second.first.second ) cerr <<"scales wrong way round!\n";
    internal_wgt /= Sud( scale, cit->first );
    if(internal_wgt > 1.1 || internal_wgt < 0.)cerr<<"\n\nbig internal weight of "<< internal_wgt
						   <<"\nnum scale = "
						   <<cit->second.first.second / GeV
						   <<"\nden scale = "
						   <<cit->second.second.second /GeV
						   <<"\n\n";
  }
 
  double alphaWgt = 1.;
  //need to add the alphaS weight
  //the alphaS ratio evaluated at all nodal values
  for(map<HardBranchingPtr,double>::const_iterator cit=_theNodes.begin();
      cit != _theNodes.end(); ++cit ) {
    alphaWgt *= _alphaS->value( cit->second * _s ) / _alphaSMG;
  }
  (*_hSud) += SudWgt;
  (*_halphaS) += alphaWgt;
  if( SudWgt > 1.1 ) {
    cerr<<"\n\nweight exceeded 1 in PowhegHandler::reweight() !!! \n";
    cerr<<"  alpha wgt = "<<alphaWgt
    	<<"\n  sudWgt = "<<SudWgt<<"\n\n";
  }
  return SudWgt*alphaWgt;
}

void PowhegHandler::setBeams(HardTreePtr tree) {
  PPair beams=lastXCombPtr()->lastParticles();
  if((**tree->incoming().begin()).branchingParticle()->momentum().z()/
     beams.first->momentum().z()<0.)
    swap(beams.first,beams.second);
  set<HardBranchingPtr>::iterator it = tree->incoming().begin();
  HardBranchingPtr br=*it;
  br->beam(beams.first);
  while (!br->children().empty()) {
    for(unsigned int ix=0;ix<br->children().size();++ix) {
      if(br->children()[ix]->incoming()) {
	br = br->children()[ix];
	break;
      }
    }
    br->beam(beams.first);
  }
  ++it;
  br=*it;
  br->beam(beams.second);
  while (!br->children().empty()) {
    for(unsigned int ix=0;ix<br->children().size();++ix) {
      if(br->children()[ix]->incoming()) {
	br = br->children()[ix];
	break;
      }
    }
    br->beam(beams.second);
  }
}

bool PowhegHandler::checkTree(HardTreePtr tree) {
  set<HardBranchingPtr>::const_iterator it;
  bool reject = false;
  for(it=tree->incoming().begin();it!=tree->incoming().end();++it) {
    reject |=checkBranching(*it);
  }
  for(it=tree->branchings().begin();it!=tree->branchings().end();++it) {
    if((**it).incoming()) continue;
    reject |=checkBranching(*it);
  }
  return !reject;
}

bool PowhegHandler::checkBranching(HardBranchingPtr br) {
  static const double eps(1e-5);
  bool reject(false);
  for(vector<HardBranchingPtr>::const_iterator it=br->children().begin();
      it!=br->children().end();++it) {
    reject |=checkBranching(*it);
  }
  reject |= br->z()<-eps || br->z()>1.+eps;
  return reject;
}
