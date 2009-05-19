// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegHandler class.
//

#include "PowhegHandler.h"
#include "ThePEG/Utilities/CFileLineReader.h"
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
#include "ThePEG/Repository/UseRandom.h"


using namespace Herwig;

IBPtr PowhegHandler::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegHandler::fullclone() const {
  return new_ptr(*this);
}

void PowhegHandler::persistentOutput(PersistentOStream & os) const {
  os  << _alphaS << _sudopt << _sudname << _jetMeasureMode << _allowedInitial
      << _allowedFinal << _matrixElement << _lepton << _highestMult << _reweightOpt << _testSudakovs <<_yini 
      << _alphaSMG << _npoint <<  ounit( _max_qtilde, GeV ) <<  ounit( _max_pt_cut, GeV ) 
      <<  ounit( _min_pt_cut, GeV ) << _clusterOption << _dalitzOn << _qtildeDist;
}

void PowhegHandler::persistentInput(PersistentIStream & is, int) {
  is  >> _alphaS >> _sudopt >> _sudname >> _jetMeasureMode >> _allowedInitial
      >> _allowedFinal >> _matrixElement >> _lepton >> _highestMult >> _reweightOpt >> _testSudakovs >> _yini 
      >> _alphaSMG >> _npoint >> iunit( _max_qtilde, GeV ) >> iunit( _max_pt_cut, GeV ) 
      >> iunit( _min_pt_cut, GeV ) >> _clusterOption >> _dalitzOn >> _qtildeDist;
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
     &PowhegHandler::_jetMeasureMode, 1, false, false);
  
  static SwitchOption Durham
    (ifaceJetMeasureMode,"Durham","Durham jet algorithm", 0);
  
  static SwitchOption ShowerPt
    (ifaceJetMeasureMode,"ShowerPt","ShowerPt", 1);
  
  static SwitchOption LUCLUS
    (ifaceJetMeasureMode,"LUCLUS","LUCLUS jet algorithm", 2);

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

   static Switch<PowhegHandler,bool> interfaceHighestMultiplicity
    ("HighestMult",
     "Whether we are treating the highest mutliplicity treatment",
     &PowhegHandler::_highestMult, false, false, false);

  static SwitchOption interfaceMultHighest
    (interfaceHighestMultiplicity,
     "Highest",
     "highest multiplicity",
     true);

  static SwitchOption interfaceMultNotHighest
    (interfaceHighestMultiplicity,
     "NotHighest",
     "Not the highest multiplicity",
     false);

  static Switch<PowhegHandler, unsigned int> interfaceReweight
    ("ReweightOption",
     "Whether to switch off the sudakov reweighting",
     &PowhegHandler::_reweightOpt, 0, false, false);
  static SwitchOption interfaceReweightOff
    (interfaceReweight,
     "Off",
     "No Sudakov reweighting",
     1);
  static SwitchOption interfaceReweightOn
    (interfaceReweight,
     "On",
     "Do Sudakov reweighting",
     0);
  static SwitchOption interfaceReweightNoSud
    (interfaceReweight,
     "NoSud",
     "alphaS reweighting only",
     2);

  static Switch<PowhegHandler,bool> interfaceTestSudakov
    ("testSudakov",
     "Whether to output Sudakov test histograms",
     &PowhegHandler::_testSudakovs, false, false, false);

  static SwitchOption interfaceTestSudakovOff
    (interfaceTestSudakov,
     "Off",
     "No Sudakov testing",
     false);
  static SwitchOption interfaceTestSudakovOn
    (interfaceTestSudakov,
     "On",
     "Do Sudakov testing",
     true);


  static Reference<PowhegHandler,MEBase> interfaceMatrixElement
    ("MatrixElement",
     "The matrix element class for the core 2->2 process",
     &PowhegHandler::_matrixElement, false, false, true, true, false);

  static Parameter<PowhegHandler,unsigned int> interfaceInterpPoints
    ("InterpolatorPoints",
     "The number of points used for sudakov interpolation tables",
     &PowhegHandler::_npoint, 10, 0, 1000000,
     false, false, Interface::limited );

  static Parameter<PowhegHandler, Energy> interfaceMaxQTilde
    ("maxQTilde",
     "The maximum QTilde scale for sudakov interpolation tables",
     &PowhegHandler::_max_qtilde, GeV, 91.2*GeV, 1.*GeV, 1000000.*GeV,
     false, false, Interface::limited);

  static Parameter<PowhegHandler, Energy> interfaceMaxPtCut
    ("maxPtCut",
     "The maximum pt cut for sudakov interpolation tables",
     &PowhegHandler::_max_pt_cut, GeV, 45.6*GeV, 1.*GeV, 1000000.*GeV,
     false, false, Interface::limited);

  static Parameter<PowhegHandler, Energy> interfaceMinPtCut
    ("minPtCut",
     "The minimum pt cut for sudakov interpolation tables",
     &PowhegHandler::_min_pt_cut, GeV, 0.*GeV, 0.*GeV, 1000000.*GeV,
     false, false, Interface::limited);

  static Switch<PowhegHandler, unsigned int> ifaceClusterOption
    ("ClusterOption",
     "Choice of the clustering scheme",
     &PowhegHandler::_clusterOption, 0, false, false);
  
  static SwitchOption allHistories
    (ifaceClusterOption,"allHistories","make all histories and require angular ordering", 0);
  
  static SwitchOption jetClustered
    (ifaceClusterOption,"jetClustered", "cluster according to jet algorithm", 1);
  
  static SwitchOption ptChoice
    (ifaceClusterOption,"ptChoice", "choose ordered history with lowest total pt", 2);
  
  static SwitchOption highProbChoice
    (ifaceClusterOption,"highProbChoice", "choose ordered history with highest probability", 3);

  static Switch<PowhegHandler,bool> interfaceDalitz
    ("Dalitz",
     "Switch for dalitz analysis of hard tree clustering (3 jets only)",
     &PowhegHandler::_dalitzOn, false, false, false);
  static SwitchOption interfaceDalitzOff
    (interfaceDalitz,
     "Off",
     "Dalitz analysis off",
     false);
  static SwitchOption interfaceDalitzOn
    (interfaceDalitz,
     "On",
     "Dalitz analysis on",
     true);

  static Switch<PowhegHandler,bool> interfaceQtildeDist
    ("QtildeDist",
     "Switch for qtilde distribution hists from sudakov tables",
     &PowhegHandler::_qtildeDist, false, false, false);
  static SwitchOption interfaceQtildeDistOff
    (interfaceQtildeDist,
     "Off",
     "Qtilde analysis off",
     false);
  static SwitchOption interfaceQtildeDistOn
    (interfaceQtildeDist,
     "On",
     "Qtilde analysis on",
     true);
}

double PowhegHandler::reweightCKKW(int minMult, int maxMult) {
  // cluster the event
  _max_mult = maxMult;

  //calculate a global alphaS factor to ensure alphaS weight is
  //always less than one.
  double alphaS_max =  _alphaS->value( sqr( _fbranchings.find( long( 1 ) )->second.second ) );
  _global_alphaS_wgt = 1.;
  for( int ix = minMult; ix < maxMult; ix++ ) 
    _global_alphaS_wgt *= alphaS_max / _alphaSMG;

  if( _clusterOption == 0 || _clusterOption == 2 || _clusterOption == 3 )
    _theHardTree = doClusteringOrdered();
  else
    _theHardTree = doClustering();

  //if highest mult set veto def
  // needs fixing
//   if( _highestMult == true ) evolver()->setHighest( true );
//   else evolver()->setHighest( false );


  // return if fails
  if(!_theHardTree)
    return 0.;
  //call dalitz analysis if asked for
  if( _dalitzOn ) getDalitz();

  // compute the Sudakov weight
  double SudWgt;
  if(  _reweightOpt != 1 )
    SudWgt = _lepton ? sudakovWeight( _theHardTree ) : 1.;
  else SudWgt = 1.;
 
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
      if((**it).status()==HardBranching::Incoming) continue;
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
      if((**it).status()==HardBranching::Incoming) {
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
  //divide by constant 5 here is just a factor becasue alphas = 0.2 isn't high enough to guarrantee
  //that the alphaS weight is < 1
  if( _reweightOpt != 1 ){
    if ( SudWgt / 2. / _global_alphaS_wgt > 1 ) cerr << SudWgt / 2. / _global_alphaS_wgt <<"\n";
    return SudWgt / 2. / _global_alphaS_wgt;
  }
  
  else{
    return 1.;
  }
}
void PowhegHandler::dofinish() {
  ShowerHandler::dofinish();
  
  if( _clusterOption == 1 ){
    cout<<"\n---------\n"
	<<"proportion of unordered trees created = "
	<< ( 1. - double( _ordered_trees_created )
	     / double( _trees_created ) ) * 100.
	<<" %\n--------\n";
  }
  
  string fname = generator()->filename() + string("-") 
    + string("wgts.top");
  ofstream output(fname.c_str());
  //output ddalitz analysis if switched on
  if( _dalitzOn ){
    string dalitz_fname = generator()->filename() + string("-") + string("dalitz.top");
    ofstream dalitz( dalitz_fname.c_str() );
    dalitz<<"SET WINDOW X 2 9 Y 2 9\n";
    dalitz<<"SET FONT DUPLEX\n";
    dalitz<<"SET LIMITS X 0 1 Y 0 1\n";
    dalitz<<"TITLE BOTTOM \"X011\" \n";
    dalitz<<"CASE         \" X X\" \n";
    dalitz<<"TITLE LEFT \"X021\" \n";
    dalitz<<"CASE       \" X X\" \n";
    for(unsigned int ix = 0; ix < _dalitz_from_q1.size(); ix++ )
      dalitz<< _dalitz_from_q1[ix].first <<"\t"<< _dalitz_from_q1[ix].second <<"\n";
    dalitz << "PLOT RED\n";
    for(unsigned int ix = 0; ix < _dalitz_from_q2.size(); ix++ )
      dalitz<< _dalitz_from_q2[ix].first <<"\t"<< _dalitz_from_q2[ix].second <<"\n";
    dalitz << "PLOT BLUE\n";
    
    dalitz<<"NEW FRAME \n";
    dalitz<<"SET WINDOW X 2 9 Y 2 9\n";
    dalitz<<"SET FONT DUPLEX\n";
    dalitz<<"SET LIMITS X 0 1 Y 0 1\n";
    dalitz<<"TITLE BOTTOM \"X011\" \n";
    dalitz<<"CASE         \" X X\" \n";
    dalitz<<"TITLE LEFT \"X021\" \n";
    dalitz<<"CASE       \" X X\" \n";
    for(int ix = 0; ix < _dalitz_from_q1.size(); ix++ )
      dalitz<< _dalitz_from_q1[ix].first <<"\t"<< _dalitz_from_q1[ix].second <<"\n";
    dalitz << "PLOT RED\n";
    
    dalitz<<"NEW FRAME \n";
    dalitz<<"SET WINDOW X 2 9 Y 2 9\n";
    dalitz<<"SET FONT DUPLEX\n";
    dalitz<<"SET LIMITS X 0 1 Y 0 1\n";
    dalitz<<"TITLE BOTTOM \"X011\" \n";
    dalitz<<"CASE         \" X X\" \n";
    dalitz<<"TITLE LEFT \"X021\" \n";
    dalitz<<"CASE       \" X X\" \n";
    for(int ix = 0; ix < _dalitz_from_q2.size(); ix++ )
      dalitz<< _dalitz_from_q2[ix].first <<"\t"<< _dalitz_from_q2[ix].second <<"\n";
    dalitz << "PLOT BLUE\n";
 
    cerr<<"no q1 dalitz = "<<_dalitz_from_q1.size()<<"\n";
    cerr<<"no q2 dalitz = "<<_dalitz_from_q2.size()<<"\n";

    dalitz.close();
  }
  
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
  _trees_created = 0;
  _ordered_trees_created = 0;

  _dalitz_from_q1.clear();
  _dalitz_from_q2.clear();

  ShowerHandler::doinitrun();  
  _s = sqr( generator()->maximumCMEnergy() );

  _hSud = new_ptr(Histogram(0.,2.,100));
  _halphaS = new_ptr(Histogram(0.,2.,100));

  ofstream sudFileOutput;
  if(_sudopt==1) sudFileOutput.open(_sudname.c_str());
 
  // integrator for the outer integral
  GaussianIntegrator outer;
  // get the final-state branchings from the evolver
  if(_sudopt!=2) {
   
    for(BranchingList::const_iterator 
	  it = evolver()->splittingGenerator()->finalStateBranchings().begin();
	it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
   
      //skip sudakovs involving tops
      if( abs( it->second.second[0] ) == 6 ||
	  abs( it->second.second[1] ) == 6 ||
	  abs( it->second.second[2] ) == 6 ) continue;

      Ptr<QTildeSudakovIntegrator>::pointer integrator = 
	new_ptr( QTildeSudakovIntegrator(it->second, _jetMeasureMode, _s ) );
    
      Energy qtildemax = _max_qtilde;
      Energy qtildemin = integrator->minimumScale();

      //initialise sudakov values on grid ij ( pt_i, qtilde_j )
      vector< double > dummy( _npoint + 1, 0. );
      vector< vector< double > > sud( _npoint + 1, dummy );
      vector< Energy > ptCut;
      vector< Energy > scale;
  
      //fill scales at start
      for( unsigned int ix = 0; ix < _npoint + 1; ++ix ){
	ptCut.push_back( _min_pt_cut + double( ix ) * ( _max_pt_cut - _min_pt_cut ) / double( _npoint - 1 ) );
	scale.push_back( qtildemin + double( ix ) * ( qtildemax - qtildemin ) / double( _npoint - 1 ) );
      }
      //fill sud integrals
      for( unsigned int ix = 0; ix < _npoint + 1; ++ix ){
	sud[ix][0] = 0.; 
	for( unsigned int jx = 1; jx < _npoint + 1; ++jx ) {
	  //the pt_cut here is pt in the jet measure variable used
	  double currentSud = integrator->value( scale[ jx ], scale[ jx - 1 ], ptCut[ix] );
	  sud[ix][jx] = ( sud[ix][ jx - 1 ]  + currentSud );
	  // cerr<<jx<<"\t";
	}
	//	cerr<<ix<<"\n";
      }

      cerr<<"\n\nextreme scales pt, qt = "<< ptCut[ _npoint  ] / GeV <<" "<<  scale[ _npoint ] / GeV <<"\n\n\n";
      //exponentiate to the Sudakov
      for( unsigned int ix = 0; ix < _npoint+1; ++ix ) {
	for( unsigned int jx = 0; jx < _npoint+1; ++jx ) {
	  sud[ix][jx] = exp( - sud[ix][jx] );
	}
      }

      Interpolator2d< double, Energy, Energy >::Ptr theInterpolator = 
	new_ptr( Interpolator2d<double,Energy,Energy>( sud, ptCut, scale ) );

      _fbranchings.insert( make_pair( it->first, make_pair( theInterpolator, qtildemin ) ) );
      
      //write current sud grid to selected output file
      if(_sudopt==1) {
	
	sudFileOutput << it->second.second[0] << "\t"
		      << it->second.second[1] << "\t"
		      << it->second.second[2] << "\n";

	//output the grid size i * j
	sudFileOutput << ptCut.size() << "\t" << scale.size() << "\n";

	for( unsigned int jx = 0;jx < scale.size(); ++jx ){
	  //write a row
	  for( unsigned int ix = 0; ix < ptCut.size(); ++ix ){
	    sudFileOutput << sud[ix][jx] << "\t";
	  }
	  sudFileOutput << "\n";
	}
	//output pts in a line
	for( unsigned int ix = 0; ix < ptCut.size(); ++ix ){
	  sudFileOutput << ptCut[ix] / GeV << "\t";
	}
	sudFileOutput << "\n";
	//output scales in a line
	for( unsigned int ix = 0; ix < scale.size(); ++ix ){
	  sudFileOutput << scale[ix] / GeV << "\t";
	}
	sudFileOutput << "\n";
      }
      
    }
    sudFileOutput.close();
  }
 
  else {
    CFileLineReader file(_sudname);
    while(file.readline()) {
      string line = file.getline();
      istringstream is;
      is.str(line);
      IdList ids(3);
      //GET NAMES AND IDS FROM FIRST LINE 
      is >> ids[0] >> ids[1] >> ids[2];
 
      file.readline();
      unsigned int isize;
      unsigned int jsize;

      //GET THE NUMBER OF POINTS IN THIS
      is.str(file.getline());
      is >> isize >> jsize;

      //initialise vectors and matrix to the correct size
      vector< double > dummy( jsize, 0. );
      vector< vector< double > > sud( isize, dummy );
      vector< Energy > pt( isize );
      vector< Energy > scale( jsize );
      
      //read in matrix of sud values sud( pt_i, scale_i )

      //read a column- different line is different qtilde
      for( unsigned int jx = 0; jx < jsize; ++jx ) {
	//read a horizontal line diferent entry ids different pt
	file.readline();
	is.str(file.getline());
	for( unsigned int ix = 0; ix < isize; ++ix )
	  is >> sud[ix][jx];
      }
      //read in pt values
      file.readline();
      is.str(file.getline());
      for( unsigned int ix = 0; ix < isize; ++ix ){
	double val;
	is >> val;
	pt[ix] = val * GeV;
      }
      file.readline();
      is.str(file.getline());
      for( unsigned int jx = 0; jx < jsize; ++jx ){
	double val;
	is >> val;
	scale[jx] = val * GeV;
      }
    
      //find branching list with matching ids
      BranchingList::const_iterator it,
	start = evolver()->splittingGenerator()->finalStateBranchings().lower_bound(ids[0]),
	end   = evolver()->splittingGenerator()->finalStateBranchings().upper_bound(ids[0]);
      for( it = start; it != end; ++it ) {
	if( it->second.second[0] == ids[0] &&
	    it->second.second[1] == ids[1] &&
	    it->second.second[2] == ids[2] ) {
	  // construct the Interpolators
	  Interpolator2d< double, Energy, Energy >::Ptr theInterpolator = 
	    new_ptr( Interpolator2d< double, Energy, Energy >( sud, pt, scale ) );
	  Energy qtildemin = scale[0];
	  _fbranchings.insert( make_pair( it->first, make_pair( theInterpolator, qtildemin ) ) );
	  break;
	}
	
	if( it == end ) {
	  cerr << "sud read error: could not find correct branching list\n";
	}
      }
    }
  }

  if( _testSudakovs ) testSudakovs();
  if( _qtildeDist ) makeQtildeDist();
}

void PowhegHandler::doinit() {
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
  case 2:
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
bool PowhegHandler::splittingAllowed( ShowerParticlePtr part_i,
				      ShowerParticlePtr part_j,
				      int qq_pairs ) {
  // g 2 q qbar or an incorrect qq type
  if ( abs ( part_i->id() ) < 7 && abs ( part_j->id() ) < 7 ) { 
    if ( abs ( part_i->id() ) != abs ( part_j->id() ) ) return false;
    if ( ( part_i->id() < 0 &&  part_j->id() < 0 ) ||
	 ( part_i->id() > 0 &&  part_j->id() > 0 ) ) return false;
    if ( qq_pairs < 2 ) return false;
  }
  return true;
}

// finds the id of the emitting particle and sudakov for the desired clustering
// also swaps order of children pointers as required (not pointers passed by reference)
SudakovPtr PowhegHandler::getSud( long & emmitter_id,
				  ShowerParticlePtr & part_i, 
				  ShowerParticlePtr & part_j ) {
  // g 2 q qbar or an incorrect qq type
  if ( abs ( part_i->id() ) < 7 && abs ( part_j->id() ) < 7 ) { 
    if ( abs ( part_i->id() ) != abs ( part_j->id() )  ) return SudakovPtr();
    if ( ( part_i->id() < 0 &&  part_j->id() < 0 ) ||
	 ( part_i->id() > 0 &&  part_j->id() > 0 ) ) return SudakovPtr();
    //if the q and qbar are the wrong way round then switch order
    if ( part_j->id() > part_i->id() ) swap( part_i, part_j );
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
      if( abs( ids[1] ) == abs( part_j->id() ) && 
	  abs( ids[2] ) == abs( part_i->id() ) ) {
	swap( part_i, part_j );
	return cit->second.first;
      }
    }
  }
  return SudakovPtr();
}
//have the & in here so that can remove pointer - is this right??
bool PowhegHandler::fillProtoTrees( map< ShowerParticlePtr, HardBranchingPtr > theParticles, 
				     ProtoTreePtr currentProtoTree ){
  if( theParticles.size() < 3 ) return true;
  //find number of qqpairs - count quarks - used
  int no_qqbar = 0;
  for( map<ShowerParticlePtr, HardBranchingPtr>::iterator itc = theParticles.begin();
       itc != theParticles.end(); itc++ ) 
    if( itc->first->id() > 0 && itc->first->id() < 7 ) no_qqbar ++;

  HardBranchingPtr currentBranching;
  ProtoTreePtr newProtoTree;
  map< ShowerParticlePtr, HardBranchingPtr > newParticles;
  for( map< ShowerParticlePtr, HardBranchingPtr >::iterator ita = theParticles.begin();
       ita != theParticles.end(); ita++ ) {
    for( map< ShowerParticlePtr, HardBranchingPtr >::iterator itb = theParticles.begin();
	 itb != ita; itb++) {
      if( ! splittingAllowed( ita->first, itb->first, no_qqbar )  ) continue;
     
      currentBranching = getCluster( make_pair( ita->first, itb->first ), 
				     theParticles );
 
      set< HardBranchingPtr > newTreeBranchings = currentProtoTree->getBranchings();
      
      newProtoTree = new_ptr( ProtoTree( newTreeBranchings ) );
      //now remove the hard branchings of clustered particles from newTree
      
      //currentProtoTree should contain which ever hardbranchings are in theParticles
      if( ! newProtoTree->removeBranching( ita->second ) )
	cerr<<"fill proto tree problem!!! can't find clustered in newProtoTree 1 \n"
	    << "couldn't find: " << ita->second <<"\n";
      if( ! newProtoTree->removeBranching( itb->second ) )
	cerr<<"fill proto tree problem!!! can't find clustered in newProtoTree 1 \n"
	    << "couldn't find: " << itb->second <<"\n";
     
      newProtoTree->addBranching( currentBranching );
  
      newParticles = theParticles;
      if( newParticles.find( ita->first ) != newParticles.end() )
	newParticles.erase( ita->first );
      else
	cerr<<"fill newParticles problem!!! can't find clustered in newParticles \n";
      if( newParticles.find( itb->first ) != newParticles.end() )
	  newParticles.erase( itb->first );
      else
	cerr<<"fill newParticles problem!!! can't find clustered in newParticles \n";
      newParticles.insert( make_pair( currentBranching->branchingParticle(),
				      currentBranching ) );
    
      if( ! repeatProtoTree( newProtoTree ) ) _proto_trees.insert( newProtoTree );
     
      //remove the current tree if it hasn't already been removed
      if( _proto_trees.find( currentProtoTree ) != _proto_trees.end() )
	_proto_trees.erase( currentProtoTree );
      //do recursion
      fillProtoTrees( newParticles, newProtoTree );
    }
  }
  return true;
}

HardBranchingPtr PowhegHandler::getCluster( pair< ShowerParticlePtr, ShowerParticlePtr > clusterPair,
					    map< ShowerParticlePtr, HardBranchingPtr > theParticles ){
  for( map< HardBranchingPtr , pair< ShowerParticlePtr, ShowerParticlePtr > >::const_iterator 
	 cit = _all_clusters.begin(); cit != _all_clusters.end(); ++cit ){
    if( ( cit->second.first == clusterPair.first && cit->second.second == clusterPair.second ) ||
	( cit->second.first == clusterPair.second && cit->second.second == clusterPair.first ) ){
      return cit->first;
    }
  }
  //branching not found create -- with sudakov and everything
  long thePartId;
  tcPDPtr particle_data;
  SudakovPtr theSudakov = getSud( thePartId, clusterPair.first, clusterPair.second ); 
  if( !theSudakov ){
    cerr << "can't find the sudakov in: \n"
	 << *clusterPair.first<<"\n"
	 << *clusterPair.second<<"\n";
  }
  Lorentz5Momentum pairMomentum = clusterPair.first->momentum() + 
    clusterPair.second->momentum();
  pairMomentum.setMass( 0. * MeV );
  particle_data = getParticleData( thePartId );
  
  //creates emitter particle
  ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, true ) );
  clustered->set5Momentum( pairMomentum );
  
  HardBranchingPtr clusteredBranch( new_ptr( HardBranching( clustered, theSudakov,
							    HardBranchingPtr(),
							    HardBranching::Outgoing ) ) );
  _all_clusters.insert( make_pair( clusteredBranch, clusterPair ) );

  //join children
  clusteredBranch->addChild( theParticles.find( clusterPair.first )->second );	    
  clusteredBranch->addChild( theParticles.find( clusterPair.second )->second );
  
  return clusteredBranch;
}

bool PowhegHandler::repeatProtoTree( ProtoTreePtr currentProtoTree ){
  //loop over all prototrees and see how many hardbranchings of curentProtoTree are found in each
  for( set< ProtoTreePtr >::const_iterator cit = _proto_trees.begin();
       cit != _proto_trees.end(); ++cit ){
    unsigned int no_matches = 0;
    for( set< HardBranchingPtr >::const_iterator ckt 
	   = currentProtoTree->getBranchings().begin(); ckt != currentProtoTree->getBranchings().end(); ckt++ ){
      if( (*cit)->getBranchings().find( *ckt ) != (*cit)->getBranchings().end() )
	no_matches++;
    }
    if( no_matches == currentProtoTree->getBranchings().size() ) return true;
  }
  return false;
}

bool PowhegHandler::simpleColConnections( ProtoTreePtr theProtoTree ){
  set< HardBranchingPtr > currentProtoTree = theProtoTree->getBranchings();
  if( currentProtoTree.size() != 2 ) {
    cerr<<"\n\nwrong size of proto tree: " 
	<< currentProtoTree.size() <<"\n\n\n";
    return false;
  }
  //colourline to join up q and qbar
  ColinePtr newline = new_ptr( ColourLine() );
  for( set< HardBranchingPtr >::iterator it = currentProtoTree.begin();
       it != currentProtoTree.end(); ++it ){
    (*it)->branchingParticle()->resetColour();
    if( (*it)->branchingParticle()->dataPtr()->iColour() 
	== PDT::Colour3 )
      newline->addColoured( (*it)->branchingParticle() );
   
    else if( (*it)->branchingParticle()->dataPtr()->iColour() 
	     == PDT::Colour3bar )
      newline->addAntiColoured( (*it)->branchingParticle() );
   
    else cerr<< "\n\n\nClustered back to gluon\n\n";
  }
  return true;
}

bool PowhegHandler::simpleColConnections( HardTreePtr theHardTree ){
  set< HardBranchingPtr > particles = theHardTree->branchings();
  //colourline to join up q and qbar
  ColinePtr newline = new_ptr( ColourLine() );
  for( set< HardBranchingPtr >::iterator it = particles.begin();
       it != particles.end(); ++it ){
    if( (*it)->status()==HardBranching::Incoming ) continue;
    (*it)->branchingParticle()->resetColour();
    if( (*it)->branchingParticle()->dataPtr()->iColour() 
	== PDT::Colour3 )
      newline->addColoured( (*it)->branchingParticle() );
    
    else if( (*it)->branchingParticle()->dataPtr()->iColour() 
	     == PDT::Colour3bar )
      newline->addAntiColoured( (*it)->branchingParticle() );
    else cerr<< "\n\n\nClustered back to gluon\n\n";
  }
  return true;
}

HardTreePtr PowhegHandler::doClusteringOrdered() {
  if(!_lepton) {
    return generalClustering();
  }
  ParticleVector theParts  = lastXCombPtr()->subProcess()->outgoing();
  
  //initialise global variables
  _all_clusters.clear();
  _proto_trees.clear();
  _hardTrees.clear();
  
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

  map <ShowerParticlePtr,HardBranchingPtr> theParticles;
  tcPDPtr particle_data;
  ShowerParticlePtr vBoson = new_ptr( ShowerParticle( *vb, 1, false, false ) );

  //is this the highest multiplcity channel
  if( theParts.size() == _max_mult ) {
    // cerr<<"highest mult channel \n";
    _highestMult = true;
  }
  else _highestMult = false;

  //loops through the FS particles and create hardBranchings
  for( unsigned int i = 0; i < theParts.size(); i++){
    ShowerParticlePtr currentParticle = 
      new_ptr( ShowerParticle( *theParts[i], 1, true, false ) );
    HardBranchingPtr currentBranching = 
      new_ptr( HardBranching( currentParticle, SudakovPtr(),
			      HardBranchingPtr(), HardBranching::Outgoing ) );     
    theParticles.insert( make_pair( currentParticle, currentBranching ) );
  }
  //create and initialise the first tree
  ProtoTreePtr initialProtoTree = new_ptr( ProtoTree() );
  for( map<ShowerParticlePtr, HardBranchingPtr>::iterator ita = theParticles.begin();
       ita != theParticles.end(); ita++ ){
    initialProtoTree->addBranching( ita->second );
  }
  _proto_trees.insert( initialProtoTree );

  fillProtoTrees( theParticles, initialProtoTree );
 
  double totalWeight = 0.;
  
  //create a hardtree from each proto tree and fill _hardTrees with angular ordered configs
  for( set< ProtoTreePtr >::const_iterator cit = _proto_trees.begin(); 
       cit != _proto_trees.end(); ++cit ){
    simpleColConnections( *cit );
    //vector boson branching
    vector<HardBranchingPtr> spaceBranchings;
    //all branchings
    vector<HardBranchingPtr> theBranchings;
    //fill theBranchings
    for( set< HardBranchingPtr >::const_iterator cjt = (*cit)->getBranchings().begin(); 
	 cjt != (*cit)->getBranchings().end(); ++cjt )
      theBranchings.push_back( *cjt );
    spaceBranchings.push_back( new_ptr( HardBranching( vBoson, SudakovPtr(),
						       HardBranchingPtr(), 
						       HardBranching::Incoming ) ) );
    theBranchings.push_back( spaceBranchings.back() );
    HardTreePtr powhegtree = new_ptr( HardTree( theBranchings,
						spaceBranchings,
						ShowerInteraction::QCD) );  
    // Calculate the shower variables
    // if momentum deconstruction fails then continue and ignore
    if( ! evolver()->showerModel()->kinematicsReconstructor()
	->deconstructDecayJets( powhegtree, evolver() , ShowerInteraction::QCD )) continue;
    
    //only insert angular ordered hardTrees
    if( powhegtree->checkHardOrdering() ) {
      //find the wgt and fill _hardTrees map
      powhegtree->findNodes();
      double treeWeight = sudakovWeight( powhegtree );
      //  treeWeight *= splittingFnWeight( powhegtree );
      _hardTrees.push_back( make_pair( powhegtree, treeWeight ) );
      totalWeight += treeWeight;
    }
  }

  if( _hardTrees.empty() )
    return HardTreePtr();

  //the hardTreePtr that is to be returned
  HardTreePtr chosen_hardTree;
  
  //choose a hardTree from shower probability
  if( _clusterOption == 0 ){
    long treeIndex;
    do{
   
    } while ( _hardTrees[ treeIndex ].second / totalWeight < UseRandom::rnd() );
    chosen_hardTree = _hardTrees[ treeIndex ].first;
  }

  //choose hardtree with lowest pt
  else if( _clusterOption == 2 ){
    //set min pt to be large
    Energy min_pt = 9999999999.*GeV;
    for(unsigned int ix = 0; ix < _hardTrees.size(); ix++ ){
      if( _hardTrees[ix].first->totalPt() < min_pt ){
	min_pt = _hardTrees[ix].first->totalPt();
	chosen_hardTree = _hardTrees[ix].first;
      }
    }
  }
  
  //choose hardtree with highest prob (cluster option = 3)
  else{
    //set min pt to be large
    double max_prob = 0.;
    for(unsigned int ix = 0; ix < _hardTrees.size(); ix++ ){
      if( _hardTrees[ ix ].second > max_prob ){
	max_prob = _hardTrees[ ix ].second;
	chosen_hardTree = _hardTrees[ix].first;
      }
      if( _hardTrees[ ix ].second == max_prob && UseRandom::rndbool() ){
	max_prob = _hardTrees[ ix ].second;
	chosen_hardTree = _hardTrees[ix].first;
      }
    }
  }
  //re-do momentum deconstruction (has been overridden by other trees otherwise)
  simpleColConnections( chosen_hardTree );
  if( ! evolver()->showerModel()->kinematicsReconstructor()
      ->deconstructDecayJets( chosen_hardTree, evolver(), ShowerInteraction::QCD ) )
    cerr<<"\n\nproblem doing momentum decon in selected tree \n\n";
  if( ! chosen_hardTree ) {
    cerr<<"PowhegHandler::problem in choosing hard tree\n";
    return HardTreePtr();
  }

  return chosen_hardTree;
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
  // and remove ones which aren't allowed
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
      deconstructDecayJets(newTree,evolver(),ShowerInteraction::QCD);
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
					      particle->isFinalState() ?
					      HardBranching::Outgoing : 
					      HardBranching::Incoming));
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
  // incoming lines and spacelike into the hard process
  for(it=incoming.begin();it!=incoming.end();++it) {
    spacelike.push_back((**it).convert());
    HardBranchingPtr br(spacelike.back());
    while (!br->children().empty()) {
      for(unsigned int ix=0;ix<br->children().size();++ix) {
	if(br->children()[ix]->status()==HardBranching::Incoming) {
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
  HardTreePtr newTree = new_ptr(HardTree(branchings,spacelike,
					 ShowerInteraction::QCD));
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
    if((**it).status()==HardBranching::Incoming) {
      meMomenta.push_back((**it).branchingParticle()->momentum());
      mePartonData.push_back((**it).branchingParticle()->dataPtr());
      particles.push_back((**it).branchingParticle());
    }
  }
  for(it=tree->branchings().begin();it!=tree->branchings().end();++it) {
    if(!(**it).status()==HardBranching::Incoming) {
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

double PowhegHandler::Sud( Energy scale, long id, Energy pt_cut ){
  //upper limit on scale 
  double sudwgt = 1.;
  Energy scale_cut = _max_qtilde;
  multimap< long, pair< Interpolator2d< double, Energy, Energy >::Ptr, Energy > >::const_iterator cjt;
  for( cjt =  _fbranchings.lower_bound( abs( id ) );
       cjt != _fbranchings.upper_bound( abs( id ) );
       ++cjt ) {
    //check to see if we are within qtilde limits before calling interpolator
   
    //pt_cut called here is in the durham or luclus jet measure
    
    if( scale < scale_cut && scale > cjt->second.second ) 
      sudwgt *= (* cjt->second.first )( pt_cut, scale );
    else if ( scale < cjt->second.second )
      sudwgt *= (* cjt->second.first )( pt_cut, cjt->second.second);
    else 
      sudwgt *= (* cjt->second.first )( pt_cut, scale_cut);
    
    //check for errors in sud value
    if( sudwgt == 0 ){ cerr<<"zero sud at scale = "<< scale / GeV
			   <<", id = "<< id 
			   <<", pt_cut = " << pt_cut / GeV 
			   <<", scale_cut = " << scale_cut / GeV
			   <<"\n";
    }
    if( sudwgt > 1.2 ){ cerr<<"large sud = "<< sudwgt<<" at scale = "<< scale / GeV
			   <<", id = "<< id 
			   <<", pt_cut = " << pt_cut / GeV 
			   <<", scale_cut = " << scale_cut / GeV
			   <<"\n";
    }
    if( isnan( sudwgt ) ){ cerr<<"nan sud at scale = "<< scale / GeV
			      <<", id = "<< id 
			       << ", pt_cut = " << pt_cut / GeV <<"\n";
    }
    if( isinf( sudwgt ) ){ cerr<<"inf sud at scale = "<< scale / GeV
			       <<", id = "<< id 
			       << ", pt_cut = " << pt_cut / GeV <<"\n";
    }
  }
  return sudwgt;
}

double PowhegHandler::splittingFnWeight( HardTreePtr theTree ){
  double splitFnWgt = 1.;
  for( map< HardBranchingPtr, Energy>::const_iterator cit = theTree->getNodes().begin();
       cit != theTree->getNodes().end(); ++cit ) {
    vector< long > ids;
    ids.push_back( cit->first->branchingParticle()->id() );
    if( ! cit->first->children().empty() ) {
      ids.push_back( cit->first->children()[0]->branchingParticle()->id() );
      ids.push_back( cit->first->children()[1]->branchingParticle()->id() );
    }
    else cerr<< "splittingFnWeight(): node with no children found\n";
    double z = cit->first->children()[0]->z();
    Energy2 t = z * ( 1. - z ) * sqr( cit->second );
    splitFnWgt *= cit->first->sudakov()->splittingFn()->P( z, t, ids, true );
  }
  return splitFnWgt;
}

double PowhegHandler::sudakovWeight( HardTreePtr theTree ) {
  double SudWgt = 1.;
  //ktcut for sudakovs
  Energy kt_cut;
  if( ! _highestMult ) kt_cut = sqrt( _yini * _s );
  else {
    //does this lowest pt method return something in the luc/dur jet measure?
    kt_cut = theTree->lowestPt( _jetMeasureMode, _s );
    if( kt_cut > _max_pt_cut )
      kt_cut = _max_pt_cut;
  }

  //external line weight
  for( map< ShowerParticlePtr, HardBranchingPtr >::const_iterator cit = 
	 theTree->getExternals().begin();
       cit != theTree->getExternals().end(); ++cit ) {
    
    Energy scale= sqrt( _s );
    if( cit->second ){
      if( cit->first == cit->second->children()[0]->branchingParticle() )
	scale =  cit->second->scale() * cit->second->children()[0]->z();
      else if(  cit->first == cit->second->children()[1]->branchingParticle() )
	scale =  cit->second->scale() * cit->second->children()[1]->z();
      else cerr<<"could not find child in external HardBranching \n";
    }
    else{
      scale = cit->first->evolutionScale();
    }
    SudWgt *= Sud( scale, cit->first->id(), kt_cut );
   
  }
  if(SudWgt > 1.1) cerr<<"\n\n\nsudakov from externals > 1!!\n\n\n";
  //intermediate line wgts
  for( map< long, pair< Energy, Energy > >::const_iterator cit 
	 = theTree->getInternals().begin();
       cit != theTree->getInternals().end(); ++cit ) {
    
    Energy scale =  cit->second.first;
    
    double internal_wgt = Sud( scale, cit->first, kt_cut );
    scale =  cit->second.second;

    internal_wgt /= Sud( scale, cit->first, kt_cut );
    if(internal_wgt > 1.1 || internal_wgt < 0.)cerr<<"\n\nbig internal weight of "<< internal_wgt
						   <<"\nnum scale = "
						   <<cit->second.first / GeV
						   <<"\nden scale = "
						   <<cit->second.second /GeV
						   <<"\n\n";
    SudWgt *= internal_wgt;
  }
  if(SudWgt > 1.1 ) cerr<<"sud wgt is "<<SudWgt<<"\n";
 
  double alphaWgt = 1.;
  //alphaS weight
  for( map< HardBranchingPtr, Energy >::const_iterator cit = theTree->getNodes().begin(); 
       cit != theTree->getNodes().end(); ++cit ) {
    if( ! cit->first->children().empty() ){
      alphaWgt *= _alphaS->value( sqr( cit->first->children()[0]->pT() ) ) / _alphaSMG;
    }
      else cerr << "sudakovWeight(): node with no children \n"; 
  }
 
  if( SudWgt > 1.1 ) {
    cerr<<"\n\nweight exceeded 1 in PowhegHandler::reweight() !!! \n";
    cerr<<"  alpha wgt = "<< alphaWgt
   	<<"\n  sudWgt = "<< SudWgt<<"\n\n";
  }
  if( _reweightOpt == 0 )
    return SudWgt*alphaWgt;
  else
    return SudWgt;
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
      if(br->children()[ix]->status()==HardBranching::Incoming) {
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
      if(br->children()[ix]->status()==HardBranching::Incoming) {
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
    if((**it).status()==HardBranching::Incoming) continue;
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

void PowhegHandler::testSudakovs(){
  ofstream sudTestOut;
  Energy deltaQt = 0.5 * GeV;
  //vector of pts to evaluate at and the colours they should be on plot
  vector< pair<Energy, string> > thePts;
  thePts.push_back( make_pair( 0.*GeV, string("BLACK" ) ) );
  thePts.push_back( make_pair( sqrt(_s*0.001), string("RED" ) ) );
  thePts.push_back( make_pair( sqrt(_s*0.005 ), string("BLUE" ) ) );
  thePts.push_back( make_pair( sqrt(_s*0.01 ), string("GREEN" ) ) );
  thePts.push_back( make_pair( sqrt(_s*0.05 ), string("CYAN" ) ) );

  sudTestOut.open( "sudTest.top" );


  multimap< long, pair< Interpolator2d< double, Energy, Energy >::Ptr, Energy > >::const_iterator cjt;
  for( cjt =  _fbranchings.begin();
       cjt != _fbranchings.end();
       ++cjt ) {
 
    //loop over pts
    sudTestOut <<"NEW FRAME \nSET WINDOW X 1.6 8 Y 3.5 9\nSET FONT DUPLEX\n"
	       <<"TITLE TOP \"Sud Test "<< cjt->first
	       <<": BLACK:y_ms=0, RED:y_ms=0.001, BLUE:y_ms=0.005, GREEN:y_ms=0.01, CYAN:y_ms=0.05\" \n"
	       <<"CASE      \"\" \nTITLE LEFT \"Sud(qtilde)\" \nCASE \" \" \n"
	       <<"SET ORDER X Y DX \nTITLE BOTTOM \"qtilde / GeV \" \n";
    for( vector< pair< Energy, string > >::const_iterator cit = thePts.begin();
	 cit != thePts.end(); ++cit ){ 
      Energy qtilde = cjt->second.second;
      while( qtilde < _max_qtilde ){
	double sud_val = (*cjt->second.first)( cit->first, qtilde );
	sudTestOut << qtilde / GeV <<"\t"<< sud_val <<"\t"<< deltaQt / 2. / GeV<<"\n";
	qtilde += deltaQt;
      }
      sudTestOut <<"HIST "<< cit->second <<"\n";
    }
  }
  sudTestOut.close();
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

  //is this the highest multiplcity channel
  if( theParts.size() == _max_mult ) _highestMult = true;
  else _highestMult = false;

  map <ShowerParticlePtr,HardBranchingPtr> theParticles;
  tcPDPtr particle_data;
  ShowerParticlePtr vBoson = new_ptr( ShowerParticle( *vb, 1, false, false ) );
  //loops through the FS particles and create naon branchings
  for( unsigned int i = 0; i < theParts.size(); i++){
    ShowerParticlePtr currentParticle = 
      new_ptr( ShowerParticle( *theParts[i], 1, true, false ) );
    theParticles.insert(make_pair(currentParticle, 
				  new_ptr( HardBranching( currentParticle, SudakovPtr(),
							  HardBranchingPtr(),
							  HardBranching::Outgoing ) ) ) );
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

    //get number of qqbar pairs
    int qq_pairs = 0;
    for( map<ShowerParticlePtr, HardBranchingPtr>::iterator ita = theParticles.begin();
	 ita != theParticles.end() ; ita++ ) {
      if( ita->first->id() > 0 && ita->first->id() < 7 ) qq_pairs++;
    }

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
    SudakovPtr theSudakov = getSud( thePartId,
				    clusterPair.first, clusterPair.second ); 
    if( !theSudakov ){
      cerr << "can't find the sudakov in: \n";
      cerr << *clusterPair.first<<"\n"
	   << *clusterPair.second<<"\n";
    }
    Lorentz5Momentum pairMomentum = clusterPair.first->momentum() + 
      clusterPair.second->momentum();
    pairMomentum.setMass(0.*MeV);
    particle_data = getParticleData( thePartId );
    
    //creates emitter particle
    ShowerParticlePtr clustered = new_ptr( ShowerParticle( particle_data, true ) );
    clustered->set5Momentum( pairMomentum );
  
    HardBranchingPtr clusteredBranch( new_ptr( HardBranching( clustered, theSudakov,
							      HardBranchingPtr(), 
							      HardBranching::Outgoing) ) );
    fixColours( clustered, clusterPair.first, clusterPair.second );
    theParticles.insert( make_pair( clustered, clusteredBranch ) );

    //add children
    clusteredBranch->addChild(theParticles.find(clusterPair.first )->second);
    clusteredBranch->addChild(theParticles.find(clusterPair.second)->second);
    theParticles.erase( clusterPair.first  );
    theParticles.erase( clusterPair.second );
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
						     HardBranching::Incoming ) ) );
  theBranchings.push_back( spaceBranchings.back() );
  HardTreePtr powhegtree = new_ptr( HardTree( theBranchings,
					      spaceBranchings,
					      ShowerInteraction::QCD) );

  // Calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->
    deconstructDecayJets( powhegtree, evolver(),ShowerInteraction::QCD );

  //keep track of proportion of trees that are ordered
  _trees_created++;
  if( powhegtree->checkHardOrdering() )
    _ordered_trees_created++;
  
  powhegtree->findNodes();

  return powhegtree;
}

void PowhegHandler::getDalitz(){
  if( _theHardTree->getExternals().size() == 3 ){
    //   cerr<<"\n\ndoing dalitz:\n";
    Energy total_energy = 0. * GeV;
    double x1_dal=0.;
    double x2_dal=0.;
    long parent_id = 0;
    for( map< ShowerParticlePtr, HardBranchingPtr >::const_iterator cit = 
	   _theHardTree->getExternals().begin();
	 cit != _theHardTree->getExternals().end(); ++cit ) {
      total_energy += cit->first->momentum().t();
      
      if( abs( cit->first->id() ) < 7 ){
	if( cit->first->id() > 0 ) 
	  x1_dal = 2. * cit->first->momentum().t() / GeV;
	else  x2_dal = 2. * cit->first->momentum().t() / GeV;
      } 
    }
    x1_dal /= total_energy / GeV;
    x2_dal /= total_energy / GeV;
    
    for( set< HardBranchingPtr >::const_iterator cit = 
	   _theHardTree->branchings().begin();
	 cit != _theHardTree->branchings().end(); ++cit ) {
      if( (*cit)->status()==HardBranching::Incoming ) continue;
      if( (*cit)->branchingParticle()->id() > 0 ){ 
	if( (*cit)->children().size() == 2 && (*cit)->children()[0] 
	    && (*cit)->children()[1] ) {
	  parent_id = (*cit)->branchingParticle()->id();
	}
      }
      else {
	
	if( (*cit)->children().size() == 2 && (*cit)->children()[0] 
	    && (*cit)->children()[0] ) {
	  parent_id = (*cit)->branchingParticle()->id();
	}
      }
    }
    if( parent_id > 0 ) 
      _dalitz_from_q1.push_back( make_pair(x1_dal, x2_dal) );
    else if( parent_id < 0 ) 
      _dalitz_from_q2.push_back( make_pair(x1_dal, x2_dal) );
  }
}



void PowhegHandler::makeQtildeDist(){
  //currently this only makes plots for up quarks

  ofstream sudTestOut;
  Energy deltaQt = 0.5 * GeV;
  //vector of pts to evaluate at and the colours they should be on plot
  vector< Energy > thePts;
  thePts.push_back( sqrt( 0.05 * _s ) );
  thePts.push_back( sqrt( 0.01 * _s ) );
  thePts.push_back( sqrt( 0.005 * _s ) );
  thePts.push_back( sqrt( 0.001 * _s ) );

  //iniialise the histograms
  vector< HistogramPtr > theHists;  
  for(int ix = 0; ix < thePts.size(); ix ++ )
    theHists.push_back( new_ptr(Histogram(0.,100.,100) ) );

  sudTestOut.open( "QtildeDist.top" );
  
  //for each pt value a vector of qtilde and sud_bar(qtilde) can be obtained from 2d interpolator
  //from this can calculate sud(qtilde) and make a 1d interpolator going both ways from the vectors
  
  //iterator to the 2d interpolator of the up quark
  multimap< long, pair< Interpolator2d< double, Energy, Energy >::Ptr, Energy > >::const_iterator cjt = 
    _fbranchings.find( long( 1 ) );

  //fill the vector of qtilde values once and for all
  //how to know exactly what the minimum value is - this came from integrator
  vector< Energy > scale;
  Energy qtildemin = cjt->second.second;
  Energy qtildemax = sqrt( _s );
  for( unsigned int ix = 0; ix < _npoint; ix++ )
    scale.push_back( qtildemin + double( ix ) * ( qtildemax - qtildemin ) / double( _npoint - 1 ) );

  //create a vector of pairs of 1d interpolators
  vector < pair< Interpolator< double, Energy >::Ptr, Interpolator< Energy, double >::Ptr >  > sud_interp;

  for( unsigned int jx = 0; jx < thePts.size(); jx++ ){
    double sud_min =  (* cjt->second.first )( thePts[jx], qtildemax );

    vector<double> sud;
    for( unsigned int ix = 0; ix < scale.size(); ix ++)
      sud.push_back( sud_min /  (* cjt->second.first )( thePts[jx], scale[ix] ) );
    // construct the Interpolators                                                                                                                         
    Interpolator<double,Energy>::Ptr
      intq = new_ptr(Interpolator<double,Energy>(sud,scale,3));
    Interpolator<Energy,double>::Ptr
      ints = new_ptr(Interpolator<Energy,double>(scale,sud,3));
    sud_interp.push_back(  make_pair( intq, ints ) );
    //make the interpolators both ways and add to the vectors
  }
  //fill histograms with the qtilde distributions
  for( unsigned int ix = 0; ix < sud_interp.size(); ix ++ ){
    for(unsigned int jx = 0; jx < 100000000; jx ++ ){
      double sud_wgt = UseRandom::rnd();
      Energy solved_scale = (* sud_interp[ix].second )( sud_wgt );
      (* theHists[ix] ) += solved_scale / GeV;
    }
  }
  //output the histograms
  for( unsigned int ix = 0; ix < theHists.size(); ix ++ ){
    using namespace HistogramOptions;
    theHists[ix]->topdrawOutput(sudTestOut,Frame,
			    "RED",
			    "qtilde distribution",
			    "",
			    "",
			    "",
			    "qtilde / GeV",
			    "");
  }
  sudTestOut.close();
}
