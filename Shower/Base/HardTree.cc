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
  : _interaction(type),
    _branchings(branchings.begin(),branchings.end()),
    _spacelike (spacelike .begin(),spacelike .end())
{}

bool HardTree::findNodes() {
  //clear all containers to be filled
  _theExternals.clear();
  _theNodes.clear();
  _theInternals.clear();

  _lowestPt = HardBranchingPtr();

  //fix forward relations based on backChild relations
  //(these are correctly set during clustering)
  fixFwdBranchings();
  //find all branchings that initiate a timelike shower
  //also add nodes, externals and intermediates from the spacelike line
  set< HardBranchingPtr > FS_initiators;
  for( set<HardBranchingPtr>::const_iterator it = this->branchings().begin();
       it != this->branchings().end(); ++it)  {
    if( ! (*it)->branchingParticle()->coloured() ) continue;
    if( (*it)->status() ==HardBranching::Outgoing ){
      //remove any parent ptr that might be set
      (*it)->parent( HardBranchingPtr() );
      FS_initiators.insert( *it );
      continue;
    }
    HardBranchingPtr spacelike = *it;
    while( spacelike->parent() ) {
      assert( spacelike == spacelike->parent()->children()[0] );
      spacelike = spacelike->parent();
      if( ! _lowestPt || _lowestPt->children()[0]->pT() 
	  > spacelike->children()[0]->pT() )
	_lowestPt = spacelike;
      _theIntermediates.insert( make_pair( spacelike, 
					   spacelike->children()[0] ) );
      //nb for bkwd branchings the child contains the splitting variables
      _theNodes.insert( make_pair( spacelike->children()[0], spacelike->children()[0]->scale() ) );
      FS_initiators.insert( spacelike->children()[0] );
    }
    _theExternals.insert( make_pair( (*it)->branchingParticle(), *it ) );
  }
  //fix timelike parent relations once and for all
  //fill timelike shower externals nodes and intermediates
  for( set<HardBranchingPtr>::const_iterator it = FS_initiators.begin();
       it != FS_initiators.end(); ++it){
    fixParents( *it );
    fillNodes( *it );
  }
  return true;
}

bool HardTree::fillNodes( HardBranchingPtr branch ){
  if( ! branch->children().empty() ) {
    if( ! _lowestPt || _lowestPt->children()[0]->pT() > branch->children()[0]->pT() )
      _lowestPt = branch;
    
    _theIntermediates.insert( make_pair( branch, branch->children()[0] ) );
    _theNodes.insert( make_pair( branch, branch->scale() ) );
    fillNodes( branch->children()[0] );
    fillNodes( branch->children()[1] ); 
  }
  else 
    _theExternals.insert( make_pair( branch->branchingParticle(), branch ) );
  return true;
}

bool HardTree::connect(ShowerTreePtr shower) {
  _particles.clear();
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
  if( particles().size() == progenitors.size() ) return true;
  else{
    cerr<<"hardTree connect:: size of particles and progenitors does not match \n";
    return false;
  }
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

void HardTree::fillHardScales( HardBranchingPtr branch, vector< pair< Energy, double > > & currentLine ){
  if( branch->children().empty() ) return;
  else{
    //child[0] continues currentline child[1] creates a new line
    //copy contents of old line into newline
    vector< pair< Energy, double > > newHardLine = currentLine;
    currentLine.push_back( make_pair( branch->scale(),
				      branch->children()[0]->z() ) );
    fillHardScales( branch->children()[0], currentLine );
    newHardLine.push_back( make_pair( branch->scale(),
				      branch->children()[1]->z() ) );
    fillHardScales( branch->children()[1], newHardLine );
    _hard_line_scales.push_back( newHardLine );
  }
}
//this fn needs to be adapted for IS lines
bool HardTree::checkHardOrdering( ) {
  //this function also caculates sum of pts of all branchings
  _total_pt = 0. * GeV;
  _hard_line_scales.clear();
  //create timelike proto lines from the outgoing and hardbranchings (the ones in hard process)
  vector< pair< HardBranchingPtr, vector< pair< Energy, double > > > >  proto_lines;
  for( set<HardBranchingPtr>::const_iterator it = this->branchings().begin();
       it != this->branchings().end(); ++it)  {
    if( ! (*it)->branchingParticle()->coloured() ) continue;
    if( ! (*it)->status()==HardBranching::Incoming && ! (*it)->children().empty() ) { 
      vector< pair< Energy, double > > new_hard_line1;
      new_hard_line1.push_back( make_pair( (*it)->scale(), (*it)->children()[0]->z() ) );
      vector< pair< Energy, double > > new_hard_line2;
      new_hard_line2.push_back( make_pair( (*it)->scale(), (*it)->children()[1]->z() ) );
      proto_lines.push_back( make_pair( (*it)->children()[0], new_hard_line1 ) ); 
      proto_lines.push_back( make_pair( (*it)->children()[1], new_hard_line2 ) ); 
      //pts of children are equal so just add once 
      Energy branchingPt = (*it)->children()[0]->pT();
      //included factor to favour ISFS branchings between partons in same z direction
      if( (*it)->status() == HardBranching::Incoming && 
	  ( ( (*it)->branchingParticle()->momentum().z() > ZERO )
	    == ( (*it)->children()[1]->branchingParticle()->momentum().z() > ZERO ) ) )
	branchingPt *= 0.9;  
      _total_pt += branchingPt;
    }
    else if( (*it)->parent() ) {
      //trace all spacelike branchings back following parents
      HardBranchingPtr spacelike = *it;
      vector< pair< Energy, double > > space_like_line;
      while( spacelike->parent() ) {
	_total_pt += spacelike->pT();
	//create a protoline from the time like child of parent
	vector< pair< Energy, double > > time_like_line = space_like_line;
	if( spacelike->parent()->children().empty() ) 
	  cerr<<"HardTree::checkHardOrdering() connection problem\n";
	//timelike child is always first child
	time_like_line.push_back( make_pair( spacelike->parent()->scale(), 
					     spacelike->parent()->children()[1]->z() ) );
	proto_lines.push_back( make_pair( spacelike->parent()->children()[1], time_like_line ) );
	//add the parent to the space_like_line
	space_like_line.push_back( make_pair( spacelike->parent()->scale(), 1. ) );
	spacelike = spacelike->parent();
      }
      _hard_line_scales.push_back( space_like_line );
    }
  }
  //recursively fill all timelike lines from proto_lines
  for( unsigned int ix = 0; ix < proto_lines.size(); ix++ ) {
    fillHardScales( proto_lines[ix].first, proto_lines[ix].second );
    _hard_line_scales.push_back( proto_lines[ix].second );
  }
  //go down each line (outwards from hard sub process) checking angular ordering condition
  for( unsigned int ix = 0; ix < _hard_line_scales.size(); ix++ ){
    for( unsigned int jx = 0; jx < _hard_line_scales[ix].size(); jx++ ){
      if( jx == 0 ) 
	continue;
      //angular ordering condition: z_1*q_1 > q2
      //this should also work for spacelike lines since for those z was set to 1
      if( _hard_line_scales[ix][jx].first  > _hard_line_scales[ix][jx - 1].first * _hard_line_scales[ix][ jx - 1 ].second )
	return false;
    }
  }
  return true;
}

Energy HardTree::lowestPt( int jetMeasureMode, Energy2 s ){
  
  //check to see we found the _lowest pt correctly
  if( !_lowestPt ) {
    cerr<<"null lowestpt from tree:\n";
    return 0.*GeV;
  }
  if(  _lowestPt->children().size() != 2 ) {
    cerr<< "HardTree: wrong no. children = " << _lowestPt->children().size() << " \n";
    return 0.*GeV;
  }
  if( ! _lowestPt->children()[0] ) {
    cerr<< "HardTree: null child \n";
    return 0.*GeV;
  }

  Energy kt_softest = _lowestPt->children()[0]->pT();
  if( jetMeasureMode != 1 ){
    Energy pt = kt_softest;
    double z = _lowestPt->children()[0]->z();
   
    Energy2 m0 = sqr( _lowestPt->branchingParticle()->nominalMass() );
    Energy2 m1 = sqr( _lowestPt->children()[0]->branchingParticle()->nominalMass() );
    Energy2 m2 = sqr( _lowestPt->children()[1]->branchingParticle()->nominalMass() );

    double lambda = sqrt( 1. - 4.*m0/s );
    double beta1 = 2.*( m1 - sqr(z)*m0 + sqr(pt) )
      / z / lambda / ( lambda + 1. ) / s;
    double beta2 = 2.*( m2 - sqr( 1. - z )*m0 + sqr(pt) )
      / ( 1. - z ) / lambda / ( lambda + 1. ) / s;

    Energy E1 = sqrt(s)/2.*( z + lambda*beta1 );
    Energy E2 = sqrt(s)/2.*( (1.-z) + lambda*beta2 );
    Energy Z1 = sqrt(s)/2.*lambda*( z - beta1 );
    Energy Z2 = sqrt(s)/2.*lambda*( (1.-z) - beta2 );

    double costheta = ( Z1*Z2 - sqr(pt) )
      / sqrt( sqr(Z1)+sqr(pt) ) / sqrt( sqr(Z2)+sqr(pt) );

    Energy2 kt_measure;
    if(  jetMeasureMode == 0 )
      kt_measure = 2.*min( sqr(E1), sqr(E2) )*( 1. - costheta );
    else if( jetMeasureMode == 2 )
      kt_measure = 2.*sqr(E1)*sqr(E2)/sqr(E1+E2)*( 1. - costheta );
    
    kt_softest = sqrt( kt_measure );
  }
  return kt_softest;
}

bool HardTree::fixFwdBranchings(){
  //loop over the spacelike hardbranchings in hard process
  set< HardBranchingPtr >::iterator it;
  for( it = _branchings.begin(); it != _branchings.end(); ++it ){
    HardBranchingPtr current = *it;
    if( current->status() ==HardBranching::Outgoing ) continue;
    if( !current->branchingParticle()->coloured() ) continue;
    current->clearChildren();
    current->sudakov( SudakovPtr() );
    if( current->backChildren().empty() ) continue;
    vector< HardBranchingPtr > backChildren = current->backChildren();
    while( ! backChildren.empty() ){
      if( backChildren.size() != 2 ) {
	cerr << "fixFwdBranchings:: wrong number of back childrem\n";
	continue;
      }
      if( !backChildren[0] || !backChildren[1] ) continue;
      //remove any exiting children
      backChildren[0]->clearChildren();
      backChildren[0]->addChild( current );
      backChildren[0]->addChild( backChildren[1] );
      current->parent( backChildren[0] );
      backChildren[1]->parent ( backChildren[0] );
      if( !current->backSudakov() ){
	cerr<<"fixFwdBranchings: problem finding backSudakov \n";
	continue;
      }
      backChildren[0]->sudakov( current->backSudakov() );
      //continue along incoming line (always the first child)
      current = backChildren[0];
      backChildren = backChildren[0]->backChildren();
    } 
  }
  return true;
}

void HardTree::removeBackChildren(){
  //this function will only work once the forward relations have been set - do this first
  fixFwdBranchings();
  for( set<HardBranchingPtr>::iterator it = _branchings.begin();
       it != _branchings.end(); ++it ){
    HardBranchingPtr currentSpacelike = *it;
    while( currentSpacelike->parent() ){
      currentSpacelike->clearBackChildren();
      currentSpacelike = currentSpacelike->parent();
    }
  }
}

bool HardTree::fixParents( HardBranchingPtr branch ){
  if( branch->children().empty() ) return true;
  branch->children()[0]->parent( branch );
  fixParents(  branch->children()[0] );
  branch->children()[1]->parent( branch );
  fixParents(  branch->children()[1] );
  return true;
}

//check that spacelike lines areordered with x_frac getting smaller towards hard subprocess
bool HardTree::checkXOrdering( ) {
  //  return true;

  double x_frac;
  for( set< HardBranchingPtr >::const_iterator cit = incoming().begin(); 
       cit != incoming().end(); ++cit ) {
    //find the pdf for this line
    HardBranchingPtr currentParticle = *cit;
    if( currentParticle->status() == HardBranching::Outgoing || 
	!currentParticle->branchingParticle()->coloured() ) continue;
    Lorentz5Momentum beamMom = currentParticle->beam()->momentum();
    Lorentz5Momentum otherBeam = beamMom;
    otherBeam.setZ( -otherBeam.z() );

    // x_frac =  currentParticle->branchingParticle()->momentum() 
    //  * otherBeam / ( beamMom*otherBeam );
    x_frac = currentParticle->x_frac();

    //trace from incoming to hard process checking x_{i+1} < x_{i}
    while( !currentParticle->children().empty() ){
      currentParticle = currentParticle->children()[0];
      double last_xfrac = x_frac;
      //is this the correct momentum
      x_frac = currentParticle->x_frac();
      //  x_frac =  currentParticle->branchingParticle()->momentum() 
      //	* otherBeam / ( beamMom*otherBeam );
      if( x_frac  > last_xfrac ) {
	return false;
      }
      //  if( currentParticle->branchingParticle()->momentum().t() < ZERO )
      //	return false;
    }
  } 
  return true;
}
