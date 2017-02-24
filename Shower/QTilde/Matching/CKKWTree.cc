// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWTree class.
//

#include "CKKWTree.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig/Utilities/Maths.h"

using namespace Herwig;

CKKWTree::CKKWTree(vector<HardBranchingPtr> branchings,
		   vector<HardBranchingPtr> spacelike,
		   ShowerInteraction type) : HardTree(branchings,spacelike,type),
						   lowestpTMomentum_(ZERO),
						   totalpT_(ZERO), ordered_(false)
{}

bool CKKWTree::checkXOrdering() {
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

//calculates the highest merging scale in the jet measure
//does this based on the hardtree momenta rather than the shower variables
//cutting events based on this is equivalent to the mad graph cuts
Energy CKKWTree::lowestPtMomentum( int jetMeasureMode, int cutOption ){
  Energy lowestpTMomentum_ = 9999999. * GeV;
  //vector of timelike initiators
  vector< HardBranchingPtr > FSsHat_hower_initiators;
  //add FS shower initiators from FS progenitors
  for( set<HardBranchingPtr>::const_iterator it = branchings().begin();
       it != branchings().end(); ++it)  {
    if( ! (*it)->branchingParticle()->coloured() || (*it)->status() == HardBranching::Incoming  ) continue;
    FSsHat_hower_initiators.push_back( *it );
  }
  //add FS shower initiators from IS branchings
  //trace back all incoming lines
  for( set< HardBranchingPtr >::const_iterator cit = incoming().begin(); 
       cit != incoming().end(); ++cit ) {
    if( ! (*cit)->branchingParticle()->coloured() ) continue;
    HardBranchingPtr currentParticle = *cit;
    while( !currentParticle->children().empty() ){
      //get scale for hardbranching - order of partons here IS first??
      //is this correct for ISFS branching should be between widtilde{ij} and j creating i
      //I think it was wrong before so have changed it -depended only on the second FS parton though!!
      Energy kt_measure = hadronJetMeasure( currentParticle->branchingParticle()->momentum(), 
					    currentParticle->children()[1]->branchingParticle()->momentum(),
					    false );
      //if cutOption == 2 then only set this if the only partons 
      //involved are endpoints (no parent for IS and no children for FS)
      if( !( cutOption == 2 && !externalBranching( currentParticle, currentParticle->children()[1] ) )
	  && kt_measure < lowestpTMomentum_ ) 
	lowestpTMomentum_ = kt_measure;				  
      FSsHat_hower_initiators.push_back(  currentParticle->children()[1] );
      currentParticle = currentParticle->children()[0];
    }
  }
  //loop over the time like initiators
  for( vector< HardBranchingPtr >::const_iterator cit = FSsHat_hower_initiators.begin(); 
       cit != FSsHat_hower_initiators.end(); ++cit ){
    getLowestJetMeasure( *cit, jetMeasureMode, cutOption );
  }
  return lowestpTMomentum_;
}

void CKKWTree::getLowestJetMeasure( HardBranchingPtr branch, int jetMeasureMode, int cutOption ){
  //if branching has children then find the jet measure from them
  if( ! branch->children().empty() ){
    Energy kt_measure = ZERO;
    if( jetMeasureMode == 0 || jetMeasureMode == 0 ){
      kt_measure = getJetMeasure( branch->children()[0]->branchingParticle()->momentum(),
				  branch->children()[1]->branchingParticle()->momentum(),
				  jetMeasureMode );
    }
    else if( jetMeasureMode == 3 ) {
      kt_measure = hadronJetMeasure( branch->children()[0]->branchingParticle()->momentum(),
				     branch->children()[1]->branchingParticle()->momentum(),
				     true );
    }
    if( !( cutOption == 2 && !externalBranching( branch->children()[0], branch->children()[1] ) )
	&& kt_measure < lowestpTMomentum_ ) lowestpTMomentum_ = kt_measure;
    //do the same for children
    getLowestJetMeasure( branch->children()[0], jetMeasureMode, cutOption );
    getLowestJetMeasure( branch->children()[1], jetMeasureMode, cutOption );
  }  
}

Energy CKKWTree::hadronJetMeasure( const Lorentz5Momentum & p1,
				   const Lorentz5Momentum & p2,
				   bool final ) {
  Energy kt_measure(ZERO);
  //FSFS case
  if( final ) {
    double deltay   = p1.rapidity() - p2.rapidity();
    double deltaphi = Herwig::Math::angleMinusPiToPi(p1.phi() - p2.phi());
    double deltaR = sqr( deltay ) + sqr( deltaphi );
    kt_measure = sqrt( min( p1.perp2(), p2.perp2() ) * deltaR );
  }
  //in the case of ISFS the merge scale is given by the pt 
  //of the FS parton w.r.t incoming hadron (assumed second argument)
  else {
    kt_measure = p2.perp();
  }
  return kt_measure;
}

Energy CKKWTree::getJetMeasure( const Lorentz5Momentum & p1,
				const Lorentz5Momentum & p2,
				int jetMeasureMode ) {
  Energy kt_measure(ZERO);
  double costheta = p1.vect().cosTheta( p2.vect() );
  switch( jetMeasureMode ){
  case 0:
    kt_measure = sqrt( 2. * sqr( max(p1.e(), p2.e() ) ) * ( 1. - costheta ) );
    break;
  case 2:
    kt_measure = sqrt( 2. * sqr( p1.e() * p2.e() / ( p1.e() + p2.e() ) )
		       * ( 1. - costheta ) );
    break;
  default:
    assert(false);
  }
  return kt_measure;
}

bool CKKWTree::externalBranching( HardBranchingPtr a, HardBranchingPtr b ){
  if( a->status() == HardBranching::Incoming &&  a->parent()           ) return false; 
  if( a->status() == HardBranching::Outgoing && !a->children().empty() ) return false; 
  if( b->status() == HardBranching::Incoming &&  b->parent()           ) return false; 
  if( b->status() == HardBranching::Outgoing && !b->children().empty() ) return false; 
  return true;
}

bool CKKWTree::checkHardOrdering( ) {
  //this function also caculates sum of pts of all branchings
  totalpT_ = 0. * GeV;
  hardLineScales_.clear();
  //create timelike proto lines from the outgoing and hardbranchings (the ones in hard process)
  vector< pair< HardBranchingPtr, vector< pair< Energy, double > > > >  proto_lines;
  for( set<HardBranchingPtr>::const_iterator it = this->branchings().begin();
       it != this->branchings().end(); ++it)  {
    if( ! (*it)->branchingParticle()->coloured() ) continue;
    if(   (*it)->status()!=HardBranching::Incoming && ! (*it)->children().empty() ) { 
      vector< pair< Energy, double > > new_hard_line1;
      new_hard_line1.push_back( make_pair( (*it)->scale(), (*it)->children()[0]->z() ) );
      vector< pair< Energy, double > > new_hard_line2;
      new_hard_line2.push_back( make_pair( (*it)->scale(), (*it)->children()[1]->z() ) );
      proto_lines.push_back( make_pair( (*it)->children()[0], new_hard_line1 ) ); 
      proto_lines.push_back( make_pair( (*it)->children()[1], new_hard_line2 ) ); 
      //pts of children are equal so just add once 
      Energy branchingPt = (*it)->children()[0]->pT();
      //included factor to favour ISFS branchings between partons in same z direction
      //hard coded factor here should make a parameter
      if( (*it)->status() == HardBranching::Incoming && 
	  ( ( (*it)->branchingParticle()->momentum().z() > ZERO )
	    == ( (*it)->children()[1]->branchingParticle()->momentum().z() > ZERO ) ) )
	branchingPt *= 0.9;  
      totalpT_ += branchingPt;
    }
    else if( (*it)->parent() ) {
      //trace all spacelike branchings back following parents
      HardBranchingPtr spacelike = *it;
      vector< pair< Energy, double > > space_like_line;
      while( spacelike->parent() ) {
	Energy branchingPt = spacelike->pT();
	assert( spacelike->parent()->children().size() == 2 );
	HardBranchingPtr emitted = spacelike->parent()->children()[0] == spacelike ?
	  spacelike->parent()->children()[1] : spacelike->parent()->children()[0];
	// included factor to favour ISFS branchings between partons in same z direction
	//hard coded factor here should make a parameter
	if( ( emitted            ->branchingParticle()->momentum().z() > ZERO ) ==
	    ( spacelike->parent()->branchingParticle()->momentum().z() > ZERO ) )
	  branchingPt *= 0.9;  
 	totalpT_ += branchingPt;
	//create a protoline from the time like child of parent
	vector< pair< Energy, double > > time_like_line = space_like_line;
	//timelike child is always first child
	time_like_line.push_back( make_pair( spacelike->parent()->scale(), 
					     spacelike->parent()->children()[1]->z() ) );
	proto_lines.push_back( make_pair( spacelike->parent()->children()[1], time_like_line ) );
	//add the parent to the space_like_line
	space_like_line.push_back( make_pair( spacelike->parent()->scale(), 1. ) );
	spacelike = spacelike->parent();
      }
      hardLineScales_.push_back( space_like_line );
    }
  }
  //recursively fill all timelike lines from proto_lines
  for( unsigned int ix = 0; ix < proto_lines.size(); ix++ ) {
    fillHardScales( proto_lines[ix].first, proto_lines[ix].second );
    hardLineScales_.push_back( proto_lines[ix].second );
  }
  //go down each line (outwards from hard sub process) checking angular ordering condition
  for( unsigned int ix = 0; ix < hardLineScales_.size(); ix++ ){
    for( unsigned int jx = 0; jx < hardLineScales_[ix].size(); jx++ ){
      if( jx == 0 ) 
	continue;
      //angular ordering condition: z_1*q_1 > q2
      //this should also work for spacelike lines since for those z was set to 1
      if( hardLineScales_[ix][jx].first  > hardLineScales_[ix][jx - 1].first * hardLineScales_[ix][ jx - 1 ].second )
	return false;
    }
  }
  return true;
}

void CKKWTree::fillHardScales( HardBranchingPtr branch, vector< pair< Energy, double > > & currentLine ){
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
    hardLineScales_.push_back( newHardLine );
  }
}

void CKKWTree::findNodes() {
  //clear all containers to be filled
  nodes_.clear();
  lowestpT_ = HardBranchingPtr();
  // find all branchings that initiate a timelike shower
  // also add nodes from the spacelike line
  for( set<HardBranchingPtr>::const_iterator it = this->branchings().begin();
       it != this->branchings().end(); ++it)  {
    if( ! (*it)->branchingParticle()->coloured() ) continue;
    if( (*it)->status() ==HardBranching::Outgoing ){
      //remove any parent ptr that might be set
      (*it)->parent( HardBranchingPtr() );
      fillNodes( *it );
      continue;
    }
    HardBranchingPtr spacelike = *it;
    while( spacelike->parent() ) {
      assert( spacelike == spacelike->parent()->children()[0] );
      spacelike = spacelike->parent();
      if( ! lowestpT_ || lowestpT_->children()[0]->pT() 
	  > spacelike->children()[0]->pT() )
	lowestpT_ = spacelike;
      //nb for bkwd branchings the child contains the splitting variables
      nodes_.insert( make_pair( spacelike->children()[0], spacelike->children()[0]->scale() ) );
      fillNodes( spacelike->children()[0] );
    }
  }
  ordered_ = checkHardOrdering() && checkXOrdering();
}

bool CKKWTree::fillNodes( HardBranchingPtr branch ){
  if( ! branch->children().empty() ) {
    if( ! lowestpT_ || lowestpT_->children()[0]->pT() > branch->children()[0]->pT() )
      lowestpT_ = branch;
    
    nodes_.insert( make_pair( branch, branch->scale() ) );
    fillNodes( branch->children()[0] );
    fillNodes( branch->children()[1] ); 
  }
  return true;
}

Energy CKKWTree::lowestPt( int jetMeasureMode, Energy2 s ){
  // perform various checks
  assert( lowestpT_ && lowestpT_->children().size() == 2 && 
	  lowestpT_->children()[0] );


  Energy ktsHat_oftest = lowestpT_->children()[0]->pT();
  if( jetMeasureMode == 0 || jetMeasureMode == 2 ){
    Energy pt = ktsHat_oftest;
    double z = lowestpT_->children()[0]->z();
   
    Energy2 m0 = sqr( lowestpT_->branchingParticle()->nominalMass() );
    Energy2 m1 = sqr( lowestpT_->children()[0]->branchingParticle()->nominalMass() );
    Energy2 m2 = sqr( lowestpT_->children()[1]->branchingParticle()->nominalMass() );

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

    Energy2 kt_measure(ZERO);
    if(  jetMeasureMode == 0 )
      kt_measure = 2.*min( sqr(E1), sqr(E2) )*( 1. - costheta );
    else if( jetMeasureMode == 2 )
      kt_measure = 2.*sqr(E1)*sqr(E2)/sqr(E1+E2)*( 1. - costheta );
    else
      assert(false);
    ktsHat_oftest = sqrt( kt_measure );
  }
  else if( jetMeasureMode == 3 && lowestpT_->status() != HardBranching::Incoming ){
    Energy2 m1 = sqr( lowestpT_->children()[0]->branchingParticle()->nominalMass() );
    Energy2 m2 = sqr( lowestpT_->children()[1]->branchingParticle()->nominalMass() );

    Energy pt = ktsHat_oftest;
    double z = lowestpT_->z();
    
    double beta1 = 2.*( m1 + sqr(pt) ) / z  / s;
    double beta2 = 2.*( m2 + sqr(pt) ) / ( 1. - z ) / s;
      
    //delta phi is always pi for first emission (qt_i = +-pt)
    double deltaR = sqr(  log( z / beta1 ) - log( (1-z) / beta2 ) ) / 4. 
      + sqr( Constants::pi );
    ktsHat_oftest = pt*sqrt( deltaR );
  }
  return ktsHat_oftest;
}
