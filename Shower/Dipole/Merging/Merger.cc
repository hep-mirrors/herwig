  // -*- C++ -*-
  //
  // Merger.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright ( C ) 2002-2007 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL , see COPYING for details.
  // Please respect the MCnet academic guidelines , see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined , non-templated member
  // functions of the Merger class.
  //

#include "Merger.h"
#include "Node.h"
#include "MergingFactory.h"
  // other includes when needed below.
using namespace Herwig;

Merger::Merger() : MergerBase() ,
  Unlopsweights( true ) , theCMWScheme( true ) ,
  isUnitarized( true ) , isNLOUnitarized( true ) ,
  defMERegionByJetAlg( false ) , theOpenInitialStateZ( false ) ,
  theChooseHistory( 0 ) , theN0( 0 ) , theOnlyN( -1 ) ,
  theCurrentMaxLegs( -1 )  , theGamma( 1. ) ,
  ee_ycut( -1 ) , pp_dcut( -1 ) , theSmearing( 0. ) ,
  theIRSafePT( 1.*GeV ) ,
  theMergePt( 4.*GeV ) , theCentralMergePt( 4.*GeV ) ,
  theMergingJetFinder() , theLargeNBasis() ,
  theDipoleShowerHandler() , theTreeFactory()
{
  FFLTK = new_ptr( FFLightTildeKinematics() );
  FILTK = new_ptr( FILightTildeKinematics() );
  IFLTK = new_ptr( IFLightTildeKinematics() );
  IILTK = new_ptr( IILightTildeKinematics() );
  FFMTK = new_ptr( FFMassiveTildeKinematics() );
  FIMTK = new_ptr( FIMassiveTildeKinematics() );
  IFMTK = new_ptr( IFMassiveTildeKinematics() );
}

Merger::~Merger() {}

IBPtr Merger::clone() const {
  return new_ptr( *this );
}

IBPtr Merger::fullclone() const {
  return new_ptr( *this );
}

pair<PVector , PVector> getInOut( NodePtr Node ){
  PVector incoming;
  const auto me=Node->nodeME();
  for( auto const & i : {0 , 1} )
  incoming.push_back(
         me->mePartonData()[i]->produceParticle(
         me->lastMEMomenta()[i] ) );
  PVector outgoing;
  for ( size_t i = 2;i<Node->nodeME()->mePartonData().size();i++ ){
    PPtr p  = me->mePartonData()[i]->produceParticle(
              me->lastMEMomenta()[i] );
    outgoing.push_back( p );
  }
  return make_pair( incoming , outgoing );
}

CrossSection Merger::MergingDSigDRBornGamma( NodePtr Node ){
  
  double weightCL = 0.;
  weight = 1.;
  
  if ( !Node->children().empty() ) {
    auto const inoutpair = getInOut( Node );
      // Here we make sure that clustering and splitting are invertible.
    NodePtr rndCh =  Node->randomChild();
      // Check if point is part of the ME region.
    if( !matrixElementRegion( inoutpair.first ,
                              inoutpair.second ,
                              rndCh->pT() ,
                              theMergePt ) )weight *= 0.;
  }
  
  NodePtr Born =  Node-> getHistory( true , DSH()->hardScaleFactor() );
  NodePtr CLNode;
  NodePtr BornCL;
  
  
  if( !Node->children().empty() ){
    if ( UseRandom::rndbool() ){
      CLNode =  Node->randomChild();
      bool inhist = CLNode->isInHistoryOf( Born );
      weight *= inhist?( -2.*int( Node->children().size() ) ):0.;
      projected = true;
      weightCL = 2.*int( Node->children().size() );
      BornCL = CLNode-> getHistory( false , DSH()->hardScaleFactor() );
    }else{
      weight = 2.;
      projected = false;
    }
  }else{
    weight = 1.;
    projected = false;
  }
  
  
  if ( treefactory()->onlymulti() != -1&&
      treefactory()->onlymulti() != int( Node->legsize()-(projected ? 1 : 0) ) )
  return ZERO;
  
  
  if( !Node->allAbove( mergePt()*(1.-1e-6)  ) )weight = 0.;
  if( CLNode&&!CLNode->allAbove( mergePt()*(1.-1e-6) ) )weightCL = 0.;
  if ( !Born->xcomb()->willPassCuts() ){
    return ZERO;
  }
  
  CrossSection res = ZERO;
  bool maxMulti = Node->legsize() == int( maxLegsLO() );
  
  
  if( weight != 0. ){
    Energy startscale = CKKW_StartScale( Born );
    fillHistory( startscale , Born , Node );
    Node->runningPt( fillProjector( (projected ? 1 : 0) ) );
    weight *= history.back().weight*alphaReweight()*pdfReweight();
    if( weight == 0.&&weightCL == 0. )return ZERO;
    
    res += weight*TreedSigDR( startscale , Node , ( !maxMulti&&!projected )?theGamma:1. );
  }
  
  if( CLNode&&theGamma != 1. ){
    Energy startscale = CKKW_StartScale( BornCL );
    fillHistory( startscale , BornCL , CLNode );
    Node->runningPt( fillProjector( projected ? 1 : 0 ) );
    weightCL *= history.back().weight*alphaReweight()*pdfReweight();
    CrossSection diff = ZERO;
    Node->nodeME()->factory()->setAlphaParameter( 1. );
    diff -=  weightCL*CLNode->dipole()->dSigHatDR( sqr( startscale*theCurrentME->renFac() ) );
    Node->nodeME()->factory()->setAlphaParameter( theGamma );
    
    string alp = ( CLNode->dipole()->aboveAlpha()?"Above":"Below" );
    
    diff += weightCL*CLNode->dipole()->dSigHatDR( sqr( startscale*theCurrentME->renFac() ) );
    Node->nodeME()->factory()->setAlphaParameter( 1. );
    
    res += diff;
  }
  
  
  return res;
}



double decideClustering(const NodePtr sub,const NodePtr head,bool& pro){
  if( sub !=  head ){
      // at least one history step -> unitarisation
    if ( UseRandom::rndbool() ){
      pro = true;
      return -2.;
    }else{
      pro = false;
      return 2.;
    }
  }else{
      // no ordered history -> no projection
    pro = false;
    return 1.;
  }
}



CrossSection Merger::MergingDSigDRBornStandard( NodePtr Node ){
      // Check if phase space poing is in ME region
  if ( !Node->children().empty() ) {
    auto const & inoutpair = getInOut( Node );
    NodePtr rndCh =  Node->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inoutpair.first , inoutpair.second , rndCh->pT() , theMergePt ) )return ZERO;
  }
  
    // get the history for the process
  NodePtr Born =  Node-> getHistory( true , DSH()->hardScaleFactor() );
  
    // decide if to cluster
  weight = decideClustering(Born,Node,projected);
  
  if ( treefactory()->onlymulti() != -1&&
      treefactory()->onlymulti() != int( Node->legsize()-( projected ? 1 : 0 ) ) )
      return ZERO;
    // Check for cuts on the production proces.
  if ( !Born->xcomb()->willPassCuts() )return ZERO;
    // calculate the staring scale
  Energy startscale = CKKW_StartScale( Born );
    // fill history with caluclation of sudakov supression
  fillHistory( startscale , Born , Node );
    // fill the projector -> return the scale of the last splitting
  Node->runningPt( fillProjector( projected ? 1 : 0 ) );
    // the weight has three components to get shower history weight
  weight *= history.back().weight*alphaReweight()*pdfReweight();
  if( weight == 0. )return ZERO;
    //calculate the cross section
  return weight*TreedSigDR( startscale , Node , 1. );
}


  /// calculate the history weighted born cross section
#include "ThePEG/Handlers/SamplerBase.h"
CrossSection Merger::MergingDSigDRBornCheapME( NodePtr Node ){
    // Check if phase space poing is in ME region
  if ( !Node->children().empty() ) {
    auto const & inoutpair = getInOut( Node );
    NodePtr rndCh =  Node->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inoutpair.first , inoutpair.second , rndCh->pT() , theMergePt ) )return ZERO;
  }
  
  CrossSection meweight=(SamplerBase::runLevel() == SamplerBase::RunMode)?TreedSigDR( 60.*GeV , Node , 1. ):1.*nanobarn;
  
  if (SamplerBase::runLevel() == SamplerBase::RunMode && theHighMeWeightMap.count(Node) ) {
    if (abs(meweight)>abs(theHighMeWeightMap[Node])) {
      theHighMeWeightMap[Node]=abs(meweight);
    }else{
      if (meweight/theHighMeWeightMap[Node]<UseRandom::rnd()) {
        return ZERO;
      }
    }
  }else{
    theHighMeWeightMap.insert(make_pair(Node,abs(meweight)));
  }
  
  
  
    // get the history for the process
  NodePtr Born =  Node-> getHistory( true , DSH()->hardScaleFactor() );
  
    // decide if to cluster
  weight = decideClustering(Born,Node,projected);
  
  if ( treefactory()->onlymulti() != -1&&
      treefactory()->onlymulti() != int( Node->legsize()-( projected ? 1 : 0 ) ) )
  return ZERO;
    // Check for cuts on the production proces.
  if ( !Born->xcomb()->willPassCuts() )return ZERO;
    // calculate the staring scale
  Energy startscale = CKKW_StartScale( Born );
  
    // fill history with caluclation of sudakov supression
  fillHistory( startscale , Born , Node );
    // fill the projector -> return the scale of the last splitting
  Node->runningPt( fillProjector( projected ? 1 : 0 ) );
    // the weight has three components to get shower history weight
  weight *= history.back().weight*alphaReweight()*pdfReweight();
  if( weight == 0. )return ZERO;
    //calculate the cross section
  return weight*TreedSigDR( startscale , Node , 1. )*theHighMeWeightMap[Node]/meweight;
}




CrossSection Merger::MergingDSigDRVirtualStandard( NodePtr Node ){
    // Check if phase space poing is in ME region
  if ( !Node->children().empty() ) {
    auto inoutpair = getInOut( Node );
    NodePtr rndCh =  Node->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inoutpair.first , inoutpair.second , rndCh->pT() , theMergePt ) )return ZERO;
  }
  
  NodePtr Born =  Node-> getHistory( true , DSH()->hardScaleFactor() );

  weight = decideClustering(Born,Node,projected);
  
  if ( !Born->xcomb()->willPassCuts() )return ZERO;

  Energy startscale = CKKW_StartScale( Born );
  fillHistory( startscale , Born , Node );
  Node->runningPt( fillProjector( projected ? 1 : 0 ) );
  
  double ww1 = history.back().weight;
  double ww2 = alphaReweight();
  double ww3 = pdfReweight();
  
  
  weight *= ww1*ww2*ww3;
  if( weight == 0. )return ZERO;
  
  CrossSection matrixElement = LoopdSigDR( startscale , Node );
  
  CrossSection Bornweight = Node->nodeME()->dSigHatDRB();
  
  double w1 = -sumpdfReweightUnlops();
  double w2 = -sumalphaReweightUnlops();
  double w3 = -sumfillHistoryUnlops();
  
  CrossSection unlopsweight  = ( w1+w2+w3 )
  *Bornweight
  *SM().alphaS()/( 2.*ThePEG::Constants::pi );
  
  if ( Node->legsize() == 5&&Debug::level > 2 ) {
    Energy minPT = Constants::MaxEnergy;
    for( auto const & no :Node->children() )minPT = min( no->pT() , minPT );
    
    generator()->log() << "\nVIRT " << minPT/GeV << " " << weight << " " << w1;
    generator()->log() << " " << w2;
    generator()->log() << " " << w3;
    generator()->log() << " " << ( matrixElement/nanobarn ) << " " << ww1 << " " << ww2 << " " << ww3 << " " << Born->pT()/GeV << " " << Born->nodeME()->mePartonData()[3]->mass()/GeV << " " << ( Bornweight*SM().alphaS()/( 2.*ThePEG::Constants::pi )/nanobarn );
  }
  
  
  return weight*
  as( startscale*DSH()->renFac() )/SM().alphaS()*
  ( matrixElement+unlopsweight );
}


CrossSection Merger::MergingDSigDRRealStandard( NodePtr Node ){
  bool allAbove = Node->allAbove( mergePt() );
    //TODO: VW Abgas Skandal
  if( !Node->allAbove( ( Debug::level > 2?0.01:1. )*theIRSafePT ) )return ZERO;
  if ( allAbove )return MergingDSigDRRealAllAbove( Node );
  if ( UseRandom::rndbool() )
  return 2.*MergingDSigDRRealBelowSubReal( Node );
  return 2.*MergingDSigDRRealBelowSubInt( Node );
}

CrossSection Merger::MergingDSigDRRealAllAbove( NodePtr Node ){
  if ( Node->children().empty() ) {
    throw Exception()
    << "Real emission contribution without underlying born."
    << "These are finite contibutions already handled in LO merging."
    << Exception::abortnow;
  }
  
  
    //If all dipoles pts are above , we cluster to this dipole.
  NodePtr CLNode =  Node->randomChild();
    // Check if phase space poing is in ME region--> else rm PSP
  if ( !CLNode->children().empty() ) {
    auto inoutpair = getInOut( CLNode );
    NodePtr rndCh =  CLNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inoutpair.first , inoutpair.second , rndCh->pT() , theMergePt ) )return ZERO;
  }
  
    // first find the history for the acctual Node
  NodePtr Born =  Node-> getHistory( true , DSH()->hardScaleFactor() );
    // check if the randomly choosen CLNode is part of the history.
    // If CLNode is not part of the history , dont calculate the Real contribution
    // else multiply the real contribution with N ( number of children ).
    // this returns the sudakov suppression according to the clustering of the born parts.
  bool inhist = CLNode->isInHistoryOf( Born );
    // get the history for the clustered Node.
  Born = CLNode-> getHistory( false , DSH()->hardScaleFactor() );

  weight = decideClustering(Born,CLNode,projected);

  if ( !CLNode->allAbove( mergePt() ) )return ZERO;
  if ( !Born->xcomb()->willPassCuts() )return ZERO;
  Energy startscale = CKKW_StartScale( Born );
  fillHistory( startscale , Born , CLNode );
  Node->runningPt( fillProjector( projected ? 2 : 1 ) );
  weight *= history.back().weight*alphaReweight()*pdfReweight();
  if( weight == 0. )return ZERO;
  
  CrossSection me = ( inhist?TreedSigDR( startscale , Node ):ZERO );
  CrossSection dip = CLNode->calcDip( startscale*theCurrentME->renFac() );
  
  
  CrossSection res =  weight*as( startscale*DSH()->renFac() )/SM().alphaS()*
  ( double )Node->children().size()*( me -dip );
  if ( Node->legsize() == 6&&Debug::level > 2 ) {
    Energy minPT = Constants::MaxEnergy;
    for( auto const & no :Node->children() )minPT = min( no->pT() , minPT );
    
    generator()->log() << "\nRAA " << minPT/GeV << " " << weight << " " << ( me-dip )/nanobarn << " " << me/nanobarn << " " << dip/nanobarn;
  }
  
  return res;
}

CrossSection Merger::MergingDSigDRRealBelowSubReal( NodePtr Node ){
  NodePtr HistNode = Node->randomChild();
  if ( !HistNode->children().empty() ) {
    auto inoutpair = getInOut( HistNode );
    NodePtr rndCh =  HistNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inoutpair.first , inoutpair.second , rndCh->pT() , theMergePt ) )return ZERO;
  }
  
  NodePtr Born = HistNode-> getHistory( false , DSH()->hardScaleFactor() );
  
  weight = decideClustering(Born,HistNode,projected);

  if ( !Born->xcomb()->willPassCuts() )return ZERO;
  
  Energy startscale = CKKW_StartScale( Born );
  fillHistory( startscale , Born , HistNode );
  Node->runningPt( fillProjector( projected ? 1 : 0 ) );
  weight *= history.back().weight*alphaReweight()*pdfReweight();
  
  if( weight == 0. )return ZERO;
  
  CrossSection sumPS = ZERO;
  
  for( auto const & child : Node->children() ){
    if ( child->allAbove( mergePt() ) ){
      if( ( child )->pT()>mergePt()/3. )//TODO: this is a dynamical cutoff( see below )
      sumPS += child->calcPs( startscale*theCurrentME->renFac() );
      else
      sumPS += child->calcDip( startscale*theCurrentME->renFac() );
    }else{
      assert( child->pT()>mergePt() );
    }
  }
  
  CrossSection me = TreedSigDR( startscale , Node );
  
  if ( Node->legsize() == 6&&Debug::level > 2 ) {
    Energy minPT = Constants::MaxEnergy;
    for( auto const & no :Node->children() )minPT = min( no->pT() , minPT );
    
    generator()->log() << "\nRBSR " << minPT/GeV << " " << weight << " " << ( me-sumPS )/nanobarn << " " << me/nanobarn << " " << sumPS/nanobarn;
  }
    //Here we subtract the PS ( and below the dynamical cutoff the Dip )
  return weight*as( startscale*DSH()->renFac() )/SM().alphaS()*
  ( me-sumPS );
}



CrossSection Merger::MergingDSigDRRealBelowSubInt( NodePtr Node ){
  
  if( Node->children().empty() )return ZERO;
  NodePtr CLNode =  Node->randomChild();
  if( CLNode->pT()<mergePt()/3. )return ZERO;//TODO: this is a dynamical cutoff( see above )
  
  if ( !CLNode->children().empty() ) {
    auto inoutpair = getInOut( CLNode );
    NodePtr rndCh =  CLNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inoutpair.first , inoutpair.second , rndCh->pT() , theMergePt ) )return ZERO;
  }
  
  
  NodePtr Born = CLNode-> getHistory( false , DSH()->hardScaleFactor() );
  
  weight = decideClustering(Born,CLNode,projected);

  if ( !CLNode->allAbove( mergePt() ) )return ZERO;
  
  if ( !Born->xcomb()->willPassCuts() )return ZERO;
  
  Energy startscale = CKKW_StartScale( Born );
  
  fillHistory( startscale , Born , CLNode );
  
  Node->runningPt( fillProjector( projected ? 2 : 1 ) );
  
  weight *= history.back().weight*alphaReweight()*pdfReweight();
  
  if( weight == 0. )return ZERO;
  
  
  pair<CrossSection , CrossSection> DipAndPs =
  CLNode->calcDipandPS( startscale*theCurrentME->renFac() );
  
  if ( Node->legsize() == 6&&Debug::level > 2 ) {
    Energy minPT = Constants::MaxEnergy;
    for( auto const & no :Node->children() )minPT = min( no->pT() , minPT );
    
    
    generator()->log() << "\nRBSI " << CLNode->pT()/GeV << " " << weight << " " << ( DipAndPs.second-DipAndPs.first )/nanobarn << " " << DipAndPs.second/nanobarn << " " << DipAndPs.first/nanobarn;
  }
    //Here we add the PS and subtrac the Dip ( only above the dynamical cutoff )
  return weight*as( startscale*DSH()->renFac() )/SM().alphaS()*
  ( double )Node->children().size()*( DipAndPs.second-DipAndPs.first );
}










CrossSection Merger::TreedSigDR( Energy startscale , NodePtr Node , double diffAlpha ){
  assert( Node->deepHead() == Node );
  
  Node->nodeME()->setScale( sqr( startscale ) , sqr( startscale ) );
  CrossSection res = Node->nodeME()->dSigHatDRB();
  if ( diffAlpha != 1. ) {
    res += Node->nodeME()->dSigHatDRAlphaDiff( diffAlpha );
  }
  if( std::isnan( double( res/nanobarn ) ) ){
    generator()->logWarning(Exception()
                            << "Merger: TreedSigDR is nan"
                            << Exception::warning);
    res = ZERO;};
  return res;
}

CrossSection Merger::LoopdSigDR( Energy startscale , NodePtr Node ){
   // The deephead should be calculated here.
  assert( Node->deepHead() == Node );
    //Node->nodeME()->setXComb( Node->xcomb() );
  Node->nodeME()->setScale( sqr( startscale ) , sqr( startscale ) );
  Node->nodeME()->doOneLoopNoBorn();
  CrossSection res = Node->nodeME()->dSigHatDRV()+Node->nodeME()->dSigHatDRI();
  Node->nodeME()->noOneLoopNoBorn();
  return res;
}

Energy Merger::fillProjector( int pjs ){
    // in the shower handler the scale is multiplied by DSH()->hardScaleFactor()
    // so here we need to devide this
  double xiQSh = history.begin()->node->legsize() == N0()?DSH()->hardScaleFactor():1.;
  if( pjs == 0 ){
    return ( history.size() == 1?1.:( 1./xiQSh ) )*history.back().scale;
  }
  for( auto const & hs : history )
  if ( isProjectorStage( hs.node , pjs )&&pjs !=  0 ){
    history.begin()->node->deepHead()->xcomb()->lastProjector( hs.node->xcomb() );
    return ( hs.node == history[0].node?1.:( 1./xiQSh ) )*hs.scale;
  }
  
  throw Exception() << "Could not fill projector." << Exception::abortnow;
  return ZERO;
}

double Merger::pdfReweight(){
  double res = 1.;
  double facfactor = ( history[0].node->legsize() == N0()?theCurrentME->renFac():DSH()->facFac() );
  for( int side : {0 , 1} ){
    if( history[0].node->xcomb()->mePartonData()[side]->coloured() ){
      for ( auto const & hs : history ){
          //pdf-ratio only to the last step
        if( !hs.node->parent() )continue;
        if ( hs.node == history.back().node )continue;
        if( !hs.node->dipole() ){
          throw Exception() << "\nMerger: pdfReweight: history step has no dipol. "<< Exception::abortnow;
          return 0.;
        }
        res *= pdfratio( hs.node , facfactor*hs.scale , DSH()->facFac()*( hs.node->pT() ) , side );
        facfactor = DSH()->facFac();
      }
      res *= pdfratio( history.back().node , facfactor*history.back().scale  , history[0].scale*theCurrentME->renFac() , side );
    }
  }
  if ( std::isnan( res ) )
     generator()->logWarning(Exception()
                          << "Merger: pdfReweight is nan."
                          << Exception::warning);

  return res;
}

double Merger::alphaReweight(){
  double res = 1.;
  
  Energy Q_R = ( history[0].node->legsize() == N0()?
                theCurrentME->renFac():
                DSH()->renFac() )*
                history[0].scale;
  
  const auto pi=Constants::pi;
  const auto Q_qed=history[0].node->nodeME()->factory()->scaleChoice()->renormalizationScaleQED();
  const auto Oew=history[0].node->nodeME()->orderInAlphaEW();
  const auto Oqcd=history[0].node->nodeME()->orderInAlphaS();
  const auto glfac=3.*( 67./18.-1./6.*pi*pi );
  
  
  res *= pow( as( Q_R ) / SM().alphaS() , Oqcd );
  res *= pow( SM().alphaEMME( Q_qed )/ SM().alphaEMMZ() , Oew );
  
  if ( !( history[0].node->children().empty() ) ){
    res *= pow( ( theCMWScheme?( 1.+ ( glfac-5./9.*Nf( Q_R ) )*as( Q_R )/2./pi ):1. ) , int( history[0].node->legsize()-N0() ) );
  }
  
  
  
  for ( auto const & hs : history )
  if ( hs.node!= history.back().node ){
    Energy q_i = DSH()->renFac()* hs.node->pT();
    res *= as( q_i )/ SM().alphaS()
    *( theCMWScheme?( 1.+ ( glfac-5./9.*Nf( q_i ) )*as( q_i )/2./pi ):1. );
  }
  
  if ( std::isnan( res ) )
      generator()->logWarning(Exception()
               << "Merger: alphaReweight is nan. "<< Exception::warning);
  return res;
}

void Merger::fillHistory( Energy scale , NodePtr Begin , NodePtr EndNode ){
  
  history.clear();
  double sudakov0_n = 1.;
  
  history.push_back( {Begin , sudakov0_n , scale} );
  
  
  double xiQSh = history.begin()->node->legsize() == N0()?DSH()->hardScaleFactor():1.;
  
  scale *= xiQSh;
  if ( Begin->parent()||!isUnitarized ) {
    while ( Begin->parent() && ( Begin !=  EndNode ) ) {
      if ( !dosudakov( Begin , scale , Begin->pT() , sudakov0_n ) ){
        history.push_back( { Begin->parent() , 0. , scale } );
      }
      
      if ( std::isnan( sudakov0_n ) )
      generator()->logWarning(Exception()
                              << "Merger: sudakov"<<scale/GeV<<" "<<Begin->pT()/GeV<<"0_n is nan. "
                              << Exception::warning);
      
      scale = Begin->pT();
      
      history.push_back( { Begin->parent() , sudakov0_n , Begin->pT() } );
      Begin = Begin->parent();
    }
    
    Energy notunirunning = scale;
    
    if ( !isUnitarized&&N()+N0() > int( Begin->deepHead()->legsize() ) ) {
      if ( !dosudakov( Begin , notunirunning , mergePt() , sudakov0_n ) ){
        history.back().weight = 0.;
      }else{
        history.back().weight = sudakov0_n;
      }
    }
  }
  if( history.size() == 1 )scale /= DSH()->hardScaleFactor();
}




double Merger::sumpdfReweightUnlops(){
  double res = 0.;
  Energy beam1Scale = history[0].scale*( history[0].node->legsize() == N0()?theCurrentME->renFac():DSH()->facFac() );
  Energy beam2Scale = history[0].scale*( history[0].node->legsize() == N0()?theCurrentME->renFac():DSH()->facFac() );
  
  for ( auto const & hs : history ){
      //pdf expansion only to the last step
    if( !hs.node->parent() )continue;
    if( hs.node->xcomb()->mePartonData()[0]->coloured()&&hs.node->nodeME()->lastX1()>0.99 )return 0.;
    if( hs.node->xcomb()->mePartonData()[1]->coloured()&&hs.node->nodeME()->lastX2()>0.99 )return 0.;
    
    if( hs.node->nodeME()->lastX1()<0.00001 )return 0.;
    if( hs.node->nodeME()->lastX2()<0.00001 )return 0.;
    
    if ( hs.node->dipole()->bornEmitter() == 0 ){
      res += pdfUnlops( hs.node , 0 , 
                       beam1Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX1() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam1Scale = ( hs.node->pT() )*DSH()->facFac();
    }
    if ( hs.node->dipole()->bornEmitter() == 1 ){
      res += pdfUnlops( hs.node , 1 , 
                       beam2Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX2() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam2Scale = ( hs.node->pT() )*DSH()->facFac();
    }
    if ( hs.node->dipole()->bornSpectator() == 0 &&
         hs.node->dipole()->bornEmitter() >1 ){
      res += pdfUnlops( hs.node , 0 , 
                       beam1Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX1() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam1Scale = ( hs.node->pT() )*DSH()->facFac();
    }
    if ( hs.node->dipole()->bornSpectator() == 1 &&
         hs.node->dipole()->bornEmitter() >1 ){
      res += pdfUnlops( hs.node , 1 , 
                       beam2Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX2() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam2Scale = ( hs.node->pT() )*DSH()->facFac();
    }
  }
  
  if ( history[0].node->deepHead()->xcomb()->mePartonData()[0]->coloured() ){
    res += pdfUnlops( history.back().node , 0 , 
                     beam1Scale , 
                     history[0].scale*theCurrentME->renFac() , 
                     ( history.back() ).node->nodeME()->lastX1() , 
                     Nf( history[0].scale ) , 
                     history[0].scale );
    
  }
  if ( history[0].node->deepHead()->xcomb()->mePartonData()[1]->coloured() ) {
    res += pdfUnlops( history.back().node , 1 , 
                     beam2Scale , 
                     history[0].scale*theCurrentME->renFac() , 
                     ( history.back() ).node->nodeME()->lastX2() , 
                     Nf( history[0].scale ) , 
                     history[0].scale );
  }
  return res;
}





#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"
double Merger::pdfUnlops( NodePtr node  , int side , Energy  running , Energy next , double x , int nlp , Energy fixedScale )  {
  
  tcPDPtr particle , parton;
  tcPDFPtr pdf;
  if ( side == 0 ) {
    particle = node->nodeME()->lastParticles().first->dataPtr();
    parton = node->nodeME()->lastPartons().first->dataPtr();
    pdf  = node->xcomb()->partonBins().first->pdf();
  }else{
    assert( side == 1 );
    particle = node->nodeME()->lastParticles().second->dataPtr();
    parton = node->nodeME()->lastPartons().second->dataPtr();
    pdf  = node->xcomb()->partonBins().second->pdf();
  }
  
    //copied from PKOperator
  double res = 0.;
  int number = 10;
  for ( int nr = 0;nr<number;nr++ ){
    double restmp = 0;
    
    using namespace RandomHelpers;
    double r = UseRandom::rnd();
    double eps = 1e-3;
    
    pair<double , double> zw  =
    generate( ( piecewise() , 
               flat( 0.0 , x ) , 
               match( inverse( 0.0 , x , 1.0 ) + inverse( 1.0+eps , x , 1.0 ) ) ) , r );
    
    double z = zw.first;
    double mapz = zw.second;
    double PDFxparton = pdf->xfx( particle , parton , sqr( fixedScale ) , x )/x;
    double CA = SM().Nc();
    double CF = ( SM().Nc()*SM().Nc()-1.0 )/( 2.*SM().Nc() );
    
    
    if ( abs( parton->id() ) < 7 ) {
      
      double PDFxByzgluon = pdf->xfx( particle , getParticleData( ParticleID::g ) , sqr( fixedScale ) , x/z )*z/x;
      double PDFxByzparton = pdf->xfx( particle , parton , sqr( fixedScale ) , x/z )*z/x;
      assert( abs( parton->id() ) < 7 );
      
      restmp += CF*( 3./2.+2.*log( 1.-x ) ) * PDFxparton;
      if ( z > x ) {
        restmp += 0.5 * ( sqr( z ) + sqr( 1.-z ) ) * PDFxByzgluon / z;
        restmp += CF*2.*( PDFxByzparton - z*PDFxparton )/( z*( 1.-z ) );
        restmp -=  CF*PDFxByzparton * ( 1.+z )/z;
      }
    }else{
      
      assert( parton->id() == ParticleID::g );
      double PDFxByzgluon = pdf->xfx( particle , getParticleData( ParticleID::g ) , sqr( fixedScale ) , x/z )*z/x;
      if ( z > x ){
        double factor = CF * ( 1. + sqr( 1.-z ) ) / sqr( z );
        
        
        for ( int f = -nlp; f <=  nlp; ++f ) {
          if ( f == 0 )
          continue;
          restmp += pdf->xfx( particle , getParticleData( f ) , sqr( fixedScale ) , x/z )*z/x*factor;
        }
      }
      
      restmp += ( ( 11./6. ) * CA - ( 1./3. )*Nf( history[0].scale ) + 2.*CA*log( 1.-x ) ) *PDFxparton;
      if ( z > x ) {
        restmp += 2. * CA * ( PDFxByzgluon - z*PDFxparton ) / ( z*( 1.-z ) );
        restmp += 2.* CA *( ( 1.-z )/z - 1. + z*( 1.-z ) ) * PDFxByzgluon / z;
      }
      
    }
    if ( PDFxparton<1e-8 )restmp =  0.;
    res += 1*restmp*log( sqr( running/next ) )/PDFxparton*mapz;
    
  }
  return res/number;
}



double Merger::sumalphaReweightUnlops(){
  double res = 0.;
  if ( !( history[0].node->children().empty() ) ){
    res += alphasUnlops( history[0].scale*DSH()->renFac() ,
                        history[0].scale*theCurrentME->renFac() )*
                        int( history[0].node->legsize()-N0() );
  }
    // dsig is calculated with fixed alpha_s
  for ( auto const & hs : history ){
      //expansion only to the last step
    if( !hs.node->parent() )continue;
    res += alphasUnlops( hs.node->pT()*DSH()->renFac()  , history[0].scale );
  }
  return res;
}

double Merger::sumfillHistoryUnlops(){
  double res = 0.;
  double xiQSh = history[0].node->legsize() == N0()?DSH()->hardScaleFactor():1.;
  for ( auto const & hs : history ){
    if( !hs.node->parent() )continue;
    doUNLOPS( hs.node , ( hs.node == history[0].node?xiQSh:1. )*hs.scale , hs.node->pT() , history[0].scale , res );
  }
  return res;
}


MergingFactoryPtr Merger::treefactory(){return theTreeFactory;}


void Merger::doinit(){
  if ( !DSH()->hardScaleIsMuF() ) {
    throw Exception()
    << "Merger: Merging is currently only sensible if we are using the hardScale as MuF."
    << Exception::abortnow;
  }
}


CrossSection Merger::MergingDSigDR() {
  
  assert(history.empty());
  
  if ( theFirstNodeMap.find( theCurrentME ) == theFirstNodeMap.end() ) {
    throw Exception()
    << "Merger: The current matrix element could not be found."
    << Exception::abortnow;
  }
  
  
  NodePtr Node = theFirstNodeMap[theCurrentME];
    //DSH()->setCurrentHandler();
  DSH()->eventHandler( generator()->eventHandler() );
  
  CrossSection res = ZERO;
  if( Node->deepHead()->subtractedReal() ){
    res = MergingDSigDRRealStandard( Node );
    theCurrentMaxLegs = maxLegsNLO();
  }else if( Node->deepHead()->virtualContribution() ){
    res = MergingDSigDRVirtualStandard( Node );
    theCurrentMaxLegs = maxLegsNLO();
  }else if (true){
    res = MergingDSigDRBornCheapME( Node );
    theCurrentMaxLegs = maxLegsLO();
  }else if( theGamma!= 1. ){
    res = MergingDSigDRBornGamma( Node );
    theCurrentMaxLegs = maxLegsLO();
  }else{
    res = MergingDSigDRBornStandard( Node );
    theCurrentMaxLegs = maxLegsLO();
  }
  
  auto lxc=theCurrentME->lastXCombPtr();
    //lxc->lastCentralScale( sqr( Node->runningPt() ) );
  lxc->lastShowerScale( sqr( Node->runningPt() ) );
  
  auto lp=theCurrentME->lastXCombPtr()->lastProjector();
  if( lp ){
      //lp->lastCentralScale( sqr( Node->runningPt() ) );
    lp->lastShowerScale( sqr( Node->runningPt() ) );
  }
  
  if ( res == ZERO ){
    history.clear();
    return res;
  }
  
  
  cleanup( Node );
  lxc->subProcess( SubProPtr() );
  
  history.clear();
  
  if( std::isnan( double( res/nanobarn ) )|| !std::isfinite( double( res/nanobarn ) ) ){
    generator()->logWarning(Exception()
       << "Merger weight is " << res/nanobarn<< " -> setting to 0"
       << Exception::warning);
    return ZERO;
  }
  
  return res;
  
}


#include "Herwig/PDF/HwRemDecayer.h"
void Merger::CKKW_PrepareSudakov( NodePtr Born , Energy running ){
    //cleanup( Born );
  tSubProPtr sub = Born->xcomb()->construct();
  DSH()->resetPDFs( make_pair( Born->xcomb()->partonBins().first->pdf() , 
                              Born->xcomb()->partonBins().second->pdf() ) , 
                   Born->xcomb()->partonBins() );
  DSH()->setCurrentHandler();
  
  DSH()->currentHandler()->generator()->currentEventHandler( Born->deepHead()->xcomb()->eventHandlerPtr() );
  
  DSH()->currentHandler()->remnantDecayer()->setHadronContent( Born->deepHead()->xcomb()->lastParticles() );
  DSH()->eventRecord().clear();
  DSH()->eventRecord().slimprepare( sub , dynamic_ptr_cast<tStdXCombPtr>( Born->xcomb() ) , DSH()->pdfs() , Born->deepHead()->xcomb()->lastParticles() );
  DSH()->hardScales( sqr( running ) );
}


Energy Merger::CKKW_StartScale( NodePtr Born ){
  Energy res = generator()->maximumCMEnergy();;
  if( !Born->children().empty() ){
    for ( int i = 0;i<Born->legsize();i++ ){
      if ( !Born->nodeME()->mePartonData()[i]->coloured() )continue;
      for ( int j = i+1;j<Born->legsize();j++ ){
        if ( i == j||!Born->nodeME()->mePartonData()[j]->coloured() )continue;
        res =  min( res , sqrt( 2.*Born->nodeME()->lastMEMomenta()[i]*Born->nodeME()->lastMEMomenta()[j] ) );
      }
    }
  }else{
    Born->nodeME()->factory()->scaleChoice()->setXComb( Born->xcomb() );
    res =  sqrt( Born->nodeME()->factory()->scaleChoice()->renormalizationScale() );
  }
  Born->nodeME()->factory()->scaleChoice()->setXComb( Born->xcomb() );
  res = max( res , sqrt( Born->nodeME()->factory()->scaleChoice()->renormalizationScale() ) );
  return res;
}




double Merger::alphasUnlops( Energy next , Energy fixedScale )  {
  double betaZero = ( 11./6. )*SM().Nc() - ( 1./3. )*Nf( history[0].scale );
  return ( betaZero*log( sqr( fixedScale/next ) ) )+
  ( theCMWScheme?( ( ( 3.*( 67./18.-1./6.*Constants::pi*Constants::pi )-5./9.*Nf( history[0].scale ) ) ) ):0. );
}


double Merger::pdfratio( NodePtr  Born , Energy  nominator_scale , Energy denominator_scale , int side ){
  
  StdXCombPtr bXc = Born->xcomb();
  if( !bXc->mePartonData()[side]->coloured() )
  throw Exception()
  << "Merger: pdf-ratio required for non-coloured particle."
  << Exception::abortnow;
  
  
  double from , to;
  if ( side == 0 ){
    if ( denominator_scale == nominator_scale ) {
      return 1.;
    }
    from = Born->nodeME()->pdf1( sqr( denominator_scale ) );
    to   = Born->nodeME()->pdf1( sqr( nominator_scale ) );
    if ( ( to < 1e-8||from < 1e-8 )&&( to/from>10000000. ) ){
      generator()->logWarning(Exception()
                              << "Merger: pdfratio to = " << to << " from = " << from
                              << Exception::warning);
      return 0.;
    }
  }
  else{
    if ( denominator_scale == nominator_scale ) {
      return 1.;
    }
    from = Born->nodeME()->pdf2( sqr( denominator_scale ) );
    to = Born->nodeME()->pdf2( sqr( nominator_scale ) );
    if ( ( to < 1e-8||from < 1e-8 )&&( to/from>10000000. ) ){
      generator()->logWarning(Exception()
                              << "Merger: pdfratio to = " << to << " from = " << from
                              << Exception::warning);
      return 0.;}
  }
  return to/from;
}



bool Merger::dosudakov( NodePtr Born , Energy running , Energy next , double& sudakov0_n ) {
  CKKW_PrepareSudakov( Born , running );
  for( DipoleChain const & chain : DSH()->eventRecord().chains() ){
    for( Dipole const & dip : chain.dipoles() ){
      sudakov0_n *= singlesudakov( dip , next , running , make_pair( true , false ) );
      sudakov0_n *= singlesudakov( dip , next , running , make_pair( false , true ) );
      if ( sudakov0_n == 0.0 ){
          cleanup( Born );
        return false;
      }
    }
  }
    cleanup( Born );
  return true;
}

bool Merger::doUNLOPS( NodePtr Born , Energy  running , Energy next , Energy fixedScale , double& UNLOPS ) {
  CKKW_PrepareSudakov( Born , running );
  for( DipoleChain const & chain : DSH()->eventRecord().chains() ){
    for( Dipole const & dip : chain.dipoles() ){
      UNLOPS += singleUNLOPS( dip , next , running , fixedScale , make_pair( true , false ) );;
      UNLOPS += singleUNLOPS( dip , next , running , fixedScale , make_pair( false , true ) );
    }
  }
  cleanup( Born );
  return true;
}



bool Merger::isProjectorStage( NodePtr  Born , int pjs ){
  return ( pjs == int( ( Born->deepHead()->legsize() - Born->legsize() ) ) );
}

void Merger::cleanup( NodePtr Born ) {
  DSH()->eventRecord().clear();
  if( !Born->xcomb()->subProcess() )return;
  ParticleVector vecfirst = Born->xcomb()->subProcess()->incoming().first->children();
  for( auto const & particle : vecfirst )
  Born->xcomb()->subProcess()->incoming().first->abandonChild( particle );
  
  ParticleVector vecsecond = Born->xcomb()->subProcess()->incoming().second->children();
  for( auto const & particle : vecsecond )
  Born->xcomb()->subProcess()->incoming().second->abandonChild( particle );
  Born->xcomb()->subProcess( SubProPtr() );
}

double Merger::singlesudakov( Dipole dip  , Energy next , Energy running , pair<bool , bool> conf ){
  
  double res = 1.;
  tPPtr emitter = dip.emitter( conf );
  tPPtr spectator = dip.spectator( conf );
  DipoleSplittingInfo candidate( dip.index( conf ) , conf , 
                                dip.emitterX( conf ) , 
                                dip.spectatorX( conf ) , 
                                emitter , spectator );
  
  
  if ( DSH()->generators().find( candidate.index() ) == DSH()->generators().end() ) DSH()->getGenerators( candidate.index() );
  
  auto const & gens = DSH()->generators().equal_range( candidate.index() );
  
  for ( auto   gen = gens.first; gen !=  gens.second; ++gen ) {
    if ( !( gen->first == candidate.index() ) )
    continue;
    
    Energy dScale  = 	gen->second->splittingKinematics()->dipoleScale( emitter->momentum() , spectator->momentum() );
    candidate.scale( dScale );
    candidate.continuesEvolving();
    Energy ptMax = gen->second->splittingKinematics()->ptMax( candidate.scale() , candidate.emitterX() , candidate.spectatorX() , 
                                                             candidate.index() , *gen->second->splittingKernel() );
    
    candidate.hardPt( min( running , ptMax ) );
    
    if ( candidate.hardPt()>next ){
      res *= gen->second->sudakov( candidate , next );
    }
  }
  
  return res;
}


double Merger::singleUNLOPS( Dipole dip  , Energy next , Energy running , Energy fixedScale , pair<bool , bool> conf ){
  
  double res = 0.;
  tPPtr emitter = dip.emitter( conf );
  tPPtr spectator = dip.spectator( conf );
  DipoleSplittingInfo candidate( dip.index( conf ) , 
                                 conf , dip.emitterX( conf ) , 
                                 dip.spectatorX( conf ) , 
                                 emitter , spectator );
  
  if ( DSH()->generators().find( candidate.index() ) ==
       DSH()->generators().end() )
       DSH()->getGenerators( candidate.index() );
  
  auto const & gens = DSH()->generators().equal_range( candidate.index() );
  
  for ( auto gen = gens.first; gen !=  gens.second; ++gen ) {
    if ( !( gen->first == candidate.index() ) )
    continue;
    Energy dScale = gen->second->splittingKinematics()->dipoleScale(
                          emitter->momentum() , spectator->momentum() );
    candidate.scale( dScale );
    candidate.continuesEvolving();
    Energy ptMax = gen->second->
                   splittingKinematics()->ptMax(
                      candidate.scale() , candidate.emitterX() , 
                      candidate.spectatorX() , candidate.index() , 
                      *gen->second->splittingKernel() );
    
    candidate.hardPt( min( running , ptMax ) );
    if ( candidate.hardPt()>next ){
      res += gen->second->unlops( candidate , next , fixedScale );
    }
  }
	 
	 return res;
}




void Merger::firstNodeMap( MatchboxMEBasePtr a , NodePtr b ){theFirstNodeMap.insert( make_pair( a , b ) );}



map<MatchboxMEBasePtr , NodePtr> Merger::firstNodeMap(){return theFirstNodeMap;}




void Merger::setXComb( MatchboxMEBasePtr me , tStdXCombPtr xc ){
  theFirstNodeMap[me]->setXComb( xc );
}
void Merger::setKinematics( MatchboxMEBasePtr me ){
  theFirstNodeMap[me]->setKinematics();
}
void Merger::clearKinematics( MatchboxMEBasePtr me ){
  theFirstNodeMap[me]->clearKinematics();
}
bool Merger::generateKinematics( MatchboxMEBasePtr me , const double * r ){
  
  theFirstNodeMap[me]->firstgenerateKinematics( r , 0 );
  
  if ( theFirstNodeMap[me]->cutStage() == 0 ){
    
    bool inAlphaPS = false;
    NodePtrVec children = theFirstNodeMap[me]->children();
    for ( NodePtr const & child: children ) {
      treefactory()->setAlphaParameter( theGamma );
      inAlphaPS |=  theGamma!= 1.?child->dipole()->aboveAlpha():false;
      treefactory()->setAlphaParameter( 1. );
    }
    
    SafeClusterMap temp = theFirstNodeMap[me]->clusterSafe();
    for ( auto const & cs: temp ) {
      if ( !cs.second.first&&!inAlphaPS )return false;
    }
  }
  if ( theFirstNodeMap[me]->cutStage() == 1 ){
    SafeClusterMap temp = theFirstNodeMap[me]->clusterSafe();
    for ( auto const & sc: temp ) {
      if ( !sc.second.first && !sc.second.second )return false;
    }
  }
  return true;
  
}



void Merger::fillProjectors( MatchboxMEBasePtr me ){
  for ( auto const & propair: theFirstNodeMap[me]->Projector() ) {
    me->lastXCombPtr()->projectors().insert( propair.first , 
                                            propair.second->xcomb() );
  }
}
pair<bool , bool> Merger::clusterSafe( MatchboxMEBasePtr me , int emit , int emis , int spec ){
  return theFirstNodeMap[me]->clusterSafe().find({emit, emis, spec})->second;
  
}

  //#include "fastjet/ClusterSequence.hh"
bool Merger::matrixElementRegion( PVector incoming , PVector outgoing , Energy winnerScale , Energy cutscale ){
  
    //generator()->log() << "\nparticles s" << particles.size() << " " << particles[0] << " " << particles[1] << flush;
  /*
   if ( defMERegionByJetAlg && !particles[0]->coloured()&& !particles[1]->coloured() ) {
   assert( false );
   vector<fastjet::PseudoJet> input_particles;
   for( size_t em = 2; em < particles.size();em++ ){
   input_particles.push_back( fastjet::PseudoJet( em->momentum().x()/GeV , 
   em->momentum().y()/GeV , 
   em->momentum().z()/GeV , 
   em->momentum().e()/GeV ) );
   }
   fastjet::JetDefinition jet_def( fastjet::ee_kt_algorithm );
   fastjet::ClusterSequence clust_seq( input_particles , jet_def );
   size_t n = particles.size()-2;
   vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets_ycut( ee_ycut );
   return n == exclusive_jets.size();
   }
   
   
   if ( defMERegionByJetAlg ) {
   assert( false );
   size_t noncol = 0;
   vector<fastjet::PseudoJet> input_particles;
   for( size_t em = 2; em < particles.size();em++ ){
   if ( em->coloured() )
   input_particles.push_back( fastjet::PseudoJet( em->momentum().x()/GeV , 
   em->momentum().y()/GeV , 
   em->momentum().z()/GeV , 
   em->momentum().e()/GeV ) );
   else
   noncol++;
   }
   double Rparam = 1.0;
   fastjet::Strategy strategy = fastjet::Best;
   fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
   fastjet::JetDefinition jet_def( fastjet::kt_algorithm , Rparam , recomb_scheme , strategy );
   
   fastjet::ClusterSequence clust_seq( input_particles , jet_def );
   size_t n = particles.size()-2-noncol;
   vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets( pp_dcut );
   
   // generator()->log() << "\nn =  " << n << " jets =  " << exclusive_jets.size();
   //for ( size_t i = 0; i<exclusive_jets.size(); i++ ) {
   //generator()->log() << "\nn =  " << n << " jetspt =  " << exclusive_jets[i].perp();
   //}
   
   return n == exclusive_jets.size();
   }
   
   */
  
  
  
  
  Energy ptx = Constants::MaxEnergy;
  bool foundwinnerpt = false;
  using namespace boost;
    //FF
  
  for( auto const & em : outgoing ){ if ( ! em->coloured() ) continue;
    for( auto const & emm : outgoing ){ if ( !emm->coloured() ) continue; if ( em == emm ) continue;
      for( auto const & spe : outgoing ){ if ( !spe->coloured() ) continue; if ( em == spe||emm == spe ) continue;
        
        if ( !( em->id() == -emm->id()||emm->id()>6 ) )continue;
        
        Energy pt = ZERO;
        if ( em->momentum().m()<= 0.001*GeV&&
            emm->momentum().m()<= 0.001*GeV&&
            spe->momentum().m()<= 0.001*GeV ) {
          pt = FFLTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }else{
          pt = FFMTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }
        
        if ( abs( pt-winnerScale )<0.001*GeV ) {
          foundwinnerpt = true;
        }
        ptx  = min( ptx , pt );
      }
    }
  }
  
    //FI
  for( auto const & spe : incoming ){ if ( ! spe->coloured() ) continue;
    for( auto const & emm : outgoing ){ if ( ! emm->coloured() ) continue;
      for( auto const & em : outgoing ){ if ( ! em->coloured() ) continue; if ( em == emm ) continue;
        
        if ( !( em->id() == -emm->id() || emm->id()>6 ) )continue;
        

        Energy pt = ZERO;
        if ( em->momentum().m()<= 0.001*GeV&&
            emm->momentum().m()<= 0.001*GeV&&
            spe->momentum().m()<= 0.001*GeV ) {
          pt = FILTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }else{
          pt = FIMTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }
        
        
        if ( abs( pt-winnerScale )<0.001*GeV ) {
          foundwinnerpt = true;
        }
        
        if( pt > 0.*GeV )
        ptx  = min( ptx , pt );
      }
    }
  }
  
    //IF
  for( auto const & em : incoming ){ if ( ! em->coloured() ) continue;
    for( auto const & emm : outgoing ){ if ( ! emm->coloured() ) continue;
      for( auto const & spe : outgoing ){ if ( ! spe->coloured() ) continue; if ( emm == spe ) continue;
        
        if ( !( em->id()>6|| em->id() == emm->id() ||emm->id()>6 ) )continue;
        
        Energy pt = ZERO;
        
        if ( em->momentum().m()<= 0.001*GeV&&
             emm->momentum().m()<= 0.001*GeV&&
             spe->momentum().m()<= 0.001*GeV ) {
            //massless
          pt = IFLTK->lastPt( em->momentum() , emm->momentum() , spe->momentum()  );
        }else{
            //massiv
          pt = IFMTK->lastPt( em->momentum() , emm->momentum() , spe->momentum()  );
        }
        
        
        if ( abs( pt-winnerScale )<0.01*GeV ) {
          foundwinnerpt = true;
        }
        ptx  = min( ptx , pt );
      }
    }
  }
  
    //II
  for( auto const & em : incoming ){ if ( ! em->coloured() ) continue;
    for( auto const & spe : incoming ){ if ( ! spe->coloured() ) continue; if ( em == spe ) continue;
      for( auto const & emm : outgoing ){ if ( ! emm->coloured() ) continue;
        if ( !( em->id()>6||em->id() == emm->id() ||emm->id()>6 ) )continue;
        
        Energy  pt = IILTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        
        if ( abs( pt-winnerScale )<0.01*GeV ) {
          foundwinnerpt = true;
        }
        ptx  = min( ptx , pt );
      }
    }
  }
  
  if( !foundwinnerpt ){
    generator()->logWarning( Exception()
                            << "Merger: Could not find winner with pt."
                            << "Run with -d2 to get phase space points. "
                            << Exception::warning );
    
    if( Debug::level > 2 ) {
      generator()->log() << "\nWarning: could not find winner with pt: " << winnerScale/GeV;
      for( auto const & emm : incoming ){
        generator()->log() << "\n" << emm->id() << " " << emm->momentum()/GeV << " " << emm->momentum().m()/GeV << flush;
      }
      for( auto const & emm : outgoing ){
        generator()->log() <<"\n" << emm->id() << " " << emm->momentum()/GeV << " " << emm->momentum().m()/GeV << flush;
      }
    }
    
  }
  
  return ( ptx>cutscale ) ;
  
}

int Merger::M()const{return theTreeFactory->M();}

int Merger::N()const{return theTreeFactory->N();}


  // If needed , insert default implementations of virtual function defined
  // in the InterfacedBase class here ( using ThePEG-interfaced-impl in Emacs ).

#include "ThePEG/Persistency/PersistentOStream.h"
void Merger::persistentOutput( PersistentOStream & os ) const {
  os << Unlopsweights << theCMWScheme << projected <<
  isUnitarized << isNLOUnitarized << defMERegionByJetAlg <<
  theOpenInitialStateZ << theChooseHistory <<
  theN0 << theOnlyN   << weight <<
  weightCB << theGamma << ee_ycut << pp_dcut <<
  theSmearing <<  ounit( theIRSafePT , GeV ) <<
  ounit( theMergePt , GeV ) <<  ounit( theCentralMergePt , GeV ) <<
  theMergingJetFinder << theLargeNBasis << FFLTK << FILTK <<
  IFLTK << IILTK << FFMTK << FIMTK << IFMTK <<
  theDipoleShowerHandler << theTreeFactory << theFirstNodeMap;
  
}
#include "ThePEG/Persistency/PersistentIStream.h"
void Merger::persistentInput( PersistentIStream & is , int ) {
  is >> Unlopsweights >> theCMWScheme >> projected >>
  isUnitarized >> isNLOUnitarized >>
  defMERegionByJetAlg >> theOpenInitialStateZ >>
  theChooseHistory >> theN0 >> theOnlyN >>
  weight >> weightCB >>
  theGamma >> ee_ycut >> pp_dcut >>
  theSmearing >>  iunit( theIRSafePT , GeV ) >>
  iunit( theMergePt , GeV ) >>
  iunit( theCentralMergePt , GeV ) >>
  theMergingJetFinder >> theLargeNBasis >>
  FFLTK >> FILTK >> IFLTK >>
  IILTK >> FFMTK >> FIMTK >>
  IFMTK >> theDipoleShowerHandler >>
  theTreeFactory >> theFirstNodeMap ;
}

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct ( the class and its base class ) , and that the constructor
  // arguments are correct ( the class name and the name of the dynamically
  // loadable library where the class implementation can be found ).
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<Merger , Herwig::MergerBase>
describeHerwigMerger( "Herwig::Merger" , "HwDipoleShower.so" );

  //ClassDescription<Merger> Merger::initMerger;
  // Definition of the static class description member.

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"

void Merger::Init() {
  
  static ClassDocumentation<Merger> documentation
  ( "Merger." );
  
  
  
  static Reference<Merger , DipoleShowerHandler> interfaceShowerHandler
  ( "DipoleShowerHandler" , 
   "" , 
   &Merger::theDipoleShowerHandler , false , false , true , true , false );
  
  
  
  static Switch<Merger , bool>
  interfaceUnlopsweights ( "Unlopsweights" , "" , &Merger::Unlopsweights , false , false , false );
  static SwitchOption interfaceUnlopsweightsYes
  ( interfaceUnlopsweights , "Yes" , "" , true );
  static SwitchOption interfaceUnlopsweightsNo
  ( interfaceUnlopsweights , "No" , "" , false );
  
  static Switch<Merger , bool>
  interfacetheCMWScheme ( "CMWScheme" , "" , &Merger::theCMWScheme , false , false , false );
  static SwitchOption interfacetheCMWSchemeYes
  ( interfacetheCMWScheme , "Yes" , "" , true );
  static SwitchOption interfacetheCMWSchemeNo
  ( interfacetheCMWScheme , "No" , "" , false );
  
  static Parameter<Merger , Energy> interfaceMergerScale
  ( "MergingScale" , 
   "The MergingScale." , 
   &Merger::theCentralMergePt , GeV , 20.0*GeV , 0.0*GeV , 0*GeV , 
   false , false , Interface::lowerlim );
  
  static Reference<Merger , MergingFactory> interfaceMergerHelper
  ( "MergingFactory" , 
   "" , 
   &Merger::theTreeFactory , false , false , true , true , false );
  
  static Parameter<Merger , double> interfacedcut
  ( "pp_dcut" , 
   "The MergingScale." , 
   &Merger::pp_dcut , 5.0 , 0.0 , 0 , 
   false , false , Interface::lowerlim );
  
  static Parameter<Merger , double> interfaceycut
  ( "ee_ycut" , 
   "The MergingScale." , 
   &Merger::ee_ycut , 0.2 , 0.0 , 0 , 
   false , false , Interface::lowerlim );
  
  static Parameter<Merger , double> interfacegamma
  ( "gamma" , 
   "gamma parameter." , 
   &Merger::theGamma , 1.0 , 0.0 , 0 , 
   false , false , Interface::lowerlim );
  
  static Reference<Merger , JetFinder> interfaceMergingJetFinder
  ( "MergingJetFinder" , 
   "A reference to the jet finder used in Merging." , 
   &Merger::theMergingJetFinder , false , false , true , false , false );
  
  
  
  static Reference<Merger , ColourBasis> interfaceLargeNBasis
  ( "LargeNBasis" , 
   "Set the large-N colour basis implementation." , 
   &Merger::theLargeNBasis , false , false , true , true , false );
  
  
  
  
  
  
  static Switch<Merger , bool>
  interfacedefMERegionByJetAlg
  ( "MERegionByJetAlg" , "" , &Merger::defMERegionByJetAlg , false , false , false );
  
  static SwitchOption interfacedefMERegionByJetAlgYes
  ( interfacedefMERegionByJetAlg , "Yes" , "" , true );
  static SwitchOption interfacedefMERegionByJetAlgNo
  ( interfacedefMERegionByJetAlg , "No" , "" , false );
  
  
  static Switch<Merger , bool>
  interfaceOpenInitialSateZ
  ( "OpenInitialStateZ" , "" , &Merger::theOpenInitialStateZ , false , false , false );
  
  static SwitchOption interfaceOpenInitialSateZYes
  ( interfaceOpenInitialSateZ , "Yes" , "" , true );
  static SwitchOption interfaceOpenInitialSateZNo
  ( interfaceOpenInitialSateZ , "No" , "" , false );
  
  
  
  
  
  
  
  
  static Parameter<Merger , Energy>
  interfaceIRSafePT
  ( "IRSafePT" , "Set the pt for which a matrixelement should be treated as IR-safe." , 
   
   &Merger::theIRSafePT , 
   GeV , 0.0 * GeV , ZERO , Constants::MaxEnergy , true , false , Interface::limited );
  interfaceIRSafePT.setHasDefault( false );
  
  
  
  static Parameter<Merger , double> interfacemergePtsmearing( "MergingScaleSmearing" , "Set the percentage the merging pt should be smeared." , 
                                                            &Merger::theSmearing , 0. , 0. , 
                                                            0.0 , 0.5 , true , false , Interface::limited );
  
  
  
  static Parameter<Merger , int> interfacechooseHistory( "chooseHistory" , "different ways of choosing the history" , &Merger::theChooseHistory , 3 , -1 , 0 , 
                                                       false , false , Interface::lowerlim );
  
  
  
  
  
  
}
















