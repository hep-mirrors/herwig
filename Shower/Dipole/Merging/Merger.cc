  // -*- C++ -*-
  //
  // Merger.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL , see COPYING for details.
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

IBPtr Merger::clone() const {
  return new_ptr( *this );
}

IBPtr Merger::fullclone() const {
  return new_ptr( *this );
}





namespace {
  double decideClustering(const NodePtr sub,const NodePtr head,bool& pro){
    if( sub !=  head ){// at least one history step -> unitarisation
      if ( UseRandom::rndbool() ){ pro = true;    return -2.;   }
      else{                        pro = false;   return  2.;   }
    }                  // no ordered history        -> no projection
    else{                          pro = false;   return  1.;   }
  }
}

CrossSection Merger::MergingDSigDRBornStandard( ){
    // get the history for the process
  const NodePtr productionNode =
         currentNode()-> getHistory( true, DSH()->hardScaleFactor() );
    // decide if to cluster
  weight = decideClustering(productionNode, currentNode(), projected);
    // check if we only want to calculate the current multiplicity.
  if(notOnlyMulti())                                      return ZERO;
    // Check for cuts on the production proces.
  if ( !productionNode->xcomb()->willPassCuts() )         return ZERO;
    // calculate the staring scale for the production node
  Energy startscale = CKKW_StartScale( productionNode );
    // fill history with caluclation of sudakov supression
  fillHistory( startscale , productionNode , currentNode() );
    // fill the projector -> return the scale of the last splitting
  currentNode()->runningPt( fillProjector( projected ? 1 : 0 ) );
    // the weight has three components to get shower history weight
  weight *= history.back().weight*  // Sudakov suppression
            alphaReweight()*        // alpha_s reweight
            pdfReweight();          // pdf reweight
    // If weight is zero return.
  if( weight == ZERO )                                    return ZERO;
    //calculate the cross section
  return weight*TreedSigDR( startscale ,  1. );
}




CrossSection Merger::MergingDSigDRVirtualStandard( ){
   // get the history for the process
  const NodePtr productionNode =  currentNode()-> getHistory( true , DSH()->hardScaleFactor() );
   // decide if to cluster
  weight = decideClustering(productionNode,currentNode(),projected);
    // Check for cuts on the production proces.
  if ( !productionNode->xcomb()->willPassCuts() )return ZERO;
    // calculate the staring scale
  Energy startscale = CKKW_StartScale( productionNode );
    // fill history with caluclation of sudakov supression
  fillHistory( startscale , productionNode , currentNode() );
    // fill the projector -> return the scale of the last splitting
  currentNode()->runningPt( fillProjector( projected ? 1 : 0 ) );
    // the weight has three components to get shower history weight
  double ww1 = history.back().weight;
  double ww2 = alphaReweight(true);
  double ww3 = pdfReweight();
  weight *= ww1*ww2*ww3;
    // If weight is zero return.
  if( weight == 0. )return ZERO;
    // calculate the cross section for virtual contribution
    // and various insertion operators.
  CrossSection matrixElement = LoopdSigDR( startscale );
    // Now calculate the expansion of the shower history.
    // first: the born contibution:
  CrossSection bornWeight =  currentME()->dSigHatDRB();
    // second: expansion of pdf ,alpha_s-ratio and sudakov suppression.
  double w1 = -sumPdfReweightExpansion();
  double w2 = -sumAlphaSReweightExpansion();
  double w3 = -sumFillHistoryExpansion();
    // put together the expansion weights.
  CrossSection expansionweight  = 
	bornWeight*SM().alphaS()/( 2.*ThePEG::Constants::pi );



  if (theShowerExpansionWeights == 0){
	expansionweight  *=0.;
  }else if ( theShowerExpansionWeights == 1 ){
	expansionweight  *=w1+w2+w3;
  }else if ( theShowerExpansionWeights == 2 ){
  	expansionweight  *=w1+w2+w3*pow(as( startscale*DSH()->renFac() )/SM().alphaS(),currentME()->orderInAlphaS())/ww2;
  }else if ( theShowerExpansionWeights == 3 ){
        expansionweight  *=(w1+w2+w3)*pow(as( startscale*DSH()->renFac() )/SM().alphaS(),currentME()->orderInAlphaS())/ww2;
  }else if ( theShowerExpansionWeights == 4 ){
        expansionweight  *= w1+w3+w2*pow(as( startscale*DSH()->renFac() )/SM().alphaS(),currentME()->orderInAlphaS())/ww2;
  }else assert(false && theShowerExpansionWeights);

    // [ DEBUG ]
  if ( currentNode()->legsize() == 5 && Debug::level > 2 )
    debugVirt(weight,w1,w2,w3,matrixElement,ww1,ww2,ww3,productionNode,bornWeight);
    // return with correction that ME was calculated with fixed alpha_s
  return weight* as( startscale*DSH()->renFac() )/
         SM().alphaS()* ( matrixElement+expansionweight );
}




CrossSection Merger::MergingDSigDRRealStandard(){
  if ( currentNode()->children().empty() ) {
    throw Exception()
    << "Real emission contribution without underlying born."
    << "These are finite contibutions already handled in LO merging."
    << Exception::abortnow;
  }
    // check for IR Safe Cutoff
  if( !currentNode()->allAbove( theIRSafePT ) )return ZERO;
  
  auto inOutPair =  currentNode()->getInOut();
  NodePtr randomChild =  currentNode()->randomChild();
  bool meRegion =matrixElementRegion(
                            inOutPair.first ,
                            inOutPair.second ,
                            randomChild->pT() ,
                            theMergePt );
  
  
  if ( meRegion )return MergingDSigDRRealAllAbove(  );
  if ( UseRandom::rndbool() )
  return 2.*MergingDSigDRRealBelowSubReal(  );
  return 2.*MergingDSigDRRealBelowSubInt(  );
}

CrossSection Merger::MergingDSigDRRealAllAbove( ){
    //If all dipoles pts are above , we cluster to this dipole.
  NodePtr CLNode =  currentNode()->randomChild();
    // Check if phase space poing is in ME region--> else rm PSP
  if ( !CLNode->children().empty() ) {
    auto inOutPair =  CLNode->getInOut();
    NodePtr randomChild =  CLNode->randomChild();
    if( !matrixElementRegion( inOutPair.first ,
                              inOutPair.second ,
                              randomChild->pT() ,
                              theMergePt ) )return ZERO;
  }
  
    // first find the history for the acctual Node
  NodePtr productionNode =  currentNode()-> getHistory( true , DSH()->hardScaleFactor() );
    // If CLNode is not part of the history , dont calculate the Real contribution
    // else multiply the real contribution with N ( number of children ).
    // this returns the sudakov suppression according to
    // the clustering of the born parts.
  bool inhist = CLNode->isInHistoryOf( productionNode );
  if(productionNode== currentNode())assert(!inhist);
    // get the history for the clustered Node.
  productionNode = CLNode-> getHistory( false , DSH()->hardScaleFactor() );
    // decide if to cluster
  weight = decideClustering(productionNode,CLNode,projected);
    // Check for cuts on the production process.
  if ( !productionNode->xcomb()->willPassCuts() )return ZERO;
    // calculate the staring scale
  Energy startscale = CKKW_StartScale( productionNode );
    // fill history with caluclation of sudakov supression
  fillHistory( startscale , productionNode , CLNode );
    // fill the projector -> return the scale of the last splitting
  currentNode()->runningPt( fillProjector( projected ? 2 : 1 ) );
    // the weight has three components to get shower history weight
  weight *= history.back().weight*alphaReweight(true)*pdfReweight();
  if( weight == 0. )return ZERO;
    // The inhist flag produces the correct cluster density.
  CrossSection me = ( inhist?TreedSigDR( startscale ):ZERO );
    // calculate the dipole
  CrossSection dip = CLNode->calcDip( startscale* currentME()->renFac() );
  
  CrossSection res =  	weight*as( startscale*DSH()->renFac() )/SM().alphaS()*
 			currentNode()->children().size()*
  		       	( me - dip );
    // [ DEBUG ]
  if ( currentNode()->legsize() == 6&&Debug::level > 2 )
     debugReal("RAA",weight,me,dip);
  
  return res;
}

CrossSection Merger::MergingDSigDRRealBelowSubReal(  ){
   // Choose a random child to produce history from. 
  NodePtr HistNode = currentNode()->randomChild();
   // Check that this subleading point is in ME region.
  if ( !HistNode->children().empty() ) {
    auto inOutPair = HistNode->getInOut();
    NodePtr randomChild =  HistNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inOutPair.first , inOutPair.second , randomChild->pT() , theMergePt ) )return ZERO;
  }
  
   // As the HistNode is now in ME region, we can cluster according to LO merging.
   // In this clustering we do not require a orering for the last step to the currentNode.
  const NodePtr productionNode = HistNode-> getHistory( false , DSH()->hardScaleFactor() );
   // If the real emission contribution should be unitarised, we decide here if we cluster. 
   // This applies to NLO corrections for the first emission w.r.t. the production process the first time.
  weight = decideClustering(productionNode,HistNode,projected);
   // The production node needs to fulfill the cut criterion of the production process.
  if ( !productionNode->xcomb()->willPassCuts() && !currentNode()->xcomb()->willPassCuts() )return ZERO;
   // Calculate the starting scale w.r.t. the production process.
  Energy startscale = CKKW_StartScale( productionNode );
   // If the production process does not fullfill the cut criterion use the real emission point (DEBUG)
  if (!productionNode->xcomb()->willPassCuts()) startscale = CKKW_StartScale( currentNode() );
   // DEBUG trial
  //currentNode()->xcomb()->lastProjector( productionNode->xcomb());
   // Calculate the sudakov weights starting from the production node to the histNode
  fillHistory( startscale , productionNode , HistNode );
   // Set the running Pt of the process as it is used in the vetoed parton shower.
  currentNode()->runningPt( fillProjector( projected ? 1 : 0 ) );
  currentNode()->runningPt(max(HistNode->pT(),theMergePt));
   // Calculate the alpha_S ratios and pdf ratios.
   weight *= history.back().weight*alphaReweight(true)*pdfReweight();
  
  if( weight == 0. )return ZERO;

   // Start calculation of subtraction contribution. 
  CrossSection sumPS = ZERO;

   // Iterate over all subtraction contributions.
  for( auto const & child : currentNode()->children() ){
    if ( child->allAbove( mergePt() ) && child->xcomb()->willPassCuts() ){
      Energy relevantScale=child->children().empty()?CKKW_StartScale( child ):child->maxChildPt();
      if( ( child )->pT()>mergePt()/theRealSubtractionRatio ){
	if( child ->pT()<relevantScale&&(child->inShowerPS(relevantScale))){  //DEBUG: CKKW_StartScale( child);????
          sumPS += child->calcPs( startscale* currentME()->renFac() );
	}
      }else{
	if( child ->pT()<relevantScale){
    	  sumPS += child->calcDip( startscale* currentME()->renFac() );
        }
      }
     }
  }
  
  CrossSection me = ZERO;
  if(currentNode()->xcomb()->willPassCuts()){
    me = TreedSigDR( startscale );
  } 


   // [ DEBUG ]
  if ( currentNode()->legsize() == 6&&Debug::level > 2 )
      debugReal("RBSR",weight,me,sumPS);
   //Here we subtract the PS ( and below the dynamical cutoff the Dip )  
  return weight*as( startscale*DSH()->renFac() )/SM().alphaS()*
  ( me-sumPS );
}



CrossSection Merger::MergingDSigDRRealBelowSubInt( ){
  
  if( currentNode()->children().empty() )return ZERO;
  NodePtr CLNode =  currentNode()->randomChild();
  if( CLNode->pT()<mergePt()/theRealSubtractionRatio )return ZERO;
  
  if ( !CLNode->children().empty() ) {
    auto inOutPair = CLNode->getInOut( );
    NodePtr randomChild =  CLNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if( !matrixElementRegion( inOutPair.first , inOutPair.second , randomChild->pT() , theMergePt ) )return ZERO;
  }
  
  
  const NodePtr productionNode = CLNode-> getHistory( false , DSH()->hardScaleFactor() );
  
  weight = decideClustering(productionNode,CLNode,projected);

  if ( !CLNode->allAbove( mergePt() ) )return ZERO;
  
  if ( !productionNode->xcomb()->willPassCuts() )return ZERO;
  
  Energy startscale = CKKW_StartScale( productionNode );
  
  fillHistory( startscale , productionNode , CLNode );
  
  currentNode()->runningPt( fillProjector( projected ? 2 : 1 ) );
  
  weight *= history.back().weight*alphaReweight(true)*pdfReweight();

  if( weight == 0. )return ZERO;
    
  pair<CrossSection , CrossSection> DipAndPs =
  CLNode->calcDipandPS( startscale* currentME()->renFac() );
    // [ DEBUG ]
  if ( currentNode()->legsize() == 6&&Debug::level > 2 )
     debugReal("RBSI",weight,DipAndPs.second,DipAndPs.first);
  
  Energy relevantScale=CLNode->children().empty()?CKKW_StartScale( CLNode ):CLNode->maxChildPt();
  bool calcPScontribution=CLNode->pT()<relevantScale&&(CLNode->inShowerPS(relevantScale));
    //Here we add the PS and subtrac the Dip ( only above the dynamical cutoff )
  return weight*as( startscale*DSH()->renFac() )/SM().alphaS()*
  currentNode()->children().size()*( (calcPScontribution?DipAndPs.second:ZERO)-DipAndPs.first );
}


CrossSection Merger::MergingDSigDRBornGamma( ){
  
  double weightCL = 0.;
  weight = 1.;
  
  if ( !currentNode()->children().empty() ) {
    auto const inOutPair = currentNode()->getInOut();
      // Here we make sure that clustering and splitting are invertible.
    NodePtr randomChild =  currentNode()->randomChild();
      // Check if point is part of the ME region.
    if( !matrixElementRegion( inOutPair.first ,
                              inOutPair.second ,
                              randomChild->pT() ,
                              theMergePt ) )weight *= 0.;
  }
  
  const NodePtr productionNode = currentNode()->getHistory( true , DSH()->hardScaleFactor() );
  NodePtr CLNode;
  NodePtr BornCL;
  
  
  if( !currentNode()->children().empty() ){
    if ( UseRandom::rndbool() ){
      CLNode =  currentNode()->randomChild();
      bool inhist = CLNode->isInHistoryOf( productionNode );
      weight *= inhist?( -2.*int( currentNode()->children().size() ) ):0.;
      projected = true;
      weightCL = 2.*int( currentNode()->children().size() );
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
       treefactory()->onlymulti() !=
         int( currentNode()->legsize()-(projected ? 1 : 0) ) )
    return ZERO;
  
  
  if( !currentNode()->allAbove( mergePt()*(1.-1e-6)  ) )weight = 0.;
  if( CLNode&&!CLNode->allAbove( mergePt()*(1.-1e-6) ) )weightCL = 0.;
  if ( !productionNode->xcomb()->willPassCuts() ){
    return ZERO;
  }
  
  CrossSection res = ZERO;
  bool maxMulti = currentNode()->legsize() == int( maxLegsLO() );
  
  
  if( weight != 0. ){
    Energy startscale = CKKW_StartScale( productionNode );
    fillHistory( startscale , productionNode , currentNode() );
    currentNode()->runningPt( fillProjector( (projected ? 1 : 0) ) );
    weight *= history.back().weight*alphaReweight()*pdfReweight();
    if( weight == 0.&&weightCL == 0. )return ZERO;
    
    res += weight*TreedSigDR( startscale , ( !maxMulti&&!projected )?theGamma:1. );
  }
  
  if( CLNode&&theGamma != 1. ){
    Energy startscale = CKKW_StartScale( BornCL );
    fillHistory( startscale , BornCL , CLNode );
    currentNode()->runningPt( fillProjector( projected ? 1 : 0 ) );
    weightCL *= history.back().weight*alphaReweight()*pdfReweight();
    CrossSection diff = ZERO;
     currentME()->factory()->setAlphaParameter( 1. );
    diff -=  weightCL*CLNode->dipole()->dSigHatDR( sqr( startscale* currentME()->renFac() ) );
     currentME()->factory()->setAlphaParameter( theGamma );
    
    string alp = ( CLNode->dipole()->aboveAlpha()?"Above":"Below" );
    
    diff += weightCL*CLNode->dipole()->dSigHatDR( sqr( startscale* currentME()->renFac() ) );
     currentME()->factory()->setAlphaParameter( 1. );
    
    res += diff;
  }
  return res;
}




CrossSection Merger::TreedSigDR( Energy startscale , double diffAlpha ){
  currentME()->setScale( sqr( startscale ) , sqr( startscale ) );
  CrossSection res =  currentME()->dSigHatDRB();
  /*bool useDipolesForME=false;
  if (useDipolesForME && !currentNode()->children().empty()){
	res=ZERO;
	for (auto const & child : currentNode()->children() )
	    res-=child->dipole()->dSigHatDR(sqr( startscale ));	
  }
  */
  if ( projected && emitDipoleMEDiff ) {
    CrossSection resDip=ZERO;
    for (auto const & child : currentNode()->children() )
      resDip-=child->dipole()->dSigHatDR(sqr( startscale ));
    setEmissionProbability(1.-min(resDip/res,res/resDip));
  }else{
    setEmissionProbability(0.);	
  }

  if ( diffAlpha != 1. ) {
    res +=  currentME()->dSigHatDRAlphaDiff( diffAlpha );
  }
  if( std::isnan( double( res/nanobarn ) ) ){
    generator()->logWarning(Exception()
                            << "Merger: TreedSigDR is nan"
                            << Exception::warning);
    res = ZERO;};
  return res;
}

CrossSection Merger::LoopdSigDR( Energy startscale ){
   currentME()->setScale( sqr( startscale ) , sqr( startscale ) );
   currentME()->doOneLoopNoBorn();
   CrossSection res =  currentME()->dSigHatDRV()+
                       currentME()->dSigHatDRI();
   currentME()->noOneLoopNoBorn();
  return res;
}

Energy Merger::fillProjector( int pjs ){
    // in the shower handler the scale is multiplied
    // by DSH()->hardScaleFactor() so here we need
    // to devide by the factor.
  double xiQSh = history.begin()->node->legsize() == N0()?DSH()->hardScaleFactor():1.;
  if( pjs == 0 ){
     return ( history.size() == 1?1.:( 1./xiQSh ) )*history.back().scale;
  }
  for( auto const & hs : history )
  if ( isProjectorStage( hs.node , pjs )&&pjs !=  0 ){
    currentNode()->xcomb()->lastProjector( hs.node->xcomb() );
    return ( hs.node == history[0].node?1.:( 1./xiQSh ) )*hs.scale;
  }
  
  throw Exception() << "Could not fill projector." << Exception::abortnow;
  return ZERO;
}

double Merger::pdfReweight(){ // TODO factorization scale inside
  
  double res = 1.;
  // consider both sides.
  for( int side : {0 , 1} ){
    // The side scale defines the scale that the leg is changig thou the history.
    // We start reweighting at the seed process.
    // We only need to calculate the pdf if the emission whould change the leg,
    // otherwise the leg remains the same.
    // If no emission is prduced from this leg the pdf ratios from this leg are always 1.
    Energy sidescale=history[0].scale*(
		       history[0].node->legsize() == N0()  ?
                       currentME()->facFac():
                       DSH()->facFac());
    bool sidechanged=false;
    // only if the incoming parton is coloured.
    if( history[0].node->xcomb()->mePartonData()[side]->coloured() ){
      // go though the history.
      for ( auto const & hs : history ){
          //pdf-ratio only to the last step
        if ( !hs.node->parent() ) continue;
        if ( hs.node == history.back().node ) continue;
        if ( !hs.node->dipole() ){
          throw Exception()
          << "\nMerger: pdfReweight: history step has no dipol. "
          << Exception::abortnow;
          return 0.;
        }
	// if the emitter is the side the emission changes the momentum fraction.
	if(!(hs.node->dipole()->bornEmitter()==side||
	// II dipoles dont change the momentum fraction of the spectator only FI
	   ( hs.node->dipole()->bornSpectator()==side &&  hs.node->dipole()->bornEmitter()>1)))
		continue;
        const bool fromIsME = false;
        const bool toIsME   = history[0].node == hs.node;
        res *= pdfratio( hs.node,
              		 //numerator
	                 sidescale,
			 // denominator
                         DSH()->facFac()*hs.node->pT(),
                         side,
                         fromIsME,
                         toIsME);
	sidescale= DSH()->facFac()*hs.node->pT();
        sidechanged=true;
      }

      const bool fromIsME = true;
      const bool toIsME   = !sidechanged && history[0].node->legsize() == N0();

      res *= pdfratio( history.back().node,
                       sidescale,
                       history[0].scale * currentME()->facFac(),
                       side,
                       fromIsME,
                       toIsME );
    }
  }
  if ( std::isnan( res ) )
     generator()->logWarning(Exception()
                          << "Merger: pdfReweight is nan."
                          << Exception::warning);

  return res;
}


double Merger::cmwAlphaS(Energy q)const{
  using Constants::pi;
    // No cmw-scheme
  if (theCMWScheme==0) 
    return as( q );
    // Linear cmw-scheme
  else if(theCMWScheme==1){
    double als=as( q );
    return als *
           (1.+(3.*(67./18.-1./6.*sqr(pi))
                 -5./9.*Nf(q))* als/2./pi);
  }
    // cmw-scheme as factor in argument.
  else if(theCMWScheme==2){
    double kg=exp(-(67.-3.*sqr(pi)-10/3*Nf(q))
                  /(     2.     *(33.-2.*Nf(q))));
    //Note factor 2 since we here dealing with Energy
    return as(max(kg*q,1_GeV));
  }else{
    throw Exception()
    << "This CMW-Scheme is not implemented."
    << Exception::abortnow;
  }
  return -1;
}




double Merger::alphaReweight(bool nocmw){
  double res = 1.;
  
  Energy Q_R = ( history[0].node->legsize() == N0()?
                 currentME()->renFac():
                DSH()->renFac() )*
                history[0].scale;
  
  using Constants::pi;
  const auto Q_qed=history[0].node->nodeME()->factory()->scaleChoice()->renormalizationScaleQED();
  const auto Oew=history[0].node->nodeME()->orderInAlphaEW();
  const auto Oqcd=history[0].node->nodeME()->orderInAlphaS();

  if (!history[0].node->children().empty()) {
    assert(Oqcd!=0);
  }
  
  res *= pow( SM().alphaEMME( Q_qed )/ SM().alphaEMMZ() , Oew );
  
  res *= pow( (nocmw?as(Q_R):cmwAlphaS(Q_R)) / SM().alphaS() , Oqcd );

  
  
  
  for ( auto const & hs : history )
  if ( hs.node!= history.back().node ){
    Energy q_i = DSH()->renFac()* hs.node->pT();
    res *= cmwAlphaS(q_i) / SM().alphaS();
  }
  
  if ( std::isnan( res ) )
      generator()->logWarning(Exception()
               << "Merger: alphaReweight is nan. "<< Exception::warning);
  return res;
}

void Merger::fillHistory( Energy scale , NodePtr begin , NodePtr endNode ){
  
  history.clear();
  double sudakov0_n = 1.;
  
  history.push_back( {begin , sudakov0_n , scale} );
  
  
  double xiQSh = history.begin()->node->legsize() == N0()?
                 DSH()->hardScaleFactor():1.;
  
  scale *= xiQSh;
  if ( begin->parent()||!isUnitarized ) {
    while ( begin->parent() && ( begin !=  endNode ) ) {
      if ( !dosudakov( begin , scale , begin->pT() , sudakov0_n ) ){
        history.push_back( { begin->parent() , 0. , scale } );
      }
      
      if ( std::isnan( sudakov0_n ) )
             generator()->logWarning(Exception()
                  << "Merger: sudakov"<<scale/GeV<<" "
                  <<begin->pT()/GeV<<"0_n is nan. "
                  << Exception::warning);
      
      scale = begin->pT();
      history.push_back( { begin->parent() , sudakov0_n , begin->pT() } );
      begin = begin->parent();
    }
    
    Energy notunirunning = scale;
    
    if ( !isUnitarized&&N()+N0() > int( currentNode()->legsize() ) ) {
      if ( !dosudakov( begin , notunirunning , mergePt() , sudakov0_n ) ){
        history.back().weight = 0.;
      }else{
        history.back().weight = sudakov0_n;
      }
    }
  }
  if( history.size() == 1 )scale /= DSH()->hardScaleFactor();
}




double Merger::sumPdfReweightExpansion()const{
  double res = 0.;
  Energy beam1Scale = history[0].scale*
                    ( history[0].node->legsize() == N0()?
                      currentME()->facFac():
                      DSH()->facFac() );
  Energy beam2Scale = history[0].scale*
                    ( history[0].node->legsize() == N0()?
                      currentME()->facFac():
                      DSH()->facFac() );
  
  for ( auto const & hs : history ){
      //pdf expansion only to the last step
    if( !hs.node->parent() )continue;
    if( hs.node->xcomb()->mePartonData()[0]->coloured()&&
        hs.node->nodeME()->lastX1()>0.99 )return 0.;
    if( hs.node->xcomb()->mePartonData()[1]->coloured()&&
        hs.node->nodeME()->lastX2()>0.99 )return 0.;
    
    if( hs.node->nodeME()->lastX1()<0.00001 )return 0.;
    if( hs.node->nodeME()->lastX2()<0.00001 )return 0.;
    
    if ( hs.node->dipole()->bornEmitter() == 0 ){
      res += pdfExpansion( hs.node , 0 ,
                       beam1Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX1() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam1Scale = ( hs.node->pT() )*DSH()->facFac();
    }
    else if ( hs.node->dipole()->bornEmitter() == 1 ){
      res += pdfExpansion( hs.node , 1 ,
                       beam2Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX2() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam2Scale = ( hs.node->pT() )*DSH()->facFac();
    }
    // if we're here we know hs.node->dipole()->bornEmitter() > 1
    // works only in collinear scheme
    else if ( hs.node->dipole()->bornSpectator() == 0 ){
      res += pdfExpansion( hs.node , 0 ,
                       beam1Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX1() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam1Scale = ( hs.node->pT() )*DSH()->facFac();
    }
    else if ( hs.node->dipole()->bornSpectator() == 1 ){
      res += pdfExpansion( hs.node , 1 ,
                       beam2Scale , 
                       ( hs.node->pT() ) , 
                       hs.node->nodeME()->lastX2() , 
                       Nf( history[0].scale ) , 
                       history[0].scale );
      beam2Scale = ( hs.node->pT() )*DSH()->facFac();
    }
  }
  
  if ( currentNode()->xcomb()->mePartonData()[0]->coloured() ){
    res += pdfExpansion( history.back().node , 0 ,
                     beam1Scale , 
                     history[0].scale* currentME()->facFac() , 
                     ( history.back() ).node->nodeME()->lastX1() , 
                     Nf( history[0].scale ) , 
                     history[0].scale );
    
  }
  if ( currentNode()->xcomb()->mePartonData()[1]->coloured() ) {
    res += pdfExpansion( history.back().node , 1 ,
                     beam2Scale , 
                     history[0].scale* currentME()->facFac() , 
                     ( history.back() ).node->nodeME()->lastX2() , 
                     Nf( history[0].scale ) , 
                     history[0].scale );
  }
  return res;
}





#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"
double Merger::pdfExpansion( NodePtr node  ,
                             int side , Energy  running ,
                             Energy next , double x ,
                             int nlp , Energy fixedScale ) const {
  
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
               match( inverse( 0.0 , x , 1.0 ) +
                      inverse( 1.0+eps , x , 1.0 ) ) ) , r );
    
    double z = zw.first;
    double mapz = zw.second;
    double PDFxparton = pdf->xfx( particle , parton , sqr( fixedScale ) , x )/x;
    static double CA = SM().Nc();
    static double CF = ( SM().Nc()*SM().Nc()-1.0 )/( 2.*SM().Nc() );
    
    
    if ( abs( parton->id() ) < 7 ) {
      
      double PDFxByzgluon = pdf->xfx( particle ,
                                     getParticleData( ParticleID::g ) ,
                                     sqr( fixedScale ) , x/z )*z/x;
      double PDFxByzparton = pdf->xfx( particle , parton ,
                                      sqr( fixedScale ) , x/z )*z/x;
      assert( abs( parton->id() ) < 7 );
      
      restmp += CF*( 3./2.+2.*log( 1.-x ) ) * PDFxparton;
      if ( z > x ) {
        restmp += 0.5 * ( sqr( z ) + sqr( 1.-z ) ) * PDFxByzgluon / z;
        restmp += CF*2.*( PDFxByzparton - z*PDFxparton )/( z*( 1.-z ) );
        restmp -=  CF*PDFxByzparton * ( 1.+z )/z;
      }
    }else{
      
      assert( parton->id() == ParticleID::g );
      double PDFxByzgluon = pdf->xfx( particle ,
                                      getParticleData( ParticleID::g ) ,
                                      sqr( fixedScale ) , x/z )*z/x;
      // Pqg
      if ( z > x ){
        double factor = CF * ( 1. + sqr( 1.-z ) ) / sqr( z );
        for ( int f = -nlp; f <=  nlp; ++f ) {
          if ( f == 0 )
          continue;
          restmp += pdf->xfx( particle , getParticleData( f ) ,
                              sqr( fixedScale ) , x/z )*z/x*factor;
        }
      }
      // Pgg 
      restmp += ( ( 11./6. ) * CA
                  - ( 1./3. )*Nf( history[0].scale )
                  + 2.*CA*log( 1.-x ) ) *PDFxparton;
      if ( z > x ) {
        restmp += 2. * CA * ( PDFxByzgluon - z*PDFxparton ) / ( z*( 1.-z ) );
        restmp += 2.* CA *( ( 1.-z )/z - 1. + z*( 1.-z ) ) * PDFxByzgluon / z;
      }
      
    }
    if ( PDFxparton<1e-8 )
	restmp =  0.;
    else
    	res += 1*restmp*log( sqr( running/next ) )/PDFxparton*mapz;
    
  }
  return res/number;
}



double Merger::sumAlphaSReweightExpansion()const{
  double res = 0.;
  const auto Oqcd=history[0].node->nodeME()->orderInAlphaS();
  res += alphasExpansion( history[0].scale* DSH()->renFac() ,
                            history[0].scale* currentME()->renFac() )*
                            Oqcd;

    // dsig is calculated with fixed alpha_s
  for ( auto const & hs : history ){
      //expansion only to the last step
    if( !hs.node->parent() )continue;
    res += alphasExpansion( hs.node->pT()*DSH()->renFac()  ,currentME()->renFac()*history[0].scale );
  }
  return res;
}

double Merger::sumFillHistoryExpansion(){
  double res = 0.;
  double xiQSh = history[0].node->legsize() == N0()?DSH()->hardScaleFactor():1.;
  for ( auto const & hs : history ){
    if( !hs.node->parent() )continue;
    doHistExpansion( hs.node , ( hs.node == history[0].node?xiQSh:1. )*hs.scale ,
                     hs.node->pT() , history[0].scale , res );
  }
  return res;
}


MergingFactoryPtr Merger::treefactory()const{return theTreeFactory;}


void Merger::doinit(){
  if ( !DSH()->hardScaleIsMuF() ) {
    throw Exception()
    << "Merger: Merging is currently only sensible "
    << "if we are using the hardScale as MuF."
    << Exception::abortnow;
  }
  DSH()->init();
  
  if ( !( DSH()->showerPhaseSpaceOption() == 0 ||
          DSH()->showerPhaseSpaceOption() == 1) ){
    throw InitException() << "Merger::doinit(): Choice of shower phase space cannot be handled by the merging";
  }
}


namespace{
void setXCombScales(StdXCombPtr xc,Energy2 scale){
    xc->lastShowerScale                       ( scale );
    xc->partonBinInstances().first->scale     ( scale );
    xc->partonBinInstances().second->scale    ( scale );
}
}

CrossSection Merger::MergingDSigDR() {
  
  history.clear();
  assert(currentNode()==theFirstNodeMap[ currentME()]);

  if(DSH()->doesSplitHardProcess()){
	throw Exception()
  	<< "Merger: The splithardprocess option is currently not supported."
  	<< Exception::abortnow;	
  }


      //get the PDF's (from ShowerHandler.cc)
  if (!DSH()->getPDFA()||!DSH()->firstPDF().particle()){
     tSubProPtr sub = currentNode()->xcomb()->construct();
     const auto pb=currentNode()->xcomb()->partonBins();
     tcPDFPtr first  = DSH()->getPDFA() ?
          tcPDFPtr(DSH()->getPDFA()) : DSH()->firstPDF().pdf();
     tcPDFPtr second = DSH()->getPDFB() ?
          tcPDFPtr(DSH()->getPDFB()) : DSH()->secondPDF().pdf();
     DSH()->resetPDFs( {first,second },pb );
  }

 
 
  DSH()->eventHandler( generator()->eventHandler() );
  
  CrossSection res = ZERO;
  if( currentNode()->subtractedReal() ){
    res = MergingDSigDRRealStandard();
    theCurrentMaxLegs = maxLegsNLO();
  }else if( currentNode()->virtualContribution() ){
    res = MergingDSigDRVirtualStandard();
    theCurrentMaxLegs = maxLegsNLO();
  }else if( theGamma!= 1. ){
    res = MergingDSigDRBornGamma();
    theCurrentMaxLegs = maxLegsLO();
  }else{
    res = MergingDSigDRBornStandard();
    theCurrentMaxLegs = maxLegsLO();
  }
  
  auto lxc= currentME()->lastXCombPtr();
  setXCombScales(lxc,sqr( currentNode()->runningPt()));


  auto lp= currentME()->lastXCombPtr()->lastProjector();
  
  if( lp )
    setXCombScales( lp, sqr( currentNode()->runningPt()));
  
  if ( res == ZERO ){
    history.clear();
    return ZERO;
  }
  
  
  cleanup( currentNode() );
  lxc->subProcess( SubProPtr() );
  
  history.clear();
  
  if( !std::isfinite( double( res/nanobarn ) ) ){
    generator()->logWarning(Exception()
       << "Merger weight is " << res/nanobarn<< " -> setting to 0"
       << Exception::warning);
    return ZERO;
  }
  
  return res;
  
}




#include "Herwig/PDF/HwRemDecayer.h"
void Merger::CKKW_PrepareSudakov( NodePtr node , Energy running ){
  tSubProPtr sub = node->xcomb()->construct();
  const auto pb=node->xcomb()->partonBins();
  DSH()->setCurrentHandler();
  DSH()->currentHandler()->generator()->currentEventHandler(currentNode()->xcomb()->eventHandlerPtr() );
  DSH()->currentHandler()->remnantDecayer()->setHadronContent( currentNode()->xcomb()->lastParticles() );
  DSH()->eventRecord().clear();
  DSH()->eventRecord().slimprepare( sub , dynamic_ptr_cast<tStdXCombPtr>( node->xcomb() ) , DSH()->pdfs() ,
				    currentNode()->xcomb()->lastParticles(), DSH()->offShellPartons() );
  DSH()->hardScales( sqr( running ) );
}


Energy Merger::CKKW_StartScale( NodePtr node ) const {
  Energy res = generator()->maximumCMEnergy();
  const int N=node->legsize();
  const auto & data=node->nodeME()->mePartonData();
  const auto & momenta=node->nodeME()->lastMEMomenta();
  if( !node->children().empty() ){
    for ( int i = 0; i<N ; i++ ){ if ( !data[i]->coloured() )continue;
      for ( int j = 2; j<N ; j++ ){ if ( i == j||!data[j]->coloured() )continue;
	for ( int k = 0; k<N ; k++ ){ if ( !data[k]->coloured() || i == k || j == k )continue;
	  if(i<2){
	    if(k<2) {
		res =  min( res , IILTK->lastPt( momenta[i] , momenta[j] , momenta[k] ));
	    }
	    else    {

                res =  min( res ,
                (data[k]->mass()+data[j]->mass()+data[i]->mass()>ZERO)?
                IFMTK->lastPt( momenta[i] , momenta[j] , momenta[k] ):
                IFLTK->lastPt( momenta[i] , momenta[j] , momenta[k] ));
            }
	  }else{
	    if(k<2) {
                res =  min( res ,
                (data[k]->mass()+data[j]->mass()+data[i]->mass()>ZERO)?
                FIMTK->lastPt( momenta[i] , momenta[j] , momenta[k] ):
                FILTK->lastPt( momenta[i] , momenta[j] , momenta[k] ));

	    }
	    else    {
		res =  min( res , 
		(data[k]->mass()+data[j]->mass()+data[i]->mass()>ZERO)?
		FFMTK->lastPt( momenta[i] , momenta[j] , momenta[k] ):
		FFLTK->lastPt( momenta[i] , momenta[j] , momenta[k] ));
	    }
	  }
	}
      }
    }
  }else{
    node->nodeME()->factory()->scaleChoice()->setXComb( node->xcomb() );
    res =  sqrt( node->nodeME()->factory()->scaleChoice()->renormalizationScale() );
  }
  node->nodeME()->factory()->scaleChoice()->setXComb( node->xcomb() );
  res = max( res , sqrt( node->nodeME()->factory()->scaleChoice()->renormalizationScale() ) );
  return res;
}




double Merger::alphasExpansion( Energy next , Energy fixedScale ) const {
  double betaZero = ( 11./6. )*SM().Nc() - ( 1./3. )*Nf( history[0].scale );
  double K=3.*( 67./18.-1./6.*sqr(Constants::pi) ) -5./9.*Nf( history[0].scale);
  return ( betaZero*log( sqr( fixedScale/next ) ) )+( theCMWScheme>0?K:0. );
}


double Merger::pdfratio( NodePtr  node ,
                         Energy  numerator_scale ,
                         Energy denominator_scale ,
                         int side ,
                         bool fromIsME,
                         bool toIsME ){
  StdXCombPtr bXc = node->xcomb();
  if( !bXc->mePartonData()[side]->coloured() )
  throw Exception()
  << "Merger: pdf-ratio required for non-coloured particle."
  << Exception::abortnow;
  
  
  double from = 1.; 
  double to   = 1.;
  if ( side == 0 ){
    if ( denominator_scale == numerator_scale && fromIsME==toIsME ) {
      return 1.;
    }


    if (fromIsME) {
      from = node->nodeME()->pdf1( sqr( denominator_scale ) );
    }else{
      from = DSH()->firstPDF().xfx(node->xcomb()->lastPartons().first->dataPtr(),
                                  sqr( denominator_scale ),
                                  node->xcomb()->lastX1())/node->xcomb()->lastX1();
    }


    if (toIsME) {
      to = node->nodeME()->pdf1( sqr( numerator_scale ) );
    }else{
      to = DSH()->firstPDF().xfx(node->xcomb()->lastPartons().first->dataPtr(),
                                  sqr( numerator_scale ),
                                  node->xcomb()->lastX1())/node->xcomb()->lastX1();
    }
    
    if ( ( to < 1e-8||from < 1e-8 )&&( to/from>10000000. ) ){
      generator()->logWarning(Exception()
                              << "Merger: pdfratio to = " << to << " from = " << from
                              << Exception::warning);
      return 0.;
    }
  }
  else{
    if ( denominator_scale == numerator_scale && fromIsME==toIsME ) {
      return 1.;
    }
    
    if (fromIsME) {
      from = node->nodeME()->pdf2( sqr( denominator_scale ) );
    }else{
      from =DSH()->secondPDF().xfx(node->xcomb()->lastPartons().second->dataPtr(),
                                  sqr( denominator_scale ),
                                  node->xcomb()->lastX2())/node->xcomb()->lastX2();
    }


    if (toIsME) {
      to = node->nodeME()->pdf2( sqr( numerator_scale ) );
    }else{
      to = DSH()->secondPDF().xfx(node->xcomb()->lastPartons().second->dataPtr(),
                                sqr( numerator_scale ),
                                node->xcomb()->lastX2())/node->xcomb()->lastX2();
    }
    
    if ( ( to < 1e-8||from < 1e-8 )&&( to/from>10000000. ) ){
      generator()->logWarning(Exception()
                              << "Merger: pdfratio to = "
                              << to << " from = " << from
                              << Exception::warning);
      return 0.;}
  }
  return to/from;
}



bool Merger::dosudakov( NodePtr node , Energy running , Energy next , double& sudakov0_n ) {
  CKKW_PrepareSudakov( node , running );
  for( DipoleChain const & chain : DSH()->eventRecord().chains() ){
    for( Dipole const & dip : chain.dipoles() ){
      sudakov0_n *= singlesudakov( dip , next , running , { true , false } );
      sudakov0_n *= singlesudakov( dip , next , running , { false , true } );
      if ( sudakov0_n == 0.0 ){
          cleanup( node );
        return false;
      }
    }
  }
    cleanup( node );
  return true;
}

bool Merger::doHistExpansion( NodePtr node , Energy  running , Energy next , Energy fixedScale , double& histExpansion ) {
  CKKW_PrepareSudakov( node , running );
  for( DipoleChain const & chain : DSH()->eventRecord().chains() ){
    for( Dipole const & dip : chain.dipoles() ){
      histExpansion += singleHistExpansion( dip , next , running , fixedScale , { true , false } );;
      histExpansion += singleHistExpansion( dip , next , running , fixedScale , { false , true } );
    }
  }
  cleanup( node );
  return true;
}



bool Merger::isProjectorStage( NodePtr  node , int pjs )const{
  return ( pjs == int( ( currentNode()->legsize() - node->legsize() ) ) );
}

void Merger::cleanup( NodePtr node ) {
  DSH()->eventRecord().clear();
  if( !node->xcomb()->subProcess() )return;
  ParticleVector vecfirst = node->xcomb()->subProcess()->incoming().first->children();
  for( auto const & particle : vecfirst )
  node->xcomb()->subProcess()->incoming().first->abandonChild( particle );
  
  ParticleVector vecsecond = node->xcomb()->subProcess()->incoming().second->children();
  for( auto const & particle : vecsecond )
  node->xcomb()->subProcess()->incoming().second->abandonChild( particle );
  node->xcomb()->subProcess( SubProPtr() );
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

    Energy ptMax = gen->second->splittingKinematics()->ptMax( 
                  candidate.scale() , 
                  candidate.emitterX() , 
                  candidate.spectatorX() ,
                  candidate , 
                  *gen->second->splittingKernel() );
    
    candidate.hardPt( min( running , ptMax ) );
    
    if ( candidate.hardPt()>next ){
      res *= gen->second->sudakov( candidate , next );
    }
  }
  
  return res;
}


double Merger::singleHistExpansion( Dipole dip  , Energy next , Energy running , Energy fixedScale , pair<bool , bool> conf ){
  
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
                      candidate.scale() ,
                      candidate.emitterX() ,
                      candidate.spectatorX() ,
                      candidate ,
                      *gen->second->splittingKernel() );
    
    candidate.hardPt( min( running , ptMax ) );
    if ( candidate.hardPt()>next ){
      res += gen->second->sudakovExpansion( candidate , next , fixedScale );
    }
  }
	 
	 return res;
}




void Merger::firstNodeMap( MatchboxMEBasePtr a , NodePtr b ){
  theFirstNodeMap.insert( { a , b } );
}



map<MatchboxMEBasePtr , NodePtr> Merger::firstNodeMap()const{return theFirstNodeMap;}




void Merger::setXComb( tStdXCombPtr xc ){
  currentNode()->setXComb( xc );
}
void Merger::setKinematics(  ){
  currentNode()->setKinematics();
}
void Merger::clearKinematics(  ){
  currentNode()->clearKinematics();
}

void Merger::flushCaches(){
  if (currentNode()&&currentNode()->xcomb()->lastParticles().first) {
    currentNode()->flushCaches();
  }
}

bool Merger::generateKinematics( const double * r ){
  return currentNode()->firstgenerateKinematics( r , ! currentNode()->subtractedReal() );
}

bool Merger::matrixElementRegion( PVector incoming , PVector outgoing , Energy winnerScale , Energy cutscale )const{
  

  Energy ptx = Constants::MaxEnergy;
  bool foundwinnerpt = false;
  using namespace boost;
    //FF
  
  for( auto const & em : outgoing ){ if ( ! em->coloured() ) continue;
    for( auto const & emm : outgoing ){ if ( !emm->coloured() ) continue; if ( em == emm ) continue;
      for( auto const & spe : outgoing ){ if ( !spe->coloured() ) continue; if ( em == spe||emm == spe ) continue;
        
        if ( !( em->id() == -emm->id()||emm->id()>6 ) )continue;
        
        Energy pt = ZERO;
        if ( em->momentum().m()<= 10_MeV &&
            emm->momentum().m()<= 10_MeV &&
            spe->momentum().m()<= 10_MeV  ) {
          pt = FFLTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }else{
          pt = FFMTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }
        
        if ( abs( pt-winnerScale ) < 10_MeV  ) {
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
        if ( em->momentum().m()<= 10_MeV &&
            emm->momentum().m()<= 10_MeV &&
            spe->momentum().m()<= 10_MeV  ) {
          pt = FILTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }else{
          pt = FIMTK->lastPt( em->momentum() , emm->momentum() , spe->momentum() );
        }
        
        
        if ( abs( pt-winnerScale )<10_MeV  ) {
          foundwinnerpt = true;
        }
        
        if( pt > ZERO )
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
        
        if ( em->momentum().m()<= 10_MeV &&
             emm->momentum().m()<= 10_MeV &&
             spe->momentum().m()<= 10_MeV  ) {
            //massless
          pt = IFLTK->lastPt( em->momentum() , emm->momentum() , spe->momentum()  );
        }else{
            //massiv
          pt = IFMTK->lastPt( em->momentum() , emm->momentum() , spe->momentum()  );
        }
        
        
        if ( abs( pt-winnerScale )< 10_MeV ) {
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
        
        if ( abs( pt-winnerScale )< 10_MeV ) {
          foundwinnerpt = true;
        }
        ptx  = min( ptx , pt );
      }
    }
  }
  
  if( !foundwinnerpt ){
    generator()->logWarning( Exception()
                            << "Merger: Could not find winner with pt."
                            << "Run with -d3 to get phase space points. "
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

bool Merger::notOnlyMulti() const {
  return ( treefactory()->onlymulti() != -1&&
          treefactory()->onlymulti() !=
          int( currentNode()->legsize()-( projected ? 1 : 0 ) ) );
}

int Merger::M()const{return theTreeFactory->M();}

int Merger::N()const{return theTreeFactory->N();}

void Merger::debugVirt(double weight,double w1,double w2,double w3,
                       CrossSection matrixElement,double ww1,double ww2,
                       double ww3, NodePtr node,CrossSection bornWeight)const{
  Energy minPT = Constants::MaxEnergy;
  for( auto const & no :currentNode()->children() )minPT = min( no->pT() , minPT );
  
  generator()->log() << "\nVIRT " << minPT/GeV << " " << weight << " " << w1;
  generator()->log() << " " << w2;
  generator()->log() << " " << w3;
  generator()->log() << " " << ( matrixElement/nanobarn ) << " " << ww1 << " " << ww2 << " " << ww3 << " " << node->pT()/GeV << " " << node->nodeME()->mePartonData()[3]->mass()/GeV << " " << ( bornWeight*SM().alphaS()/( 2.*ThePEG::Constants::pi )/nanobarn );
}



void Merger::debugReal(string realstr, double weight,
                       CrossSection me, CrossSection dip)const {
  Energy minPT = Constants::MaxEnergy;
  for( auto const & no :currentNode()->children() )minPT = min( no->pT() , minPT );
  
  generator()->log() << "\n"<<realstr <<" "<< minPT/GeV << " " << weight << " " << ( me-dip )/nanobarn << " " << me/nanobarn << " " << dip/nanobarn;
  
}

  // If needed , insert default implementations of virtual function defined
  // in the InterfacedBase class here ( using ThePEG-interfaced-impl in Emacs ).

#include "ThePEG/Persistency/PersistentOStream.h"
void Merger::persistentOutput( PersistentOStream & os ) const {
  os << theShowerExpansionWeights << theCMWScheme << projected <<
  isUnitarized << isNLOUnitarized  <<
  theChooseHistory <<
  theN0 << theOnlyN   << weight <<
  theGamma << theSmearing <<  ounit( theIRSafePT , GeV ) <<
  ounit( theMergePt , GeV ) <<  ounit( theCentralMergePt , GeV )
  << theLargeNBasis << FFLTK << FILTK <<
  IFLTK << IILTK << FFMTK << FIMTK << IFMTK <<
  theDipoleShowerHandler << theTreeFactory << 
  theFirstNodeMap<<theRealSubtractionRatio << emitDipoleMEDiff;
  
}
#include "ThePEG/Persistency/PersistentIStream.h"
void Merger::persistentInput( PersistentIStream & is , int ) {
  is >> theShowerExpansionWeights >> theCMWScheme >> projected >>
  isUnitarized >> isNLOUnitarized >>
  theChooseHistory >> theN0 >> theOnlyN >>
  weight >> 
  theGamma >>   theSmearing >>  iunit( theIRSafePT , GeV ) >>
  iunit( theMergePt , GeV ) >>
  iunit( theCentralMergePt , GeV ) >> theLargeNBasis >>
  FFLTK >> FILTK >> IFLTK >>
  IILTK >> FFMTK >> FIMTK >>
  IFMTK >> theDipoleShowerHandler >>  theTreeFactory >> 
  theFirstNodeMap >> theRealSubtractionRatio >> emitDipoleMEDiff;
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
  ( "The Merger class takes care of merging multiple LO & NLO cross sections." );
  
 ////////////////////////////////////////////////// 
  
  static Switch<Merger , unsigned int>
  interfaceShowerExpansionWeights
  ( "ShowerExpansionWeights" ,
   "Calculate the expansions of the shower history to be NLO accurate, in different schemes." ,
   &Merger::theShowerExpansionWeights ,
   3 , false , false );
  
  static SwitchOption
  interfaceShowerExpansionWeightsScheme0
  ( interfaceShowerExpansionWeights ,
   "NoWeights" ,
   "Switch off the expansion." ,
   0 );
  static SwitchOption
  interfaceShowerExpansionWeightsScheme1
  ( interfaceShowerExpansionWeights ,
   "FlatAndHistoryReweighted" ,
   "Sum alphaS expansion and weight with same Historyweight as LO." ,
   1 );
  static SwitchOption
  interfaceShowerExpansionWeightsScheme2
  ( interfaceShowerExpansionWeights ,
   "NoAlphaSReweightForSudakovExpansion" ,
   "Switch off the expansion." ,
   2 );
  static SwitchOption
  interfaceShowerExpansionWeightsScheme3
  ( interfaceShowerExpansionWeights ,
   "NoAlphaSReweightForAllExpansions" ,
   "Switch off the expansion." ,
   3 );
  static SwitchOption
  interfaceShowerExpansionWeightsScheme4
  ( interfaceShowerExpansionWeights ,
   "NoAlphaSReweightForAlphaSExpansion" ,
   "Switch off the expansion." ,
   4 ); 

 /////////////////////////////////////////////////////////////////////
 
  static Switch<Merger , unsigned int>
  interfacetheCMWScheme
  ( "CMWScheme" ,
   "Use CMW-Scheme to calculate the alpha_s for the shower expressions." ,
   &Merger::theCMWScheme ,
   0 ,
   false , false );
  static SwitchOption interfacetheCMWSchemeNo
  (interfacetheCMWScheme,
   "No",
   "No CMW-Scheme",
   0);
  static SwitchOption interfacetheCMWSchemeLinear
  (interfacetheCMWScheme,
   "Linear",
   "Linear CMW multiplication: alpha_s(q) -> alpha_s(q)(1+K_g*alpha_s(q)/2pi )",
   1);
  static SwitchOption interfacetheCMWSchemeFactor
  (interfacetheCMWScheme,
   "Factor",
   "Use factor in alpha_s argument: alpha_s(q) -> alpha_s(fac*q) with fac=exp(-(67-3pi^2-10/3*Nf)/(33-2Nf)) ",
   2);
  
  static Parameter<Merger , Energy> interfaceMergerScale
  ( "MergingScale" ,  "The MergingScale." , 
   &Merger::theCentralMergePt ,
   GeV , 20.0*GeV , 0.0*GeV , 0*GeV ,
   false , false , Interface::lowerlim );
  
  static Parameter<Merger , double> interfacegamma
  ( "gamma" , 
   "gamma parameter." , 
   &Merger::theGamma , 1.0 , 0.0 , 0 , 
   false , false , Interface::lowerlim );
   interfacegamma.rank(-1); 
  
  static Switch<Merger , bool>
  interfaceemitDipoleMEDiff
  ( "emitDipoleMEDiff" , 
  "Allow emissions of the unitarisation contribution with prob. 1-min(B_n/sum Dip_n,sum Dip_n/B_n)" , 
  &Merger::emitDipoleMEDiff , false , false , false );
  static SwitchOption interfaceemitDipoleMEDiffYes
  ( interfaceemitDipoleMEDiff , "Yes" , "" , true );
  static SwitchOption interfaceemitDipoleMEDiffNo
  ( interfaceemitDipoleMEDiff , "No" , "" , false );
  interfaceemitDipoleMEDiff.rank(-1);

  static Parameter<Merger , Energy>  interfaceIRSafePT
  ( "IRSafePT" ,
  "Set the pt for which a matrixelement should be treated as IR-safe." ,
  &Merger::theIRSafePT ,
  GeV , 0.0 * GeV , ZERO , Constants::MaxEnergy ,
  true , false , Interface::limited );
  interfaceIRSafePT.setHasDefault( false );
  
  
  static Parameter<Merger , double>  interfacemergePtsmearing(
  "MergingScaleSmearing" ,
  "Set the percentage the merging pt should be smeared." ,
  &Merger::theSmearing ,
  0. , 0. ,  0.0 , 0.5 ,
  true , false , Interface::limited );
  
 
  static Parameter<Merger , int>  interfacechooseHistory(
  "chooseHistory" ,
  "Various   ways of choosing the history weights: 0(default):dipole xs, 1:dipole/born, 2:flat, 3: 1/pt_dip" ,
  &Merger::theChooseHistory ,
  3 , -1 , 0 ,
  false , false , Interface::lowerlim );
  
  
  static Parameter<Merger , double>  interfacetheRealSubtractionRatio(
  "RealSubtractionRatio" ,
  "If the pt of the Dipole divided by the merging scale is lower than the ratio, the dipole is used to subtract the real emission. Otherwise the shower approximation is used, " ,
  &Merger::theRealSubtractionRatio , 0. , 0. ,
  1. , 10.0 ,
  true , false , Interface::limited );

}
















