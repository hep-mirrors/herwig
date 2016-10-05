  // -*- C++ -*-
  //
  // Merger.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2007 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the Merger class.
  //

#include "Merger.h"
#include "Node.h"
#include "MFactory.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "Herwig/Shower/Dipole/DipoleShowerHandler.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/PDF/HwRemDecayer.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"


#include "fastjet/ClusterSequence.hh"


using namespace Herwig;



Merger::Merger()
: MergerBase() {
  
  Unlopsweights=true;
  theCMWScheme=true;
  theSmearing=0.;
  theDipoleShowerHandler=
  Ptr<DipoleShowerHandler>::ptr();
  FFLTK=new_ptr(FFLightTildeKinematics());
  FILTK=new_ptr(FILightTildeKinematics());
  IFLTK=new_ptr(IFLightTildeKinematics());
  IILTK=new_ptr(IILightTildeKinematics());
  FFMTK=new_ptr(FFMassiveTildeKinematics());
  FIMTK=new_ptr(FIMassiveTildeKinematics());
  IFMTK=new_ptr(IFMassiveTildeKinematics());
  isUnitarized=true;
  isNLOUnitarized=true;
  defMERegionByJetAlg=false;
  theOpenInitialStateZ=false;
  theChooseHistory=8;
  xiRenME=1.;
  xiFacME=1.;
  xiRenSh=1.;
  xiFacSh=1.;
  xiQSh=1.;
  theGamma=1.;
  ee_ycut=1000000;
  pp_dcut=1000000;
  therenormscale=0.*GeV;
  theIRSafePT=1.*GeV;
  theMergePt=3.*GeV;
  theCentralMergePt=3.*GeV;
}


Merger::~Merger() {}

IBPtr Merger::clone() const {
  return new_ptr(*this);
}

IBPtr Merger::fullclone() const {
  return new_ptr(*this);
}

pair<PVector,PVector> getInOut(NPtr Node){
  PVector incoming;
  for(auto const & i : {0,1})
    incoming.push_back(Node->nodeME()->mePartonData()[i]->produceParticle(Node->nodeME()->lastMEMomenta()[i]));
    PVector outgoing;
  for (size_t i=2;i<Node->nodeME()->mePartonData().size();i++){
    Ptr<ThePEG::Particle>::ptr p =Node->nodeME()->mePartonData()[i]->produceParticle(Node->nodeME()->lastMEMomenta()[i]);
    outgoing.push_back(p);
  }
  return make_pair(incoming,outgoing);
}


CrossSection Merger::MergingDSigDRBornGamma(NPtr Node){
  
  double weightCL=0.;
  weight=1.;


 
  if (!Node->children().empty()) {
    auto const inoutpair=getInOut(Node);
    NPtr rndCh= Node->randomChild(); // Here we make sure that clustering and splitting are invertible
    if(!matrixElementRegion(inoutpair.first,inoutpair.second, rndCh->pT(), theMergePt))weight*= 0.;
  }
  
  NPtr Born= Node-> getHistory(true,xiQSh);
  NPtr CLNode;
  NPtr BornCL;
  
  
  if( !Node->children().empty()){
    if (UseRandom::rnd()<.5){
      CLNode= Node->randomChild();
      bool inhist=CLNode->isInHistoryOf(Born);
      weight*=inhist?(-2.*int(Node->children().size())):0.;
      projected=true;Node->nodeME()->projectorStage(1);
      weightCL=2.*int(Node->children().size());
      BornCL=CLNode-> getHistory(false,xiQSh);
    }else{
      weight=2.;projected=false;Node->nodeME()->projectorStage(0);
    }
  }else{
    weight=1.;projected=false;Node->nodeME()->projectorStage(0);
  }

 
  if (treefactory()->onlymulti()!=-1&&
      treefactory()->onlymulti()!=int(Node->legsize()-Node->nodeME()->projectorStage()))
    return ZERO;
  
  
  if(!Node->allAbove(mergePt()-0.1*GeV))weight=0.;
  if(CLNode&&!CLNode->allAbove(mergePt()-0.1*GeV))weightCL=0.;
  if (!Born->xcomb()->willPassCuts()){
    
    if (Born->parent()) {
        //cout<<"\n parent not pass";
    }
    return ZERO;
    
  }
  
  CrossSection res=ZERO;
  bool maxMulti=Node->legsize() == int(maxLegsLO());
  
  
  if(weight!=0.){
    Energy startscale=CKKW_StartScale(Born);
    fillHistory( startscale,  Born, Node);
    Node->runningPt(fillProjector());
    weight*=history.back().weight*alphaReweight()*pdfReweight();
    if(weight==0.&&weightCL==0.)return ZERO;
    
    res+=weight*TreedSigDR(startscale,Node,(!maxMulti&&!projected)?theGamma:1.);
  }
  
  if(CLNode&&theGamma!=1.){
    Energy startscale=CKKW_StartScale(BornCL);
    fillHistory( startscale,  BornCL, CLNode);
    Node->runningPt(fillProjector());
    weightCL*=history.back().weight*alphaReweight()*pdfReweight();
    CrossSection diff=ZERO;
    Node->nodeME()->factory()->setAlphaParameter(1.);
    diff-=weightCL*CLNode->dipole()->dSigHatDR(sqr(startscale*xiFacME));
    Node->nodeME()->factory()->setAlphaParameter(theGamma);
    
    string alp=(CLNode->dipole()->aboveAlpha()?"Above":"Below");
    
    diff+=weightCL*CLNode->dipole()->dSigHatDR(sqr(startscale*xiFacME));
    Node->nodeME()->factory()->setAlphaParameter(1.);
    
    res+=diff;
  }
  
  
  return res;
}





CrossSection Merger::MergingDSigDRBornStandard(NPtr Node){
    // Check if phase space poing is in ME region
  if (!Node->children().empty()) {
    auto const & inoutpair=getInOut(Node);
    NPtr rndCh= Node->randomChild(); // Here we make sure that clustering and splitting are invertible
    if(!matrixElementRegion(inoutpair.first,inoutpair.second, rndCh->pT(), theMergePt))return ZERO;
  }
  
    // get the history for the process
  NPtr Born= Node-> getHistory(true,xiQSh);
  if( Born!= Node){
      // at least one history step -> unitarisation
    if (UseRandom::rnd()<.5){
      weight=-2.;projected=true;Node->nodeME()->projectorStage(1);
    }else{
      weight=2.;projected=false;Node->nodeME()->projectorStage(0);
    }
  }else{
      // no ordered history -> no projection
    weight=1.;projected=false;Node->nodeME()->projectorStage(0);
  }
  
  if (treefactory()->onlymulti()!=-1&&
      treefactory()->onlymulti()!=int(Node->legsize()-(projected?1:0)))
    return ZERO;
  
  assert(Node->allAbove(mergePt()-0.1*GeV));
  
  if (!Born->xcomb()->willPassCuts())return ZERO;
    // calculate the staring scale
  Energy startscale=CKKW_StartScale(Born);
    // fill history with caluclation of sudakov supression
  fillHistory( startscale,  Born, Node);
    // fill the projector, the argument gets set to current scale
  Node->runningPt(fillProjector());
  assert(Node->runningPt()!=ZERO);
    // the weight has three components to get shower history weight
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.)return ZERO;
    //calculate the cross section
  return weight*TreedSigDR(startscale,Node,1.);
}



CrossSection Merger::MergingDSigDRVirtualStandard(NPtr Node ){
  
    // Check if phase space poing is in ME region
  if (!Node->children().empty()) {
    auto inoutpair=getInOut(Node);
    NPtr rndCh= Node->randomChild(); // Here we make sure that clustering and splitting are invertible
    if(!matrixElementRegion(inoutpair.first,inoutpair.second, rndCh->pT(), theMergePt))return ZERO;
  }

  
  NPtr Born= Node-> getHistory(true,xiQSh);
  if( Born!= Node){
    if (UseRandom::rnd()<.5){
      weight=-2.;projected=true;Node->nodeME()->projectorStage(1);
    }else{
      weight=2.;projected=false;Node->nodeME()->projectorStage(0);
    }
  }else{
    weight=1.;projected=false;Node->nodeME()->projectorStage(0);
  }
  if (!Born->xcomb()->willPassCuts())return ZERO;
  Energy startscale=CKKW_StartScale(Born);
  fillHistory( startscale,  Born, Node);
  Node->runningPt(fillProjector());
  
  double ww1=history.back().weight;
  double ww2=alphaReweight();
  double ww3=pdfReweight();
  
  
  weight*=ww1*ww2*ww3;
  if(weight==0.)return ZERO;
  
  CrossSection matrixElement=LoopdSigDR(startscale,Node);
  
  CrossSection Bornweight=Node->nodeME()->dSigHatDRB();
  
  double w1=-sumpdfReweightUnlops();
  double w2=-sumalphaReweightUnlops();
  double w3=-sumfillHistoryUnlops();
  
  CrossSection unlopsweight =(w1+w2+w3)
                       *Bornweight
                       *SM().alphaS()/(2.*ThePEG::Constants::pi);
  
  if (Node->legsize()==5&&Debug::level > 2) {
    Energy minPT=1000*GeV;
    for(auto const & no :Node->children() )minPT=min(no->pT(),minPT);
    
    cout<<"\nVIRT "<<minPT/GeV<<" "<<weight<<" "<<w1;
    cout<<" "<<w2;
    cout<<" "<<w3;
    cout<<" "<<(matrixElement/nanobarn)<<" "<<ww1<<" "<<ww2<<" "<<ww3<<" "<<Born->pT()/GeV<<" "<<Born->nodeME()->mePartonData()[3]->mass()/GeV<<" "<<(Bornweight*SM().alphaS()/(2.*ThePEG::Constants::pi)/nanobarn);
  }
  
  
  return weight*
         as(startscale*xiRenSh)/SM().alphaS()*
         (matrixElement+unlopsweight);
}


CrossSection Merger::MergingDSigDRRealStandard(NPtr Node){
  bool allAbove=Node->allAbove(mergePt());
    //TODO: VW Abgas Skandal
  if(!Node->allAbove((Debug::level > 2?0.01:1.)*theIRSafePT))return ZERO;
  if (allAbove)return MergingDSigDRRealAllAbove(Node);
  if (UseRandom::rnd()<.5)
    return 2.*MergingDSigDRRealBelowSubReal( Node );
  return 2.*MergingDSigDRRealBelowSubInt( Node);
}

CrossSection Merger::MergingDSigDRRealAllAbove(NPtr Node){
  
  assert(!Node->children().empty());
    //If all dipoles pts are above, we cluster to this dipole.
  NPtr CLNode= Node->randomChild();
    // Check if phase space poing is in ME region--> else rm PSP
  if (!CLNode->children().empty()) {
    auto inoutpair=getInOut(CLNode);
    NPtr rndCh= CLNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if(!matrixElementRegion(inoutpair.first,inoutpair.second, rndCh->pT(), theMergePt))return ZERO;
  }

    // first find the history for the acctual Node
  NPtr Born= Node-> getHistory(true,xiQSh);
    // check if the randomly choosen CLNode is part of the history.
    // If CLNode is not part of the history, dont calculate the Real contribution
    // else multiply the real contribution with N (number of children).
    // this returns the sudakov suppression according to the clustering of the born parts.
  bool inhist=CLNode->isInHistoryOf(Born);
    // get the history for the clustered Node.
  Born=CLNode-> getHistory(false,xiQSh);
  if( Born!= CLNode){
    if (UseRandom::rnd()<.5){
      weight=-2.; projected=true;  Node->nodeME()->projectorStage(2);}
    else{
      weight= 2.; projected=false; Node->nodeME()->projectorStage(1);}
  }else{
    weight=1.;  projected=false; Node->nodeME()->projectorStage(1);
  }
  if (!CLNode->allAbove(mergePt()))return ZERO;
  if (!Born->xcomb()->willPassCuts())return ZERO;
  Energy startscale=CKKW_StartScale(Born);
  fillHistory( startscale,  Born, CLNode);
  Node->runningPt(fillProjector());
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.)return ZERO;
  
  CrossSection me=(inhist?TreedSigDR(startscale,Node):ZERO);
  CrossSection dip=CLNode->calcDip(startscale*xiFacME);
  
  
  CrossSection res= weight*as(startscale*xiRenSh)/SM().alphaS()*
  (double)Node->children().size()*(me -dip);
  if (Node->legsize()==6&&Debug::level > 2) {
  Energy minPT=1000*GeV;
  for(auto const & no :Node->children() )minPT=min(no->pT(),minPT);
  
  cout<<"\nRAA "<<minPT/GeV<<" "<<weight<<" "<<(me-dip)/nanobarn<< " "<<me/nanobarn<<" "<<dip/nanobarn;
  }
    // cout<<"\nCLNode->dipole()->dSigHatDR(sqr(startscale*xiFacME))/nanobarn "<<CLNode->dipole()->dSigHatDR(sqr(startscale*xiFacME))/nanobarn;
  return res;
}

CrossSection Merger::MergingDSigDRRealBelowSubReal(NPtr Node){
  NPtr HistNode=Node->randomChild();
  if (!HistNode->children().empty()) {
    auto inoutpair=getInOut(HistNode);
    NPtr rndCh= HistNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if(!matrixElementRegion(inoutpair.first,inoutpair.second, rndCh->pT(), theMergePt))return ZERO;
  }
  
  NPtr Born=HistNode-> getHistory(false,xiQSh);
  
  if( Born!= HistNode){
    if (UseRandom::rnd()<.5){
      weight=-2.; projected=true;  Node->nodeME()->projectorStage(1);}
    else{
      weight= 2.; projected=false; Node->nodeME()->projectorStage(0);}
  }else{
    weight=1.;  projected=false; Node->nodeME()->projectorStage(0);
  }
  if (!Born->xcomb()->willPassCuts())return ZERO;
  
  Energy startscale=CKKW_StartScale(Born);
  fillHistory( startscale,  Born, HistNode);
  Node->runningPt(fillProjector());
  weight*=history.back().weight*alphaReweight()*pdfReweight();
    //cout<<"\n"<<history.back().weight<<" "<<alphaReweight()<<" "<<pdfReweight();
  if(weight==0.)return ZERO;
  
  CrossSection sumPS=ZERO;

  for(auto const & child : Node->children()){
    if (child->allAbove(mergePt())){
      if((child)->pT()>mergePt()/3.)//TODO: this is a dynamical cutoff(see below)
        sumPS+=child->calcPs(startscale*xiFacME);
      else
        sumPS+=child->calcDip(startscale*xiFacME);
    }else{
      assert(child->pT()>mergePt());
    }
  }
  
  CrossSection me=TreedSigDR(startscale,Node);
  
  if (Node->legsize()==6&&Debug::level > 2) {
    Energy minPT=1000*GeV;
    for(auto const & no :Node->children() )minPT=min(no->pT(),minPT);
  
  cout<<"\nRBSR "<<minPT/GeV<< " " <<weight<<" "<<(me-sumPS)/nanobarn<<" "<<me/nanobarn<<" "<<sumPS/nanobarn;
  }
    //Here we subtract the PS (and below the dynamical cutoff the Dip)
  return weight*as(startscale*xiRenSh)/SM().alphaS()*
  (me-sumPS);
}



CrossSection Merger::MergingDSigDRRealBelowSubInt(NPtr Node){
  if(Node->children().empty())return ZERO;
  NPtr CLNode= Node->randomChild();
  if(CLNode->pT()<mergePt()/3.)return ZERO;//TODO: this is a dynamical cutoff(see above)

  if (!CLNode->children().empty()) {
    auto inoutpair=getInOut(CLNode);
    NPtr rndCh= CLNode->randomChild(); // Here we make sure that clustering and splitting are invertible
    if(!matrixElementRegion(inoutpair.first,inoutpair.second, rndCh->pT(), theMergePt))return ZERO;
  }
  
  
  NPtr Born=CLNode-> getHistory(false,xiQSh);
  if( Born!= CLNode){
    if (UseRandom::rnd()<.5){
      weight=-2.; projected=true;  Node->nodeME()->projectorStage(2);}
    else{
      weight= 2.; projected=false; Node->nodeME()->projectorStage(1);}
  }else{
    weight=1.;  projected=false; Node->nodeME()->projectorStage(1);
  }
  if (!CLNode->allAbove(mergePt()))return ZERO;
  if (!Born->xcomb()->willPassCuts())return ZERO;
  Energy startscale=CKKW_StartScale(Born);
  fillHistory( startscale,  Born, CLNode);
  Node->runningPt(fillProjector());
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.)return ZERO;
  
  

  
  pair<CrossSection,CrossSection> DipAndPs=
  CLNode->calcDipandPS(startscale*xiFacME);
  
  if (Node->legsize()==6&&Debug::level > 2) {
    Energy minPT=1000*GeV;
    for(auto const & no :Node->children() )minPT=min(no->pT(),minPT);
  

    cout<<"\nRBSI "<<CLNode->pT()/GeV<<" "<<weight<<" "<<(DipAndPs.second-DipAndPs.first)/nanobarn<<" "<<DipAndPs.second/nanobarn<<" "<<DipAndPs.first/nanobarn;
  }
    //Here we add the PS and subtrac the Dip (only above the dynamical cutoff)
  return weight*as(startscale*xiRenSh)/SM().alphaS()*
  (double)Node->children().size()*(DipAndPs.second-DipAndPs.first);
}










CrossSection Merger::TreedSigDR(Energy startscale,NPtr Node,double diffAlpha){
  
  NPtr DeepHead=Node;//->deepHead();
  renormscale(startscale);
  DeepHead->nodeME()->factory()->scaleChoice()->setXComb(DeepHead->xcomb());
  DeepHead->nodeME()->setScale(sqr(startscale),sqr(startscale));
  theCalculateInNode=false;
  CrossSection res=DeepHead->nodeME()->dSigHatDRB();
  if (diffAlpha!=1.) {
    res+=DeepHead->nodeME()->dSigHatDRAlphaDiff(diffAlpha);
  }
  renormscale(0.0*GeV);
  theCalculateInNode=true;
  if(std::isnan(double(res/nanobarn))){cout<<"\nTreedSigDR is nan";res=ZERO;};
  return res;
}

CrossSection Merger::LoopdSigDR(Energy startscale,NPtr Node){
    // The deephead should be calculated here.
  NPtr DeepHead=Node;//->deepHead();
  renormscale(startscale);
  DeepHead->nodeME()->setXComb(DeepHead->xcomb());
  DeepHead->nodeME()->factory()->scaleChoice()->setXComb(DeepHead->xcomb());
  DeepHead->nodeME()->setScale(sqr(startscale),sqr(startscale));
  theCalculateInNode=false;
  DeepHead->nodeME()->doOneLoopNoBorn();
  CrossSection res=DeepHead->nodeME()->dSigHatDRV()+DeepHead->nodeME()->dSigHatDRI();
  DeepHead->nodeME()->noOneLoopNoBorn();
  renormscale(0.0*GeV);
  theCalculateInNode=true;
  return res;
}

Energy Merger::fillProjector(){
    // in the shower handler the scale is multiplied by xiQSh
    // so here we need to devide this
  double xiQShfactor=history.begin()->node->legsize()==N0()?xiQSh:1.;
  if(history.begin()->node->deepHead()->nodeME()->projectorStage() == 0){
    return (history.size()==1?1.:(1./xiQShfactor))*history.back().scale;
  }
  for(auto const & hs : history)
    if (projectorStage(hs.node)&&history.begin()->node->deepHead()->nodeME()->projectorStage() != 0){
      history.begin()->node->deepHead()->xcomb()->lastProjector(hs.node->xcomb());
      return (hs.node==history[0].node?1.:(1./xiQShfactor))*hs.scale;
    }
  
  assert(false);
  return ZERO;
}

double Merger::pdfReweight(){
  double res=1.;
  double facfactor=(history[0].node->legsize()==N0()?xiFacME:xiFacSh);
  for(int side : {0,1}){
    if(history[0].node->xcomb()->mePartonData()[side]->coloured()){
      for (auto const & hs : history){
          //pdf-ratio only to the last step
        if(!hs.node->parent())continue;
        if (hs.node==history.back().node)continue;
        if(!hs.node->dipole()){
          cout<<"\nthis should not happen";
          return 0.;
        }
        res *= pdfratio(hs.node, facfactor*hs.scale,xiFacSh*(hs.node->pT()), side);
        facfactor=xiFacSh;
      }
      res*=pdfratio(history.back().node,facfactor*history.back().scale ,history[0].scale*xiFacME, side);
    }
  }
  if (std::isnan(res))cout<<"\nWarning:pdfReweight is nan.";
  return res;
}

double Merger::alphaReweight(){
  double res=1.;
  Energy Q_R=(history[0].node->legsize()==N0()?xiRenME:xiRenSh)*history[0].scale;
  res *= pow(as(Q_R) / SM().alphaS(), history[0].node->nodeME()->orderInAlphaS());
  res *= pow(SM().alphaEMME(history[0].node->nodeME()->factory()->scaleChoice()->renormalizationScaleQED())/ SM().alphaEMMZ(), history[0].node->nodeME()->orderInAlphaEW());
  
  
  
  if (!(history[0].node->children().empty())){
    res *=pow((theCMWScheme?(1.+((3.*(67./18.-1./6.*Constants::pi*Constants::pi)-5./9.*Nf(Q_R))*as(Q_R))/2./Constants::pi):1.),int(history[0].node->legsize()-N0()));
  }

  
  for (auto const & hs : history)
    if (hs.node!=history.back().node){
      Energy q_i=xiRenSh* hs.node->pT();
      res *= as(q_i)/ SM().alphaS()
      *(theCMWScheme?(1.+((3.*(67./18.-1./6.*Constants::pi*Constants::pi)-5./9.*Nf(q_i))*as(q_i))/2./Constants::pi):1.);
      if (std::isnan(res))cout<<"\nWarning:alphaReweight is nan:q_i "<<q_i/GeV<<" "<<Nf(q_i);
    }
  
  if (std::isnan(res))cout<<"\nWarning:alphaReweight is nan.";
  return res;
}

void Merger::fillHistory(Energy scale, NPtr Begin, NPtr EndNode){
  
  history.clear();
  double sudakov0_n=1.;
  history.push_back(HistoryStep(Begin,sudakov0_n,scale));
  
  
  double xiQShfactor=history.begin()->node->legsize()==N0()?xiQSh:1.;
  
  scale*=xiQShfactor;
  if (Begin->parent()||!isUnitarized) {
    while (Begin->parent() && (Begin != EndNode)) {
      if (!dosudakov(Begin,scale, Begin->pT(),sudakov0_n)){
        history.push_back(HistoryStep(Begin->parent(),0.,scale));
      }
      scale=Begin->pT();
      
      if (std::isnan(sudakov0_n))cout<<"\nWarning:sudakov0_n is nan.";
      history.push_back(HistoryStep(Begin->parent(),sudakov0_n,Begin->pT()));
      Begin = Begin->parent();
    }
    
    Energy notunirunning=scale;
    
    if (!isUnitarized&&N()+N0() > int(Begin->deepHead()->legsize())) {
      if (!dosudakov(Begin,notunirunning,mergePt(),sudakov0_n)){
        history.back().weight=0.;
      }else{
        history.back().weight=sudakov0_n;
      }
    }
  }
  if(  history.size()==1)scale/=xiQSh;
}




double Merger::sumpdfReweightUnlops(){
  double res=0.;
  Energy beam1Scale=history[0].scale*(history[0].node->legsize()==N0()?xiFacME:xiFacSh);
  Energy beam2Scale=history[0].scale*(history[0].node->legsize()==N0()?xiFacME:xiFacSh);
  
  for (auto const & hs : history){
      //pdf expansion only to the last step
    if(!hs.node->parent())continue;
    if(hs.node->xcomb()->mePartonData()[0]->coloured()&&hs.node->nodeME()->lastX1()>0.99)return 0.;
    if(hs.node->xcomb()->mePartonData()[1]->coloured()&&hs.node->nodeME()->lastX2()>0.99)return 0.;
    
    if(hs.node->nodeME()->lastX1()<0.00001)return 0.;
    if(hs.node->nodeME()->lastX2()<0.00001)return 0.;
    
    if (hs.node->dipole()->bornEmitter() == 0 ){
      res +=pdfUnlops(hs.node,0,
                      beam1Scale,
                      (hs.node->pT()),
                      hs.node->nodeME()->lastX1(),
                      Nf(history[0].scale),
                      history[0].scale);
      beam1Scale=(hs.node->pT())*xiFacSh;
    }
    if (hs.node->dipole()->bornEmitter() == 1 ){
      res +=pdfUnlops(hs.node,1,
                      beam2Scale,
                      (hs.node->pT()),
                      hs.node->nodeME()->lastX2(),
                      Nf(history[0].scale),
                      history[0].scale);
      beam2Scale=(hs.node->pT())*xiFacSh;
    }
    if (hs.node->dipole()->bornSpectator() == 0 &&hs.node->dipole()->bornEmitter() >1){//
      res +=pdfUnlops(hs.node,0,
                      beam1Scale,
                      (hs.node->pT()),
                      hs.node->nodeME()->lastX1(),
                      Nf(history[0].scale),
                      history[0].scale);
        //pdfratio(hs.node, beam1Scale, sqrt(hs.node->pT()), 1);
      beam1Scale=(hs.node->pT())*xiFacSh;
    }
    if (hs.node->dipole()->bornSpectator() == 1 &&hs.node->dipole()->bornEmitter() >1){//
      res +=pdfUnlops(hs.node,1,
                      beam2Scale,
                      (hs.node->pT()),
                      hs.node->nodeME()->lastX2(),
                      Nf(history[0].scale),
                      history[0].scale);
        //pdfratio(hs.node, beam2Scale , sqrt(hs.node->pT()), 2);
      beam2Scale=(hs.node->pT())*xiFacSh;
    }
  }
  
  if (history[0].node->deepHead()->xcomb()->mePartonData()[0]->coloured()){
    res +=pdfUnlops(history.back().node,0,
                    beam1Scale,
                    history[0].scale*xiFacME,
                    (history.back()).node->nodeME()->lastX1(),
                    Nf(history[0].scale),
                    history[0].scale);
    
  }
  if (history[0].node->deepHead()->xcomb()->mePartonData()[1]->coloured()) {
    res +=pdfUnlops(history.back().node,1,
                    beam2Scale,
                    history[0].scale*xiFacME,
                    (history.back()).node->nodeME()->lastX2(),
                    Nf(history[0].scale),
                    history[0].scale);
  }
  return res;
}






double Merger::pdfUnlops(NPtr node ,int side,Energy  running, Energy next,double x,int nlp,Energy fixedScale)  {
  
  tcPDPtr particle,parton;
  tcPDFPtr pdf;
  if (side==0) {
    particle = node->nodeME()->lastParticles().first->dataPtr();
    parton=node->nodeME()->lastPartons().first->dataPtr();
    pdf =node->xcomb()->partonBins().first->pdf();
  }else{
    assert(side==1);
    particle = node->nodeME()->lastParticles().second->dataPtr();
    parton=node->nodeME()->lastPartons().second->dataPtr();
    pdf =node->xcomb()->partonBins().second->pdf();
  }
  
    //copied from PKOperator
  double res=0.;
  int number=10;//*int((max(running/next,next/running)));
  for (int nr=0;nr<number;nr++){
    double restmp=0;
    
    using namespace RandomHelpers;
    double r = UseRandom::rnd();
    double eps = 1e-3;
    
    pair<double,double> zw =
    generate((piecewise(),
              flat(0.0,x),
              match(inverse(0.0,x,1.0) + inverse(1.0+eps,x,1.0))),r);
    
    double z = zw.first;
    double mapz = zw.second;
    double PDFxparton=pdf->xfx(particle,parton,sqr(fixedScale),x)/x;
    double CA = SM().Nc();
    double CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    
    
    if (abs(parton->id()) < 7) {
      
      double PDFxByzgluon=pdf->xfx(particle,getParticleData(ParticleID::g),sqr(fixedScale),x/z)*z/x;
      double PDFxByzparton=pdf->xfx(particle,parton,sqr(fixedScale),x/z)*z/x;
      assert(abs(parton->id()) < 7);
      
      restmp += CF*(3./2.+2.*log(1.-x)) * PDFxparton;
      if ( z > x ) {
        restmp+= 0.5 * ( sqr(z) + sqr(1.-z) ) * PDFxByzgluon / z;
        restmp += CF*2.*(PDFxByzparton - z*PDFxparton)/(z*(1.-z));
        restmp -= CF*PDFxByzparton * (1.+z)/z;
      }
    }else{
      
      assert(parton->id() == ParticleID::g);
      double PDFxByzgluon=pdf->xfx(particle,getParticleData(ParticleID::g),sqr(fixedScale),x/z)*z/x;
      if ( z > x ){
        double factor = CF * ( 1. + sqr(1.-z) ) / sqr(z);
        
        
        for ( int f = -nlp; f <= nlp; ++f ) {
          if ( f == 0 )
            continue;
          restmp += pdf->xfx(particle,getParticleData(f),sqr(fixedScale),x/z)*z/x*factor;
        }
      }
      
      restmp += ( (11./6.) * CA - (1./3.)*Nf(history[0].scale) + 2.*CA*log(1.-x) ) *PDFxparton;
      if ( z > x ) {
        restmp += 2. * CA * ( PDFxByzgluon - z*PDFxparton ) / (z*(1.-z));
        restmp += 2.* CA *( (1.-z)/z - 1. + z*(1.-z) ) * PDFxByzgluon / z;
      }
      
    }
    if (PDFxparton<1e-8)restmp= 0.;
    res+=1*restmp*log(sqr(running/next))/PDFxparton*mapz;
    
  }
  return res/number;
}



double Merger::sumalphaReweightUnlops(){
  double res=0.;
  if (!(history[0].node->children().empty())){
    res +=alphasUnlops(history[0].scale*xiRenSh,
                       history[0].scale*xiRenME)*int(history[0].node->legsize()-N0());
    assert(int(history[0].node->legsize()-N0()>0));
  }
    // dsig is calculated with fixed alpha_s
  for (auto const & hs : history){
      //expansion only to the last step
    if(!hs.node->parent())continue;
    res +=alphasUnlops(hs.node->pT()*xiRenSh ,history[0].scale);
  }
  return res;
}

double Merger::sumfillHistoryUnlops(){
  double res=0.;
  double xiQShfactor=history[0].node->legsize()==N0()?xiQSh:1.;
  for (auto const & hs : history){
    if(!hs.node->parent())continue;
    doUNLOPS(hs.node,(hs.node == history[0].node?xiQShfactor:1.)*hs.scale , hs.node->pT() , history[0].scale, res);
  }
  return res;
}


Ptr<MFactory>::ptr Merger::treefactory(){return theTreeFactory;}


void Merger::doinit(){
  assert(DSH()->hardScaleIsMuF());
  
}


CrossSection Merger::MergingDSigDR() {
  
  history.clear();
  
  
  if (theFirstNodeMap.find(theCurrentME)==theFirstNodeMap.end()) {
    cout<<"\nnot in map:"<<theFirstNodeMap.size();
    assert(false);
  }
  
  
  NPtr Node = theFirstNodeMap[theCurrentME];  assert(Node);
  
  
  xiRenME=theCurrentME->renormalizationScaleFactor();
  xiFacME=theCurrentME->factorizationScaleFactor();
  xiRenSh=DSH()->renormalizationScaleFactor();
  xiFacSh=DSH()->factorizationScaleFactor();
  xiQSh=DSH()->hardScaleFactor();
  
  DSH()->setCurrentHandler();
  
  DSH()->eventHandler(generator()->eventHandler());
  assert(DSH()->hardScaleIsMuF());
  
  
  CrossSection res=ZERO;
  if(Node->deepHead()->subtractedReal()){
    res=MergingDSigDRRealStandard(Node);
    theCurrentMaxLegs=maxLegsNLO();
  }else if(Node->deepHead()->virtualContribution()){
    res=MergingDSigDRVirtualStandard(Node);
    theCurrentMaxLegs=maxLegsNLO();
  }else if(theGamma!=1.){
    res=MergingDSigDRBornGamma(Node);
    theCurrentMaxLegs=maxLegsLO();
  }else{
    res=MergingDSigDRBornStandard(Node);
    theCurrentMaxLegs=maxLegsLO();
  }
  
  
  theCurrentME->lastXCombPtr()->lastCentralScale(sqr(Node->runningPt()));
  theCurrentME->lastXCombPtr()->lastShowerScale(sqr(Node->runningPt()));
  if(theCurrentME->lastXCombPtr()->lastProjector()){
    theCurrentME->lastXCombPtr()->lastProjector()->lastCentralScale(sqr(Node->runningPt()));
    theCurrentME->lastXCombPtr()->lastProjector()->lastShowerScale(sqr(Node->runningPt()));
  }
  
  renormscale(0.0*GeV);
  if (res == ZERO){
    history.clear();
    return res;
  }
  
  
  cleanup(Node);
  DSH()->eventRecord().clear();
  theCurrentME->lastXCombPtr()->subProcess(SubProPtr());
  
  history.clear();
  
  if(std::isnan(double(res/nanobarn))|| !std::isfinite(double(res/nanobarn))){
    cout<<"Warning merger weight is "<<res/nanobarn<<" -> setting to 0";
    return ZERO;
  }
  
  return res;
  
}



void Merger::CKKW_PrepareSudakov(NPtr Born,Energy running){
    //cleanup(Born);
  tSubProPtr sub = Born->xcomb()->construct();
  DSH()->resetPDFs(make_pair(Born->xcomb()->partonBins().first->pdf(),
                             Born->xcomb()->partonBins().second->pdf()),
                   Born->xcomb()->partonBins());
  DSH()->setCurrentHandler();
  
  DSH()->currentHandler()->generator()->currentEventHandler(Born->deepHead()->xcomb()->eventHandlerPtr());
  
  DSH()->currentHandler()->remnantDecayer()->setHadronContent(Born->deepHead()->xcomb()->lastParticles());
  DSH()->eventRecord().clear();
  DSH()->eventRecord().slimprepare(sub, dynamic_ptr_cast<tStdXCombPtr>(Born->xcomb()), DSH()->pdfs(), Born->deepHead()->xcomb()->lastParticles());
  DSH()->hardScales(sqr(running));
}


Energy Merger::CKKW_StartScale(NPtr Born){
  Energy res=generator()->maximumCMEnergy();;
  if(!Born->children().empty()){
    for (int i=0;i<Born->legsize();i++){
      if (!Born->nodeME()->mePartonData()[i]->coloured())continue;
      for (int j=i+1;j<Born->legsize();j++){
        if (i==j||!Born->nodeME()->mePartonData()[j]->coloured())continue;
        res= min(res,sqrt(2.*Born->nodeME()->lastMEMomenta()[i]*Born->nodeME()->lastMEMomenta()[j]));
      }
    }
  }else{
    Born->nodeME()->factory()->scaleChoice()->setXComb(Born->xcomb());
    res= sqrt(Born->nodeME()->factory()->scaleChoice()->renormalizationScale());
  }
  Born->nodeME()->factory()->scaleChoice()->setXComb(Born->xcomb());
  res=max(res, sqrt(Born->nodeME()->factory()->scaleChoice()->renormalizationScale()));
  return res;
}




double Merger::alphasUnlops( Energy next,Energy fixedScale)  {
  double betaZero =  (11./6.)*SM().Nc() - (1./3.)*Nf(history[0].scale);
  return (betaZero*log(sqr(fixedScale/next)))+
  (theCMWScheme?(((3.*(67./18.-1./6.*Constants::pi*Constants::pi)-5./9.*Nf(history[0].scale)))):0.);
}


double Merger::pdfratio(NPtr  Born,Energy  nominator_scale, Energy denominator_scale,int side){
  
  StdXCombPtr bXc = Born->xcomb();
  assert (bXc->mePartonData()[side]->coloured());
  double from,to;
  if (side==0){
    if (denominator_scale==nominator_scale) {
      return 1.;
    }
    from=Born->nodeME()->pdf1(sqr(denominator_scale));
    to  =Born->nodeME()->pdf1(sqr( nominator_scale ));
    if ( (to < 1e-8||from < 1e-8)&&(to/from>10000000.)){cout<<"\npdfratio to="<<to<<" from="<<from;return 0.;}
  }
  else{
    if (denominator_scale==nominator_scale) {
      return 1.;
    }
    from=Born->nodeME()->pdf2(sqr(denominator_scale));
    to=Born->nodeME()->pdf2( sqr( nominator_scale ));
    if ( (to < 1e-8||from < 1e-8)&&(to/from>10000000.)){cout<<"\npdfratio to="<<to<<" from="<<from;return 0.;}
  }
  return to/from;
}



bool Merger::dosudakov(NPtr Born,Energy running, Energy next, double& sudakov0_n) {
  CKKW_PrepareSudakov(Born, running);
  for(DipoleChain const & chain : DSH()->eventRecord().chains()){
    for(Dipole const & dip : chain.dipoles()){
      sudakov0_n*=singlesudakov( dip, next,running,make_pair(true,false) );
      sudakov0_n*=singlesudakov( dip, next,running,make_pair(false,true) );
      if (sudakov0_n==0.0){
        cleanup(Born);
        return false;
      }
    }
  }
  cleanup(Born);
  return true;
}

bool Merger::doUNLOPS(NPtr Born,Energy  running, Energy next,Energy fixedScale, double& UNLOPS) {
  CKKW_PrepareSudakov(Born, running);
  for(DipoleChain const & chain : DSH()->eventRecord().chains()){
    for(Dipole const & dip : chain.dipoles()){
      UNLOPS+=singleUNLOPS( dip, next,running,fixedScale,make_pair(true,false) );;
      UNLOPS+=singleUNLOPS( dip, next,running,fixedScale,make_pair(false,true) );
    }
  }
  cleanup(Born);
  return true;
}



bool Merger::projectorStage(NPtr  Born){
  
  return	(Born->deepHead()->nodeME()->projectorStage() ==
             int((Born->deepHead()->legsize()
                  - Born->legsize())));
}

void Merger::cleanup(NPtr Born) {
  DSH()->eventRecord().clear();
  if(!Born->xcomb()->subProcess())return;
  ParticleVector vecfirst = Born->xcomb()->subProcess()->incoming().first->children();
  for(auto const & particle : vecfirst)
    Born->xcomb()->subProcess()->incoming().first->abandonChild(particle);
  
  ParticleVector vecsecond = Born->xcomb()->subProcess()->incoming().second->children();
  for(auto const & particle : vecsecond)
    Born->xcomb()->subProcess()->incoming().second->abandonChild(particle);
  Born->xcomb()->subProcess(SubProPtr());
}

double Merger::singlesudakov(Dipole dip ,Energy next,Energy running,pair<bool,bool> conf ){
  
  double res=1.;
  tPPtr emitter = dip.emitter(conf);
  tPPtr spectator = dip.spectator(conf);
  DipoleSplittingInfo candidate(dip.index(conf),conf,
                                dip.emitterX(conf),
                                dip.spectatorX(conf),
                                emitter,spectator);
  
  
  if ( DSH()->generators().find(candidate.index()) == DSH()->generators().end() ) DSH()->getGenerators(candidate.index());
  
  auto const & gens = DSH()->generators().equal_range(candidate.index());
  
  for ( auto   gen = gens.first; gen != gens.second; ++gen ) {
    if ( !(gen->first == candidate.index()) )
      continue;
    
    Energy dScale =	gen->second->splittingKinematics()->dipoleScale(emitter->momentum(),spectator->momentum());
    candidate.scale(dScale);
    candidate.continuesEvolving();
    Energy ptMax=gen->second->splittingKinematics()->ptMax(candidate.scale(),candidate.emitterX(), candidate.spectatorX(),
                                                             candidate.index(),*gen->second->splittingKernel());
    
    candidate.hardPt(min(running,ptMax));
    
    if (candidate.hardPt()>next){
      res*=gen->second->sudakov(candidate,next);
    }
  }
  
  return res;
}


double Merger::singleUNLOPS(Dipole dip ,Energy next,Energy running,Energy fixedScale,pair<bool,bool> conf ){
  
  double res=0.;
  tPPtr emitter = dip.emitter(conf);
  tPPtr spectator = dip.spectator(conf);
  DipoleSplittingInfo candidate(dip.index(conf),conf,dip.emitterX(conf),dip.spectatorX(conf),emitter,spectator);
  
  if ( DSH()->generators().find(candidate.index()) == DSH()->generators().end() ) DSH()->getGenerators(candidate.index());
  
   auto const & gens = DSH()->generators().equal_range(candidate.index());
  
  for ( auto gen = gens.first; gen != gens.second; ++gen ) {
    if ( !(gen->first == candidate.index()) )
      continue;
    Energy dScale =	gen->second->splittingKinematics()->dipoleScale(emitter->momentum(),spectator->momentum());
    candidate.scale(dScale);
    candidate.continuesEvolving();
    Energy ptMax=gen->second->splittingKinematics()->ptMax(candidate.scale(),candidate.emitterX(), candidate.spectatorX(),
                                                           candidate.index(),*gen->second->splittingKernel());
    
    candidate.hardPt(min(running,ptMax));
    if (candidate.hardPt()>next){
      res+=gen->second->unlops(candidate,next,fixedScale);
    }
  }
	 
	 return res;
}




void Merger::firstNodeMap(Ptr<MatchboxMEBase>::ptr a,NPtr b){theFirstNodeMap.insert(make_pair(a,b));}



map<Ptr<MatchboxMEBase>::ptr,NPtr> Merger::firstNodeMap(){return theFirstNodeMap;}




void Merger::setXComb(Ptr<MatchboxMEBase>::ptr me,tStdXCombPtr xc,int st){
  theFirstNodeMap[me]->setXComb(xc, st);
}
void Merger::setKinematics(Ptr<MatchboxMEBase>::ptr me){
  theFirstNodeMap[me]->setKinematics();
}
void Merger::clearKinematics(Ptr<MatchboxMEBase>::ptr me){
  theFirstNodeMap[me]->clearKinematics();
}
bool Merger::generateKinematics(Ptr<MatchboxMEBase>::ptr me,const double * r){
  
  theFirstNodeMap[me]->firstgenerateKinematics(r, 0);
  
  if (theFirstNodeMap[me]->cutStage()==0 ){
    
    bool inAlphaPS=false;
    NPtrVec children=theFirstNodeMap[me]->children();
    for (Ptr<Node>::ptr const & child: children) {
      treefactory()->setAlphaParameter(theGamma);
      inAlphaPS|=theGamma!=1.?child->dipole()->aboveAlpha():false;
      treefactory()->setAlphaParameter(1.);
    }
    
    SafeClusterMap temp=theFirstNodeMap[me]->clusterSafe();
    for (auto const & cs: temp) {
      if (!cs.second.first&&!inAlphaPS)return false;
    }
  }
  if (theFirstNodeMap[me]->cutStage()==1 ){
    SafeClusterMap temp=theFirstNodeMap[me]->clusterSafe();
    for (auto const & sc: temp) {
      if (!sc.second.first && !sc.second.second)return false;
    }
  }
  return true;
  
}
bool Merger::calculateInNode() const{
  return theCalculateInNode;
}


void Merger::fillProjectors(Ptr<MatchboxMEBase>::ptr me){
  for (auto const & propair: theFirstNodeMap[me]->Projector()) {
    me->lastXCombPtr()->projectors().insert(propair.first,
                                            propair.second->xcomb());
  }
}
pair<bool,bool> Merger::clusterSafe(Ptr<MatchboxMEBase>::ptr me,int emit,int emis,int spec){
  return theFirstNodeMap[me]->clusterSafe().find(make_pair(make_pair(emit,emis),spec))->second;
  
}


bool Merger::matrixElementRegion(PVector incoming,PVector outgoing,Energy winnerScale,Energy cutscale){
  
    //cout<<"\nparticles s"<<particles.size()<<" "<<particles[0]<<" "<<particles[1]<<flush;
  /*
   if (defMERegionByJetAlg && !particles[0]->coloured()&& !particles[1]->coloured()) {
   assert(false);
   vector<fastjet::PseudoJet> input_particles;
   for(size_t em=2; em < particles.size();em++){
   input_particles.push_back(fastjet::PseudoJet(em->momentum().x()/GeV,
   em->momentum().y()/GeV,
   em->momentum().z()/GeV,
   em->momentum().e()/GeV));
   }
   fastjet::JetDefinition jet_def(fastjet::ee_kt_algorithm);
   fastjet::ClusterSequence clust_seq(input_particles, jet_def);
   size_t n = particles.size()-2;
   vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets_ycut(ee_ycut);
   return n==exclusive_jets.size();
   }
   
   
   if (defMERegionByJetAlg ) {
   assert(false);
   size_t noncol=0;
   vector<fastjet::PseudoJet> input_particles;
   for(size_t em=2; em < particles.size();em++){
   if (em->coloured())
   input_particles.push_back(fastjet::PseudoJet(em->momentum().x()/GeV,
   em->momentum().y()/GeV,
   em->momentum().z()/GeV,
   em->momentum().e()/GeV));
   else
   noncol++;
   }
   double Rparam = 1.0;
   fastjet::Strategy strategy = fastjet::Best;
   fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
   fastjet::JetDefinition jet_def(fastjet::kt_algorithm, Rparam, recomb_scheme, strategy);
   
   fastjet::ClusterSequence clust_seq(input_particles, jet_def);
   size_t n = particles.size()-2-noncol;
   vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets(pp_dcut);
   
   // cout<<"\nn= "<<n<<" jets= "<<exclusive_jets.size();
   //for (size_t i=0; i<exclusive_jets.size(); i++) {
   //cout<<"\nn= "<<n<<" jetspt= "<<exclusive_jets[i].perp();
   //}
   
   return n==exclusive_jets.size();
   }
   
   */
  
    //assert(false);
  
  
  
  Energy ptx=1000000.*GeV;
  bool foundwinnerpt=false;
  using namespace boost;
    //FF
  
  for(auto const & em : outgoing){        if (! em->coloured()) continue;
    for(auto const & emm : outgoing){     if (!emm->coloured()) continue; if (em==emm) continue;
      for(auto const & spe : outgoing){   if (!spe->coloured()) continue; if (em==spe||emm==spe) continue;
        
        if (!(em->id()==-emm->id()||emm->id()>6))continue;

        Lorentz5Momentum emittermom = em->momentum();
        Lorentz5Momentum emissionmom = emm->momentum();
        Lorentz5Momentum spectatormom = spe->momentum();
        Energy pt=0*GeV;
        if (emittermom.m()<=0.001*GeV&&emissionmom.m()<=0.001*GeV&&spectatormom.m()<=0.001*GeV) {
          pt=FFLTK->lastPt(emittermom,emissionmom,spectatormom);
        }else{
          pt=FFMTK->lastPt(emittermom,emissionmom,spectatormom);
        }
        
        if (abs(pt-winnerScale)<0.001*GeV) {
          foundwinnerpt=true;
        }
        ptx =min(ptx,pt);
      }
    }
  }
  
    //FI
  for(auto const & spe : incoming){          if (! spe->coloured()) continue;
    for(auto const & emm : outgoing){        if (! emm->coloured()) continue;
      for(auto const & em : outgoing){       if (! em->coloured()) continue; if (em==emm) continue;
        
        if (!(em->id()==-emm->id()||emm->id()>6))continue;
        
        Lorentz5Momentum emittermom = em->momentum();
        Lorentz5Momentum emissionmom = emm->momentum();
        Lorentz5Momentum spectatormom = spe->momentum();
        Energy pt=0*GeV;
        if (emittermom.m()<=0.001*GeV&&emissionmom.m()<=0.001*GeV&&spectatormom.m()<=0.001*GeV) {
          pt=FILTK->lastPt(emittermom,emissionmom,spectatormom);
        }else{
          pt=FIMTK->lastPt(emittermom,emissionmom,spectatormom);
        }
        
        
        if (abs(pt-winnerScale)<0.001*GeV) {
          foundwinnerpt=true;
        }
        
        if(pt>0.*GeV)
          ptx =min(ptx,pt);
      }
    }
  }
  
    //IF
  for(auto const & em : incoming){          if (! em->coloured()) continue;
    for(auto const & emm : outgoing){        if (! emm->coloured()) continue;
      for(auto const & spe : outgoing){       if (! spe->coloured()) continue; if (emm==spe) continue;
        
        if (!(em->id()>6|| em->id()==emm->id() ||emm->id()>6))continue;
        
        Lorentz5Momentum emittermom = em->momentum();
        Lorentz5Momentum emissionmom = emm->momentum();
        Lorentz5Momentum spectatormom = spe->momentum();
        Energy pt=0*GeV;
        
        if (emittermom.m()<=0.001*GeV&&emissionmom.m()<=0.001*GeV&&spectatormom.m()<=0.001*GeV) {
          pt=IFLTK->lastPt(emittermom,emissionmom,spectatormom);
        }else{
          pt=IFMTK->lastPt(emittermom,emissionmom,spectatormom);
        }
        
        
        if (abs(pt-winnerScale)<0.01*GeV) {
          foundwinnerpt=true;
        }
        ptx =min(ptx,pt);
      }
    }
  }
  
    //II
  for(auto const & em : incoming){          if (! em->coloured()) continue;
    for(auto const & spe : incoming){       if (! spe->coloured()) continue; if (em==spe) continue;
      for(auto const & emm : outgoing){        if (! emm->coloured()) continue;
      
        
        if (!(em->id()>6||em->id()==emm->id() ||emm->id()>6))continue;
        
        Lorentz5Momentum emittermom = em->momentum();
        Lorentz5Momentum emissionmom = emm->momentum();
        Lorentz5Momentum spectatormom = spe->momentum();
        
        Energy  pt=IILTK->lastPt(emittermom,emissionmom,spectatormom);
        
        if (abs(pt-winnerScale)<0.01*GeV) {
          foundwinnerpt=true;
        }
        ptx =min(ptx, pt);
      }
    }
  }
  
  if(!foundwinnerpt){
    cout<<"\nWarning: could not find winner with pt: "<<winnerScale/GeV;
    for(auto const & emm : incoming){
      cout<<"\n"<<emm->id()<<" "<<emm->momentum()/GeV<<" "<<emm->momentum().m()/GeV<<flush;
    }
    for(auto const & emm : outgoing){
      cout<<"\n"<<emm->id()<<" "<<emm->momentum()/GeV<<" "<<emm->momentum().m()/GeV<<flush;
    }
  
  }
  
  return (ptx>cutscale) ;
  
}





int Merger::M()const{return theTreeFactory->M();}

int Merger::N()const{return theTreeFactory->N();}











  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Merger::persistentOutput(PersistentOStream & os) const {
  
  os <<  Unlopsweights<<
  theCMWScheme<<   projected<<
  isUnitarized<<   isNLOUnitarized<<
  defMERegionByJetAlg<<theOpenInitialStateZ<<
  theChooseHistory<<
  theN0<<    theOnlyN<<
  xiRenME<<     xiFacME<<
  xiRenSh<<     xiFacSh<<
  xiQSh<<     weight<<weightCB<<theGamma<<ee_ycut<<pp_dcut<<theSmearing<<ounit(therenormscale, GeV)<<ounit(theIRSafePT, GeV)<<ounit(theMergePt, GeV)<<ounit(theCentralMergePt, GeV)<<theMergingJetFinder
  
  <<theLargeNBasis
  
  
  << FFLTK
  << FILTK
  << IFLTK
  << IILTK
  << FFMTK
  << FIMTK
  << IFMTK
  
  <<theDipoleShowerHandler<< theTreeFactory<<theFirstNodeMap ;
  
}

void Merger::persistentInput(PersistentIStream & is, int) {
  is >>  Unlopsweights>>
  theCMWScheme>>   projected>>
  isUnitarized>>   isNLOUnitarized>>
  defMERegionByJetAlg>>theOpenInitialStateZ>>
  theChooseHistory>>
  theN0>>    theOnlyN>>
  xiRenME>>     xiFacME>>
  xiRenSh>>     xiFacSh>>
  xiQSh>>
  weight>>weightCB>>theGamma>>ee_ycut>>pp_dcut>>theSmearing>>iunit(therenormscale, GeV)>>iunit(theIRSafePT, GeV)>>iunit(theMergePt, GeV)>>iunit(theCentralMergePt, GeV)>>theMergingJetFinder>>theLargeNBasis>>
  
  
  FFLTK
  >> FILTK
  >> IFLTK
  >> IILTK
  >> FFMTK
  >> FIMTK
  >> IFMTK>>
  
  theDipoleShowerHandler>> theTreeFactory>>theFirstNodeMap ;
}

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<Merger, Herwig::MergerBase>
describeHerwigMerger("Herwig::Merger", "HwDipoleShower.so");

  //ClassDescription<Merger> Merger::initMerger;
  // Definition of the static class description member.

void Merger::Init() {
  
  static ClassDocumentation<Merger> documentation
  ("Merger.");
  
  
  
  static Reference<Merger,DipoleShowerHandler> interfaceShowerHandler
  ("DipoleShowerHandler",
   "",
   &Merger::theDipoleShowerHandler, false, false, true, true, false);
  
  
  
  static Switch<Merger,bool>
  interfaceUnlopsweights ("Unlopsweights","",&Merger::Unlopsweights, false, false, false);
  static SwitchOption interfaceUnlopsweightsYes
  (interfaceUnlopsweights,"Yes","",true);
  static SwitchOption interfaceUnlopsweightsNo
  (interfaceUnlopsweights,"No","",false);
  
  static Switch<Merger,bool>
  interfacetheCMWScheme ("CMWScheme","",&Merger::theCMWScheme, false, false, false);
  static SwitchOption interfacetheCMWSchemeYes
  (interfacetheCMWScheme,"Yes","",true);
  static SwitchOption interfacetheCMWSchemeNo
  (interfacetheCMWScheme,"No","",false);
  
  
  
  
  
  
  
  
  static Parameter<Merger,Energy> interfaceMergerScale
  ("MergingScale",
   "The MergingScale.",
   &Merger::theCentralMergePt, GeV, 20.0*GeV, 0.0*GeV, 0*GeV,
   false, false, Interface::lowerlim);
  
  
  
  
  
  
  
  static Reference<Merger,MFactory> interfaceMergerHelper
  ("MFactory",
   "",
   &Merger::theTreeFactory, false, false, true, true, false);
  
  
  
  
  
  static Parameter<Merger,double> interfacedcut
  ("pp_dcut",
   "The MergingScale.",
   &Merger::pp_dcut, 5.0, 0.0, 0,
   false, false, Interface::lowerlim);
  
  static Parameter<Merger,double> interfaceycut
  ("ee_ycut",
   "The MergingScale.",
   &Merger::ee_ycut, 0.2, 0.0, 0,
   false, false, Interface::lowerlim);
  
  static Parameter<Merger,double> interfacegamma
  ("gamma",
   "gamma parameter.",
   &Merger::theGamma, 1.0, 0.0, 0,
   false, false, Interface::lowerlim);
  
  
  static Reference<Merger,JetFinder> interfaceMergingJetFinder
  ("MergingJetFinder",
   "A reference to the jet finder used in Merging.",
   &Merger::theMergingJetFinder, false, false, true, false, false);
  
  
  
  static Reference<Merger,ColourBasis> interfaceLargeNBasis
  ("LargeNBasis",
   "Set the large-N colour basis implementation.",
   &Merger::theLargeNBasis, false, false, true, true, false);
  
  
  
  
  
  
  static Switch<Merger,bool>
  interfacedefMERegionByJetAlg
  ("MERegionByJetAlg","",&Merger::defMERegionByJetAlg, false, false, false);
  
  static SwitchOption interfacedefMERegionByJetAlgYes
  (interfacedefMERegionByJetAlg,"Yes","",true);
  static SwitchOption interfacedefMERegionByJetAlgNo
  (interfacedefMERegionByJetAlg,"No","",false);
  
  
  static Switch<Merger,bool>
  interfaceOpenInitialSateZ
  ("OpenInitialStateZ","",&Merger::theOpenInitialStateZ, false, false, false);
  
  static SwitchOption interfaceOpenInitialSateZYes
  (interfaceOpenInitialSateZ,"Yes","",true);
  static SwitchOption interfaceOpenInitialSateZNo
  (interfaceOpenInitialSateZ,"No","",false);
  
  
  
  
  
  
  
  
  static Parameter<Merger, Energy>
  interfaceIRSafePT
  ("IRSafePT", "Set the pt for which a matrixelement should be treated as IR-safe.",
   
   &Merger::theIRSafePT,
   GeV, 0.0 * GeV, ZERO, Constants::MaxEnergy, true, false, Interface::limited);
  interfaceIRSafePT.setHasDefault(false);
  
  
  
  static Parameter<Merger, double> interfacemergePtsmearing("MergingScaleSmearing", "Set the percentage the merging pt should be smeared.",
                                                            &Merger::theSmearing, 0., 0.,
                                                            0.0, 0.5, true, false, Interface::limited);
  
  
  
  static Parameter<Merger, int> interfacechooseHistory("chooseHistory", "different ways of choosing the history", &Merger::theChooseHistory, 3, -1, 0,
                                                       false, false, Interface::lowerlim);
  
  
  
  
  
  
}
















