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
  
  minusL=false;
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










CrossSection Merger::MergingDSigDRBornGamma(NPtr Node){
  
  double weightCL=0.;
  weight=1.;


 
  if (!Node->children().empty()) {
    PVector particles;
    for (size_t i=0;i<Node->nodeME()->mePartonData().size();i++){
      Ptr<ThePEG::Particle>::ptr p =Node->nodeME()->mePartonData()[i]->produceParticle(Node->nodeME()->lastMEMomenta()[i]);
      particles.push_back(p);
    }
    if(!matrixElementRegion(particles, Node->children()[0]->dipol()->lastPt(), theMergePt))weight*= 0.;
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
  Energy projectedscale=startscale;
  fillHistory( startscale,  Born, Node);
  if (!fillProjector(projectedscale))weight=0.;
  Node->runningPt(projectedscale);
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.&&weightCL==0.)return ZERO;
 
  res+=weight*TreedSigDR(startscale,Node,(!maxMulti&&!projected)?theGamma:1.);
  }
  
  if(CLNode&&theGamma!=1.){
    Energy startscale=CKKW_StartScale(BornCL);
    Energy projectedscale=startscale;
    fillHistory( startscale,  BornCL, CLNode);
    if (!fillProjector(projectedscale))return ZERO;
    Node->runningPt(projectedscale);
    weightCL*=history.back().weight*alphaReweight()*pdfReweight();
    CrossSection diff=ZERO;
    Node->nodeME()->factory()->setAlphaParameter(1.);
    diff-=weightCL*CLNode->dipol()->dSigHatDR(sqr(startscale*xiFacME));
    Node->nodeME()->factory()->setAlphaParameter(theGamma);
    
    string alp=(CLNode->dipol()->aboveAlpha()?"Above":"Below");
    
    diff+=weightCL*CLNode->dipol()->dSigHatDR(sqr(startscale*xiFacME));
    Node->nodeME()->factory()->setAlphaParameter(1.);
    
    res+=diff;
  }
  
  
  return res;
}





CrossSection Merger::MergingDSigDRBornStandard(NPtr Node){
  

  
  if (!Node->children().empty()) {
    PVector particles;
    for (size_t i=0;i<Node->nodeME()->mePartonData().size();i++){
      Ptr<ThePEG::Particle>::ptr p =Node->nodeME()->mePartonData()[i]->produceParticle(Node->nodeME()->lastMEMomenta()[i]);
      particles.push_back(p);
    }
    if(!matrixElementRegion(particles, Node->children()[0]->dipol()->lastPt(), theMergePt))return ZERO;
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
  
  if (treefactory()->onlymulti()!=-1&&
      treefactory()->onlymulti()!=int(Node->xcomb()->mePartonData().size()-Node->nodeME()->projectorStage()))
    return ZERO;
  
  assert(Node->allAbove(mergePt()-0.1*GeV));
  

  
  if (!Born->xcomb()->willPassCuts()){
    return ZERO;
  }
  
  
  Energy startscale=CKKW_StartScale(Born);
  Energy projectedscale=startscale;
  fillHistory( startscale,  Born, Node);
  if (!fillProjector(projectedscale))return ZERO;
  Node->runningPt(projectedscale);
  double w1=history.back().weight;
  double w2=alphaReweight();
  double w3=pdfReweight();
  
  if(std::isnan(w1)){cout<<"\nhistory weight is nan";w1=0.;};
  if(std::isnan(w2)){cout<<"\nalphaReweight weight is nan";w1=0.;};
  if(std::isnan(w3)){cout<<"\npdfReweight weight is nan";w1=0.;};
  weight*=w1*w2*w3;
  if(weight==0.)return ZERO;
  bool maxMulti=Node->xcomb()->meMomenta().size() == maxLegsLO();
  Node->vetoPt((projected&&maxMulti)?mergePt():history.back().scale);
  CrossSection w4=TreedSigDR(startscale,Node,1.);
  if(std::isnan(double(w4/nanobarn))){cout<<"\nTreedSigDR is nan";w1=0.;};
  return weight*w4;
}

 

CrossSection Merger::MergingDSigDRVirtualStandard(NPtr Node ){
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
  Energy projectedscale=startscale;
  fillHistory( startscale,  Born, Node);
  //if (history.size()==1&&Node->children().size()!=0) {
  //  cout<<"\n1-->"<<startscale/GeV<<" "<<weight;
  //}
  if (!fillProjector(projectedscale))return ZERO;
  Node->runningPt(projectedscale);
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.)return ZERO;
  bool maxMulti=Node->xcomb()->meMomenta().size() == maxLegsNLO();
  Node->vetoPt((projected&&maxMulti)?mergePt():history.back().scale);
  
  CrossSection matrixElement=LoopdSigDR(startscale,Node);
  
  CrossSection Bornweight=Node->nodeME()->dSigHatDRB();
  
 
  
  CrossSection unlopsweight =(-sumpdfReweightUnlops()
                        -sumalphaReweightUnlops()
                        -sumfillHistoryUnlops())
                       *Bornweight
                       *SM().alphaS()/(2.*ThePEG::Constants::pi);


  return weight*
         as(startscale*xiRenSh)/SM().alphaS()*
         (matrixElement+unlopsweight);
}


CrossSection Merger::MergingDSigDRRealStandard(NPtr Node){
  bool allAbove=Node->allAbove(mergePt());
  if(!Node->allAbove(theIRSafePT))return ZERO;
  if (allAbove)return MergingDSigDRRealAllAbove(Node);
  if (UseRandom::rnd()<.5)
    return 2.*MergingDSigDRRealBelowSubReal( Node );
  return 2.*MergingDSigDRRealBelowSubInt( Node);
}

CrossSection Merger::MergingDSigDRRealAllAbove(NPtr Node){
  NPtr Born= Node-> getHistory(true,xiQSh);
  NPtr CLNode= Node->randomChild();
  bool inhist=CLNode->isInHistoryOf(Born);
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
  Energy projectedscale=startscale;
  fillHistory( startscale,  Born, CLNode);
  if (!fillProjector(projectedscale))return ZERO;
  Node->runningPt(projectedscale);
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.)return ZERO;
  bool maxMulti=CLNode->xcomb()->meMomenta().size() == maxLegsNLO();
  Node->vetoPt((projected&&maxMulti)?mergePt():history.back().scale);
  
  CrossSection res= weight*as(startscale*xiRenSh)/SM().alphaS()*
         (double)Node->children().size()*
         ((inhist?TreedSigDR(startscale,Node):ZERO)
          +CLNode->dipol()->dSigHatDR(sqr(startscale*xiFacME)));
    // cout<<"\nCLNode->dipol()->dSigHatDR(sqr(startscale*xiFacME))/nanobarn "<<CLNode->dipol()->dSigHatDR(sqr(startscale*xiFacME))/nanobarn;
  return res;
}

CrossSection Merger::MergingDSigDRRealBelowSubReal(NPtr Node){
  NPtrVec children=Node->children();
  Selector<NPtr> HistNodeSel;
  Energy minScale=generator()->maximumCMEnergy();
  for (NPtrVec::iterator child = children.begin();
       child != children.end(); child++){
    if ((*child)->dipol()->lastPt()<minScale) {
      minScale=(*child)->dipol()->lastPt();
      HistNodeSel.insert(1.,*child);
    }
  }
  NPtr HistNode=HistNodeSel.select(UseRandom::rnd());
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
  Energy projectedscale=startscale;
  fillHistory( startscale,  Born, HistNode);
  if (!fillProjector(projectedscale))return ZERO;
  Node->runningPt(projectedscale);
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.)return ZERO;
  bool maxMulti=HistNode->xcomb()->meMomenta().size() == maxLegsNLO();
  Node->vetoPt((projected&&maxMulti)?mergePt():history.back().scale);
  
  CrossSection sumPS=ZERO;
  for (NPtrVec::iterator child = children.begin();
       child != children.end(); child++){
    if ((*child)->allAbove(mergePt()))
      sumPS-=(*child)->calcPs(startscale*xiFacME);
  }
  CrossSection me=TreedSigDR(startscale,Node);
  //cout<<"\nSubReal "<<Node->miniPt()/GeV<<" "<<me<<" "<<sumPS<<" "<<me/sumPS;
  
  return weight*as(startscale*xiRenSh)/SM().alphaS()*
  (me-sumPS);
}



CrossSection Merger::MergingDSigDRRealBelowSubInt(NPtr Node){
  NPtr CLNode= Node->randomChild();
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
  Energy projectedscale=startscale;
  fillHistory( startscale,  Born, CLNode);
  if (!fillProjector(projectedscale))return ZERO;
  Node->runningPt(projectedscale);
  weight*=history.back().weight*alphaReweight()*pdfReweight();
  if(weight==0.)return ZERO;
  bool maxMulti=CLNode->xcomb()->meMomenta().size() == maxLegsNLO();
  Node->vetoPt((projected&&maxMulti)?mergePt():history.back().scale);
  
  pair<CrossSection,CrossSection> DipAndPs=
               CLNode->calcDipandPS(startscale*xiFacME);
  
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
  return res;
}

CrossSection Merger::LoopdSigDR(Energy startscale,NPtr Node){
    // The deephead should be calculated here.
  NPtr DeepHead=Node;//->deepHead();
  renormscale(startscale);
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

bool Merger::fillProjector(Energy& prerunning){
  // in the shower handler the scale is multiplied by xiQSh
  // so here we need to devide this
  
  double xiQShfactor=history.begin()->node->legsize()==N0()?xiQSh:1.;
  if(history.begin()->node->deepHead()->nodeME()->projectorStage() == 0){
    prerunning=(history.size()==1?1.:(1./xiQShfactor))*history.back().scale;
    return true;
  }
  for (History::iterator it=history.begin();it!=history.end();it++){
    if (projectorStage((*it).node)&&history.begin()->node->deepHead()->nodeME()->projectorStage() != 0){
      history.begin()->node->deepHead()->xcomb()->lastProjector((*it).node->xcomb());
      prerunning=(it==history.begin()?1.:(1./xiQShfactor))*(*it).scale;
      return true;
    }
  }
  return false;
}

double Merger::pdfReweight(){
  
  double res=1.;
  double facfactor=(history[0].node->legsize()==N0()?xiFacME:xiFacSh);
  for(int side=0;side!=2;side++){
    if(history[0].node->xcomb()->mePartonData()[side]->coloured()){
      History::iterator it=history.begin();
      for (;it+1!=history.end();it++){
        res *= pdfratio((*it).node, facfactor*(*it).scale,xiFacSh*((*it).node->dipol()->lastPt()), side);
	facfactor=xiFacSh;
      }
      res*=pdfratio(history.back().node,facfactor*history.back().scale ,history[0].scale*xiFacME, side);
    }
  }
  return res;
}

double Merger::alphaReweight(){
  double res=1.;
  Energy Q_R=(history[0].node->legsize()==N0()?xiRenME:xiRenSh)*history[0].scale;
  res *= pow(as(Q_R) / SM().alphaS(), history[0].node->nodeME()->orderInAlphaS());
  res *= pow(history[0].node->deepHead()->xcomb()->eventHandler().SM().alphaEMPtr()->value(history[0].node->nodeME()->factory()->scaleChoice()->renormalizationScaleQED())/ SM().alphaEMMZ(), history[0].node->nodeME()->orderInAlphaEW());
 

  if (!(history[0].node->children().empty())){ 
    res *=pow((theCMWScheme?(1.+((3.*(67./18.-1./6.*Constants::pi*Constants::pi)-5./9.*Nf(Q_R))*as(Q_R))/2./Constants::pi):1.),int(history[0].node->legsize()-N0()));
  }


 
  for (History::iterator it=history.begin();(it+1)!=history.end();it++){
    if ((*it).node->parent()){
      Energy q_i=xiRenSh* (*it).node->dipol()->lastPt();
      res *= as(q_i)/ SM().alphaS()
      *(theCMWScheme?(1.+((3.*(67./18.-1./6.*Constants::pi*Constants::pi)-5./9.*Nf(q_i))*as(q_i))/2./Constants::pi):1.);
    }
  }
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
      if (!dosudakov(Begin,scale, Begin->dipol()->lastPt(),sudakov0_n)){
        history.push_back(HistoryStep(Begin->parent(),0.,scale));
      }
      scale=Begin->dipol()->lastPt();
      history.push_back(HistoryStep(Begin->parent(),sudakov0_n,Begin->dipol()->lastPt()));
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
  History::iterator it=history.begin();
  
  
  
  
  for (;it!=history.end();it++){
    History::iterator ittmp=it;
    ittmp++;
    if((*it).node->xcomb()->mePartonData()[0]->coloured()&&(*it).node->nodeME()->lastX1()>0.99)return 0.;
    if((*it).node->xcomb()->mePartonData()[1]->coloured()&&(*it).node->nodeME()->lastX2()>0.99)return 0.;
    
    if (ittmp!=history.end()){
      if((*it).node->nodeME()->lastX1()<0.00001)return 0.;
      if((*it).node->nodeME()->lastX2()<0.00001)return 0.;
      
      if ((*it).node->dipol()->bornEmitter() == 0 ){
        res +=pdfUnlops((*it).node->nodeME()->lastParticles().first->dataPtr(),
                        (*it).node->nodeME()->lastPartons().first->dataPtr(),
                        (*it).node->xcomb()->partonBins().first->pdf(),
                        beam1Scale,
                        ((*it).node->dipol()->lastPt()),
                        (*it).node->nodeME()->lastX1(),
                        Nf(history[0].scale),
                        history[0].scale);
        beam1Scale=((*it).node->dipol()->lastPt())*xiFacSh;
      }
      if ((*it).node->dipol()->bornEmitter() == 1 ){
        res +=pdfUnlops((*it).node->nodeME()->lastParticles().second->dataPtr(),
                        (*it).node->nodeME()->lastPartons().second->dataPtr(),
                        (*it).node->xcomb()->partonBins().second->pdf(),
                        beam2Scale,
                        ((*it).node->dipol()->lastPt()),
                        (*it).node->nodeME()->lastX2(),
                        Nf(history[0].scale),
                        history[0].scale);
        beam2Scale=((*it).node->dipol()->lastPt())*xiFacSh;
      }
      if ((*it).node->dipol()->bornSpectator() == 0 &&(*it).node->dipol()->bornEmitter() >1){//
        res +=pdfUnlops((*it).node->nodeME()->lastParticles().first->dataPtr(),
                        (*it).node->nodeME()->lastPartons().first->dataPtr(),
                        (*it).node->xcomb()->partonBins().first->pdf(),
                        beam1Scale,
                        ((*it).node->dipol()->lastPt()),
                        (*it).node->nodeME()->lastX1(),
                        Nf(history[0].scale),
                        history[0].scale);
          //pdfratio((*it).node, beam1Scale, sqrt((*it).node->dipol()->lastPt()), 1);
        beam1Scale=((*it).node->dipol()->lastPt())*xiFacSh;
      }
      if ((*it).node->dipol()->bornSpectator() == 1 &&(*it).node->dipol()->bornEmitter() >1){//
        res +=pdfUnlops((*it).node->nodeME()->lastParticles().second->dataPtr(),
                        (*it).node->nodeME()->lastPartons().second->dataPtr(),
                        (*it).node->xcomb()->partonBins().second->pdf(),
                        beam2Scale,
                        ((*it).node->dipol()->lastPt()),
                        (*it).node->nodeME()->lastX2(),
                        Nf(history[0].scale),
                        history[0].scale);
          //pdfratio((*it).node, beam2Scale , sqrt((*it).node->dipol()->lastPt()), 2);
        beam2Scale=((*it).node->dipol()->lastPt())*xiFacSh;
      }
    }
  }
  if (history[0].node->deepHead()->xcomb()->mePartonData()[0]->coloured()){
    res +=pdfUnlops(history.back().node->nodeME()->lastParticles().first->dataPtr(),
                    history.back().node->nodeME()->lastPartons().first->dataPtr(),
                    history.back().node->xcomb()->partonBins().first->pdf(),
                    beam1Scale,
                    history[0].scale*xiFacME,
                    (history.back()).node->nodeME()->lastX1(),
                    Nf(history[0].scale),
                    history[0].scale);
    
  }
  if (history[0].node->deepHead()->xcomb()->mePartonData()[1]->coloured()) {
    res +=pdfUnlops(history.back().node->nodeME()->lastParticles().second->dataPtr(),
                    history.back().node->nodeME()->lastPartons().second->dataPtr(),
                    history.back().node->xcomb()->partonBins().second->pdf(),
                    beam2Scale,
                    history[0].scale*xiFacME,
                    (history.back()).node->nodeME()->lastX2(),
                    Nf(history[0].scale),
                    history[0].scale);
    
      //history[0].node->deepHead()->nodeME()->pdf2(sqr(beam2Scale))/history[0].node->deepHead()->nodeME()->pdf2(sqr(history[0].scale));
  }
  return res;
}






double Merger::pdfUnlops(tcPDPtr particle,tcPDPtr parton,tcPDFPtr pdf,Energy  running, Energy next,double x,int nlp,Energy fixedScale)  {
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
  for (History::iterator it=history.begin();(it+1)!=history.end();it++){
    assert((*it).node->parent());
    res +=alphasUnlops((*it).node->dipol()->lastPt()*xiRenSh ,history[0].scale);
  }
  return res;
}

double Merger::sumfillHistoryUnlops(){
  double res=0.;
  double xiQShfactor=history.begin()->node->legsize()==N0()?xiQSh:1.;
  for (History::iterator it = history.begin(); (it+1) != history.end();it++){
    doUNLOPS((*it).node,(it == history.begin()?xiQShfactor:1.)*(*it).scale , (*it).node->dipol()->lastPt() , history[0].scale, res);
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
  
  
  //cout<<"\nrunning "<<Node->runningPt()/GeV<<" "<<Node->legsize(); 
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











  //----------------------------------Reviewed--------------------------------------//




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
  for ( list<DipoleChain>::iterator chain = DSH()->eventRecord().chains().begin() ;
       chain != DSH()->eventRecord().chains().end() ; chain++ ) {
    for ( list<Dipole>::iterator dip = (*chain).dipoles().begin() ; dip != (*chain).dipoles().end() ; ++dip ) {
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
  for ( list<DipoleChain>::iterator chain = DSH()->eventRecord().chains().begin() ;
       chain != DSH()->eventRecord().chains().end() ; chain++ ) {
    for ( list<Dipole>::iterator dip = (*chain).dipoles().begin() ; dip != (*chain).dipoles().end() ; ++dip ) {
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
  for ( ParticleVector::const_iterator it = vecfirst.begin() ; it != vecfirst.end() ; it++ )
    Born->xcomb()->subProcess()->incoming().first->abandonChild(*it);
  ParticleVector vecsecond = Born->xcomb()->subProcess()->incoming().second->children();
  for ( ParticleVector::const_iterator it = vecsecond.begin() ; it != vecsecond.end() ; it++ )
    Born->xcomb()->subProcess()->incoming().second->abandonChild(*it);
  Born->xcomb()->subProcess(SubProPtr());
}

double Merger::singlesudakov(list<Dipole>::iterator dip ,Energy next,Energy running,pair<bool,bool> conf ){
  
  double res=1.;
  tPPtr emitter = dip->emitter(conf);
  tPPtr spectator = dip->spectator(conf);
  DipoleSplittingInfo candidate((*dip).index(conf),conf,(*dip).emitterX(conf),(*dip).spectatorX(conf),emitter,spectator);
  
  
  if ( DSH()->generators().find(candidate.index()) == DSH()->generators().end() ) DSH()->getGenerators(candidate.index());
  
  pair<GeneratorMap2::iterator,GeneratorMap2::iterator> gens = DSH()->generators().equal_range(candidate.index());
  
  for ( GeneratorMap2::iterator gen = gens.first; gen != gens.second; ++gen ) {
    if ( !(gen->first == candidate.index()) )
      continue;
    
    Energy dScale =	gen->second->splittingKinematics()->dipoleScale(emitter->momentum(),spectator->momentum());
    candidate.scale(dScale);
    candidate.continuesEvolving();
    Energy ptMax=(*gen).second->splittingKinematics()->ptMax(candidate.scale(),candidate.emitterX(), candidate.spectatorX(),
                                                             candidate.index(),*gen->second->splittingKernel());
    
    candidate.hardPt(min(running,ptMax));
    
    if (candidate.hardPt()>next){
      res*=gen->second->sudakov(candidate,next);
    }
  }
  
  return res;
}


double Merger::singleUNLOPS(list<Dipole>::iterator dip ,Energy next,Energy running,Energy fixedScale,pair<bool,bool> conf ){
  
  double res=0.;
  tPPtr emitter = dip->emitter(conf);
  tPPtr spectator = dip->spectator(conf);
  DipoleSplittingInfo candidate((*dip).index(conf),conf,(*dip).emitterX(conf),(*dip).spectatorX(conf),emitter,spectator);
  
  if ( DSH()->generators().find(candidate.index()) == DSH()->generators().end() ) DSH()->getGenerators(candidate.index());
  
  pair<GeneratorMap2::iterator,GeneratorMap2::iterator> gens = DSH()->generators().equal_range(candidate.index());
  
  for ( GeneratorMap2::iterator gen = gens.first; gen != gens.second; ++gen ) {
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
  
  theFirstNodeMap[me]->firstgenerateKinematics(r, 0,me->lastXCombPtr()->lastSHat());
  
    if (theFirstNodeMap[me]->cutStage()==0 ){
      
      bool inAlphaPS=false;
      NPtrVec children=theFirstNodeMap[me]->children();
      for (NPtrVec::iterator child = children.begin();
           child != children.end(); child++){
       
        treefactory()->setAlphaParameter(theGamma);
        inAlphaPS|=theGamma!=1.?(*child)->dipol()->aboveAlpha():false;
        treefactory()->setAlphaParameter(1.);
      }
      
      SafeClusterMap temp=theFirstNodeMap[me]->clusterSafe();
      for(SafeClusterMap::iterator
          it=temp.begin();
          it!=temp.end();++it){
        if (!it->second.first&&!inAlphaPS)return false;
    		}
    }
    if (theFirstNodeMap[me]->cutStage()==1 ){
      SafeClusterMap temp=theFirstNodeMap[me]->clusterSafe();
      for(SafeClusterMap::iterator
          it=temp.begin();
          it!=temp.end();++it){
      		if (!it->second.first && !it->second.second)return false;
      }
    }
  return true;
  
}
bool Merger::calculateInNode() const{
    return theCalculateInNode;
}


void Merger::fillProjectors(Ptr<MatchboxMEBase>::ptr me){
  for (unsigned int i = 0; i < (theFirstNodeMap[me]->Projector()).size(); ++i) {
    me->lastXCombPtr()->projectors().insert(
                                        (theFirstNodeMap[me]->Projector())[i].first,
                                        (theFirstNodeMap[me]->Projector())[i].second->xcomb());
  }
}
pair<bool,bool> Merger::clusterSafe(Ptr<MatchboxMEBase>::ptr me,int emit,int emis,int spec){
  return theFirstNodeMap[me]->clusterSafe().find(make_pair(make_pair(emit,emis),spec))->second;
  
}








bool Merger::matrixElementRegion(PVector particles,Energy winnerScale,Energy cutscale){
  
    //cout<<"\nparticles s"<<particles.size()<<" "<<particles[0]<<" "<<particles[1]<<flush;
  /*
   if (defMERegionByJetAlg && !particles[0]->coloured()&& !particles[1]->coloured()) {
   assert(false);
   vector<fastjet::PseudoJet> input_particles;
   for(size_t em=2; em < particles.size();em++){
   input_particles.push_back(fastjet::PseudoJet(particles[em]->momentum().x()/GeV,
   particles[em]->momentum().y()/GeV,
   particles[em]->momentum().z()/GeV,
   particles[em]->momentum().e()/GeV));
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
   if (particles[em]->coloured())
   input_particles.push_back(fastjet::PseudoJet(particles[em]->momentum().x()/GeV,
   particles[em]->momentum().y()/GeV,
   particles[em]->momentum().z()/GeV,
   particles[em]->momentum().e()/GeV));
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
  
    //FF
  for(size_t em=2; em < particles.size();em++){
    if (!particles[em]->coloured()) continue;
    for(size_t emm=2; emm < particles.size();emm++){
      if (!particles[emm]->coloured()) continue;
      if (em==emm) continue;
      for(size_t spe=2; spe < particles.size();spe++){
        if (!particles[spe]->coloured()) continue;
        if (em==spe||emm==spe) continue;
        if (!(particles[em]->id()==-particles[emm]->id()||particles[emm]->id()>6))continue;
          //  assert(false);
        Lorentz5Momentum emittermom = particles[em]->momentum();
        Lorentz5Momentum emissionmom = particles[emm]->momentum();
        Lorentz5Momentum spectatormom = particles[spe]->momentum();
        Energy pt=0*GeV;
        if (emittermom.m()==0.001*GeV&&emissionmom.m()==0.001*GeV&&spectatormom.m()==0.001*GeV) {
          pt=FFLTK->lastPt(emittermom,emissionmom,spectatormom);
        }else{
          pt=FFMTK->lastPt(emittermom,emissionmom,spectatormom);
        }
        
          //cout<<"\npt "<<pt/GeV<<" "<<winnerScale/GeV;
        if (abs(pt-winnerScale)<0.001*GeV) {
          foundwinnerpt=true;
        }
          //          if(scale * sqrt(y*z*(1.-z))<optVeto&&winnerScale>optVeto)cout<<"\nFF "<<(scale * sqrt(y*z*(1.-z))/GeV);
        ptx =min(ptx,pt);
      }
    }
  }
  
    //FI
  for(size_t spe=0; spe < 2;spe++){
    if (!particles[spe]->coloured()) continue;
    for(size_t em=2; em < particles.size();em++){
      if (!particles[em]->coloured()) continue;
      for(size_t emm=2; emm < particles.size();emm++){
        if (!particles[emm]->coloured()) continue;
        if (em==emm) continue;
        if (!(particles[em]->id()==-particles[emm]->id()||particles[emm]->id()>6))continue;
          //  assert(false);
        Lorentz5Momentum emittermom = particles[em]->momentum();
        Lorentz5Momentum emissionmom = particles[emm]->momentum();
        Lorentz5Momentum spectatormom = particles[spe]->momentum();
        Energy pt=0*GeV;
        if (emittermom.m()==0.001*GeV&&emissionmom.m()==0.001*GeV&&spectatormom.m()==0.001*GeV) {
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
  for(size_t em=0; em < 2;em++){
    if (!particles[em]->coloured()) continue;
    for(size_t emm=2; emm < particles.size();emm++){
      if (!particles[emm]->coloured()) continue;
      for(size_t spe=2; spe < particles.size();spe++){
        if (!particles[spe]->coloured()) continue;
        
        if (emm==spe) continue;
        if (!(particles[em]->id()>6|| particles[em]->id()==particles[emm]->id() ||particles[emm]->id()>6))continue;
          // assert(false);
        Lorentz5Momentum emittermom = particles[em]->momentum();
        Lorentz5Momentum emissionmom = particles[emm]->momentum();
        Lorentz5Momentum spectatormom = particles[spe]->momentum();
        
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
  for(size_t em=0; em < 2;em++){
    if (!particles[em]->coloured()) continue;
    for(size_t spe=0; spe < 2;spe++){
      if (!particles[spe]->coloured()) continue;
      for(size_t emm=2; emm < particles.size();emm++){
        if (!particles[emm]->coloured()) continue;
        if (em==spe) continue;
        if (!(particles[em]->id()>6||
              particles[em]->id()==particles[emm]->id() ||
              particles[emm]->id()>6))continue;
          //assert(false);
        Lorentz5Momentum emittermom = particles[em]->momentum();
        Lorentz5Momentum emissionmom = particles[emm]->momentum();
        Lorentz5Momentum spectatormom = particles[spe]->momentum();
        Energy  pt=IILTK->lastPt(emittermom,emissionmom,spectatormom);
       if (abs(pt-winnerScale)<0.01*GeV) {
          foundwinnerpt=true;
        }
        ptx =min(ptx, pt);
      }
    }
  }
  
    // cout<<"\n"<<cutscale/GeV<< " "<<foundwinnerpt<<" "<<ptx/GeV;
  
  
  
  if(!foundwinnerpt){
    cout<<"\nWarning: could not find winner with pt: "<<winnerScale/GeV;
    for(size_t emm=0; emm < particles.size();emm++){
    cout<<"\n"<<particles[emm]->id()<<" "<<particles[emm]->momentum()/GeV<<" "<<particles[emm]->momentum().m()/GeV<<flush;
    }
  }
  
  return (ptx>cutscale) ;
  
}





int Merger::M()const{return theTreeFactory->M();}

int Merger::N()const{return theTreeFactory->N();}











  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Merger::persistentOutput(PersistentOStream & os) const {
  
  os << minusL<<  Unlopsweights<<
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
  is >> minusL>>  Unlopsweights>>
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
  interfaceminusL ("minusL","",&Merger::minusL, false, false, false);
  static SwitchOption interfaceminusLYes
  (interfaceminusL,"Yes","",true);
  static SwitchOption interfaceminusLNo
  (interfaceminusL,"No","",false);
  
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
















