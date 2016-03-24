  // -*- C++ -*-
  //
  // Merging.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2007 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the Merging class.
  //

#include "Merging.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "Herwig/DipoleShower/DipoleShowerHandler.h"
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

#include "Herwig/MatrixElement/Matchbox/Mergeing/MergeFactory.h"


using namespace Herwig;



Merging::Merging()
: HandlerBase() {
  StartingBorn0=CNPtr();
  StartingBorn1=CNPtr();
  StartingBorn2=CNPtr();
  CNPtr CalcBorn0=CNPtr();
  CNPtr CalcBorn1=CNPtr();
  CNPtr CalcBorn2=CNPtr();
  theNf=5;
  minusL=false;
  Unlopsweights=true;
  theKImproved=true;
  MergingScale=20.*GeV;
}


Merging::~Merging() {}

IBPtr Merging::clone() const {
  return new_ptr(*this);
}

IBPtr Merging::fullclone() const {
  return new_ptr(*this);
}




bool Merging::sudakov(CNPtr Born, Energy & running, Energy next) {
  if(running<next)return false;
  Born->VetoedShower(false);
  Born->deepHead()->VetoedShower(false);
  tSubProPtr sub = Born->xcomb()->construct();
  theDipoleShowerHandler->resetPDFs(make_pair(Born->xcomb()->partonBins().first->pdf(), Born->xcomb()->partonBins().second->pdf()), Born->xcomb()->partonBins());
  theDipoleShowerHandler->setCurrentHandler();
  theDipoleShowerHandler->generator()->currentEventHandler(Born->deepHead()->xcomb()->eventHandlerPtr());
  theDipoleShowerHandler->remnantDecayer()->setHadronContent(Born->deepHead()->xcomb()->lastParticles());
  theDipoleShowerHandler->eventRecord().clear();
  theDipoleShowerHandler->eventRecord().prepare(sub, dynamic_ptr_cast<tStdXCombPtr>(Born->xcomb()), theDipoleShowerHandler->pdfs(), Born->deepHead()->xcomb()->lastParticles());
  Born->xcomb()->lastCentralScale(running*running);
  theDipoleShowerHandler->hardScales(running*running);
  
    //didRadiate = false;//TODO?
    //didRealign = false;//TODO?
  unsigned int nEmitted = 0;
  theDipoleShowerHandler->setNEmissions(1);
  try {
    theDipoleShowerHandler->doCascade(nEmitted);
  } catch (...) {
    theDipoleShowerHandler->setNEmissions(1);
    return false;
  }
  theDipoleShowerHandler->setNEmissions(1);
  if ( nEmitted == 1 ) {
    Energy after = 0. * GeV;
    for ( list<DipoleChain>::iterator chain = theDipoleShowerHandler->eventRecord().chains().begin() ; chain != theDipoleShowerHandler->eventRecord().chains().end() ; chain++ ) {
      for ( list<Dipole>::iterator dip = (*chain).dipoles().begin() ; dip != (*chain).dipoles().end() ; ++dip ) {
        after = max(dip->leftScale(), max(dip->rightScale(), after));
      }
    }
    for ( list<DipoleChain>::iterator chain = theDipoleShowerHandler->eventRecord().doneChains().begin() ; chain != theDipoleShowerHandler->eventRecord().doneChains().end() ; chain++ ) {
      for ( list<Dipole>::iterator dip = (*chain).dipoles().begin() ; dip != (*chain).dipoles().end() ; ++dip ) {
        after = max(dip->leftScale(), max(dip->rightScale(), after));
      }
    }
    if ( after > next ) {
      
      running = -1. * GeV;
      return false;
    }
  }
  running = next;
  return true;
}

double Merging::singlesudakov(list<Dipole>::iterator dip ,Energy next,Energy running,pair<bool,bool> conf ,bool fast){
  DipoleSplittingInfo candidate;
  double tmp=1.;
  tPPtr emitter = dip->emitter(conf);
  tPPtr spectator = dip->spectator(conf);
  candidate.index((*dip).index(conf));
  candidate.configuration(conf);
  candidate.emitterX((*dip).emitterX(conf));
  candidate.spectatorX((*dip).spectatorX(conf));
  candidate.emitter(emitter);
  candidate.spectator(spectator);
  if ( theDipoleShowerHandler->generators().find(candidate.index()) == theDipoleShowerHandler->generators().end() ) theDipoleShowerHandler->getGenerators(candidate.index());
  
  pair<GeneratorMap2::iterator,GeneratorMap2::iterator> gens = theDipoleShowerHandler->generators().equal_range(candidate.index());
  
  for ( GeneratorMap2::iterator gen = gens.first; gen != gens.second; ++gen ) {
      //TODO
      //assert(false);
       if ( !(gen->first == candidate.index()) )
        continue;
    
      //TODO!!!!
    Energy dScale =	gen->second->splittingKinematics()->dipoleScale(emitter->momentum(),spectator->momentum());
    assert (! isnan(dScale/GeV ) );
    candidate.scale(dScale);
    candidate.continuesEvolving();
    Energy ptMax=(*gen).second->splittingKinematics()->ptMax(candidate.scale(),candidate.emitterX(), candidate.spectatorX(),
                                                             candidate.index(),*gen->second->splittingKernel());
    
    candidate.hardPt(min(running,ptMax));
    
    if (candidate.hardPt()>next){
      tmp*=gen->second->sudakov(candidate,next,fast);
    }
  }
  
	 return tmp;
}

double Merging::singleUNLOPS(list<Dipole>::iterator dip ,Energy next,Energy running,Energy fixedScale,pair<bool,bool> conf ){
  DipoleSplittingInfo candidate;
  double res=0.;
  tPPtr emitter = dip->emitter(conf);
  tPPtr spectator = dip->spectator(conf);
  candidate.index((*dip).index(conf));
  candidate.configuration(conf);
  candidate.emitterX((*dip).emitterX(conf));
  candidate.spectatorX((*dip).spectatorX(conf));
  candidate.emitter(emitter);
  candidate.spectator(spectator);
  if ( theDipoleShowerHandler->generators().find(candidate.index()) == theDipoleShowerHandler->generators().end() ) theDipoleShowerHandler->getGenerators(candidate.index());
  
  pair<GeneratorMap2::iterator,GeneratorMap2::iterator> gens = theDipoleShowerHandler->generators().equal_range(candidate.index());
  
  for ( GeneratorMap2::iterator gen = gens.first; gen != gens.second; ++gen ) {
    if ( !(gen->first == candidate.index()) )
      continue;
    Energy dScale =	gen->second->splittingKinematics()->dipoleScale(emitter->momentum(),spectator->momentum());
    assert (! isnan(dScale/GeV ) );
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


bool Merging::dosudakov(Energy & running, Energy next, double& sudakov0_n,bool fast) {
  double suda=1.;
  for ( list<DipoleChain>::iterator chain = theDipoleShowerHandler->eventRecord().chains().begin() ; chain != theDipoleShowerHandler->eventRecord().chains().end() ; chain++ ) {
    for ( list<Dipole>::iterator dip = (*chain).dipoles().begin() ; dip != (*chain).dipoles().end() ; ++dip ) {
      suda*=singlesudakov( dip, next,running,make_pair(true,false), fast );
      suda*=singlesudakov( dip, next,running,make_pair(false,true), fast );
    }
  }
  running=next;
  sudakov0_n*=suda;
  if (suda==0.0)
    return false;
  return true;
}

bool Merging::doUNLOPS(Energy  running, Energy next,Energy fixedScale, double& UNLOPS) {
  double unlops=0.;
  
  for ( list<DipoleChain>::iterator chain = theDipoleShowerHandler->eventRecord().chains().begin() ; chain != theDipoleShowerHandler->eventRecord().chains().end() ; chain++ ) {
    for ( list<Dipole>::iterator dip = (*chain).dipoles().begin() ; dip != (*chain).dipoles().end() ; ++dip ) {
      double tmp =singleUNLOPS( dip, next,running,fixedScale,make_pair(true,false) );
      unlops+=tmp;
      tmp=singleUNLOPS( dip, next,running,fixedScale,make_pair(false,true) );
      unlops+=tmp;
    }
  }
  UNLOPS+=unlops;
  return true;
}


double Merging::pdfUnlops(tcPDPtr particle,tcPDPtr parton,tcPDFPtr pdf,Energy  running, Energy next,double x,int nlp,Energy fixedScale)  {
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
      
      restmp += ( (11./6.) * CA - (1./3.)*theNf + 2.*CA*log(1.-x) ) *PDFxparton;
      if ( z > x ) {
        restmp += 2. * CA * ( PDFxByzgluon - z*PDFxparton ) / (z*(1.-z));
        restmp += 2.* CA *( (1.-z)/z - 1. + z*(1.-z) ) * PDFxByzgluon / z;
      }
      
    }
    if (PDFxparton<1e-8)restmp= 0.;
    res+=1*restmp*log(sqr(running/next))/PDFxparton*mapz;
    
  }
  return res/number; }








bool Merging::dosudakovold(CNPtr Born, Energy & running, Energy next, double& sudakov0_n) {
  double suda=0.0;
  Energy scale=running;
  cleanup(Born);
  double sudakovtries=5.;
  for (int step=0;step<sudakovtries;step++){
    if (sudakov(Born,scale, next))
      suda+=1./sudakovtries;
    cleanup(Born);
    scale=running;
  }
  running=next;
  sudakov0_n*=suda;
  cleanup(Born);
  if (suda==0.0)
    return false;
  return true;
}

void Merging::cleanup(CNPtr Born) {
  theDipoleShowerHandler->eventRecord().clear();
  if(!Born->xcomb()->subProcess())return;
  ParticleVector vecfirst = Born->xcomb()->subProcess()->incoming().first->children();
  for ( ParticleVector::const_iterator it = vecfirst.begin() ; it != vecfirst.end() ; it++ )
    Born->xcomb()->subProcess()->incoming().first->abandonChild(*it);
  ParticleVector vecsecond = Born->xcomb()->subProcess()->incoming().second->children();
  for ( ParticleVector::const_iterator it = vecsecond.begin() ; it != vecsecond.end() ; it++ )
    Born->xcomb()->subProcess()->incoming().second->abandonChild(*it);
  Born->xcomb()->subProcess(SubProPtr());
}

double Merging::pdfratio(CNPtr  Born,Energy & nominator_scale, Energy denominator_scale,int side){
  StdXCombPtr bXc = Born->xcomb();
  
  
    //  Ptr<StandardEventHandler>::ptr eh =
    //dynamic_ptr_cast<Ptr<StandardEventHandler>::ptr>(theDipoleShowerHandler->eventHandler());
  
    //if (eh->didEstimate()) {
    //cout<<"\n---> den "<<denominator_scale/GeV<<" nom "<<nominator_scale/GeV;
    //}
  
  if (side==1){
    if (!bXc->mePartonData()[0]->coloured())return 1.;
    double from=Born->nodeME()->pdf1(sqr(denominator_scale));
    double to=Born->nodeME()->pdf1(sqr( nominator_scale ));
    if ( to < 1e-8||from < 1e-8)return 0.;
    return to/from;
  }
  if (side==2){
    if (!bXc->mePartonData()[1]->coloured())return 1.;
    double from=Born->nodeME()->pdf2(sqr(denominator_scale));
    double to=Born->nodeME()->pdf2( sqr( nominator_scale ));
    if ( to < 1e-8||from < 1e-8)return 0.;
    return to/from;
  }
  return 0.0;
}

bool Merging::projectorStage(CNPtr  Born){
  Born->deepHead()->nodeME()->mePartonData().size();
  return	(Born->deepHead()->nodeME()->projectorStage() == int((Born->deepHead()->nodeME()->mePartonData().size() - Born->nodeME()->mePartonData().size())));
}


Energy Merging::CKKW_StartScale(CNPtr Born){
  Energy res=generator()->maximumCMEnergy();;
  if(!Born->children().empty()){
      // assert(false);
    for (size_t i=0;i<Born->nodeME()->mePartonData().size();i++){
      if (!Born->nodeME()->mePartonData()[i]->coloured())continue;
      for (size_t j=i+1;j<Born->nodeME()->mePartonData().size();j++){
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


void Merging::CKKW_PrepareSudakov(CNPtr Born,Energy running){
  cleanup(Born);
  tSubProPtr sub = Born->xcomb()->construct();
  theDipoleShowerHandler->resetPDFs(make_pair(Born->xcomb()->partonBins().first->pdf(), Born->xcomb()->partonBins().second->pdf()), Born->xcomb()->partonBins());
  theDipoleShowerHandler->setCurrentHandler();

  theDipoleShowerHandler->currentHandler()->generator()->currentEventHandler(Born->deepHead()->xcomb()->eventHandlerPtr());

  theDipoleShowerHandler->currentHandler()->remnantDecayer()->setHadronContent(Born->deepHead()->xcomb()->lastParticles());
  theDipoleShowerHandler->eventRecord().clear();
  theDipoleShowerHandler->eventRecord().prepare(sub, dynamic_ptr_cast<tStdXCombPtr>(Born->xcomb()), theDipoleShowerHandler->pdfs(), Born->deepHead()->xcomb()->lastParticles());
  theDipoleShowerHandler->hardScales(running*running);
}




double Merging::reweightCKKWBorn3(CNPtr CalcBorn2,bool fast){
  
  
  if(fast||!StartingBorn0){
    if( CalcBorn2->xcomb()->meMomenta().size()-2 == theMaxLegsLO){
      double anteil0=3.;
      double anteil1=1.;
      double anteil2=1.;
      double anteilsum=(anteil0+anteil1+anteil2);
      double stage=UseRandom::rnd()*anteilsum;//needs random number from sampler
      if (stage<anteil0) {
          //         Stage  Multi
        projectorWeight2Born0=0.;
        projectorWeight1Born1=0.;
        projectorWeight2Born1=0.;
        projectorWeight0Born2=anteilsum/anteil0;
        projectorWeight1Born2=0.;
        
        projectorWeight1Real1=0.;
        projectorWeight2Real1=0.;
        projectorWeight0Real2=anteilsum/anteil0;;//Theta<
        projectorWeight1Real2=0.;
        projectorWeight2Real2=0.;
        
        CalcBorn2->nodeME()->projectorStage(0);
      }else if (stage<(anteil0+anteil1)) {

        projectorWeight2Born0=0.;
        projectorWeight1Born1=anteilsum/anteil1;
        projectorWeight2Born1=0.;
        projectorWeight0Born2=0.;
        projectorWeight1Born2=-anteilsum/anteil1;
        
        projectorWeight1Real1=anteilsum/anteil1;;//Theta<  theta1=false;
        projectorWeight2Real1=0.;
        projectorWeight0Real2=0.;
        projectorWeight1Real2=anteilsum/anteil1;//Theta>  theta2=true;
        projectorWeight2Real2=0.;
        
        CalcBorn2->nodeME()->projectorStage(1);
      }else if (stage<=anteilsum) {
        
        projectorWeight2Born0=anteilsum/anteil2;;
        projectorWeight1Born1=0.;
        projectorWeight2Born1=-anteilsum/anteil2;
        projectorWeight0Born2=0.;
        projectorWeight1Born2=0.;
        
        projectorWeight1Real1=0.;
        projectorWeight2Real1=anteilsum/anteil2;//Theta>  theta1=true;
        projectorWeight0Real2=0.;
        projectorWeight1Real2=0.;
        projectorWeight2Real2=-anteilsum/anteil2;//Theta> theta2=false;
        
        CalcBorn2->nodeME()->projectorStage(2);
      }else{
        assert(false);
      }

    }
    else{
      assert(false);

      
      projectorWeight2Born0=1.;
      projectorWeight1Born1=0.;
      projectorWeight2Born1=-1.;
      projectorWeight0Born2=0.;
      projectorWeight1Born2=0.;
      
      projectorWeight1Real1=0.;
      projectorWeight2Real1=1.;//Theta>  theta1=true;
      projectorWeight0Real2=0.;
      projectorWeight1Real2=0.;
      projectorWeight2Real2=-1.;//Theta>   theta2=false;
      CalcBorn2->nodeME()->projectorStage(2);
    }
  }
  
  if (fast) {
    assert(!StartingBorn0);
    CalcBorn0=CNPtr();
    CalcBorn1=CNPtr();
  }
  
  double weightB0K0 =projectorWeight2Born0;  
  double weightB1K1 =projectorWeight1Born1;
  double weightB1K0 =projectorWeight2Born1;
  double weightB2K2 =projectorWeight0Born2;
  double weightB2K1 =projectorWeight1Born2;
  
  double weightR1K1 =projectorWeight1Real1;
  double weightR1K0 =projectorWeight2Real1;
  double weightR2K2 =projectorWeight0Real2;
  double weightR2K1 =projectorWeight1Real2;
  double weightR2K0 =projectorWeight2Real2;
  double headR1=1.;
  double headR2=1.;
  
  bool theta1=true;
  bool theta2=true;
  
  bool safetheta1=true;
  bool safetheta2=true;
  
  bool inhist=false;
  bool inhist2=false;
  
  CNPtr Born2;
  CNPtr Born1;
  CNPtr Born0;
  assert(!CalcBorn2->children().empty());
  
    // weightB0K0 =0.;
    // weightB1K1 =0.;
    // weightB1K0 =0.;
    // weightB2K2 =0.;
    // weightB2K1 =0.;
  
   weightR1K1 =0.;
   weightR1K0 =0.;
   weightR2K2 =0.;
   weightR2K1 =0.;
   weightR2K0 =0.;
  
  
  
  

  int randomIndex = 0;
  if(!CalcBorn1){
      //Set up the history by random choice
    randomIndex = (int)(UseRandom::rnd() *  CalcBorn2->children().size());
    CalcBorn1=CalcBorn2->children()[randomIndex];
    randomIndex = (int)(UseRandom::rnd() *  CalcBorn1->children().size());
    CalcBorn0=CalcBorn1->children()[randomIndex];
  }
  assert(CalcBorn1&&CalcBorn0);

  vector<Ptr<ClusterNode>::ptr> Children;
  
  Children = CalcBorn2->children();
  for (vector<Ptr<ClusterNode>::ptr>::iterator it = Children.begin(); it != Children.end(); it++) {
    if((*it)->dipol()->lastPt()<CalcBorn2->mergePt()){
      
      weightB2K2 =0.;
      weightB2K1 =0.;
      weightR2K1 =0.;
      
      theta2=false;
      if((*it)->dipol()->lastPt()<1.*GeV)safetheta2=false;
    }
  }
  Children = CalcBorn1->children();
  for (vector<Ptr<ClusterNode>::ptr>::iterator it = Children.begin(); it != Children.end(); it++) {
    if((*it)->dipol()->lastPt()<CalcBorn2->mergePt()){

       weightB1K1 =0.;
       weightB1K0 =0.;

       weightR1K0 =0.;
       weightR2K2 =0.;
       weightR2K1 =0.;
       weightR2K0 =0.;
      theta1=false;
      if((*it)->dipol()->lastPt()<1.*GeV)safetheta1=false;
    }
  }

  Children = CalcBorn0->children();
  for (vector<Ptr<ClusterNode>::ptr>::iterator it = Children.begin(); it != Children.end(); it++) {
    if((*it)->dipol()->lastPt()<CalcBorn2->mergePt()){
      weightB0K0 =0.;

      weightR1K1 =0.;
      weightR1K0 =0.;
        //weightR2K2 =0.;
        //weightR2K1 =0.;
        //weightR2K0 =0.;
      assert(false);//More than 2 NLO
    }
  }
  

  
  
  
  
  
  if (StartingBorn0) {
    assert(!fast&&StartingBorn1&&StartingBorn2);
    Born2=StartingBorn2;
    Born1=StartingBorn1;
    Born0=StartingBorn0;
    while (Born2->parent()) {
      inhist2|=(Born2==CalcBorn1);
      Born2=Born2->parent();
    }
    Born2=StartingBorn2;
    
    while (Born1->parent()) {
      inhist|=(Born1==CalcBorn0);
      Born1=Born1->parent();
    }
    Born0 = StartingBorn0;
    Born1 = StartingBorn1;
    Born2 = StartingBorn2;
  }else{
    Born2 = CalcBorn2->getLongestHistory_simple(true,xiQSh);
    Born1 = CalcBorn1->getLongestHistory_simple(false,xiQSh);
    Born0 = CalcBorn0->getLongestHistory_simple(false,xiQSh);
    StartingBorn2=Born2;
    StartingBorn1=Born1;
    StartingBorn0=Born0;
    while (Born2->parent()) {
      inhist2|=(Born2==CalcBorn1);
      Born2=Born2->parent();
    }
    Born2=StartingBorn2;
    while (Born1->parent()) {
     inhist|=(Born1==CalcBorn0);
     Born1=Born1->parent();
    }
    Born1=StartingBorn1;
    if(!fast){
      StartingBorn0=CNPtr();
      StartingBorn1=CNPtr();
      StartingBorn2=CNPtr();
    }
  }

  
  
  
  
   /**
    *  Setup the combination factors.
    *
    *
    *
    **/
  
  weightB0K0*=1.* CalcBorn1->children().size()/CalcBorn0->numberOfSplittings()*
                  CalcBorn2->children().size()/CalcBorn1->numberOfSplittings()*
                  CalcBorn1->dipol()->jacobianMerging(CalcBorn1->xcomb()->lastSHat(),
                                                      CalcBorn2->xcomb()->lastSHat(),
                                                      CalcBorn1->xcomb()->meMomenta().size())*
                  CalcBorn0->dipol()->jacobianMerging(CalcBorn0->xcomb()->lastSHat(),
                                                      CalcBorn1->xcomb()->lastSHat(),
                                                      CalcBorn0->xcomb()->meMomenta().size());
  
  
  
  if (CalcBorn1==Born1) { //Noclustering found --> finite Born: Do not calculate the real.(But the Dipole!!! (TODO: headR1))
      double factor=1.*CalcBorn2->children().size()/CalcBorn1->numberOfSplittings()*
                      CalcBorn1->dipol()->jacobianMerging(CalcBorn1->xcomb()->lastSHat(),
                                                          CalcBorn2->xcomb()->lastSHat(),
                                                          CalcBorn1->xcomb()->meMomenta().size());;
    weightB1K1 *=factor;
      weightB1K0 *=0.;
      weightR1K1 *=factor; //Just for the Dipole
      if(theta1)headR1*=0.;
      // if(weightB1K1!=0.)cout<<"\nCalcBorn1==Born1"<<weightB1K1;
  }else{
    

    
    
    if (!inhist){
      weightB1K0 =0.;
      weightB1K1 =0.;
      weightB1K0 =0.;
      weightB2K1 =0.;
    
      weightR1K1 =0.;
      weightR1K0 =0.;
      weightR2K2 =0.;
      weightR2K1 =0.;
      weightR2K0 =0.;
    }else{
      double factor=1.*CalcBorn1->children().size()*                                 //Since we choose CalcBorn1
                       CalcBorn2->children().size()/CalcBorn1->numberOfSplittings()*
                       CalcBorn1->dipol()->jacobianMerging(CalcBorn1->xcomb()->lastSHat(),
                                          CalcBorn2->xcomb()->lastSHat(),
                                          CalcBorn1->xcomb()->meMomenta().size());
      
      weightB1K1 *=factor;
      weightB1K0 *=factor;
      
      weightR1K1 *=factor;
      weightR1K0 *=factor;
    }
  }
  
  if (CalcBorn2==Born2) {
      weightB2K2 *=1.;
      weightB2K1 *=0.;
      if(theta2)headR2*=0.;
  }else{
    if (!inhist2){
      weightB2K2 =0.;
      weightB2K1 =0.;
      weightR2K2 =0.;
      weightR2K1 =0.;
      weightR2K0 =0.;
    }else{
      weightB2K2 *=1.*CalcBorn2->children().size();
      weightB2K1 *=1.*CalcBorn2->children().size();
      weightR2K2 *=1.*CalcBorn2->children().size();
      weightR2K1 *=1.*CalcBorn2->children().size();
      weightR2K0 *=1.*CalcBorn2->children().size();
    }
  }

  
  
  
  
  /**
   *
   *    Cuts!
   *    Need to be done!
   *
   **/
  
  
  if (CalcBorn2->nodeME()->projectorStage()==1&&!Born1->xcomb()->willPassCuts()){
    weightB0K0 =0.;
    weightB1K1 =0.;
    weightB1K0 =0.;
    weightB2K1 =0.;
    weightB2K2 =0.;

    weightR1K1 =0.;
    weightR1K0 =0.;
    weightR2K2 =0.;
    weightR2K1 =0.;
    weightR2K0 =0.;
      //assert(false);
  }
  if (CalcBorn2->nodeME()->projectorStage()==0&&!Born0->xcomb()->willPassCuts()){
    weightB0K0 =0.;
    weightB1K1 =0.;
    weightB1K0 =0.;
    weightB2K2 =0.;
    weightB2K1 =0.;
    weightR1K1 =0.;
    weightR1K0 =0.;
    weightR2K2 =0.;
    weightR2K1 =0.;
    weightR2K0 =0.;
      //assert(false);
  }
  if (CalcBorn2->nodeME()->projectorStage()==2&&!Born0->xcomb()->willPassCuts()){
    weightB0K0 =0.;
    weightB1K1 =0.;
    weightB1K0 =0.;
    weightB2K2 =0.;
    weightB2K1 =0.;
    weightR1K1 =0.;
    weightR1K0 =0.;
    weightR2K2 =0.;
    weightR2K1 =0.;
    weightR2K0 =0.;
      //assert(false);
  }

if(weightB0K0 == 0.&&
   weightB1K1 == 0.&&
   weightB1K0 == 0.&&
   weightB2K2 == 0.&&
   weightB2K1 == 0.&&
   weightR1K1 == 0.&&
   weightR1K0 == 0.&&
   weightR2K2 == 0.&&
   weightR2K1 == 0.&&
   weightR2K0 == 0.)return 0.;

  
  Energy startscaleB00=CKKW_StartScale(Born0);
  Energy startscaleB10=startscaleB00;
  Energy startscaleB20=startscaleB00;
  Energy startscaleB01=CKKW_StartScale(Born1);
  Energy startscaleB11=startscaleB01;
  Energy startscaleB21=startscaleB01;
  Energy startscaleB02=CKKW_StartScale(Born2);
  Energy startscaleB12=startscaleB02;
  Energy startscaleB22=startscaleB02;
  Energy startscaleR11=startscaleB00;
  Energy startscaleR21=startscaleB00;
  Energy startscaleR02=startscaleB01;
  Energy startscaleR12=startscaleB01;
  Energy startscaleR22=startscaleB01;

  
  double unlopsweightNLO1=0.;
  double unlopsweightNLO2=0.;
  
  
  
  
  Energy running;
  Energy prerunning;
  if(CalcBorn2->nodeME()->projectorStage()==0){
    
    /**
     * Relevant history weights:
     * weightB2K2
     * weightR2K2
     **/
    if(weightB2K2!=0.){
    running=startscaleB02;
    fillHistory( running,  Born2, CalcBorn2,fast);
    weightB2K2*=history.back().weight*alphaReweight()*pdfReweight();
    prerunning=CalcBorn1->dipol()->lastPt();
    if (!fillProjector(running)){cout<<"\n"<<"could not find4";return 0.;}//Actualy no projector
    prerunning=running;
    CalcBorn2->runningPt(prerunning);
    }
    if (!theta2&&theta1&&weightR2K2!=0.) {                    //subcorrections
      running=startscaleR02;
      assert(CalcBorn1->dipol()->clustersafe());
      fillHistory( running, Born1 , CalcBorn1,fast);
      prerunning=running;
      CalcBorn2->runningPt(prerunning);
      weightR2K2 *=history.back().weight*alphaReweight()*pdfReweight();
    }else{
      weightR2K2=0.;
    }
    
  }else if(CalcBorn2->nodeME()->projectorStage()==1){
    
    /**
     * Relevant history weights:
     * weightB2K1
     * weightB1K1
     * weightR1K1
     * weightR2K1
     **/
    running=startscaleB12;
    if(weightB2K1!=0.){
      fillHistory( running,  Born2, CalcBorn2,fast);
      weightB2K1*=history.back().weight*alphaReweight()*pdfReweight();
    }
    if(weightB1K1!=0.){
    prerunning=running;
    running=startscaleB11;
    fillHistory( running,  Born1, CalcBorn1,fast);
    prerunning=running;
    weightB1K1*=history.back().weight*alphaReweight()*pdfReweight();
    }
      //cout<<"\n.-.-.-.-.-1 "<<weightB1K1<<flush;
    if (theta1&&weightB1K1!=0.&&!fast)
      unlopsweightNLO2-=sumpdfReweightUnlops()+sumalphaReweightUnlops()+sumfillHistoryUnlops();
    
    if (!fillProjector(running)){cout<<"\n"<<"could not find3";return 0.;}
    prerunning=running;
    CalcBorn2->runningPt(prerunning);
    if (theta2&&theta1) {
      if(weightR2K1!=0.)weightR2K1*=history.back().weight*alphaReweight()*pdfReweight();;
    }
    
    if (!theta1&&weightR1K1!=0.) { //subcorrections
      running=startscaleR11;
      fillHistory( running, Born0 , CalcBorn0,fast);
      weightR1K1*=history.back().weight*alphaReweight()*pdfReweight();
    }
    
  }else{
    assert(CalcBorn1&&CalcBorn2->nodeME()->projectorStage()==2);
      /**
       * Relevant history weights:
       * weightB1K0
       * weightB0K0
       * weightR2K0
       * weightR1K0
       **/
    if(weightB1K0!=0.){
      running=startscaleB21;
      fillHistory( running,  Born1, CalcBorn1,fast);
      
      weightB1K0*=history.back().weight*alphaReweight()*pdfReweight();
        //cout<<"\n.-.-.-.-.-2"<<flush;
      if (theta1&&inhist&&!fast)
        unlopsweightNLO2-=sumpdfReweightUnlops()+sumalphaReweightUnlops()+sumfillHistoryUnlops();
    }else{
      assert(weightB1K0==0.);
    }

    if (weightB0K0!=0.) {
      running=startscaleB20;
      fillHistory( running,  Born0, CalcBorn0,fast);
      prerunning=running;
      weightB0K0*=history.back().weight*alphaReweight()*pdfReweight();
      //cout<<"\n.-.-.-.-.-3 "<<weightB0K0<<" "<<startscaleB20<<flush;
      if (fast)
      unlopsweightNLO1-=sumpdfReweightUnlops()+sumalphaReweightUnlops()+sumfillHistoryUnlops();
      if (!fillProjector(running)){cout<<"\n"<<"could not find1";return 0.;}
      prerunning=running;
      CalcBorn2->runningPt(prerunning);
    }
    weightR1K0*=history.back().weight*alphaReweight()*pdfReweight();
    if(weightR2K0!=0.){
      assert(theta1);
      running=startscaleR22;
      fillHistory( running,  Born1, CalcBorn1,fast);
      weightR2K0*=history.back().weight*alphaReweight()*pdfReweight();
    }
  }
  
  
  
  if(weightB0K0 == 0.&&
     weightB1K0 == 0.&&
     weightB1K1 == 0.&&
     weightB2K2 == 0.&&
     weightB2K1 == 0.&&
     weightR1K1 == 0.&&
     weightR1K0 == 0.&&
     weightR2K2 == 0.&&
     weightR2K1 == 0.&&
     weightR2K0 == 0.)return 0.;
  
  
  if(CalcBorn2->N()==CalcBorn2->nodeME()->lastMEMomenta().size()&&CalcBorn2->nodeME()->projectorStage() == 0){
      // cout<<"\n"<< CalcBorn2->runningPt()/GeV<<" "<< prerunning/GeV<<" "<<CalcBorn2->nodeME()->projectorStage();
    CalcBorn2->vetoPt(prerunning);
      //cout<<"\ncalc1 dip "<<CalcBorn1->dipol()->lastPt()/GeV<<" "<<CalcBorn0->dipol()->lastPt()/GeV;
  }else{
    CalcBorn2->vetoPt(CalcBorn2->mergePt());
      //cout<<"\n"<< CalcBorn2->runningPt()/GeV<<" "<< prerunning/GeV<<" "<<CalcBorn2->nodeME()->projectorStage();
  }
  
  
  double resLO0=0.;
  double resLO1=0.;
  double resLO2=0.;
  double resNLO1R=0.;
  double resNLO2R=0.;
  double resNLO1L=0.;
  double resNLO2L=0.;
  CNPtr selectedNode;
  Energy minpt;
  int nDipoles;
  
  

  if(CalcBorn2->nodeME()->projectorStage()==0){//-----------------------------------------------------------------------------------------
    /**   
     * Born with two extra emissions,
     * weighted with the full history.
     * weiht is always 0 if random choice of history 
     * is not the one from the getHistory call.
     * If right history is choosen reweight with number of posibilities.
     **/
    
    if( weightB2K2!=0.)      resLO2+=1.*weightB2K2* matrixElementWeight(startscaleB02,CalcBorn2);
    
    if( weightR2K2 !=0.) {
      /**
       * Real emission weighted with history up to the born state.
       * Multiply with extra Aslpha_s factor
       * Only calculate if Emission is below the merging scale.
       **/
      
      double real=matrixElementWeight(startscaleR02,CalcBorn2)*headR2;
      
      double subtraction=0.;

      if (safetheta2)  CalcBorn2->psBelowMergeingScale(selectedNode, subtraction, minpt, nDipoles);
      else             CalcBorn2->dipBelowMergeingScale(selectedNode, subtraction, minpt, nDipoles);
      
      if (CalcBorn2->xcomb()->mePartonData()[0]->coloured())
        subtraction*=CalcBorn2->nodeME()->pdf1(sqr(startscaleR02*xiFacME))/CalcBorn2->nodeME()->pdf1(sqr(10.*GeV));
      if (CalcBorn2->xcomb()->mePartonData()[1]->coloured())
        subtraction *=CalcBorn2->nodeME()->pdf2(sqr(startscaleR02*xiFacME))/CalcBorn2->nodeME()->pdf2(sqr(10.*GeV));
      
      
      resNLO2R+=(real-subtraction)*weightR2K2*theDipoleShowerHandler->as(startscaleR02*xiRenSh) / SM().alphaS();
    }
    
    
    
  }else
  if(CalcBorn2->nodeME()->projectorStage()==1){//-----------------------------------------------------------------------------------------
   
    if (weightB1K1!=0.)
      resLO1+= weightB1K1* matrixElementWeight(startscaleB11,CalcBorn1);
    
    if (weightB2K1!=0.)
      resLO2+= weightB2K1* matrixElementWeight(startscaleB12,CalcBorn2);
    
        //NLO 2
    if (weightB1K1!=0.) {
    double NLOME=matrixElementWeightWithLoops(startscaleB11,CalcBorn1,fast);
    double Bornweight=CalcBorn1->nodeME()->lastBorndSigHatDR()* SM().alphaS()/(2.*ThePEG::Constants::pi);
    
    //cout<<"\n "<<matrixElementWeight(startscaleB11,CalcBorn1)<<" "<<NLOME<<" "<<CalcBorn1->nodeME()->lastBorndSigHatDR()<<" unlops: "<<unlopsweightNLO2;
    resNLO2L+=1.*
            weightB1K1*
            (NLOME+unlopsweightNLO2*Bornweight)*
            theDipoleShowerHandler->as(startscaleB11*xiRenSh) / SM().alphaS();
    }
    if (weightR2K1!=0.) {
      assert(theta1&&theta2);
      double real= matrixElementWeight(startscaleR12,CalcBorn2)*headR2;
      
      double subtraction=1.*CalcBorn1->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
      
      if (CalcBorn2->xcomb()->mePartonData()[0]->coloured()){
        subtraction*=CalcBorn2->nodeME()->pdf1(sqr(startscaleR12*xiFacME))/CalcBorn2->nodeME()->pdf1(sqr(10.*GeV));
      }
      if (CalcBorn2->xcomb()->mePartonData()[1]->coloured()){
        subtraction *=CalcBorn2->nodeME()->pdf2(sqr(startscaleR12*xiFacME))/CalcBorn2->nodeME()->pdf2(sqr(10.*GeV));
      }
      resNLO2R+=(real-subtraction)*
                weightR2K1*
                theDipoleShowerHandler->as(startscaleR12*xiRenSh)/SM().alphaS();;
    }else{
        //TODO
        //diffPsDipBelowMergeingScale()
    }
  
    if (!theta1&&weightR1K1!=0.) {
      double real=matrixElementWeight(startscaleR11,CalcBorn1);
      double subtraction=0.;
      if (safetheta1)CalcBorn1->psBelowMergeingScale(selectedNode, subtraction, minpt, nDipoles);
      else           CalcBorn1->dipBelowMergeingScale(selectedNode, subtraction, minpt, nDipoles);
      if (CalcBorn2->xcomb()->mePartonData()[0]->coloured()){
        subtraction *=CalcBorn1->nodeME()->pdf1(sqr(startscaleR11*xiFacME))/CalcBorn1->nodeME()->pdf1(sqr(10.*GeV));
      }
      if (CalcBorn2->xcomb()->mePartonData()[1]->coloured()){
        subtraction *=CalcBorn1->nodeME()->pdf2(sqr(startscaleR11*xiFacME))/CalcBorn1->nodeME()->pdf2(sqr(10.*GeV));
      }

      resNLO1R+=(real-subtraction)*
                weightR1K1*
                theDipoleShowerHandler->as(startscaleR11*xiRenSh) / SM().alphaS();
    }
  }else
    if(CalcBorn2->nodeME()->projectorStage()==2){  //-----------------------------------------------------------------------------------------
      if(weightB1K0!=0.) resLO1+=1.* weightB1K0*matrixElementWeight(startscaleB21,CalcBorn1);
      
    if(weightB1K0!=0.){
      double NLOME=matrixElementWeightWithLoops(startscaleB11,CalcBorn1,fast);
      double Bornweight=CalcBorn1->nodeME()->lastBorndSigHatDR()* SM().alphaS()/(2.*ThePEG::Constants::pi);;
    resNLO2L+=1.*weightB1K0*
                (NLOME+unlopsweightNLO2*Bornweight)*
                theDipoleShowerHandler->as(startscaleB11*xiRenSh) / SM().alphaS();
    }
      
      

    if (weightB0K0!=0.)  resLO0+= weightB0K0*matrixElementWeight(startscaleB20,CalcBorn0);

    if (weightB0K0!=0.){
    
      double NLOME= matrixElementWeightWithLoops(startscaleB00,CalcBorn0,fast);
      double Bornweight=CalcBorn0->nodeME()->lastBorndSigHatDR()* SM().alphaS()/(2.*ThePEG::Constants::pi);;
    
      resNLO1L+= weightB0K0*
                 (NLOME+unlopsweightNLO1*Bornweight)*
                 theDipoleShowerHandler->as(startscaleB20*xiRenSh) / SM().alphaS();
    }
      
      
    
    if (weightR1K0!=0.) {
      double real=matrixElementWeight(startscaleR21,CalcBorn1);
      
      double subtraction=CalcBorn0->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
      
      if (CalcBorn2->xcomb()->mePartonData()[0]->coloured()){
        subtraction*=CalcBorn1->nodeME()->pdf1(sqr(startscaleR21*xiFacME))/
                    CalcBorn1->nodeME()->pdf1(sqr(10.*GeV));
      }
      if (CalcBorn2->xcomb()->mePartonData()[1]->coloured()){
        subtraction *=CalcBorn1->nodeME()->pdf2(sqr(startscaleR21*xiFacME))/
                     CalcBorn1->nodeME()->pdf2(sqr(10.*GeV));
      }
      resNLO1R+=(real-subtraction)*
                weightR1K0*
                theDipoleShowerHandler->as(startscaleR21*xiRenSh) / SM().alphaS();
    }
    
    if (theta2&&weightR2K0!=0.) {
      assert(theta2);
      double real=matrixElementWeight(startscaleR22,CalcBorn2);
      double subtration=CalcBorn1->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
      
      if (CalcBorn2->xcomb()->mePartonData()[0]->coloured()){
        subtration*=CalcBorn2->nodeME()->pdf1(sqr(startscaleR22*xiFacME))/
                    CalcBorn2->nodeME()->pdf1(sqr(10.*GeV));
      }
      if (CalcBorn2->xcomb()->mePartonData()[1]->coloured()){
        subtration*=CalcBorn2->nodeME()->pdf2(sqr(startscaleR22*xiFacME))/
                    CalcBorn2->nodeME()->pdf2(sqr(10.*GeV));
      }
      resNLO2R+=(real-subtration)*
               weightR2K0*
               theDipoleShowerHandler->as(startscaleR22*xiRenSh) / SM().alphaS();;

    }else{
      if (!theta2&&weightR2K0!=0.) {
       
        double real=matrixElementWeight(startscaleR22,CalcBorn2);
        double subtraction=0.;
        if (safetheta2)CalcBorn2->psBelowMergeingScale(selectedNode, subtraction, minpt, nDipoles);
        else           CalcBorn2->dipBelowMergeingScale(selectedNode, subtraction, minpt, nDipoles);
        if (CalcBorn2->xcomb()->mePartonData()[0]->coloured()){
          subtraction*=CalcBorn2->nodeME()->pdf1(sqr(startscaleR22*xiFacME))/
                       CalcBorn2->nodeME()->pdf1(sqr(10.*GeV));
        }
        if (CalcBorn2->xcomb()->mePartonData()[1]->coloured()){
          subtraction *=CalcBorn2->nodeME()->pdf2(sqr(startscaleR22*xiFacME))/
                        CalcBorn2->nodeME()->pdf2(sqr(10.*GeV));
        }
        
        resNLO2R+=(real-subtraction)*
                 weightR2K0*
                 theDipoleShowerHandler->as(startscaleR22*xiRenSh) / SM().alphaS();
      }
    }
  }

  
    //cout<<"\nres "<<resLO0<<" + "<<resLO1<<" + "<<resLO2<<" resNLO1(L+R) "<<resNLO1L<<" + "<<resNLO1R<<" resNLO2 "<<resNLO2L<<" + "<<resNLO2R<<" "<<CalcBorn2->nodeME()->projectorStage();
  return resLO0+resLO1+resLO2+resNLO1L+resNLO2L+resNLO1R+resNLO2R;
}

































double Merging::reweightCKKWBorn2(CNPtr Node,bool fast){
  if(fast||!StartingBorn1){
    if( Node->xcomb()->meMomenta().size()-2 == theMaxLegsLO){
      projectorWeight0=Node->setProjectorStage();
     // assert(false);
    }
    else{
      projectorWeight0=-1.;
      projectorWeight1=1.;
      Node->nodeME()->projectorStage(1);
     // cout<<"\n"<<Node->nodeME()->name()<<" "<<Node->nodeME()->projectorStage();
    }
    
  }
  //TEST
 // if(abs(projectorWeight)<1.3)return 0.;
  if (fast) {
    assert(!StartingBorn1);
    CalcBorn=CNPtr();
  }
  
  double weight = projectorWeight0;
  double weightCB=projectorWeight1;
  CNPtr Born;
  CNPtr BornCB;
  assert(!Node->children().empty());
  bool inhist=false;
  int randomIndex =0;
  if(Node->nodeME()->projectorStage()==1&&!CalcBorn){
    randomIndex = (int)(UseRandom::rnd() *  Node->children().size());
    CalcBorn=Node->children()[randomIndex];
    vector<Ptr<ClusterNode>::ptr> temp2 = CalcBorn->children();
    for (vector<Ptr<ClusterNode>::ptr>::iterator it = temp2.begin(); it != temp2.end(); it++) {
       if((*it)->dipol()->lastPt()<Node->mergePt())weightCB= 0.;
    }
    
    
  }
  
  if (StartingBorn1) {
    assert(!fast);
    Born=StartingBorn1;
    while (Born->parent()) {
      inhist|=(Born==CalcBorn);
      Born=Born->parent();
    }
    Born=StartingBorn1;
    BornCB =StartingBorn0;
  }else{
    Born   = Node->    getLongestHistory_simple(true,xiQSh);
    if(CalcBorn)BornCB = CalcBorn->getLongestHistory_simple(false,xiQSh);
    
    if(Node->nodeME()->projectorStage()!=1){
      vector<Ptr<ClusterNode>::ptr> temp2 = Node->children();
      for (vector<Ptr<ClusterNode>::ptr>::iterator it = temp2.begin(); it != temp2.end(); it++) {
        if((*it)->dipol()->lastPt()<Node->mergePt())return 0.;
      }
    }
    if(fast){
      StartingBorn1=Born;
      if(CalcBorn)StartingBorn0=BornCB;
      while (Born->parent()) {
        inhist|=(Born==CalcBorn);
        Born=Born->parent();
      }
      Born=StartingBorn1;
    }else{
      StartingBorn1=CNPtr();
      StartingBorn0=CNPtr();
    }
  }
  
  if(!inhist&&Node->nodeME()->projectorStage()==1){
    weight=0.;
  }
  
  
  if (BornCB&&Node->nodeME()->projectorStage()==1&&!BornCB->xcomb()->willPassCuts())weightCB=0.;
  if (!Born->xcomb()->willPassCuts())weight=0.;
  
  if(weight==0.&&weightCB==0.)return 0.;
  
  
  Energy startscale=0.*GeV;
  Energy startscaleCB=0.*GeV;
  if (BornCB&&Node->nodeME()->projectorStage()==1) {
    startscale=CKKW_StartScale(Born);
    startscaleCB=CKKW_StartScale(BornCB);
  }else{
    startscale=CKKW_StartScale(Born);
  }

  
 
  Energy running=startscale;
  Energy runningCB=startscaleCB;
  Energy prerunning;
  if(inhist||Node->nodeME()->projectorStage()!=1){
    fillHistory( running,  Born, Node,fast);
    weight*=history.back().weight;
    if (weight==0.&&Node->nodeME()->projectorStage()!=1) return 0.;
    prerunning=running;
    if (!fillProjector(prerunning))return 0.;
    weight*=alphaReweight();
    weight*=pdfReweight();
    Node->runningPt(prerunning);
    if(CalcBorn){
      prerunning=runningCB;
      fillHistory( runningCB,  BornCB, CalcBorn,fast);
      weightCB*=history.back().weight;
      weightCB*=alphaReweight()*pdfReweight();
    }
  }else{
    assert(CalcBorn);
    prerunning=runningCB;
      //cout<<"\nprerun "<<prerunning/GeV<<" "<<BornCB<<" "<<CalcBorn;
    fillHistory( runningCB,  BornCB, CalcBorn,fast);
    if (!fillProjector(prerunning))return 0.;
    weightCB*=history.back().weight;
    weightCB*=alphaReweight()*pdfReweight();
  }

  
  if(weight==0.&&weightCB==0.)return 0.;
  
  
  if(Node->N()==Node->nodeME()->lastMEMomenta().size()&&Node->nodeME()->projectorStage() == 0){
    Node->vetoPt(prerunning);
  }else{
    Node->vetoPt(Node->mergePt());
  }
  
  
  double res=0; 
    //theNf=0.;
  

  
 
  
  
  if(inhist&&Node->nodeME()->projectorStage()==1){

  //  cout<<"\ninhist "<<weight<<" "<<weightCB<<" "<<Node->children().size()<<" "<<CalcBorn->numberOfSplittings();
    
    double gluemitter=1.;
    if (CalcBorn->dipol()->bornEmitter()>2&&
        CalcBorn->xcomb()->mePartonData()[CalcBorn->dipol()->bornEmitter()]->id()
        ==21) {
      gluemitter=1.;//2.;
    }
    
    
    //TEST
    //double xx1=abs(projectorWeight)>1.1?1.:0.;
    //double xx2=abs(projectorWeight)<1.1?1.:0.;
    if (CalcBorn->xcomb()->meMomenta().size()==5||true)
      res= gluemitter*Node->children().size()/CalcBorn->numberOfSplittings()*
       weightCB*matrixElementWeight(startscale,CalcBorn)*
      CalcBorn->dipol()->jacobianMerging(CalcBorn->xcomb()->lastSHat(),Node->xcomb()->lastSHat(),CalcBorn->xcomb()->meMomenta().size());
    if (Node->xcomb()->meMomenta().size()==5||true)
      res+= Node->children().size()*weight*matrixElementWeight(startscale,Node);
    
    }else if(!inhist&&Node->nodeME()->projectorStage()==1){
    
    double gluemitter=1.;
    if (CalcBorn->dipol()->bornEmitter()>2&&
        CalcBorn->xcomb()->mePartonData()[CalcBorn->dipol()->bornEmitter()]->id()
        ==21) {
      gluemitter=1.;//2.;
    }
    if (CalcBorn->xcomb()->meMomenta().size()==5||true)
    res=gluemitter*Node->children().size()/CalcBorn->numberOfSplittings()*
        weightCB*
        CalcBorn->dipol()->jacobianMerging(CalcBorn->xcomb()->lastSHat(),Node->xcomb()->lastSHat(),CalcBorn->xcomb()->meMomenta().size())*
        matrixElementWeight(startscale,CalcBorn);
  }else{
   // cout<<"\nelse "<<weight<<" "<<weightCB;
   if (Node->xcomb()->meMomenta().size()==5||true)
    res=weight*matrixElementWeight(startscale,Node);
  }
  
    return res;//*fastweight/abs(fastweight);  }
}




double Merging::reweightCKKWBorn(CNPtr Node,bool fast){
    //cout<<"\nreweightCKKWBorn"<<flush;
  if(fast||!StartingBorn1) projectorWeight0=Node->setProjectorStage();
  
  double weight = projectorWeight0;
  CNPtr Born;
  
  if (StartingBorn1) {
    assert(!fast);
    Born = StartingBorn1;
  }else{
    Born = Node->getLongestHistory_simple(true,xiQSh);
    if(fast){
      StartingBorn1=Born;
    }else{
      StartingBorn1=CNPtr();
    }
  }
  
    //cout<<"\nwillPassCuts"<<flush;
  if(!Born->xcomb()->willPassCuts())return 0.;
  
  Energy startscale=CKKW_StartScale(Born);
  Energy running=startscale;
  
  fillHistory( running,  Born, Node,fast);
  
  weight*=history.back().weight;
    //cout<<"\nweight"<<weight<<flush;
  if (weight==0.) return 0.;
  
  Energy prerunning=running;
  if (!fillProjector(prerunning))return 0.;
  

    //  history reweighting
  weight*=alphaReweight();
  weight*=pdfReweight();
  
  Node->runningPt(prerunning);
  
  if(Node->N()==Node->nodeME()->lastMEMomenta().size()&&Node->nodeME()->projectorStage() == 0){
    Node->vetoPt(prerunning);
  }else{
    Node->vetoPt(Node->mergePt());
  }
    // cout<<"\n->"<<weight<<flush;
  double res=weight*matrixElementWeight(startscale,Node);
  return res;
}


double Merging::reweightCKKWVirt(CNPtr Node){
  
  double clusterweight= Node->setProjectorStage();
  
  CNPtr Born = Node->getLongestHistory_simple(true,xiQSh);
  
  if(!Born->xcomb()->willPassCuts())return 0.;
  
  Energy startscale=CKKW_StartScale(Born);  Energy running=startscale;
  
  fillHistory( running,  Born, Node);
  
  double sudaweight=history.back().weight;
  if (sudaweight==0.) return 0.;
  
    /////////////////////// Test ///////////////////////////
  if (!Node->NLOunitarized()){
    
    if(clusterweight<0.){
      cout<<"See: ClusterNode::setProjectorStage()";
      assert(false);
      return 0.;
    }else if((Node->M()) > Node->nodeME()->lastMEMomenta().size()){
      CKKW_PrepareSudakov(Node,running);
      Energy scale= running;
      if(minusL){
        if (!dosudakovold(Node,  scale, Node->mergePt(),clusterweight)){
        }
      }else{
        if (!dosudakov(scale, Node->mergePt(),clusterweight)){
        }
      }
      cleanup(Node);
    }
  }
    /////////////////////// Test ///////////////////////////
  
  Energy prerunning=running;
  if (!fillProjector(prerunning))return 0.;
  
    // alphas reweight of history
  double alphaweight=alphaReweight();
  
    // one additional for virt
  double extraalphaweight = theDipoleShowerHandler->as(startscale*xiRenSh) / SM().alphaS();
  
    // pdf reweight of history
  double pdfweight=pdfReweight();
  
    /////////////////////////// set scales for shower /////////////////////////////
  
  Node->runningPt(prerunning);
  
  if(Node->M()==Node->nodeME()->lastMEMomenta().size()&&Node->nodeME()->projectorStage() == 0){
    double smearing=(1.+(-1.+2.*UseRandom::rnd())*Node->smear());
    Node->runningPt(prerunning*smearing);
    Node->vetoPt(prerunning*smearing);
  }else{
    Node->vetoPt(Node->mergePt());
  }
  
  double matrixElement=matrixElementWeight(startscale,Node);
  
  
  double Bornweight=Node->nodeME()->lastBorndSigHatDR();
  
  
  double unlopsweight=0.;
  
  if(Unlopsweights){
    
      //- \sum_i \alpha_s(\mu) \partial_{\alpha_s} \\Delta |_{\alpha_s=0}
    double sumpdf=sumpdfReweightUnlops();
    unlopsweight-=sumpdf;
    
    double sumpartialalphs=sumalphaReweightUnlops();
      //assert(sumpartialalphs>=0.);
    unlopsweight-=sumpartialalphs;
    
      //- \sum_i \alpha_s(\mu) \partial_{\alpha_s} \\Delta |_{\alpha_s=0}
    double sumpartialsuda=sumfillHistoryUnlops();
    
    unlopsweight-=sumpartialsuda;
    
    unlopsweight*=Bornweight* theDipoleShowerHandler->as(startscale*xiRenSh)/(2.*ThePEG::Constants::pi);
    
  }
  
  Ptr<MergeFactory>::ptr MFactory = dynamic_ptr_cast<Ptr<MergeFactory>::ptr>(Node->nodeME()->factory());
  
  assert(MFactory);
  
  if (MFactory->onlyUnlopsweights()) {
    matrixElement=0.;
  }
  
  return matrixElement*
  sudaweight*
  alphaweight*
  extraalphaweight*
  pdfweight*
  clusterweight
  
  +unlopsweight*
  sudaweight*
  alphaweight*
  pdfweight*
  clusterweight;
  
}


double Merging::reweightCKKWReal(CNPtr Node){
  
  double weight=1.;
  bool calcHead= Node->headContribution(xiQSh);
  bool calchead2=true;
  double clusterweight= Node->setProjectorStage();
  
  /*
   
   
   R-\sum_i D_i
   =  R \prod_j \theta_j +  R (1-\prod_j\theta_j)   -\sum_i D_i\theta_i-\sum_i D_i(1-\theta_i)
   
   All above  (cluster these)
   =  R \prod_j \theta_j - \sum_i D_i\theta_i
   
   
   at least one below (should not be clustered)
   +  R (1-\prod_j\theta_j)-\sum_i D_i(1-\theta_i)
   
   
   */
  
  
  Node->finiteDipoles(UseRandom::rnd()<0.5);
  weight*= 2.;
  
  
  
  
    // \prod_i \theta_i
  
  bool prodthetai=true;
  bool safeprodthetai=true;
  double numdipcalc=0.;
  vector<CNPtr> children=Node->children();
  for (vector<CNPtr>::iterator it2 = children.begin(); it2 != children.end(); it2++)
    if ((*it2)->dipol()->clustersafe()) {
      numdipcalc+=1.;
      prodthetai&=(((*it2)->dipol()->lastPt()>(*it2)->deepHead()->mergePt()));
      safeprodthetai&=(((*it2)->dipol()->lastPt()>1.*GeV));
    }
  
  
  
  CNPtr selectedNode;double sumDipoles(0.);Energy minpt=100000.*GeV;int nDipoles(0);
  
  string type="";
  
  if(prodthetai&&!Node->finiteDipoles()) return 0.;
  
  if(!prodthetai&&!Node->finiteDipoles()) {
    
    type= "\n R (1-\\prod theta_i) - sum_i PS_i (1-\\theta_i)  --> dont cluster with probability r<\\Delta^q_minpt";
    
    if (safeprodthetai) {
      type+= " PS_i = PS_i";
      if(!(Node->psBelowMergeingScale(selectedNode, sumDipoles, minpt, nDipoles))) return 0.;
    }else{
      type+= " PS_i = D_i";
      if(!(Node->dipBelowMergeingScale(selectedNode, sumDipoles, minpt, nDipoles))) return 0.;
    }
  }
  
  if (prodthetai&&Node->finiteDipoles()) {
    
    type="\n   n_D(w_i/sum wj  R- D_i) u(\\tilde \\phi_i) -->cluster  ";
    if(!(Node->DipolesAboveMergeingScale(selectedNode, sumDipoles, minpt, nDipoles))) return 0.;
    
    
    CNPtr BornCalcHead = Node->getLongestHistory_simple(true,xiQSh);
    
    
    while(BornCalcHead->parent()&&BornCalcHead->parent()->parent()){
      BornCalcHead=BornCalcHead->parent();
    }
    if(selectedNode!=BornCalcHead){
      calcHead=false;
    }
    
  }
  
  if (!prodthetai&&Node->finiteDipoles()) {
    type="\n     nD ( PS_i -D_i) u(\\tilde \\phi_i)--> cluster" ;
    if (!safeprodthetai) {
      return 0.;
    }
    if(!(Node->diffPsDipBelowMergeingScale(selectedNode, sumDipoles, minpt, nDipoles))) return 0.;
    
  }
  
  assert(nDipoles>0);
  
  
  CNPtr Born = selectedNode->getLongestHistory_simple(false,xiQSh);
  
  if(!Born->xcomb()->willPassCuts())return 0.;
  
  Energy startscale=CKKW_StartScale(Born);  Energy running=startscale;
  
  fillHistory( running,  Born, selectedNode);
  
  weight*=history.back().weight;
  if (weight==0.) return 0.;
  
  Energy prerunning=running;
  
  
    /////////////////////// Test ///////////////////////////
  if (!Node->NLOunitarized()){
    if(clusterweight<0.){
      cout<<"See: ClusterNode::setProjectorStage()";
      assert(false);
      return 0.;
    }else if((selectedNode->deepHead()->M()) > selectedNode->nodeME()->lastMEMomenta().size()){
        // cout<<"\nreal"<<selectedNode->nodeME()->lastMEMomenta().size();
      CKKW_PrepareSudakov(selectedNode,running);
      Energy scale= running;
      if(minusL){
        if (!dosudakovold(selectedNode,  scale, selectedNode->mergePt(),clusterweight)){
        }
      }else{
        if (!dosudakov(scale, selectedNode->mergePt(),clusterweight)){
        }
      }
      cleanup(selectedNode);
    }
  }
    /////////////////////// Test ///////////////////////////
  
  
  
  bool docluster=true;
  if(!prodthetai&&!Node->finiteDipoles()&&selectedNode->inShowerPS(running)&&false) {
    
      //AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    CKKW_PrepareSudakov(selectedNode,running);
    
    double probabilityNotToCluster=1.;
    Energy tmpscale=running;
    
    if(minusL&&false){
      dosudakovold(selectedNode,  tmpscale, minpt,probabilityNotToCluster);
    }else{
      dosudakov(tmpscale, minpt,probabilityNotToCluster);
      
    }
    cleanup(selectedNode);
  }
  
  if(safeprodthetai&&!prodthetai&&!Node->finiteDipoles()&&Node->nodeME()->projectorStage()==1){
    Node->nodeME()->projectorStage(Node->nodeME()->projectorStage()-1);
    docluster=false;
  }
  
  
  if (!fillProjector(prerunning)){
    return 0.;
  }
  
  
  if (!docluster) {
    prerunning=minpt;
  }
  
    // alphas reweight of history
  weight*=alphaReweight();
  
  
    // one additional for real emission
  weight *= theDipoleShowerHandler->as(startscale*xiRenSh) / SM().alphaS();
  
  
    // pdf reweight of history
  weight*=pdfReweight();
  
  
    /////////////////////////// set scales for shower /////////////////////////////
  
  Node->runningPt(prerunning);
  
  if((Node->M()+1)==Node->nodeME()->lastMEMomenta().size()
     &&
     ((Node->nodeME()->projectorStage() == 1)||
      ((Node->nodeME()->projectorStage() == 0)&&!Node->finiteDipoles()))//Not so clear, maybe all "not clustered" reals need to be full showered
     ){
    double smearing=(1.+(-1.+2.*UseRandom::rnd())*Node->smear());
    Node->runningPt(prerunning*smearing);
    Node->vetoPt(prerunning*smearing);
  }else{
    Node->vetoPt(Node->mergePt());
  }
  
  
  double Da=-1.*selectedNode->dipol()->dSigHatDR(sqr(startscale*xiFacME))/nanobarn;
  
  vector<CNPtr> tmp2=selectedNode->children();
  for (vector<CNPtr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++) assert(((*it2)->dipol()->lastPt()>(*it2)->deepHead()->mergePt()));
  
  if (Node->xcomb()->mePartonData()[0]->coloured()){
    sumDipoles*=Node->nodeME()->pdf1(sqr(startscale*xiFacME))/Node->nodeME()->pdf1(sqr(10.*GeV));
  }
  if (Node->xcomb()->mePartonData()[1]->coloured()){
    sumDipoles *=Node->nodeME()->pdf2(sqr(startscale*xiFacME))/Node->nodeME()->pdf2(sqr(10.*GeV));
  }
  
  
  
    //       R (1-\prod theta_i (=0.) ) - sum_i PS_i (1-\theta_i) dont cluster
  if(!prodthetai&&!Node->finiteDipoles())
    weight *=matrixElementWeight(startscale,Node)*(calchead2?1.:0.)-sumDipoles;
  
    //   (R-n_D D_i) u(\tilde \phi_i)
  if (prodthetai&&Node->finiteDipoles())
    weight *= nDipoles*(matrixElementWeight(startscale,Node)*(calcHead?1.:0.)-Da);
  
    //   nD ( PS_i -D_i) u(\tilde \phi_i)
  if (!prodthetai&&Node->finiteDipoles())
    weight *= -1*nDipoles*sumDipoles;
  
  
  return weight*clusterweight;
  
}


double Merging::matrixElementWeight(Energy startscale,CNPtr Node){
  double res;
    // The deephead should be calculated here.
  CNPtr DeepHead=Node;//->deepHead();
  DeepHead->renormscale(startscale);
  DeepHead->nodeME()->factory()->scaleChoice()->setXComb(DeepHead->xcomb());
  DeepHead->nodeME()->setScale();
  DeepHead->calculateInNode(false);
  res=DeepHead->nodeME()->dSigHatDR()/nanobarn;
  DeepHead->calculateInNode(true);
  DeepHead->renormscale(0.0*GeV);
  DeepHead->calculateInNode(true);
  return res;
}

double Merging::matrixElementWeightWithLoops(Energy startscale,CNPtr Node,bool fast){
  double res=0.;
  return res;
    // The deephead should be calculated here.
  CNPtr DeepHead=Node;//->deepHead();
  DeepHead->renormscale(startscale);
  DeepHead->nodeME()->factory()->scaleChoice()->setXComb(DeepHead->xcomb());
  DeepHead->nodeME()->setScale();
  DeepHead->calculateInNode(false);
  DeepHead->nodeME()->doOneLoopNoBorn();
  res=DeepHead->nodeME()->dSigHatDR(fast)/nanobarn;
  DeepHead->nodeME()->noOneLoopNoBorn();
  DeepHead->calculateInNode(true);
  DeepHead->renormscale(0.0*GeV);
  DeepHead->calculateInNode(true);
  return res;
}

bool Merging::fillProjector(Energy& prerunning){
  if(history.begin()->node->deepHead()->nodeME()->projectorStage() == 0){
    prerunning=(history.size()==1?xiQSh:1.)*prerunning;
    return true;
  }
  for (Hist::iterator it=history.begin();it!=history.end();it++){
    if (projectorStage((*it).node)&&history.begin()->node->deepHead()->nodeME()->projectorStage() != 0){
      history.begin()->node->deepHead()->xcomb()->lastProjector((*it).node->xcomb());
      prerunning=(it==history.begin()?xiQSh:1.)*(*it).scale;
      return true;
    }
  }
  return false;
}

double Merging::pdfReweight(){
  /*
   
   Example:
   
   Delta^0_3 B_3(x'''_1,x''_2,Q_R,Q_F)= SUDA*ALPHA*PDF* f'''_1(x'''_1,Q_F)*f''_2(x''_2,Q_F)*\tilde{B}_3
   
   B_0(x_1,x_2;Q_F)  -->  B_1(x'_1,x_2;q_1) -->  B_2(x'_1,x''_2;q_2)  -->  B_3(x'''_1,x''_2;q_3)
   
   Clustering will change the PDF weights:
   f_1(x_1,Q_F)      -->  f'_1(x'_1,q_1)/f_1(x_1,q_1)                 -->  f'''_1(x'''_1,q_3)/f'_1(x'_1,q_3)
   f_2(x_2,Q_F)                             -->  f''_2(x''_1,q_2)/f_2(x_1,q_2)
   
   
   PDF  =  1/f'''_1(x'''_1,Q_F)*1/f''_2(x''_2,Q_F)
   *  f'''_1(x'''_1,q_3)/f'_1(x'_1,q_3)
   *  f''_2(x''_1,q_2)/f_2(x_1,q_2)
   *  f'_1(x'_1,q_1)/f_1(x_1,q_1)
   *  f_1(x_1,Q_F)*f_2(x_2,Q_F)
   
   =
   f'''_1(x'''_1,q_3)/f'''_1(x'''_1,Q_F)      \__  last splitting scale at the leg.
   * f''_2(x''_1,q_2)/f''_2(x''_2,Q_F)          /
   
   * f'_1(x'_1,q_1)/f'_1(x'_1,q_3)   --> third
   
   * f_1(x_1,Q_F)/f_1(x_1,q_1)       --> first
   * f_2(x_2,Q_F)/f_2(x_1,q_2)       --> second
   
   */
  
  double res=1.;
  Energy beam1Scale=history[0].scale*xiFacME;
  Energy beam2Scale=beam1Scale;
    //  Ptr<StandardEventHandler>::ptr eh =
    //  dynamic_ptr_cast<Ptr<StandardEventHandler>::ptr>(theDipoleShowerHandler->eventHandler());
  

  
  Hist::iterator it=history.begin();
  for (;it!=history.end();it++){
    Hist::iterator ittmp=it;
    ittmp++;
    if((*it).node->xcomb()->mePartonData()[0]->coloured()&&(*it).node->nodeME()->lastX1()>0.99){ return 0.;}
    if((*it).node->xcomb()->mePartonData()[1]->coloured()&&(*it).node->nodeME()->lastX2()>0.99){ return 0.;}
    
    if (ittmp!=history.end()){
      
      if((*it).node->nodeME()->lastX1()<0.00001){return 0.;}
      if((*it).node->nodeME()->lastX2()<0.00001){return 0.;}
      
      
      if ((*it).node->dipol()->bornEmitter() == 0 ){
        res *= pdfratio((*it).node, beam1Scale,xiFacSh*((*it).node->dipol()->lastPt()), 1);
        beam1Scale=xiFacSh*((*it).node->dipol()->lastPt());
      }
      if ((*it).node->dipol()->bornEmitter() == 1 ){
        res *= pdfratio((*it).node,beam2Scale, xiFacSh*((*it).node->dipol()->lastPt()), 2);
        beam2Scale=xiFacSh*((*it).node->dipol()->lastPt());
      }
      if ((*it).node->dipol()->bornSpectator() == 0 &&(*it).node->dipol()->bornEmitter() >1){//
        res *= pdfratio((*it).node, beam1Scale,xiFacSh* ((*it).node->dipol()->lastPt()), 1);
        beam1Scale=xiFacSh*((*it).node->dipol()->lastPt());
      }
      if ((*it).node->dipol()->bornSpectator() == 1 &&(*it).node->dipol()->bornEmitter() >1){//
        res *= pdfratio((*it).node, beam2Scale ,xiFacSh* ((*it).node->dipol()->lastPt()), 2);
        beam2Scale=xiFacSh*((*it).node->dipol()->lastPt());
      }
    }
  }
  if (history[0].node->deepHead()->xcomb()->mePartonData()[0]->coloured()){
    if(history[0].node->deepHead()->nodeME()->pdf1(sqr(history.back().scale))< 1e-8){return 0.;}
    res*=history[0].node->deepHead()->nodeME()->pdf1(sqr(beam1Scale))/
    history[0].node->deepHead()->nodeME()->pdf1(sqr(history[0].scale*xiFacME));
  }
  if (history[0].node->deepHead()->xcomb()->mePartonData()[1]->coloured()) {
    if(history[0].node->deepHead()->nodeME()->pdf2(sqr(history.back().scale))< 1e-8){return 0.;}
    res*=history[0].node->deepHead()->nodeME()->pdf2(sqr(beam2Scale))/
    history[0].node->deepHead()->nodeME()->pdf2(sqr(history[0].scale*xiFacME));
  }
  return res;
}

double Merging::alphaReweight(){
  /*
   
   Note: dsig is calculated with fixed alpha_s SM().alphaS()
   
   Example:
   
   Delta^0_3 B_3(Q_R,Q_F)= SUDA*ALPHA*PDF* as_fix^(3+n) *\tilde{B}_3
   
   ALPHA= as(Q_R)^n/as_fix^(n)
   * as(q_1)/as_fix
   * as(q_2)/as_fix
   * as(q_3)/as_fix
   */
  double res=1.;
  Energy Q_R=xiRenME*history[0].scale;
  
  res *= pow(as(Q_R) / SM().alphaS(), history[0].node->nodeME()->orderInAlphaS());
  res *= pow(history[0].node->deepHead()->xcomb()->eventHandler().SM().alphaEMPtr()->value(history[0].node->nodeME()->factory()->scaleChoice()->renormalizationScaleQED())/ SM().alphaEMMZ(), history[0].node->nodeME()->orderInAlphaEW());
  
  
  for (Hist::iterator it=history.begin();it!=history.end();it++){
    Hist::iterator ittmp=it;
    ittmp++;
    if ((*it).node->parent()&&ittmp!=history.end()){
      Energy q_i=xiRenSh* (*it).node->dipol()->lastPt();
      
      res *= as(q_i)/ SM().alphaS()
      *(theKImproved?(1.+((3.*(67./18.-1./6.*Constants::pi*Constants::pi)-5./9.*5.)*as(q_i))/2./Constants::pi):1.);
    }
  }
  return res;
}

void Merging::fillHistory(Energy& scale, CNPtr Begin, CNPtr EndNode,bool fast){
  
  history.clear();
  double sudakov0_n=1.;
  
  HistStep st;
  st.node=Begin;
  st.scale=scale;
  st.weight=sudakov0_n;
  
  history.push_back(st);
  Begin->nodeME()->factory()->scaleChoice()->setXComb(Begin->xcomb());
  scale*=xiQSh;
  if (Begin->parent()||!Begin->deepHead()->unitarized()) {
    while (Begin->parent() && (Begin != EndNode)) {
      
      CKKW_PrepareSudakov(Begin,scale);
      
      if(minusL){
        if (!dosudakovold(Begin,  scale, Begin->dipol()->lastPt(),sudakov0_n)){
          HistStep st;
          st.node=Begin->parent();
          st.scale=scale;
          st.weight=0.;
          history.push_back(st);
        }
      }else{
        if (!dosudakov(scale, Begin->dipol()->lastPt(),sudakov0_n,fast)){
          HistStep st;
          st.node=Begin->parent();
          st.scale=scale;
          st.weight=0.;
          history.push_back(st);
        }
      }
      cleanup(Begin);
      Begin = Begin->parent();
      HistStep st;
      st.node=Begin;
      st.scale=scale;
      st.weight=sudakov0_n;
      history.push_back(st);
    }
    Energy notunirunning=scale;
    if (!Begin->deepHead()->unitarized()&&Begin->deepHead()->N() > Begin->deepHead()->nodeME()->lastMEMomenta().size()) {
      CKKW_PrepareSudakov(Begin,scale);
      if (!dosudakovold(Begin,  notunirunning,Begin->deepHead()->mergePt(),sudakov0_n)){
        history.back().weight=0.;
      }else{
        history.back().weight=sudakov0_n;
      }
    }
    cleanup(history[0].node);
  }  
  if(  history.size()==1)scale/=xiQSh;
}




double Merging::sumpdfReweightUnlops(){
  double res=0.;
  Energy beam1Scale=history[0].scale*xiFacME;
  Energy beam2Scale=history[0].scale*xiFacME;
  Hist::iterator it=history.begin();
  
  
  
  
  for (;it!=history.end();it++){
    Hist::iterator ittmp=it;
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
                        theNf,
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
                        theNf,
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
                        theNf,
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
                        theNf,
                        history[0].scale);
          //pdfratio((*it).node, beam2Scale , sqrt((*it).node->dipol()->lastPt()), 2);
        beam2Scale=((*it).node->dipol()->lastPt())*xiFacSh;
      }
    }
  }
  if (history[0].node->deepHead()->xcomb()->mePartonData()[0]->coloured()){
    res +=pdfUnlops(history[0].node->deepHead()->nodeME()->lastParticles().first->dataPtr(),
                    history[0].node->deepHead()->nodeME()->lastPartons().first->dataPtr(),
                    history[0].node->deepHead()->xcomb()->partonBins().first->pdf(),
                    beam1Scale,
                    history[0].scale*xiFacME,
                    (history.back()).node->nodeME()->lastX1(),
                    theNf,
                    history[0].scale);
    
  }
  if (history[0].node->deepHead()->xcomb()->mePartonData()[1]->coloured()) {
    res +=pdfUnlops(history[0].node->deepHead()->nodeME()->lastParticles().second->dataPtr(),
                    history[0].node->deepHead()->nodeME()->lastPartons().second->dataPtr(),
                    history[0].node->deepHead()->xcomb()->partonBins().second->pdf(),
                    beam2Scale,
                    history[0].scale*xiFacME,
                    (history.back()).node->nodeME()->lastX2(),
                    theNf,
                    history[0].scale);
    
      //history[0].node->deepHead()->nodeME()->pdf2(sqr(beam2Scale))/history[0].node->deepHead()->nodeME()->pdf2(sqr(history[0].scale));
  }
  return res;
}

double Merging::sumalphaReweightUnlops(){
  double res=0.;
  
  
  if (!(history[0].node->children().empty())){
    res +=alphasUnlops(history[0].scale  ,
                       history[0].scale);
  }
  
  
    // dsig is calculated with fixed alpha_s
  for (Hist::iterator it=history.begin();it!=history.end();it++){
    Hist::iterator ittmp=it;
    ittmp++;
    if ((*it).node->parent()&&ittmp!=history.end()){//TODO???????
      res +=alphasUnlops((*it).node->dipol()->lastPt()*xiRenSh ,history[0].scale*xiRenME);
        //as((*it).node->dipol()->lastPt())/ SM().alphaS();
      //cout<<"\nsumalphaReweightUnlops: "<<(*it).node->dipol()->lastPt()/GeV;
    }
  }
  return res;
}

double Merging::sumfillHistoryUnlops(){
  
  double res=0.;
  
  history[0].node->nodeME()->factory()->scaleChoice()->setXComb(history[0].node->xcomb());
  
  for (Hist::iterator it = history.begin(); it != history.end();it++){
    if ((*it).node!=history.back().node) {
      CKKW_PrepareSudakov((*it).node, (it == history.begin()?xiQSh:1.)*(*it).scale);
      doUNLOPS((it == history.begin()?xiQSh:1.)*(*it).scale , (*it).node->dipol()->lastPt() , history[0].scale, res);
      cleanup((*it).node);
    }
    
  }
  
  return res;
}
















bool Merging::reweightCKKWSingle(Ptr<MatchboxXComb>::ptr SX, double & res,bool fast) {
  
  
  Ptr<StandardEventHandler>::ptr eH =
  dynamic_ptr_cast<Ptr<StandardEventHandler>::ptr>(generator()->eventHandler());
  
  
    //cout<<"\nfast ";//<<fast<<" HS "<<history.size()<<" stnode "<<StartingBorn<<" "<<eH->didEstimate();
  
  if (!eH->didEstimate()||fast) {
    history.clear();
    StartingBorn0=CNPtr();
    StartingBorn1=CNPtr();
    StartingBorn2=CNPtr();
  }
  
  theNf=5;//TODO
  
  if (!SX) return true;
  assert(SX->eventHandlerPtr());
  Ptr<MatchboxMEBase>::ptr ME = dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(SX->matchboxME());
  if (!ME) return true;
  CNPtr Node = dynamic_ptr_cast<CNPtr>(ME->firstNode());
  if (!Node) return true;
  CNPtr MENode = Node;
  Ptr<AlphaEMBase>::transient_pointer alphaEM = SX->eventHandler().SM().alphaEMPtr();
  
  assert(theDipoleShowerHandler->hardScaleIsMuF());
  
  xiRenME=ME->renormalizationScaleFactor();
  xiFacME=ME->factorizationScaleFactor();
  xiRenSh=theDipoleShowerHandler->renormalizationScaleFactor();
  xiFacSh=theDipoleShowerHandler->factorizationScaleFactor();
  xiQSh=theDipoleShowerHandler->hardScaleFactor();
  
  if(Node->deepHead()->subtractedReal()){
    res*=reweightCKKWReal(Node);
  }else if(ME->oneLoopNoBorn()){
    res*=reweightCKKWVirt(Node);
  }else{
    res*=reweightCKKWBorn3(Node,fast);
  }
  
  
  
  SX->lastCentralScale(sqr(Node->runningPt()));
  if(SX->lastProjector())
    SX->lastProjector()->lastCentralScale(sqr(Node->runningPt()));
  
  
  Node->renormscale(0.0*GeV);
  if (res == 0.){
    history.clear();
    StartingBorn0=CNPtr();
    StartingBorn1=CNPtr();
    StartingBorn2=CNPtr();
    return false;
  }
  
  
  cleanup(MENode);
  cleanup(Node);
  theDipoleShowerHandler->eventRecord().clear();
  SX->subProcess(SubProPtr());
  Node->VetoedShower(true);
  
  if (!fast) {
    history.clear();
    StartingBorn0=CNPtr();
    StartingBorn1=CNPtr();
    StartingBorn2=CNPtr();
  }
  
  return true;
  
}









double Merging::alphasUnlops( Energy next,Energy fixedScale)  {
  double betaZero =  (11./6.)*SM().Nc() - (1./3.)*theNf;
  return (betaZero*log(sqr(fixedScale/next)))+(theKImproved?(((3.*(67./18.-1./6.*Constants::pi*Constants::pi)-5./9.*theNf))):0.);
}







  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Merging::persistentOutput(PersistentOStream & os) const {
  os <<xiRenME <<xiFacME <<xiRenSh 
     <<xiFacSh <<xiQSh <<theNf 
     <<minusL <<Unlopsweights 
     << theKImproved << ounit(MergingScale,GeV)
     <<theMaxLegsLO <<theMaxLegsNLO
     <<theDipoleShowerHandler ;
}

void Merging::persistentInput(PersistentIStream & is, int) {
    is >>xiRenME >>xiFacME >>xiRenSh 
     >>xiFacSh >>xiQSh >>theNf 
     >>minusL >>Unlopsweights 
     >> theKImproved >> iunit(MergingScale,GeV)
     >>theMaxLegsLO >>theMaxLegsNLO
     >>theDipoleShowerHandler ;
}

ClassDescription<Merging> Merging::initMerging;
  // Definition of the static class description member.

void Merging::Init() {
  
  static ClassDocumentation<Merging> documentation
  ("Merging generates intrinsic pt for massless "
   "incoming partons in a shower independent way.");
  
  
  
  static Reference<Merging,DipoleShowerHandler> interfaceShowerHandler
  ("DipoleShowerHandler",
   "",
   &Merging::theDipoleShowerHandler, false, false, true, true, false);
  
  
  
  static Switch<Merging,bool>
  interfaceminusL ("minusL","",&Merging::minusL, false, false, false);
  static SwitchOption interfaceminusLYes
  (interfaceminusL,"Yes","",true);
  static SwitchOption interfaceminusLNo
  (interfaceminusL,"No","",false);
  
  static Switch<Merging,bool>
  interfaceUnlopsweights ("Unlopsweights","",&Merging::Unlopsweights, false, false, false);
  static SwitchOption interfaceUnlopsweightsYes
  (interfaceUnlopsweights,"Yes","",true);
  static SwitchOption interfaceUnlopsweightsNo
  (interfaceUnlopsweights,"No","",false);
  
  static Switch<Merging,bool>
  interfacetheKImproved ("KImproved","",&Merging::theKImproved, false, false, false);
  static SwitchOption interfacetheKImprovedYes
  (interfacetheKImproved,"Yes","",true);
  static SwitchOption interfacetheKImprovedNo
  (interfacetheKImproved,"No","",false);
  

  static Parameter<Merging,Energy> interfaceMergingScale
  ("MergingScale",
   "The MergingScale.",
   &Merging::MergingScale, GeV, 20.0*GeV, 0.0*GeV, 0*GeV,
   false, false, Interface::lowerlim);
  
  
  
  static Parameter<Merging,unsigned int> interfaceMaxLegsLO
  ("MaxLegsLO",
   ".",
   &Merging::theMaxLegsLO, 0, 0, 0,
   false, false, Interface::lowerlim);
  
  static Parameter<Merging,unsigned int> interfaceMaxLegsNLO
  ("MaxLegsNLO",
   ".",
   &Merging::theMaxLegsNLO, 0, 0, 0,
   false, false, Interface::lowerlim);
  


  
    
  
}

