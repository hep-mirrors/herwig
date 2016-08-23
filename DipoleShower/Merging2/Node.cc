  // -*- C++ -*-
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the Node class.
  //

#include "Node.h"
#include "MFactory.h"
#include "Merger.h"


#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Handlers/StdXCombGroup.h"
#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

Node::Node() {
  theclusteredto=0;
}

bool NodeDebug=false;

Node::Node(Ptr<MatchboxMEBase>::ptr nodeME, int deepprostage, int cutstage, bool nFOH) {
  thenodeMEPtr = nodeME;
  theDeepProStage = deepprostage;
  theCutStage = cutstage;
  theNeedFullOrderedHistory = nFOH;
  needsVetoedShower=false;
  theSubtractedReal=false;
  theVirtualContribution=false;
  thefiniteDipoles=false;
  thesubCorrection=false;
  theclusteredto=0;
  theOnlyN=-1;
  theNumberOfSplittings=0;
  
}



Node::Node(Ptr<Node>::ptr deephead, Ptr<Node>::ptr head, Ptr<SubtractionDipole>::ptr dipol, Ptr<MatchboxMEBase>::ptr nodeME,
                         int deepprostage, int cutstage) {
  theDeepHead = deephead;
  theparent = head;
  thedipol = dipol;
  thenodeMEPtr = nodeME;
  theDeepProStage = deepprostage;
  theCutStage = cutstage;
  theNeedFullOrderedHistory = deephead->needFullOrderedHistory();
  needsVetoedShower=false;
  theSubtractedReal=false;
  theVirtualContribution=false;
  thefiniteDipoles=false;
  thesubCorrection=false;
  theclusteredto=0;
  theOnlyN=-1;
  theNumberOfSplittings=0;
}

Node::~Node() { }



Energy Node::theMergePt(2.*GeV);
Energy Node::theCentralMergePt(2.*GeV);
double Node::smearing(0.) ;
Energy Node::theIRsafePt(1000000.*GeV);
double Node::theIRCsafeRatio(1000);
bool   Node::isUnitarized(true);
bool   Node::isNLOUnitarized(true);
int    Node::theChooseHistory(3);

unsigned int Node::theN0(0);
int  Node::theOnlyN(-1);
unsigned int Node::theN(0); 
unsigned int Node::theM(0);



Ptr<SubtractionDipole>::ptr Node::dipol() const {
  return thedipol;
}

/** returns the matrix element pointer */

const Ptr<MatchboxMEBase>::ptr Node::nodeME() const {
  return thenodeMEPtr;
}

/** access the matrix element pointer */

Ptr<MatchboxMEBase>::ptr Node::nodeME() {
  return thenodeMEPtr;
}


Ptr<Node>::ptr  Node::randomChild() {
  return thechildren[(int)(UseRandom::rnd() *  thechildren.size())];
}

bool Node::allAbove(Energy pt){
  for (vector<Ptr<Node>::ptr>::iterator it = thechildren.begin(); it != thechildren.end(); it++) {
    if((*it)->dipol()->lastPt()<pt)return false;
  }
  return true;
}



bool Node::isInHistoryOf(Ptr<Node>::ptr other){
  while (other->parent()) {
    if(other==this)return true;
    other=other->parent();
  }
  return false;
}



void Node::flushCaches() {
  this->theProjectors.clear();
  for ( unsigned int i = 0 ; i < thechildren.size() ; ++i ) {
    
    thechildren[i]->dipol()->underlyingBornME()->flushCaches();
    for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator r = thechildren[i]->dipol()->reweights().begin() ; r != thechildren[i]->dipol()->reweights().end() ;
         ++r ) {
      (**r).flushCaches();
    }
    
    if ( thechildren[i]->xcomb() ) thechildren[i]->xcomb()->clean();
    thechildren[i]->nodeME()->flushCaches();
    thechildren[i]->flushCaches();
  }
}


double Node::smear(){
  return  treefactory()->smear();
}



void Node::setKinematics() {
  for ( unsigned int i = 0 ; i < thechildren.size() ; ++i ) {
    thechildren[i]->dipol()->setXComb(thechildren[i]->xcomb());
    thechildren[i]->dipol()->setKinematics();
    thechildren[i]->nodeME()->setKinematics();
    thechildren[i]->setKinematics();
  }
}

void Node::clearKinematics() {
  for ( unsigned int i = 0 ; i < thechildren.size() ; ++i ) {
    thechildren[i]->dipol()->setXComb(thechildren[i]->xcomb());
    thechildren[i]->nodeME()->clearKinematics();
    thechildren[i]->dipol()->clearKinematics();
    thechildren[i]->clearKinematics();
  }
}

bool Node::generateKinematics(const double *r, int stage, Energy2 shat) {
  isOrdered = false;
  isthissafe=true;
  for ( unsigned int i = 0 ; i < thechildren.size() ; ++i ) {
    thechildren[i]->dipol()->setXComb(thechildren[i]->xcomb());
    if ( !thechildren[i]->dipol()->generateKinematics(r) ) cout << "stop";
    
      // 	if(!thechildren[i]->xcomb()->willPassCuts()){
      // 	  return false;
      // 	}
    if ( dipol()->lastPt() < thechildren[i]->dipol()->lastPt() && parent()->ordered() ) {
      isOrdered = true;
      deepHead()->setOrderedSteps(stage + 1);
    }
    
    thechildren[i]->generateKinematics(r, stage + 1, thechildren[i]->xcomb()->lastSHat());
    isthissafe = (isthissafe && thechildren[i]->dipol()->lastPt() >=(1+stage*deepHead()->treefactory()->stairfactor())*theDeepHead->mergePt());
  }
  return isthissafe;
}

void Node::firstgenerateKinematics(const double *r, int stage, Energy2 shat) {
  
    //TODO: Always good?
    //  if(thechildren.size()>0)
  flushCaches();
  
  
    //Set here the new merge Pt for the next phase space point.( Smearing!!!)
  
  mergePt(centralMergePt()*(1.+0.*(-1.+2.*UseRandom::rnd())*smear()));
  isOrdered = true;
  
    // if we consider the hard scale to play a role in the steps we must change here to 0.
  setOrderedSteps(1);
  this->orderedSteps();
  
  clustersafer.clear();
  for ( unsigned int i = 0 ; i < thechildren.size() ; ++i ) {
    bool ifirst = true;
    bool isecond = true;
    thechildren[i]->dipol()->setXComb(thechildren[i]->xcomb());
    
    if ( !thechildren[i]->dipol()->generateKinematics(r) ) cout << "stop";
    //cout<<"\n child shat "<<thechildren[i]->xcomb()->lastSHat()/GeV2<<" "<<shat/GeV2;
    //thechildren[i]->xcomb()->lastSHat(shat);
    isecond = thechildren[i]->generateKinematics(r, stage + 1, thechildren[i]->xcomb()->lastSHat());
    
      //TODO rethink for NLO
    ifirst = (thechildren[i]->dipol()->lastPt() >= theDeepHead->mergePt());
      // 	if(!thechildren[i]->xcomb()->willPassCuts()){
      // 	  ifirst = false;
      // 	  isecond = false;
      // 	}
    pair<pair<int, int>, int> EmitEmisSpec = make_pair(make_pair(thechildren[i]->dipol()->realEmitter(), thechildren[i]->dipol()->realEmission()),
                                                       thechildren[i]->dipol()->realSpectator());
    clustersafer.insert(make_pair(EmitEmisSpec, make_pair(ifirst, isecond)));
    
  }
  
}

void Node::setXComb(tStdXCombPtr xc, int proStage) {
  if ( !parent() ) this->xcomb(xc);
  for ( unsigned int i = 0 ; i < thechildren.size() ; ++i ) {
    if ( !thechildren[i]->xcomb() ) {
      thechildren[i]->xcomb(thechildren[i]->dipol()->makeBornXComb(xc));
      assert(thechildren[i]->xcomb());
      thechildren[i]->xcomb()->head(xc);
      if ( !thechildren[i]->dipol()->lastXCombPtr() ) {
        thechildren[i]->dipol()->setXComb(thechildren[i]->xcomb());
      }
      thechildren[i]->setXComb(thechildren[i]->xcomb(), (proStage - 1));
      
    } else {
      if ( !(thechildren[i]->dipol()->lastXCombPtr()->lastScale() == thechildren[i]->xcomb()->lastScale()) ) {
        thechildren[i]->dipol()->setXComb(thechildren[i]->xcomb());
      }
      if ( thechildren[i]->xcomb()->head() != xc ) thechildren[i]->xcomb()->head(xc);
        //Here maybe problems.
      thechildren[i]->setXComb(thechildren[i]->xcomb(), (proStage - 1));
    }
  }
}

void Node::birth(vector<Ptr<MatchboxMEBase>::ptr> vec) {

  vector<Ptr<SubtractionDipole>::ptr> dipoles = thenodeMEPtr->getDipoles(DipoleRepository::dipoles(thenodeMEPtr->factory()->dipoleSet()), vec,true);
  
  for ( unsigned int j = 0 ; j < dipoles.size() ; ++j ) {
    dipoles[j]->doSubtraction();
    Ptr<Node>::ptr node = new_ptr(
                                         Node(theDeepHead, this, dipoles[j], dipoles[j]->underlyingBornME(), theDeepHead->DeepProStage(), theDeepHead->cutStage()));
    thechildren.push_back(node);
    
  }
}


vector<Ptr<Node>::ptr> Node::getNextOrderedNodes(bool normal,double hardScaleFactor) {
    //cout<<"\ngetNextOrderedNodes";
  vector<Ptr<Node>::ptr> temp = children();
  vector<Ptr<Node>::ptr> res;
  for ( vector<Ptr<Node>::ptr>::const_iterator it = temp.begin() ; it != temp.end() ; ++it ) {
    if((*it)->deepHead()->mergePt()>(*it)->dipol()->lastPt()){
      res.clear();
        // if any of the nodes is below the merging scale return empty vector
      return res;
      continue;
    }
    if (parent()&& normal){
      if ( (*it)->dipol()->lastPt() < dipol()->lastPt() ){
        continue;
      }
    }
    if ( (*it)->children().size() != 0 ){
      
      vector<Ptr<Node>::ptr> tempdown = (*it)->children();
      for ( vector<Ptr<Node>::ptr>::const_iterator itd = tempdown.begin() ; itd != tempdown.end() ; ++itd ) {
        if( (*itd)->dipol()->lastPt() > (*it)->dipol()->lastPt()&&(*it)->inShowerPS((*itd)->dipol()->lastPt()) ){
          res.push_back(*it);
          break;
        }
      }
      
    }
    else {
      (*it)->nodeME()->factory()->scaleChoice()->setXComb((*it)->xcomb());
        //cout<<"\n"<<sqrt((*it)->nodeME()->factory()->scaleChoice()->renormalizationScale())/GeV<<" "<<(*it)->dipol()->lastPt()/GeV<<" "<<(*it)->dipol()->name();
      if ( sqr(hardScaleFactor)*(*it)->nodeME()->factory()->scaleChoice()->renormalizationScale()
          >= (*it)->dipol()->lastPt() * (*it)->dipol()->lastPt()&&
          (*it)->inShowerPS(hardScaleFactor*sqrt((*it)->nodeME()->factory()->scaleChoice()->renormalizationScale()))) {
        res.push_back(*it);
      }
    }
  }
  return res;
}

bool Node::inShowerPS(Energy hardpT){
  //return true;
  assert(deepHead()->treefactory()->largeNBasis());
  
  double x=0.;
  double z_=dipol()->lastZ();
  
    // if (dipol()->lastPt()>50*GeV) {
    //  return false;
    // }
  
  string type;
    // II
  if( dipol()->bornEmitter()<2&&dipol()->bornSpectator()<2){
    type="II";
    x =dipol()->bornEmitter()==0?xcomb()->lastX1():xcomb()->lastX2();
    double ratio = sqr(dipol()->lastPt()/dipol()->lastDipoleScale());
    double x__ = z_*(1.-z_)/(1.-z_+ratio);
    double v = ratio*z_ /(1.-z_+ratio);
    if (dipol()->lastPt()>(1.-x) * dipol()->lastDipoleScale()/ (2.*sqrt(x)))return false;
    assert(v< 1.-x__&&x > 0. && x < 1. && v > 0.);
      //return true;
    
  }
    // IF
  if( dipol()->bornEmitter()<2&&dipol()->bornSpectator()>=2){
    type="IF";
    x =dipol()->bornEmitter()==0?xcomb()->lastX1():xcomb()->lastX2();
    if (dipol()->lastPt()>dipol()->lastDipoleScale()* sqrt((1.- x)/x) /2.)return false;
      //return true;
  }
    // FI
  if( dipol()->bornEmitter()>=2&&dipol()->bornSpectator()<2){
    type="FI";
    double lastx=dipol()->bornSpectator()==0?xcomb()->lastX1():xcomb()->lastX2();
    if (dipol()->lastPt()>dipol()->lastDipoleScale()*sqrt((1.-lastx)/lastx) /2.)return false;
  }
    // FF
  if( dipol()->bornEmitter()>=2&&dipol()->bornSpectator()>=2){
    type="FF";
    if (dipol()->lastPt()>dipol()->lastDipoleScale()/2.)return false;
      // cout<<"\n "<<dipol()->bornEmitter()<<" "<<dipol()->bornSpectator()<<" "<<dipol()->lastPt()/GeV<<" "<<hardpT/GeV<<" " << sqrt(1-sqr(dipol()->lastPt()/hardpT))<<" "<<z_;
    
    
  }
  
    //return true;
  double kappa=sqr(dipol()->lastPt()/hardpT);
    //kappa=sqr(dipol()->lastPt()/dipol()->lastDipoleScale());
  double s = sqrt(1.-kappa);
  pair<double,double> zbounds= make_pair(0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s));
    //if(zbounds.first<z_&&z_<zbounds.second)cout<<"\n "<<type<<" "<<dipol()->lastPt()/GeV<<" "<<hardpT/GeV;
    //else cout<<"\n XX "<<type<<" "<<dipol()->lastPt()/GeV<<" "<<hardpT/GeV;;
  
  if (false&&dipol()->lastBornR()<dipol()->lastRealR()&&(zbounds.first<z_&&z_<zbounds.second)) {
    cout<<"\n"<<dipol()->name()<<" pt "<<dipol()->lastPt()/GeV<<" bR "<<dipol()->lastBornR()<<" rR "<<dipol()->lastRealR();
  }
  
  if (theChooseHistory==-1) {
    
    return (zbounds.first<z_&&z_<zbounds.second)&&dipol()->lastBornR()>dipol()->lastRealR();
  }
    //dipol()->lastBornR()>dipol()->lastRealR()&&
  return (zbounds.first<z_&&z_<zbounds.second);
}




Energy Node::newHardPt(){
  assert(deepHead()->treefactory()->largeNBasis());
  
  if(dipol()->underlyingBornME()->largeNColourCorrelatedME2(
                                                            make_pair(dipol()->bornEmitter(),dipol()->bornSpectator()),
                                                            deepHead()->treefactory()->largeNBasis())==0.)return 1000000000.*GeV;
  
  double x=0.;
  double z_=dipol()->lastZ();
    // II
  if( dipol()->bornEmitter()<2&&dipol()->bornSpectator()<2){
    x =xcomb()->lastX1()*xcomb()->lastX2();
  }
    // IF
  if( dipol()->bornEmitter()<2&&dipol()->bornSpectator()>=2){
    x =dipol()->bornEmitter()==0?xcomb()->lastX1():xcomb()->lastX2();
  }
  
  double chi=1-sqr(2*z_-1-x)/sqr(1-x);
  
  return dipol()->lastPt()/sqrt(abs(chi));
}




Ptr<Node>::ptr Node::getLongestHistory_simple(bool normal,double hardScaleFactor) {
  Ptr<Node>::ptr res=this;
  vector<Ptr<Node>::ptr> temp = getNextOrderedNodes(normal,hardScaleFactor);
  Energy minpt=100000.*GeV;
  Selector<Ptr<Node>::ptr> subprosel;
  while (temp.size()!=0){
    minpt=100000.*GeV;
    subprosel.clear();
    for (vector<Ptr<Node>::ptr>::iterator it=temp.begin();it!=temp.end();it++){
     
      
         if( (*it)->dipol()->underlyingBornME()->largeNColourCorrelatedME2(
                                                                       make_pair((*it)->dipol()->bornEmitter(),(*it)->dipol()->bornSpectator()),
                                                                       (*it)->deepHead()->treefactory()->largeNBasis())!=0.
         ){
        Ptr<AlphaSBase>::transient_pointer alphaS = (*it)->xcomb()->eventHandler().SM().alphaSPtr();
        
                  if((*it)->nodeME()->dSigHatDR()/nanobarn!=0.){
           subprosel.insert((abs((*it)->dipol()->dSigHatDR() /
                                (*it)->nodeME()->dSigHatDR()*alphaS->value(
                                                                           (*it)->dipol()->lastPt()*(*it)->dipol()->lastPt()))), (*it));
              minpt=min(minpt,(*it)->dipol()->lastPt());
          }
        
 
             // subprosel.insert(1., (*it));
        
        
        
        
        
        }
    }
      //    cout<<"\n"<<subprosel.size();
    if (subprosel.empty())
      return res;
    
    res = subprosel.select(UseRandom::rnd());
    
    temp = res->getNextOrderedNodes(true,hardScaleFactor);
  }
  
  if (res->parent()&&minpt>50.*GeV) {
    Energy mmmm=minpt;
    Ptr<Node>::ptr xxxx;
    vector<Ptr<Node>::ptr> tempch =children();
    
    for (vector<Ptr<Node>::ptr>::iterator it=tempch.begin();it!=tempch.end();it++){
      if (mmmm>(*it)->dipol()->lastPt()) {
        
      
      mmmm=min(mmmm,(*it)->dipol()->lastPt());
      xxxx=*it;
      }
    }
    if(false&&mmmm*3.<minpt){
        // return res->parent();
      cout<<"\n\ncluster to "<<minpt/GeV<<" dipole:  "<<res->dipol()->name()<<" "<<res->dipol()->lastPt()/GeV;
      cout<<"\n"<<res->xcomb()->meMomenta()[2].perp()/GeV<<"  mom "<<res->xcomb()->meMomenta()[2]/GeV;
      cout<<"\n"<<res->xcomb()->meMomenta()[3].perp()/GeV<<"  mom "<<res->xcomb()->meMomenta()[3]/GeV;
      
      cout<<"\n but "<<mmmm/GeV<<" with "<<xxxx->dipol()->name()<<" "<<xxxx->dipol()->lastPt()/GeV;
      cout<<"\n"<<xxxx->xcomb()->meMomenta()[2].perp()/GeV<<"  mom "<<xxxx->xcomb()->meMomenta()[2]/GeV;
      cout<<"\n"<<xxxx->xcomb()->meMomenta()[3].perp()/GeV<<"  mom "<<xxxx->xcomb()->meMomenta()[3]/GeV;
      
      cout<<"\norig";
      cout<<"\n"<<xcomb()->meMomenta()[2].perp()/GeV<<"  mom "<<xcomb()->meMomenta()[2]/GeV;
      cout<<"\n"<<xcomb()->meMomenta()[3].perp()/GeV<<"  mom "<<xcomb()->meMomenta()[3]/GeV;
      cout<<"\n"<<xcomb()->meMomenta()[4].perp()/GeV<<"  mom "<<xcomb()->meMomenta()[4]/GeV;
      
      cout<<"\n cl bornR "<<res->dipol()->lastBornR()<<" realR "<<res->dipol()->lastRealR();
      cout<<"\n sm bornR "<<xxxx->dipol()->lastBornR()<<" realR "<<xxxx->dipol()->lastRealR();
      
      
      double deta2 = sqr(xcomb()->meMomenta()[2].eta() - xcomb()->meMomenta()[3].eta());
      double dphi = abs(xcomb()->meMomenta()[2].phi() - xcomb()->meMomenta()[3].phi());
      if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
      double dr = sqrt(deta2 + sqr(dphi));
      cout<<"\nR_23 = "<<dr<<" "<<sqrt(deta2)<<" "<<dphi;
      
      deta2 = sqr(xcomb()->meMomenta()[2].eta() - xcomb()->meMomenta()[4].eta());
      dphi = abs(xcomb()->meMomenta()[2].phi() - xcomb()->meMomenta()[4].phi());
      if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
      dr = sqrt(deta2 + sqr(dphi));
      cout<<"\nR_24 = "<<dr<<" "<<sqrt(deta2)<<" "<<dphi;
      
      deta2 = sqr(xcomb()->meMomenta()[3].eta() - xcomb()->meMomenta()[4].eta());
      dphi = abs(xcomb()->meMomenta()[3].phi() - xcomb()->meMomenta()[4].phi());
      if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
      dr = sqrt(deta2 + sqr(dphi));
      cout<<"\nR_34 = "<<dr<<" "<<sqrt(deta2)<<" "<<dphi;
      
      
      
    }
  }
  
  
  return res;
}































Ptr<Node>::ptr Node::getLongestHistory(bool normal) {
  Ptr<Node>::ptr res;
  vector<Ptr<Node>::ptr> temp = getNextOrderedNodes(normal);
  vector< vector <Ptr<Node>::ptr> > Historys;
  vector< vector <Ptr<Node>::ptr> > newHistorys;
  for ( vector<Ptr<Node>::ptr>::iterator it = temp.begin() ; it != temp.end() ; it++ ) {
    vector<Ptr<Node>::ptr> x;
    x.push_back(this);
    x.push_back(*it);
    Historys.push_back(x);
    newHistorys.push_back(x);
  }
  while(newHistorys.size()!=0){
    vector< vector <Ptr<Node>::ptr> > tempnewhist=newHistorys;
    newHistorys.clear();
    for (vector< vector <Ptr<Node>::ptr> >::iterator newhist =
         tempnewhist.begin() ; newhist != tempnewhist.end() ; newhist++ ) {
      temp =newhist->back()->getNextOrderedNodes(normal);
        //       cout<<"\nts  "<<temp.size();
      for ( vector<Ptr<Node>::ptr>::iterator it = temp.begin() ; it != temp.end() ; it++ ) {
        vector<Ptr<Node>::ptr> x=*newhist;
        x.push_back(*it);
        Historys.push_back(x);
        newHistorys.push_back(x);
      }
    }
  }
  if(Historys.empty())return this;
  
  size_t longestsize=Historys.back().size();
  
  
  map<Ptr<Node>::ptr, vector<Ptr<Node>::ptr> > Tree;
    //   cout<<"\n----";
    //         parent              considered children
  for (vector< vector <Ptr<Node>::ptr> >::iterator newhist =
       Historys.end() ; newhist != Historys.begin() ;  ) {
    newhist--;
    if (newhist->size()<longestsize){
      if(!Tree.empty())break;
    }
    
    longestsize=newhist->size();
    
    
    Ptr<Node>::ptr Born=newhist->back();
    if (Born->dipol()->lastPt()<Born->deepHead()->mergePt())continue;
    Energy Scale=100000000.*GeV;
    if(!Born->children().empty()){
      vector<Ptr<Node>::ptr> temp2 = Born->children();
      for (size_t i=0;i<Born->nodeME()->mePartonData().size();i++){
        if (!Born->nodeME()->mePartonData()[i]->coloured())continue;
        for (size_t j=i+1;j<Born->nodeME()->mePartonData().size();j++){
          if (i==j||!Born->nodeME()->mePartonData()[j]->coloured())continue;
          
          
          assert(false);//This needs to be consistent with the DipoleShower Startscale if the history is not fully ordered!
          Scale=min(Scale,sqrt(2.*Born->nodeME()->lastMEMomenta()[i]*Born->nodeME()->lastMEMomenta()[j]));
        }
      }
      
    }else{
      if(NodeDebug)cout<<"\nlong hist children";
      Born->nodeME()->factory()->scaleChoice()->setXComb(Born->xcomb());
      
      Scale = sqrt(Born->nodeME()->factory()->scaleChoice()->renormalizationScale());
    }
    if (Scale==100000000.*GeV)continue;
    
    while(Born!=this){
      
        //         cout<<"\n"<<Born->dipol()->bornEmitter()<<" "<<Born->dipol()->bornSpectator()<<" "<<Born->dipol()->underlyingBornME()->largeNColourCorrelatedME2(
        // 	  make_pair(Born->dipol()->bornEmitter(),Born->dipol()->bornSpectator()),
        // 	 Born->deepHead()->treefactory()->largeNBasis());
      if(Born->dipol()->underlyingBornME()->largeNColourCorrelatedME2(
                                                                      make_pair(Born->dipol()->bornEmitter(),Born->dipol()->bornSpectator()),
                                                                      Born->deepHead()->treefactory()->largeNBasis())==0.)break;
      if(Born->dipol()->lastPt()>Scale) break;
      if (Born->dipol()->dSigHatDR() ==0.*nanobarn) break;
      if(!Born->inShowerPS(Scale))  break;
      if (Born->dipol()->lastPt()<Born->deepHead()->mergePt())break;
      if(Born->children().empty()&& !Born->xcomb()->willPassCuts() )break;
      Scale=Born->dipol()->lastPt();
      Born=Born->parent();
      
    }
    if (Born!=this)continue;
    
    
      //This history works:
    Born=newhist->back();
    
    while (Born!=this) {
      Tree[Born->parent()].push_back(Born);
      Born = Born->parent();
    }
    
    
  }
  
  res = this;
  
  
  while (Tree.find(res) != Tree.end()) {
    Ptr<Node>::ptr subpro;
    Selector<Ptr<Node>::ptr> subprosel;
      //      cout<<"\n\ntree "<< Tree[res].size();
    for ( vector<Ptr<Node>::ptr>::iterator it = Tree[res].begin() ; it != Tree[res].end() ; it++ ) {
      assert((*it)->deepHead()->mergePt()<(*it)->dipol()->lastPt());
        //       assert((*it)->xcomb()->willPassCuts() );
        // assert((*it)->dipol()->dSigHatDR() !=0.*nanobarn);
      
      
        //case0 : flat
      if((*it)->deepHead()->chooseHistory()==0){
        subprosel.insert(1., (*it));
      }
      
        //case1 : always take one of the smallest pt
      else if((*it)->deepHead()->chooseHistory()==1){
        if(!subpro) subpro=(*it);
        if((*it)->dipol()->lastPt()<subpro->dipol()->lastPt()){
          if(NodeDebug)cout<<"\n             "<<(*it)->dipol()->lastPt()/GeV<<" < "<<subpro->dipol()->lastPt()/GeV;
          subpro=(*it);
          
          subprosel.clear();
          subprosel.insert(1., (*it));
        }else if ((*it)->dipol()->lastPt()==subpro->dipol()->lastPt()  ) subprosel.insert(1., (*it));
      }
      
        //case2 :  take correspondingly to 1/pt^2
      else if((*it)->deepHead()->chooseHistory()==2){
        subprosel.insert(1./(*it)->dipol()->lastPt()/(*it)->dipol()->lastPt()*GeV2, (*it));
      }
      
        //case3 : correspondingly to the dipoleweight
      else if((*it)->deepHead()->chooseHistory()==3){
        subprosel.insert(abs((*it)->dipol()->dSigHatDR() / nanobarn), (*it));
      }
        //case3 : correspondingly to the dipoleweight
      else if((*it)->deepHead()->chooseHistory()==5){
        
        Ptr<AlphaSBase>::transient_pointer alphaS = (*it)->xcomb()->eventHandler().SM().alphaSPtr();
        subprosel.insert((abs((*it)->dipol()->dSigHatDR() / nanobarn*alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2))), (*it));
      }
      else if((*it)->deepHead()->chooseHistory()==8){
        
        Ptr<AlphaSBase>::transient_pointer alphaS = (*it)->xcomb()->eventHandler().SM().alphaSPtr();
          //           cout<<"\n"<<(*it)->dipol()->dSigHatDR()/nanobarn<<" "<<(*it)->nodeME()->dSigHatDR()/nanobarn<<abs((*it)->dipol()->dSigHatDR() / (*it)->nodeME()->dSigHatDR()*alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2));
        if((*it)->nodeME()->dSigHatDR()/nanobarn!=0.)
          subprosel.insert((abs((*it)->dipol()->dSigHatDR() / (*it)->nodeME()->dSigHatDR()*alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2))), (*it));
      }
      
      else if((*it)->deepHead()->chooseHistory()==9){
        
        Ptr<AlphaSBase>::transient_pointer alphaS = (*it)->xcomb()->eventHandler().SM().alphaSPtr();
        
        subprosel.insert((abs((*it)->dipol()->dSigHatDR() / (*it)->nodeME()->dSigHatDR()*alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2))), (*it));
          // 		 cout<<"\n"<<(abs((*it)->dipol()->dSigHatDR() / nanobarn))<<" -> "
          // 		 <<(abs((*it)->dipol()->dSigHatDR() / nanobarn*alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2)))<<" ";
        
      }
      
        //case3 : correspondingly to the dipoleweight
      else if((*it)->deepHead()->chooseHistory()==6){
        Ptr<AlphaSBase>::transient_pointer alphaS = (*it)->xcomb()->eventHandler().SM().alphaSPtr();
        subprosel.insert((alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2)), (*it));
          // 		 cout<<"\n"<<(abs((*it)->dipol()->dSigHatDR() / nanobarn))<<" -> "
          // 		 <<(abs((*it)->dipol()->dSigHatDR() / nanobarn*alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2)))<<" ";
        
      }
      
      else if((*it)->deepHead()->chooseHistory()==7){
        Ptr<Node>::ptr xx= (*it);
        Energy sum_pt=0*GeV;
        while (xx->parent()){
          sum_pt=xx->dipol()->lastPt();
          xx=xx->parent();
        }
        
        subprosel.insert(1/sum_pt*GeV, (*it));
          // 		 cout<<"\n"<<(abs((*it)->dipol()->dSigHatDR() / nanobarn))<<" -> "
          // 		 <<(abs((*it)->dipol()->dSigHatDR() / nanobarn*alphaS->value((*it)->dipol()->lastPt()*(*it)->dipol()->lastPt())/((*it)->parent()->parent()?alphaS->value((*it)->parent()->dipol()->lastPt()*(*it)->parent()->dipol()->lastPt()):0.2)))<<" ";
        
      }
      
      
        //case4 : correspondingly to the dipoleweight*alpha_s of the last splitting
      else if((*it)->deepHead()->chooseHistory()==4){
        if(!subpro) subpro=(*it);
        if(
           sqrt(2.*(*it)->parent()->nodeME()->lastMEMomenta()[(*it)->dipol()->realEmitter()]*
                (*it)->parent()->nodeME()->lastMEMomenta()[(*it)->dipol()->realEmission()])
           <
           sqrt(2.*subpro->parent()->nodeME()->lastMEMomenta()[subpro->dipol()->realEmitter()]*
                subpro->parent()->nodeME()->lastMEMomenta()[subpro->dipol()->realEmission()])
           ){
          subpro=(*it);
          subprosel.clear();
          subprosel.insert(1., (*it));
        }else if (sqrt(2.*(*it)->parent()->nodeME()->lastMEMomenta()[(*it)->dipol()->realEmitter()]*
                       (*it)->parent()->nodeME()->lastMEMomenta()[(*it)->dipol()->realEmission()])
                  ==
                  sqrt(2.*subpro->parent()->nodeME()->lastMEMomenta()[subpro->dipol()->realEmitter()]*
                       subpro->parent()->nodeME()->lastMEMomenta()[subpro->dipol()->realEmission()])  ) subprosel.insert(1., (*it));
      }
      else{
        assert(false);
      }
      
      
    }
    
    if ( !subprosel.empty() ) {
        // 	  cout<<"\nsubprosel.select "<<subprosel.size()<<flush;
      subpro = subprosel.select(UseRandom::rnd());
        // 	  cout<<"\nc"<<flush;
      subprosel.clear();
    }
      //     if(alllongestsize!=longestsize)cout<<"\n.."<<subpro->dipol()->name()<<" ";
    if ( subpro ){
      res = subpro;
        //       if(res->parent())if(res->parent()->parent())
        //       cout<<"  "<<(1.*res->clusternumber())/(1.*res->parent()->clusternumber());
      res->addclusternumber();
    }
    
    else return res;
  }
  
  return res;
}



double Node::setProjectorStage(bool fast){
  
  nodeME()->projectorStage(0);
  
  if(deepHead()->onlyN()!=-1){
    
    if(!deepHead()->subtractedReal()){
      if(nodeME()->oneLoopNoBorn()&&!NLOunitarized()){
        assert ( onlyN()==int(nodeME()->lastMEMomenta().size()));
        return 1.;
      }else{
          // if(children().empty()&&unitarized())return 1.;
        if(int(nodeME()->lastMEMomenta().size())-onlyN()==0){
          nodeME()->projectorStage(0);
          return 1.;
        }
        if(int(nodeME()->lastMEMomenta().size())-onlyN()==1){
          nodeME()->projectorStage(1);
          return -1.;
        }
        assert(false);
      }
      
      
      
    }else{
      if((N0()+1) == nodeME()->lastMEMomenta().size() || !NLOunitarized()){
        nodeME()->projectorStage(1);
        return 1.;
      }else{
        assert((N0()+1) < nodeME()->lastMEMomenta().size());
        if(int(nodeME()->lastMEMomenta().size())-onlyN()==1){
          nodeME()->projectorStage(1);
          return 1.;
        }
        if(int(nodeME()->lastMEMomenta().size())-onlyN()==2){
          nodeME()->projectorStage(2);
          return -1.;
        }
      }
      assert(false);
    }
    assert(false);
  }
  
  
  
  
  
  double res=1.;
  nodeME()->projectorStage(0);
  if(!deepHead()->subtractedReal()){
    if(nodeME()->oneLoopNoBorn()&&!NLOunitarized())return res;
    if(!children().empty()&&unitarized()){
      res*=2.;
      if (UseRandom::rnd()<0.5){
        nodeME()->projectorStage(1);
        res*=-1.;
      }else{
        nodeME()->projectorStage(0);
      }
    }
  }else{
    if(!NLOunitarized()){
      nodeME()->projectorStage(1);
      return res;
    }
    if((N0()+1) < nodeME()->lastMEMomenta().size()){
      res*=2.;
      if (UseRandom::rnd()<0.5){
        nodeME()->projectorStage(2);
        res*=-1.;
      }else{
        nodeME()->projectorStage(1);
      }
    }else{
      nodeME()->projectorStage(1);
    }
  }
  return res;
}



bool Node::headContribution(double hardScaleFactor){
  bool allabove=true;
  vector<Ptr<Node>::ptr> temp2 = children();
    //   set<int> onebelow;
  for (vector<Ptr<Node>::ptr>::iterator it = temp2.begin(); it != temp2.end(); it++) {
      //     if ((*it)->dipol()->lastPt()<mergePt()){
      //          if(onebelow.empty()){
      // 	   onebelow.insert((*it)->dipol()->realEmitter());
      // 	   onebelow.insert((*it)->dipol()->realEmission());
      // 	 }
      // 	 else{
      // 	   if(onebelow.find((*it)->dipol()->realEmitter())==onebelow.end()&&onebelow.find((*it)->dipol()->realEmission())==onebelow.end()){
      // 	    // return false;
      // 	    }
      // 	 }
      //     }
    allabove&=(*it)->dipol()->lastPt()>mergePt();
  }
  if(allabove){
    Ptr<Node>::ptr tmpBorn = getLongestHistory_simple(true,hardScaleFactor);
    if(!tmpBorn->parent())return false;
  }
  return true;
}







bool Node::DipolesAboveMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number){
  sum=0.;
  Selector<Ptr<Node>::ptr> first_subpro;
  vector<Ptr<Node>::ptr> tmp=children();
  for (vector<Ptr<Node>::ptr>::iterator it = tmp.begin(); it != tmp.end(); it++) {
    double Di=-1.* (*it)->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
    
    if ((Di!=0.)&&(*it)->xcomb()->willPassCuts()) {
      
        //vector<Ptr<Node>::ptr> tmp2=(*it)->children();
        //for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++) assert(((*it2)->dipol()->lastPt()>deepHead()->mergePt()));
      
      assert((*it)->dipol()->clustersafe());
      minpt=min(minpt,(*it)->dipol()->lastPt());
      assert((*it)->dipol()->clustersafe());
        //assert (((*it)->dipol()->lastPt()>deepHead()->mergePt()));
      sum+=Di;
      first_subpro.insert(1., (*it));
      
    }
  }
  number=int(first_subpro.size());
  if (number!=0) {
    selectedNode=first_subpro.select(UseRandom::rnd());
    return true;
  }
  return false;
}



double Node::calcPsMinusDip(Energy scale){
  return -1.* dipol()->dipMinusPs(sqr(scale),deepHead()->treefactory()->largeNBasis())/nanobarn;
}

double Node::calcPs(Energy scale){
  return dipol()->ps(sqr(scale),deepHead()->treefactory()->largeNBasis())/nanobarn;
}





bool Node::diffPsDipBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number){
  
  Selector<Ptr<Node>::ptr> first_subpro;
  vector<Ptr<Node>::ptr> tmp=children();
  for (vector<Ptr<Node>::ptr>::iterator it = tmp.begin(); it != tmp.end(); it++) {
    vector<Ptr<Node>::ptr> tmp2=(*it)->children();
    bool isInSomePS=false;
    Energy maxx=0.*GeV;
    Energy minn=100000*GeV;
    for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++){
      isInSomePS|=(*it)->inShowerPS((*it2)->dipol()->lastPt());
      maxx=max(maxx,(*it)->dipol()->lastPt());
      minn=min(minn,(*it)->dipol()->lastPt());
    }
    double Di=0.;
    if(isInSomePS||(tmp2.empty())){
      
        //cout<<"\n "<<(*it)->dipol()->lastPt()/GeV;
      
      Di=-1.* (*it)->dipol()->dipMinusPs(sqr(10.*GeV),deepHead()->treefactory()->largeNBasis())/nanobarn;
    }else{
      Di=-1.* (*it)->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
        // -1.*-1.* (*it)->dipol()->dipMinusPs(sqr(10.*GeV),deepHead()->treefactory()->largeNBasis())/nanobarn;
    }
    
    
    if ((Di!=0.)&&(*it)->xcomb()->willPassCuts()){//&&((*it)->dipol()->lastPt()<deepHead()->mergePt())) {
      vector<Ptr<Node>::ptr> tmp2=(*it)->children();
      for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++) assert(((*it2)->dipol()->lastPt()>deepHead()->mergePt()));
      
      
      
      
      if (Di!=0.) {
        first_subpro.insert(1., (*it));
        minpt=min(minpt,(*it)->dipol()->lastPt());
      }
    }
  }
  number=int(first_subpro.size());
  if (number!=0) {
    selectedNode=first_subpro.select(UseRandom::rnd());
    
    
    
    vector<Ptr<Node>::ptr> tmp2=selectedNode->children();
    bool isInSomePS=false;
    for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++){
      isInSomePS|=selectedNode->inShowerPS((*it2)->dipol()->lastPt());
    }
    
    
    if((isInSomePS||(tmp2.empty()))){//selectedNode->dipol()->lastPt()<deepHead()->mergePt()&&
      sum=-1.* selectedNode->dipol()->dipMinusPs(sqr(10.*GeV),deepHead()->treefactory()->largeNBasis())/nanobarn;
    }else{
      sum=-1.* selectedNode->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
    }
    
      //    sum=-1.*-1.* selectedNode->dipol()->dipMinusPs(sqr(10.*GeV),deepHead()->treefactory()->largeNBasis())/nanobarn;
    
    
    
    return true;
  }
  return false;
}






bool Node::psBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number){
  sum=0.;
  Selector<Ptr<Node>::ptr> first_subpro;
  vector<Ptr<Node>::ptr> tmp=children();
  for (vector<Ptr<Node>::ptr>::iterator it = tmp.begin(); it != tmp.end(); it++) {
    double Di=-1.* (*it)->dipol()->ps(sqr(10.*GeV),deepHead()->treefactory()->largeNBasis())/nanobarn;
    if ((Di!=0)&&(*it)->xcomb()->willPassCuts()){//&&((*it)->dipol()->lastPt()<deepHead()->mergePt())) {
      vector<Ptr<Node>::ptr> tmp2=(*it)->children();
      
      bool isInSomePS=false;
      
      for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++){
        assert(((*it2)->dipol()->lastPt()>deepHead()->mergePt())||((*it)->dipol()->lastPt()>deepHead()->mergePt()));
        isInSomePS|=(*it)->inShowerPS((*it2)->dipol()->lastPt());
        
      }
        //cout<<"\nisInSomePS "<< isInSomePS<< " "<< tmp2.size()<<" "<<maxx/GeV;
      
      if(!isInSomePS&&!(tmp2.empty()))continue;
      
      sum+=Di;
      if (Di!=0) {
        first_subpro.insert(1., (*it));
        minpt=min(minpt,(*it)->dipol()->lastPt());
      }
    }
  }
  number=int(first_subpro.size());
  if (number!=0) {
    selectedNode=first_subpro.select(UseRandom::rnd());
    return true;
  }
  return false;
}


bool Node::dipBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number){
  sum=0.;
  Selector<Ptr<Node>::ptr> first_subpro;
  vector<Ptr<Node>::ptr> tmp=children();
  for (vector<Ptr<Node>::ptr>::iterator it = tmp.begin(); it != tmp.end(); it++) {
    if (
	true
       //(*it)->dipol()->clustersafe() &&
       // (*it)->xcomb()->willPassCuts()
       ) {
      
      bool calcdip=true;
      vector<Ptr<Node>::ptr> tmp2=(*it)->children();
      for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++){
        assert(((*it2)->dipol()->lastPt()>deepHead()->mergePt())||((*it)->dipol()->lastPt()>deepHead()->mergePt()));
	if(((*it2)->dipol()->lastPt()<deepHead()->mergePt())){
	   if((*it)->dipol()->lastPt()<deepHead()->mergePt()){
  	     sum=0;
	     return false;
	   }
           calcdip=false;
        }
	   
	}
      
      if(calcdip){
        double Di=-1.* (*it)->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
      
        sum+=Di;
      
      minpt=min(minpt,(*it)->dipol()->lastPt());
      if (Di!=0) {
        first_subpro.insert(1., (*it));
      }
      }
    }
  }
  number=int(first_subpro.size());
  if (number!=0) {
    selectedNode=first_subpro.select(UseRandom::rnd());
    return true;
  }
  return false;
}














Ptr<MFactory>::ptr Node::treefactory(){return theTreeFactory;}

void Node::treefactory(Ptr<MFactory>::ptr x){theTreeFactory=x;}

IBPtr Node::clone() const {
  return new_ptr(*this);
}

IBPtr Node::fullclone() const {
  return new_ptr(*this);
}

  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void Node::persistentOutput(PersistentOStream & os) const {
  os << theheadxcomb <<
  thexcomb <<
  thenodeMEPtr <<
  thedipol <<
  thechildren <<
  theparent <<
  theProjectors <<
  theDeepHead <<
  theCutStage <<
  ounit(theMergePt, GeV) <<
  ounit(theCentralMergePt, GeV) <<
  smearing <<
  isthissafe <<
  ounit(theIRsafePt, GeV) <<
  theIRCsafeRatio <<
  theDeepProStage <<
  clustersafer <<
  isOrdered <<
  theOrderdSteps <<
  theNeedFullOrderedHistory <<
  theN0 <<
  theOnlyN <<
  theN <<
  theM <<
  theSudakovSteps <<
  isUnitarized <<
  isNLOUnitarized <<
  ounit(theVetoPt, GeV) <<
  ounit(theRunningPt, GeV) <<
  needsVetoedShower <<
  theCalculateInNode <<
  theNumberOfSplittings <<
  theSubtractedReal <<
  theVirtualContribution <<
  thefiniteDipoles <<
  thesubCorrection <<
  theclusteredto <<
  theTreeFactory <<
  theChooseHistory;
  
    // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void Node::persistentInput(PersistentIStream & is, int) {
  
  is >> theheadxcomb>>
  thexcomb>>
  thenodeMEPtr>>
  thedipol>>
  thechildren>>
  theparent>>
  theProjectors>>
  theDeepHead>>
  theCutStage>>
  iunit(theMergePt, GeV)>>
  iunit(theCentralMergePt, GeV)>>
  smearing>>
  isthissafe>>
  iunit(theIRsafePt, GeV)>>
  theIRCsafeRatio>>
  theDeepProStage>>
  clustersafer>>
  isOrdered>>
  theOrderdSteps>>
  theNeedFullOrderedHistory>>
  theN0>>
  theOnlyN>>
  theN>>
  theM>>
  theSudakovSteps>>
  isUnitarized>>
  isNLOUnitarized>>
  iunit(theVetoPt, GeV)>>
  iunit(theRunningPt, GeV)>>
  needsVetoedShower>>
  theCalculateInNode>>
  theNumberOfSplittings>>
  theSubtractedReal>>
  theVirtualContribution>>
  thefiniteDipoles>>
  thesubCorrection>>
  theclusteredto>>  theTreeFactory>>
  theChooseHistory;
  
    // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<Node,Interfaced> describeHerwigNode("Herwig::Node", "HwDipoleShower.so");

void Node::Init() {
  
  static ClassDocumentation<Node> documentation("There is no documentation for the Node class");
  
}

