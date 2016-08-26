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
}

bool NodeDebug=false;

Node::Node(Ptr<MatchboxMEBase>::ptr nodeME, int deepprostage, int cutstage) {
  nodeME->maxMultCKKW(1);
  nodeME->minMultCKKW(0);
  theDeepHead = this;
  thenodeMEPtr = nodeME;
  theDeepProStage = deepprostage;
  theCutStage = cutstage;
  theSubtractedReal=false;
  theVirtualContribution=false;
  
}



Node::Node(Ptr<Node>::ptr deephead, Ptr<Node>::ptr head, Ptr<SubtractionDipole>::ptr dipol, Ptr<MatchboxMEBase>::ptr nodeME,
                         int deepprostage, int cutstage) {
  theDeepHead = deephead;
  theparent = head;
  thedipol = dipol;
  thenodeMEPtr = nodeME;
  theDeepProStage = deepprostage;
  theCutStage = cutstage;
  theSubtractedReal=false;
  theVirtualContribution=false;
}

Node::~Node() { }

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

Energy Node::miniPt() const{
   Energy res=1000000000*GeV;
   for (vector<Ptr<Node>::ptr>::const_iterator it = thechildren.begin(); it != thechildren.end(); it++) {
    res=min(res,(*it)->dipol()->lastPt());
  }
  return res;
   	
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
    if ( dipol()->lastPt() < thechildren[i]->dipol()->lastPt() && parent()->ordered() ) {
      isOrdered = true;
      deepHead()->setOrderedSteps(stage + 1);
    }
    thechildren[i]->generateKinematics(r, stage + 1, thechildren[i]->xcomb()->lastSHat());
    isthissafe = (isthissafe && thechildren[i]->dipol()->lastPt() >=deepHead()->MH()->mergePt());
  }
  return isthissafe;
}

void Node::firstgenerateKinematics(const double *r, int stage, Energy2 shat) {
  flushCaches();
  
  
  
    //Set here the new merge Pt for the next phase space point.( Smearing!!!)
  
  MH()->smeareMergePt();
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
    
    isecond = thechildren[i]->generateKinematics(r, stage + 1, thechildren[i]->xcomb()->lastSHat());
    ifirst = (thechildren[i]->dipol()->lastPt() >= deepHead()->MH()->mergePt());

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
      thechildren[i]->setXComb(thechildren[i]->xcomb(), (proStage - 1));
    }
  }
}

void Node::birth(vector<Ptr<MatchboxMEBase>::ptr> vec) {

  vector<Ptr<SubtractionDipole>::ptr> dipoles =
           nodeME()->getDipoles(DipoleRepository::dipoles(
                        nodeME()->factory()->dipoleSet()), vec,true);
  
  for ( unsigned int j = 0 ; j < dipoles.size() ; ++j ) {
    dipoles[j]->doSubtraction();
    Ptr<Node>::ptr node = new_ptr(Node(theDeepHead,this, dipoles[j],
                                       dipoles[j]->underlyingBornME(),
                                       theDeepHead->DeepProStage(),
                                       theDeepHead->cutStage()));
    thechildren.push_back(node);
    
  }
}


vector<Ptr<Node>::ptr> Node::getNextOrderedNodes(bool normal,double hardScaleFactor) {

  vector<Ptr<Node>::ptr> temp = children();
  vector<Ptr<Node>::ptr> res;
  for ( vector<Ptr<Node>::ptr>::const_iterator it = temp.begin() ; it != temp.end() ; ++it ) {
    if(deepHead()->MH()->mergePt()>(*it)->dipol()->lastPt()){
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

  assert(deepHead()->MH()->largeNBasis());
  
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
    if (deepHead()->MH()->openInitialStateZ()) {
      return true;
    }
  }
    // IF
  if( dipol()->bornEmitter()<2&&dipol()->bornSpectator()>=2){
    type="IF";
    x =dipol()->bornEmitter()==0?xcomb()->lastX1():xcomb()->lastX2();
    if (dipol()->lastPt()>dipol()->lastDipoleScale()* sqrt((1.- x)/x) /2.)return false;
    if (deepHead()->MH()->openInitialStateZ()) {
      return true;
    }
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
  }
  
  double kappa=sqr(dipol()->lastPt()/hardpT);

  double s = sqrt(1.-kappa);
  pair<double,double> zbounds= make_pair(0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s));
 
  return (zbounds.first<z_&&z_<zbounds.second);
}






Ptr<Node>::ptr Node::getHistory(bool normal,double hardScaleFactor) {
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
                                  deepHead()->MH()->largeNBasis())!=0.
         ){
         
          if((*it)->nodeME()->dSigHatDR()/nanobarn!=0.){
             subprosel.insert((abs((*it)->dipol()->dSigHatDR() /
              (*it)->nodeME()->dSigHatDR()*deepHead()->MH()->as((*it)->dipol()->lastPt()))), (*it));
              minpt=min(minpt,(*it)->dipol()->lastPt());
          }
        
           //TODO choosehistories
        
        
        }
    }
    if (subprosel.empty())
      return res;
    
    res = subprosel.select(UseRandom::rnd());
    temp = res->getNextOrderedNodes(true,hardScaleFactor);
  }
  return res;
}



bool Node::headContribution(double hardScaleFactor){
  bool allabove=true;
  vector<Ptr<Node>::ptr> temp2 = children();
  for (vector<Ptr<Node>::ptr>::iterator it = temp2.begin(); it != temp2.end(); it++) {
    allabove&=(*it)->dipol()->lastPt()>deepHead()->MH()->mergePt();
  }
  if(allabove){
    Ptr<Node>::ptr tmpBorn = getHistory(true,hardScaleFactor);
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
      assert((*it)->dipol()->clustersafe());
      minpt=min(minpt,(*it)->dipol()->lastPt());
      assert((*it)->dipol()->clustersafe());
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



pair<double,double> Node::calcDipandPS(Energy scale){
  return  dipol()->dipandPs(sqr(scale),deepHead()->MH()->largeNBasis());
}

double Node::calcPs(Energy scale){
  return dipol()->ps(sqr(scale),deepHead()->MH()->largeNBasis())/nanobarn;
}




/*
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
    
      Di=-1.* (*it)->dipol()->dipMinusPs(sqr(10.*GeV),deepHead()->MH()->largeNBasis())/nanobarn;
    }else{
      Di=-1.* (*it)->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
    }
    
    
    if ((Di!=0.)&&(*it)->xcomb()->willPassCuts()){//&&((*it)->dipol()->lastPt()<MH()->mergePt())) {
      vector<Ptr<Node>::ptr> tmp2=(*it)->children();
      for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++) assert(((*it2)->dipol()->lastPt()>deepHead()->MH()->mergePt()));
      
      
      
      
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
    
    
    if((isInSomePS||(tmp2.empty()))){//selectedNode->dipol()->lastPt()<MH()->mergePt()&&
      sum=-1.* selectedNode->dipol()->dipMinusPs(sqr(10.*GeV),deepHead()->MH()->largeNBasis())/nanobarn;
    }else{
      sum=-1.* selectedNode->dipol()->dSigHatDR(sqr(10.*GeV))/nanobarn;
    }
    return true;
  }
  return false;
}

*/

/*


bool Node::psBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number){
  sum=0.;
  Selector<Ptr<Node>::ptr> first_subpro;
  vector<Ptr<Node>::ptr> tmp=children();
  for (vector<Ptr<Node>::ptr>::iterator it = tmp.begin(); it != tmp.end(); it++) {
    double Di=-1.* (*it)->dipol()->ps(sqr(10.*GeV),deepHead()->MH()->largeNBasis())/nanobarn;
    if ((Di!=0)&&(*it)->xcomb()->willPassCuts()){//&&((*it)->dipol()->lastPt()<MH()->mergePt())) {
      vector<Ptr<Node>::ptr> tmp2=(*it)->children();
      
      bool isInSomePS=false;
      
      for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++){
        assert(((*it2)->dipol()->lastPt()>deepHead()->MH()->mergePt())||((*it)->dipol()->lastPt()>deepHead()->MH()->mergePt()));
        isInSomePS|=(*it)->inShowerPS((*it2)->dipol()->lastPt());
        
      }
      
      
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
*/

/*

bool Node::dipBelowMergeingScale(Ptr<Node>::ptr& selectedNode,double & sum,Energy& minpt,int& number){
  sum=0.;
  Selector<Ptr<Node>::ptr> first_subpro;
  vector<Ptr<Node>::ptr> tmp=children();
  for (vector<Ptr<Node>::ptr>::iterator it = tmp.begin(); it != tmp.end(); it++) {
    if (true
       //(*it)->dipol()->clustersafe() &&
       // (*it)->xcomb()->willPassCuts()
       ) {
      
      bool calcdip=true;
      vector<Ptr<Node>::ptr> tmp2=(*it)->children();
      for (vector<Ptr<Node>::ptr>::iterator it2 = tmp2.begin(); it2 != tmp2.end(); it2++){
        assert(((*it2)->dipol()->lastPt()>deepHead()->MH()->mergePt())||((*it)->dipol()->lastPt()>deepHead()->MH()->mergePt()));
	if(((*it2)->dipol()->lastPt()<deepHead()->MH()->mergePt())){
	   if((*it)->dipol()->lastPt()<deepHead()->MH()->mergePt()){
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
*/


IBPtr Node::clone() const {
  return new_ptr(*this);
}

IBPtr Node::fullclone() const {
  return new_ptr(*this);
}

  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void Node::persistentOutput(PersistentOStream & os) const {
  os <<
  theheadxcomb             <<  thexcomb          <<  thenodeMEPtr <<
  thedipol                 <<  thechildren   <<theparent <<
            theDeepHead       <<  theCutStage <<
  theDeepProStage          <<  clustersafer      <<  ounit(theVetoPt, GeV) <<
  ounit(theRunningPt, GeV) <<  theSubtractedReal <<  theVirtualContribution<<theMergingHelper;
  
  
  
    // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void Node::persistentInput(PersistentIStream & is, int) {
  
  is >>
  theheadxcomb             >>  thexcomb          >>  thenodeMEPtr >>
  thedipol                 >>  thechildren       >>  theparent >>  theDeepHead       >>  theCutStage >>
  theDeepProStage          >>  clustersafer      >>  iunit(theVetoPt, GeV) >>
  iunit(theRunningPt, GeV) >>  theSubtractedReal >>  theVirtualContribution
  >>theMergingHelper;
  
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

