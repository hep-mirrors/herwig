  // -*- C++ -*-
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the Node class.
  //

#include "Node.h"
#include "MergingFactory.h"
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

Node::Node(Ptr<MatchboxMEBase>::ptr nodeME, int cutstage,Ptr<Merger>::ptr mh)
  :Interfaced(),
 thenodeMEPtr(nodeME),
 thedipol(),
 thechildren(),
 theparent(),
 theDeepHead(this),
 theCutStage(cutstage),
 isOrdered(true),
 theSubtractedReal(false),
 theVirtualContribution(false),
 theMergingHelper(mh)
{
  nodeME->maxMultCKKW(1);
  nodeME->minMultCKKW(0);
}



Node::Node(Ptr<Node>::ptr deephead,
           Ptr<Node>::ptr head,
           Ptr<SubtractionDipole>::ptr dipol,
           Ptr<MatchboxMEBase>::ptr nodeME,
           int cutstage)
:Interfaced(),thenodeMEPtr(nodeME),
thedipol(dipol),
theparent(head),
theDeepHead(deephead),
theCutStage(cutstage),
isOrdered(true),
theSubtractedReal(false),
theVirtualContribution(false),
theMergingHelper() //The subnodes have no merging helper
{
}

Node::~Node() { }

Ptr<SubtractionDipole>::ptr Node::dipole() const {
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

int Node::legsize() const {return nodeME()->legsize();}

Ptr<Node>::ptr  Node::randomChild() {
  return thechildren[(int)(UseRandom::rnd() *  thechildren.size())];
}

bool Node::allAbove(Energy pt){
  for (Ptr<Node>::ptr child : thechildren)
    if(child->pT()<pt)return false;
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
  for ( auto const & ch: thechildren) {
    ch->dipole()->underlyingBornME()->flushCaches();
    for (Ptr<MatchboxReweightBase>::ptr r : ch->dipole()->reweights())
      r->flushCaches();
    if ( ch->xcomb() ) ch->xcomb()->clean();
    ch->nodeME()->flushCaches();
    ch->flushCaches();
  }
}

void Node::setKinematics() {
  for (auto const & ch: thechildren) {
    ch->dipole()->setXComb(ch->xcomb());
    ch->dipole()->setKinematics();
    ch->nodeME()->setKinematics();
    ch->setKinematics();
  }
}

void Node::clearKinematics() {
  for (auto const & ch: thechildren) {
    ch->dipole()->setXComb(ch->xcomb());
    ch->nodeME()->clearKinematics();
    ch->dipole()->clearKinematics();
    ch->clearKinematics();
  }
}

bool Node::generateKinematics(const double *r, int stage, Energy2 ) {
  bool isthissafe=true;
  for (auto const & ch: thechildren) {
    ch->dipole()->setXComb(ch->xcomb());
    if ( !ch->dipole()->generateKinematics(r) ) cout << "stop";
    ch->generateKinematics(r, stage + 1, ch->xcomb()->lastSHat());
    isthissafe = (isthissafe && ch->pT() >=deepHead()->MH()->mergePt());
  }
  return isthissafe;
}

void Node::firstgenerateKinematics(const double *r, int stage) {
  flushCaches();
  
  MH()->smeareMergePt();
    //Set here the new merge Pt for the next phase space point.( Smearing!!!)
  clustersafer.clear();
  for (auto const & ch: thechildren) {
    bool ifirst = true;
    bool isecond = true;
    ch->dipole()->setXComb(ch->xcomb());
    
    if ( !ch->dipole()->generateKinematics(r) ) cout << "stop";
    
    isecond = ch->generateKinematics(r, stage + 1, ch->xcomb()->lastSHat());
    ifirst = (ch->pT() >= deepHead()->MH()->mergePt());
    
    pair<pair<int, int>, int> EmitEmisSpec =
    make_pair(make_pair(ch->dipole()->realEmitter(),
                        ch->dipole()->realEmission()),
              ch->dipole()->realSpectator());
    clustersafer.insert(make_pair(EmitEmisSpec, make_pair(ifirst, isecond)));
    
  }
}



void Node::setXComb(tStdXCombPtr xc) {
  if ( !parent() ) this->xcomb(xc);
  for (auto const & ch: thechildren) {
    if ( !ch->xcomb() ) {
      ch->xcomb(ch->dipole()->makeBornXComb(xc));
      ch->xcomb()->head(xc);
      if ( !ch->dipole()->lastXCombPtr() ) {
        ch->dipole()->setXComb(ch->xcomb());
      }
      ch->setXComb(ch->xcomb());
      
    } else {
      if ( !(ch->dipole()->lastXCombPtr()->lastScale() == ch->xcomb()->lastScale()) ) {
        ch->dipole()->setXComb(ch->xcomb());
      }
      if ( ch->xcomb()->head() != xc ) ch->xcomb()->head(xc);
      ch->setXComb(ch->xcomb());
    }
  }
}

void Node::birth(vector<Ptr<MatchboxMEBase>::ptr> vec) {

  vector<Ptr<SubtractionDipole>::ptr> dipoles =
           nodeME()->getDipoles(DipoleRepository::dipoles(
                        nodeME()->factory()->dipoleSet()), vec,true);
  
  for ( auto const & dip : dipoles ) {
    dip->doSubtraction();
    Ptr<Node>::ptr node = new_ptr(Node(theDeepHead,
                                       this,
                                       dip,
                                       dip->underlyingBornME(),
                                       theDeepHead->cutStage()));
    thechildren.push_back(node);
    
  }
}


vector<Ptr<Node>::ptr> Node::getNextOrderedNodes(bool normal,double hardScaleFactor) {

  vector<Ptr<Node>::ptr> temp = children();
  vector<Ptr<Node>::ptr> res;

  for (Ptr<Node>::ptr  const & child : children()) {
    if(deepHead()->MH()->mergePt()>child->pT()){
      res.clear();
      return res;
    }
  }
  
  for (Ptr<Node>::ptr const & child: children()) {

    if (parent()&& normal){
      if ( child->pT() < pT() ){
        continue;
      }
    }
    if ( child->children().size() != 0 ){
      
      for (Ptr<Node>::ptr itChild:  child->children()) {
        if( itChild->pT() > child->pT()&&child->inShowerPS(itChild->pT()) ){
          res.push_back(child);
          break;
        }
      }
    }
    else {
      child->nodeME()->factory()->scaleChoice()->setXComb(child->xcomb());
      if ( sqr(hardScaleFactor)*child->nodeME()->factory()->scaleChoice()->renormalizationScale()
          >= sqr(child->pT()) &&
          child->inShowerPS(hardScaleFactor*sqrt(child->nodeME()->factory()->scaleChoice()->renormalizationScale()))) {
        res.push_back(child);
      }
    }
  }
  return res;
}

bool Node::inShowerPS(Energy hardpT){
    //Here we decide if the current phase space point can be reached from the underlying Node.
  double z_=dipole()->lastZ();
    // II
  if( dipole()->bornEmitter()<2&&dipole()->bornSpectator()<2&&deepHead()->MH()->openInitialStateZ()) return true;
    // IF
  if( dipole()->bornEmitter()<2&&dipole()->bornSpectator()>=2&&deepHead()->MH()->openInitialStateZ())
    return true;
  
  pair<double,double> zbounds=
  dipole()->tildeKinematics()->zBounds(pT(),hardpT);
  
  return (zbounds.first<z_&&z_<zbounds.second);
}






Ptr<Node>::ptr Node::getHistory(bool normal,double hardScaleFactor) {
  Ptr<Node>::ptr res=this;
    //cout<<"\nstart get next"<<flush;
  vector<Ptr<Node>::ptr> temp = getNextOrderedNodes(normal,hardScaleFactor);

  Energy minpt=100000.*GeV;
  Selector<Ptr<Node>::ptr> subprosel;
  while (temp.size()!=0){
    minpt=100000.*GeV;
    subprosel.clear();
    for (Ptr<Node>::ptr const &  child : temp) {
      if( child->dipole()->underlyingBornME()->largeNColourCorrelatedME2(
                                                                         make_pair(child->dipole()->bornEmitter(),child->dipole()->bornSpectator()),
                                                                         deepHead()->MH()->largeNBasis())!=0.
         ){
        
        double weight=abs(child->dipole()->dSigHatDR()/nanobarn);
        if(weight!=0.){
          subprosel.insert(weight , child);
          minpt=min(minpt,child->pT());
        }
        
        /*
         if((*it)->nodeME()->dSigHatDR()/nanobarn!=0.){
         subprosel.insert((abs((*it)->dipole()->dSigHatDR() /
         (*it)->nodeME()->dSigHatDR()*deepHead()->MH()->as((*it)->pT()))), (*it));
         minpt=min(minpt,(*it)->pT());
         }
         */
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










pair<CrossSection,CrossSection> Node::calcDipandPS(Energy scale){
  return  dipole()->dipandPs(sqr(scale),deepHead()->MH()->largeNBasis());
}

CrossSection Node::calcPs(Energy scale){
  return dipole()->ps(sqr(scale),deepHead()->MH()->largeNBasis());
}

CrossSection Node::calcDip(Energy scale){
  return dipole()->dip(sqr(scale));
}


IBPtr Node::clone() const {
  return new_ptr(*this);
}

IBPtr Node::fullclone() const {
  return new_ptr(*this);
}


void Node::persistentOutput(PersistentOStream & os) const {
  
  os <<
  thexcomb<<
  thenodeMEPtr<<
  thedipol<<
  thechildren<<
  theparent<<
  theProjectors<<
  theDeepHead<<
  theCutStage<<
  clustersafer<<
  ounit(theRunningPt,GeV)<<
  theSubtractedReal<<
  theVirtualContribution<<
  theMergingHelper;
}

void Node::persistentInput(PersistentIStream & is, int) {
  
  is >>
  thexcomb>>
  thenodeMEPtr>>
  thedipol>>
  thechildren>>
  theparent>>
  theProjectors>>
  theDeepHead>>
  theCutStage>>
  clustersafer>>
  iunit(theRunningPt,GeV)>>
  theSubtractedReal>>
  theVirtualContribution>>
  theMergingHelper;
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

