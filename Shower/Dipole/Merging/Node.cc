  // -*- C++ -*-
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the Node class.
  //

#include "Node.h"
#include "MergingFactory.h"
#include "Merger.h"

using namespace Herwig;

Node::Node(MatchboxMEBasePtr nodeME, int cutstage, MergerPtr mh)
  :Interfaced(), 
 thenodeMEPtr(nodeME), 
 thedipol(), 
 theparent(),
 theCutStage(cutstage), 
 //isOrdered(true), 
 theSubtractedReal(false), 
 theVirtualContribution(false), 
 theMergingHelper(mh)
{
  nodeME->maxMultCKKW(1);
  nodeME->minMultCKKW(0);
}



Node::Node(NodePtr deephead, 
           NodePtr head, 
           SubtractionDipolePtr dipol, 
           MatchboxMEBasePtr nodeME, 
           int cutstage)
:Interfaced(), thenodeMEPtr(nodeME), 
thedipol(dipol), 
theparent(head), 
theDeepHead(deephead), 
theCutStage(cutstage), 
//isOrdered(true), 
theSubtractedReal(false), 
theVirtualContribution(false),
 theMergingHelper() //The subnodes have no merging helper
{
}

Node::~Node() { }

SubtractionDipolePtr Node::dipole() const {
  return thedipol;
}

/** returns the matrix element pointer */

const MatchboxMEBasePtr Node::nodeME() const {
  return thenodeMEPtr;
}

/** access the matrix element pointer */

MatchboxMEBasePtr Node::nodeME() {
  return thenodeMEPtr;
}

pair<PVector , PVector> Node::getInOut( ){
    PVector in;
    const auto me= nodeME();
    const auto pd=me->mePartonData();
    for( auto i : {0 , 1} )
     in.push_back(pd[i]->produceParticle( me->lastMEMomenta()[i] ) );
    PVector out;
    for ( size_t i = 2;i< pd.size();i++ ){
      PPtr p  = pd[i]->produceParticle( me->lastMEMomenta()[i] );
      out.push_back( p );
    }
    return { in , out };
}

int Node::legsize() const {return nodeME()->legsize();}

NodePtr  Node::randomChild() {
  return thechildren[UseRandom::irnd(thechildren.size())];
}

bool Node::allAbove(Energy pt) {
  for (NodePtr child : thechildren)
    if ( child->pT() < pt )
    	return false;
  return true;
}

Energy Node::maxChildPt(){
  Energy maxi=-1*GeV;
  for (NodePtr child : thechildren)maxi=max(child->pT(),maxi);
  return maxi;
}


bool Node::isInHistoryOf(NodePtr other) {
  while (other->parent()) {
    if (other == this) 
    	return true;
    other = other->parent();
  }
  return false;
}


void Node::flushCaches() {
  if (didflush) return;
  didflush=true;
  for ( auto const & ch: thechildren) {
    ch->xcomb()->clean();
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

bool Node::generateKinematics(const double *r, bool directCut) {
  didflush=false;
  // If there are no children to the child process we are done.
  if(children().empty()) return true;
  
  assert(parent());

  if ( ! directCut && pT() < deepHead()->MH()->mergePt()) {
      // Real emission:
      // If there are children to the child process,
      // we now require that all subsequent children
      // with pt < merging scale are in their ME region.
      // Since the possible children of the real emission
      // contribution are now in their ME region,
      // it is clear that a second clustering is possible
      // -- modulo phase space restrictions.
      // Therefore the real emission contribution
      // are unitarised and the cross section is
      // hardly modified.
    auto inOutPair = getInOut();
    NodePtr rc =  randomChild();
    rc->dipole()->setXComb(rc->xcomb());
    if(!rc->dipole()->generateKinematics(r))assert(false);
    
      // If not in ME -> return false
    if(!deepHead()->MH()->matrixElementRegion( inOutPair.first ,
                            inOutPair.second ,
                            rc->pT() ,
                            deepHead()->MH()->mergePt() ) )return false;
  }
  for (auto & ch : children() ) {
    ch->dipole()->setXComb(ch->xcomb());
    if ( !ch->dipole()->generateKinematics(r) ) { assert(false); }
    ch->generateKinematics( r, true);
  }
  return true;
}

bool Node::firstgenerateKinematics(const double *r, bool directCut) {
  didflush=false;
    // This is called form the merging helper for the first node. So:
  assert(!parent());
  assert(xcomb());
  
  
  if(MH()->treefactory()->nonQCDCuts()){
  tcPDVector outdata(xcomb()->mePartonData().begin()+2,
                     xcomb()->mePartonData().end());
  vector<LorentzMomentum> outmomenta(xcomb()->meMomenta().begin()+2,
                                     xcomb()->meMomenta().end());
  

  if ( !MH()->treefactory()->nonQCDCuts()->passCuts(outdata,outmomenta,
                         xcomb()->mePartonData()[0],
                         xcomb()->mePartonData()[1]) )
    return false;
  }
  
   ///// This should not be needed!!!
   ///// ( Warning inMerger::matrixElementRegion gets triggered.)
  flushCaches();
    //Set here the new merge Pt for the next phase space point.( Smearing!!!)
  MH()->smearMergePt();
    // If there are no children to this node, we are done here:
  if (children().empty())
     return true;
    // directCut is for born and for virtual contributions.
    // if directCut is true, then cut on the first ME region.
    // call recursiv generate kinematics for subsequent nodes.
  if ( directCut ){
    
    auto inOutPair = getInOut();

    NodePtr rc = randomChild();
    rc->dipole()->setXComb(rc->xcomb());
    if ( !rc->dipole()->generateKinematics(r) ) { return false; }
    rc->nodeME()->setXComb(rc->xcomb());
    
    if(MH()->gamma() == 1.){
      if(!MH()->matrixElementRegion( inOutPair.first ,
                                     inOutPair.second ,
                                     rc->pT() ,
                                    MH()->mergePt() ) ){
        return false;
      }
      
    }else{
        
          // Different  treatment if gamma is not 1.
          // Since the dipoles need to be calculated always
          // their alpha region is touched.
          // AlphaRegion != MERegion !!!
        bool inAlphaPS = false;
        for (auto const & ch: thechildren) {
          ch->dipole()->setXComb(ch->xcomb());
          if ( !ch->dipole()->generateKinematics(r) )
            return false;
          MH()->treefactory()->setAlphaParameter( MH()->gamma() );
          inAlphaPS |= ch->dipole()->aboveAlpha();
          MH()->treefactory()->setAlphaParameter( 1. );
        }
        NodePtr rc = randomChild();
        if(!inAlphaPS&&
           !MH()->matrixElementRegion( inOutPair.first ,
                                      inOutPair.second ,
                                      rc->pT() ,
                                      MH()->mergePt()    ) )
             return false;
             
           
        
      
    }
  }
  for (auto const & ch: thechildren) {
    ch->dipole()->setXComb(ch->xcomb());
    if ( !ch->dipole()->generateKinematics(r) ) { cout<<"\nCould not generate dipole kinematics";;return false; }
    if( ! ch->generateKinematics(r,directCut) )return false;
  }
  
  
  return true;
}


StdXCombPtr Node::xcomb() const {
  assert(thexcomb);
  return thexcomb;
}

StdXCombPtr Node::xcomb(){
  if(thexcomb)return thexcomb;
  assert(parent());
  thexcomb=dipole()->makeBornXComb(parent()->xcomb());
  xcomb()->head(parent()->xcomb());
  dipole()->setXComb(thexcomb);
  return thexcomb;
}


void Node::setXComb(tStdXCombPtr xc) {
  assert ( !parent() );
  thexcomb=xc;
  assert(thexcomb->lastParticles().first);
}

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
void Node::birth(const vector<MatchboxMEBasePtr> & vec) {
    // produce the children
  vector<SubtractionDipolePtr> dipoles  = 
           nodeME()->getDipoles(DipoleRepository::dipoles(
                        nodeME()->factory()->dipoleSet()), vec, true);
  
  for ( auto const & dip : dipoles ) {
    dip->doSubtraction();
    NodePtr node = new_ptr(Node(theDeepHead, 
                                       this, 
                                       dip, 
                                       dip->underlyingBornME(), 
                                       theDeepHead->cutStage()));
    thechildren.push_back(node);
    
  }
}


vector<NodePtr> Node::getNextOrderedNodes(bool normal, double hardScaleFactor) const {

  vector<NodePtr> temp = children();
  vector<NodePtr> res;

  for (NodePtr  const & child : children()) {
    if(deepHead()->MH()->mergePt()>child->pT()) {
      res.clear();
      return res;
    }
  }
  
  for (NodePtr const & child: children()) {

    if (parent()&& normal) {
      if ( child->pT() < pT() ) {
        continue;
      }
    }
    if ( child->children().size() != 0 ) {
      
      for (NodePtr itChild:  child->children()) {
        if( itChild->pT() > child->pT()&&child->inShowerPS(itChild->pT()) ) {
          res.push_back(child);
          break;
        }
      }
    }
    else {
      const auto sc=child->nodeME()->factory()->scaleChoice();
      sc->setXComb(child->xcomb());
      if ( sqr(hardScaleFactor)*
           sc->renormalizationScale() >= sqr(child->pT()) &&
          child->inShowerPS(hardScaleFactor*sqrt(sc->renormalizationScale()))) {
        res.push_back(child);
      }
    }
  }
  return res;
}

bool Node::inShowerPS(Energy hardpT)const {
    // Here we decide if the current phase space
    // point can be reached from the underlying Node.
  
  // Full phase space available -> Tilde Kinematic is always fine.
  if(deepHead()->MH()->openZBoundaries()==1)
      return true;
  
  double z_ = dipole()->lastZ();
  // restrict according to hard scale
  if(deepHead()->MH()->openZBoundaries()==0){
      pair<double, double> zbounds =
      dipole()->tildeKinematics()->zBounds(pT(), hardpT);
      return (zbounds.first<z_&&z_<zbounds.second);
  }
  assert(false);
  return false;
}






NodePtr Node::getHistory(bool normal, double hardScaleFactor) {
  NodePtr res = this;

  vector<NodePtr> temp = getNextOrderedNodes(normal, hardScaleFactor);

  Energy minpt = Constants::MaxEnergy;
  Selector<NodePtr> subprosel;
  while (temp.size() != 0) {
    minpt = Constants::MaxEnergy;
    subprosel.clear();    
    for (NodePtr const &  child : temp) {
      assert(deepHead()->MH()->largeNBasis());
      if( child->dipole()->underlyingBornME()->largeNColourCorrelatedME2(
                             {child->dipole()->bornEmitter(),
                              child->dipole()->bornSpectator()},
                              deepHead()->MH()->largeNBasis()) != 0.
        ) {
        
        double weight = 1.;
	if ( deepHead()->MH()->chooseHistory() == 0 )
	   weight = abs(child->dipole()->dSigHatDR()/nanobarn);
	else if ( deepHead()->MH()->chooseHistory() == 1  )
	   weight = abs(child->dipole()->dSigHatDR()/child->nodeME()->dSigHatDRB());
	else if ( deepHead()->MH()->chooseHistory() == 2  )
	   weight = 1.;
        else if ( deepHead()->MH()->chooseHistory() == 3  )
	   weight = 1_GeV/child->pT();
	else 
	   assert(false);
        if(weight != 0.) {
          subprosel.insert(weight , child);
          minpt = min(minpt, child->pT());
        }
      }
    }
    if (subprosel.empty())
      return res;
    
    res = subprosel.select(UseRandom::rnd());
    temp = res->getNextOrderedNodes(true, hardScaleFactor);
  }
  return res;
}


pair<CrossSection, CrossSection> Node::calcDipandPS(Energy scale)const {
  return  dipole()->dipandPs(sqr(scale), deepHead()->MH()->largeNBasis());
}

CrossSection Node::calcPs(Energy scale)const {
  return dipole()->ps(sqr(scale), deepHead()->MH()->largeNBasis());
}

CrossSection Node::calcDip(Energy scale)const {
  return dipole()->dip(sqr(scale));
}


IBPtr Node::clone() const {
  return new_ptr(*this);
}

IBPtr Node::fullclone() const {
  return new_ptr(*this);
}

#include "ThePEG/Persistency/PersistentOStream.h"
void Node::persistentOutput(PersistentOStream & os) const {
  
  os <<
  thexcomb<<
  thenodeMEPtr<<
  thedipol<<
  thechildren<<
  theparent<<
  theProjector<<
  theDeepHead<<
  theCutStage<<
  ounit(theRunningPt, GeV)<<
  theSubtractedReal<<
  theVirtualContribution<<
  theMergingHelper;
}

#include "ThePEG/Persistency/PersistentIStream.h"
void Node::persistentInput(PersistentIStream & is, int) {
  
  is >>
  thexcomb>>
  thenodeMEPtr>>
  thedipol>>
  thechildren>>
  theparent>>
  theProjector>>
  theDeepHead>>
  theCutStage>>
  iunit(theRunningPt, GeV)>>
  theSubtractedReal>>
  theVirtualContribution>>
  theMergingHelper;
}



  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<Node, Interfaced>
 describeHerwigNode("Herwig::Node", "HwDipoleShower.so");

void Node::Init() {
  
  static ClassDocumentation<Node>
  documentation("There is no documentation for the Node class");
  
}

