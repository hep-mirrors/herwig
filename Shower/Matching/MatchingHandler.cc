// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchingHandler class.
//

#include "MatchingHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/PDF/PDF.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/PDF/HwRemDecayer.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

namespace {
struct ParticleOrdering {
  bool operator()(tcPDPtr p1, tcPDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};
}

MatchingHandler::MatchingHandler(bool reWeight) 
  : reWeight_(reWeight), rejectNonAO_( true ), rejectNOHistories_( true ),
    includeDecays_(false)
{}

void MatchingHandler::persistentOutput(PersistentOStream & os) const {
  os << alphaS_ << matrixElement_ << HWmatrixElement_ << includeDecays_
     << partonExtractor_ << cuts_ << rejectNonAO_ << rejectNOHistories_ 
     << allowedInitial_ << allowedFinal_;
}

void MatchingHandler::persistentInput(PersistentIStream & is, int) {
  is >> alphaS_ >> matrixElement_ >> HWmatrixElement_ >> includeDecays_
     >> partonExtractor_ >> cuts_ >> rejectNonAO_ >> rejectNOHistories_
     >> allowedInitial_ >> allowedFinal_;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchingHandler,ShowerHandler>
describeHerwigMatchingHandler("Herwig::MatchingHandler", "HwMatching.so");

void MatchingHandler::Init() {

  static ClassDocumentation<MatchingHandler> documentation
    ("The MatchingHandler class is the base class implementating"
     " many of the features needed for matching.");

  static Reference<MatchingHandler,MEBase> interfaceMatrixElement
    ("MatrixElement",
     "The matrix element class for the core 2->2 process",
     &MatchingHandler::matrixElement_, false, false, true, false, false);

  static Reference<MatchingHandler,PartonExtractor> interfacePartonExtractor
    ("PartonExtractor",
     "The PartonExtractor object used to construct remnants. If no object is "
     "provided the LesHouchesEventHandler object must provide one instead.",
     &MatchingHandler::partonExtractor_, true, false, true, false, false);

  static Reference<MatchingHandler,Cuts> interfaceCuts
    ("Cuts",
     "The Cuts object to be used for this reader. Note that these "
     "must not be looser cuts than those used in the actual generation. "
     "If no object is provided the LesHouchesEventHandler object must "
     "provide one instead.",
     &MatchingHandler::cuts_, true, false, true, false, false);

  static Reference<MatchingHandler,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &MatchingHandler::alphaS_, false, false, true, true, false);

 static Switch<MatchingHandler,bool> interfaceReject
    ("RejectNonOrdered",
     "Whether to reject non angular ordered cluster histories",
     &MatchingHandler::rejectNonAO_, true, false, false);
 static SwitchOption interfaceRejectYes
    (interfaceReject,
     "Reject",
     "Reject all non angular ordered events",
     true);
  static SwitchOption interfaceRejectNo
    (interfaceReject,
     "Select",
     "Select a history for non angular ordered events",
     false);

  static Switch<MatchingHandler,bool> interfaceRejectNoHist
    ("RejectNoHistories",
     "Whether to reject events with no shower interpretation",
     &MatchingHandler::rejectNOHistories_, true, false, false);
  static SwitchOption interfaceRejectNoHistYes
    (interfaceRejectNoHist,
     "Reject",
     "Reject all events with no shower interpretation",
     true);
  static SwitchOption interfaceRejectNoHistNo
    (interfaceRejectNoHist,
     "Shower",
     "Shower events with no shower interpretation directly",
     false);

  static Switch<MatchingHandler,bool> interfaceDecayingParticles
    ("IncludeDecayingParticles",
     "Separate production of decay of unstable particles",
     &MatchingHandler::includeDecays_, false, false, false);
  static SwitchOption interfaceDecayingParticlesYes
    (interfaceDecayingParticles,
     "Yes",
     "Separate them",
     true);
  static SwitchOption interfaceDecayingParticlesNo
    (interfaceDecayingParticles,
     "No",
     "Don't separate them",
     false);

}

void MatchingHandler::doinit() {
  ShowerHandler::doinit();
  HWmatrixElement_ = dynamic_ptr_cast<HwMEBasePtr>(matrixElement_);
  // check 
  if(reWeight_ && !alphaS_) {
    throw Exception() << "ShowerAlpha must be set in MatchingHandler if "
		      << "reweighting events"
		      << Exception::runerror;
  }
  // extract the allowed branchings
  // final-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->finalStateBranchings().begin();
      it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
    pair<long,long> prod(make_pair(it->second.second[1],it->second.second[2]));
    allowedFinal_.insert(make_pair(prod,it->second));
    swap(prod.first,prod.second);
    allowedFinal_.insert(make_pair(prod,it->second));
  }
  // initial-state
  for(BranchingList::const_iterator 
	it = evolver()->splittingGenerator()->initialStateBranchings().begin();
      it != evolver()->splittingGenerator()->initialStateBranchings().end(); ++it) {
    allowedInitial_.insert(make_pair(it->second.second[0],it->second));
  }
}

void MatchingHandler::fillProtoTrees( ProtoTreePtr currentProtoTree ) {
  if(currentProtoTree->branchings().size()==4) return;
  for( set<tProtoBranchingPtr>::const_iterator 
	 ita = currentProtoTree->branchings().begin();
       ita!=currentProtoTree->branchings().end();++ita) {
    for( set<tProtoBranchingPtr>::const_iterator 
	   itb = currentProtoTree->branchings().begin();
	 itb!=ita;++itb) {
      // can't merge two incoming branchings
      if( (**ita).status() == HardBranching::Incoming && 
	  (**itb).status() == HardBranching::Incoming ) continue;
      // get a new branching for this pair
      ProtoBranchingPtr currentBranching = getCluster(*ita,*itb);
      if( ! currentBranching ) continue;
      // branching allowed so make a new Tree out of these branchings
      set< tProtoBranchingPtr > newTreeBranchings = currentProtoTree->branchings();
      newTreeBranchings.erase(*ita);
      newTreeBranchings.erase(*itb);
      newTreeBranchings.insert(currentBranching);
      ProtoTreePtr newProtoTree = new_ptr( ProtoTree( newTreeBranchings ) );
      // remove duplicate trees
      if( ! repeatProtoTree( newProtoTree ) ) protoTrees().insert( newProtoTree );
      // remove the current tree if it hasn't already been removed
      if( protoTrees().find( currentProtoTree ) != protoTrees().end() )
	protoTrees().erase( currentProtoTree );
      // do recursion
      fillProtoTrees( newProtoTree );
    }
  }
}

tProtoBranchingPtr MatchingHandler::getCluster( tProtoBranchingPtr b1,
					    tProtoBranchingPtr b2 ) {
  //look for the clustered pair in protoBranchings_
  for(set<ProtoBranchingPtr>::const_iterator cit = protoBranchings().begin();
      cit != protoBranchings().end(); ++cit) {
    // both outgoing
    if(b1->status()==HardBranching::Outgoing &&
       b2->status()==HardBranching::Outgoing) {
      if((**cit).status()!=HardBranching::Outgoing||
	 (**cit).children().empty()) continue;
      if( ( b1 == (**cit).children()[0] && b2 == (**cit).children()[1] ) ||
	  ( b1 == (**cit).children()[1] && b2 == (**cit).children()[0] ) )
	return *cit;
    }
    // first incoming
    else if(b1->status()==HardBranching::Incoming) {
      if((**cit).backChildren().empty() ) continue;
      if(b1!=(**cit).backChildren()[0]) continue;
      if(b2==(**cit).backChildren()[1]) return *cit;
    }
    // second incoming
    else if(b2->status()==HardBranching::Incoming) {
      if((**cit).backChildren().empty() ) continue;
      if(b2!=(**cit).backChildren()[0]) continue;
      if(b1==(**cit).backChildren()[1]) return *cit;
    }
  }
  // is branching incoming or outgoing
  bool incoming = b1->status()==HardBranching::Incoming ||
                  b2->status()==HardBranching::Incoming;
  // get the branching
  BranchingElement theBranching;
  if( !incoming ) theBranching =   allowedFinalStateBranching( b1, b2 );
  else            theBranching = allowedInitialStateBranching( b1, b2 );

  //if branching is not allowed return null ProtoBrancing
  if( !theBranching.first ) return ProtoBranchingPtr();

  // get the PArticleData object for the new branching
  tcPDPtr particle_data = incoming ?
    getParticleData( theBranching.second[1] ) : getParticleData( theBranching.second[0] );

  // create clustered ProtoBranching
  ProtoBranchingPtr clusteredBranch;
  // outgoing
  if( !incoming ){
    Lorentz5Momentum pairMomentum = b1->momentum() + b2->momentum();
    pairMomentum.setMass(ZERO);
    clusteredBranch = new_ptr(ProtoBranching(particle_data,HardBranching::Outgoing,
					     pairMomentum, theBranching.first));
  }
  // incoming
  else {
    Lorentz5Momentum pairMomentum = b1->momentum() - b2->momentum();
    pairMomentum.setMass( ZERO );
    // check for CC
    if( particle_data->CC() &&
	( b1->id() != theBranching.second[0] ||
	  b2->id() != theBranching.second[2] ) ) {
      particle_data = particle_data->CC();
    }
    clusteredBranch = new_ptr(ProtoBranching(particle_data,HardBranching::Incoming,
					     pairMomentum,theBranching.first));
  }
  protoBranchings().insert(clusteredBranch);
  //set children relations 
  // outgoing
  if( !incoming ){
    clusteredBranch->addChild( b1 );	    
    clusteredBranch->addChild( b2 );  
  }
  else {
    clusteredBranch->addBackChild( b1 );	    
    clusteredBranch->addBackChild( b2 );
  }
  return clusteredBranch;
}

BranchingElement MatchingHandler::
allowedFinalStateBranching( tProtoBranchingPtr & b1, tProtoBranchingPtr & b2) {
  // check with normal ID's
  pair< long, long > ptest = make_pair( b1->id(), b2->id() );
  map< pair< long, long >, pair< SudakovPtr, IdList > >::const_iterator 
    split = allowedFinal_.find(ptest);
  if( split != allowedFinal_.end() ) {
    if(  split->second.second[1] != ptest.first ) swap( b1, b2 );
    return split->second;
  }
  // check with CC
  if( b1->particle()->CC() ) ptest.first  *= -1;
  if( b2->particle()->CC() ) ptest.second *= -1;
  split = allowedFinal_.find( ptest );
  if( split != allowedFinal_.end() ) {
    // cc the idlist only be for qbar g clusterings
    BranchingElement ccBranch = split->second;
    if( getParticleData( ccBranch.second[0] )->CC() ) ccBranch.second[0] *= -1;
    if( getParticleData( ccBranch.second[1] )->CC() ) ccBranch.second[1] *= -1;
    if( getParticleData( ccBranch.second[2] )->CC() ) ccBranch.second[2] *= -1;
    if( split->second.second[1] !=  ptest.first ) swap( b1, b2);
    return ccBranch;
  }
  // not found found null pointer
  return make_pair( SudakovPtr(), IdList() );
}

BranchingElement MatchingHandler::
allowedInitialStateBranching( tProtoBranchingPtr & b1,
			      tProtoBranchingPtr & b2) {
  if(b2->status()==HardBranching::Incoming) swap(b1,b2);
  // is initial parton an antiparticle
  bool cc = b1->id() < 0;
  //gives range of allowedInitial_ with matching first abs( id )
  pair< multimap< long, pair< SudakovPtr, IdList > >::const_iterator,
    multimap< long, pair< SudakovPtr, IdList > >::const_iterator >
    location = allowedInitial_.equal_range( abs( b1->id() ) );
  //iterates over this range
  for( multimap< long, pair< SudakovPtr, IdList> >::const_iterator it = location.first;
       it != location.second; ++it ) {
    //test id for second particle in pair
    long idtest = it->second.second[2];
    //if it is antiparticle *= -1
    if( cc && getParticleData( idtest )->CC() ) idtest *= -1;
    // does second id match the test
    if( idtest == b2->id() ) return it->second;
    //if the the IS parton is a gluon and charge conjugate of second parton mathes accept
    if( idtest == -b2->id() &&
        ! b1->particle()->CC() ) return it->second;
  }
  // not found found null pointer
  return make_pair(SudakovPtr(),IdList());
}


bool MatchingHandler::repeatProtoTree( ProtoTreePtr currentProtoTree ) {
  // loop over all prototrees and see 
  // how many ProtoBranchings of current ProtoTree are found in each
  for( set< ProtoTreePtr >::const_iterator cit = protoTrees().begin();
       cit != protoTrees().end(); ++cit ) {
    unsigned int no_matches = 0;
    for( set< tProtoBranchingPtr >::const_iterator ckt 
	   = currentProtoTree->branchings().begin(); 
	 ckt != currentProtoTree->branchings().end(); ckt++ ) {
      if( (*cit)->branchings().find( *ckt ) != (*cit)->branchings().end() )
	++no_matches;
    }
    // return true if all match
    if( no_matches == currentProtoTree->branchings().size() )
      return true;
  }
  return false;
}

double MatchingHandler::getDiagram(PotentialTree & tree) {
  if(tree.diagrams().empty()) {
    set<HardBranchingPtr>::const_iterator cit;
    tcPDPair incoming;
    multiset<tcPDPtr,ParticleOrdering> outgoing;  
    //get the incoming and outgoing partons involved in hard process
    for( cit = tree.tree()->branchings().begin(); 
	 cit != tree.tree()->branchings().end(); ++cit ){ 
      if( (*cit)->status() ==HardBranching::Incoming) {
	HardBranchingPtr parent = *cit;
	while(parent->parent()) parent = parent->parent();
	if( parent->branchingParticle()->momentum().z()>ZERO )
	  incoming.first  = (*cit)->branchingParticle()->dataPtr();
	else
	  incoming.second = (*cit)->branchingParticle()->dataPtr();
      }
      else {
	outgoing.insert( (*cit)->branchingParticle()->dataPtr() );
      }
    }
    if(!incoming.first || !incoming.second) return 0.;
    pair<string,string> tag;
    tag.first  = incoming.first  ->PDGName() + "," + incoming.second->PDGName() + "->";
    tag.second = incoming.second ->PDGName() + "," + incoming.first ->PDGName() + "->";

    string tag_out;
    for ( multiset<tcPDPtr,ParticleOrdering>::iterator i = outgoing.begin();
	  i != outgoing.end(); ++i ) {
      if ( i != outgoing.begin() ) tag_out += ",";
      tag_out += (**i).PDGName();
    }
    tag.first  += tag_out;
    tag.second += tag_out;
    // find the diagrams
    for ( int i = 0, N = matrixElement()->diagrams().size(); i < N; ++i ) {
      string temp = matrixElement()->diagrams()[i]->getTag();
      if ( temp == tag.first || temp == tag.second )
	tree.diagrams().push_back(matrixElement()->diagrams()[i]);
    }
  }
  if(tree.diagrams().empty()) return 0.;
  // construct a set of on-shell momenta for the hard collison
  vector<Lorentz5Momentum> meMomenta;
  vector<tcPDPtr> mePartonData;
  PVector particles;
  set<HardBranchingPtr>::const_iterator it;
  // incoming particles
  for( it = tree.tree()->branchings().begin();it != tree.tree()->branchings().end(); ++it ) {
    if( (**it).status() == HardBranching::Incoming ) {
      meMomenta.push_back( (**it).branchingParticle()->momentum() );
      mePartonData.push_back( (**it).branchingParticle()->dataPtr() );
      particles.push_back( (**it).branchingParticle() );
    }
  }
  assert(particles.size()==2);
  for( it = tree.tree()->branchings().begin(); it != tree.tree()->branchings().end(); ++it ) {
    if( (**it).status() == HardBranching::Outgoing ) {
      meMomenta.push_back( (**it).branchingParticle()->momentum() );
      mePartonData.push_back( (**it).branchingParticle()->dataPtr() );
      particles.push_back( (**it).branchingParticle() );
    }
  }
  const cPDVector partons = tree.diagrams()[0]->partons();
  // order of the incoming partons
  if(mePartonData[0] != partons[0]) {
    swap( mePartonData[0], mePartonData[1] );
    swap( meMomenta[0], meMomenta[1] );
    swap( particles[0], particles[1] );
  }
  // order of the outgoing partons
  for(unsigned int ix=2;ix<partons.size();++ix) {
    for(unsigned int iy=ix;iy<meMomenta.size();++iy) {
      if(partons[ix]==mePartonData[iy]) {
	if(ix!=iy) {
	  swap(mePartonData[ix],mePartonData[iy]);
	  swap(meMomenta[ix],meMomenta[iy]);
	  swap(particles[ix],particles[iy]);
	}
	break;
      }
    }
  }
  // boost to the CMF frame
  Lorentz5Momentum prest(meMomenta[0]+meMomenta[1]);
  LorentzRotation R(-prest.boostVector());
  // and then to put beams along the axis
  Lorentz5Momentum ptest = R*meMomenta[0];
  Axis axis(ptest.vect().unit());
  if(axis.perp2()>0.) {
    R.rotateZ(-axis.phi());
    R.rotateY(-acos(axis.z()));
  }
  for( unsigned int ix = 0; ix < meMomenta.size(); ++ix )
    meMomenta[ix].transform(R);
  // now rescale to put on shell
  Energy Ebeam = 0.5 * ( max(meMomenta[0].e(),abs(meMomenta[0].z())) + 
			 max(meMomenta[1].e(),abs(meMomenta[1].z())) );
  for( unsigned int i = 0; i < 2; ++i ) {
    meMomenta[i].setZ( meMomenta[i].z() / abs(meMomenta[i].z()) * Ebeam  );
    meMomenta[i].setT( Ebeam );
  }
  Energy2 s = 4.0 * sqr(Ebeam);
  Energy m1 = mePartonData[2]->mass();
  Energy m2 = mePartonData[3]->mass();
  // \todo need to improve this
  if(m1+m2>sqrt(s)) return 0.;
  double lambda = 0.25/Ebeam/meMomenta[2].rho() * 
    sqrt( ( s - sqr(m1+m2) ) * ( s - sqr(m1-m2) ) );
  for( unsigned int i = 2; i < meMomenta.size(); ++i ) {
    meMomenta[i] *= lambda;
    meMomenta[i].setMass(mePartonData[i]->mass());
    meMomenta[i].rescaleEnergy();
  }
  // incoming pair
  PPair in( mePartonData[0]->produceParticle( meMomenta[0] ),
	    mePartonData[1]->produceParticle( meMomenta[1] ) );
  // outgoing
  PVector out;
  for(unsigned int ix=2;ix<meMomenta.size();++ix) {
    out.push_back(mePartonData[ix]->produceParticle(meMomenta[ix]));
  }
  // call the matrix element to initialize
  matrixElement()->setKinematics(in,out);
  if(HWMatrixElement()) {
    vector<Lorentz5Momentum> momenta;
    cPDVector data;
    data.push_back(in. first->dataPtr());
    momenta.push_back(in. first->momentum());
    data.push_back(in.second->dataPtr());
    momenta.push_back(in.second->momentum());
    for(unsigned int ix=0;ix<out.size();++ix) {
      data.push_back(out[ix]->dataPtr());
      momenta.push_back(out[ix]->momentum());
    }
    HWMatrixElement()->rescaleMomenta(momenta,data);
  }
  if(!cuts()->scale(matrixElement()->scale())) {
    return 0.;
  }
  matrixElement()->dSigHatDR();
  // select the diagram
  if(!tree.diagram())
    tree.diagram(tree.diagrams()[matrixElement()->diagram(tree.diagrams())]);
  // get the colour structure
  if(!tree.colourLines()) {
    Selector<const ColourLines *> sel = matrixElement()->colourGeometries(tree.diagram());
    tree.colourLines(sel.select(rnd()));
  }
  PVector slike;
  tPVector ret;
  slike.push_back(particles[0]);
  Ptr<Tree2toNDiagram>::pointer diagram2 = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::pointer>(tree.diagram());
  for ( int i = 1; i < diagram2->nSpace() - 1; ++i )
    slike.push_back(diagram2->allPartons()[i]->produceParticle());
  slike.push_back(particles[1]);
  ret = tPVector(slike.begin(), slike.end());
  int io = particles.size();
  PVector tlike(diagram2->allPartons().size() - diagram2->nSpace());
  for ( int i = diagram2->allPartons().size() - 1; i >=  diagram2->nSpace(); --i ) {
    int it = i - diagram2->nSpace();
    pair<int,int> ch = diagram2->children(i);
    bool iso = ch.first < 0;
    if ( iso ) {
      tlike[it] = particles[--io];
    } 
    else {
      Lorentz5Momentum p = tlike[ch.first - diagram2->nSpace()]->momentum() +
 	tlike[ch.second - diagram2->nSpace()]->momentum();
      tlike[it] = diagram2->allPartons()[i]->produceParticle(p);
    }
  }
  ret.insert( ret.end(), tlike.begin(), tlike.end() );
  tree.colourLines()->connect(ret);
  for( unsigned int ix = 0; ix < ret.size(); ++ix ) {
    PVector::iterator it = find( particles.begin(), particles.end(),ret[ix] );
    if( it == particles.end() ) {
      ColinePtr line = ret[ix]->colourLine();
      if(line) line->removeColoured(ret[ix]);
      line = ret[ix]->antiColourLine();
      if(line) line->removeAntiColoured(ret[ix]);
    }
  }
  // now the colours of the rest of the particles
  for( set<HardBranchingPtr>::const_iterator it = tree.tree()->branchings().begin();
       it!=tree.tree()->branchings().end(); ++it ) (**it).fixColours();
  // now make the colour connections in the tree
  ShowerParticleVector branchingParticles;
  map<ShowerParticlePtr,HardBranchingPtr> branchingMap;
  for( set< HardBranchingPtr >::iterator it = tree.tree()->branchings().begin();
       it != tree.tree()->branchings().end(); ++it ) {
    branchingParticles.push_back((**it).branchingParticle());
    branchingMap.insert(make_pair((**it).branchingParticle(),*it));
  }
  // find the colour partners
  evolver()->showerModel()->partnerFinder()
    ->setInitialEvolutionScales(branchingParticles,false,
				ShowerInteraction::QCD,true);
  for(unsigned int ix=0;ix<branchingParticles.size();++ix) {
    if(branchingParticles[ix]->partner()) {
      HardBranchingPtr partner = branchingMap[branchingParticles[ix]->partner()];
      branchingMap[branchingParticles[ix]]->colourPartner(partner);
    }
  }
  double weight = 1.;
  // calculate the weight if needed
  if(reWeight_) {
    if(matrixElement()->orderInAlphaS()>0) {
      weight = pow(alphaS_->value( max(matrixElement()->scale(),sqr(pdfScale_)) )
		   / alphaSMG_, int(matrixElement()->orderInAlphaS()));
    }
  }
  // return the weight
  return weight;
}

bool MatchingHandler::updateSubProcess() {
  assert( hardTree().diagram() && hardTree().tree() );
  PPair beams = lastXCombPtr()->lastParticles();
  // remove children of beams
  PVector beam_children = beams.first->children();
  for( unsigned int ix = 0; ix != beam_children.size(); ix++ )
    beams.first->abandonChild( beam_children[ix] );
  beam_children = beams.second->children();
  for( unsigned int ix = 0; ix != beam_children.size(); ix++ )
    beams.second->abandonChild( beam_children[ix] );
  if( (**hardTree().tree()->incoming().begin()).branchingParticle()->momentum().z() /
      beams.first->momentum().z() < 0.)
    swap( beams.first, beams.second );
  Ptr<Tree2toNDiagram>::pointer diagram = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::pointer>(hardTree().diagram());
  assert(diagram);
  set<HardBranchingPtr>::const_iterator it; 
  map< ColinePtr, ColinePtr> colourMap;
  PPair incoming;
  // loop over the branchings and sort out incoming particles
  for( it = hardTree().tree()->branchings().begin();
       it != hardTree().tree()->branchings().end(); ++it) {
    if( (*it)->status() == HardBranching::Outgoing ) continue;
    PPtr newParticle = new_ptr( Particle( (**it).branchingParticle()->dataPtr() ) );
    newParticle->set5Momentum( (**it).showerMomentum() );
    if( (**it).branchingParticle()->colourLine() ) {
      map< ColinePtr, ColinePtr>::iterator loc 
	= colourMap.find( (**it).branchingParticle()->colourLine() );
      if( loc != colourMap.end() ) {
	loc->second->addColoured( newParticle );
      }
      else {
	ColinePtr newLine = new_ptr( ColourLine() );
	colourMap[ (**it).branchingParticle()->colourLine() ] = newLine;
	newLine->addColoured( newParticle );
      }
    }
    if( (**it).branchingParticle()->antiColourLine() ) {
      map< ColinePtr, ColinePtr> ::iterator loc 
	= colourMap.find( (**it).branchingParticle()->antiColourLine() );
      if( loc != colourMap.end() ) {
	loc->second->addAntiColoured( newParticle );
      }
      else {
	ColinePtr newLine = new_ptr( ColourLine() );
	colourMap[ (**it).branchingParticle()->antiColourLine() ] = newLine;
	newLine->addAntiColoured( newParticle );
      }
    }
    if( lastXCombPtr()->subProcess()->incoming().first->momentum().z() /
	newParticle->momentum().z() > 0. )
      incoming.first  = newParticle;
    else
      incoming.second = newParticle;
  }
  bool mirror = 
    incoming.first ->dataPtr() != diagram->partons()[0] &&
    incoming.second->dataPtr() != diagram->partons()[1];
  // create the new subprocess
  SubProPtr newSubProcess =
    new_ptr( SubProcess( incoming, lastXCombPtr()->subProcess()->collision(),
			 lastXCombPtr()->subProcess()->handler() ) );
  // add the spacelike intermediates
  PVector slike;
  slike.push_back( !mirror ? incoming.first : incoming.second);
  for ( int i = 1; i < diagram->nSpace() - 1; ++i )
    slike.push_back(diagram->allPartons()[i]->produceParticle());
  slike.push_back( !mirror ? incoming.second : incoming.first);
  tPVector ret = tPVector(slike.begin(), slike.end());
  for ( size_t i = 1; i < slike.size() - 1; ++i ) {
    slike[i-1]->addChild(slike[i]);
    newSubProcess->addIntermediate(slike[ mirror ? i: slike.size() - 1 - i], false);
  }
  // get the parton bins from the parton extractor if first time
  static bool first = true;
  if(first) {
    first = false;
    Energy e1 = lastXCombPtr()->lastParticles().first ->momentum().e();
    Energy e2 = lastXCombPtr()->lastParticles().second->momentum().e();
    Energy emax = 2.0*sqrt(e1*e2);
    cPDPair inData = make_pair(lastXCombPtr()->lastParticles().first ->dataPtr(),
			       lastXCombPtr()->lastParticles().second->dataPtr());
    cuts()->initialize(sqr(emax),0.5*log(e1/e2));
    partonBins(partonExtractor()->getPartons(emax, inData, *cuts()));
  }
  // get the parton bins for this event
  tcPBPair sel;
  for ( int i = 0, N = partonBins().size(); i < N; ++i ) {
    tcPBPtr bin = partonBins()[i].first;
    tPPtr p = incoming.first;
    while ( bin && p ) {
      if ( p->dataPtr() != bin->parton() ) break;
      bin = bin->incoming();
      p = p != lastXCombPtr()->lastParticles().first ?
	lastXCombPtr()->lastParticles().first : PPtr();
    }
    if ( bin || p ) continue;
    bin = partonBins()[i].second;
    p = incoming.second;
    while ( bin && p ) {
      if ( p->dataPtr() != bin->parton() ) break;
      bin = bin->incoming();
      p = p != lastXCombPtr()->lastParticles().second ?
	lastXCombPtr()->lastParticles().second : PPtr();
    }
    if ( bin || p ) continue;
    sel = partonBins()[i];
    break;
  }
  if ( !sel.first || !sel.second ) Throw<Exception>()
				     << "Could not find appropriate "
				     << "PartonBin objects for event in "
				     << "MatchingHandler " << Exception::runerror;
  // create the new parton bin instances
  Direction<0> dir(true);
  PBIPair partonBinInstances;
  // temporary mother/daugther settings
  // to get parton bin instances correct
  lastXCombPtr()->lastParticles().first ->addChild(incoming.first );
  lastXCombPtr()->lastParticles().second->addChild(incoming.second);
  // make the parton bin instances
  partonBinInstances.first =
    new_ptr(PartonBinInstance(incoming.first, sel.first,
			      lastXCombPtr()->partonBinInstances().first->scale()));
  dir.reverse();
  partonBinInstances.second =
    new_ptr(PartonBinInstance(incoming.second, sel.second,
			      lastXCombPtr()->partonBinInstances().second->scale()));
  // remove temporary mother/daugther settings
  lastXCombPtr()->lastParticles().first ->abandonChild(incoming.first );
  lastXCombPtr()->lastParticles().second->abandonChild(incoming.second);
  // set the parton bin instances
  lastXCombPtr()->setPartonBinInstances(partonBinInstances,
					lastXCombPtr()->lastScale());
  // momenta of the time like partons
  vector<Lorentz5Momentum> meMomenta;
  vector<tcPDPtr> mePartonData;
  vector<HardBranchingPtr> branchings;
  for( it = hardTree().tree()->branchings().begin(); 
       it != hardTree().tree()->branchings().end(); ++it ) {
    if( (**it).status() == HardBranching::Outgoing ) {
      meMomenta.push_back( (**it).showerMomentum() );
      mePartonData.push_back( (**it).branchingParticle()->dataPtr() );
      branchings.push_back(*it);
    }
  }
  // order of the outgoing partons
  for(int ix=2;ix<int(diagram->partons().size());++ix) {
    for(int iy=ix-2;iy<int(meMomenta.size());++iy) {
      if(diagram->partons()[ix]==mePartonData[iy]) {
	if(ix!=iy) {
	  swap(mePartonData[ix-2],mePartonData[iy]);
	  swap(meMomenta   [ix-2],meMomenta   [iy]);
	  swap(branchings  [ix-2],branchings  [iy]);
	}
	break;
      }
    }
  }
  // time like particles
  int io = meMomenta.size();
  PVector tlike(diagram->allPartons().size() - diagram->nSpace());
  vector<HardBranchingPtr> tBranchings;
  tPVector out;
  for ( int i = diagram->allPartons().size() - 1; i >=  diagram->nSpace(); --i ) {
    int it = i - diagram->nSpace();
    pair<int,int> ch = diagram->children(i);
    bool iso = ch.first < 0;
    if ( iso ) {
      tlike[it] = diagram->allPartons()[i]->produceParticle(meMomenta[--io]);
    } 
    else {
      Lorentz5Momentum p = tlike[ch.first - diagram->nSpace()]->momentum() +
	tlike[ch.second - diagram->nSpace()]->momentum();
      tlike[it] = diagram->allPartons()[i]->produceParticle(p);
    }
    if ( diagram->parent(i) < diagram->nSpace() ) {
      slike[diagram->parent(i)]->addChild(tlike[it]);
      if ( diagram->parent(i) == diagram->nSpace() - 2 )
	slike[diagram->parent(i) + 1]->addChild(tlike[it]);
    }
    if ( !iso ) {
      tlike[it]->addChild(tlike[ch.first - diagram->nSpace()]);
      tlike[it]->addChild(tlike[ch.second - diagram->nSpace()]);
    }
    if ( iso )  {
      out.push_back(tlike[it]);
      tBranchings.push_back(branchings[io]);
    }
    else        
      newSubProcess->addIntermediate(tlike[it], false);
  }
  ret.insert(ret.end(), tlike.begin(), tlike.end());
  // select the colour structure now   
  const ColourLines & cl = matrixElement()->selectColourGeometry(hardTree().diagram());
  cl.connect(ret);
  // add the particles
  for ( int i = 0, N = out.size(); i < N; ++i ) {
    tPPtr particle = out[mirror ? i: out.size() - i - 1];
    // if not including decays add as outgoing
    if(!includeDecays_) {
      newSubProcess->addOutgoing(particle, false);
      continue;
    }
    HardBranchingPtr branching = tBranchings[mirror ? i: out.size() - i - 1];
    // move to end of chain for time-like branchings
    while(!branching->children().empty()) {
      bool found=false;
      for(unsigned int ix=0;ix<branching->children().size();++ix) {
	if(branching->children()[ix]->branchingParticle()->id()==
	   branching->branchingParticle()->id()) {
	  found = true;
	  branching = branching->children()[ix];
	  break;
	}
      }
      if(!found) break;
    }
    // check if from decay
    map<PPtr,ParticleVector>::const_iterator pit = parent(branching->branchingParticle());
    // if not add as outgoing
    if(pit==decayingParticles_.end()) {
      newSubProcess->addOutgoing(particle, false);
      continue;
    }
    LorentzRotation decayBoost = branching->showerBoost();
    // final boost if FSR
    Lorentz5Momentum pShower = decayBoost*(pit->first->momentum());
    decayBoost = LorentzRotation(-pShower.boostVector())*decayBoost;
    decayBoost = LorentzRotation(particle->momentum().boostVector())*decayBoost;
    // add decayed particle as intermediate
    newSubProcess->addIntermediate(particle,false);
    // add decay products
    addDecayProducts(newSubProcess,particle,pit,decayBoost);
  }
  for ( PVector::size_type i = 0; i < slike.size() - 2; ++i ) {
    pair<int,int> ch = diagram->children(i);
    slike[ch.first]->set5Momentum(slike[i]->momentum() -
				  tlike[ch.second - diagram->nSpace()]->momentum());
  }
  // set the subprocess
  lastXCombPtr()->subProcess( newSubProcess );
  decayingParticles_.clear();
  return true;
}

void MatchingHandler::findDecayingParticles(PPtr parent) {
  ParticleVector decayProducts;
  for(unsigned int ix=0;ix<parent->children().size();++ix) {
      decayProducts.push_back(parent->children()[ix]);
    if(!parent->children()[ix]->children().empty()) 
      findDecayingParticles(parent->children()[ix]);
  }
  if(decayingParticles_.find(parent)==decayingParticles_.end())
    decayingParticles_[parent] = decayProducts;
}

PotentialTree MatchingHandler::doClustering() {
  noShowerHists_ = false;
  // clear storage of the protoTrees
  protoBranchings().clear();
  protoTrees().clear();
  nonOrderedTrees_.clear();
  hardTrees_.clear();
  assert( matrixElement() );
  // get particles from the XComb object
  ParticleVector outgoing  = lastXCombPtr()->subProcess()->outgoing();
  PPair incoming = lastXCombPtr()->subProcess()->incoming();
  // storage of decayed particles
  decayingParticles_.clear();
  // all outgoing particles as a set for checking
  set<PPtr> outgoingset(outgoing.begin(),outgoing.end());
  // loop through the FS particles and create ProtoBranchings
  for( unsigned int i = 0; i < outgoing.size(); ++i) {
    tPPtr parent = outgoing[i]->parents()[0];
    bool decayProd = decayProduct(parent,lastXCombPtr()->subProcess());
    if(!decayProd||!includeDecays_) {
      ProtoBranchingPtr currentBranching =
	new_ptr(ProtoBranching(outgoing[i]->dataPtr(),HardBranching::Outgoing,
			       outgoing[i]->momentum(),tSudakovPtr()));
      protoBranchings().insert(currentBranching);
    }
    else {
      bool isHard = true;
      PPtr newParent = findParent(parent,isHard,outgoingset,false,
				  lastXCombPtr()->subProcess());
      assert(newParent);
      if(decayingParticles_.find(newParent)==decayingParticles_.end()) {
	ProtoBranchingPtr currentBranching =
	  new_ptr(ProtoBranching(newParent->dataPtr(),HardBranching::Outgoing,
				 newParent->momentum(),tSudakovPtr()));
	protoBranchings().insert(currentBranching);
	findDecayingParticles(newParent);
      }
    }
  }
  // add IS hardBranchings
  ProtoBranchingPtr currentBranching = 
    new_ptr(ProtoBranching(incoming.first ->dataPtr(),HardBranching::Incoming,
			   incoming.first ->momentum(),tSudakovPtr()));
  protoBranchings().insert(currentBranching);
  currentBranching =
    new_ptr(ProtoBranching(incoming.second->dataPtr(),HardBranching::Incoming,
			   incoming.second->momentum(),tSudakovPtr()));
  protoBranchings().insert(currentBranching);
  //create and initialise the first tree
  ProtoTreePtr initialProtoTree = new_ptr( ProtoTree() );
  for(set<ProtoBranchingPtr>::const_iterator it=protoBranchings().begin();
      it!=protoBranchings().end();++it) {
    initialProtoTree->addBranching(*it);
  }
  //fill _proto_trees with all possible trees
  protoTrees().insert(initialProtoTree );
  fillProtoTrees( initialProtoTree );
  double totalWeight = 0., nonOrderedWeight = 0.;
  // create a HardTree from each ProtoTree and fill hardTrees()
  for( set< ProtoTreePtr >::const_iterator cit = protoTrees().begin(); 
       cit != protoTrees().end(); ++cit ) {
    PotentialTree newTree;
    newTree.tree((**cit).createHardTree());
    // check the created CKKWTree corresponds to an allowed LO configuration
    // (does matrix element have a corresponding diagram)
    double meWgt = getDiagram( newTree );
    if( !newTree.diagram() ) continue;
    // set the beam particles
    PPair beams = lastXCombPtr()->lastParticles();
    // remove children of beams
    PVector beam_children = beams.first->children();
    if( (**newTree.tree()->incoming().begin()).branchingParticle()->momentum().z() /
	beams.first->momentum().z() < 0.)
      swap( beams.first, beams.second );
    set<HardBranchingPtr>::iterator it = newTree.tree()->incoming().begin();
    HardBranchingPtr br = *it;
    br->beam( beams.first );
    while ( !br->children().empty() ) {
      for(unsigned int ix = 0; ix < br->children().size(); ++ix ) {
	if( br->children()[ix]->status() == HardBranching::Incoming ) {
	  br = br->children()[ix];
	  break;
	}
      }
      br->beam( beams.first );
    }
    ++it;
    br = *it;
    br->beam( beams.second );
    while ( !br->children().empty() ) {
      for( unsigned int ix = 0; ix < br->children().size(); ++ix ) {
	if( br->children()[ix]->status() == HardBranching::Incoming ) {
	  br = br->children()[ix];
	  break;
	}
      }
      br->beam( beams.second );
    }
    // do inverse momentum reconstruction
    if( !evolver()->showerModel()->kinematicsReconstructor()
    	->deconstructHardJets( newTree.tree(), evolver(), ShowerInteraction::QCD ) )
      continue;
    newTree.tree()->findNodes();
    if( newTree.tree()->ordered() ) {
      // find the wgt and fill hardTrees() map
      double treeWeight = 1.;
      if(reWeight_) treeWeight = meWgt*sudakovWeight( newTree.tree() );
      newTree.weight(treeWeight);
      hardTrees_.push_back( make_pair( newTree, treeWeight ) );
      totalWeight += treeWeight;
    }
    else {
      nonOrderedTrees_.push_back( make_pair( newTree, 1. ) );
      nonOrderedWeight += 1.;
    }
  }
  // if rejecting non-ordered trees return
  if( hardTrees_.empty() && rejectNonAO_  ) return PotentialTree();
  // select the tree
  PotentialTree chosen_hardTree=chooseHardTree(totalWeight,
					       nonOrderedWeight);
  protoBranchings().clear();
  protoTrees().clear();
  nonOrderedTrees_.clear();
  hardTrees_.clear();
  if(! chosen_hardTree.tree() ) {
    noShowerHists_ = true;
    return PotentialTree();
  }
  else
    return chosen_hardTree;
}

void MatchingHandler::initialiseMatching(int minMult, int maxMult) {
  // check multiplicity of FS partons
  int nOut = lastXCombPtr()->subProcess()->outgoing().size();
  // is it the lowest  multiplicity
  lowestMult_  = nOut == minMult;
  // or the highest
  highestMult_ = nOut == maxMult;
  // centre-of-mass energy
  sHat(lastXCombPtr()->lastSHat());
  // scale for the PDFs 
  pdfScale(sqrt(lastXCombPtr()->lastScale()));
  // alphaS value used by the ME generate
  alphaSMG(lastXCombPtr()->lastAlphaS());
  // create a hard tree by clustering the event
  hardTree(doClustering());
}

map<PPtr,ParticleVector>::const_iterator MatchingHandler::parent(PPtr parent) {
  long id = parent->id();
  for(map<PPtr,ParticleVector>::const_iterator it = decayingParticles_.begin();
      it!=decayingParticles_.end();++it) {
    if(id!=it->first->id()) continue;
    Energy2 diff = 
      sqr(it->first->momentum().x()-parent->momentum().x()) +
      sqr(it->first->momentum().y()-parent->momentum().y()) +
      sqr(it->first->momentum().z()-parent->momentum().z()) +
      sqr(it->first->momentum().t()-parent->momentum().t());
    Energy2 sum = 
      sqr(it->first->momentum().x()+parent->momentum().x()) +
      sqr(it->first->momentum().y()+parent->momentum().y()) +
      sqr(it->first->momentum().z()+parent->momentum().z()) +
      sqr(it->first->momentum().t()+parent->momentum().t());
    double ratio = diff/sum;
    if(ratio<1e-10) return it;
  }
  return decayingParticles_.end();
}

void MatchingHandler::addDecayProducts(SubProPtr subProcess, PPtr parent,
				       map<PPtr,ParticleVector>::const_iterator decay,
				       const LorentzRotation & boost) const {
  // map colours of the parent
  map<ColinePtr,ColinePtr> cmap;
  if(parent->colourLine())
    cmap.insert(make_pair(decay->first->    colourLine(),
			  parent->    colourLine()));
  if(parent->antiColourLine())
    cmap.insert(make_pair(decay->first->antiColourLine(),
			  parent->antiColourLine()));
  // add the decay products
  for(unsigned int ix=0;ix<decay->second.size();++ix) {
    Lorentz5Momentum pnew = boost*decay->second[ix]->momentum();
    PPtr newParticle = decay->second[ix]->dataPtr()->produceParticle(pnew);
    parent->addChild(newParticle);
    if(decay->second[ix]->colourLine()) {
      if(cmap.find(decay->second[ix]->colourLine())==cmap.end()) {
	ColinePtr newLine(new_ptr(ColourLine()));
	cmap.insert(make_pair(decay->second[ix]->colourLine(),newLine));
      }
      cmap[decay->second[ix]->colourLine()]->addColoured(newParticle);
    }
    if(decay->second[ix]->antiColourLine()) {
      if(cmap.find(decay->second[ix]->antiColourLine())==cmap.end()) {
	ColinePtr newLine(new_ptr(ColourLine()));
	cmap.insert(make_pair(decay->second[ix]->antiColourLine(),newLine));
      }
      cmap[decay->second[ix]->antiColourLine()]->addAntiColoured(newParticle);
    }
    map<PPtr,ParticleVector>::const_iterator pit = 
      decayingParticles_.find(decay->second[ix]);
    if(pit!=decayingParticles_.end()) {
      subProcess->addIntermediate(newParticle, false);
      addDecayProducts(subProcess,newParticle,pit,boost);
    }
    else {
      subProcess->addOutgoing(newParticle, false);
    }
  }
}
