  // -*- C++ -*-
  //
  // MergeboxFactory.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the MergeboxFactory class.
  //


#include "MergingFactory.h"
#include "Node.h"
#include "ThePEG/Repository/Repository.h"

#include "ThePEG/Utilities/ColourOutput.h"

using namespace Herwig;
using std::ostream_iterator;

IBPtr MergingFactory::clone() const {
  return new_ptr(*this);
}

IBPtr MergingFactory::fullclone() const {
  return new_ptr(*this);
}

void MergingFactory::doinit(){
  MatchboxFactory::doinit();
  if (subProcessGroups()) {
    throw InitException() << "There are no subprocess groups in merging!";
  }
}

void MergingFactory::productionMode() {
  if(M()<0)
  for ( vector<Ptr<MatchboxAmplitude>::ptr>::iterator amp
         = amplitudes().begin(); amp != amplitudes().end(); ++amp ) {
        Repository::clog() << "One-loop contributions from '"
        << (**amp).name()
        << "' are not required and will be disabled.\n"
        << flush;
        (**amp).disableOneLoop();
  }
  MatchboxFactory::productionMode();
}



void MergingFactory::fillMEsMap() {
  
  olpProcesses().clear();
  
  assert( getProcesses().size() == 1 );
  processMap[0] = getProcesses()[0];
  
  if ( MH()->M() >= 0 ) 
    setHighestVirt(processMap[0].size()+MH()->M());
  
  MH()->N0(processMap[0].size());
  for ( int i = 1 ; i <= MH()->N() ; ++i ) {
    processMap[i] = processMap[i - 1];
    processMap[i].push_back("j");
  }
  
  for ( int i = 0 ; i <= MH()->N() ; ++i ) {
  	const bool below_maxNLO = i < MH()->M() + 1;
    vector<MatchboxMEBasePtr> ames 
    	= makeMEs(processMap[i], orderInAlphaS() + i, below_maxNLO );
    copy(ames.begin(), ames.end(), back_inserter(pureMEsMap()[i]));
  }
}


#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
void MergingFactory::prepare_BV(int i) {
    // check if we have virtual contributions
  bool haveVirtuals = true;

  for ( auto born : pureMEsMap()[i]) {
    
    prepareME(born);
     
    if ( born->isOLPTree() ) {
      int id = orderOLPProcess(born->subProcess(),
                               born->matchboxAmplitude(),
                               ProcessType::treeME2);
      born->olpProcess(ProcessType::treeME2,id);
      id = orderOLPProcess(born->subProcess(),
                           born->matchboxAmplitude(),
                           ProcessType::colourCorrelatedME2);
      born->olpProcess(ProcessType::colourCorrelatedME2,id);
      
      bool haveGluon = false;
      for ( const auto & p : born->subProcess().legs )
      if ( p->id() == 21 ) {
        haveGluon = true;
        break;
      }
      if ( haveGluon ) {
        id = orderOLPProcess(born->subProcess(),
                             born->matchboxAmplitude(),
                             ProcessType::spinColourCorrelatedME2);
        born->olpProcess(ProcessType::spinColourCorrelatedME2,id);
      }
    }
    if ( born->isOLPLoop() &&  i <= MH()->M() ) {
      int id = orderOLPProcess(born->subProcess(),
                               born->matchboxAmplitude(),
                               ProcessType::oneLoopInterference);
      born->olpProcess(ProcessType::oneLoopInterference,id);
      if ( !born->onlyOneLoop() && born->needsOLPCorrelators() ) {
        id = orderOLPProcess(born->subProcess(),
                             born->matchboxAmplitude(),
                             ProcessType::colourCorrelatedME2);
        born->olpProcess(ProcessType::colourCorrelatedME2,id);
      }
    }
    haveVirtuals &= born->haveOneLoop();
  }
   
    // check for consistent conventions on virtuals, if we are to include MH()->M()
  
  if (!(i > MH()->M()||haveVirtuals))
  	throw InitException() 
	<< MH()->M()
  	<< " NLO corrections requested,\n"
  	<< "but no virtual contributions are found.";
}

void MergingFactory::prepare_R(int i) {
  for ( auto real : pureMEsMap()[i])
    prepareME(real);
}


#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
void MergingFactory::getVirtuals(MatchboxMEBasePtr nlo, bool clone){
  
  const auto & partons = nlo->diagrams().front()->partons();

  for ( auto I : DipoleRepository::insertionIOperators(dipoleSet()) )
  if ( I->apply(partons) ){
  	auto myI = I;
  	if ( clone ) myI = I->cloneMe();
        nlo->virtuals().push_back(myI);
  }
  
  for ( auto PK : DipoleRepository::insertionPKOperators(dipoleSet()) )
  if ( PK->apply(partons) ){
  	auto myPK = PK;
  	if ( clone ) myPK = PK->cloneMe();
        nlo->virtuals().push_back(myPK);
  }
  
}



void MergingFactory::pushB(MatchboxMEBasePtr born, int i) {
  MatchboxMEBasePtr bornme = born->cloneMe();
  bornme->maxMultCKKW(1);
  bornme->minMultCKKW(0);
  
  
  string pname = fullName() + "/" + bornme->name() + ".Born";
  if ( !(generator()->preinitRegister(bornme, pname)) ) 
  	throw InitException() 
  	<< "Born ME "<< pname << " already existing.";
  
  if (MH()->gamma()!=1.)
    getVirtuals(bornme,false);
  
  
  NodePtr clusternode = new_ptr(Node(bornme, 0, MH()));
  
  clusternode->deepHead(clusternode);
  MH()->firstNodeMap(bornme,clusternode);
  bornme->merger(MH());

  
  vector<NodePtr> current = {{clusternode}};
  vector<NodePtr> children;

  unsigned int k = 1;
  while ( ! thePureMEsMap[i - k].empty() ) {
    for ( auto tmp : current ){//j
      tmp->birth(thePureMEsMap[i - k]);
      for ( auto tmpchild : tmp->children() ) {//m
        children.push_back(tmpchild);
      }
    }
    current = children;
    children.clear();
    ++k;
  }
  
  if ( MH()->N() > i )
	  bornme->needsCorrelations();
  else 
  	bornme->needsNoCorrelations();
  
  bornme->cloneDependencies();
  MEs().push_back(bornme);
  
}




void MergingFactory::pushV(MatchboxMEBasePtr born, int i) {
  
  MatchboxMEBasePtr nlo = born->cloneMe();
  
  string pname = fullName() + "/" + nlo->name() + ".Virtual";
  if ( !(generator()->preinitRegister(nlo, pname)) ) 
  	throw InitException() 
  		<< "Virtual ME "<< pname << " already existing.";
    ////////////////////////////////////NLO///////////////////////////
  nlo->virtuals().clear();
  
  getVirtuals(nlo , false);
  
  if ( nlo->virtuals().empty() )
        throw InitException() 
      		<< "No insertion operators have been found for "
              	<< born->name() << ".\n";
  
  nlo->doOneLoopNoBorn();
    ////////////////////////////////////NLO///////////////////////////
  NodePtr clusternode = new_ptr(Node(nlo, 0,MH()));
  
  clusternode->deepHead(clusternode);
  clusternode->virtualContribution(true);
  MH()->firstNodeMap(nlo,clusternode);
  nlo->merger(MH());
  
  vector<NodePtr> current = {{clusternode}};
  vector<NodePtr> children;
  unsigned int k = 1;
  while ( ! thePureMEsMap[i - k].empty() ) {
    for ( auto tmp : current ){
      tmp->birth(thePureMEsMap[i - k]);
      for ( auto tmpchild : tmp->children())
      children.push_back(tmpchild);
    }
    current = children;
    children.clear();
    ++k;
  }
  
  if ( nlo->isOLPLoop() ) {
    int id = orderOLPProcess(nlo->subProcess(),
                             born->matchboxAmplitude(),
                             ProcessType::oneLoopInterference);
    nlo->olpProcess(ProcessType::oneLoopInterference,id);
    if ( !nlo->onlyOneLoop() && nlo->needsOLPCorrelators() ) {
      id = orderOLPProcess(nlo->subProcess(),
                           born->matchboxAmplitude(),
                           ProcessType::colourCorrelatedME2);
      nlo->olpProcess(ProcessType::colourCorrelatedME2,id);
    }
  }
  
  nlo->needsCorrelations();
  nlo->cloneDependencies();
  MEs().push_back(nlo);
  
}

void MergingFactory::pushR(MatchboxMEBasePtr born, int i) {
  MatchboxMEBasePtr bornme = born->cloneMe();
  
  string pname = fullName() + "/" + bornme->name() + ".Real";
  if ( !(generator()->preinitRegister(bornme, pname)) ) 
  	throw InitException() 
  	<< "Subtracted ME " << pname << " already existing.";
  
  
  NodePtr clusternode = new_ptr(Node(bornme, 1, MH()));
  clusternode->deepHead(clusternode);
  clusternode->subtractedReal(true);
  MH()->firstNodeMap(bornme,clusternode);
  bornme->merger(MH());
  
  
  
  vector<NodePtr> current = {{clusternode}};
  vector<NodePtr> children;
 
  unsigned int k = 1;
  
  while ( ! thePureMEsMap[i - k].empty() ) {
    for ( auto tmp : current ){
      tmp->birth(thePureMEsMap[i - k]);
      for ( auto tmpchild : tmp->children()) 
      	children.push_back(tmpchild);
    }
    current = children;
    children.clear();
    ++k;
  }
  if(clusternode->children().empty()){
	// This is a finite real contribution. 
	// This process is included in the LO merging.
 	return;
  }
  
  if ( MH()->N() > i )   bornme->needsCorrelations();
  else                   bornme->needsNoCorrelations();
  
  bornme->cloneDependencies(pname);
  
  MEs().push_back(bornme);
}

// MergingFactory should never order OLPs here,
// they're done elsewhere. 
void MergingFactory::orderOLPs() {}

#include "ThePEG/Utilities/StringUtils.h"
vector<string> MergingFactory::parseProcess(string in) {
  vector<string> process = StringUtils::split(in);
  if ( process.size() < 3 )
    throw Exception() 
	<< "MatchboxFactory: Invalid process."<< Exception::runerror;

  for ( string & p : process) {
    p = StringUtils::stripws(p);
  }
  theN = 0;
  bool prodprocess = true;
  vector<string> result;
  for ( const string & p : process ) {
    if ( p == "->" )
    	continue;
    
    if (p=="[") {
      prodprocess = false;
    } else if (p=="]") {
      prodprocess = false;
      // TODO what if there's stuff after the bracket?
      assert( p == process.back() );
      break;
    } else if (p=="[j") {
      prodprocess = false;
      ++theN;
    } else if (p=="j" && !prodprocess) {
      ++theN;
      prodprocess = false;
    } else if (p=="j]") {
      ++theN;
      prodprocess = false;
      // TODO what if there's stuff after the bracket?
      assert( p == process.back() );
      break;
    } else if ( prodprocess ) {
      result.push_back(p);
    } else {
      throw InitException()
      	<< "Unknown particle class \"" << p << '"' 
      	<< " in the process definition merging bracket.\n"
      	<< "Only jets (\"j\") are supported at the moment.";
    }
    
  }
  return result;
}




#include "Herwig/Utilities/Progress.h"
void MergingFactory::setup() {
  
  useMe();
  
  DipoleShowerHandlerPtr dsh=dynamic_ptr_cast<DipoleShowerHandlerPtr>(this->CKKWHandler());
  if(! dsh  )throw InitException() << "The showerhandler for the MergingFactory must be the DipoleShower. ";
  
  dsh->setMerger(MH());
  MH()->setFactory(this);
  MH()->setDipoleShower(dsh); 

 
  if(!ransetup){
    
    generator()->log() <<"\nStarting merging setup.\n\n";

    olpProcesses().clear();
    externalAmplitudes().clear();

    // We set the couplings in the ME to be fixed
    // and reweight in the history weight for this.
    setFixedCouplings(true);
    setFixedQEDCouplings(true);
    
    // rebind the particle data objects, can't use rebind() function
    for ( auto & g : particleGroups()) {
      for ( auto & p : g.second) {
        p = getParticleData(p->id());
      }
    }
    
    const PDVector& partons = particleGroups()["j"];
    unsigned int nl = 0;
    for ( const auto p : partons ) {
    	const Energy mass = p->hardProcessMass();
    	const long pid = p->id();

      if ( abs(pid) < 7 && mass == ZERO )
        ++nl;
      if ( pid > 0 && pid < 7 && mass == ZERO )
        nLightJetVec( pid );
      if ( pid > 0 && pid < 7 && mass != ZERO )
        nHeavyJetVec( pid );
    }
    nLight(nl/2);
    
    const PDVector& partonsInP = particleGroups()["p"];
    for ( const auto pip : partonsInP )
      if ( pip->id() > 0 && pip->id() < 7 && pip->hardProcessMass() == ZERO )
        nLightProtonVec( pip->id() );
    
      // fill the amplitudes
    if ( !amplitudes().empty() )  fillMEsMap();
    
      // Use the colour basis of the first element of amplitudes
      // to set the large N colour basis for the MergingHelper
    assert(!amplitudes().empty() );
    if ( !amplitudes()[0]->colourBasis() )
        throw Exception() << "MergingFactory::setup(): Expecting a colour basis object."
        << Exception::runerror;
    auto largeNBasis =
    amplitudes()[0]->colourBasis()->cloneMe();
    largeNBasis->clear();
    largeNBasis->doLargeN();
    MH()->largeNBasis(largeNBasis);
    
      // prepare the Born and virtual matrix elements
    for ( int i = 0 ; i <= max(0, MH()->N()) ; ++i ) prepare_BV(i);
    
      // prepare the real emission matrix elements
    for ( int i = 0 ; i <= MH()->N() ; ++i )  prepare_R(i);
    
    
    if (MH()->N()<=MH()->M()) {
      throw InitException() << "Merging: The number of NLOs need to be"
                            << "\nsmaller than the number of LO processes.\n";
    }
    
      // Order the external Amplitudes.
    orderOLPs();
    
      // start creating matrix elements
    MEs().clear();
    
      // count the subprocesses
    size_t numb = 0;
    size_t numv = 0;
    size_t numr = 0;
    
    for (int i = 0; i <= max(0, MH()->N()) ; ++i ) {
      for ( auto born : thePureMEsMap[i] ) {
        if (bornContributions()
	   ) {
              numb++;
        }
      }
    }
    for (int i = 0 ; i <=max(0, MH()->N()); ++i )
      for ( auto virt : thePureMEsMap[i] )
        if ( virtualContributions() && i <= MH()->M()) {
            numv++;
        }
    
    for (int i = 1; i <= max(0, MH()->N()) ; ++i )
      for ( auto real : thePureMEsMap[i] )
        if (realContributions() && i <= MH()->M() + 1   ){
            numr++;
    }
    
    
    
    
    
    if(int(numb+numv+numr) < theChunk){
      throw InitException() << "You try to chunk (Chunk="<<theChunk<<") the number of "<<
      "subprocesses ("<<int(numb+numv+numr)<<") into more parts than there are subprocesses.";
    }
    if(theChunk<theChunkPart)
      throw InitException() <<" The ChunkPart is larger than the Chunk.";
    
    if(theChunk > 0 && theChunkPart == 0 )
      throw InitException() <<" Set the ChunkPart ( = "<<theChunkPart<<" ) to i in [1,..,N] with N=Chunk.";
    
    if ( theChunkPart != 0 )
      generator()->log() << ANSI::yellow
      << "\n\nWarning: \nYou split up the runs into theChunks. This is no standard feature."
      << "\nYou are now responsible to make sure to run all theChunk parts."
      << "\nThis setup run is for: " << theChunkPart << "/" << theChunk<<"\n\n ";
    
    
    
    generator()->log() << ANSI::red << "Preparing Merging: ";
    generator()->log() << ANSI::green << numb << " x Born " << ANSI::red;
    if (MH()->M()>-1) {
      generator()->log() << ANSI::yellow << numv << " x Virtual ";
      generator()->log() << ANSI::blue << numr << " x Real " << ANSI::red << flush;
    }
    
    int countchunk=0;
    
    progress_display progressBar{ numb+numv+numr, generator()->log() };
    
    // insert the born contributions to ME vector
    for (int i = 0; i <= max(0, MH()->N()) ; ++i )
      for ( auto born : thePureMEsMap[i])
        if (bornContributions() ){
          countchunk++; // theChunkPart is in [1,...,theChunk]
          if ( theChunk == 0 || theChunkPart == countchunk ) pushB( born , i );
          if ( countchunk == theChunk ) countchunk=0;
          generator()->log() << ANSI::green;
          ++progressBar;
          generator()->log() << ANSI::reset;
        }
    // insert the virtual contributions to ME vector
    for (int i = 0 ; i <=max(0, MH()->N()); ++i )
      for ( auto virt : thePureMEsMap[i])
        if ( virtualContributions() && i <= MH()->M()){
          countchunk++;// theChunkPart is in [1,...,theChunk]
          if ( theChunk == 0 || theChunkPart == countchunk ) pushV(virt, i);
          if ( countchunk == theChunk ) countchunk=0;
          generator()->log() << ANSI::yellow;
          ++progressBar;
        }
    // insert the real contributions to ME vector
    for (int i = 1; i <= max(0, MH()->N()) ; ++i )
      for ( auto real : thePureMEsMap[i] )
        if (realContributions()&& i <= MH()->M() + 1 ){
          countchunk++;// theChunkPart is in [1,...,theChunk]
          if ( theChunk == 0 || theChunkPart == countchunk ) pushR(real, i);
          if ( countchunk == theChunk ) countchunk=0;
          generator()->log() << ANSI::blue;
          ++progressBar;
        }
    
    generator()->log() << ANSI::reset;
    
    if ( !externalAmplitudes().empty() ) {
      generator()->log() << "Initializing external amplitudes." << endl;
      progress_display progressBar{ externalAmplitudes().size(), generator()->log() };
      for ( const auto ext : externalAmplitudes() ) {
        if ( ! ext->initializeExternal() ) {
          throw InitException()
  << "error: failed to initialize amplitude '" << ext->name() << "'\n";
        }
        ++progressBar;
      }
      generator()->log()
      << "---------------------------------------------------" << endl;
    }
    
    
    if ( !olpProcesses().empty() ) {
      generator()->log() << "Initializing one-loop provider(s)." << endl;
      map<Ptr<MatchboxAmplitude>::tptr, map<pair<Process, int>, int> > olps;
      for (const auto oit : olpProcesses()) {
        olps[oit.first] = oit.second;
      }
      progress_display progressBar{ olps.size(), generator()->log() };
      for ( const auto olpit : olps ) {
        if ( !olpit.first->startOLP(olpit.second) ) {
          throw InitException() << "error: failed to start OLP for amplitude '" << olpit.first->name() << "'\n";
        }
        ++progressBar;
      }
      generator()->log()
      << "---------------------------------------------------\n" << flush;
    }
    
    
    generator()->log() <<"\nGenerated "<<MEs().size()<<" Subprocesses.\n"<<flush;
    generator()->log()
    << "---------------------------------------------------\n" << flush;
    
    generator()->log() <<"\n\n" << ANSI::red
                        <<"Note: Due to the unitarization of the higher  "
                       <<"\nmultiplicities, the individual cross sections "
                       <<"\ngiven in the integration and run step are not"
                       <<"\nmeaningful without merging." << ANSI::reset << endl;
    ransetup=true;
    
  }
  
  
}

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

void MergingFactory::persistentOutput(PersistentOStream & os) const {
  
  
  os
  << theonlymulti
  << ransetup
  << processMap << theMergingHelper <<theM<<theN<<theNonQCDCuts<<theChunk<<theChunkPart;
}

void MergingFactory::persistentInput(PersistentIStream & is, int) {
  is
  >> theonlymulti
  >> ransetup
  >> processMap >> theMergingHelper >>theM>>theN>>theNonQCDCuts>>theChunk>>theChunkPart;
}


#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"



void MergingFactory::Init() {
  
  static Parameter<MergingFactory, int> interfaceonlymulti("onlymulti", 
	 "Calculate only the ME with k additional partons.", 
         &MergingFactory::theonlymulti, -1, -1, 0,
         false, false, Interface::lowerlim);
  
    
  static Switch<MergingFactory, bool> interface_Unitarized("Unitarized", 
	 "Unitarize the cross section (default is unitarised. NLO merging must be unitarised).", 
	 &MergingFactory::unitarized, true, false, false);
  static SwitchOption interface_UnitarizedYes(interface_Unitarized, "Yes", 
	 "Switch on the unitarized cross section.", true);
  static SwitchOption interface_UnitarizedNo(interface_Unitarized, "No", 
	 "Switch off the unitarized cross section.", false);
  

  static Reference<MergingFactory,Merger> interfaceMergingHelper("MergingHelper",
   	 "Pointer to the Merging Helper.",
   	 &MergingFactory::theMergingHelper, false, false, true, true, false);
  
  
    static Parameter<MergingFactory, int> interfaceaddNLOLegs("NLOProcesses",
         "Set the number of virtual corrections to consider. 0 is default for no virtual correction.", 
	 &MergingFactory::theM, 0, 0, 0, false, false, Interface::lowerlim);
  
  static Reference<MergingFactory, Cuts > interfaceNonQcdCuts("NonQCDCuts",
        "Cut on non-QCD modified observables. Be carefull!",
        &MergingFactory::theNonQCDCuts, false, false, true, true, false);
  
  
  
  static Parameter<MergingFactory, int> interfacetheChunk("Chunk",
         "Cut the number of subprocesses into n theChunks.",
         &MergingFactory::theChunk, -1, -1, 0,
         false, false, Interface::lowerlim);
  
  static Parameter<MergingFactory, int> interfacetheChunkPart("ChunkPart",
         "If theChunk is larger then 0, set this parameter to the n'th part. Make sure to add the ChunksParts afterwards.",
         &MergingFactory::theChunkPart, -1, -1, 0,
         false, false, Interface::lowerlim);
 
}

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<MergingFactory, Herwig::MatchboxFactory> 
describeHerwigMergingFactory("Herwig::MergingFactory", "HwDipoleShower.so");
