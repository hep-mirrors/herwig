  // -*- C++ -*-
  //
  // MergeboxFactory.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2012 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the MergeboxFactory class.
  //


#include "MergingFactory.h"
#include "Node.h"

using namespace Herwig;
using std::ostream_iterator;

MergingFactory::MergingFactory() :MatchboxFactory(),
theonlyNLOParts(false),
theonlyvirtualNLOParts(false),
theonlyrealNLOParts(false),
theonlyUnlopsweights(false),
theunitarizeNLOParts(true),
calc_born(true),
calc_virtual(true),
calc_real(true),
unitarized(true),
NLOunitarized(true),
ransetup(false),
theonlymulti(-1),
theonlysub(-1),
divideSub(-1),
divideSubNumber(-1)
{}

MergingFactory::~MergingFactory(){}

IBPtr MergingFactory::clone() const {
  return new_ptr(*this);
}

IBPtr MergingFactory::fullclone() const {
  return new_ptr(*this);
}

void MergingFactory::doinit(){
  MatchboxFactory::doinit();
  if (subProcessGroups()) {
    throw InitException() << "There are no subprocess groups in merging!!!";
  }
  
  if(!(!(divideSub!=-1&&divideSubNumber==-1)||!(divideSub==-1&&divideSubNumber!=-1))){
    throw InitException() << "dividing the subprocesses is not performed correct.";
  }
  
  MH()->largeNBasis()->factory(this);

}


void MergingFactory::fill_amplitudes() {
  
  olpProcesses().clear();
  
  processMap[0] = getProcesses()[0];
  
  if(MH()->M()>=0) setHighestVirt(processMap[0].size()+MH()->M());
  
  MH()->N0(processMap[0].size());
  for ( int i = 1 ; i <= MH()->N() ; ++i ) {
    processMap[i] = processMap[i - 1];
    processMap[i].push_back("j");
  }
  
  for ( int i = 0 ; i <= MH()->N() ; ++i ) {
    vector<MatchboxMEBasePtr> ames = makeMEs(processMap[i], orderInAlphaS() + i,i<MH()->M()+1);
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
      for ( const auto p : born->subProcess().legs )
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
    if ( born->isOLPLoop() &&  i <= MH()->M()) {
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
  
    // check the additional insertion operators
  if ( !theVirtualsMap[i].empty() ) haveVirtuals = true;
  
  for ( auto & virt : theVirtualsMap[i] ) virt->factory(this);
  
    // check for consistent conventions on virtuals, if we are to include MH()->M()
  
  if(!(i > MH()->M()||haveVirtuals))throw InitException() << MH()->M()<<" NLO corrections requested,\n but no virtual contributions are found.";
  
  for ( auto & virt : DipoleRepository::insertionIOperators(dipoleSet()) )
  virt->factory(this);
  
  for ( auto & virt : DipoleRepository::insertionPKOperators(dipoleSet()) )
  virt->factory(this);
}

void MergingFactory::prepare_R(int i) {
  for ( auto real : pureMEsMap()[i])
  prepareME(real);
}


#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
void MergingFactory::getVirtuals(MatchboxMEBasePtr nlo,int i){
  
  for ( auto OV : theVirtualsMap[i] )
  if ( OV->apply(nlo->diagrams().front()->partons()) ){
    nlo->virtuals().push_back(OV);
  }
  
  for ( auto I : DipoleRepository::insertionIOperators(dipoleSet()) )
  if ( I->apply(nlo->diagrams().front()->partons()) ){
    nlo->virtuals().push_back(I);
  }
  
  for ( auto PK : DipoleRepository::insertionPKOperators(dipoleSet()) )
  if ( PK->apply(nlo->diagrams().front()->partons()) ){
    nlo->virtuals().push_back(PK);
  }
  
}



void MergingFactory::pushB(MatchboxMEBasePtr born, int i) {
  MatchboxMEBasePtr bornme = born->cloneMe();
  bornme->maxMultCKKW(1);
  bornme->minMultCKKW(0);
  
  
  string pname = fullName() + "/" + bornme->name() + ".Born";
  if ( !(generator()->preinitRegister(bornme, pname)) ) throw InitException() << "Born ME " << pname << " already existing.";
  
  
  getVirtuals(bornme,i);
  
  
  NodePtr clusternode = new_ptr(Node(bornme, 0,MH()));
  
  clusternode->deepHead(clusternode);
  MH()->firstNodeMap(bornme,clusternode);
  bornme->factory(this);
  bornme->merger(MH());
  
  
  
  
  vector<NodePtr> temp;
  vector<NodePtr> temp1;
  temp.push_back(clusternode);
  unsigned int k = 1;
  while (thePureMEsMap[i - k].size() != 0) {
    for ( auto tmp : temp){//j
      tmp->birth(thePureMEsMap[i - k]);
      for ( auto tmpchild : tmp->children()){//m
        if ( i <= MH()->M() ) {
          
          for ( auto I : DipoleRepository::insertionIOperators(dipoleSet()) )
          if ( I->apply(tmpchild->nodeME()->diagrams().front()->partons()) ){
            Ptr<MatchboxInsertionOperator>::ptr myI = I->cloneMe();
            tmpchild->nodeME()->virtuals().push_back(myI);
          }
          for ( auto I : DipoleRepository::insertionPKOperators(dipoleSet()) )
          if ( I->apply(tmpchild->nodeME()->diagrams().front()->partons()) ){
            Ptr<MatchboxInsertionOperator>::ptr myI = I->cloneMe();
            tmpchild->nodeME()->virtuals().push_back(myI);
          }
          
          tmpchild->nodeME()->noOneLoop();
        }
        temp1.push_back(tmpchild);
      }
    }
    temp = temp1;
    temp1.clear();
    k++;
  }
  
  if(MH()->N()>i)
  bornme->needsCorrelations();
  else bornme->needsNoCorrelations();
  
  bornme->cloneDependencies();
  MEs().push_back(bornme);
  
}




void MergingFactory::pushV(MatchboxMEBasePtr born, int i) {
  
  MatchboxMEBasePtr nlo = born->cloneMe();
  
  string pname = fullName() + "/" + nlo->name() + ".Virtual";
  if ( !(generator()->preinitRegister(nlo, pname)) ) throw InitException() << "Virtual ME " << pname << " already existing.";
    ////////////////////////////////////NLO///////////////////////////
  nlo->virtuals().clear();
  
  getVirtuals(nlo , i);
  
  if ( nlo->virtuals().empty() )
        throw InitException() << "No insertion operators have been found for "
              << born->name() << "\n";
  
  nlo->doOneLoopNoBorn();
    ////////////////////////////////////NLO///////////////////////////
  NodePtr clusternode = new_ptr(Node(nlo, 0,MH()));
  
  clusternode->deepHead(clusternode);
  clusternode->virtualContribution(true);
  MH()->firstNodeMap(nlo,clusternode);
  nlo->merger(MH());
  
  vector<NodePtr> temp;
  vector<NodePtr> temp1;
  temp.push_back(clusternode);
  unsigned int k = 1;
  while (thePureMEsMap[i - k].size() != 0) {
    for ( auto tmp : temp ){
      tmp->birth(thePureMEsMap[i - k]);
      for ( auto tmpchild : tmp->children())
      temp1.push_back(tmpchild);
    }
    temp = temp1;
    k++;
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

void MergingFactory::pushProR(MatchboxMEBasePtr born, int i) {
  MatchboxMEBasePtr bornme = born->cloneMe();
  
  string pname = fullName() + "/" + bornme->name() + ".Real";
  if ( !(generator()->preinitRegister(bornme, pname)) ) throw InitException() << "Subtracted ME " << pname << " already existing.";
  
  
  NodePtr clusternode = new_ptr(Node(bornme, 1,MH()));
  clusternode->deepHead(clusternode);
  clusternode->subtractedReal(true);
  MH()->firstNodeMap(bornme,clusternode);
  bornme->merger(MH());
  
  
  
  vector<NodePtr> temp;
  vector<NodePtr> temp1;
  temp.push_back(clusternode);
  
  unsigned int k = 1;
  
  while (thePureMEsMap[i - k].size() != 0) {
    for ( auto tmp : temp ){
      tmp->birth(thePureMEsMap[i - k]);
      for ( auto tmpchild : tmp->children()) temp1.push_back(tmpchild);
    }
    temp = temp1;
    k++;
  }
  
  
  
  if(MH()->N()>i)
  bornme->needsCorrelations();
  else bornme->needsNoCorrelations();
  
  bornme->cloneDependencies(pname);
  
  MEs().push_back(bornme);
}

void MergingFactory::orderOLPs() {
  
}

#include "ThePEG/Utilities/StringUtils.h"
vector<string> MergingFactory::parseProcess(string in) {
  vector<string> process = StringUtils::split(in);
  if ( process.size() < 3 )
  throw Exception() << "MatchboxFactory: Invalid process."
  << Exception::runerror;
  for ( string & p : process) {
    p = StringUtils::stripws(p);
  }
  theN=0;
  bool prodprocess=true;
  vector<string> pprocess;
  for (string const p : process) {
    if ( p == "->" )
    continue;
    
    if (p=="[") {
      prodprocess=false;
    }else if (p=="]") {
      prodprocess=false;
      break;
    }else if (p=="[j") {
      prodprocess=false;
      theN++;
    }else if (p=="j"&&!prodprocess) {
      theN++;
      prodprocess=false;
    }else if (p=="j]") {
      theN++;
      prodprocess=false;
      break;
    }else if(prodprocess){
      pprocess.push_back(p);
    }else{
      throw InitException()
      <<p<<" in the brackets of the process definition is not recognized."
      <<   "\n Only j as jets are recognized.";
    }
    
  }
  return pprocess;
}




#include <boost/progress.hpp>
void MergingFactory::setup() {
  
  useMe();
  
  if(!ransetup){
    
    generator()->log() <<"\nStarting merging setup.\n\n";

    olpProcesses().clear();
    externalAmplitudes().clear();
    setFixedCouplings(true);
    setFixedQEDCouplings(true);
    
      // rebind the particle data objects
    for ( auto & g : particleGroups()) {
      for ( auto & p : g.second) {
        p = getParticleData(p->id());
      }
    }
    
    const PDVector& partons = particleGroups()["j"];
    unsigned int nl = 0;
    for ( const auto p : partons ) {
      if ( abs(p->id()) < 7 && p->hardProcessMass() == ZERO )
      ++nl;
      if ( p->id() > 0 && p->id() < 7 && p->hardProcessMass() == ZERO )
      nLightJetVec( p->id() );
      if ( p->id() > 0 && p->id() < 7 && p->hardProcessMass() != ZERO )
      nHeavyJetVec( p->id() );
    }
    nLight(nl/2);
    
    const PDVector& partonsInP = particleGroups()["p"];
    for ( const auto pip : partonsInP )
      if ( pip->id() > 0 && pip->id() < 7 && pip->hardProcessMass() == ZERO )
        nLightProtonVec( pip->id() );
    
    for ( auto & amp: amplitudes() ) amp->factory(this);
    
      //fill the amplitudes
    if ( !amplitudes().empty() )  fill_amplitudes();
    
      // prepare the Born and virtual matrix elements
    for ( int i = 0 ; i <= max(0, MH()->N()) ; ++i ) prepare_BV(i);
    
      // prepare the real emission matrix elements
    for ( int i = 0 ; i <= MH()->N() ; ++i )  prepare_R(i);
    
    
    orderOLPs();
    
      // start creating matrix elements
    MEs().clear();
    
    int onlysubcounter=0;
      // count the subprocesses
    size_t numb = 0;
    size_t numv = 0;
    size_t numr = 0;
    
    for (int i = 0; i <= max(0, MH()->N()) ; ++i ){
      for ( auto born : thePureMEsMap[i]){
        if (calc_born && !theonlyUnlopsweights){
          if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)
             ||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber)){
              numb++;
          }
        }
      }
    }
    for (int i = 0 ; i <=max(0, MH()->N()); ++i )
    for ( auto virt : thePureMEsMap[i])
    if ( calc_virtual && i <= MH()->M()){
      if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)
         ||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
      numv++;
    }
    
    for (int i = 0; i <= max(0, MH()->N()) ; ++i )
    for ( auto real : thePureMEsMap[i] )
    if (calc_real&& i <= MH()->M() + 1 && i > 0 && !theonlyvirtualNLOParts&&!theonlyUnlopsweights){
      if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)
         ||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
      numr++;
    }
    
    
    onlysubcounter=0;
    
    generator()->log() << "\e[31mPreparing Merging:" << flush;
    generator()->log() << "\e[32m "<<numb<<" x Born \e[31m" << flush;
    if(MH()->M()>-1){
      generator()->log() << "\e[93m"<<numv<<" x Virtual " << flush;
      generator()->log() << "\e[34m"<<numr<<" x Real \e[31m" << flush;
    }
    
    boost::progress_display * progressBar = new boost::progress_display(numb+numv+numr,generator()->log());
      for (int i = 0; i <= max(0, MH()->N()) ; ++i ){
        for ( auto born : thePureMEsMap[i]){
          if (calc_born && !theonlyUnlopsweights){
            if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)
               ||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber)){
              pushB(born, i);
              onlysubcounter++;
              generator()->log() << "\e[32m";
              ++(*progressBar);
              generator()->log() << "\e[0m";
            }
          }
        }
      }
    
    
  
   
    for (int i = 0 ; i <=max(0, MH()->N()); ++i )
    for ( auto virt : thePureMEsMap[i])
    if ( calc_virtual && i <= MH()->M()){
      if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)
         ||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
      pushV(virt, i);
      onlysubcounter++;
       generator()->log() << "\e[93m";
      ++(*progressBar);
    }
    
    for (int i = 0; i <= max(0, MH()->N()) ; ++i )
    for ( auto real : thePureMEsMap[i] )
    if (calc_real&& i <= MH()->M() + 1 && i > 0 && !theonlyvirtualNLOParts&&!theonlyUnlopsweights){
      if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)
         ||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
      pushProR(real, i);
      onlysubcounter++;
      generator()->log() << "\e[34m";
      ++(*progressBar);
    }
    generator()->log() << "\e[0m";
    delete progressBar;
    
    if ( !externalAmplitudes().empty() ) {
      generator()->log() << "Initializing external amplitudes.\n" << flush;
      boost::progress_display * progressBar =
      new boost::progress_display(externalAmplitudes().size(),generator()->log());
      for ( const auto ext : externalAmplitudes()) {
        if ( ! ext->initializeExternal() ) {
          throw InitException()
  << "error: failed to initialize amplitude '" << ext->name() << "'\n";
        }
        ++(*progressBar);
      }
      delete progressBar;
      generator()->log()
      << "---------------------------------------------------\n"
  << flush;
    }
    
    
    if ( !olpProcesses().empty() ) {
      generator()->log() << "Initializing one-loop provider(s).\n" << flush;
      map<Ptr<MatchboxAmplitude>::tptr, map<pair<Process, int>, int> > olps;
      for (const auto oit : olpProcesses()) {
        olps[oit.first] = oit.second;
      }
      boost::progress_display * progressBar = new boost::progress_display(olps.size(), generator()->log());
      for ( const auto olpit : olps ) {
        if ( !olpit.first->startOLP(olpit.second) ) {
          throw InitException() << "error: failed to start OLP for amplitude '" << olpit.first->name() << "'\n";
        }
        ++(*progressBar);
      }
      delete progressBar;
      generator()->log()
      << "---------------------------------------------------\n" << flush;
    }
    
    
    generator()->log() <<"\nGenerated "<<MEs().size()<<" Subprocesses.\n"<<flush;
    if(theonlysub!=-1)generator()->log() <<" ( "<<theonlysub<<"/"<<onlysubcounter<<" )"<<flush;
    generator()->log()
    << "---------------------------------------------------\n" << flush;
    
    generator()->log() <<"\n\n\e[31m"
                        <<"Note: Due to the unitarization of the higher  "
                       <<"\nmultiplicities, the individual cross sections "
                       <<"\ngiven in the integration and run step are not"
                       <<"\nmeaningful without merging.\e[0m\n"<<flush;
    ransetup=true;
    
  }
  
  
}

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

void MergingFactory::persistentOutput(PersistentOStream & os) const {
  
  
  os
  << calc_born  << calc_virtual << calc_real
  << theonlyUnlopsweights  << theonlymulti
  << divideSub  << divideSubNumber
  << theonlysub  << ransetup
  << processMap << theMergingHelper <<theM<<theN;
}

void MergingFactory::persistentInput(PersistentIStream & is, int) {
  is
  >> calc_born >> calc_virtual >> calc_real
  >> theonlyUnlopsweights >> theonlymulti
  >> divideSub >> divideSubNumber
  >> theonlysub >> ransetup
  >> processMap >> theMergingHelper >>theM>>theN;
}


#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"



void MergingFactory::Init() {
  
  
  
  static Parameter<MergingFactory, int> interfaceonlymulti("onlymulti", "calculate only the ME with k additional partons.", &MergingFactory::theonlymulti, -1, -1, 0,
                                                           false, false, Interface::lowerlim);
  
  
  static Parameter<MergingFactory, int> interfaceonlysub("onlysub", "calculate only one subProcess. this is for building grids.", &MergingFactory::theonlysub, -1, -1, 0,
                                                         false, false, Interface::lowerlim);
  
  
  
  
  
  
  
  
  
  static Parameter<MergingFactory, int> interfacedivideSub("divideSub", "calculate only one subProcess. this is for building grids.", &MergingFactory::divideSub, -1, -1, 0,
                                                           false, false, Interface::lowerlim);
  
  
  static Parameter<MergingFactory, int> interfacedivideSubNumber("divideSubNumber", "calculate only one subProcess. this is for building grids.", &MergingFactory::divideSubNumber, -1, -1, 0,
                                                                 false, false, Interface::lowerlim);
  
  
  
  
  
  static Switch<MergingFactory, bool> interface_calc_born("calc_born", "[debug] Switch on or off the born contribution.", &MergingFactory::calc_born, true,
                                                          false, false);
  static SwitchOption interface_calc_bornOn(interface_calc_born, "On", "Switch on calculation of born.", true);
  static SwitchOption interface_calc_bornOff(interface_calc_born, "Off", "Switch off calculation of born.", false);
  
  static Switch<MergingFactory, bool> interface_calc_virtual("calc_virtual", "[debug] Switch on or off the virtual contribution.",
                                                             &MergingFactory::calc_virtual, true, false, false);
  static SwitchOption interface_calc_virtualOn(interface_calc_virtual, "On", "Switch on calculation of virtual.", true);
  static SwitchOption interface_calc_virtualOff(interface_calc_virtual, "Off", "Switch off calculation of virtual.", false);
  
  static Switch<MergingFactory, bool> interface_calc_real("calc_real", "[debug] Switch on or off the real contribution.", &MergingFactory::calc_real, true,
                                                          false, false);
  static SwitchOption interface_calc_realOn(interface_calc_real, "On", "Switch on calculation of real.", true);
  static SwitchOption interface_calc_realOff(interface_calc_real, "Off", "Switch off calculation of real.", false);
  
  
  
  
  static Switch<MergingFactory, bool> interface_theonlyNLOParts("onlyNLOParts", "Switch on or off the onlyNLOParts.", &MergingFactory::theonlyNLOParts, true, false,
                                                                false);
  static SwitchOption interface_theonlyNLOPartsOn(interface_theonlyNLOParts, "On", "Switch on the theonlyNLOParts.", true);
  static SwitchOption interface_theonlyNLOPartsOff(interface_theonlyNLOParts, "Off", "Switch off the theonlyNLOParts.", false);
  
  static Switch<MergingFactory, bool> interface_theonlyvirtualNLOParts("onlyvirtualNLOParts", "Switch on or off the onlyvirtualNLOParts.", &MergingFactory::theonlyvirtualNLOParts, true, false,
                                                                       false);
  static SwitchOption interface_theonlyvirtualNLOPartsOn(interface_theonlyvirtualNLOParts, "On", "Switch on the theonlyvirtualNLOParts.", true);
  static SwitchOption interface_theonlyvirtualNLOPartsOff(interface_theonlyvirtualNLOParts, "Off", "Switch off the theonlyvirtualNLOParts.", false);
  
  static Switch<MergingFactory, bool> interface_theonlyrealNLOParts("onlyrealNLOParts", "Switch on or off the onlyrealNLOParts.", &MergingFactory::theonlyrealNLOParts, true, false,
                                                                    false);
  static SwitchOption interface_theonlyrealNLOPartsOn(interface_theonlyrealNLOParts, "On", "Switch on the theonlyrealNLOParts.", true);
  static SwitchOption interface_theonlyrealNLOPartsOff(interface_theonlyrealNLOParts, "Off", "Switch off the theonlyrealNLOParts.", false);
  
  static Switch<MergingFactory, bool> interface_theunitarizeNLOParts("unitarizeNLOParts", "Switch on or off the unitarizeNLOParts.", &MergingFactory::theunitarizeNLOParts, true, false,
                                                                     false);
  static SwitchOption interface_theunitarizeNLOPartsOn(interface_theunitarizeNLOParts, "On", "Switch on the unitarizeNLOParts.", true);
  static SwitchOption interface_theunitarizeNLOPartsOff(interface_theunitarizeNLOParts, "Off", "Switch off the unitarizeNLOParts.", false);
  
  
  static Switch<MergingFactory, bool> interface_theonlyUnlopsweights("onlyUnlopsweights", "Switch on or off the onlyUnlopsweights.", &MergingFactory::theonlyUnlopsweights, true, false,
                                                                     false);
  static SwitchOption interface_theonlyUnlopsweightsOn(interface_theonlyUnlopsweights, "On", "Switch on the onlyUnlopsweights.", true);
  static SwitchOption interface_theonlyUnlopsweightsOff(interface_theonlyUnlopsweights, "Off", "Switch off the onlyUnlopsweights.", false);
  
  
  
  static Switch<MergingFactory, bool> interface_Unitarized("Unitarized", "Unitarize the cross section.", &MergingFactory::unitarized, true, false, false);
  static SwitchOption interface_UnitarizedOn(interface_Unitarized, "On", "Switch on the unitarized cross section.", true);
  static SwitchOption interface_UnitarizedOff(interface_Unitarized, "Off", "Switch off the unitarized cross section.", false);
  
  
  
  
  static Switch<MergingFactory, bool> interface_NLOUnitarized("NLOUnitarized", "Unitarize the cross section.", &MergingFactory::NLOunitarized, true, false, false);
  static SwitchOption interface_NLOUnitarizedOn(interface_NLOUnitarized, "On", "Switch on the unitarized NLO cross section.", true);
  static SwitchOption interface_NLOUnitarizedOff(interface_NLOUnitarized, "Off", "Switch off the unitarized NLO cross section.", false);
  
  static Reference<MergingFactory,Merger> interfaceMergingHelper
  ("MergingHelper",
   "",
   &MergingFactory::theMergingHelper, false, false, true, true, false);
  
  
  
  static Parameter<MergingFactory, int> interfaceaddNLOLegs("NLOProcesses",
                                                            "Set the number of virtual corrections to consider. 0 is default for no virtual correction.", &MergingFactory::theM, 0, 0, 0, false, false,
                                                            Interface::lowerlim);
  
  
  
}

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<MergingFactory, Herwig::MatchboxFactory> describeHerwigMergingFactory("Herwig::MergingFactory", "HwDipoleShower.so");
