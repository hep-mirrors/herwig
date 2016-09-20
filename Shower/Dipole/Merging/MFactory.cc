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
#include "MFactory.h"
#include "Node.h"

#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Utility/DiagramDrawer.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

#include "ThePEG/MatrixElement/MEBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SU2Helper.h"

#include "ThePEG/Cuts/JetFinder.h"
#include "fastjet/ClusterSequence.hh"

#include <boost/progress.hpp>

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;
using std::ostream_iterator;

MFactory::MFactory() :MatchboxFactory(),
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
 theonlyk(-1),
 theonlymulti(-1),
 theonlyabove(-1),
 theonlysub(-1),
 divideSub(-1),
 divideSubNumber(-1)
{}

MFactory::~MFactory(){}

IBPtr MFactory::clone() const {
  return new_ptr(*this);
}

IBPtr MFactory::fullclone() const {
  return new_ptr(*this);
}




void MFactory::fill_amplitudes() {
  
  olpProcesses().clear();
  
  
  processMap[0] = getProcesses()[0];
  if(MH()->M()>=0)
    setHighestVirt(processMap[0].size()+MH()->M());
 
  
  
  MH()->N0(processMap[0].size());
  for ( int i = 1 ; i <= MH()->N() ; ++i ) {
    processMap[i] = processMap[i - 1];
    processMap[i].push_back("j");
  }
  
  for ( int i = 0 ; i <= MH()->N() ; ++i ) {
    vector<Ptr<MatchboxMEBase>::ptr> ames = makeMEs(processMap[i], orderInAlphaS() + i,i<MH()->M()+1);
    copy(ames.begin(), ames.end(), back_inserter(pureMEsMap()[i]));
  }
}

void MFactory::prepare_BV(int i) {
    // check if we have virtual contributions
  bool haveVirtuals = true;
  for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born = pureMEsMap()[i].begin() ; born != pureMEsMap()[i].end() ; ++born ) {
    
    prepareME(*born);
    
    
    if ( (*born)->isOLPTree() ) {
      int id = orderOLPProcess((*born)->subProcess(),
                               (*born)->matchboxAmplitude(),
                               ProcessType::treeME2);
      (*born)->olpProcess(ProcessType::treeME2,id);
      id = orderOLPProcess((*born)->subProcess(),
                           (*born)->matchboxAmplitude(),
                           ProcessType::colourCorrelatedME2);
      (*born)->olpProcess(ProcessType::colourCorrelatedME2,id);
      
      bool haveGluon = false;
      for ( PDVector::const_iterator p = (*born)->subProcess().legs.begin();
           p != (*born)->subProcess().legs.end(); ++p )
        if ( (**p).id() == 21 ) {
          haveGluon = true;
          break;
        }
      if ( haveGluon ) {
        id = orderOLPProcess((*born)->subProcess(),
                             (*born)->matchboxAmplitude(),
                             ProcessType::spinColourCorrelatedME2);
        (*born)->olpProcess(ProcessType::spinColourCorrelatedME2,id);
      }
    }
    if ( (*born)->isOLPLoop() &&  i <= MH()->M()) {
      int id = orderOLPProcess((*born)->subProcess(),
                               (*born)->matchboxAmplitude(),
                               ProcessType::oneLoopInterference);
      (*born)->olpProcess(ProcessType::oneLoopInterference,id);
      if ( !(*born)->onlyOneLoop() && (*born)->needsOLPCorrelators() ) {
        id = orderOLPProcess((*born)->subProcess(),
                             (*born)->matchboxAmplitude(),
                             ProcessType::colourCorrelatedME2);
        (*born)->olpProcess(ProcessType::colourCorrelatedME2,id);
      }
    }
    
    
    
    
    
    
    haveVirtuals &= (**born).haveOneLoop();
  }
  
    // check the additional insertion operators
  if ( !theVirtualsMap[i].empty() ) haveVirtuals = true;
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = theVirtualsMap[i].begin() ; virt != theVirtualsMap[i].end() ; ++virt ) {
    (**virt).factory(this);
  }
  
  
    // check for consistent conventions on virtuals, if we are to include MH()->M()
  if ( i <= MH()->M() ) {
    if ( !haveVirtuals ) {
      throw InitException() << "Could not find amplitudes for all virtual contributions needed.\n";
    }
  }
  
    // prepare dipole insertion operators
    // Need MH()->M() for alpha-improvement
    //if ( i <= MH()->M() ) {
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
       = DipoleRepository::insertionIOperators(dipoleSet()).begin();
       virt != DipoleRepository::insertionIOperators(dipoleSet()).end(); ++virt ) {
    (**virt).factory(this);
  }
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
       = DipoleRepository::insertionPKOperators(dipoleSet()).begin();
       virt != DipoleRepository::insertionPKOperators(dipoleSet()).end(); ++virt ) {
    
    (**virt).factory(this);
  }
  
    //}
  
}

void MFactory::prepare_R(int i) {
  for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator
       real = thePureMEsMap[i].begin() ;
       real!= thePureMEsMap[i].end() ; ++real )
    prepareME(*real);
}

void MFactory::pushB(Ptr<MatchboxMEBase>::ptr born, int i) {
  Ptr<MatchboxMEBase>::ptr bornme = (*born).cloneMe();
  bornme->maxMultCKKW(1);
  bornme->minMultCKKW(0);
  
  
  string pname = fullName() + "/" + bornme->name() + ".Born";
  if ( !(generator()->preinitRegister(bornme, pname)) ) throw InitException() << "Born ME " << pname << " already existing.";
  
  bornme->virtuals().clear();
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = theVirtualsMap[i].begin() ; virt != theVirtualsMap[i].end() ; ++virt ) {
    if ( (**virt).apply((*born).diagrams().front()->partons()) ) bornme->virtuals().push_back(*virt);
  }
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionIOperators(dipoleSet()).begin() ;
       virt != DipoleRepository::insertionIOperators(dipoleSet()).end() ; ++virt ) {
    if ( (**virt).apply((*born).diagrams().front()->partons()) ) bornme->virtuals().push_back(*virt);
  }
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionPKOperators(dipoleSet()).begin() ;
       virt != DipoleRepository::insertionPKOperators(dipoleSet()).end() ; ++virt ) {
    if ( (**virt).apply((*born).diagrams().front()->partons()) ) bornme->virtuals().push_back(*virt);
  }
  
  
  Ptr<Node>::ptr clusternode = new_ptr(Node(bornme, i, 0));
  MH()->firstNodeMap(bornme,clusternode);
  bornme->factory(this);
  bornme->merger(MH());
  clusternode->MH(MH());
  clusternode->deepHead(clusternode);
  
  
  
  vector<Ptr<Node>::ptr> temp;
  vector<Ptr<Node>::ptr> temp1;
  temp.push_back(clusternode);
  unsigned int k = 1;
  while (thePureMEsMap[i - k].size() != 0) {
    for ( unsigned int j = 0 ; j < temp.size() ; j++ ) {
      temp[j]->birth(thePureMEsMap[i - k]);
      for ( unsigned int m = 0 ; m < temp[j]->children().size() ; m++ ) {
        if ( i <= MH()->M() ) {
          for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionIOperators(dipoleSet()).begin() ;
               virt != DipoleRepository::insertionIOperators(dipoleSet()).end() ; ++virt ) {
            if ( (**virt).apply((*(temp[j]->children()[m]->nodeME())).diagrams().front()->partons()) ){
              Ptr<MatchboxInsertionOperator>::ptr myIOP = (**virt).cloneMe();
              temp[j]->children()[m]->nodeME()->virtuals().push_back(myIOP);
            }
          }
          for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionPKOperators(dipoleSet()).begin() ;
               virt != DipoleRepository::insertionPKOperators(dipoleSet()).end() ; ++virt ) {
            if ( (**virt).apply((*(temp[j]->children()[m]->nodeME())).diagrams().front()->partons()) ){
              Ptr<MatchboxInsertionOperator>::ptr myIOP = (**virt).cloneMe();
              temp[j]->children()[m]->nodeME()->virtuals().push_back(myIOP);
            }
          }
          temp[j]->children()[m]->nodeME()->noOneLoop();
        }
        temp1.push_back(temp[j]->children()[m]);
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
  if ( theonlyk == -1  || (theonlyk == i-1 &&unitarized) || theonlyk == i) {
    MEs().push_back(bornme);
  }
}





void MFactory::pushV(Ptr<MatchboxMEBase>::ptr born, int i) {
  
  Ptr<MatchboxMEBase>::ptr nlo = (*born).cloneMe();
  
  string pname = fullName() + "/" + nlo->name() + ".Virtual";
  if ( !(generator()->preinitRegister(nlo, pname)) ) throw InitException() << "Virtual ME " << pname << " already existing.";
    ////////////////////////////////////NLO///////////////////////////
  nlo->virtuals().clear();
  if ( !nlo->onlyOneLoop() ) {
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = theVirtualsMap[i].begin() ; virt != theVirtualsMap[i].end() ; ++virt ) {
      if ( (**virt).apply((*born).diagrams().front()->partons()) ) nlo->virtuals().push_back(*virt);
    }
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionIOperators(dipoleSet()).begin() ;
         virt != DipoleRepository::insertionIOperators(dipoleSet()).end() ; ++virt ) {
      if ( (**virt).apply((*born).diagrams().front()->partons()) ) nlo->virtuals().push_back(*virt);
    }
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionPKOperators(dipoleSet()).begin() ;
         virt != DipoleRepository::insertionPKOperators(dipoleSet()).end() ; ++virt ) {
      if ( (**virt).apply((*born).diagrams().front()->partons()) ) nlo->virtuals().push_back(*virt);
    }
    if ( nlo->virtuals().empty() ) throw InitException() << "No insertion operators have been found for " << (*born).name() << "\n";
  }
  nlo->doOneLoopNoBorn();
    ////////////////////////////////////NLO///////////////////////////
  Ptr<Node>::ptr clusternode = new_ptr(Node(nlo, i, 0));
  clusternode->virtualContribution(true);
  MH()->firstNodeMap(nlo,clusternode);
  nlo->merger(MH());
  clusternode->deepHead(clusternode);
  clusternode->MH(MH());
  
  vector<Ptr<Node>::ptr> temp;
  vector<Ptr<Node>::ptr> temp1;
  temp.push_back(clusternode);
  unsigned int k = 1;
  while (thePureMEsMap[i - k].size() != 0) {
    for ( unsigned int j = 0 ; j < temp.size() ; j++ ) {
      temp[j]->birth(thePureMEsMap[i - k]);
      for ( unsigned int m = 0 ; m < temp[j]->children().size() ; ++m ) {
        temp1.push_back(temp[j]->children()[m]);
      }
    }
    temp = temp1;
    k++;
  }
  
  if ( nlo->isOLPLoop() ) {
    int id = orderOLPProcess(nlo->subProcess(),
                             (*born).matchboxAmplitude(),
                             ProcessType::oneLoopInterference);
    nlo->olpProcess(ProcessType::oneLoopInterference,id);
    if ( !nlo->onlyOneLoop() && nlo->needsOLPCorrelators() ) {
      id = orderOLPProcess(nlo->subProcess(),
                           (*born).matchboxAmplitude(),
                           ProcessType::colourCorrelatedME2);
      nlo->olpProcess(ProcessType::colourCorrelatedME2,id);
    }
  }
  
  nlo->needsCorrelations();
  nlo->cloneDependencies();
  if ( theonlyk == -1 || theonlyk == i - 1|| theonlyk == i - 2  ) {
    MEs().push_back(nlo);
  }
}

void MFactory::pushProR(Ptr<MatchboxMEBase>::ptr born, int i) {
  Ptr<MatchboxMEBase>::ptr bornme = (*born).cloneMe();
  
  string pname = fullName() + "/" + bornme->name() + ".Real";
  if ( !(generator()->preinitRegister(bornme, pname)) ) throw InitException() << "Subtracted ME " << pname << " already existing.";
  
  
  Ptr<Node>::ptr clusternode = new_ptr(Node(bornme, i - 2, 1));
  clusternode->subtractedReal(true);
  MH()->firstNodeMap(bornme,clusternode);
  bornme->merger(MH());
  clusternode->deepHead(clusternode);
  clusternode->MH(MH());
  
  vector<Ptr<Node>::ptr> temp;
  vector<Ptr<Node>::ptr> temp1;
  temp.push_back(clusternode);
  
  unsigned int k = 1;
  while (thePureMEsMap[i - k].size() != 0) {
    for ( unsigned int j = 0 ; j < temp.size() ; j++ ) {
      temp[j]->birth(thePureMEsMap[i - k]);
      for ( unsigned int m = 0 ; m < temp[j]->children().size() ; ++m ) {
        temp1.push_back(temp[j]->children()[m]);
      }
    }
    temp = temp1;
    k++;
  }
  
  if(MH()->N()>i)
    bornme->needsCorrelations();
  else bornme->needsNoCorrelations();
  
  bornme->cloneDependencies(pname);
  
  if ( theonlyk == -1 || theonlyk == i - 1 || theonlyk == i - 2 ) {
    MEs().push_back(bornme);
  }
}

void MFactory::orderOLPs() {
  
}


vector<string> MFactory::parseProcess(string in) {
  vector<string> process = StringUtils::split(in);
  if ( process.size() < 3 )
    throw Exception() << "MatchboxFactory: Invalid process."
    << Exception::runerror;
  for ( vector<string>::iterator p = process.begin();
       p != process.end(); ++p ) {
    *p = StringUtils::stripws(*p);
  }
  theN=0;
  bool prodprocess=true;
  vector<string> pprocess;
  for ( vector<string>::const_iterator p = process.begin();
       p != process.end(); ++p ) {
    if ( *p == "->" )
      continue;

    if (*p=="[") {
      prodprocess=false;
    }else if (*p=="]") {
      prodprocess=false;
      break;
    }else if (*p=="[j") {
      prodprocess=false;
      theN++;
    }else if (*p=="j"&&!prodprocess) {
      theN++;
      prodprocess=false;
    }else if (*p=="j]") {
      theN++;
      prodprocess=false;
      break;
    }else if(prodprocess){
      pprocess.push_back(*p);
    }else{
      cout<<"\nWarning: "<<*p<<" in the brackets of the process definition is not recognized.\n Only j as jets are recognized.";
    }
    
  }
  return pprocess;
}





void MFactory::setup() {
  
  useMe();
  
  if(!ransetup){
    
    olpProcesses().clear();
    externalAmplitudes().clear();
    setFixedCouplings(true);
    setFixedQEDCouplings(true);
    
    
    
    const PDVector& partons = particleGroups()["j"];
    unsigned int nl = 0;
    
    
      // rebind the particle data objects
    for ( map<string, PDVector>::iterator g = particleGroups().begin() ; g != particleGroups().end() ; ++g ) {
      for ( PDVector::iterator p = g->second.begin() ; p != g->second.end() ; ++p ) {
        *p = getParticleData((**p).id());
      }
    }
    
    for ( PDVector::const_iterator p = partons.begin();
         p != partons.end(); ++p ) {
      if ( abs((**p).id()) < 7 && (**p).hardProcessMass() == ZERO )
        ++nl;
      if ( (**p).id() > 0 && (**p).id() < 7 && (**p).hardProcessMass() == ZERO )
        nLightJetVec( (**p).id() );
      if ( (**p).id() > 0 && (**p).id() < 7 && (**p).hardProcessMass() != ZERO )
        nHeavyJetVec( (**p).id() );
    }
    nLight(nl/2);
    
    const PDVector& partonsInP = particleGroups()["p"];
    for ( PDVector::const_iterator pip = partonsInP.begin();
         pip != partonsInP.end(); ++pip ) {
      if ( (**pip).id() > 0 && (**pip).id() < 7 && (**pip).hardProcessMass() == ZERO )
        nLightProtonVec( (**pip).id() );
    }
    
    
    for ( vector<Ptr<MatchboxAmplitude>::ptr>::iterator amp
         = amplitudes().begin(); amp != amplitudes().end(); ++amp )
      (**amp).factory(this);
    
    MH()->largeNBasis()->factory(this);
    
    assert(!(divideSub!=-1&&divideSubNumber==-1)||!(divideSub==-1&&divideSubNumber!=-1));
    assert(!subProcessGroups());
    
      //fill the amplitudes
    if ( !amplitudes().empty() ) {
      fill_amplitudes();
    }
    
      // prepare the Born and virtual matrix elements
    for ( int i = 0 ; i <= max(0, MH()->N()) ; ++i ) prepare_BV(i);
    
      // prepare the real emission matrix elements
    for ( int i = 0 ; i <= MH()->N() ; ++i )  prepare_R(i);
    
    
    orderOLPs();
    
      // start creating matrix elements
    MEs().clear();
    
    int onlysubcounter=0;
    int i = theonlyabove ;
    
    
    if (calc_born) {
      for (; i <= max(0, MH()->N()) ; ++i ) {
        for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born = thePureMEsMap[i].begin() ; born != thePureMEsMap[i].end() ; ++born ) {
          if (!theonlyUnlopsweights&&( theonlyk == -1  || theonlyk == i -1|| theonlyk == i )){
            if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
              pushB(*born, i);
            onlysubcounter++;
          }
        }
      }
    }
    
    
    
    if (calc_virtual) {
      i = theonlyabove ;
      for (; i <=max(0, MH()->N()); ++i ) {
        for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born = thePureMEsMap[i].begin() ; born != thePureMEsMap[i].end() ; ++born ) {
          if ( i <= MH()->M()  &&( theonlyk == -1  || theonlyk == i -1|| theonlyk == i )){
            if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
              pushV(*born, i);
            onlysubcounter++;
          }
        }}
    }
    if (calc_real) {
      i = theonlyabove ;
      for (; i <= max(0, MH()->N()) ; ++i ) {
        for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born = thePureMEsMap[i].begin() ; born != thePureMEsMap[i].end() ; ++born ) {
          if ( i <= MH()->M() + 1 && i > 0 && !theonlyvirtualNLOParts&&!theonlyUnlopsweights&&( theonlyk == -1  || theonlyk == i -2|| theonlyk == i-1 )){
            if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
              pushProR(*born, i);
            onlysubcounter++;
          }
        }
      }
    }
    
    
    if ( !externalAmplitudes().empty() ) {
      generator()->log() << "Initializing external amplitudes.\n" << flush;
      boost::progress_display * progressBar =
      new boost::progress_display(externalAmplitudes().size(),generator()->log());
      for ( set<Ptr<MatchboxAmplitude>::tptr>::const_iterator ext =
           externalAmplitudes().begin(); ext != externalAmplitudes().end(); ++ext ) {
        if ( !(**ext).initializeExternal() ) {
          throw InitException()
          << "error: failed to initialize amplitude '" << (**ext).name() << "'\n";
        }
        ++(*progressBar);
      }
      delete progressBar;
      generator()->log() << "--------------------------------------------------------------------------------\n"
      << flush;
    }
    
    
    if ( !olpProcesses().empty() ) {
      generator()->log() << "Initializing one-loop provider(s).\n" << flush;
      map<Ptr<MatchboxAmplitude>::tptr, map<pair<Process, int>, int> > olps;
      for ( map<Ptr<MatchboxAmplitude>::tptr, map<pair<Process, int>, int> >::const_iterator oit = olpProcesses().begin() ; oit != olpProcesses().end() ;
           ++oit ) {
        olps[oit->first] = oit->second;
      }
      boost::progress_display * progressBar = new boost::progress_display(olps.size(), generator()->log());
      for ( map<Ptr<MatchboxAmplitude>::tptr, map<pair<Process, int>, int> >::const_iterator olpit = olps.begin() ; olpit != olps.end() ; ++olpit ) {
        if ( !olpit->first->startOLP(olpit->second) ) {
          throw InitException() << "error: failed to start OLP for amplitude '" << olpit->first->name() << "'\n";
        }
        ++(*progressBar);
      }
      delete progressBar;
      generator()->log() << "--------------------------------------------------------------------------------\n" << flush;
    }
    
    
    cout<<"\n Generated "<<MEs().size()<<" Subprocesses."<<flush;
    if(theonlysub!=-1)cout<<" ( "<<theonlysub<<"/"<<onlysubcounter<<" )"<<flush;
    cout<<"\n"<<flush;
    generator()->log() << "process setup finished.\n" << flush;
    
    
    
    ransetup=true;
    
  }
  
  
}

void MFactory::persistentOutput(PersistentOStream & os) const {
  
  
  os
  << calc_born            << calc_virtual       << calc_real
  << theonlyUnlopsweights << theonlyk           << theonlymulti
  << divideSub            << divideSubNumber    << theonlyabove
  << theonlysub           << theSubtractionData << ransetup
  << processMap           << theMergingHelper   <<theM<<theN;
}

void MFactory::persistentInput(PersistentIStream & is, int) {
  is
  >> calc_born            >> calc_virtual       >> calc_real
  >> theonlyUnlopsweights >> theonlyk           >> theonlymulti
  >> divideSub            >> divideSubNumber    >> theonlyabove
  >> theonlysub           >> theSubtractionData >> ransetup
  >> processMap           >> theMergingHelper   >>theM>>theN;
}

void MFactory::Init() {
  
  
  
  static Parameter<MFactory, int> interfaceonlyk("onlyk", "calculate only the ME with k additional partons.", &MFactory::theonlyk, -1, -1, 0,
                                                 false, false, Interface::lowerlim);
  static Parameter<MFactory, int> interfaceonlymulti("onlymulti", "calculate only the ME with k additional partons.", &MFactory::theonlymulti, -1, -1, 0,
                                                     false, false, Interface::lowerlim);
  static Parameter<MFactory, int> interfaceonlyabove("onlyabove", "calculate only the ME with more than k additional partons.", &MFactory::theonlyabove, -1, -1, 0,
                                                     false, false, Interface::lowerlim);
  
  static Parameter<MFactory, int> interfaceonlysub("onlysub", "calculate only one subProcess. this is for building grids.", &MFactory::theonlysub, -1, -1, 0,
                                                   false, false, Interface::lowerlim);
  
  
  
  
  
  
  
  
  
  static Parameter<MFactory, int> interfacedivideSub("divideSub", "calculate only one subProcess. this is for building grids.", &MFactory::divideSub, -1, -1, 0,
                                                     false, false, Interface::lowerlim);
  
  
  static Parameter<MFactory, int> interfacedivideSubNumber("divideSubNumber", "calculate only one subProcess. this is for building grids.", &MFactory::divideSubNumber, -1, -1, 0,
                                                           false, false, Interface::lowerlim);
  
  
  
  
  
  static Switch<MFactory, bool> interface_calc_born("calc_born", "[debug] Switch on or off the born contribution.", &MFactory::calc_born, true,
                                                    false, false);
  static SwitchOption interface_calc_bornOn(interface_calc_born, "On", "Switch on calculation of born.", true);
  static SwitchOption interface_calc_bornOff(interface_calc_born, "Off", "Switch off calculation of born.", false);
  
  static Switch<MFactory, bool> interface_calc_virtual("calc_virtual", "[debug] Switch on or off the virtual contribution.",
                                                       &MFactory::calc_virtual, true, false, false);
  static SwitchOption interface_calc_virtualOn(interface_calc_virtual, "On", "Switch on calculation of virtual.", true);
  static SwitchOption interface_calc_virtualOff(interface_calc_virtual, "Off", "Switch off calculation of virtual.", false);
  
  static Switch<MFactory, bool> interface_calc_real("calc_real", "[debug] Switch on or off the real contribution.", &MFactory::calc_real, true,
                                                    false, false);
  static SwitchOption interface_calc_realOn(interface_calc_real, "On", "Switch on calculation of real.", true);
  static SwitchOption interface_calc_realOff(interface_calc_real, "Off", "Switch off calculation of real.", false);
  
  
  
  
  static Switch<MFactory, bool> interface_theonlyNLOParts("onlyNLOParts", "Switch on or off the onlyNLOParts.", &MFactory::theonlyNLOParts, true, false,
                                                          false);
  static SwitchOption interface_theonlyNLOPartsOn(interface_theonlyNLOParts, "On", "Switch on the theonlyNLOParts.", true);
  static SwitchOption interface_theonlyNLOPartsOff(interface_theonlyNLOParts, "Off", "Switch off the theonlyNLOParts.", false);
  
  
    //  theonlyvirtualNLOParts;
    //  theonlyrealNLOParts;
  
  static Switch<MFactory, bool> interface_theonlyvirtualNLOParts("onlyvirtualNLOParts", "Switch on or off the onlyvirtualNLOParts.", &MFactory::theonlyvirtualNLOParts, true, false,
                                                                 false);
  static SwitchOption interface_theonlyvirtualNLOPartsOn(interface_theonlyvirtualNLOParts, "On", "Switch on the theonlyvirtualNLOParts.", true);
  static SwitchOption interface_theonlyvirtualNLOPartsOff(interface_theonlyvirtualNLOParts, "Off", "Switch off the theonlyvirtualNLOParts.", false);
  
  static Switch<MFactory, bool> interface_theonlyrealNLOParts("onlyrealNLOParts", "Switch on or off the onlyrealNLOParts.", &MFactory::theonlyrealNLOParts, true, false,
                                                              false);
  static SwitchOption interface_theonlyrealNLOPartsOn(interface_theonlyrealNLOParts, "On", "Switch on the theonlyrealNLOParts.", true);
  static SwitchOption interface_theonlyrealNLOPartsOff(interface_theonlyrealNLOParts, "Off", "Switch off the theonlyrealNLOParts.", false);
  
  static Switch<MFactory, bool> interface_theunitarizeNLOParts("unitarizeNLOParts", "Switch on or off the unitarizeNLOParts.", &MFactory::theunitarizeNLOParts, true, false,
                                                               false);
  static SwitchOption interface_theunitarizeNLOPartsOn(interface_theunitarizeNLOParts, "On", "Switch on the unitarizeNLOParts.", true);
  static SwitchOption interface_theunitarizeNLOPartsOff(interface_theunitarizeNLOParts, "Off", "Switch off the unitarizeNLOParts.", false);
  
  
  static Switch<MFactory, bool> interface_theonlyUnlopsweights("onlyUnlopsweights", "Switch on or off the onlyUnlopsweights.", &MFactory::theonlyUnlopsweights, true, false,
                                                               false);
  static SwitchOption interface_theonlyUnlopsweightsOn(interface_theonlyUnlopsweights, "On", "Switch on the onlyUnlopsweights.", true);
  static SwitchOption interface_theonlyUnlopsweightsOff(interface_theonlyUnlopsweights, "Off", "Switch off the onlyUnlopsweights.", false);
  
  
  
  static Switch<MFactory, bool> interface_Unitarized("Unitarized", "Unitarize the cross section.", &MFactory::unitarized, true, false, false);
  static SwitchOption interface_UnitarizedOn(interface_Unitarized, "On", "Switch on the unitarized cross section.", true);
  static SwitchOption interface_UnitarizedOff(interface_Unitarized, "Off", "Switch off the unitarized cross section.", false);
  
  
  
  
  static Switch<MFactory, bool> interface_NLOUnitarized("NLOUnitarized", "Unitarize the cross section.", &MFactory::NLOunitarized, true, false, false);
  static SwitchOption interface_NLOUnitarizedOn(interface_NLOUnitarized, "On", "Switch on the unitarized NLO cross section.", true);
  static SwitchOption interface_NLOUnitarizedOff(interface_NLOUnitarized, "Off", "Switch off the unitarized NLO cross section.", false);
  
  
  
  
  static Parameter<MFactory, string> interfaceSubtractionData("SubtractionData", "Prefix for subtraction check data.",
                                                              &MFactory::theSubtractionData, "", false, false);
  
  
  
  
  static Reference<MFactory,Merger> interfaceMergingHelper
  ("MergingHelper",
   "",
   &MFactory::theMergingHelper, false, false, true, true, false);
  
  
  
  static Parameter<MFactory, int> interfaceaddNLOLegs("NLOProcesses",
                                                    "Set the number of virtual corrections to consider. 0 is default for no virtual correction.", &MFactory::theM, 0, 0, 0, false, false,
                                                    Interface::lowerlim);
  
  
  
}

  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<MFactory, Herwig::MatchboxFactory> describeHerwigMFactory("Herwig::MFactory", "HwDipoleShower.so");
