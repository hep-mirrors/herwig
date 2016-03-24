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
  #include "MergeFactory.h"
  #include "ClusterNode.h"

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

  #include <boost/progress.hpp>

  #include <iterator>
  using std::ostream_iterator;

  using namespace Herwig;
  using std::ostream_iterator;

  MergeFactory::MergeFactory() :MatchboxFactory(),
  theonlyNLOParts(false),
  theonlyvirtualNLOParts(false),
  theonlyrealNLOParts(false),
  theunitarizeNLOParts(true),
  calc_born(true),
  calc_virtual(true),
  calc_real(true),
  calc_projected_virtual(true),
  calc_inclusive_real(true),
  calc_double_inclusive(true),
  needFOH(false),
  unitarized(true),
  NLOunitarized(true),
  theonlyUnlopsweights(false),
  theonlyk(-1),
  theonlymulti(-1),
  divideSub(-1),
  divideSubNumber(-1),
  theonlyabove(0),
  theonlysub(-1),
  theChooseHistory(3),
  theN(1),
  theM(-1),
  theMergePTsmearing(0.),
  theIRSafeRATIO(100),
  theStairFactor(0.),
  theMergePT(2.*GeV),
  theNLOMergePT(2.*GeV),
  theIRSafePT(1000000.0 * GeV),ransetup(false){
  }

  MergeFactory::~MergeFactory() {
  }

  IBPtr MergeFactory::clone() const {
    return new_ptr(*this);
  }

  IBPtr MergeFactory::fullclone() const {
    return new_ptr(*this);
  }

  void MergeFactory::fill_amplitudes() {
    
    olpProcesses().clear();
    
    if ( particleGroups().find("j") == particleGroups().end() ) throw InitException() << "Could not find a jet particle group named 'j'";
    
      // rebind the particle data objects
    for ( map<string, PDVector>::iterator g = particleGroups().begin() ; g != particleGroups().end() ; ++g ) {
      for ( PDVector::iterator p = g->second.begin() ; p != g->second.end() ; ++p ) {
  #ifndef NDEBUG
        long checkid = (**p).id();
  #endif
        *p = getParticleData((**p).id());
        assert((**p).id() == checkid);
      }
    }
    
    const PDVector& partons = particleGroups()["j"];
    
    unsigned int nl = 0;
    
    for ( PDVector::const_iterator p = partons.begin() ; p != partons.end() ; ++p )
      if ( abs((**p).id()) < 6 ) ++nl;
    
    nLight(nl / 2);
    processMap[0] = getProcesses()[0];
    if(theM>=0&&!theonlyrealNLOParts)
      setHighestVirt(processMap[0].size()+theM);
    
    for ( int i = 1 ; i <= theN ; ++i ) {
      processMap[i] = processMap[i - 1];
      if (particleGroups().find("jm") != particleGroups().end()){
        for(vector<string>::iterator it= processMap[i].begin();it!= processMap[i].end();it++){
          if(*it=="j")*it="jm";
          if(*it=="p")*it="pm";
        }
        processMap[i].push_back("jm");
      }else{
        processMap[i].push_back("j");
      }
            
    }
    
    for ( int i = 0 ; i <= theN ; ++i ) {
      vector<Ptr<MatchboxMEBase>::ptr> ames = makeMEs(processMap[i], orderInAlphaS() + i,i<theM+1);
      copy(ames.begin(), ames.end(), back_inserter(pureMEsMap()[i]));
    }
  }

  void MergeFactory::prepare_BV(int i) {
      // check if we have virtual contributions
    bool haveVirtuals = true;
    bool virtualsAreDR = false;
    bool virtualsAreDRbar = false;
    bool virtualsAreCDR = false;
    bool virtualsAreCS = false;
    bool virtualsAreBDK = false;
    bool virtualsAreExpanded = false;
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
      
      haveVirtuals &= (**born).haveOneLoop();
      if ( (**born).haveOneLoop() ) {
        virtualsAreDR |= (**born).isDR();
        virtualsAreCDR |= !(**born).isDR();
        virtualsAreDRbar |= (**born).isDRbar();
        virtualsAreCS |= (**born).isCS();
        virtualsAreBDK |= (**born).isBDK();
        virtualsAreExpanded |= (**born).isExpanded();
      }
    }
    
      // check the additional insertion operators
    if ( !theVirtualsMap[i].empty() ) haveVirtuals = true;
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = theVirtualsMap[i].begin() ; virt != theVirtualsMap[i].end() ; ++virt ) {
      (**virt).factory(this);
      virtualsAreDR |= (**virt).isDR();
      virtualsAreDRbar |= (**virt).isDRbar();
      virtualsAreCDR |= !(**virt).isDR();
      virtualsAreCS |= (**virt).isCS();
      virtualsAreBDK |= (**virt).isBDK();
      virtualsAreExpanded |= (**virt).isExpanded();
    }
    
    
      // check for consistent conventions on virtuals, if we are to include them
    if ( i <= theM ) {
      if ( virtualsAreDR && virtualsAreCDR ) {
        throw InitException() << "Virtual corrections use inconsistent regularization schemes.\n";
      }
      if ( (virtualsAreCS && virtualsAreBDK) || (virtualsAreCS && virtualsAreExpanded) || (virtualsAreBDK && virtualsAreExpanded)
          || (!virtualsAreCS && !virtualsAreBDK && !virtualsAreExpanded) ) {
        throw InitException() << "Virtual corrections use inconsistent conventions on finite terms.\n";
      }
      if ( !haveVirtuals ) {
        throw InitException() << "Could not find amplitudes for all virtual contributions needed.\n";
      }
    }
    
      // prepare dipole insertion operators
    if ( i <= theM ) {
      for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
           = DipoleRepository::insertionIOperators(dipoleSet()).begin();
           virt != DipoleRepository::insertionIOperators(dipoleSet()).end(); ++virt ) {
        (**virt).factory(this);
        if ( virtualsAreDRbar )
          (**virt).useDRbar();
        if ( virtualsAreDR )
          (**virt).useDR();
        else
          (**virt).useCDR();
        if ( virtualsAreCS )
          (**virt).useCS();
        if ( virtualsAreBDK )
          (**virt).useBDK();
        if ( virtualsAreExpanded )
          (**virt).useExpanded();
      }
      for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
           = DipoleRepository::insertionPKOperators(dipoleSet()).begin();
           virt != DipoleRepository::insertionPKOperators(dipoleSet()).end(); ++virt ) {
        
        (**virt).factory(this);
        if ( virtualsAreDRbar )
          (**virt).useDRbar();
        if ( virtualsAreDR )
          (**virt).useDR();
        else
          (**virt).useCDR();
        if ( virtualsAreCS )
          (**virt).useCS();
        if ( virtualsAreBDK )
          (**virt).useBDK();
        if ( virtualsAreExpanded )
          (**virt).useExpanded();
      }
      
    }
    thevirtualsAreExpandedMap[i] = virtualsAreExpanded;
    
  }

  void MergeFactory::prepare_R(int i) {
    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator real = thePureMEsMap[i].begin() ; real != thePureMEsMap[i].end() ; ++real ) {
      prepareME(*real);
    }
  }

  void MergeFactory::pushB(Ptr<MatchboxMEBase>::ptr born, int i, bool projected) {
    Ptr<MatchboxMEBase>::ptr bornme = (*born).cloneMe();
    bornme->maxMultCKKW(1);
    bornme->minMultCKKW(0);

    
    string pname = fullName() + "/" + bornme->name() + ".Born";
    if ( projected ) pname += "pro";
    if ( !(generator()->preinitRegister(bornme, pname)) ) throw InitException() << "Born ME " << pname << " already existing.";
    
    Ptr<ClusterNode>::ptr clusternode = new_ptr(ClusterNode(bornme, i, 1, needFOH));
    clusternode->mergePt(theMergePT);
    clusternode->centralMergePt(theMergePT);
    clusternode->N(theN + getProcesses()[0].size());clusternode->N0( getProcesses()[0].size());
    clusternode->M(theM + getProcesses()[0].size());
    if(theonlyk!=-1)clusternode->onlyN(theonlyk + getProcesses()[0].size());
    clusternode->unitarized(unitarized);
    clusternode->NLOunitarized(NLOunitarized);
    clusternode->calculateInNode(true);
    clusternode->treefactory(this);
    bornme->firstNode(clusternode);
    if ( projected ) bornme->projectorStage(1);
    else bornme->projectorStage(0);
    bornme->firstNode()->deepHead(clusternode);
    bornme->firstNode()->chooseHistory(theChooseHistory);
    
    vector<Ptr<ClusterNode>::ptr> temp;
    vector<Ptr<ClusterNode>::ptr> temp1;
    temp.push_back(bornme->firstNode());
    unsigned int k = 1;
    while (thePureMEsMap[i - k].size() != 0) {
      for ( unsigned int j = 0 ; j < temp.size() ; j++ ) {
        temp[j]->birth(thePureMEsMap[i - k]);
        for ( unsigned int m = 0 ; m < temp[j]->children().size() ; ++m ) {
          temp[j]->children()[m]->numberOfSplittings(temp[j]->children()[m]->nodeME()->numberOfSplittings(DipoleRepository::dipoles(dipoleSet()),thePureMEsMap[i - k+1]));

          
          
          for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionIOperators(dipoleSet()).begin() ;
               virt != DipoleRepository::insertionIOperators(dipoleSet()).end() ; ++virt ) {
            if ( (**virt).apply((*(temp[j]->children()[m]->nodeME())).diagrams().front()->partons()) ){
              Ptr<MatchboxInsertionOperator>::ptr myIOP = (**virt).cloneMe();
              ostringstream pname;
              pname <<  temp[j]->children()[m]->nodeME()->fullName()  << "/" << (**virt).name();
              if ( ! (generator()->preinitRegister(myIOP,pname.str()) ) )
                throw Exception() << "MatchboxMEBase::cloneDependencies(): Insertion operator " << pname.str() << " already existing." << Exception::runerror;
              temp[j]->children()[m]->nodeME()->virtuals().push_back(myIOP);
            }
          }
          for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt = DipoleRepository::insertionPKOperators(dipoleSet()).begin() ;
               virt != DipoleRepository::insertionPKOperators(dipoleSet()).end() ; ++virt ) {
            if ( (**virt).apply((*(temp[j]->children()[m]->nodeME())).diagrams().front()->partons()) ){
              Ptr<MatchboxInsertionOperator>::ptr myIOP = (**virt).cloneMe();
              ostringstream pname;
              pname <<  temp[j]->children()[m]->nodeME()->fullName()  << "/" << (**virt).name();
              if ( ! (generator()->preinitRegister(myIOP,pname.str()) ) )
                throw Exception() << "MatchboxMEBase::cloneDependencies(): Insertion operator " << pname.str() << " already existing." << Exception::runerror;
              temp[j]->children()[m]->nodeME()->virtuals().push_back(myIOP);
            }
          }
          
          
          temp[j]->children()[m]->nodeME()->noOneLoop();
          
          temp1.push_back(temp[j]->children()[m]);
        }
      }
      temp = temp1;
      k++;
    }
    
    if(theN>i)
         bornme->needsCorrelations();
    else bornme->needsNoCorrelations();
    
    
    bornme->cloneDependencies();
    if ( theonlyk == -1  || (theonlyk == i-1 &&unitarized) || theonlyk == i) {
      MEs().push_back(bornme);
    }
  }





  void MergeFactory::pushV(Ptr<MatchboxMEBase>::ptr born, int i, bool projected) {
    
    Ptr<MatchboxMEBase>::ptr nlo = (*born).cloneMe();
    nlo->maxMultCKKW(1);
    nlo->minMultCKKW(0);
    int pro = projected ? 1 : 0;
    
    string pname = fullName() + "/" + nlo->name() + ".Virtual";
    if ( projected ) pname += "pro";
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
      if ( checkPoles() ) {
        if ( !thevirtualsAreExpandedMap[i] ) {
          throw InitException() << "Cannot check epsilon poles if virtuals are not in `expanded' convention.\n";
        }
      }
    }
    nlo->doOneLoopNoBorn();
      ////////////////////////////////////NLO///////////////////////////
    Ptr<ClusterNode>::ptr clusternode = new_ptr(ClusterNode(nlo, i, 0, needFOH));
    clusternode->mergePt(theMergePT);
    clusternode->centralMergePt(theMergePT);
    clusternode->unitarized(unitarized);
    clusternode->NLOunitarized(NLOunitarized);
    clusternode->N(theN + getProcesses()[0].size());clusternode->N0( getProcesses()[0].size());
    clusternode->M(theM + getProcesses()[0].size());
    if(theonlyk!=-1)clusternode->onlyN(theonlyk + getProcesses()[0].size());
    clusternode->calculateInNode(true);
    clusternode->virtualContribution(true);
    clusternode->treefactory(this);
    nlo->firstNode(clusternode);
    if ( projected ) nlo->projectorStage(1);
    else nlo->projectorStage(0);
    nlo->firstNode()->deepHead(clusternode);
    nlo->firstNode()->chooseHistory(theChooseHistory);
    vector<Ptr<ClusterNode>::ptr> temp;
    vector<Ptr<ClusterNode>::ptr> temp1;
    temp.push_back(nlo->firstNode());
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
    if ( theonlyk == -1 || theonlyk == i - pro|| theonlyk == i - pro-1  ) {
      MEs().push_back(nlo);
    }
  }

  void MergeFactory::pushProR(Ptr<MatchboxMEBase>::ptr born, int i, bool projected,bool finiteDipoles) {
    Ptr<MatchboxMEBase>::ptr bornme = (*born).cloneMe();
    bornme->maxMultCKKW(1);
    bornme->minMultCKKW(0);
      //Ptr<SubtractedME>::ptr sub = new_ptr(SubtractedME());
    int pro = projected ? 2 : 1;
    
    string pname = fullName() + "/" + bornme->name() + (!finiteDipoles?".proReal":"finDI");
    if ( projected ) pname += "pro";
    
    if ( !(generator()->preinitRegister(bornme, pname)) ) throw InitException() << "Subtracted ME " << pname << " already existing.";
    
    if ( projected ) {
      bornme->projectorStage(2);
    } else {
      bornme->projectorStage(1);
    }
    
    Ptr<ClusterNode>::ptr clusternode = new_ptr(ClusterNode(bornme, i - 2, 1, needFOH));
    clusternode->mergePt(theMergePT);
    clusternode->centralMergePt(theMergePT);
    clusternode->N(theN + getProcesses()[0].size());clusternode->N0( getProcesses()[0].size());
    clusternode->M(theM + getProcesses()[0].size());
    if(theonlyk!=-1)clusternode->onlyN(theonlyk + getProcesses()[0].size());
    clusternode->calculateInNode(true);
    clusternode->subtractedReal(true);
    clusternode->finiteDipoles(finiteDipoles);
    clusternode->unitarized(unitarized);
    clusternode->NLOunitarized(NLOunitarized);
    clusternode->treefactory(this);
    bornme->firstNode(clusternode);
    if ( projected ) bornme->projectorStage(2);
    else bornme->projectorStage(1);
    bornme->firstNode()->deepHead(clusternode);
    bornme->firstNode()->chooseHistory(theChooseHistory);
    vector<Ptr<ClusterNode>::ptr> temp;
    vector<Ptr<ClusterNode>::ptr> temp1;
    temp.push_back(bornme->firstNode());
    unsigned int k = 1;
    while (thePureMEsMap[i - k].size() != 0) {
      for ( unsigned int j = 0 ; j < temp.size() ; j++ ) {
        temp[j]->birth(thePureMEsMap[i - k]);
        for ( unsigned int m = 0 ; m < temp[j]->children().size() ; ++m ) {
          
            // 	if ( temp[j]->children()[m]->nodeME()->isOLPTree() ) {
            // 	int id = orderOLPProcess(temp[j]->children()[m]->nodeME()->subProcess(),
            // 				 temp[j]->children()[m]->nodeME()->matchboxAmplitude(),
            // 				 ProcessType::treeME2);
            // 	temp[j]->children()[m]->nodeME()->olpProcess(ProcessType::treeME2,id);
            // 	id = orderOLPProcess(temp[j]->children()[m]->nodeME()->subProcess(),
            // 				       temp[j]->children()[m]->nodeME()->matchboxAmplitude(),
            // 				       ProcessType::colourCorrelatedME2);
            // 	temp[j]->children()[m]->nodeME()->olpProcess(ProcessType::colourCorrelatedME2,id);
            //
            // 	bool haveGluon = false;
            // 	for ( PDVector::const_iterator p = temp[j]->children()[m]->nodeME()->subProcess().legs.begin();
            // 	      p != temp[j]->children()[m]->nodeME()->subProcess().legs.end(); ++p )
            // 	  if ( (**p).id() == 21 ) {
            // 	    haveGluon = true;
            // 	    break;
            // 	  }
            // 	if ( haveGluon ) {
            // 	  id = orderOLPProcess(temp[j]->children()[m]->nodeME()->subProcess(),
            // 			       temp[j]->children()[m]->nodeME()->matchboxAmplitude(),
            // 			       ProcessType::spinColourCorrelatedME2);
            // 	  temp[j]->children()[m]->nodeME()->olpProcess(ProcessType::spinColourCorrelatedME2,id);
            // 	}
            //
            //
            //     }
          
          
          temp1.push_back(temp[j]->children()[m]);
        }
      }
      temp = temp1;
      k++;
    }
    
    
      //     if ( bornme->isOLPTree() ) {
      // 	int id = orderOLPProcess(bornme->subProcess(),
      // 				 (*born).matchboxAmplitude(),
      // 				 ProcessType::treeME2);
      // 	bornme->olpProcess(ProcessType::treeME2,id);
      // 	id = orderOLPProcess(bornme->subProcess(),
      // 				       (*born).matchboxAmplitude(),
      // 				       ProcessType::colourCorrelatedME2);
      // 	bornme->olpProcess(ProcessType::colourCorrelatedME2,id);
      //
      // 	bool haveGluon = false;
      // 	for ( PDVector::const_iterator p = (*born).subProcess().legs.begin();
      // 	      p != (*born).subProcess().legs.end(); ++p )
      // 	  if ( (**p).id() == 21 ) {
      // 	    haveGluon = true;
      // 	    break;
      // 	  }
      // 	if ( haveGluon ) {
      // 	  id = orderOLPProcess((*born).subProcess(),
      // 			       (*born).matchboxAmplitude(),
      // 			       ProcessType::spinColourCorrelatedME2);
      // 	  (*born).olpProcess(ProcessType::spinColourCorrelatedME2,id);
      // 	}
      //
      //
      //     }
    
    
    
    
    
    
    
    if(theN>i)
      bornme->needsCorrelations();
    else bornme->needsNoCorrelations();
    
    bornme->cloneDependencies(pname);
    
      ///////////////////////////////REAL/////////////////////////
    
      //sub->setBorns(thePureMEsMap[i - 1]);
      //sub->head(bornme);
      //sub->factory(this);
      //sub->dependent().clear();
      //sub->getDipoles();
      //if ( sub->dependent().empty() ) {
      //	cout << "???TreeMerge????" << flush;
      //	abort();
      //}
    
      ///////////////////////////////REAL/////////////////////////
    
      //if ( projected ) {
      
      //  Ptr<ReweightNumber>::ptr rw1 = new_ptr(ReweightNumber(1.));
        //sub->head()->addReweighter(rw1);
        //bornme->addReweighter(rw1);
        //} else {
        //Ptr<ReweightNumber>::ptr rw1 = new_ptr(ReweightNumber(1.));
        //sub->head()->addReweighter(rw1);
        //bornme->addReweighter(rw1);
        // }
    if ( theonlyk == -1 || theonlyk == i - pro || theonlyk == i - pro-1 ) {
      MEs().push_back(bornme);
    }
  }

  void MergeFactory::orderOLPs() {
    
  }

  void MergeFactory::setup() {
   
    useMe();
 
    if(!ransetup){

    olpProcesses().clear();
    externalAmplitudes().clear();
    setFixedCouplings(true);
    setFixedQEDCouplings(true);
    
    
    
    const PDVector& partons = particleGroups()["j"];
    unsigned int nl = 0;
    
    for ( PDVector::const_iterator p = partons.begin();
         p != partons.end(); ++p ) {
      if ( abs((**p).id()) < 7 && (**p).mass() == ZERO )
        ++nl;
      if ( (**p).id() > 0 && (**p).id() < 7 && (**p).mass() == ZERO )
        nLightJetVec( (**p).id() );
      if ( (**p).id() > 0 && (**p).id() < 7 && (**p).mass() != ZERO )
        nHeavyJetVec( (**p).id() );
    }
    nLight(nl/2);
    
    const PDVector& partonsInP = particleGroups()["p"];
    for ( PDVector::const_iterator pip = partonsInP.begin();
         pip != partonsInP.end(); ++pip ) {
      if ( (**pip).id() > 0 && (**pip).id() < 7 && (**pip).mass() == ZERO )
        nLightProtonVec( (**pip).id() );
    }
    
    
    for ( vector<Ptr<MatchboxAmplitude>::ptr>::iterator amp
         = amplitudes().begin(); amp != amplitudes().end(); ++amp )
      (**amp).factory(this);
    largeNBasis()->factory(this);
    
      // 	assert(!(theonlyk!=-1&&theonlysub!=-1));
      //         assert(!(theonlyk!=-1&&divideSub!=-1));
    assert(!(divideSub!=-1&&divideSubNumber==-1)||!(divideSub==-1&&divideSubNumber!=-1));
    assert(!subProcessGroups());
    
      //fill the amplitudes
    if ( !amplitudes().empty() ) {
      fill_amplitudes();
    }
    
      // prepare the Born and virtual matrix elements
    for ( int i = 0 ; i <= max(0, theN) ; ++i ) {
      prepare_BV(i);
    }
      // prepare the real emission matrix elements
    for ( int i = 0 ; i <= theN ; ++i ) {
      prepare_R(i);
    }
    
    orderOLPs();
    
      // start creating matrix elements
    MEs().clear();
    int onlysubcounter=0;
    int i = theonlyabove ;
    
    cout<<"\nstart filling LO "<<i<<" "<<theonlymulti<<flush;
    
    
    for (; i <= max(0, theN) ; ++i ) {
      if(i==theonlymulti||theonlymulti==-1)
        for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born = thePureMEsMap[i].begin() ; born != thePureMEsMap[i].end() ; ++born ) {
          if (!theonlyNLOParts&&!theonlyUnlopsweights&&( theonlyk == -1  || theonlyk == i -1|| theonlyk == i )){
            if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
              pushB(*born, i, false);
            onlysubcounter++;
          }
        }
    }
    
    cout<<"\nstart filling NLO"<<flush;
      if (false) {
        
      
    i = theonlyabove ;
    for (; i <=max(0, theN); ++i ) {
      if(i==theonlymulti||theonlymulti==-1)
        for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born = thePureMEsMap[i].begin() ; born != thePureMEsMap[i].end() ; ++born ) {
          
          if ( i <= theM  && !theonlyrealNLOParts&&( theonlyk == -1  || theonlyk == i -1|| theonlyk == i )){
            if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
              pushV(*born, i, false);
            onlysubcounter++;
          }
        }}
    i = theonlyabove ;
    for (; i <= max(0, theN) ; ++i ) {
      if((i)==theonlymulti+1||theonlymulti==-1)
        for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born = thePureMEsMap[i].begin() ; born != thePureMEsMap[i].end() ; ++born ) {
          if ( i <= theM + 1 && i > 0 && !theonlyvirtualNLOParts&&!theonlyUnlopsweights&&( theonlyk == -1  || theonlyk == i -2|| theonlyk == i-1 )){
            if(((theonlysub==-1||theonlysub==onlysubcounter)&&divideSub==-1)||(divideSub!=-1&&onlysubcounter%divideSub==divideSubNumber))
              pushProR(*born, i, false,false);
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

  void MergeFactory::persistentOutput(PersistentOStream & os) const {
    os<< theonlyNLOParts
    << theonlyvirtualNLOParts
    << theonlyrealNLOParts
    << theunitarizeNLOParts
    << calc_born
    << calc_virtual
    << calc_real
    << calc_projected_virtual
    << calc_inclusive_real
    << calc_double_inclusive
    << needFOH
    << unitarized
    << NLOunitarized
    << theonlyUnlopsweights
    << theonlyk
    << theonlymulti
    << divideSub
    << divideSubNumber
    << theonlyabove
    << theonlysub
    << theChooseHistory
    << theN
    << theM
    << theMergePTsmearing
    << theIRSafeRATIO
    << theStairFactor
    << ounit(theMergePT,GeV)
    << ounit(theNLOMergePT,GeV)
    << ounit(theIRSafePT,GeV)
    << theSubtractionData
    << theLargeNBasis
    << ransetup;
  }

  void MergeFactory::persistentInput(PersistentIStream & is, int) {
    is>> theonlyNLOParts
    >> theonlyvirtualNLOParts
    >> theonlyrealNLOParts
    >> theunitarizeNLOParts
    >> calc_born
    >> calc_virtual
    >> calc_real
    >> calc_projected_virtual
    >> calc_inclusive_real
    >> calc_double_inclusive
    >> needFOH
    >> unitarized
    >> NLOunitarized
    >> theonlyUnlopsweights
    >> theonlyk
    >> theonlymulti
    >> divideSub
    >> divideSubNumber
    >> theonlyabove
    >> theonlysub
    >> theChooseHistory
    >> theN
    >> theM
    >> theMergePTsmearing
    >> theIRSafeRATIO
    >> theStairFactor
    >> iunit(theMergePT,GeV)
    >> iunit(theNLOMergePT,GeV)
    >> iunit(theIRSafePT,GeV)
    >> theSubtractionData
    >> theLargeNBasis
    >> ransetup;
  }

  void MergeFactory::Init() {
    
    static Parameter<MergeFactory, int> interfaceadditionalN("additionalN", "Set the number of additional jets to consider.", &MergeFactory::theN, 0, 0,
                                                             0, false, false, Interface::lowerlim);
    
    static Parameter<MergeFactory, int> interfacevirtualM("virtualM",
                                                          "Set the number of virtual corrections to consider. -1 is default for no virtual correction.", &MergeFactory::theM, -1, -1, 0, false, false,
                                                          Interface::lowerlim);
    
    static Parameter<MergeFactory, int> interfaceonlyk("onlyk", "calculate only the ME with k additional partons.", &MergeFactory::theonlyk, -1, -1, 0,
                                                       false, false, Interface::lowerlim);
    static Parameter<MergeFactory, int> interfaceonlymulti("onlymulti", "calculate only the ME with k additional partons.", &MergeFactory::theonlymulti, -1, -1, 0,
                                                           false, false, Interface::lowerlim);
          static Parameter<MergeFactory, int> interfaceonlyabove("onlyabove", "calculate only the ME with more than k additional partons.", &MergeFactory::theonlyabove, -1, -1, 0,
                                                                 false, false, Interface::lowerlim);
    
    static Parameter<MergeFactory, int> interfaceonlysub("onlysub", "calculate only one subProcess. this is for building grids.", &MergeFactory::theonlysub, -1, -1, 0,
                                                         false, false, Interface::lowerlim);
    
    
    
    
    
    
    
    
    
    static Parameter<MergeFactory, int> interfacedivideSub("divideSub", "calculate only one subProcess. this is for building grids.", &MergeFactory::divideSub, -1, -1, 0,
                                                           false, false, Interface::lowerlim);
    
    
    static Parameter<MergeFactory, int> interfacedivideSubNumber("divideSubNumber", "calculate only one subProcess. this is for building grids.", &MergeFactory::divideSubNumber, -1, -1, 0,
                                                                 false, false, Interface::lowerlim);
    
    
    
    
    
    static Parameter<MergeFactory, int> interfacechooseHistory("chooseHistory", "different ways of choosing the history", &MergeFactory::theChooseHistory, 3, -1, 0,
                                                               false, false, Interface::lowerlim);
    
    
    static Switch<MergeFactory, bool> interface_calc_born("calc_born", "[debug] Switch on or off the born contribution.", &MergeFactory::calc_born, true,
                                                          false, false);
    static SwitchOption interface_calc_bornOn(interface_calc_born, "On", "Switch on calculation of born.", true);
    static SwitchOption interface_calc_bornOff(interface_calc_born, "Off", "Switch off calculation of born.", false);
    
    static Switch<MergeFactory, bool> interface_calc_virtual("calc_virtual", "[debug] Switch on or off the virtual contribution.",
                                                             &MergeFactory::calc_virtual, true, false, false);
    static SwitchOption interface_calc_virtualOn(interface_calc_virtual, "On", "Switch on calculation of virtual.", true);
    static SwitchOption interface_calc_virtualOff(interface_calc_virtual, "Off", "Switch off calculation of virtual.", false);
    
    static Switch<MergeFactory, bool> interface_calc_real("calc_real", "[debug] Switch on or off the real contribution.", &MergeFactory::calc_real, true,
                                                          false, false);
    static SwitchOption interface_calc_realOn(interface_calc_real, "On", "Switch on calculation of real.", true);
    static SwitchOption interface_calc_realOff(interface_calc_real, "Off", "Switch off calculation of real.", false);
    
    static Switch<MergeFactory, bool> interface_calc_projected_virtual("calc_projected_virtual",
                                                                       "[debug] Switch on or off the projected_virtual contribution.", &MergeFactory::calc_projected_virtual, true, false, false);
    static SwitchOption interface_calc_projected_virtualOn(interface_calc_projected_virtual, "On", "Switch on calculation of projected_virtual.", true);
    static SwitchOption interface_calc_projected_virtualOff(interface_calc_projected_virtual, "Off", "Switch off calculation of projected_virtual.", false);
    
    static Switch<MergeFactory, bool> interface_calc_inclusive_real("calc_inclusive_real", "[debug] Switch on or off the real contribution.",
                                                                    &MergeFactory::calc_inclusive_real, true, false, false);
    static SwitchOption interface_calc_inclusive_realOn(interface_calc_inclusive_real, "On", "Switch on calculation of inclusive_real.", true);
    static SwitchOption interface_calc_inclusive_realOff(interface_calc_inclusive_real, "Off", "Switch off calculation of inclusive_real.", false);
    
    static Switch<MergeFactory, bool> interface_calc_double_inclusive("calc_double_inclusive", "[debug] Switch on or off the real contribution.",
                                                                      &MergeFactory::calc_double_inclusive, true, false, false);
    static SwitchOption interface_calc_double_inclusiveOn(interface_calc_double_inclusive, "On", "Switch on calculation of double_inclusive.", true);
    static SwitchOption interface_calc_double_inclusiveOff(interface_calc_double_inclusive, "Off", "Switch off calculation of double_inclusive.", false);
    
    static Switch<MergeFactory, bool> interface_needFOH("needFOH", "Switch on or off the full ordered histories.", &MergeFactory::needFOH, true, false,
                                                        false);
    static SwitchOption interface_needFOHOn(interface_needFOH, "On", "Switch on the full ordered histories.", true);
    static SwitchOption interface_needFOHOff(interface_needFOH, "Off", "Switch off the full ordered histories.", false);
    
    static Switch<MergeFactory, bool> interface_theonlyNLOParts("onlyNLOParts", "Switch on or off the onlyNLOParts.", &MergeFactory::theonlyNLOParts, true, false,
                                                                false);
    static SwitchOption interface_theonlyNLOPartsOn(interface_theonlyNLOParts, "On", "Switch on the theonlyNLOParts.", true);
    static SwitchOption interface_theonlyNLOPartsOff(interface_theonlyNLOParts, "Off", "Switch off the theonlyNLOParts.", false);
    
    
      //  theonlyvirtualNLOParts;
      //  theonlyrealNLOParts;
    
    static Switch<MergeFactory, bool> interface_theonlyvirtualNLOParts("onlyvirtualNLOParts", "Switch on or off the onlyvirtualNLOParts.", &MergeFactory::theonlyvirtualNLOParts, true, false,
                                                                       false);
    static SwitchOption interface_theonlyvirtualNLOPartsOn(interface_theonlyvirtualNLOParts, "On", "Switch on the theonlyvirtualNLOParts.", true);
    static SwitchOption interface_theonlyvirtualNLOPartsOff(interface_theonlyvirtualNLOParts, "Off", "Switch off the theonlyvirtualNLOParts.", false);
    
    static Switch<MergeFactory, bool> interface_theonlyrealNLOParts("onlyrealNLOParts", "Switch on or off the onlyrealNLOParts.", &MergeFactory::theonlyrealNLOParts, true, false,
                                                                    false);
    static SwitchOption interface_theonlyrealNLOPartsOn(interface_theonlyrealNLOParts, "On", "Switch on the theonlyrealNLOParts.", true);
    static SwitchOption interface_theonlyrealNLOPartsOff(interface_theonlyrealNLOParts, "Off", "Switch off the theonlyrealNLOParts.", false);
    
    static Switch<MergeFactory, bool> interface_theunitarizeNLOParts("unitarizeNLOParts", "Switch on or off the unitarizeNLOParts.", &MergeFactory::theunitarizeNLOParts, true, false,
                                                                     false);
    static SwitchOption interface_theunitarizeNLOPartsOn(interface_theunitarizeNLOParts, "On", "Switch on the unitarizeNLOParts.", true);
    static SwitchOption interface_theunitarizeNLOPartsOff(interface_theunitarizeNLOParts, "Off", "Switch off the unitarizeNLOParts.", false);
    
    
    static Switch<MergeFactory, bool> interface_theonlyUnlopsweights("onlyUnlopsweights", "Switch on or off the onlyUnlopsweights.", &MergeFactory::theonlyUnlopsweights, true, false,
                                                                     false);
    static SwitchOption interface_theonlyUnlopsweightsOn(interface_theonlyUnlopsweights, "On", "Switch on the onlyUnlopsweights.", true);
    static SwitchOption interface_theonlyUnlopsweightsOff(interface_theonlyUnlopsweights, "Off", "Switch off the onlyUnlopsweights.", false);
    
    
    
    static Switch<MergeFactory, bool> interface_Unitarized("Unitarized", "Unitarize the cross section.", &MergeFactory::unitarized, true, false, false);
    static SwitchOption interface_UnitarizedOn(interface_Unitarized, "On", "Switch on the unitarized cross section.", true);
    static SwitchOption interface_UnitarizedOff(interface_Unitarized, "Off", "Switch off the unitarized cross section.", false);
    
    static Switch<MergeFactory, bool> interface_NLOUnitarized("NLOUnitarized", "Unitarize the cross section.", &MergeFactory::NLOunitarized, true, false, false);
    static SwitchOption interface_NLOUnitarizedOn(interface_NLOUnitarized, "On", "Switch on the unitarized NLO cross section.", true);
    static SwitchOption interface_NLOUnitarizedOff(interface_NLOUnitarized, "Off", "Switch off the unitarized NLO cross section.", false);
    
    static Parameter<MergeFactory, Energy> interfacemergept("mergept", "Set the pt to be merged to consider.", &MergeFactory::theMergePT, GeV, 2.0 * GeV,
                                                            ZERO, Constants::MaxEnergy, true, false, Interface::limited);
    
    
    
    static Parameter<MergeFactory, double> interfacemergeptsmearing("mergeptsmearing", "Set the percentage the merging pt should be smeared.",
                                                                    &MergeFactory::theMergePTsmearing, 0., 0.,
                                                                    0.0, 0.5, true, false, Interface::limited);	
    
    static Parameter<MergeFactory, double> interfacestairfactor("stairfactor", "Set the pt to be merged to consider.", &MergeFactory::theStairFactor, 1, 1,
                                                                1, 4, true, false, Interface::limited);
    interfacestairfactor.setHasDefault(false);
    
    
    
    interfacemergept.setHasDefault(false);
    static Parameter<MergeFactory, Energy> interfaceNLOmergept("NLOmergept", "NLO:Set the pt to be merged to consider.", &MergeFactory::theNLOMergePT, GeV, 2.0 * GeV,
                                                               ZERO, Constants::MaxEnergy, true, false, Interface::limited);
    interfaceNLOmergept.setHasDefault(false);
    
    static Parameter<MergeFactory, Energy> interfaceIRSafePT("IRSafePT", "Set the pt for which a matrix element should be treated as IR-safe.",
                                                             &MergeFactory::theIRSafePT, GeV, 0.0 * GeV, ZERO, Constants::MaxEnergy, true, false, Interface::limited);
    interfaceIRSafePT.setHasDefault(false);
    
    static Parameter<MergeFactory, double> interfaceIRSafeRATIO("IRSafeRATIO", "Set the RATIO for which a matrix element should be treated as IR-safe.",
                                                                &MergeFactory::theIRSafeRATIO, 100.0, 0.0, ZERO, false, false, Interface::lowerlim);
    interfaceIRSafeRATIO.setHasDefault(false);
    
    /*	static Command<MergeFactory> interfaceProcess("Process", "Set the process to consider.", &MergeFactory::doProcess, false);*/
    
    static Parameter<MergeFactory, string> interfaceSubtractionData("SubtractionData", "Prefix for subtraction check data.",
                                                                    &MergeFactory::theSubtractionData, "", false, false);
    
    
    static Reference<MergeFactory,ColourBasis> interfaceLargeNBasis
    ("LargeNBasis",
     "Set the large-N colour basis implementation.",
     &MergeFactory::theLargeNBasis, false, false, true, true, false);
    
    
  }

    // *** Attention *** The following static variable is needed for the type
    // description system in ThePEG. Please check that the template arguments
    // are correct (the class and its base class), and that the constructor
    // arguments are correct (the class name and the name of the dynamically
    // loadable library where the class implementation can be found).
  DescribeClass<MergeFactory, MatchboxFactory> describeHerwigMergeFactory("Herwig::MergeFactory", "Herwig.so");
