// -*- C++ -*-
//
// MatchboxFactory.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxFactory class.
//

#include "MatchboxFactory.h"
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

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"

using namespace Herwig;
using std::ostream_iterator;

MatchboxFactory::MatchboxFactory() 
  : SubProcessHandler(), theNLight(0),
    theOrderInAlphaS(0), theOrderInAlphaEW(0),
    theBornContributions(true), theVirtualContributions(true),
    theRealContributions(true), theSubProcessGroups(false),
    theFactorizationScaleFactor(1.0), theRenormalizationScaleFactor(1.0),
    theFixedCouplings(false), theVetoScales(false),
    theVerbose(false), theSubtractionData("") {}

MatchboxFactory::~MatchboxFactory() {}

IBPtr MatchboxFactory::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxFactory::fullclone() const {
  return new_ptr(*this);
}

void MatchboxFactory::prepareME(Ptr<MatchboxMEBase>::ptr me) const {

  Ptr<MatchboxAmplitude>::ptr amp =
    dynamic_ptr_cast<Ptr<MatchboxAmplitude>::ptr>((*me).amplitude());
  me->matchboxAmplitude(amp);

  if ( diagramGenerator() && !me->diagramGenerator() )
    me->diagramGenerator(diagramGenerator());

  if ( me->nLight() == 0 )
    me->nLight(nLight());

  if ( phasespace() && !me->phasespace() )
    me->phasespace(phasespace());

  if ( scaleChoice() && !me->scaleChoice() )
    me->scaleChoice(scaleChoice());

  if ( me->factorizationScaleFactor() == 1.0 )
    me->factorizationScaleFactor(factorizationScaleFactor());

  if ( me->renormalizationScaleFactor() == 1.0 )
    me->renormalizationScaleFactor(renormalizationScaleFactor());

  if ( fixedCouplings() )
    me->setFixedCouplings();

  if ( cache() && !me->cache() )
    me->cache(cache());

  if ( verbose() )
    me->setVerbose();

}

struct ProjectQN {
  inline pair<int,pair<int,int> > operator()(PDPtr p) const {
    return
      pair<int,pair<int,int> >((*p).iSpin(),pair<int,int>((*p).iCharge(),(*p).iColour()));
  }
};

string pid(const vector<pair<int,pair<int,int> > >& key) {
  ostringstream res;
  for ( vector<pair<int,pair<int,int> > >::const_iterator k =
	  key.begin(); k != key.end(); ++k )
    res << k->first << k->second.first
	<< k->second.second;
  return res.str();
}

vector<Ptr<MatchboxMEBase>::ptr> MatchboxFactory::
makeMEs(const vector<string>& proc, unsigned int orderas) const {

  typedef vector<pair<int,pair<int,int> > > QNKey;

  map<Ptr<MatchboxAmplitude>::ptr,map<QNKey,vector<PDVector> > > ampProcs;
  set<PDVector> processes = makeSubProcesses(proc);

  for ( vector<Ptr<MatchboxAmplitude>::ptr>::const_iterator amp
	  = amplitudes().begin(); amp != amplitudes().end(); ++amp ) {
    if ( (**amp).orderInGs() != orderas ||
	 (**amp).orderInGem() != orderInAlphaEW() )
      continue;
    for ( set<PDVector>::const_iterator p = processes.begin();
	  p != processes.end(); ++p ) {
      if ( !(**amp).canHandle(*p) )
	continue;
      QNKey key;
      transform(p->begin(),p->end(),back_inserter(key),ProjectQN());
      ampProcs[*amp][key].push_back(*p);
    }
  }

  vector<Ptr<MatchboxMEBase>::ptr> res;
  for ( map<Ptr<MatchboxAmplitude>::ptr,map<QNKey,vector<PDVector> > >::const_iterator
	  ap = ampProcs.begin(); ap != ampProcs.end(); ++ap ) {
    for ( map<QNKey,vector<PDVector> >::const_iterator m = ap->second.begin();
	  m != ap->second.end(); ++m ) {
      Ptr<MatchboxMEBase>::ptr me = new_ptr(MatchboxMEBase());
      me->subProcesses() = m->second;
      me->amplitude(ap->first);
      string pname = "ME" + ap->first->name() + pid(m->first);
      if ( ! (generator()->preinitRegister(me,pname) ) )
	throw InitException() << "Matrix element " << pname << " already existing.";
      res.push_back(me);
    }
  }

  return res;

}

void MatchboxFactory::setup() {

  if ( !amplitudes().empty() ) {

    if ( particleGroups().find("j") == particleGroups().end() )
      throw InitException() << "Could not find a jet particle group named 'j'";

    // rebind the particle data objects
    for ( map<string,PDVector>::iterator g = particleGroups().begin();
	  g != particleGroups().end(); ++g )
      for ( PDVector::iterator p = g->second.begin();
	    p != g->second.end(); ++p ) {
	long checkid = (**p).id();
	*p = getParticleData((**p).id());
	assert((**p).id() == checkid);
      }

    nLight(particleGroups()["j"].size());

    bornMEs() = makeMEs(process,orderInAlphaS());

    if ( realContributions() ) {
      vector<string> rproc = process;
      rproc.push_back("j");
      realEmissionMEs() = makeMEs(rproc,orderInAlphaS()+1);
    }

  }

  // check if we have virtual contributions
  bool haveVirtuals = true;

  // check DR conventions of virtual contributions
  bool virtualsAreDR = false;
  bool virtualsAreCDR = false;

  // check finite term conventions of virtual contributions
  bool virtualsAreCS = false;
  bool virtualsAreStandard = false;

  // check and prepare the Born and virtual matrix elements
  for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	  = bornMEs().begin(); born != bornMEs().end(); ++born ) {
    prepareME(*born);
    haveVirtuals &= (**born).haveOneLoop();
    if ( (**born).haveOneLoop() ) {
      virtualsAreDR |= (**born).isDR();
      virtualsAreCDR |= !(**born).isDR();
      virtualsAreCS |= (**born).isCS();
      virtualsAreStandard |= !(**born).isCS();
    }
  }

  // check the additional insertion operators
  if ( !virtuals().empty() )
    haveVirtuals = true;    
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
	  = virtuals().begin(); virt != virtuals().end(); ++virt ) {
    virtualsAreDR |= (**virt).isDR();
    virtualsAreCDR |= !(**virt).isDR();
    virtualsAreCS |= (**virt).isCS();
    virtualsAreStandard |= !(**virt).isCS();
  }

  // check for consistent conventions on virtuals, if we are to include them
  if ( virtualContributions() ) {
    if ( virtualsAreDR && virtualsAreCDR ) {
      throw InitException() << "Virtual corrections use inconsistent regularization schemes.\n";
    }
    if ( virtualsAreCS && virtualsAreStandard ) {
      throw InitException() << "Virtual corrections use inconsistent conventions on finite terms.\n";
    }
    if ( !haveVirtuals ) {
      throw InitException() << "Could not find amplitudes for all virtual contributions needed.\n";
    }
  }

  // prepare dipole insertion operators
  if ( virtualContributions() ) {
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
	    = DipoleRepository::insertionOperators().begin(); 
	  virt != DipoleRepository::insertionOperators().end(); ++virt ) {
      if ( virtualsAreDR )
	(**virt).useDR();
      else
	(**virt).useCDR();
      if ( virtualsAreCS )
	(**virt).useCS();
      else
	(**virt).useNonCS();
    }
  }

  // prepare the real emission matrix elements
  if ( realContributions() ) {
    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator real
	    = realEmissionMEs().begin(); real != realEmissionMEs().end(); ++real ) {
      prepareME(*real);
    }
  }

  // start creating matrix elements
  MEs().clear();

  // setup born and virtual contributions

  if ( !bornContributions() && virtualContributions() ) {
    throw InitException() << "Virtual corrections without Born contributions not yet supported.\n";
  }

  if ( bornContributions() && !virtualContributions() ) {
    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	    = bornMEs().begin(); born != bornMEs().end(); ++born ) {

      if ( (**born).onlyOneLoop() )
	continue;

      Ptr<MatchboxMEBase>::ptr bornme = (**born).cloneMe();
      string pname = fullName() + "/" + (**born).name();
      if ( ! (generator()->preinitRegister(bornme,pname) ) )
	throw InitException() << "Matrix element " << pname << " already existing.";
      bornme->cloneDependencies();
      MEs().push_back(bornme);

    }
  }

  if ( bornContributions() && virtualContributions() ) {

    bornVirtualMEs().clear();

    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	    = bornMEs().begin(); born != bornMEs().end(); ++born ) {

      Ptr<MatchboxNLOME>::ptr nlo = new_ptr(MatchboxNLOME());
      string pname = fullName() + "/" + (**born).name();
      if ( ! (generator()->preinitRegister(nlo,pname) ) )
	throw InitException() << "NLO ME " << pname << " already existing.";

      nlo->matrixElement(*born);
      nlo->virtuals().clear();

      if ( !nlo->matrixElement()->onlyOneLoop() ) {
	for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
		= virtuals().begin(); virt != virtuals().end(); ++virt ) {
	  if ( (**virt).apply((**born).diagrams().front()->partons()) )
	    nlo->virtuals().push_back(*virt);
	}
	for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
		= DipoleRepository::insertionOperators().begin(); 
	      virt != DipoleRepository::insertionOperators().end(); ++virt ) {
	  if ( (**virt).apply((**born).diagrams().front()->partons()) )
	    nlo->virtuals().push_back(*virt);
	}
	if ( nlo->virtuals().empty() )
	  throw InitException() << "No insertion operators have been found for "
				<< (**born).name() << "\n";
      }

      nlo->cloneDependencies();

      bornVirtualMEs().push_back(nlo);
      MEs().push_back(nlo);

    }

  }

  if ( realContributions() ) {

    if ( theSubtractionData != "" )
      if ( theSubtractionData[theSubtractionData.size()-1] != '/' )
	theSubtractionData += "/";

    subtractedMEs().clear();    

    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator real
	    = realEmissionMEs().begin(); real != realEmissionMEs().end(); ++real ) {

      Ptr<SubtractedME>::ptr sub = new_ptr(SubtractedME());
      string pname = fullName() + "/" + (**real).name();
      if ( ! (generator()->preinitRegister(sub,pname) ) )
	throw InitException() << "Subtracted ME " << pname << " already existing.";

      sub->borns() = bornMEs();
      sub->head(*real);

      sub->allDipoles().clear();
      sub->dependent().clear();

      sub->getDipoles();

      if ( verbose() )
	sub->setVerbose();

      if ( subProcessGroups() )
	sub->setSubProcessGroups();

      if ( vetoScales() )
	sub->doVetoScales();

      if ( subtractionData() != "" )
	sub->subtractionData(subtractionData());

      subtractedMEs().push_back(sub);

      MEs().push_back(sub);

    }

  }

}

void MatchboxFactory::print(ostream& os) const {

  os << "--- MatchboxFactory setup -----------------------------------------------------------\n";

  if ( !amplitudes().empty() ) {
    os << " generated Born matrix elements:\n";
    for ( vector<Ptr<MatchboxMEBase>::ptr>::const_iterator m = bornMEs().begin();
	  m != bornMEs().end(); ++m ) {
      os << " '" << (**m).name() << "'\n";
    }
    os << flush;
    os << " generated real emission matrix elements:\n";
    for ( vector<Ptr<MatchboxMEBase>::ptr>::const_iterator m = realEmissionMEs().begin();
	  m != realEmissionMEs().end(); ++m ) {
      os << " '" << (**m).name() << "'\n";
    }
    os << flush;
  }

  os << " generated Born+virtual matrix elements:\n";

  for ( vector<Ptr<MatchboxNLOME>::ptr>::const_iterator bv
	  = bornVirtualMEs().begin(); bv != bornVirtualMEs().end(); ++bv ) {
    (**bv).print(os);
  }

  os << " generated subtracted matrix elements:\n";

  for ( vector<Ptr<SubtractedME>::ptr>::const_iterator sub
	  = subtractedMEs().begin(); sub != subtractedMEs().end(); ++sub ) {
    os << " '" << (**sub).name() << "'\n";
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void MatchboxFactory::doinit() {
  setup();
  if ( theVerbose )
    print(Repository::clog());
  SubProcessHandler::doinit();
}


void MatchboxFactory::persistentOutput(PersistentOStream & os) const {
  os << theDiagramGenerator 
     << theNLight << theOrderInAlphaS << theOrderInAlphaEW 
     << theBornContributions << theVirtualContributions
     << theRealContributions << theSubProcessGroups
     << thePhasespace << theScaleChoice
     << theFactorizationScaleFactor << theRenormalizationScaleFactor
     << theFixedCouplings << theVetoScales
     << theAmplitudes << theCache
     << theBornMEs << theVirtuals << theRealEmissionMEs
     << theBornVirtualMEs << theSubtractedMEs
     << theVerbose << theSubtractionData
     << theParticleGroups << process;
}

void MatchboxFactory::persistentInput(PersistentIStream & is, int) {
  is >> theDiagramGenerator 
     >> theNLight >> theOrderInAlphaS >> theOrderInAlphaEW 
     >> theBornContributions >> theVirtualContributions
     >> theRealContributions >> theSubProcessGroups
     >> thePhasespace >> theScaleChoice
     >> theFactorizationScaleFactor >> theRenormalizationScaleFactor
     >> theFixedCouplings >> theVetoScales
     >> theAmplitudes >> theCache
     >> theBornMEs >> theVirtuals >> theRealEmissionMEs
     >> theBornVirtualMEs >> theSubtractedMEs
     >> theVerbose >> theSubtractionData
     >> theParticleGroups >> process;
}

string MatchboxFactory::startParticleGroup(string name) {
  particleGroupName = StringUtils::stripws(name);
  particleGroup.clear();
  return "";
}

string MatchboxFactory::endParticleGroup(string) {
  if ( particleGroup.empty() )
    throw InitException() << "Empty particle group.";
  particleGroups()[particleGroupName] = particleGroup;
  particleGroup.clear();
  return "";
}

string MatchboxFactory::doProcess(string in) {
  process = StringUtils::split(in);
  if ( process.size() < 3 )
    throw InitException() << "Invalid process.";
  for ( vector<string>::iterator p = process.begin();
	p != process.end(); ++p ) {
    *p = StringUtils::stripws(*p);
  }
  return "";
}

struct SortPID {
  inline bool operator()(PDPtr a, PDPtr b) const {
    return a->id() < b->id();
  }
};

set<PDVector> MatchboxFactory::
makeSubProcesses(const vector<string>& proc) const {

  if ( proc.empty() )
    throw InitException() << "No process specified.";

  vector<PDVector> allProcs(1);
  size_t pos = 0;
  typedef map<string,PDVector>::const_iterator GroupIterator;

  while ( pos < proc.size() ) {

    GroupIterator git =
      particleGroups().find(proc[pos]);

    if ( git == particleGroups().end() ) {
      throw InitException() << "particle group '"
			    << proc[pos] << "' not defined.";
    }

    vector<PDVector> mine;

    for ( vector<PDVector>::const_iterator i = allProcs.begin();
	  i != allProcs.end(); ++i ) {
      for ( PDVector::const_iterator p = git->second.begin();
	    p != git->second.end(); ++p ) {
	PDVector v = *i;
	v.push_back(*p);
	mine.push_back(v);
      }
    }

    allProcs = mine;
    ++pos;

  }

  set<PDVector> allCheckedProcs;
  for ( vector<PDVector>::const_iterator p = allProcs.begin();
	p != allProcs.end(); ++p ) {
    int charge = -(*p)[0]->iCharge() -(*p)[1]->iCharge();
    for ( size_t k = 2; k < (*p).size(); ++k )
      charge += (*p)[k]->iCharge();
    if ( charge != 0 )
      continue;
    PDVector pr = *p;
    sort(pr.begin()+2,pr.end(),SortPID());
    allCheckedProcs.insert(pr);
  }

  return allCheckedProcs;

}

void MatchboxFactory::Init() {

  static ClassDocumentation<MatchboxFactory> documentation
    ("MatchboxFactory",
     "NLO QCD corrections have been calculated "
     "using Matchbox \\cite{Platzer:2011bc}",
     "%\\cite{Platzer:2011bc}\n"
     "\\bibitem{Platzer:2011bc}\n"
     "S.~Platzer and S.~Gieseke,\n"
     "``Dipole Showers and Automated NLO Matching in Herwig++,''\n"
     "arXiv:1109.6256 [hep-ph].\n"
     "%%CITATION = ARXIV:1109.6256;%%");


  static Reference<MatchboxFactory,Tree2toNGenerator> interfaceDiagramGenerator
    ("DiagramGenerator",
     "Set the diagram generator.",
     &MatchboxFactory::theDiagramGenerator, false, false, true, true, false);


  static Parameter<MatchboxFactory,unsigned int> interfaceOrderInAlphaS
    ("OrderInAlphaS",
     "The order in alpha_s to consider.",
     &MatchboxFactory::theOrderInAlphaS, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxFactory,unsigned int> interfaceOrderInAlphaEW
    ("OrderInAlphaEW",
     "The order in alpha_EW",
     &MatchboxFactory::theOrderInAlphaEW, 2, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<MatchboxFactory,bool> interfaceBornContributions
    ("BornContributions",
     "Switch on or off the Born contributions.",
     &MatchboxFactory::theBornContributions, true, false, false);
  static SwitchOption interfaceBornContributionsOn
    (interfaceBornContributions,
     "On",
     "Switch on Born contributions.",
     true);
  static SwitchOption interfaceBornContributionsOff
    (interfaceBornContributions,
     "Off",
     "Switch off Born contributions.",
     false);

  static Switch<MatchboxFactory,bool> interfaceVirtualContributions
    ("VirtualContributions",
     "Switch on or off the virtual contributions.",
     &MatchboxFactory::theVirtualContributions, true, false, false);
  static SwitchOption interfaceVirtualContributionsOn
    (interfaceVirtualContributions,
     "On",
     "Switch on virtual contributions.",
     true);
  static SwitchOption interfaceVirtualContributionsOff
    (interfaceVirtualContributions,
     "Off",
     "Switch off virtual contributions.",
     false);

  static Switch<MatchboxFactory,bool> interfaceRealContributions
    ("RealContributions",
     "Switch on or off the real contributions.",
     &MatchboxFactory::theRealContributions, true, false, false);
  static SwitchOption interfaceRealContributionsOn
    (interfaceRealContributions,
     "On",
     "Switch on real contributions.",
     true);
  static SwitchOption interfaceRealContributionsOff
    (interfaceRealContributions,
     "Off",
     "Switch off real contributions.",
     false);

  static Switch<MatchboxFactory,bool> interfaceSubProcessGroups
    ("SubProcessGroups",
     "Switch on or off production of sub-process groups.",
     &MatchboxFactory::theSubProcessGroups, false, false, false);
  static SwitchOption interfaceSubProcessGroupsOn
    (interfaceSubProcessGroups,
     "On",
     "On",
     true);
  static SwitchOption interfaceSubProcessGroupsOff
    (interfaceSubProcessGroups,
     "Off",
     "Off",
     false);

  static Reference<MatchboxFactory,MatchboxPhasespace> interfacePhasespace
    ("Phasespace",
     "Set the phasespace generator.",
     &MatchboxFactory::thePhasespace, false, false, true, true, false);

  static Reference<MatchboxFactory,MatchboxScaleChoice> interfaceScaleChoice
    ("ScaleChoice",
     "Set the scale choice object.",
     &MatchboxFactory::theScaleChoice, false, false, true, true, false);

  static Parameter<MatchboxFactory,double> interfaceFactorizationScaleFactor
    ("FactorizationScaleFactor",
     "The factorization scale factor.",
     &MatchboxFactory::theFactorizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxFactory,double> interfaceRenormalizationScaleFactor
    ("RenormalizationScaleFactor",
     "The renormalization scale factor.",
     &MatchboxFactory::theRenormalizationScaleFactor, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Switch<MatchboxFactory,bool> interfaceFixedCouplings
    ("FixedCouplings",
     "Switch on or off fixed couplings.",
     &MatchboxFactory::theFixedCouplings, true, false, false);
  static SwitchOption interfaceFixedCouplingsOn
    (interfaceFixedCouplings,
     "On",
     "On",
     true);
  static SwitchOption interfaceFixedCouplingsOff
    (interfaceFixedCouplings,
     "Off",
     "Off",
     false);

  static Switch<MatchboxFactory,bool> interfaceVetoScales
    ("VetoScales",
     "Switch on or setting veto scales.",
     &MatchboxFactory::theVetoScales, true, false, false);
  static SwitchOption interfaceVetoScalesOn
    (interfaceVetoScales,
     "On",
     "On",
     true);
  static SwitchOption interfaceVetoScalesOff
    (interfaceVetoScales,
     "Off",
     "Off",
     false);

  static RefVector<MatchboxFactory,MatchboxAmplitude> interfaceAmplitudes
    ("Amplitudes",
     "The amplitude objects.",
     &MatchboxFactory::theAmplitudes, -1, false, false, true, true, false);

  static Reference<MatchboxFactory,MatchboxMECache> interfaceCache
    ("Cache",
     "Set the matrix element cache object.",
     &MatchboxFactory::theCache, false, false, true, true, false);

  static RefVector<MatchboxFactory,MatchboxMEBase> interfaceBornMEs
    ("BornMEs",
     "The Born matrix elements to be used",
     &MatchboxFactory::theBornMEs, -1, false, false, true, true, false);


  static RefVector<MatchboxFactory,MatchboxInsertionOperator> interfaceVirtuals
    ("Virtuals",
     "The virtual corrections to include",
     &MatchboxFactory::theVirtuals, -1, false, false, true, true, false);

  static RefVector<MatchboxFactory,MatchboxMEBase> interfaceRealEmissionMEs
    ("RealEmissionMEs",
     "The RealEmission matrix elements to be used",
     &MatchboxFactory::theRealEmissionMEs, -1, false, false, true, true, false);


  static RefVector<MatchboxFactory,MatchboxNLOME> interfaceBornVirtuals
    ("BornVirtualMEs",
     "The generated Born/virtual contributions",
     &MatchboxFactory::theBornVirtualMEs, -1, false, true, true, true, false);

  static RefVector<MatchboxFactory,SubtractedME> interfaceSubtractedMEs
    ("SubtractedMEs",
     "The generated Born/virtual contributions",
     &MatchboxFactory::theSubtractedMEs, -1, false, true, true, true, false);

  static Switch<MatchboxFactory,bool> interfaceVerbose
    ("Verbose",
     "Print full infomation on each evaluated phase space point.",
     &MatchboxFactory::theVerbose, false, false, false);
  static SwitchOption interfaceVerboseOn
    (interfaceVerbose,
     "On",
     "On",
     true);
  static SwitchOption interfaceVerboseOff
    (interfaceVerbose,
     "Off",
     "Off",
     false);

  static Parameter<MatchboxFactory,string> interfaceSubtractionData
    ("SubtractionData",
     "Prefix for subtraction check data.",
     &MatchboxFactory::theSubtractionData, "",
     false, false);


  static RefVector<MatchboxFactory,ParticleData> interfaceParticleGroup
    ("ParticleGroup",
     "The particle group just started.",
     &MatchboxFactory::particleGroup, -1, false, false, true, false, false);

  static Command<MatchboxFactory> interfaceStartParticleGroup
    ("StartParticleGroup",
     "Start a particle group.",
     &MatchboxFactory::startParticleGroup, false);

  static Command<MatchboxFactory> interfaceEndParticleGroup
    ("EndParticleGroup",
     "End a particle group.",
     &MatchboxFactory::endParticleGroup, false);


  static Command<MatchboxFactory> interfaceProcess
    ("Process",
     "Set the process to consider.",
     &MatchboxFactory::doProcess, false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxFactory,SubProcessHandler>
describeHerwigMatchboxFactory("Herwig::MatchboxFactory", "HwMatchbox.so");
