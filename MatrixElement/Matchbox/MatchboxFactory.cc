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
#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

#include <boost/progress.hpp>

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;
using std::ostream_iterator;

MatchboxFactory::MatchboxFactory() 
  : SubProcessHandler(), theNLight(0),
    theOrderInAlphaS(0), theOrderInAlphaEW(0),
    theBornContributions(true), theVirtualContributions(true),
    theRealContributions(true), theIndependentVirtuals(false),
    theSubProcessGroups(false), theInclusive(false),
    theFactorizationScaleFactor(1.0), theRenormalizationScaleFactor(1.0),
    theFixedCouplings(false), theFixedQEDCouplings(false), theVetoScales(false),
    theVerbose(false), theInitVerbose(false), theSubtractionData(""), theCheckPoles(false) {}

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

  if ( processData() && !me->processData() )
    me->processData(processData());

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

  if ( fixedQEDCouplings() )
    me->setFixedQEDCouplings();

  if ( cache() && !me->cache() )
    me->cache(cache());

  if ( verbose() )
    me->setVerbose();

}

struct LegIndex {

  int spin;
  int charge;
  int colour;

  int isSameAs;
  int isSameFamilyAs;

  inline bool operator==(const LegIndex& other) const {
    return 
      spin == other.spin &&
      charge == other.charge &&
      colour == other.colour &&
      isSameAs == other.isSameAs &&
      isSameFamilyAs == other.isSameFamilyAs;
  }

  inline bool operator<(const LegIndex& other) const {

    if ( spin != other.spin )
      return spin < other.spin;

    if ( charge != other.charge )
      return charge < other.charge;

    if ( colour != other.colour )
      return colour < other.colour;

    if ( isSameAs != other.isSameAs )
      return isSameAs < other.isSameAs;

    if ( isSameFamilyAs != other.isSameFamilyAs )
      return isSameFamilyAs < other.isSameFamilyAs;

    return false;

  }

};

vector<LegIndex> makeIndex(const PDVector& proc) {

  map<long,int> idMap;
  map<int,int> familyIdMap;
  int lastId = 0;
  int lastFamilyId = 0;

  vector<LegIndex> res;

  for ( PDVector::const_iterator p = proc.begin();
	p != proc.end(); ++p ) {
    int id;
    if ( idMap.find((**p).id()) != idMap.end() ) {
      id = idMap[(**p).id()];
    } else {
      id = lastId;
      idMap[(**p).id()] = lastId;
      ++lastId;
    }
    int familyId;
    if ( familyIdMap.find(SU2Helper::family(*p)) != familyIdMap.end() ) {
      familyId = familyIdMap[SU2Helper::family(*p)];
    } else {
      familyId = lastFamilyId;
      familyIdMap[SU2Helper::family(*p)] = lastFamilyId;
      ++lastFamilyId;
    }
    LegIndex idx;
    idx.spin = (**p).iSpin();
    idx.charge = (**p).iCharge();
    idx.colour = (**p).iColour();
    idx.isSameAs = id;
    idx.isSameFamilyAs = familyId;
    res.push_back(idx);
  }

  return res;

}

string pid(const vector<LegIndex>& key) {
  ostringstream res;
  for ( vector<LegIndex>::const_iterator k =
	  key.begin(); k != key.end(); ++k )
    res << k->spin << k->charge
	<< k->colour << k->isSameAs
	<< k->isSameFamilyAs;
  return res.str();
}

vector<Ptr<MatchboxMEBase>::ptr> MatchboxFactory::
makeMEs(const vector<string>& proc, unsigned int orderas) const {

  typedef vector<LegIndex> QNKey;

  generator()->log() << "determining subprocesses for ";
  copy(proc.begin(),proc.end(),ostream_iterator<string>(generator()->log()," "));
  generator()->log() << flush;

  map<Ptr<MatchboxAmplitude>::ptr,map<QNKey,vector<PDVector> > > ampProcs;
  set<PDVector> processes = makeSubProcesses(proc);

  vector<Ptr<MatchboxAmplitude>::ptr> matchAmplitudes;

  for ( vector<Ptr<MatchboxAmplitude>::ptr>::const_iterator amp
	  = amplitudes().begin(); amp != amplitudes().end(); ++amp ) {
    (**amp).orderInGs(orderas);
    (**amp).orderInGem(orderInAlphaEW());
    if ( (**amp).orderInGs() != orderas ||
	 (**amp).orderInGem() != orderInAlphaEW() ) {
      continue;
    }
    matchAmplitudes.push_back(*amp);
  }

  size_t combinations = processes.size()*matchAmplitudes.size();
  size_t procCount = 0;

  boost::progress_display * progressBar = 
    new boost::progress_display(combinations,generator()->log());

  for ( vector<Ptr<MatchboxAmplitude>::ptr>::const_iterator amp
	  = matchAmplitudes.begin(); amp != matchAmplitudes.end(); ++amp ) {
    (**amp).orderInGs(orderas);
    (**amp).orderInGem(orderInAlphaEW());
    for ( set<PDVector>::const_iterator p = processes.begin();
	  p != processes.end(); ++p ) {
      ++(*progressBar);
      if ( !(**amp).canHandle(*p) )
	continue;
      QNKey key = makeIndex(*p);
      ++procCount;
      ampProcs[*amp][key].push_back(*p);
    }
  }

  delete progressBar;
  generator()->log() << flush;

  vector<Ptr<MatchboxMEBase>::ptr> res;
  for ( map<Ptr<MatchboxAmplitude>::ptr,map<QNKey,vector<PDVector> > >::const_iterator
	  ap = ampProcs.begin(); ap != ampProcs.end(); ++ap ) {
  for ( map<QNKey,vector<PDVector> >::const_iterator m = ap->second.begin();
	  m != ap->second.end(); ++m ) {
      Ptr<MatchboxMEBase>::ptr me = ap->first->makeME(m->second);
      me->subProcesses() = m->second;
      me->amplitude(ap->first);
      string pname = "ME" + ap->first->name() + pid(m->first);
      if ( ! (generator()->preinitRegister(me,pname) ) )
	throw InitException() << "Matrix element " << pname << " already existing.";
      res.push_back(me);
    }
  }

  generator()->log() << "created " << res.size()
		     << " matrix element objects for "
		     << procCount << " subprocesses.\n";
  generator()->log() << "--------------------------------------------------------------------------------\n"
		     << flush;

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
#ifndef NDEBUG
	long checkid = (**p).id();
#endif
	*p = getParticleData((**p).id());
	assert((**p).id() == checkid);
      }

    const PDVector& partons = particleGroups()["j"];
    unsigned int nl = 0;
    for ( PDVector::const_iterator p = partons.begin();
	  p != partons.end(); ++p )
      if ( abs((**p).id()) < 6 )
	++nl;
    nLight(nl/2);

    vector<Ptr<MatchboxMEBase>::ptr> ames = makeMEs(process,orderInAlphaS());
    copy(ames.begin(),ames.end(),back_inserter(bornMEs()));

    if ( realContributions() ) {
      vector<string> rproc = process;
      if ( realEmissionProcess.empty() ) {
	rproc.push_back("j");
      } else {
	rproc = realEmissionProcess;
      }
      ames = makeMEs(rproc,orderInAlphaS()+1);
      copy(ames.begin(),ames.end(),back_inserter(realEmissionMEs()));
    }

  }

  // check if we have virtual contributions
  bool haveVirtuals = true;

  // check DR conventions of virtual contributions
  bool virtualsAreDR = false;
  bool virtualsAreCDR = false;

  // check finite term conventions of virtual contributions
  bool virtualsAreCS = false;
  bool virtualsAreBDK = false;
  bool virtualsAreExpanded = false;

  // check and prepare the Born and virtual matrix elements
  for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	  = bornMEs().begin(); born != bornMEs().end(); ++born ) {
    prepareME(*born);
    haveVirtuals &= (**born).haveOneLoop();
    if ( (**born).haveOneLoop() ) {
      virtualsAreDR |= (**born).isDR();
      virtualsAreCDR |= !(**born).isDR();
      virtualsAreCS |= (**born).isCS();
      virtualsAreBDK |= (**born).isBDK();
      virtualsAreExpanded |= (**born).isExpanded();
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
    virtualsAreBDK |= (**virt).isBDK();
    virtualsAreExpanded |= (**virt).isExpanded();
  }

  // check for consistent conventions on virtuals, if we are to include them
  if ( virtualContributions() ) {
    if ( virtualsAreDR && virtualsAreCDR ) {
      throw InitException() << "Virtual corrections use inconsistent regularization schemes.\n";
    }
    if ( (virtualsAreCS && virtualsAreBDK) ||
	 (virtualsAreCS && virtualsAreExpanded) ||
	 (virtualsAreBDK && virtualsAreExpanded) ||
	 (!virtualsAreCS && !virtualsAreBDK && !virtualsAreExpanded) ) {
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
      if ( virtualsAreBDK )
	(**virt).useBDK();
      if ( virtualsAreExpanded )
	(**virt).useExpanded();
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

  if ( bornContributions() || virtualContributions() ) {
    generator()->log() << "preparing Born "
		       << (virtualContributions() ? "and virtual" : "")
		       << " matrix elements." << flush;
  }

  if ( (bornContributions() && !virtualContributions()) || independentVirtuals() ) {
    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	    = bornMEs().begin(); born != bornMEs().end(); ++born ) {

      if ( (**born).onlyOneLoop() )
	continue;

      Ptr<MatchboxMEBase>::ptr bornme = (**born).cloneMe();
      string pname = fullName() + "/" + (**born).name();
      if ( independentVirtuals() )
	pname += ".Born";
      if ( ! (generator()->preinitRegister(bornme,pname) ) )
	throw InitException() << "Matrix element " << pname << " already existing.";
      bornme->cloneDependencies();
      MEs().push_back(bornme);

    }
  }

  if ( virtualContributions() ) {

    bornVirtualMEs().clear();

    boost::progress_display * progressBar = 
      new boost::progress_display(bornMEs().size(),generator()->log());

    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	    = bornMEs().begin(); born != bornMEs().end(); ++born ) {

      Ptr<MatchboxMEBase>::ptr nlo = (**born).cloneMe();
      string pname = fullName() + "/" + (**born).name();
      if ( !independentVirtuals() )
	pname += ".BornVirtual";
      else
	pname += ".Virtual";
      if ( ! (generator()->preinitRegister(nlo,pname) ) )
	throw InitException() << "NLO ME " << pname << " already existing.";

      nlo->virtuals().clear();

      if ( !nlo->onlyOneLoop() ) {
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
	if ( checkPoles() ) {
	  if ( !virtualsAreExpanded ) {
	    throw InitException() << "Cannot check epsilon poles if virtuals are not in `expanded' convention.\n";
	  }
	  nlo->doCheckPoles();
	}
      }

      if ( !bornContributions() || independentVirtuals() ) {
	nlo->doOneLoopNoBorn();
      } else {
	nlo->doOneLoop();
      }

      nlo->cloneDependencies();

      bornVirtualMEs().push_back(nlo);
      MEs().push_back(nlo);

      ++(*progressBar);

    }

    delete progressBar;

    generator()->log() << "--------------------------------------------------------------------------------\n"
		       << flush;

  }

  theSplittingDipoles.clear();
  set<cPDVector> bornProcs;
  if ( showerApproximation() ) 
    if ( showerApproximation()->needsSplittingGenerator() ) {
      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	      = bornMEs().begin(); born != bornMEs().end(); ++born )
	for ( MEBase::DiagramVector::const_iterator d = (**born).diagrams().begin();
	      d != (**born).diagrams().end(); ++d )
	  bornProcs.insert((**d).partons());
    }

  if ( realContributions() ) {

    generator()->log() << "preparing real emission matrix elements." << flush;

    if ( theSubtractionData != "" )
      if ( theSubtractionData[theSubtractionData.size()-1] != '/' )
	theSubtractionData += "/";

    subtractedMEs().clear();    

    boost::progress_display * progressBar = 
      new boost::progress_display(realEmissionMEs().size(),generator()->log());

    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator real
	    = realEmissionMEs().begin(); real != realEmissionMEs().end(); ++real ) {

      Ptr<SubtractedME>::ptr sub = new_ptr(SubtractedME());
      string pname = fullName() + "/" + (**real).name() + ".Real";
      if ( ! (generator()->preinitRegister(sub,pname) ) )
	throw InitException() << "Subtracted ME " << pname << " already existing.";

      sub->borns() = bornMEs();
      sub->head(*real);

      sub->allDipoles().clear();
      sub->dependent().clear();

      if ( subtractionData() != "" )
	sub->subtractionData(subtractionData());

      sub->getDipoles();

      if ( sub->dependent().empty() ) {
	// finite real contribution
	MEs().push_back(sub->head());
	Ptr<MatchboxMEBase>::ptr fme = dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(sub->head());
        finiteRealMEs().push_back(fme);
	sub->head(tMEPtr());
	continue;
      }

      if ( verbose() )
	sub->setVerbose();

      if ( subProcessGroups() )
	sub->setSubProcessGroups();

      if ( inclusive() )
	sub->setInclusive();

      if ( vetoScales() )
	sub->doVetoScales();

      subtractedMEs().push_back(sub);
      MEs().push_back(sub);

      if ( showerApproximation() ) {
	sub->showerApproximation(showerApproximation());
	Ptr<SubtractedME>::ptr subv = new_ptr(*sub);
	string vname = sub->fullName() + ".vsub";
	if ( ! (generator()->preinitRegister(subv,pname) ) )
	  throw InitException() << "Subtracted ME " << vname << " already existing.";
	sub->doRealShowerSubtraction();
	subv->doVirtualShowerSubtraction();
	subtractedMEs().push_back(subv);
	MEs().push_back(subv);
	if ( showerApproximation()->needsSplittingGenerator() )
	  for ( set<cPDVector>::const_iterator p = bornProcs.begin();
		p != bornProcs.end(); ++p ) {
	    vector<Ptr<SubtractionDipole>::ptr> sdip = sub->splitDipoles(*p);
	    set<Ptr<SubtractionDipole>::ptr>& dips = theSplittingDipoles[*p];
	    copy(sdip.begin(),sdip.end(),inserter(dips,dips.begin()));
	  }
      }

      ++(*progressBar);

    }

    delete progressBar;

    generator()->log() << "--------------------------------------------------------------------------------\n"
		       << flush;

  }

  if ( !theSplittingDipoles.empty() ) {
    map<Ptr<SubtractionDipole>::ptr,Ptr<SubtractionDipole>::ptr> cloneMap;
    for ( map<cPDVector,set<Ptr<SubtractionDipole>::ptr> >::const_iterator sd = theSplittingDipoles.begin();
	  sd != theSplittingDipoles.end(); ++sd ) {
      for ( set<Ptr<SubtractionDipole>::ptr>::const_iterator d = sd->second.begin();
	    d != sd->second.end(); ++d ) {
	cloneMap[*d] = Ptr<SubtractionDipole>::ptr();
      }
    }
    for ( map<Ptr<SubtractionDipole>::ptr,Ptr<SubtractionDipole>::ptr>::iterator cd =
	    cloneMap.begin(); cd != cloneMap.end(); ++cd ) {
      Ptr<SubtractionDipole>::ptr cloned = cd->first->cloneMe();
      string dname = cd->first->fullName() + ".splitting";
      if ( ! (generator()->preinitRegister(cloned,dname)) )
	throw InitException() << "Dipole '" << dname << "' already existing.";
      cloned->cloneDependencies();
      cloned->showerApproximation(Ptr<ShowerApproximation>::tptr());
      cloned->doSplitting();
      cd->second = cloned;
    }
    for ( map<cPDVector,set<Ptr<SubtractionDipole>::ptr> >::iterator sd = theSplittingDipoles.begin();
	  sd != theSplittingDipoles.end(); ++sd ) {
      set<Ptr<SubtractionDipole>::ptr> cloned;
      for ( set<Ptr<SubtractionDipole>::ptr>::iterator d = sd->second.begin();
	    d != sd->second.end(); ++d ) {
	cloned.insert(cloneMap[*d]);
      }
      sd->second = cloned;
    }
  }

  generator()->log() << "process setup finished.\n" << flush;

}

list<MatchboxFactory::SplittingChannel> 
MatchboxFactory::getSplittingChannels(tStdXCombPtr xcptr) const {

  if ( xcptr->lastProjector() )
    xcptr = xcptr->lastProjector();

  const StandardXComb& xc = *xcptr;

  cPDVector proc = xc.mePartonData();

  map<cPDVector,set<Ptr<SubtractionDipole>::ptr> >::const_iterator splitEntries
    = splittingDipoles().find(proc);

  list<SplittingChannel> res;
  if ( splitEntries == splittingDipoles().end() )
    return res;

  const set<Ptr<SubtractionDipole>::ptr>& splitDipoles = splitEntries->second;


  for ( set<Ptr<SubtractionDipole>::ptr>::const_iterator sd =
	  splitDipoles.begin(); sd != splitDipoles.end(); ++sd ) {

    Ptr<MatchboxMEBase>::tptr bornME = 
      const_ptr_cast<Ptr<MatchboxMEBase>::tptr>((**sd).underlyingBornME());

    SplittingChannel channel;
    channel.dipole = *sd;
    channel.bornXComb =
      bornME->makeXComb(xc.maxEnergy(),xc.particles(),xc.eventHandlerPtr(),
			const_ptr_cast<tSubHdlPtr>(xc.subProcessHandler()),
			xc.pExtractor(),xc.CKKWHandler(),
			xc.partonBins(),xc.cuts(),xc.diagrams(),xc.mirror(),
			PartonPairVec());

    vector<StdXCombPtr> realXCombs = (**sd).makeRealXCombs(channel.bornXComb);

    for ( vector<StdXCombPtr>::const_iterator rxc = realXCombs.begin();
	  rxc != realXCombs.end(); ++rxc ) {
      channel.realXComb = *rxc;
      if ( showerApproximation()->needsTildeXCombs() ) {
	channel.tildeXCombs.clear();
	assert(!channel.dipole->partnerDipoles().empty());
	for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator p =
		channel.dipole->partnerDipoles().begin();
	      p != channel.dipole->partnerDipoles().end(); ++p ) {
	  StdXCombPtr txc = channel.dipole->makeBornXComb(channel.realXComb);
	  if ( txc )
	    channel.tildeXCombs.push_back(txc);
	}
      }
      res.push_back(channel);
    }

  }

  return res;

}

void MatchboxFactory::print(ostream& os) const {

  os << "--- MatchboxFactory setup -----------------------------------------------------------\n";

  if ( !amplitudes().empty() ) {
    os << " generated Born matrix elements:\n";
    for ( vector<Ptr<MatchboxMEBase>::ptr>::const_iterator m = bornMEs().begin();
	  m != bornMEs().end(); ++m ) {
      (**m).print(os);
    }
    os << flush;
    os << " generated real emission matrix elements:\n";
    for ( vector<Ptr<MatchboxMEBase>::ptr>::const_iterator m = realEmissionMEs().begin();
	  m != realEmissionMEs().end(); ++m ) {
      (**m).print(os);
    }
    os << flush;
  }

  os << " generated Born+virtual matrix elements:\n";

  for ( vector<Ptr<MatchboxMEBase>::ptr>::const_iterator bv
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
  if ( initVerbose() )
    print(Repository::clog());
  SubProcessHandler::doinit();
}


void MatchboxFactory::persistentOutput(PersistentOStream & os) const {
  os << theDiagramGenerator << theProcessData
     << theNLight << theOrderInAlphaS << theOrderInAlphaEW 
     << theBornContributions << theVirtualContributions
     << theRealContributions << theIndependentVirtuals << theSubProcessGroups << theInclusive
     << thePhasespace << theScaleChoice
     << theFactorizationScaleFactor << theRenormalizationScaleFactor
     << theFixedCouplings << theFixedQEDCouplings << theVetoScales
     << theAmplitudes << theCache
     << theBornMEs << theVirtuals << theRealEmissionMEs
     << theBornVirtualMEs << theSubtractedMEs << theFiniteRealMEs
     << theVerbose << theInitVerbose << theSubtractionData << theCheckPoles
     << theParticleGroups << process << realEmissionProcess
     << theShowerApproximation << theSplittingDipoles;
}

void MatchboxFactory::persistentInput(PersistentIStream & is, int) {
  is >> theDiagramGenerator >> theProcessData
     >> theNLight >> theOrderInAlphaS >> theOrderInAlphaEW 
     >> theBornContributions >> theVirtualContributions
     >> theRealContributions >> theIndependentVirtuals >> theSubProcessGroups >> theInclusive
     >> thePhasespace >> theScaleChoice
     >> theFactorizationScaleFactor >> theRenormalizationScaleFactor
     >> theFixedCouplings >> theFixedQEDCouplings >> theVetoScales
     >> theAmplitudes >> theCache
     >> theBornMEs >> theVirtuals >> theRealEmissionMEs
     >> theBornVirtualMEs >> theSubtractedMEs >> theFiniteRealMEs
     >> theVerbose >> theInitVerbose >> theSubtractionData >> theCheckPoles
     >> theParticleGroups >> process >> realEmissionProcess
     >> theShowerApproximation >> theSplittingDipoles;
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

string MatchboxFactory::doSingleRealProcess(string in) {
  realEmissionProcess = StringUtils::split(in);
  if ( realEmissionProcess.size() < 3 )
    throw InitException() << "Invalid process.";
  for ( vector<string>::iterator p = realEmissionProcess.begin();
	p != realEmissionProcess.end(); ++p ) {
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

  static Reference<MatchboxFactory,ProcessData> interfaceProcessData
    ("ProcessData",
     "Set the process data object to be used.",
     &MatchboxFactory::theProcessData, false, false, true, true, false);

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

  static Switch<MatchboxFactory,bool> interfaceIndependentVirtuals
    ("IndependentVirtuals",
     "Switch on or off virtual contributions as separate subprocesses.",
     &MatchboxFactory::theIndependentVirtuals, true, false, false);
  static SwitchOption interfaceIndependentVirtualsOn
    (interfaceIndependentVirtuals,
     "On",
     "Switch on virtual contributions as separate subprocesses.",
     true);
  static SwitchOption interfaceIndependentVirtualsOff
    (interfaceIndependentVirtuals,
     "Off",
     "Switch off virtual contributions as separate subprocesses.",
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

  static Switch<MatchboxFactory,bool> interfaceInclusive
    ("Inclusive",
     "Switch on or off production of inclusive cross section.",
     &MatchboxFactory::theInclusive, false, false, false);
  static SwitchOption interfaceInclusiveOn
    (interfaceInclusive,
     "On",
     "On",
     true);
  static SwitchOption interfaceInclusiveOff
    (interfaceInclusive,
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

  static Switch<MatchboxFactory,bool> interfaceFixedQEDCouplings
    ("FixedQEDCouplings",
     "Switch on or off fixed QED couplings.",
     &MatchboxFactory::theFixedQEDCouplings, true, false, false);
  static SwitchOption interfaceFixedQEDCouplingsOn
    (interfaceFixedQEDCouplings,
     "On",
     "On",
     true);
  static SwitchOption interfaceFixedQEDCouplingsOff
    (interfaceFixedQEDCouplings,
     "Off",
     "Off",
     false);

  static Switch<MatchboxFactory,bool> interfaceVetoScales
    ("VetoScales",
     "Switch on or setting veto scales.",
     &MatchboxFactory::theVetoScales, false, false, false);
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


  static RefVector<MatchboxFactory,MatchboxMEBase> interfaceBornVirtuals
    ("BornVirtualMEs",
     "The generated Born/virtual contributions",
     &MatchboxFactory::theBornVirtualMEs, -1, false, true, true, true, false);

  static RefVector<MatchboxFactory,SubtractedME> interfaceSubtractedMEs
    ("SubtractedMEs",
     "The generated subtracted real emission contributions",
     &MatchboxFactory::theSubtractedMEs, -1, false, true, true, true, false);

  static RefVector<MatchboxFactory,MatchboxMEBase> interfaceFiniteRealMEs
    ("FiniteRealMEs",
     "The generated finite real contributions",
     &MatchboxFactory::theFiniteRealMEs, -1, false, true, true, true, false);

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

  static Switch<MatchboxFactory,bool> interfaceInitVerbose
    ("InitVerbose",
     "Print setup information.",
     &MatchboxFactory::theInitVerbose, false, false, false);
  static SwitchOption interfaceInitVerboseOn
    (interfaceInitVerbose,
     "On",
     "On",
     true);
  static SwitchOption interfaceInitVerboseOff
    (interfaceInitVerbose,
     "Off",
     "Off",
     false);

  static Parameter<MatchboxFactory,string> interfaceSubtractionData
    ("SubtractionData",
     "Prefix for subtraction check data.",
     &MatchboxFactory::theSubtractionData, "",
     false, false);

  static Switch<MatchboxFactory,bool> interfaceCheckPoles
    ("CheckPoles",
     "Switch on or off checks of epsilon poles.",
     &MatchboxFactory::theCheckPoles, true, false, false);
  static SwitchOption interfaceCheckPolesOn
    (interfaceCheckPoles,
     "On",
     "On",
     true);
  static SwitchOption interfaceCheckPolesOff
    (interfaceCheckPoles,
     "Off",
     "Off",
     false);

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

  static Command<MatchboxFactory> interfaceSingleRealProcess
    ("SingleRealProcess",
     "Set the process to consider.",
     &MatchboxFactory::doSingleRealProcess, false);

  static Reference<MatchboxFactory,ShowerApproximation> interfaceShowerApproximation
    ("ShowerApproximation",
     "Set the shower approximation to be considered.",
     &MatchboxFactory::theShowerApproximation, false, false, true, true, false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxFactory,SubProcessHandler>
describeHerwigMatchboxFactory("Herwig::MatchboxFactory", "HwMatchbox.so");
