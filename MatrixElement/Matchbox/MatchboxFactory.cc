// -*- C++ -*-
//
// MatchboxFactory.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "Herwig/API/RunDirectories.h"

#include "Herwig/Utilities/Progress.h"

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;
using std::ostream_iterator;

MatchboxFactory::MatchboxFactory() 
  : SubProcessHandler(), theNLight(0),
    theOrderInAlphaS(0), theOrderInAlphaEW(0),
    theBornContributions(true), theVirtualContributions(true),
    theRealContributions(true), theIndependentVirtuals(false),
    theIndependentPKs(false),
    theFactorizationScaleFactor(1.0), theRenormalizationScaleFactor(1.0),
    theFixedCouplings(false), theFixedQEDCouplings(false), theVetoScales(false),
    theDipoleSet(0), theVerbose(false), theInitVerbose(false), 
    theSubtractionData(""), theSubtractionPlotType(1), theSubtractionScatterPlot(false),
    thePoleData(""), theRealEmissionScales(false), theAllProcesses(false),
  theMECorrectionsOnly(false), theLoopSimCorrections(false), ranSetup(false),
  theFirstPerturbativePDF(true), theSecondPerturbativePDF(true),
  inProductionMode(false), theSpinCorrelations(false),theAlphaParameter(1.),
  theEnforceChargeConservation(true), theEnforceColourConservation(false),
  theEnforceLeptonNumberConservation(false), theEnforceQuarkNumberConservation(false),
  theLeptonFlavourDiagonal(false), theQuarkFlavourDiagonal(false) {}

MatchboxFactory::~MatchboxFactory() {}

Ptr<MatchboxFactory>::tptr MatchboxFactory::theCurrentFactory = Ptr<MatchboxFactory>::tptr();

bool& MatchboxFactory::theIsMatchboxRun() {
  static bool flag = false;
  return flag;
}

IBPtr MatchboxFactory::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxFactory::fullclone() const {
  return new_ptr(*this);
}

void MatchboxFactory::prepareME(Ptr<MatchboxMEBase>::ptr me) {

  Ptr<MatchboxAmplitude>::ptr amp =
    dynamic_ptr_cast<Ptr<MatchboxAmplitude>::ptr>((*me).amplitude());
  me->matchboxAmplitude(amp);

  if ( phasespace() && !me->phasespace() )
    me->phasespace(phasespace());

  if ( scaleChoice() && !me->scaleChoice() )
    me->scaleChoice(scaleChoice());

  if ( !reweighters().empty() ) {
    for ( vector<ReweightPtr>::const_iterator rw = reweighters().begin();
	  rw != reweighters().end(); ++rw )
      me->addReweighter(*rw);
  }

  if ( !preweighters().empty() ) {
    for ( vector<ReweightPtr>::const_iterator rw = preweighters().begin();
	  rw != preweighters().end(); ++rw )
      me->addPreweighter(*rw);
  }

}

string pid(const PDVector& key) {
  ostringstream res;
  res << "[" << key[0]->PDGName() << ","
      << key[1]->PDGName() << "->";
  for ( PDVector::const_iterator k =
	  key.begin() + 2; k != key.end(); ++k )
    res << (**k).PDGName() << (k != --key.end() ? "," : "");
  res << "]";
  return res.str();
}

vector<Ptr<MatchboxMEBase>::ptr> MatchboxFactory::
makeMEs(const vector<string>& proc, unsigned int orderas, bool virt) {

  generator()->log() << "determining subprocesses for ";
  copy(proc.begin(),proc.end(),ostream_iterator<string>(generator()->log()," "));
  generator()->log() << "\n" << flush;

  map<Ptr<MatchboxAmplitude>::ptr,set<Process> > ampProcs;
  map<Process,set<Ptr<MatchboxAmplitude>::ptr> > procAmps;
  set<PDVector> processes = makeSubProcesses(proc);

  set<PDVector> colouredProcesses;
  for ( set<PDVector>::const_iterator pr = processes.begin();
	pr != processes.end(); ++pr ) {
    for ( PDVector::const_iterator pp = pr->begin();
	  pp != pr->end(); ++pp ) {
      if ( (**pp).coloured() ) {
	colouredProcesses.insert(*pr);
	break;
      }
    }
  }
  if ( colouredProcesses.size() != processes.size() &&
      (virtualContributions() || realContributions()) ) {
      // NLO not working for non coloured legs
    throw Exception()
    << "Found processes without coloured legs.\n"
    << "We currently do not support NLO corrections for those processes.\n"
    << "Please switch to a setup for LO production."
    << Exception::runerror;
  }
  
  // detect external particles with non-zero width for the hard process
  bool trouble = false;
  string troubleMaker;
  for ( set<PDVector>::const_iterator pr = processes.begin();
	pr != processes.end(); ++pr ) {
    for ( PDVector::const_iterator pp = pr->begin();
	  pp != pr->end(); ++pp ) {
      if ( (**pp).hardProcessWidth() != ZERO ) {
	trouble = true;
	troubleMaker = (**pp).PDGName();
	break;
      }
    }
  }
  if ( trouble ) {
    throw Exception()
      << "MatchboxFactory::makeMEs(): Particle '"
      << troubleMaker << "' appears as external\nprocess leg with non-zero "
      << "width to be used in the hard process calculation.\n"
      << "Please check your setup and consider setting HardProcessWidth to zero."
      << Exception::runerror;
  }

  vector<Ptr<MatchboxAmplitude>::ptr> matchAmplitudes;

  unsigned int lowestAsOrder =
    allProcesses() ? 0 : orderas;
  unsigned int highestAsOrder = orderas;

  unsigned int lowestAeOrder =
    allProcesses() ? 0 : orderInAlphaEW();
  unsigned int highestAeOrder = orderInAlphaEW();

  for ( unsigned int oas = lowestAsOrder; oas <= highestAsOrder; ++oas ) {
    for ( unsigned int oae = lowestAeOrder; oae <= highestAeOrder; ++oae ) {
      for ( vector<Ptr<MatchboxAmplitude>::ptr>::const_iterator amp
	      = amplitudes().begin(); amp != amplitudes().end(); ++amp ) {
	if ( !theSelectedAmplitudes.empty() ) {
	  if ( find(theSelectedAmplitudes.begin(),theSelectedAmplitudes.end(),*amp)
	       == theSelectedAmplitudes.end() )
	    continue;
	}
	if ( !theDeselectedAmplitudes.empty() ) {
	  if ( find(theDeselectedAmplitudes.begin(),theDeselectedAmplitudes.end(),*amp)
	       != theDeselectedAmplitudes.end() )
	    continue;
	}
	(**amp).orderInGs(oas);
	(**amp).orderInGem(oae);
	if ( (**amp).orderInGs() != oas ||
	     (**amp).orderInGem() != oae ) {
	  continue;
	}
	matchAmplitudes.push_back(*amp);
      }
    }
  }

  size_t combinations =  processes.size()*matchAmplitudes.size();
  size_t procCount = 0;

  generator()->log() << "building matrix elements." << flush;

  progress_display progressBar{ combinations, generator()->log() };

  for ( unsigned int oas = lowestAsOrder; oas <= highestAsOrder; ++oas ) {
    for ( unsigned int oae = lowestAeOrder; oae <= highestAeOrder; ++oae ) {
      for ( vector<Ptr<MatchboxAmplitude>::ptr>::const_iterator amp
	      = matchAmplitudes.begin(); amp != matchAmplitudes.end(); ++amp ) {
	(**amp).orderInGs(oas);
	(**amp).orderInGem(oae);
	for ( set<PDVector>::const_iterator p = processes.begin();
	      p != processes.end(); ++p ) {
	  ++progressBar;
	  if ( !(**amp).canHandle(*p,this,virt) )
	    continue;
	  if ( (**amp).isExternal() )
	    externalAmplitudes().insert(*amp);
	  ++procCount;
	  Process proc(*p,oas,oae);
	  ampProcs[*amp].insert(proc);
	  procAmps[proc].insert(*amp);
	}
      }
    }
  }

  generator()->log() << flush;

  bool clash = false;
  for ( map<Process,set<Ptr<MatchboxAmplitude>::ptr> >::const_iterator check = 
	  procAmps.begin(); check != procAmps.end(); ++check ) {
    if ( check->second.size() > 1 ) {
      clash = true;
      generator()->log() << "Several different amplitudes have been found for: "
			 << check->first.legs[0]->PDGName() << " "
			 << check->first.legs[1]->PDGName() << " -> ";
      for ( PDVector::const_iterator p = check->first.legs.begin() + 2;
	    p != check->first.legs.end(); ++p )
	generator()->log() << (**p).PDGName() << " ";
      generator()->log() << "at alpha_s^" << check->first.orderInAlphaS
			 << " and alpha_ew^" << check->first.orderInAlphaEW
			 << "\n";
      generator()->log() << "The following amplitudes claim responsibility:\n";
      for ( set<Ptr<MatchboxAmplitude>::ptr>::const_iterator a = check->second.begin();
	    a != check->second.end(); ++a ) {
	generator()->log() << (**a).name() << " ";
      }
      generator()->log() << "\n";
    }
  }
  if ( clash ) {
    throw Exception() << "MatchboxFactory: Ambiguous amplitude setup - please check your input files.\n"
      << "To avoid this problem use the SelectAmplitudes or DeselectAmplitudes interfaces.\n"
		      << Exception::runerror;
  }

  bool canDoSpinCorrelations = true;

  vector<Ptr<MatchboxMEBase>::ptr> res;
  for ( map<Ptr<MatchboxAmplitude>::ptr,set<Process> >::const_iterator
	  ap = ampProcs.begin(); ap != ampProcs.end(); ++ap ) {
    canDoSpinCorrelations &= ap->first->canFillRhoMatrix();
    for ( set<Process>::const_iterator m = ap->second.begin();
	  m != ap->second.end(); ++m ) {
      Ptr<MatchboxMEBase>::ptr me = ap->first->makeME(m->legs);
      me->subProcess() = *m;
      me->amplitude(ap->first);
      me->matchboxAmplitude(ap->first);
      prepareME(me);
      string pname = "ME" + pid(m->legs);
      if ( ! (generator()->preinitRegister(me,pname) ) )
	throw Exception() << "MatchboxFactory: Matrix element " << pname << " already existing."
			  << Exception::runerror;
      if ( me->diagrams().empty() )continue;
      res.push_back(me);
      if ( theFirstPerturbativePDF )
	theIncoming.insert(m->legs[0]->id());
      if ( theSecondPerturbativePDF )
	theIncoming.insert(m->legs[1]->id());
    }
  }

  if ( spinCorrelations() && !canDoSpinCorrelations ) {
    generator()->log() << "Warning: Spin correlations have been requested, but no amplitude is "
		       << "capable of performing these.\n";
    theSpinCorrelations = false;
  }

  generator()->log() << "created "
		     << procCount << " subprocesses.\n";
  generator()->log() << "---------------------------------------------------\n"
		     << flush;

  return res;

}

int MatchboxFactory::orderOLPProcess(const Process& proc,
				     Ptr<MatchboxAmplitude>::tptr amp,
				     int type) {
  map<pair<Process,int>,int>& procs =
    olpProcesses()[amp];
  map<pair<Process,int>,int>::const_iterator it =
    procs.find(make_pair(proc,type));
  if ( it != procs.end() )
    return it->second;
  int id = procs.size();
  procs[make_pair(proc,type)] = id + 1;
  return id + 1;
}

void MatchboxFactory::productionMode() {

  if ( inProductionMode )
    return;

  if ( !bornContributions() && !virtualContributions() && !realContributions() )
    throw Exception() << "MatchboxFactory: At least one cross section contribution needs to be enabled.\n"
      << "Please check your setup.\n"
      << Exception::runerror;

  bool needTrueVirtuals =
    virtualContributions() && !meCorrectionsOnly() && !loopSimCorrections();

  for ( vector<Ptr<MatchboxAmplitude>::ptr>::iterator amp
	  = amplitudes().begin(); amp != amplitudes().end(); ++amp ) {
    if ( !needTrueVirtuals && (**amp).oneLoopAmplitude() ) {
      Repository::clog() << "One-loop contributions from '"
			 << (**amp).name()
			 << "' are not required and will be disabled.\n"
			 << flush;
      (**amp).disableOneLoop();
    }
  }

  if ( subtractionData() != "" && !subProcessGroups() ) {
    throw Exception() << "MatchboxFactory: Plain NLO settings are required for subtraction checks.\n"
      << "Please check your setup.\n"
      << Exception::runerror;
  }
  if ( showerApproximation() && !virtualContributions() && !realContributions() ) {
    Repository::clog() << "Warning: Matching requested for LO run. Matching disabled.\n" << flush;
    showerApproximation(Ptr<ShowerApproximation>::tptr());
  }
  if ( showerApproximation() && (subtractionData() != "" || subProcessGroups()) ) {
    Repository::clog() << "Warning: Matching requested for plain NLO run. Matching disabled.\n" << flush;
    showerApproximation(Ptr<ShowerApproximation>::tptr());
  }

  if ( showerApproximation() ) {
    if ( spinCorrelations() && !showerApproximation()->hasSpinCorrelations() ) {
      Repository::clog() << "Warning: Spin correlations have been requested but the matching "
			 << "object is not capable of these. Spin correlations will be turned of.\n"
			 << flush;
      theSpinCorrelations = false;
    }
  }

  inProductionMode = true;

}

void MatchboxFactory::setup() {

  useMe();

  if ( !ranSetup ) {

    if ( !inProductionMode )
      throw Exception() << "MatchboxFactory: The MatchboxFactory object '"
			<< name() << "' has not been switched to production mode.\n"
			<< "Did you use 'do "
			<< name() << ":ProductionMode' before isolating the event generator?\n"
			<< Exception::runerror;

    olpProcesses().clear();
    externalAmplitudes().clear();
    theHighestVirtualsize = 0;
    theIncoming.clear();

    bool needTrueVirtuals =
      virtualContributions() && !meCorrectionsOnly() && !loopSimCorrections();

    if ( bornMEs().empty() ) {

      if ( particleGroups().find("j") == particleGroups().end() )
	throw Exception() << "MatchboxFactory: Could not find a jet particle group named 'j'"
			  << Exception::runerror;

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
	    p != partons.end(); ++p ) {
	if ( abs((**p).id()) < 7 && (**p).hardProcessMass() == ZERO )
	  ++nl;
	if ( (**p).id() > 0 && (**p).id() < 7 && (**p).hardProcessMass() == ZERO )
	  nLightJetVec( (**p).id() );
	if ( (**p).id() > 0 && (**p).id() < 7 && (**p).hardProcessMass() != ZERO )
	  nHeavyJetVec( (**p).id() );
      }
      nLight(nl/2);

      if ( particleGroups().find("p") == particleGroups().end() )
	throw Exception() << "MatchboxFactory: Could not find a hadron particle group named 'p'"
			  << Exception::runerror;

      const PDVector& partonsInP = particleGroups()["p"];
      for ( PDVector::const_iterator pip = partonsInP.begin();
	    pip != partonsInP.end(); ++pip ) {
	if ( (**pip).id() > 0 && (**pip).id() < 7 && (**pip).hardProcessMass() == ZERO )
	  nLightProtonVec( (**pip).id() );
      }

      vector<Ptr<MatchboxMEBase>::ptr> mes;
      for ( vector<vector<string> >::const_iterator p = processes.begin();
	    p != processes.end(); ++p ) {
	if( needTrueVirtuals ) {
	  theHighestVirtualsize = max(theHighestVirtualsize,(int((*p).size())));
	}
	mes = makeMEs(*p,orderInAlphaS(),needTrueVirtuals);
	copy(mes.begin(),mes.end(),back_inserter(bornMEs()));
	if ( realContributions() || meCorrectionsOnly() ||
	      (showerApproximation() && virtualContributions()) ||
	      (showerApproximation() && loopSimCorrections()) ) {
	  if ( realEmissionProcesses.empty() ) {
	    vector<string> rproc = *p;
	    rproc.push_back("j");
	    mes = makeMEs(rproc,orderInAlphaS()+1,false);
	    copy(mes.begin(),mes.end(),back_inserter(realEmissionMEs()));
	  }
	}
      }
      if ( realContributions() || meCorrectionsOnly() ||
	    (showerApproximation() && virtualContributions()) ||
	    (showerApproximation() && loopSimCorrections()) ) {
	if ( !realEmissionProcesses.empty() ) {
	  for ( vector<vector<string> >::const_iterator q =
		  realEmissionProcesses.begin(); q != realEmissionProcesses.end(); ++q ) {
	    mes = makeMEs(*q,orderInAlphaS()+1,false);
	    copy(mes.begin(),mes.end(),back_inserter(realEmissionMEs()));
	  }
	}
      }

    }

    if ( loopInducedMEs().empty() ) {

      for ( vector<vector<string> >::const_iterator p = loopInducedProcesses.begin();
	    p != loopInducedProcesses.end(); ++p ) {
	vector<Ptr<MatchboxMEBase>::ptr> mes = makeMEs(*p,orderInAlphaS(),false);
	copy(mes.begin(),mes.end(),back_inserter(loopInducedMEs()));
      }

    }

    if( bornMEs().empty() && realEmissionMEs().empty() && loopInducedMEs().empty() )
      throw Exception() << "MatchboxFactory: No matrix elements have been found.\n\
      Please check if your order of Alpha_s and Alpha_ew have the right value.\n"
			<< Exception::runerror;

    // check if we have virtual contributions
    bool haveVirtuals = true;

    // check DR conventions of virtual contributions
    bool virtualsAreDR = false;
    bool virtualsAreCDR = false;

    // check finite term conventions of virtual contributions
    bool virtualsAreCS = false;
    bool virtualsAreBDK = false;
    bool virtualsAreExpanded = false;

    // renormalization scheme
    bool virtualsAreDRbar = false;

    // check and prepare the Born and virtual matrix elements
    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	    = bornMEs().begin(); born != bornMEs().end(); ++born ) {
      prepareME(*born);
      haveVirtuals &= (**born).haveOneLoop();
      if ( needTrueVirtuals ) {
	if ( (**born).haveOneLoop() ) {
	  virtualsAreDRbar |= (**born).isDRbar();
	  virtualsAreDR |= (**born).isDR();
	  virtualsAreCDR |= !(**born).isDR();
	  virtualsAreCS |= (**born).isCS();
	  virtualsAreBDK |= (**born).isBDK();
	  virtualsAreExpanded |= (**born).isExpanded();
	}
      }
    }

    // prepare the loop induced matrix elements
    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator looped
	    = loopInducedMEs().begin(); looped != loopInducedMEs().end(); ++looped ) {
      prepareME(*looped);
    }

    if ( needTrueVirtuals ) {

      // check the additional insertion operators
      if ( !virtuals().empty() )
	haveVirtuals = true;    
      for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
	      = virtuals().begin(); virt != virtuals().end(); ++virt ) {
	virtualsAreDRbar |= (**virt).isDRbar();
	virtualsAreDR |= (**virt).isDR();
	virtualsAreCDR |= !(**virt).isDR();
	virtualsAreCS |= (**virt).isCS();
	virtualsAreBDK |= (**virt).isBDK();
	virtualsAreExpanded |= (**virt).isExpanded();
      }

      // check for consistent conventions on virtuals, if we are to include them
      if ( virtualContributions() ) {
	if ( !haveVirtuals ) {
	  throw Exception() << "MatchboxFactory: Could not find amplitudes for all virtual contributions needed.\n"
			    << Exception::runerror;
	}
	if ( virtualsAreDR && virtualsAreCDR ) {
	  throw Exception() << "MatchboxFactory: Virtual corrections use inconsistent regularization schemes.\n"
			    << Exception::runerror;
	}
	if ( (virtualsAreCS && virtualsAreBDK) ||
	     (virtualsAreCS && virtualsAreExpanded) ||
	     (virtualsAreBDK && virtualsAreExpanded) ||
	     (!virtualsAreCS && !virtualsAreBDK && !virtualsAreExpanded) ) {
	  throw Exception() << "MatchboxFactory: Virtual corrections use inconsistent conventions on finite terms.\n"
			    << Exception::runerror;
	}
      }

    }

    // prepare the real emission matrix elements
    if ( realContributions() || meCorrectionsOnly() ||
	 (showerApproximation() && virtualContributions()) ||
	 (showerApproximation() && loopSimCorrections()) ) {
      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator real
	      = realEmissionMEs().begin(); real != realEmissionMEs().end(); ++real ) {
	prepareME(*real);
      }
    }

    // start creating matrix elements
    MEs().clear();

    // setup born and virtual contributions

    if ( bornContributions() || virtualContributions() ) {
      generator()->log() << "preparing Born"
			 << (virtualContributions() ? " and virtual" : "")
			 << " matrix elements.\n" << flush;
    }

    if ( (bornContributions() && !virtualContributions()) || 
	 (bornContributions() && meCorrectionsOnly()) || 
	 (bornContributions() && virtualContributions() && independentVirtuals()) ) {
      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	      = bornMEs().begin(); born != bornMEs().end(); ++born ) {

	if ( (**born).onlyOneLoop() )
	  continue;

	Ptr<MatchboxMEBase>::ptr bornme = (**born).cloneMe();
	string pname = fullName() + "/" + (**born).name();
	if ( virtualContributions() && independentVirtuals() )
	  pname += ".Born";
	if ( ! (generator()->preinitRegister(bornme,pname) ) )
	  throw Exception() << "MatchboxFactory: Matrix element " << pname << " already existing."
			    << Exception::runerror;

	if ( bornme->isOLPTree() ) {
	  int id = orderOLPProcess(bornme->subProcess(),
				   (**born).matchboxAmplitude(),
				   ProcessType::treeME2);
	  bornme->olpProcess(ProcessType::treeME2,id);
	}

	bornme->needsNoCorrelations();

	bornme->cloneDependencies();
	MEs().push_back(bornme);

      }

    }

    if ( bornContributions() && !loopInducedMEs().empty() ) {

      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator looped
	      = loopInducedMEs().begin(); looped != loopInducedMEs().end(); ++looped ) {

	Ptr<MatchboxMEBase>::ptr loopme = (**looped).cloneMe();
	string pname = fullName() + "/" + (**looped).name() + ".LoopInduced";
	if ( ! (generator()->preinitRegister(loopme,pname) ) )
	  throw Exception() << "MatchboxFactory: Matrix element " << pname << " already existing."
			    << Exception::runerror;

	if ( loopme->isOLPTree() ) {
	  int id = orderOLPProcess(loopme->subProcess(),
				   (**looped).matchboxAmplitude(),
				   ProcessType::loopInducedME2);
	  loopme->olpProcess(ProcessType::loopInducedME2,id);
	}

	loopme->needsNoCorrelations();

	loopme->cloneDependencies();
	MEs().push_back(loopme);

      }

    }

    if ( needTrueVirtuals ) {

      bornVirtualMEs().clear();

      progress_display progressBar{ bornMEs().size(), generator()->log() };

      if ( thePoleData != "" )
	if ( thePoleData[thePoleData.size()-1] != '/' )
	  thePoleData += "/";

      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	      = bornMEs().begin(); born != bornMEs().end(); ++born ) {

	Ptr<MatchboxMEBase>::ptr nlo = (**born).cloneMe();
	string pname = fullName() + "/" + (**born).name();
	if ( !independentVirtuals() && !(!bornContributions() && virtualContributions()) )
	  pname += ".BornVirtual";
	else if ( independentPKs() && !nlo->onlyOneLoop() )
	  pname += ".VirtualVI";
        else 
	  pname += ".Virtual";
	if ( ! (generator()->preinitRegister(nlo,pname) ) )
	  throw Exception() << "MatchboxFactory: NLO ME " << pname << " already existing."
			    << Exception::runerror;

	nlo->virtuals().clear();

	if ( !nlo->onlyOneLoop() ) {
	  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
		  = virtuals().begin(); virt != virtuals().end(); ++virt ) {
	    if ( (**virt).apply((**born).diagrams().front()->partons()) )
	      nlo->virtuals().push_back(*virt);
	  }
	  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
		  = DipoleRepository::insertionIOperators(dipoleSet()).begin(); 
		virt != DipoleRepository::insertionIOperators(dipoleSet()).end(); ++virt ) {
	    if ( (**virt).apply((**born).diagrams().front()->partons()) )
	      nlo->virtuals().push_back(*virt);
	  }
          if ( !independentVirtuals() || ( independentVirtuals() && !independentPKs() ) ) { 
	    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
	            = DipoleRepository::insertionPKOperators(dipoleSet()).begin(); 
	          virt != DipoleRepository::insertionPKOperators(dipoleSet()).end(); ++virt ) {
	      if ( (**virt).apply((**born).diagrams().front()->partons()) )
	        nlo->virtuals().push_back(*virt);
	    }
          }
	  if ( nlo->virtuals().empty() )
	    throw Exception() << "MatchboxFactory: No insertion operators have been found for "
				  << (**born).name() << "\n"
			      << Exception::runerror;
	  if ( checkPoles() ) {
	    if ( !virtualsAreExpanded ) {
	      throw Exception() 
		<< "MatchboxFactory: Cannot check epsilon poles if virtuals are not in `expanded' convention.\n"
		<< Exception::runerror;
	    }
	  }
	}

	if ( !bornContributions() || independentVirtuals() ) {
	  nlo->doOneLoopNoBorn();
	} else {
	  nlo->doOneLoop();
	}

	if ( nlo->isOLPLoop() ) {
	  int id = orderOLPProcess(nlo->subProcess(),
				   (**born).matchboxAmplitude(),
				   ProcessType::oneLoopInterference);
	  nlo->olpProcess(ProcessType::oneLoopInterference,id);
	  if ( !nlo->onlyOneLoop() && nlo->needsOLPCorrelators() ) {
	    id = orderOLPProcess(nlo->subProcess(),
				 (**born).matchboxAmplitude(),
				 ProcessType::colourCorrelatedME2);
	    nlo->olpProcess(ProcessType::colourCorrelatedME2,id);	  
	  }
	}

	nlo->needsCorrelations();

	nlo->cloneDependencies();

	bornVirtualMEs().push_back(nlo);
	MEs().push_back(nlo);

	if ( independentVirtuals() && independentPKs() && !nlo->onlyOneLoop() ) {

	  Ptr<MatchboxMEBase>::ptr nlopk = (**born).cloneMe();
	  string pnamepk = fullName() + "/" + (**born).name();
	  pnamepk += ".VirtualPK";
	  if ( ! (generator()->preinitRegister(nlopk,pnamepk) ) )
	    throw Exception() << "MatchboxFactory: NLO ME " << pnamepk << " already existing."
			      << Exception::runerror;

	  nlopk->virtuals().clear();

	  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator virt
	          = DipoleRepository::insertionPKOperators(dipoleSet()).begin(); 
	        virt != DipoleRepository::insertionPKOperators(dipoleSet()).end(); ++virt ) {
	    if ( (**virt).apply((**born).diagrams().front()->partons()) )
	      nlopk->virtuals().push_back(*virt);
	  }

	  if ( !nlopk->virtuals().empty() ) {

	    nlopk->doOneLoopNoBorn();
	    nlopk->doOneLoopNoLoops();

	    if ( nlopk->isOLPLoop() ) {
	      int id = orderOLPProcess(nlopk->subProcess(),
	        		             (**born).matchboxAmplitude(),
	        		             ProcessType::treeME2);
	      nlopk->olpProcess(ProcessType::treeME2,id);	  
	      if ( nlopk->needsOLPCorrelators() ) {
	        id = orderOLPProcess(nlopk->subProcess(),
	    			   (**born).matchboxAmplitude(),
	    			   ProcessType::colourCorrelatedME2);
	        nlopk->olpProcess(ProcessType::colourCorrelatedME2,id);	  
              }
	    }

	    nlopk->needsCorrelations();

	    nlopk->cloneDependencies();

	    bornVirtualMEs().push_back(nlopk);
	    MEs().push_back(nlopk);

          }

        }

	++progressBar;

      }

      generator()->log() << "---------------------------------------------------\n"
			 << flush;

    }

    theSplittingDipoles.clear();
    set<cPDVector> bornProcs;
    if ( showerApproximation() ) {
      if ( showerApproximation()->needsSplittingGenerator() ) {
	for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
		= bornMEs().begin(); born != bornMEs().end(); ++born )
	  for ( MEBase::DiagramVector::const_iterator d = (**born).diagrams().begin();
		d != (**born).diagrams().end(); ++d )
	    bornProcs.insert((**d).partons());
      }
    }

    if ( realContributions() || meCorrectionsOnly() ||
	 (showerApproximation() && virtualContributions()) ||
	 (showerApproximation() && loopSimCorrections()) ) {

      generator()->log() << "preparing subtracted matrix elements.\n" << flush;

      if ( theSubtractionData != "" )
	if ( theSubtractionData[theSubtractionData.size()-1] != '/' )
	  theSubtractionData += "/";

      subtractedMEs().clear();    

      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	      = bornMEs().begin(); born != bornMEs().end(); ++born ) {

	if ( (**born).onlyOneLoop() )
	  continue;

	(**born).needsCorrelations();

	if ( (**born).isOLPTree() ) {
	  int id = orderOLPProcess((**born).subProcess(),
				   (**born).matchboxAmplitude(),
				   ProcessType::colourCorrelatedME2);
	  (**born).olpProcess(ProcessType::colourCorrelatedME2,id);
	  bool haveGluon = false;
	  for ( PDVector::const_iterator p = (**born).subProcess().legs.begin();
		p != (**born).subProcess().legs.end(); ++p )
	    if ( (**p).id() == 21 ) {
	      haveGluon = true;
	      break;
	    }
	  if ( haveGluon ) {
	    id = orderOLPProcess((**born).subProcess(),
				 (**born).matchboxAmplitude(),
				 ProcessType::spinColourCorrelatedME2);
	    (**born).olpProcess(ProcessType::spinColourCorrelatedME2,id);
	  }
          if ( showerApproximation() ) {
	    id = orderOLPProcess((**born).subProcess(),
				 (**born).matchboxAmplitude(),
				 ProcessType::treeME2);
	    (**born).olpProcess(ProcessType::treeME2,id);
          }
	}

      }

      progress_display progressBar{ realEmissionMEs().size(), generator()->log() };

      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator real
	      = realEmissionMEs().begin(); real != realEmissionMEs().end(); ++real ) {

	Ptr<SubtractedME>::ptr sub = new_ptr(SubtractedME());
	string pname = fullName() + "/" + (**real).name() + ".SubtractedReal";
	if ( ! (generator()->preinitRegister(sub,pname) ) )
	  throw Exception() << "MatchboxFactory: Subtracted ME " << pname << " already existing."
			    << Exception::runerror;

	(**real).needsNoCorrelations();

	if ( (**real).isOLPTree() ) {
	  int id = orderOLPProcess((**real).subProcess(),
				   (**real).matchboxAmplitude(),
				   ProcessType::treeME2);
	  (**real).olpProcess(ProcessType::treeME2,id);
	}

	sub->head(*real);

	sub->dependent().clear();

	sub->getDipoles();

	if ( sub->dependent().empty() ) {
	  // finite real contribution
	  if ( realContributions() ) {
	    Ptr<MatchboxMEBase>::ptr fme = 
	      dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(sub->head())->cloneMe();
	    string qname = fullName() + "/" + (**real).name() + ".FiniteReal";
	    if ( ! (generator()->preinitRegister(fme,qname) ) )
	      throw Exception() << "MatchboxFactory: ME " << qname << " already existing."
				<< Exception::runerror;
	    MEs().push_back(fme);	
	    finiteRealMEs().push_back(fme);
	  }
	  sub->head(tMEPtr());
	  ++progressBar;
	  continue;
	}

	if ( realEmissionScales() )
	  sub->doRealEmissionScales();

	subtractedMEs().push_back(sub);
	if ( realContributions() )
	  if ( !showerApproximation() || (showerApproximation() && showerApproximation()->hasHEvents()) )
	    MEs().push_back(sub);

	if ( showerApproximation() ) {
	  if ( virtualContributions() && !meCorrectionsOnly() && !loopSimCorrections() ) {
	    Ptr<SubtractedME>::ptr subv = new_ptr(*sub);
	    string vname = sub->fullName() + ".SubtractionIntegral";
	    if ( ! (generator()->preinitRegister(subv,vname) ) )
	      throw Exception() << "MatchboxFactory: Subtracted ME " << vname << " already existing."
				<< Exception::runerror;
	    subv->cloneDependencies(vname);
	    subv->doVirtualShowerSubtraction();
	    subtractedMEs().push_back(subv);
	    MEs().push_back(subv);
	  }
	  if ( loopSimCorrections() ) {
	    Ptr<SubtractedME>::ptr subv = new_ptr(*sub);
	    string vname = sub->fullName() + ".SubtractionIntegral";
	    if ( ! (generator()->preinitRegister(subv,vname) ) )
	      throw Exception() << "MatchboxFactory: Subtracted ME " << vname << " already existing."
				<< Exception::runerror;
	    subv->cloneDependencies(vname);
	    subv->doLoopSimSubtraction();
	    subtractedMEs().push_back(subv);
	    MEs().push_back(subv);
	  }
	  sub->doRealShowerSubtraction();
	  if ( showerApproximation()->needsSplittingGenerator() )
	    for ( set<cPDVector>::const_iterator p = bornProcs.begin();
		  p != bornProcs.end(); ++p ) {
	      vector<Ptr<SubtractionDipole>::ptr> sdip = sub->splitDipoles(*p);
	      set<Ptr<SubtractionDipole>::ptr>& dips = theSplittingDipoles[*p];
	      copy(sdip.begin(),sdip.end(),inserter(dips,dips.begin()));
	    }
	}

	++progressBar;

      }

      generator()->log() << "---------------------------------------------------\n"
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
	  throw Exception() << "MatchboxFactory: Dipole '" << dname << "' already existing."
			    << Exception::runerror;
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

    if ( !externalAmplitudes().empty() ) {
      generator()->log() << "Initializing external amplitudes.\n" << flush;
      for ( set<Ptr<MatchboxAmplitude>::tptr>::const_iterator ext =
	      externalAmplitudes().begin(); ext != externalAmplitudes().end(); ++ext ) {
	if ( !(**ext).initializeExternal() ) {
	  throw Exception()  << "Failed to initialize amplitude '" << (**ext).name() << "'\n"
			     << Exception::runerror;
	}
      }
      generator()->log() << "---------------------------------------------------\n"
			 << flush;
    }

    if ( !olpProcesses().empty() ) {
      generator()->log() << "Initializing one-loop provider(s).\n" << flush;
      map<Ptr<MatchboxAmplitude>::tptr,map<pair<Process,int>,int> > olps;
      for ( map<Ptr<MatchboxAmplitude>::tptr,map<pair<Process,int>,int> >::const_iterator
	      oit = olpProcesses().begin(); oit != olpProcesses().end(); ++oit ) {
	olps[oit->first] = oit->second;
      }
      for ( map<Ptr<MatchboxAmplitude>::tptr,map<pair<Process,int>,int> >::const_iterator
	      olpit = olps.begin(); olpit != olps.end(); ++olpit ) {
	if ( !olpit->first->startOLP(olpit->second) ) {
	  throw Exception() << "MatchboxFactory: Failed to start OLP for amplitude '" << olpit->first->name() << "'\n"
			    << Exception::runerror;
	}
      }
      generator()->log() << "---------------------------------------------------\n"
			 << flush;
    }

    generator()->log() << "Process setup finished.\n" << flush;

    ranSetup = true;

  }

}

void MatchboxFactory::SplittingChannel::print(ostream& os) const {

  os << "--- SplittingChannel setup -----------------------------------------------------\n";

  os << " Born process ";
  const StandardXComb& bxc = *bornXComb;
  os << bxc.mePartonData()[0]->PDGName() << " "
     << bxc.mePartonData()[1]->PDGName() << " -> ";
  for ( cPDVector::const_iterator p = bxc.mePartonData().begin() + 2;
	p != bxc.mePartonData().end(); ++p ) {
    os << (**p).PDGName() << " ";
  }
  os << "\n";

  os << " to real emission process ";
  const StandardXComb& rxc = *realXComb;
  os << rxc.mePartonData()[0]->PDGName() << " "
     << rxc.mePartonData()[1]->PDGName() << " -> ";
  for ( cPDVector::const_iterator p = rxc.mePartonData().begin() + 2;
	p != rxc.mePartonData().end(); ++p ) {
    os << (**p).PDGName() << " ";
  }
  os << "\n";

  os << " with dipole:\n";
  dipole->print(os);

  os << "---------------------------------------------------\n";

  os << flush;

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

  SplittingChannel channel;
  if ( !splitDipoles.empty() ) {
    Ptr<MatchboxMEBase>::tptr bornME = 
      const_ptr_cast<Ptr<MatchboxMEBase>::tptr>((**splitDipoles.begin()).underlyingBornME());
    channel.bornXComb =
      bornME->makeXComb(xc.maxEnergy(),xc.particles(),xc.eventHandlerPtr(),
			const_ptr_cast<tSubHdlPtr>(xc.subProcessHandler()),
			xc.pExtractor(),xc.CKKWHandler(),
			xc.partonBins(),xc.cuts(),xc.diagrams(),xc.mirror(),
			PartonPairVec());
  }

  for ( set<Ptr<SubtractionDipole>::ptr>::const_iterator sd =
	  splitDipoles.begin(); sd != splitDipoles.end(); ++sd ) {

    channel.dipole = *sd;

    vector<StdXCombPtr> realXCombs = (**sd).makeRealXCombs(channel.bornXComb);

    for ( vector<StdXCombPtr>::const_iterator rxc = realXCombs.begin();
	  rxc != realXCombs.end(); ++rxc ) {
      channel.realXComb = *rxc;
      if ( showerApproximation()->needsTildeXCombs() ) {
	channel.tildeXCombs.clear();
	assert(!channel.dipole->partnerDipoles().empty());
	for ( vector<Ptr<SubtractionDipole>::tptr>::const_iterator p =
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

  if ( initVerbose() ) {

    generator()->log()
      << "--- MatchboxFactory splitting channels ----------------------------------------------\n";

    const StandardXComb& bxc = *xcptr;

    generator()->log() << " hard process handled is: ";
    generator()->log() << bxc.mePartonData()[0]->PDGName() << " "
		       << bxc.mePartonData()[1]->PDGName() << " -> ";
    for ( cPDVector::const_iterator p = bxc.mePartonData().begin() + 2;
	  p != bxc.mePartonData().end(); ++p ) {
      generator()->log() << (**p).PDGName() << " ";
    }
    generator()->log() << "\n";

    for ( list<MatchboxFactory::SplittingChannel>::const_iterator sp =
	    res.begin(); sp != res.end(); ++sp ) {
      sp->print(generator()->log());
    }

    generator()->log()
      << "--------------------------------------------------------\n"
      << flush;

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

  os << "---------------------------------------------------\n";

  os << flush;

}

void MatchboxFactory::summary(ostream& os) const {
  os << "\n\n================================================================================\n"
     << " Matchbox hard process summary\n"
     << "================================================================================\n\n";

  os << " Electro-weak parameter summary:\n"
     << "---------------------------------------------------\n\n";

  os << " Electro-weak scheme : ";
  switch ( SM().ewScheme() ) {

  case 0: os << "Default"; break;
  case 1: os << "GMuScheme"; break;
  case 2: os << "alphaMZScheme"; break;
  case 3: os << "NoMass"; break;
  case 4: os << "mW"; break;
  case 5: os << "mZ"; break;
  case 6: os << "Independent"; break;
  case 7: os << "FeynRulesUFO"; break;
  default: assert(false);

  }

  os << "\n";

  os << " alphaEM is "
     << (SM().ewScheme() == 0 && !theFixedQEDCouplings ? "running" : "fixed at alphaEM(m(Z))") << "\n";

  if ( SM().ewScheme() == 0 && !theFixedQEDCouplings )
    os << " alphaEM is running at " << SM().alphaEMPtr()->nloops()
       << " loops\n\n";
  else
    os << "\n";

  os << (SM().ewScheme() != 0 ? " Tree level relations " : " Best values ")
     << "yield:\n\n"
     << " m(Z)/GeV       = "
     << getParticleData(ParticleID::Z0)->hardProcessMass()/GeV
     << "\n"
     << " g(Z)/GeV       = "
     << getParticleData(ParticleID::Z0)->hardProcessWidth()/GeV
     << "\n"
     << " m(W)/GeV       = "
     << getParticleData(ParticleID::Wplus)->hardProcessMass()/GeV
     << "\n"
     << " g(W)/GeV       = "
     << getParticleData(ParticleID::Wplus)->hardProcessWidth()/GeV
     << "\n"
     << " m(H)/GeV       = "
     << getParticleData(ParticleID::h0)->hardProcessMass()/GeV
     << "\n"
     << " g(H)/GeV       = "
     << getParticleData(ParticleID::h0)->hardProcessWidth()/GeV
     << "\n"
     << " alphaEM(m(Z))  = "
     << SM().alphaEMME(sqr(getParticleData(ParticleID::Z0)->hardProcessMass())) << "\n"
     << " sin^2(theta)   = " << SM().sin2ThetaW()
     << "\n"
     << " GeV^2 GF       = " << GeV2*SM().fermiConstant()
     << "\n\n";

  os << " Quark masses and widths are:\n"
     << "---------------------------------------------------\n\n"
     << " m(u)/GeV       = " << getParticleData(ParticleID::u)->hardProcessMass()/GeV << "\n"
     << " m(d)/GeV       = " << getParticleData(ParticleID::d)->hardProcessMass()/GeV << "\n"
     << " m(c)/GeV       = " << getParticleData(ParticleID::c)->hardProcessMass()/GeV << "\n"
     << " m(s)/GeV       = " << getParticleData(ParticleID::s)->hardProcessMass()/GeV << "\n"
     << " m(t)/GeV       = " << getParticleData(ParticleID::t)->hardProcessMass()/GeV << "\n"
     << " g(t)/GeV       = " << getParticleData(ParticleID::t)->hardProcessWidth()/GeV << "\n"
     << " m(b)/GeV       = " << getParticleData(ParticleID::b)->hardProcessMass()/GeV << "\n\n";

  os << " Lepton masses and widths are:\n"
     << "---------------------------------------------------\n\n"
     << " m(n_e)/GeV     = " << getParticleData(ParticleID::nu_e)->hardProcessMass()/GeV << "\n"
     << " m(e)/GeV       = " << getParticleData(ParticleID::eminus)->hardProcessMass()/GeV << "\n"
     << " m(n_mu)/GeV    = " << getParticleData(ParticleID::nu_mu)->hardProcessMass()/GeV << "\n"
     << " m(mu)/GeV      = " << getParticleData(ParticleID::muminus)->hardProcessMass()/GeV << "\n"
     << " m(nu_tau)/GeV  = " << getParticleData(ParticleID::nu_tau)->hardProcessMass()/GeV << "\n"
     << " m(tau)/GeV     = " << getParticleData(ParticleID::tauminus)->hardProcessMass()/GeV << "\n\n";


  os << " Strong coupling summary:\n"
     << "---------------------------------------------------\n\n";

  os << " alphaS is";
  if ( !theFixedCouplings ) {
    os << " running at " << SM().alphaSPtr()->nloops()
       << " loops with\n"
       << " alphaS(m(Z))   = " << SM().alphaSPtr()->value(sqr(getParticleData(ParticleID::Z0)->mass()))
       << "\n\n";
  } else {
    os << " fixed at "
       << SM().alphaS()
       << "\n\n";
  }

  if ( !theFixedCouplings ) {
    os << " flavour thresholds are matched at\n";
    for ( long id = 1; id <= 6; ++id ) {
      os << " m(" << id << ")/GeV       = " 
	 << (SM().alphaSPtr()->quarkMasses().empty() ?
	     getParticleData(id)->mass()/GeV : 
	     SM().alphaSPtr()->quarkMasses()[id-1]/GeV)
	 << "\n";
    }
  }

  os << "\n\n" << flush;

}


void MatchboxFactory::doinit() {
  theIsMatchboxRun() = true;
  theCurrentFactory = Ptr<MatchboxFactory>::tptr(this);
  if ( RunDirectories::empty() )
    RunDirectories::pushRunId(generator()->runName());
  setup();
  if ( theShowerApproximation )
    theShowerApproximation->init();
  if ( initVerbose() && !ranSetup )
    print(Repository::clog());
  Ptr<StandardEventHandler>::tptr eh =
    dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(generator()->eventHandler());
  assert(eh);
  if ( initVerbose() && !ranSetup ) {
    assert(standardModel());
    standardModel()->init();
    summary(Repository::clog());
  }
  SubProcessHandler::doinit();
}

void MatchboxFactory::doinitrun() {
  theIsMatchboxRun() = true;
  theCurrentFactory = Ptr<MatchboxFactory>::tptr(this);
  if ( theShowerApproximation )
    theShowerApproximation->initrun();
  Ptr<StandardEventHandler>::tptr eh =
    dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(generator()->eventHandler());
  assert(eh);
  SubProcessHandler::doinitrun();
}

const string& MatchboxFactory::buildStorage() {
  return RunDirectories::buildStorage();
}

const string& MatchboxFactory::runStorage() {
  return RunDirectories::runStorage();
}


void MatchboxFactory::persistentOutput(PersistentOStream & os) const {
  os << theDiagramGenerator << theProcessData
     << theNLight 
     << theNLightJetVec << theNHeavyJetVec << theNLightProtonVec 
     << theOrderInAlphaS << theOrderInAlphaEW 
     << theBornContributions << theVirtualContributions
     << theRealContributions << theIndependentVirtuals << theIndependentPKs
     << thePhasespace << theScaleChoice
     << theFactorizationScaleFactor << theRenormalizationScaleFactor
     << theFixedCouplings << theFixedQEDCouplings << theVetoScales
     << theAmplitudes
     << theBornMEs << theVirtuals << theRealEmissionMEs << theLoopInducedMEs
     << theBornVirtualMEs << theSubtractedMEs << theFiniteRealMEs
     << theVerbose << theInitVerbose << theSubtractionData << theSubtractionPlotType
     << theSubtractionScatterPlot << thePoleData
     << theParticleGroups << processes << loopInducedProcesses << realEmissionProcesses
     << theShowerApproximation << theSplittingDipoles
     << theRealEmissionScales << theAllProcesses
     << theOLPProcesses << theExternalAmplitudes
     << theSelectedAmplitudes << theDeselectedAmplitudes
     << theDipoleSet << theReweighters << thePreweighters
     << theMECorrectionsOnly<< theLoopSimCorrections<<theHighestVirtualsize << ranSetup
     << theIncoming << theFirstPerturbativePDF << theSecondPerturbativePDF 
     << inProductionMode << theSpinCorrelations << theAlphaParameter
     << theEnforceChargeConservation << theEnforceColourConservation
     << theEnforceLeptonNumberConservation << theEnforceQuarkNumberConservation
     << theLeptonFlavourDiagonal << theQuarkFlavourDiagonal;
}

void MatchboxFactory::persistentInput(PersistentIStream & is, int) {
  is >> theDiagramGenerator >> theProcessData
     >> theNLight 
     >> theNLightJetVec >> theNHeavyJetVec >> theNLightProtonVec 
     >> theOrderInAlphaS >> theOrderInAlphaEW 
     >> theBornContributions >> theVirtualContributions
     >> theRealContributions >> theIndependentVirtuals >> theIndependentPKs
     >> thePhasespace >> theScaleChoice
     >> theFactorizationScaleFactor >> theRenormalizationScaleFactor
     >> theFixedCouplings >> theFixedQEDCouplings >> theVetoScales
     >> theAmplitudes
     >> theBornMEs >> theVirtuals >> theRealEmissionMEs >> theLoopInducedMEs
     >> theBornVirtualMEs >> theSubtractedMEs >> theFiniteRealMEs
     >> theVerbose >> theInitVerbose >> theSubtractionData >> theSubtractionPlotType
     >> theSubtractionScatterPlot >> thePoleData
     >> theParticleGroups >> processes >> loopInducedProcesses >> realEmissionProcesses
     >> theShowerApproximation >> theSplittingDipoles
     >> theRealEmissionScales >> theAllProcesses
     >> theOLPProcesses >> theExternalAmplitudes
     >> theSelectedAmplitudes >> theDeselectedAmplitudes
     >> theDipoleSet >> theReweighters >> thePreweighters
     >> theMECorrectionsOnly>> theLoopSimCorrections>>theHighestVirtualsize >> ranSetup
     >> theIncoming >> theFirstPerturbativePDF >> theSecondPerturbativePDF
     >> inProductionMode >> theSpinCorrelations >> theAlphaParameter
     >> theEnforceChargeConservation >> theEnforceColourConservation
     >> theEnforceLeptonNumberConservation >> theEnforceQuarkNumberConservation
     >> theLeptonFlavourDiagonal >> theQuarkFlavourDiagonal;
}

string MatchboxFactory::startParticleGroup(string name) {
  particleGroupName = StringUtils::stripws(name);
  particleGroup.clear();
  return "";
}

string MatchboxFactory::endParticleGroup(string) {
  if ( particleGroup.empty() )
    throw Exception() << "MatchboxFactory: Empty particle group."
		      << Exception::runerror;
  particleGroups()[particleGroupName] = particleGroup;
  particleGroup.clear();
  return "";
}

vector<string> MatchboxFactory::parseProcess(string in) {
  vector<string> process = StringUtils::split(in);
  if ( process.size() < 3 )
    throw Exception() << "MatchboxFactory: Invalid process."
		      << Exception::runerror;
  for ( vector<string>::iterator p = process.begin();
	p != process.end(); ++p ) {
    *p = StringUtils::stripws(*p);
  }
  vector<string> pprocess;
  for ( vector<string>::const_iterator p = process.begin();
	p != process.end(); ++p ) {
    if ( *p == "->" )
      continue;
    pprocess.push_back(*p);
  }
  return pprocess;
}


string MatchboxFactory::doProcess(string in) {
  processes.push_back(parseProcess(in));
  return "";
}

string MatchboxFactory::doLoopInducedProcess(string in) {
  loopInducedProcesses.push_back(parseProcess(in));
  return "";
}

string MatchboxFactory::doSingleRealProcess(string in) {
  realEmissionProcesses.push_back(parseProcess(in));
  return "";
}

struct SortPID {
  inline bool operator()(PDPtr a, PDPtr b) const {
    return a->id() < b->id();
  }
};

//
// @TODO
//
// SP: After improving this for standard model process building this should
// actually got into a separate process builder class or something along these
// lines to have it better factored for use with BSM models.
//
//

set<PDVector> MatchboxFactory::
makeSubProcesses(const vector<string>& proc) const {

  if ( proc.empty() )
    throw Exception() << "MatchboxFactory: No process specified."
		      << Exception::runerror;

  vector<PDVector> groups;
  typedef map<string,PDVector>::const_iterator GroupIterator;
  for ( vector<string>::const_iterator gr = proc.begin();
	gr != proc.end(); ++gr ) {
    GroupIterator git = particleGroups().find(*gr);
    if ( git == particleGroups().end() ) {
      throw Exception() << "MatchboxFactory: Particle group '"
			<< *gr << "' not defined." << Exception::runerror;
    }
    groups.push_back(git->second);
  }

  vector<size_t> counts(groups.size(),0);
  PDVector proto(groups.size());

  set<PDVector> allProcs;

  while ( true ) {

    for ( size_t k = 0; k < groups.size(); ++k )
      proto[k] = groups[k][counts[k]];

    int charge = 0;
    int colour = 0;
    int nleptons = 0;
    int nquarks = 0;
    int ncolour = 0;

    int nleptonsGen[4];
    int nquarksGen[4];
    for ( size_t i = 0; i < 4; ++i ) {
      nleptonsGen[i] = 0;
      nquarksGen[i] = 0;
    }

    for ( size_t k = 0; k < proto.size(); ++k ) {
      int sign = k > 1 ? 1 : -1;
      charge += sign * proto[k]->iCharge();
      colour += sign * proto[k]->iColour();
      if ( abs(proto[k]->id()) <= 8 ) {
	int generation = (abs(proto[k]->id()) - 1)/2;
	nquarks += sign * ( proto[k]->id() < 0 ? -1 : 1);
	nquarksGen[generation] += sign * ( proto[k]->id() < 0 ? -1 : 1);
      }
      if ( abs(proto[k]->id()) > 10 &&
	   abs(proto[k]->id()) <= 18 ) {
	int generation = (abs(proto[k]->id()) - 11)/2;
	nleptons += sign * ( proto[k]->id() < 0 ? -1 : 1);
	nleptonsGen[generation] += sign * ( proto[k]->id() < 0 ? -1 : 1);
      }
      if ( proto[k]->coloured() )
	++ncolour;
    }

    bool pass = true;

    if ( theEnforceChargeConservation )
      pass &= (charge == 0);

    if ( theEnforceColourConservation )
      pass &= (colour % 8 == 0);

    if ( theEnforceLeptonNumberConservation ) {
      pass &= (nleptons == 0);
      if ( theLeptonFlavourDiagonal ) {
	for ( size_t i = 0; i < 4; ++i )
	  pass &= (nleptonsGen[i] == 0);
      }
    }

    if ( theEnforceQuarkNumberConservation ) {
      pass &= (nquarks == 0);
      if ( theQuarkFlavourDiagonal ) {
	for ( size_t i = 0; i < 4; ++i )
	  pass &= (nquarksGen[i] == 0);
      }
    }

    if ( pass ) {
      for ( int i = 0; i < 2; ++i ) {
	if ( proto[i]->coloured() &&
	     proto[i]->hardProcessMass() != ZERO )
	  throw Exception()
	    << "Inconsistent flavour scheme detected with massive incoming "
	    << proto[i]->PDGName() << ". Check your setup."
	    << Exception::runerror;
      }
      sort(proto.begin()+2,proto.end(),SortPID());
      allProcs.insert(proto);
    }

    vector<size_t>::reverse_iterator c = counts.rbegin();
    vector<PDVector>::const_reverse_iterator g = groups.rbegin();
    while ( c != counts.rend() ) {
      if ( ++(*c) == g->size() ) {
	*c = 0;
	++c; ++g;
      } else {
	break;
      }
    }
    if ( c == counts.rend() )
      break;

  }

  return allProcs;

}

void MatchboxFactory::Init() {

  static ClassDocumentation<MatchboxFactory> documentation
    ("MatchboxFactory",
     "NLO QCD corrections have been calculated "
     "using Matchbox \\cite{Platzer:2011bc}, \\cite{Matchbox:2015}",
     "%\\cite{Platzer:2011bc}\n"
     "\\bibitem{Platzer:2011bc}\n"
     "S.~Platzer and S.~Gieseke,\n"
     "``Dipole Showers and Automated NLO Matching in Herwig,''\n"
     "arXiv:1109.6256 [hep-ph].\n"
     "%%CITATION = ARXIV:1109.6256;%%\n"
     "%\\cite{Matchbox:2015}\n"
     "\\bibitem{Matchbox:2015}\n"
     "Herwig collaboration,\n"
     "``Precision LHC Event Generation with Herwig,''\n"
     "in preparation.");

  static Reference<MatchboxFactory,Tree2toNGenerator> interfaceDiagramGenerator
    ("DiagramGenerator",
     "Set the diagram generator.",
     &MatchboxFactory::theDiagramGenerator, false, false, true, true, false);
  interfaceDiagramGenerator.rank(-1);

  static Reference<MatchboxFactory,ProcessData> interfaceProcessData
    ("ProcessData",
     "Set the process data object to be used.",
     &MatchboxFactory::theProcessData, false, false, true, true, false);
  interfaceProcessData.rank(-1);

  static Parameter<MatchboxFactory,unsigned int> interfaceOrderInAlphaS
    ("OrderInAlphaS",
     "The order in alpha_s to consider.",
     &MatchboxFactory::theOrderInAlphaS, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxFactory,unsigned int> interfaceOrderInAlphaEW
    ("OrderInAlphaEW",
     "The order in alpha_EW to consider.",
     &MatchboxFactory::theOrderInAlphaEW, 2, 0, 0,
     false, false, Interface::lowerlim);

  static Switch<MatchboxFactory,bool> interfaceBornContributions
    ("BornContributions",
     "Switch on or off the Born contributions.",
     &MatchboxFactory::theBornContributions, true, false, false);
  static SwitchOption interfaceBornContributionsYes
    (interfaceBornContributions,
     "Yes",
     "Switch on Born contributions.",
     true);
  static SwitchOption interfaceBornContributionsNo
    (interfaceBornContributions,
     "No",
     "Switch off Born contributions.",
     false);

  static Switch<MatchboxFactory,bool> interfaceVirtualContributions
    ("VirtualContributions",
     "Switch on or off the virtual contributions.",
     &MatchboxFactory::theVirtualContributions, true, false, false);
  static SwitchOption interfaceVirtualContributionsYes
    (interfaceVirtualContributions,
     "Yes",
     "Switch on virtual contributions.",
     true);
  static SwitchOption interfaceVirtualContributionsNo
    (interfaceVirtualContributions,
     "No",
     "Switch off virtual contributions.",
     false);

  static Switch<MatchboxFactory,bool> interfaceRealContributions
    ("RealContributions",
     "Switch on or off the real contributions.",
     &MatchboxFactory::theRealContributions, true, false, false);
  static SwitchOption interfaceRealContributionsYes
    (interfaceRealContributions,
     "Yes",
     "Switch on real contributions.",
     true);
  static SwitchOption interfaceRealContributionsNo
    (interfaceRealContributions,
     "No",
     "Switch off real contributions.",
     false);

  static Switch<MatchboxFactory,bool> interfaceIndependentVirtuals
    ("IndependentVirtuals",
     "Switch on or off virtual contributions as separate subprocesses.",
     &MatchboxFactory::theIndependentVirtuals, true, false, false);
  static SwitchOption interfaceIndependentVirtualsYes
    (interfaceIndependentVirtuals,
     "Yes",
     "Switch on virtual contributions as separate subprocesses.",
     true);
  static SwitchOption interfaceIndependentVirtualsNo
    (interfaceIndependentVirtuals,
     "No",
     "Switch off virtual contributions as separate subprocesses.",
     false);

  static Switch<MatchboxFactory,bool> interfaceIndependentPKs
    ("IndependentPKOperators",
     "Switch on or off PK oeprators as separate subprocesses.",
     &MatchboxFactory::theIndependentPKs, true, false, false);
  static SwitchOption interfaceIndependentPKsYes
    (interfaceIndependentPKs,
     "Yes",
     "Switch on PK operators as separate subprocesses.",
     true);
  static SwitchOption interfaceIndependentPKsNo
    (interfaceIndependentPKs,
     "No",
     "Switch off PK operators as separate subprocesses.",
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
  static SwitchOption interfaceFixedCouplingsYes
    (interfaceFixedCouplings,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceFixedCouplingsNo
    (interfaceFixedCouplings,
     "No",
     "No",
     false);
  interfaceFixedCouplings.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceFixedQEDCouplings
    ("FixedQEDCouplings",
     "Switch on or off fixed QED couplings.",
     &MatchboxFactory::theFixedQEDCouplings, true, false, false);
  static SwitchOption interfaceFixedQEDCouplingsYes
    (interfaceFixedQEDCouplings,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceFixedQEDCouplingsNo
    (interfaceFixedQEDCouplings,
     "No",
     "No",
     false);
  interfaceFixedQEDCouplings.rank(-1);

  // @TDOO SP to remove this in the code as well
  /*
  static Switch<MatchboxFactory,bool> interfaceVetoScales
    ("VetoScales",
     "Switch on or setting veto scales.",
     &MatchboxFactory::theVetoScales, false, false, false);
  static SwitchOption interfaceVetoScalesYes
    (interfaceVetoScales,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceVetoScalesNo
    (interfaceVetoScales,
     "No",
     "No",
     false);
  */

  static RefVector<MatchboxFactory,MatchboxAmplitude> interfaceAmplitudes
    ("Amplitudes",
     "The amplitude objects.",
     &MatchboxFactory::theAmplitudes, -1, false, false, true, true, false);

  static Switch<MatchboxFactory,bool> interfaceVerbose
    ("Verbose",
     "Print full infomation on each evaluated phase space point.",
     &MatchboxFactory::theVerbose, false, false, false);
  static SwitchOption interfaceVerboseYes
    (interfaceVerbose,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceVerboseNo
    (interfaceVerbose,
     "No",
     "No",
     false);
  interfaceVerbose.rank(-1);
    
  static Switch<MatchboxFactory,bool> interfaceInitVerbose
    ("InitVerbose",
     "Print setup information.",
     &MatchboxFactory::theInitVerbose, false, false, false);
  static SwitchOption interfaceInitVerboseYes
    (interfaceInitVerbose,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceInitVerboseNo
    (interfaceInitVerbose,
     "No",
     "No",
     false);
  interfaceInitVerbose.rank(-1);

  static Parameter<MatchboxFactory,string> interfaceSubtractionData
    ("SubtractionData",
     "Prefix for subtraction check data.",
     &MatchboxFactory::theSubtractionData, "",
     false, false);

  static Switch<MatchboxFactory,int> interfaceSubtractionPlotType
    ("SubtractionPlotType",
     "Switch for controlling what kind of plot is generated for checking the subtraction",
     &MatchboxFactory::theSubtractionPlotType, 1, false, false);
  static SwitchOption interfaceSubtractionPlotTypeLinearRatio
    (interfaceSubtractionPlotType,
     "LinRatio",
     "Switch on the linear plot of the ratio",
     1);
  static SwitchOption interfaceSubtractionPlotTypeLogRelDiff
    (interfaceSubtractionPlotType,
     "LogRelDiff",
     "Switch on the logarithmic plot of the relative difference",
     2);
  
  static Switch<MatchboxFactory,bool> interfaceSubtractionScatterPlot
    ("SubtractionScatterPlot",
     "Switch for controlling whether subtraction data should be plotted for each phase space point individually",
     &MatchboxFactory::theSubtractionScatterPlot, false, false, false);
  static SwitchOption interfaceSubtractionScatterPlotNo
    (interfaceSubtractionScatterPlot,
     "No", "Switch off the scatter plot", false);
  static SwitchOption interfaceSubtractionScatterPlotYes
    (interfaceSubtractionScatterPlot,
     "Yes", "Switch on the scatter plot", true);

  static Parameter<MatchboxFactory,string> interfacePoleData
    ("PoleData",
     "Prefix for subtraction check data.",
     &MatchboxFactory::thePoleData, "",
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
     "Set the process(es) to consider.",
     &MatchboxFactory::doProcess, false);

  static Command<MatchboxFactory> interfaceLoopInducedProcess
    ("LoopInducedProcess",
     "Set the loop induced process(es) to consider.",
     &MatchboxFactory::doLoopInducedProcess, false);

  static Command<MatchboxFactory> interfaceSingleRealProcess
    ("SingleRealProcess",
     "Set the real emission process(es) to consider.",
     &MatchboxFactory::doSingleRealProcess, false);

  static Reference<MatchboxFactory,ShowerApproximation> interfaceShowerApproximation
    ("ShowerApproximation",
     "Set the shower approximation to be considered.",
     &MatchboxFactory::theShowerApproximation, false, false, true, true, false);

  static Switch<MatchboxFactory,bool> interfaceRealEmissionScales
    ("RealEmissionScales",
     "Switch on or off calculation of subtraction scales from real emission kinematics.",
     &MatchboxFactory::theRealEmissionScales, false, false, false);
  static SwitchOption interfaceRealEmissionScalesYes
    (interfaceRealEmissionScales,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceRealEmissionScalesNo
    (interfaceRealEmissionScales,
     "No",
     "No",
     false);
  interfaceRealEmissionScales.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceAllProcesses
    ("AllProcesses",
     "Consider all processes up to a maximum coupling order specified by the coupling order interfaces.",
     &MatchboxFactory::theAllProcesses, false, false, false);
  static SwitchOption interfaceAllProcessesYes
    (interfaceAllProcesses,
     "Yes",
     "Include all processes.",
     true);
  static SwitchOption interfaceAllProcessesNo
    (interfaceAllProcesses,
     "No",
     "Only consider processes matching the exact order in the couplings.",
     false);
  interfaceAllProcesses.rank(-1);

  static RefVector<MatchboxFactory,MatchboxAmplitude> interfaceSelectAmplitudes
    ("SelectAmplitudes",
     "The amplitude objects to be favoured in clashing responsibilities.",
     &MatchboxFactory::theSelectedAmplitudes, -1, false, false, true, true, false);

  static RefVector<MatchboxFactory,MatchboxAmplitude> interfaceDeselectAmplitudes
    ("DeselectAmplitudes",
     "The amplitude objects to be disfavoured in clashing responsibilities.",
     &MatchboxFactory::theDeselectedAmplitudes, -1, false, false, true, true, false);

  static Switch<MatchboxFactory,int> interfaceDipoleSet
    ("DipoleSet",
     "The set of subtraction terms to be considered.",
     &MatchboxFactory::theDipoleSet, 0, false, false);
  static SwitchOption interfaceDipoleSetCataniSeymour
    (interfaceDipoleSet,
     "CataniSeymour",
     "Use default Catani-Seymour dipoles.",
     0);
  interfaceDipoleSet.rank(-1);

  static RefVector<MatchboxFactory,ReweightBase> interfaceReweighters
    ("Reweighters",
     "Reweight objects for matrix elements.",
     &MatchboxFactory::theReweighters, -1, false, false, true, false, false);

  static RefVector<MatchboxFactory,ReweightBase> interfacePreweighters
    ("Preweighters",
     "Preweight objects for matrix elements.",
     &MatchboxFactory::thePreweighters, -1, false, false, true, false, false);

  static Switch<MatchboxFactory,bool> interfaceMECorrectionsOnly
    ("MECorrectionsOnly",
     "Prepare only ME corrections, but no NLO calculation.",
     &MatchboxFactory::theMECorrectionsOnly, false, false, false);
  static SwitchOption interfaceMECorrectionsOnlyYes
    (interfaceMECorrectionsOnly,
     "Yes",
     "Produce only ME corrections.",
     true);
  static SwitchOption interfaceMECorrectionsOnlyNo
    (interfaceMECorrectionsOnly,
     "No",
     "Produce full NLO.",
     false);

  static Switch<MatchboxFactory,bool> interfaceLoopSimCorrections
    ("LoopSimCorrections",
     "Prepare LoopSim corrections.",
     &MatchboxFactory::theLoopSimCorrections, false, false, false);
  static SwitchOption interfaceLoopSimCorrectionsYes
    (interfaceLoopSimCorrections,
     "Yes",
     "Produce loopsim corrections.",
     true);
  static SwitchOption interfaceLoopSimCorrectionsNo
    (interfaceLoopSimCorrections,
     "No",
     "Produce full NLO.",
     false);

  static Switch<MatchboxFactory,bool> interfaceFirstPerturbativePDF
    ("FirstPerturbativePDF",
     "",
     &MatchboxFactory::theFirstPerturbativePDF, true, false, false);
  static SwitchOption interfaceFirstPerturbativePDFYes
    (interfaceFirstPerturbativePDF,
     "Yes",
     "",
     true);
  static SwitchOption interfaceFirstPerturbativePDFNo
    (interfaceFirstPerturbativePDF,
     "No",
     "",
     false);
  interfaceFirstPerturbativePDF.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceSecondPerturbativePDF
    ("SecondPerturbativePDF",
     "",
     &MatchboxFactory::theSecondPerturbativePDF, true, false, false);
  static SwitchOption interfaceSecondPerturbativePDFYes
    (interfaceSecondPerturbativePDF,
     "Yes",
     "",
     true);
  static SwitchOption interfaceSecondPerturbativePDFNo
    (interfaceSecondPerturbativePDF,
     "No",
     "",
     false);
  interfaceSecondPerturbativePDF.rank(-1);

  static Command<MatchboxFactory> interfaceProductionMode
    ("ProductionMode",
     "Switch this factory to production mode.",
     &MatchboxFactory::doProductionMode, false);

  static Switch<MatchboxFactory,bool> interfaceSpinCorrelations
    ("SpinCorrelations",
     "Fill information for the spin correlations, if possible.",
     &MatchboxFactory::theSpinCorrelations, false, false, false);
  static SwitchOption interfaceSpinCorrelationsYes
    (interfaceSpinCorrelations,
     "Yes",
     "",
     true);
  static SwitchOption interfaceSpinCorrelationsNo
    (interfaceSpinCorrelations,
     "No",
     "",
     false);
  
  static Parameter<MatchboxFactory,double> interfaceAlphaParameter
    ("AlphaParameter",
     "Nagy-AlphaParameter.",
     &MatchboxFactory::theAlphaParameter, 1.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Switch<MatchboxFactory,bool> interfaceEnforceChargeConservation
    ("EnforceChargeConservation",
     "Enforce charge conservation while generating the hard process.",
     &MatchboxFactory::theEnforceChargeConservation, true, false, false);
  static SwitchOption interfaceEnforceChargeConservationYes
    (interfaceEnforceChargeConservation,
     "Yes",
     "Enforce charge conservation.",
     true);
  static SwitchOption interfaceEnforceChargeConservationNo
    (interfaceEnforceChargeConservation,
     "No",
     "Do not enforce charge conservation.",
     false);
  interfaceEnforceChargeConservation.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceEnforceColourConservation
    ("EnforceColourConservation",
     "Enforce colour conservation while generating the hard process.",
     &MatchboxFactory::theEnforceColourConservation, false, false, false);
  static SwitchOption interfaceEnforceColourConservationYes
    (interfaceEnforceColourConservation,
     "Yes",
     "Enforce colour conservation.",
     true);
  static SwitchOption interfaceEnforceColourConservationNo
    (interfaceEnforceColourConservation,
     "No",
     "Do not enforce colour conservation.",
     false);
  interfaceEnforceColourConservation.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceEnforceLeptonNumberConservation
    ("EnforceLeptonNumberConservation",
     "Enforce lepton number conservation while generating the hard process.",
     &MatchboxFactory::theEnforceLeptonNumberConservation, false, false, false);
  static SwitchOption interfaceEnforceLeptonNumberConservationYes
    (interfaceEnforceLeptonNumberConservation,
     "Yes",
     "Enforce lepton number conservation.",
     true);
  static SwitchOption interfaceEnforceLeptonNumberConservationNo
    (interfaceEnforceLeptonNumberConservation,
     "No",
     "Do not enforce lepton number conservation.",
     false);
  interfaceEnforceLeptonNumberConservation.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceEnforceQuarkNumberConservation
    ("EnforceQuarkNumberConservation",
     "Enforce quark number conservation while generating the hard process.",
     &MatchboxFactory::theEnforceQuarkNumberConservation, false, false, false);
  static SwitchOption interfaceEnforceQuarkNumberConservationYes
    (interfaceEnforceQuarkNumberConservation,
     "Yes",
     "Enforce quark number conservation.",
     true);
  static SwitchOption interfaceEnforceQuarkNumberConservationNo
    (interfaceEnforceQuarkNumberConservation,
     "No",
     "Do not enforce quark number conservation.",
     false);
  interfaceEnforceQuarkNumberConservation.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceLeptonFlavourDiagonal
    ("LeptonFlavourDiagonal",
     "Assume that lepton interactions are flavour diagonal while generating the hard process.",
     &MatchboxFactory::theLeptonFlavourDiagonal, false, false, false);
  static SwitchOption interfaceLeptonFlavourDiagonalYes
    (interfaceLeptonFlavourDiagonal,
     "Yes",
     "Assume that lepton interactions are flavour diagonal.",
     true);
  static SwitchOption interfaceLeptonFlavourDiagonalNo
    (interfaceLeptonFlavourDiagonal,
     "No",
     "Do not assume that lepton interactions are flavour diagonal.",
     false);
  interfaceLeptonFlavourDiagonal.rank(-1);

  static Switch<MatchboxFactory,bool> interfaceQuarkFlavourDiagonal
    ("QuarkFlavourDiagonal",
     "Assume that quark interactions are flavour diagonal while generating the hard process.",
     &MatchboxFactory::theQuarkFlavourDiagonal, false, false, false);
  static SwitchOption interfaceQuarkFlavourDiagonalYes
    (interfaceQuarkFlavourDiagonal,
     "Yes",
     "Assume that quark interactions are flavour diagonal.",
     true);
  static SwitchOption interfaceQuarkFlavourDiagonalNo
    (interfaceQuarkFlavourDiagonal,
     "No",
     "Do not assume that quark interactions are flavour diagonal.",
     false);
  interfaceQuarkFlavourDiagonal.rank(-1);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxFactory,SubProcessHandler>
describeHerwigMatchboxFactory("Herwig::MatchboxFactory", "Herwig.so");
