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
    theDipoleSet(0), theVerbose(false), theInitVerbose(false), 
    theSubtractionData(""), theSubtractionPlotType(1), theSubtractionScatterPlot(false),
    thePoleData(""), theRealEmissionScales(false), theAllProcesses(false),
    theMECorrectionsOnly(false) {}

MatchboxFactory::~MatchboxFactory() {}

MatchboxFactory*& MatchboxFactory::theCurrentFactory() {
  static MatchboxFactory* sCurrentFactory = 0;
  return sCurrentFactory;
}

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

  me->factory(this);

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
makeMEs(const vector<string>& proc, unsigned int orderas) {

  generator()->log() << "determining subprocesses for ";
  copy(proc.begin(),proc.end(),ostream_iterator<string>(generator()->log()," "));
  generator()->log() << flush;

  map<Ptr<MatchboxAmplitude>::ptr,set<Process> > ampProcs;
  map<Process,set<Ptr<MatchboxAmplitude>::ptr> > procAmps;
  set<PDVector> processes = makeSubProcesses(proc);

  bool needUnsorted = false;

  for ( vector<Ptr<MatchboxAmplitude>::ptr>::const_iterator amp
	  = amplitudes().begin(); amp != amplitudes().end(); ++amp ) {
    if ( !(**amp).sortOutgoing() ) {
      needUnsorted = true;
      break;
    }
  }

  set<PDVector> unsortedProcesses;
  if ( needUnsorted )
    unsortedProcesses = makeUnsortedSubProcesses(proc);

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

  size_t combinations =  processes.size()*matchAmplitudes.size()
    + unsortedProcesses.size()*matchAmplitudes.size();
  size_t procCount = 0;

  boost::progress_display * progressBar = 
    new boost::progress_display(combinations,generator()->log());

  for ( unsigned int oas = lowestAsOrder; oas <= highestAsOrder; ++oas ) {
    for ( unsigned int oae = lowestAeOrder; oae <= highestAeOrder; ++oae ) {
      for ( vector<Ptr<MatchboxAmplitude>::ptr>::const_iterator amp
	      = matchAmplitudes.begin(); amp != matchAmplitudes.end(); ++amp ) {
	(**amp).orderInGs(oas);
	(**amp).orderInGem(oae);
	for ( set<PDVector>::const_iterator p = processes.begin();
	      p != processes.end(); ++p ) {
	  ++(*progressBar);
	  if ( !(**amp).canHandle(*p,this) || !(**amp).sortOutgoing() )
	    continue;
	  ++procCount;
	  Process proc(*p,oas,oae);
	  ampProcs[*amp].insert(proc);
	  procAmps[proc].insert(*amp);
	}
	for ( set<PDVector>::const_iterator p = unsortedProcesses.begin();
	      p != unsortedProcesses.end(); ++p ) {
	  ++(*progressBar);
	  if ( !(**amp).canHandle(*p,this) || (**amp).sortOutgoing() )
	    continue;
	  ++procCount;
	  Process proc(*p,oas,oae);
	  ampProcs[*amp].insert(proc);
	  procAmps[proc].insert(*amp);
	}
      }
    }
  }

  delete progressBar;
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
    throw InitException()
      << "Ambiguous amplitude setup - please check your input files.\n"
      << "To avoid this problem use the SelectAmplitudes or DeselectAmplitudes interfaces.\n";
  }

  vector<Ptr<MatchboxMEBase>::ptr> res;
  for ( map<Ptr<MatchboxAmplitude>::ptr,set<Process> >::const_iterator
	  ap = ampProcs.begin(); ap != ampProcs.end(); ++ap ) {
    for ( set<Process>::const_iterator m = ap->second.begin();
	  m != ap->second.end(); ++m ) {
      Ptr<MatchboxMEBase>::ptr me = ap->first->makeME(m->legs);
      me->subProcess() = *m;
      me->amplitude(ap->first);
      me->matchboxAmplitude(ap->first);
      prepareME(me);
      string pname = "ME" + ap->first->name() + pid(m->legs);
      if ( ! (generator()->preinitRegister(me,pname) ) )
	throw InitException() << "Matrix element " << pname << " already existing.";
      if ( me->diagrams().empty() )continue;
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

void MatchboxFactory::setup() {

  olpProcesses().clear();

  if ( bornMEs().empty() ) {

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

    vector<Ptr<MatchboxMEBase>::ptr> mes;
    for ( vector<vector<string> >::const_iterator p = processes.begin();
	  p != processes.end(); ++p ) {
      mes = makeMEs(*p,orderInAlphaS());
      copy(mes.begin(),mes.end(),back_inserter(bornMEs()));
      if ( realContributions() && realEmissionMEs().empty() ) {
	if ( realEmissionProcesses.empty() ) {
	  vector<string> rproc = *p;
	  rproc.push_back("j");
	  mes = makeMEs(rproc,orderInAlphaS()+1);
	  copy(mes.begin(),mes.end(),back_inserter(realEmissionMEs()));
	}
      }
    }
    if ( realContributions() && realEmissionMEs().empty() ) {
      if ( !realEmissionProcesses.empty() ) {
	for ( vector<vector<string> >::const_iterator q =
		realEmissionProcesses.begin(); q != realEmissionProcesses.end(); ++q ) {
	  mes = makeMEs(*q,orderInAlphaS()+1);
	  copy(mes.begin(),mes.end(),back_inserter(realEmissionMEs()));
	}
      }
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

  // renormalization scheme
  bool virtualsAreDRbar = false;

  // check and prepare the Born and virtual matrix elements
  for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	  = bornMEs().begin(); born != bornMEs().end(); ++born ) {
    prepareME(*born);
    haveVirtuals &= (**born).haveOneLoop();
    if ( (**born).haveOneLoop() ) {
      virtualsAreDRbar |= (**born).isDRbar();
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
    virtualsAreDRbar |= (**virt).isDRbar();
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
	    = DipoleRepository::insertionOperators(dipoleSet()).begin(); 
	  virt != DipoleRepository::insertionOperators(dipoleSet()).end(); ++virt ) {
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
    generator()->log() << "preparing Born"
		       << (virtualContributions() ? " and virtual" : "")
		       << " matrix elements.\n" << flush;
  }

  if ( (bornContributions() && !virtualContributions()) || 
       (bornContributions() && virtualContributions() && independentVirtuals()) ) {
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

      if ( bornme->isOLPTree() ) {
	int id = orderOLPProcess(bornme->subProcess(),
				 (**born).matchboxAmplitude(),
				 ProcessType::treeME2);
	bornme->olpProcess(ProcessType::treeME2,id);
      }

      bornme->cloneDependencies();
      MEs().push_back(bornme);

    }
  }

  if ( virtualContributions() && !meCorrectionsOnly() ) {

    bornVirtualMEs().clear();

    boost::progress_display * progressBar = 
      new boost::progress_display(bornMEs().size(),generator()->log());

    if ( thePoleData != "" )
      if ( thePoleData[thePoleData.size()-1] != '/' )
	thePoleData += "/";

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
		= DipoleRepository::insertionOperators(dipoleSet()).begin(); 
	      virt != DipoleRepository::insertionOperators(dipoleSet()).end(); ++virt ) {
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
  if ( showerApproximation() && !virtualContributions() && !realContributions() ) {
    generator()->log() << "Warning: Matching requested for LO run.\n" << flush;
    showerApproximation(Ptr<ShowerApproximation>::tptr());
  }
  if ( showerApproximation() ) {
    if ( showerApproximation()->needsSplittingGenerator() ) {
      for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	      = bornMEs().begin(); born != bornMEs().end(); ++born )
	for ( MEBase::DiagramVector::const_iterator d = (**born).diagrams().begin();
	      d != (**born).diagrams().end(); ++d )
	  bornProcs.insert((**d).partons());
    }
  }

  if ( realContributions() || meCorrectionsOnly() ) {

    generator()->log() << "preparing real emission matrix elements.\n" << flush;

    if ( theSubtractionData != "" )
      if ( theSubtractionData[theSubtractionData.size()-1] != '/' )
	theSubtractionData += "/";

    subtractedMEs().clear();    

    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator born
	    = bornMEs().begin(); born != bornMEs().end(); ++born ) {

      if ( (**born).onlyOneLoop() )
	continue;

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
      }

    }

    boost::progress_display * progressBar = 
      new boost::progress_display(realEmissionMEs().size(),generator()->log());

    for ( vector<Ptr<MatchboxMEBase>::ptr>::iterator real
	    = realEmissionMEs().begin(); real != realEmissionMEs().end(); ++real ) {

      Ptr<SubtractedME>::ptr sub = new_ptr(SubtractedME());
      string pname = fullName() + "/" + (**real).name() + ".SubtractedReal";
      if ( ! (generator()->preinitRegister(sub,pname) ) )
	throw InitException() << "Subtracted ME " << pname << " already existing.";

      sub->factory(this);

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
	Ptr<MatchboxMEBase>::ptr fme = 
	  dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(sub->head())->cloneMe();
	string qname = fullName() + "/" + (**real).name() + ".FiniteReal";
	if ( ! (generator()->preinitRegister(fme,qname) ) )
	  throw InitException() << "ME " << qname << " already existing.";
	MEs().push_back(fme);	
        finiteRealMEs().push_back(fme);
	sub->head(tMEPtr());
	continue;
      }

      if ( realEmissionScales() )
	sub->doRealEmissionScales();

      subtractedMEs().push_back(sub);
      if ( !meCorrectionsOnly() )
	MEs().push_back(sub);

      if ( showerApproximation() ) {
	if ( virtualContributions() && !meCorrectionsOnly() ) {
	  Ptr<SubtractedME>::ptr subv = new_ptr(*sub);
	  string vname = sub->fullName() + ".SubtractionIntegral";
	  if ( ! (generator()->preinitRegister(subv,vname) ) )
	    throw InitException() << "Subtracted ME " << vname << " already existing.";
	  subv->cloneDependencies(vname);
	  subv->doVirtualShowerSubtraction();
	  subtractedMEs().push_back(subv);
	  MEs().push_back(subv);
	}
	if ( !meCorrectionsOnly() )
	  sub->doRealShowerSubtraction();
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

  if ( !olpProcesses().empty() ) {
    generator()->log() << "Initializing one-loop provider(s).\n" << flush;
    map<Ptr<MatchboxAmplitude>::tptr,map<pair<Process,int>,int> > olps;
    for ( map<Ptr<MatchboxAmplitude>::tptr,map<pair<Process,int>,int> >::const_iterator
	    oit = olpProcesses().begin(); oit != olpProcesses().end(); ++oit ) {
      olps[oit->first] = oit->second;
    }
    boost::progress_display * progressBar = 
      new boost::progress_display(olps.size(),generator()->log());
    for ( map<Ptr<MatchboxAmplitude>::tptr,map<pair<Process,int>,int> >::const_iterator
	    olpit = olps.begin(); olpit != olps.end(); ++olpit ) {
      if ( !olpit->first->startOLP(olpit->second) ) {
	throw InitException() 
	  << "error: failed to start OLP for amplitude '" << olpit->first->name() << "'\n";
      }
      ++(*progressBar);
    }
    delete progressBar;
    generator()->log() << "--------------------------------------------------------------------------------\n"
		       << flush;
  }

  generator()->log() << "Process setup finished.\n" << flush;

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

  os << "--------------------------------------------------------------------------------\n";

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
      << "-------------------------------------------------------------------------------------\n"
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

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void MatchboxFactory::doinit() {
  theCurrentFactory() = this;
  setup();
  if ( initVerbose() )
    print(Repository::clog());
  SubProcessHandler::doinit();
}

void MatchboxFactory::doinitrun() {
  theCurrentFactory() = this;
  SubProcessHandler::doinitrun();
}


void MatchboxFactory::persistentOutput(PersistentOStream & os) const {
  os << theDiagramGenerator << theProcessData
     << theNLight << theOrderInAlphaS << theOrderInAlphaEW 
     << theBornContributions << theVirtualContributions
     << theRealContributions << theIndependentVirtuals << theSubProcessGroups << theInclusive
     << thePhasespace << theScaleChoice
     << theFactorizationScaleFactor << theRenormalizationScaleFactor
     << theFixedCouplings << theFixedQEDCouplings << theVetoScales
     << theAmplitudes
     << theBornMEs << theVirtuals << theRealEmissionMEs
     << theBornVirtualMEs << theSubtractedMEs << theFiniteRealMEs
     << theVerbose << theInitVerbose << theSubtractionData << theSubtractionPlotType
     << theSubtractionScatterPlot << thePoleData
     << theParticleGroups << processes << realEmissionProcesses
     << theShowerApproximation << theSplittingDipoles
     << theRealEmissionScales << theAllProcesses
     << theOLPProcesses 
     << theSelectedAmplitudes << theDeselectedAmplitudes
     << theDipoleSet << theReweighters << thePreweighters
     << theMECorrectionsOnly;
}

void MatchboxFactory::persistentInput(PersistentIStream & is, int) {
  is >> theDiagramGenerator >> theProcessData
     >> theNLight >> theOrderInAlphaS >> theOrderInAlphaEW 
     >> theBornContributions >> theVirtualContributions
     >> theRealContributions >> theIndependentVirtuals >> theSubProcessGroups >> theInclusive
     >> thePhasespace >> theScaleChoice
     >> theFactorizationScaleFactor >> theRenormalizationScaleFactor
     >> theFixedCouplings >> theFixedQEDCouplings >> theVetoScales
     >> theAmplitudes
     >> theBornMEs >> theVirtuals >> theRealEmissionMEs
     >> theBornVirtualMEs >> theSubtractedMEs >> theFiniteRealMEs
     >> theVerbose >> theInitVerbose >> theSubtractionData >> theSubtractionPlotType
     >> theSubtractionScatterPlot >> thePoleData
     >> theParticleGroups >> processes >> realEmissionProcesses
     >> theShowerApproximation >> theSplittingDipoles
     >> theRealEmissionScales >> theAllProcesses
     >> theOLPProcesses
     >> theSelectedAmplitudes >> theDeselectedAmplitudes
     >> theDipoleSet >> theReweighters >> thePreweighters
     >> theMECorrectionsOnly;
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
  vector<string> process = StringUtils::split(in);
  if ( process.size() < 3 )
    throw InitException() << "Invalid process.";
  for ( vector<string>::iterator p = process.begin();
	p != process.end(); ++p ) {
    *p = StringUtils::stripws(*p);
  }
  processes.push_back(process);
  return "";
}

string MatchboxFactory::doSingleRealProcess(string in) {
  vector<string> realEmissionProcess = StringUtils::split(in);
  if ( realEmissionProcess.size() < 3 )
    throw InitException() << "Invalid process.";
  for ( vector<string>::iterator p = realEmissionProcess.begin();
	p != realEmissionProcess.end(); ++p ) {
    *p = StringUtils::stripws(*p);
  }
  realEmissionProcesses.push_back(realEmissionProcess);
  return "";
}

struct SortPID {
  inline bool operator()(PDPtr a, PDPtr b) const {
    return a->id() < b->id();
  }
};

set<PDVector> MatchboxFactory::
makeSubProcesses(const vector<string>& proc, bool sorted) const {

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
    if ( sorted )
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
  static SwitchOption interfaceSubtractionScatterPlotOff
    (interfaceSubtractionScatterPlot,
     "Off", "Switch off the scatter plot", false);
  static SwitchOption interfaceSubtractionScatterPlotOn
    (interfaceSubtractionScatterPlot,
     "On", "Switch on the scatter plot", true);

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

  static Switch<MatchboxFactory,bool> interfaceRealEmissionScales
    ("RealEmissionScales",
     "Switch on or off calculation of subtraction scales from real emission kinematics.",
     &MatchboxFactory::theRealEmissionScales, false, false, false);
  static SwitchOption interfaceRealEmissionScalesOn
    (interfaceRealEmissionScales,
     "On",
     "On",
     true);
  static SwitchOption interfaceRealEmissionScalesOff
    (interfaceRealEmissionScales,
     "Off",
     "Off",
     false);

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

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxFactory,SubProcessHandler>
describeHerwigMatchboxFactory("Herwig::MatchboxFactory", "Herwig.so");
