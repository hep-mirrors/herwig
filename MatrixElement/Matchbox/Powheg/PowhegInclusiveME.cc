// -*- C++ -*-
//
// PowhegInclusiveME.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegInclusiveME class.
//

#include "PowhegInclusiveME.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/Rebinder.h"

#include "ThePEG/Handlers/StdXCombGroup.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Powheg/PowhegInclusiveReweight.h"

using namespace Herwig;

PowhegInclusiveME::PowhegInclusiveME() 
  : MEGroup(), theVerbose(false), theMCSum(false) {}

PowhegInclusiveME::~PowhegInclusiveME() {
}

IBPtr PowhegInclusiveME::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegInclusiveME::fullclone() const {
  return new_ptr(*this);
}

MEBase::DiagramVector PowhegInclusiveME::dependentDiagrams(const cPDVector& proc,
							   tMEPtr depME) const {

  Ptr<SubtractionDipole>::tptr dipole = 
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(depME);

  return dipole->realEmissionDiagrams(proc);

}

void PowhegInclusiveME::setXComb(tStdXCombPtr xc) {

  MEGroup::setXComb(xc);
  MEVector::const_iterator me = dependent().begin();
  tStdXCombGroupPtr group = dynamic_ptr_cast<tStdXCombGroupPtr>(xc);
  assert(group);
  StdDepXCVector::const_iterator depxc = group->dependent().begin();
  for ( ; me != dependent().end(); ++me, ++depxc ) {
    theKernelMap[*me]->setXComb(*depxc);
  }

}

vector<Ptr<SubtractionDipole>::ptr> PowhegInclusiveME::dipoles() const {
  vector<Ptr<SubtractionDipole>::ptr> res;
  for ( MEVector::const_iterator k = dependent().begin();
	k != dependent().end(); ++k )
    res.push_back(dynamic_ptr_cast<Ptr<SubtractionDipole>::ptr>(*k));
  return res;
}


void PowhegInclusiveME::setup(Ptr<MatchboxNLOME>::ptr newBornVirtual,
			      const vector<Ptr<SubtractedME>::ptr>& newRealEmissions,
			      bool bornScreening) {

  head(newBornVirtual);
  theSplittingKernels.clear();
  theKernelMap.clear();

  vector<Ptr<SubtractionDipole>::ptr> dipoles;

  for ( vector<Ptr<SubtractedME>::ptr>::const_iterator realit =
	  newRealEmissions.begin(); realit != newRealEmissions.end(); ++realit ) {

    cPDVector rep = newBornVirtual->matrixElement()->diagrams().front()->partons();

    vector<Ptr<SubtractionDipole>::ptr> dips =
      (**realit).splitDipoles(rep);

    for ( vector<Ptr<SubtractionDipole>::ptr>::iterator dip =
	    dips.begin(); dip != dips.end(); ++dip ) {

      Ptr<SubtractionDipole>::ptr pdip = (**dip).cloneMe();
      string pdipname = pdip->fullName();
      ostringstream pdname;
      pdname << pdipname << ".Projection";
      if ( !(generator()->preinitRegister(pdip,pdname.str())) )
	throw InitException() << "Dipole '" << pdname.str() << "' already existing.";
      pdip->cloneDependencies();

      Ptr<PowhegInclusiveReweight>::ptr reweight =
	new_ptr(PowhegInclusiveReweight());

      ostringstream rwname;
      rwname << pdipname << ".InclusiveReweight";
      if ( !(generator()->preinitRegister(reweight,rwname.str())) )
	throw InitException() << "Reweight '" << rwname.str() << "' already existing.";

      reweight->projectionDipole(pdip);
      reweight->setup(*realit);

      pdip->doSplitting();
      pdip->addReweight(reweight);

      Ptr<SubtractionDipole>::ptr sdip = (**dip).cloneMe();
      string sdipname = sdip->fullName();

      ostringstream sdname;
      sdname << sdipname << ".SplittingDipole";
      if ( !(generator()->preinitRegister(sdip,sdname.str())) )
	throw InitException() << "Dipole '" << sdname.str() << "' already existing.";
      sdip->cloneDependencies();

      Ptr<PowhegSplittingKernel>::ptr split
	= new_ptr(PowhegSplittingKernel());

      ostringstream skname;
      skname << sdipname << ".SplittingKernel";
      if ( !(generator()->preinitRegister(split,skname.str())) )
	throw InitException() << "Splitting kernel '" << skname.str() << "' already existing.";

      split->projectionDipole(sdip);
      split->setup(*realit);
      sdip->doSplitting();

      if ( bornScreening ) {
	reweight->doBornScreening();
	split->doBornScreening();
      } else {
	reweight->noBornScreening();
	split->noBornScreening();
      }

      if ( theVerbose ) {
	sdip->print(Repository::clog());
	split->print(Repository::clog());
      }

      dipoles.push_back(pdip);
      theSplittingKernels.push_back(split);
      theKernelMap[pdip] = split;

    }

  }

  MEVector dipMEs;
  dipMEs.resize(dipoles.size());
  copy(dipoles.begin(),dipoles.end(),dipMEs.begin());

  dependent() = dipMEs;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void PowhegInclusiveME::doinit() {

  if ( theVerbose )
    print(Repository::clog());

  MEGroup::doinit();
}

void PowhegInclusiveME::print(ostream& os) const {

  os << "--- PowhegInclusiveME setup ----------------------------------------------------\n";

  os << " '" << name() << "' for Born/virtual\n '"
     << head()->name() << "':\n";

  dynamic_ptr_cast<Ptr<MatchboxNLOME>::tptr>(head())->print(os);

  os << " using the dipoles:\n";

  for ( MEVector::const_iterator d = dependent().begin();
	d != dependent().end(); ++d ) {
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*d)->name();
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*d)->reweights().front()->print(os);
  }

  os << " generated splitting kernels:\n";

  for ( vector<Ptr<PowhegSplittingKernel>::ptr>::const_iterator sp
	  = splittingKernels().begin(); sp != splittingKernels().end(); ++sp )
    (**sp).print(os);

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void PowhegInclusiveME::printLastEvent(ostream& os) const {

  os << "--- PowhegInclusiveME last event information -----------------------------------\n";

  os << " '" << name() << "' for Born/virtual\n '"
     << head()->name() << "'\n";

  os << " Born/virtual event information:\n";
  dynamic_ptr_cast<Ptr<MatchboxNLOME>::tptr>(head())->printLastEvent(os);

  os << " dipoles event information:\n";
  for ( MEVector::const_iterator d = dependent().begin();
	d != dependent().end(); ++d ) {
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*d)->printLastEvent(os);
    dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(*d)->reweights().front()->printLastEvent(os);
  }


  os << "--- end PowhegInclusiveME last event information -------------------------------\n\n\n";

  os << flush;

}

void PowhegInclusiveME::lastEventStatistics() {
  MEGroup::lastEventStatistics();
  if ( !generator() )
    return;
  /*
  if ( theVerbose )
    printLastEvent(generator()->log());
  */
}


void PowhegInclusiveME::persistentOutput(PersistentOStream & os) const {
  os << theSplittingKernels << theKernelMap << theVerbose << theMCSum;
}

void PowhegInclusiveME::persistentInput(PersistentIStream & is, int) {
  is >> theSplittingKernels >> theKernelMap >> theVerbose >> theMCSum;
}

void PowhegInclusiveME::rebind(const TranslationMap & trans) {
  map<Ptr<MEBase>::ptr,Ptr<PowhegSplittingKernel>::ptr> newKernelMap;
  for ( map<Ptr<MEBase>::ptr,Ptr<PowhegSplittingKernel>::ptr>::const_iterator mit =
	  theKernelMap.begin(); mit != theKernelMap.end(); ++mit ) {
    newKernelMap[trans.translate(mit->first)] =
      trans.translate(mit->second);
  }
  theKernelMap = newKernelMap;
  MEGroup::rebind(trans);
}

IVector PowhegInclusiveME::getReferences() {
  IVector ret = MEGroup::getReferences();
  for ( map<Ptr<MEBase>::ptr,Ptr<PowhegSplittingKernel>::ptr>::const_iterator mit =
	  theKernelMap.begin(); mit != theKernelMap.end(); ++mit ) {
    ret.push_back(mit->first);
    ret.push_back(mit->second);
  }
  return ret;
}


void PowhegInclusiveME::Init() {

  static ClassDocumentation<PowhegInclusiveME> documentation
    ("PowhegInclusiveME represents a BBar function.");

  static RefVector<PowhegInclusiveME,PowhegSplittingKernel> interfaceSplittingKernels
    ("SplittingKernels",
     "The splitting kernels to be used.",
     &PowhegInclusiveME::theSplittingKernels, -1, false, false, true, true, false);

  static Switch<PowhegInclusiveME,bool> interfaceVerbose
    ("Verbose",
     "Print full infomation on each evaluated phase space point.",
     &PowhegInclusiveME::theVerbose, false, false, false);
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

  static Switch<PowhegInclusiveME,bool> interfaceMCSum
    ("MCSum",
     "MC sum over eal emission contributions.",
     &PowhegInclusiveME::theMCSum, false, false, false);
  static SwitchOption interfaceMCSumOn
    (interfaceMCSum,
     "On",
     "On",
     true);
  static SwitchOption interfaceMCSumOff
    (interfaceMCSum,
     "Off",
     "Off",
     false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PowhegInclusiveME,MEGroup>
describeHerwigPowhegInclusiveME("Herwig::PowhegInclusiveME", "HwMatchbox.so");
