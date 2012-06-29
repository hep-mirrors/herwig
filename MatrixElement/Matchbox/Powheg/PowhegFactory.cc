// -*- C++ -*-
//
// PowhegFactory.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegFactory class.
//

#include "PowhegFactory.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/Repository.h"
#include "Herwig++/MatrixElement/Matchbox/Powheg/PowhegRealReweight.h"

using namespace Herwig;

PowhegFactory::PowhegFactory() 
  : SubProcessHandler(), theBornScreening(true), theVerbose(false) {}

PowhegFactory::~PowhegFactory() {}

IBPtr PowhegFactory::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegFactory::fullclone() const {
  return new_ptr(*this);
}

void PowhegFactory::setup() {

  theInclusiveMEs.clear();
  theRealMEs.clear();
  MEs().clear();

  if ( theMatchboxFactory ) {
    theMatchboxFactory->init();
    theBornVirtuals = theMatchboxFactory->bornVirtualMEs();
    theSubtractedMEs = theMatchboxFactory->subtractedMEs();
  }

  for ( vector<Ptr<MatchboxNLOME>::ptr>::const_iterator bv =
	  bornVirtuals().begin(); bv != bornVirtuals().end(); ++bv ) {

    Ptr<PowhegInclusiveME>::ptr bbar = new_ptr(PowhegInclusiveME());
    string pname = fullName() + "/" + (**bv).matrixElement()->name();
    if ( ! (generator()->preinitRegister(bbar,pname) ) )
      throw InitException() << "Powheg Inclusive ME " << pname << " already existing.";

    if ( theVerbose )
      bbar->beVerbose();
    else
      bbar->beQuiet();
    bbar->setup(*bv,theSubtractedMEs,theBornScreening);

    theInclusiveMEs.push_back(bbar);
    MEs().push_back(bbar);
  }

  if ( theBornScreening ) {

    for ( vector<Ptr<SubtractedME>::ptr>::const_iterator real =
	    subtractedMEs().begin(); real != subtractedMEs().end(); ++real ) {

      Ptr<MatchboxMEBase>::ptr realME = 
	dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>((**real).head());
      assert(realME);

      Ptr<MatchboxMEBase>::ptr finitereal = realME->cloneMe();

      string pname = fullName() + "/" + (**real).name();
      if ( ! (generator()->preinitRegister(finitereal,pname) ) )
	throw InitException() << "Powheg Finite Real ME " << pname << " already existing.";

      finitereal->reweights().clear();

      vector<Ptr<SubtractionDipole>::ptr> dips = (**real).dipoles();

      for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator dip
	      = dips.begin(); dip != dips.end(); ++dip ) {
	Ptr<PowhegRealReweight>::ptr reweight = new_ptr(PowhegRealReweight());

	string rwname = pname + "/" + (**dip).name() + ".FiniteRealReweight";
	if ( !(generator()->preinitRegister(reweight,rwname)) )
	  throw InitException() << "Reweight '" << rwname << "' already existing.";

	reweight->setup(*dip,*real);
	finitereal->addReweight(reweight);
      }

      theRealMEs.push_back(finitereal);
      MEs().push_back(finitereal);

    }

  }

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PowhegFactory::doinit() {
  setup();
  if ( theVerbose )
    print(Repository::clog());
  SubProcessHandler::doinit();
}

void PowhegFactory::print(ostream& os) const {

  os << "--- PowhegFactory setup --------------------------------------------------------\n";

  os << " inclusive matrix elements generated:\n";

  for ( vector<Ptr<PowhegInclusiveME>::ptr>::const_iterator bbar
	  = theInclusiveMEs.begin(); bbar != theInclusiveMEs.end(); ++bbar ) {
    os << " '" << (**bbar).name() << "'\n";
  }

  if ( theBornScreening ) {

    os << " finite real emission matrix elements generated:\n";

    for ( vector<Ptr<MatchboxMEBase>::ptr>::const_iterator real =
	    theRealMEs.begin(); real != theRealMEs.end(); ++real ) {
      os << "'" << (**real).name() << "'\n";
    }

  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void PowhegFactory::persistentOutput(PersistentOStream & os) const {
  os << theBornVirtuals << theSubtractedMEs 
     << theMatchboxFactory << theBornScreening 
     << theInclusiveMEs << theRealMEs << theVerbose;
}

void PowhegFactory::persistentInput(PersistentIStream & is, int) {
  is >> theBornVirtuals >> theSubtractedMEs 
     >> theMatchboxFactory >> theBornScreening 
     >> theInclusiveMEs >> theRealMEs >> theVerbose;
}

void PowhegFactory::Init() {

  static ClassDocumentation<PowhegFactory> documentation
    ("PowhegFactory",
     "NLO QCD corrections and POWHEG matching have been calculated "
     "using Matchbox \\cite{Platzer:2011bc}",
     "%\\cite{Platzer:2011bc}\n"
     "\\bibitem{Platzer:2011bc}\n"
     "S.~Platzer and S.~Gieseke,\n"
     "``Dipole Showers and Automated NLO Matching in Herwig++,''\n"
     "arXiv:1109.6256 [hep-ph].\n"
     "%%CITATION = ARXIV:1109.6256;%%");


  static RefVector<PowhegFactory,MatchboxNLOME> interfaceBornVirtuals
    ("BornVirtuals",
     "Born processes along with virtual corrections to consider",
     &PowhegFactory::theBornVirtuals, -1, false, false, true, true, false);


  static RefVector<PowhegFactory,SubtractedME> interfaceSubtractedMEs
    ("SubtractedMEs",
     "The subtracted real emission matrix elements to consider",
     &PowhegFactory::theSubtractedMEs, -1, false, false, true, true, false);


  static Switch<PowhegFactory,bool> interfaceBornScreening
    ("BornScreening",
     "Switch on or off Born screening",
     &PowhegFactory::theBornScreening, true, false, false);
  static SwitchOption interfaceBornScreeningOn
    (interfaceBornScreening,
     "On",
     "Perform Born screening",
     true);
  static SwitchOption interfaceBornScreeningOff
    (interfaceBornScreening,
     "Off",
     "Do not perform Born screening",
     false);

  static Reference<PowhegFactory,MatchboxFactory> interfaceMatchboxFactory
    ("MatchboxFactory",
     "An optional MatchboxFactory object to pick matrix elements from.",
     &PowhegFactory::theMatchboxFactory, false, false, true, true, false);

  static RefVector<PowhegFactory,PowhegInclusiveME> interfaceInclusiveMEs
    ("InclusiveMEs",
     "The inclusive matrix elements generated",
     &PowhegFactory::theInclusiveMEs, -1, false, true, true, true, false);

  static RefVector<PowhegFactory,MatchboxMEBase> interfaceRealMEs
    ("RealMEs",
     "The finite real matrix elements generated",
     &PowhegFactory::theRealMEs, -1, false, true, true, true, false);

  static Switch<PowhegFactory,bool> interfaceVerbose
    ("Verbose",
     "Print full infomation on each evaluated phase space point.",
     &PowhegFactory::theVerbose, false, false, false);
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

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PowhegFactory,SubProcessHandler>
describeHerwigPowhegFactory("Herwig::PowhegFactory", "HwMatchbox.so");
