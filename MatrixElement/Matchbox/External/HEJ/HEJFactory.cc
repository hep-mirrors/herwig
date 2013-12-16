// -*- C++ -*-
//
// HEJFactory.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HEJFactory class.
//

#include "HEJFactory.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "HEJMEBase.h"

using namespace Herwig;

HEJFactory::HEJFactory() 
  : theNGluons(20),
    theJetPTMin(20.*GeV), theExtremalPTMin(15.*GeV), 
    theScaleChoice(1),
    theFixedScale(40.*GeV), theJets(0) {}

HEJFactory::~HEJFactory() {
  if ( theJets ) {
    MEs().clear();
    delete theJets;
    theJets = 0;
  }
}

IBPtr HEJFactory::clone() const {
  return new_ptr(*this);
}

IBPtr HEJFactory::fullclone() const {
  return new_ptr(*this);
}

void HEJFactory::makeME(cPDPtr ini, cPDPtr inj, unsigned int n) {

  theHEJME->firstIncoming(ini);
  theHEJME->secondIncoming(inj);
  theHEJME->nGluons(n);

  ostringstream mename;
  mename << fullName() << "/"
	 << ini->PDGName() << inj->PDGName()
	 << "to"
	 << ini->PDGName() << inj->PDGName()
	 << n << "g";

  string pname = mename.str();
  Ptr<HEJMEBase>::ptr myme = theHEJME->cloneMe();
  if ( ! (generator()->preinitRegister(myme,pname) ) )
    throw InitException() << "Matrix element " << pname << " already existing.";
  myme->cloneDependencies();
  myme->factory(this);
  MEs().push_back(myme);

}

SGeneratorFlags HEJFactory::setup(tStdXCombPtr) const {

  SGeneratorFlags flags;

  if ( abs(generator()->eventHandler()->incoming().first->id()) != ParticleID::pplus ||
       abs(generator()->eventHandler()->incoming().second->id()) != ParticleID::pplus )
    throw Exception() << "can only handle proton/antiproton beams so far\n"
		      << Exception::abortnow;

  flags.beam[0] = generator()->eventHandler()->incoming().first->id() > 0 ? 1 : -1;
  flags.beam[1] = generator()->eventHandler()->incoming().second->id() > 0 ? 1 : -1;

  flags.Ebeam = generator()->maximumCMEnergy()/2./GeV;
  flags.partonptmax = flags.Ebeam;

  if ( theExtremalPTMin > theJetPTMin )
    throw Exception() << "Minimum jet pt in Cuts needs to be smaller than HEJFactory:JetPTMin"
		      << Exception::abortnow;

  flags.jetptmin = theJetPTMin/GeV;
  flags.extpartonptmin =  theExtremalPTMin/GeV;

  flags.scale = theFixedScale/GeV;
  flags.scalesetting = theScaleChoice;

  /*
  cerr << "HEJFactory settings:\n"
       << "flags.Ebeam = " << flags.Ebeam << "\n"
       << "flags.jetptmin = " << flags.jetptmin << "\n"
       << "flags.extpartonptmin = " << flags.extpartonptmin << "\n"
       << "flags.scalesetting = " << flags.scalesetting << "\n"
       << flush;
  */

  return flags;

}

void HEJFactory::makeCMultijet(tStdXCombPtr xc) {

  theJets = new CMultijet(setup(xc));

  for ( MEVector::iterator m = MEs().begin();
	m != MEs().end(); ++m )
    dynamic_cast<HEJMEBase&>(**m).jets(theJets);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void HEJFactory::doinit() {
  SubProcessHandler::doinit();

  MEs().clear();

  for ( unsigned int n = 0; n <= theNGluons; ++n ) {

    for ( PDVector::const_iterator pi = theIncomingFlavours.begin();
	  pi != theIncomingFlavours.end(); ++pi ) {
      for ( PDVector::const_iterator pj = pi;
	    pj != theIncomingFlavours.end(); ++pj ) {
	makeME(*pi,*pj,n);
	if ( *pi != *pj )
	  makeME(*pj,*pi,n);
      }
    }

  }

}

void HEJFactory::persistentOutput(PersistentOStream & os) const {
  os << theHEJME << theIncomingFlavours << theNGluons
     << ounit(theJetPTMin,GeV) << ounit(theExtremalPTMin,GeV) 
     << theScaleChoice << ounit(theFixedScale,GeV);
}

void HEJFactory::persistentInput(PersistentIStream & is, int) {
  is >> theHEJME >> theIncomingFlavours >> theNGluons
     >> iunit(theJetPTMin,GeV) >> iunit(theExtremalPTMin,GeV) 
     >> theScaleChoice >> iunit(theFixedScale,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HEJFactory,SubProcessHandler>
  describeHerwigHEJFactory("Herwig::HEJFactory", "HwMatchbox.so HwHEJ.so");

void HEJFactory::Init() {

  static ClassDocumentation<HEJFactory> documentation
    ("HEJFactory sets up HEJ matrix elements.");


  static Reference<HEJFactory,HEJMEBase> interfaceHEJME
    ("HEJME",
     "The HEJ matrix element to consider",
     &HEJFactory::theHEJME, false, false, true, false, false);


  static RefVector<HEJFactory,ParticleData> interfaceIncomingFlavours
    ("IncomingFlavours",
     "The incoming flavours to take into account.",
     &HEJFactory::theIncomingFlavours, -1, false, false, true, false, false);


  static Parameter<HEJFactory,unsigned int> interfaceNGluons
    ("NGluons",
     "The number of additional gluons to consider.",
     &HEJFactory::theNGluons, 20, 0, 0,
     false, false, Interface::lowerlim);


  static Parameter<HEJFactory,Energy> interfaceJetPTMin
    ("JetPTMin",
     "The minimum jet pt of the two hardest jets.",
     &HEJFactory::theJetPTMin, GeV, 20.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<HEJFactory,Energy> interfaceExtremalPTMin
    ("ExtremalPTMin",
     "The minimum jet pt of the extremal partons.",
     &HEJFactory::theExtremalPTMin, GeV, 15.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<HEJFactory,int> interfaceScaleChoice
    ("ScaleChoice",
     "The scale setting to use.",
     &HEJFactory::theScaleChoice, 1, false, false);
  static SwitchOption interfaceScaleChoiceFixedScale
    (interfaceScaleChoice,
     "FixedScale",
     "Use a fixed scale.",
     0);
  static SwitchOption interfaceScaleChoiceHardestPt
    (interfaceScaleChoice,
     "HardestPt",
     "Use the pt of the hardest jet.",
     1);
  static SwitchOption interfaceScaleChoiceGeometricMean
    (interfaceScaleChoice,
     "GeometricMean",
     "Use the geometric mean scale choice.",
     2);



  static Parameter<HEJFactory,Energy> interfaceFixedScale
    ("FixedScale",
     "The fixed scale value, if chosen.",
     &HEJFactory::theFixedScale, GeV, 40.0*GeV, 0*GeV, 0*GeV,
     false, false, Interface::lowerlim);


}

