// -*- C++ -*-
//
// ME2byDipoles.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ME2byDipoles class.
//

#include "ME2byDipoles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ME2byDipoles::ME2byDipoles() 
  : MatchboxReweightBase() {}

ME2byDipoles::~ME2byDipoles() {}

IBPtr ME2byDipoles::clone() const {
  return new_ptr(*this);
}

IBPtr ME2byDipoles::fullclone() const {
  return new_ptr(*this);
}

void ME2byDipoles::setup(Ptr<SubtractedME>::tptr sub) {

  theRealME = dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(sub->head());
  assert(theRealME);
  Ptr<MatchboxMEBase>::ptr nreal = theRealME->cloneMe();
  ostringstream pname;
  pname << fullName() << "/" << nreal->name();
  if ( ! (generator()->preinitRegister(nreal,pname.str()) ) )
    throw InitException() << "Matrix element " << pname.str() << " already existing.";
  nreal->cloneDependencies();
  theRealME = nreal;

  vector<Ptr<SubtractionDipole>::ptr> inDipoles = sub->dipoles();
  theDipoles.clear();

  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d =
	  inDipoles.begin(); d != inDipoles.end(); ++d ) {
    Ptr<SubtractionDipole>::ptr ndipole = (**d).cloneMe();
    ostringstream dname;
    dname << fullName() << "/" << (**d).name();
    if ( ! (generator()->preinitRegister(ndipole,dname.str())) )
      throw InitException() << "Dipole '" << dname.str() << "' already existing.";
    ndipole->cloneDependencies();
    theDipoles.push_back(ndipole);
  }

  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d =
	  theDipoles.begin(); d != theDipoles.end(); ++d )
    (**d).doSubtraction();

}

void ME2byDipoles::setup(Ptr<SubtractionDipole>::tptr dip, 
			 Ptr<SubtractedME>::tptr sub) {

  vector<Ptr<SubtractionDipole>::ptr> inDipoles = sub->dipoles();

  theDipoles.clear();

  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d =
	  inDipoles.begin(); d != inDipoles.end(); ++d ) {
    Ptr<SubtractionDipole>::ptr ndipole = (**d).cloneMe();
    ostringstream dname;
    dname << fullName() << "/" << (**d).name();
    if ( ! (generator()->preinitRegister(ndipole,dname.str())) )
      throw InitException() << "Dipole '" << dname.str() << "' already existing.";
    ndipole->cloneDependencies();
    theDipoles.push_back(ndipole);
    if ( *d == dip ) {
      projectionDipole(theDipoles.back());
    }
  }

  theRealME = Ptr<MatchboxMEBase>::ptr();

  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d =
	  theDipoles.begin(); d != theDipoles.end(); ++d )
    (**d).doSubtraction();

}

void ME2byDipoles::setXComb(tStdXCombPtr real) {

  assert(real);
  theLastXComb = real;

  if ( theRealME )
    theRealME->setXComb(theLastXComb);

  map<StdXCombPtr,vector<StdXCombPtr> >::iterator xcs
    = theXCombMap.find(theLastXComb);

  if ( xcs == theXCombMap.end() )
    getXCombs(theLastXComb);

}

double ME2byDipoles::scaledBornScreen() const {

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' evaluating Born screening\n";

  Energy scale = projectionDipole()->lastDipoleScale();
  Energy pt = projectionDipole()->lastPt();

  if ( projectionDipole()->verbose() )
    generator()->log() << "from pt/GeV = " << (pt/GeV) 
		       << " scale/GeV = " << (scale/GeV)
		       << "\n" << flush;

  return pow(pt/scale,4.);

}

double ME2byDipoles::scaledBorn(Energy2 factorizationScale) const {

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' evaluating scaled Born\n" << flush;

  projectionDipole()->underlyingBornME()->setScale();
  projectionDipole()->underlyingBornME()->getPDFWeight(factorizationScale);
  double me2 = projectionDipole()->underlyingBornME()->me2();
  double pdf = projectionDipole()->underlyingBornME()->lastXComb().lastMEPDFWeight();

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' done evaluating scaled Born\n" << flush;

  return me2 * pdf;

}

void ME2byDipoles::flushCaches() {
  if ( theRealME )
    theRealME->flushCaches();
  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator d =
	  theDipoles.begin(); d != theDipoles.end(); ++d )
    (**d).flushCaches();

}

double ME2byDipoles::evaluate(double& sratio, bool doSRatio) const {

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' ME2byDipoles evaluating\n" << flush;

  double den = 0.0;
  sratio = 0.;

  double numDip = 0.0;

  map<StdXCombPtr,vector<StdXCombPtr> >::const_iterator xcs
    = theXCombMap.find(theLastXComb);
  assert(xcs != theXCombMap.end());

  for ( vector<StdXCombPtr>::const_iterator xcit = xcs->second.begin(); xcit != xcs->second.end(); ++xcit ) {
    (**xcit).clean();
    Ptr<SubtractionDipole>::tptr dip = 
      dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>((**xcit).matrixElement());
    assert(dip);
    dip->setXComb(*xcit);
    if ( !dip->apply() )
      continue;
    tStdXCombPtr depXComb = dip->lastXCombPtr();  
    if ( !dip->generateTildeKinematics() )
      continue;
    depXComb->setIncomingPartons();

    dip->realEmissionME()->setScale();
    dip->underlyingBornME()->setScale();
    double res = dip->me2Avg(-dip->underlyingBornME()->me2());
    den += res;
    if ( doSRatio )
      if ( depXComb->willPassCuts() )
	sratio += dip->me2();

    if ( dip == projectionDipole() ) {
      numDip = res;
    }
  }

  if ( sratio != 0. ) {
    assert(abs(den) != 0.0);
    sratio /= den;
  }

  if ( theRealME ) {
    if ( !realME()->lastXCombPtr()->willPassCuts() ) {
      if ( projectionDipole()->verbose() )
	generator()->log() << "real emission did not pass the cuts\n" << flush;
      return 0.0;
    }
  }

  assert(abs(den) != 0.);
  double num = theRealME ? realME()->me2() : numDip;

  double res = num / den;

  if ( projectionDipole()->verbose() ) {
    generator()->log() << "'" << name() << "' done evaluating\n"
		       << "numerator = " << num << " denominator = "
		       << den << "\n" << flush;
  }

  return res;

}

void ME2byDipoles::getXCombs(tStdXCombPtr xc) {

  vector<StdXCombPtr> xcs;

  for ( vector<Ptr<SubtractionDipole>::ptr>::iterator dip = theDipoles.begin();
	dip != theDipoles.end(); ++dip ) {
    StdXCombPtr depxc = (**dip).makeBornXComb(xc);
    xcs.push_back(depxc);
  }

  theXCombMap[xc] = xcs;

}

void ME2byDipoles::print(ostream& os) const {

  os << "--- ME2byDipoles setup ---------------------------------------------------------\n";

  os << " '"  << name() << "'\n"
     << " real emission matrix element '" << theRealME->name() << "'\n"
     << " projection dipole: '"
     << (projectionDipole() ? projectionDipole()->name() : "")
     << "'\n";

  os << " associated dipoles are:\n";

  for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator dip =
	  theDipoles.begin(); dip != theDipoles.end(); ++dip ) {
    os << " '" << (**dip).name() << "'\n";
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void ME2byDipoles::printLastEvent(ostream& os) const {

  os << "--- ME2byDipoles last event information ----------------------------------------\n";

  os << " for ratio '" << name() << "'\n";

  os << " real emission event information:\n";

  if ( dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(theRealME) )
    dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(theRealME)->printLastEvent(os);
  else if ( dynamic_ptr_cast<Ptr<SubtractionDipole>::ptr>(theRealME) )
    dynamic_ptr_cast<Ptr<SubtractionDipole>::ptr>(theRealME)->printLastEvent(os);
  else
    os << " unknown MEBase object.\n";

  if ( projectionDipole() ) {
    os << " projection dipole event information:\n";
    projectionDipole()->printLastEvent(os);
  }

  os << " dipoles event information:\n";
  for ( vector<Ptr<SubtractionDipole>::ptr>::const_iterator d = theDipoles.begin();
	d != theDipoles.end(); ++d )
    (**d).printLastEvent(os);


  os << "--- end ME2byDipoles last event information ------------------------------------\n";

  os << flush;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ME2byDipoles::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theRealME << theProjectionDipole << theDipoles;
}

void ME2byDipoles::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theRealME >> theProjectionDipole >> theDipoles;
}

void ME2byDipoles::Init() {

  static ClassDocumentation<ME2byDipoles> documentation
    ("ME2byDipoles");


  static Reference<ME2byDipoles,MatchboxMEBase> interfaceRealME
    ("RealME",
     "The real emission matrix element.",
     &ME2byDipoles::theRealME, false, false, true, false, false);

  static Reference<ME2byDipoles,SubtractionDipole> interfaceProjectionDipole
    ("ProjectionDipole",
     "The projection dipole.",
     &ME2byDipoles::theProjectionDipole, false, false, true, false, false);


  static RefVector<ME2byDipoles,SubtractionDipole> interfaceDipoles
    ("Dipoles",
     "The dipoles associated to the real emission matrix element.",
     &ME2byDipoles::theDipoles, -1, false, false, true, false, false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<ME2byDipoles,MatchboxReweightBase>
describeME2byDipoles("Herwig::ME2byDipoles", "HwMatchbox.so");
