// -*- C++ -*-
//
// MatchboxNLOME.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxNLOME class.
//

#include "MatchboxNLOME.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

MatchboxNLOME::MatchboxNLOME() 
  : MEBase(), theNDim(-1), theCheckPoles(false) {}

MatchboxNLOME::~MatchboxNLOME() {}


IBPtr MatchboxNLOME::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxNLOME::fullclone() const {
  return new_ptr(*this);
}

void MatchboxNLOME::cloneDependencies(const std::string& prefix) {

  Ptr<MatchboxMEBase>::ptr myMatrixElement = theMatrixElement->cloneMe();
  ostringstream pname;
  pname << (prefix == "" ? fullName() : prefix) << "/" << myMatrixElement->name();
  if ( ! (generator()->preinitRegister(myMatrixElement,pname.str()) ) )
    throw InitException() << "Matrix element " << pname.str() << " already existing.";
  myMatrixElement->cloneDependencies(pname.str());
  theMatrixElement = myMatrixElement;

  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    Ptr<MatchboxInsertionOperator>::ptr myIOP = (**v).cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << (**v).name();
    if ( ! (generator()->preinitRegister(myIOP,pname.str()) ) )
      throw InitException() << "Insertion operator " << pname.str() << " already existing.";
    *v = myIOP;
    (**v).setBorn(theMatrixElement);
  }

}

bool MatchboxNLOME::generateKinematics(const double * r) {
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v )
    if ( (**v).nDimAdditional() )
      (**v).additionalKinematics(r+theMatrixElement->nDim());
  bool ret = matrixElement()->generateKinematics(r);
  jacobian(matrixElement()->lastXComb().jacobian());
  return ret;
}

void MatchboxNLOME::logPoles() const {
  double res2me = matrixElement()->oneLoopDoublePole();
  double res1me = matrixElement()->oneLoopSinglePole();
  double res2i = 0.;
  double res1i = 0.;
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    res2i += (**v).oneLoopDoublePole();
    res1i += (**v).oneLoopSinglePole();
  }
  double diff2 = abs(res2me) != 0. ? 1.-abs(res2i/res2me) : abs(res2i)-abs(res2me);
  double diff1 = abs(res1me) != 0. ? 1.-abs(res1i/res1me) : abs(res1i)-abs(res1me);
  generator()->log() 
    << "check "
    << log10(abs(diff2)) << " " << log10(abs(diff1)) << "\n"
    << flush;
}

double MatchboxNLOME::me2() const {
  if ( !matrixElement()->onlyOneLoop() && checkPoles() )
    logPoles();
  double res = !matrixElement()->onlyOneLoop() ? matrixElement()->me2() : 0.;
  if ( matrixElement()->haveOneLoop() ) {
    res += 
      matrixElement()->oneLoopInterference();
  }
  if ( !matrixElement()->onlyOneLoop() )
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	    virtuals().begin(); v != virtuals().end(); ++v )
      res += (**v).me2();
  return res;
}

CrossSection MatchboxNLOME::dSigHatDR() const {
  CrossSection res = 
    !matrixElement()->onlyOneLoop() ? matrixElement()->dSigHatDR() : ZERO;
  if ( res == ZERO && !matrixElement()->onlyOneLoop() )
    return res;
  if ( !matrixElement()->onlyOneLoop() && checkPoles() )
    logPoles();
  if ( matrixElement()->haveOneLoop() ) {
    if ( matrixElement()->onlyOneLoop() )
      matrixElement()->getPDFWeight();
    double vme =
      matrixElement()->oneLoopInterference();
    res +=
      sqr(hbarc) * vme *
      jacobian() * lastMEPDFWeight() /
      (2.*lastSHat());
  }
  if ( !matrixElement()->onlyOneLoop() )
    for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	    virtuals().begin(); v != virtuals().end(); ++v )
      res += (**v).dSigHatDR();
  lastMECrossSection(res);
  return res;
}

void MatchboxNLOME::flushCaches() {
  theMatrixElement->flushCaches();
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v )
    (**v).flushCaches();
}

void MatchboxNLOME::setXComb(tStdXCombPtr xc) {

  MEBase::setXComb(xc);
  theMatrixElement->setXComb(xc);
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v )
    (**v).setXComb(xc);
}

void MatchboxNLOME::getNDim() const {
  int maxadd = 0;
  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    if ( (**v).nDimAdditional() > 1 ) {
      throw InitException() << "at most one additional random number supported for "
			    << "virtual corrections currently";
    }
    maxadd = max(maxadd,(**v).nDimAdditional());
  }
  theNDim = theMatrixElement->nDim() + maxadd;
}

void MatchboxNLOME::print(ostream& os) const {

  os << "--- MatchboxNLOME setup --------------------------------------------------------\n";

  os << " '" << name() << "' using\n"
     << " matrix element '" << matrixElement()->name() << "' and insertion operators:\n";

  for ( vector<Ptr<MatchboxInsertionOperator>::ptr>::const_iterator v =
	  virtuals().begin(); v != virtuals().end(); ++v ) {
    os << " '" << (**v).name() << "' with " 
       << ((**v).isDR() ? "" : "C") << "DR/";
    if ( (**v).isCS() )
      os << "CS";
    if ( (**v).isBDK() )
      os << "BDK";
    if ( (**v).isExpanded() )
      os << "expanded";
    os << " conventions\n";
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void MatchboxNLOME::printLastEvent(ostream& os) const {

  os << "--- MatchboxNLOME last event information ---------------------------------------\n";

  os << " for matrix element '" << name() << "'\n";

  os << " process considered:\n ";

  int in = 0;
  for ( cPDVector::const_iterator p = mePartonData().begin();
	p != mePartonData().end(); ++p ) {
    os << (**p).PDGName() << " ";
    if ( ++in == 2 )
      os << " -> ";
  }

  os << " kinematic environment as set by the XComb " << lastXCombPtr() << ":\n"
     << " sqrt(shat)/GeV = " << sqrt(lastSHat()/GeV2)
     << " x1 = " << lastX1() << " x2 = " << lastX2() 
     << " alphaS = " << lastAlphaS() << "\n";

  os << " momenta/GeV generated from random numbers\n ";
  copy(meInfo().begin(),meInfo().end(),ostream_iterator<double>(os," "));
  os << ":\n ";

  for ( vector<Lorentz5Momentum>::const_iterator p = meMomenta().begin();
	p != meMomenta().end(); ++p ) {
    os << (*p/GeV) << "\n ";
  }

  os << "last cross section/nb calculated was:\n "
     << (lastMECrossSection()/nanobarn) << " (pdf weight " << lastMEPDFWeight() << ")\n";

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void MatchboxNLOME::doinit() {
  MEBase::doinit();
}

void MatchboxNLOME::doinitrun() {
  MEBase::doinitrun();
}

void MatchboxNLOME::persistentOutput(PersistentOStream & os) const {
  os << theMatrixElement << theVirtuals << theNDim << theCheckPoles;
}

void MatchboxNLOME::persistentInput(PersistentIStream & is, int) {
  is >> theMatrixElement >> theVirtuals >> theNDim >> theCheckPoles;
}

void MatchboxNLOME::Init() {

  static ClassDocumentation<MatchboxNLOME> documentation
    ("MatchboxNLOME");


  static Reference<MatchboxNLOME,MatchboxMEBase> interfaceBornME
    ("BornME",
     "The Born matrix element",
     &MatchboxNLOME::theMatrixElement, false, false, true, false, false);


  static RefVector<MatchboxNLOME,MatchboxInsertionOperator> interfaceVirtuals
    ("Virtuals",
     "The virtual corrections to be added.",
     &MatchboxNLOME::theVirtuals, -1, false, false, true, true, false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxNLOME,MEBase>
describeMatchboxNLOME("Herwig::MatchboxNLOME", "HwMatchbox.so");

