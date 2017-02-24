// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BlobME class.
//

#include "BlobME.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

BlobME::BlobME() 
  : theNAdditional(0) {}

BlobME::~BlobME() {}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void BlobME::persistentOutput(PersistentOStream & os) const {
  os << thePhasespace << theNAdditional;
}

void BlobME::persistentInput(PersistentIStream & is, int) {
  is >> thePhasespace >> theNAdditional;
}

void BlobME::setXComb(tStdXCombPtr xc) {
  BlobMEBase::setXComb(xc);
  thePhasespace->setXComb(xc);
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<BlobME,BlobMEBase>
  describeHerwigBlobME("Herwig::BlobME", "Herwig.so");

void BlobME::Init() {

  static ClassDocumentation<BlobME> documentation
    ("BlobME serves as a base class for special processes such as "
     "instanton or sphaleron induced ones.");

  static Parameter<BlobME,size_t> interfaceNAdditional
    ("NAdditional",
     "The number of additional objects to consider.",
     &BlobME::theNAdditional, 0, 0, 0,
     false, false, Interface::lowerlim);

  static Reference<BlobME,MatchboxPhasespace> interfacePhasespace
    ("Phasespace",
     "The phase space to use.",
     &BlobME::thePhasespace, false, false, true, false, false);

}

