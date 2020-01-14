// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KMatrix class.
//

#include "KMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KMatrix::KMatrix(FlavourInfo flavour, vector<Channels> channels,
		 vector<Energy> poles) : flavour_(flavour), channels_(channels), poles_(poles)
{}

void KMatrix::persistentOutput(PersistentOStream & os) const {
  os << ounit(poles_,GeV);
}

void KMatrix::persistentInput(PersistentIStream & is, int) {
  is >> iunit(poles_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<KMatrix,Interfaced>
describeHerwigKMatrix("Herwig::KMatrix", "Herwig.so");

void KMatrix::Init() {

  static ClassDocumentation<KMatrix> documentation
    ("The KMatrix class provides a base class for the implementation of "
     "K-matrix parameterizations in Herwig");

}

