// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVHGGVertex class.
//

#include "RPVHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVHGGVertex::RPVHGGVertex() {
  orderInGs(2);
  orderInGem(1);
}

IBPtr RPVHGGVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVHGGVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVHGGVertex::persistentOutput(PersistentOStream & os) const {
}

void RPVHGGVertex::persistentInput(PersistentIStream & is, int) {
}


// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<RPVHGGVertex,VVSLoopVertex>
describeHerwigRPVHGGVertex("Herwig::RPVHGGVertex", "HwRPV.so");

void RPVHGGVertex::Init() {

  static ClassDocumentation<RPVHGGVertex> documentation
    ("There is no documentation for the RPVHGGVertex class");

}

void RPVHGGVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			       tcPDPtr particle3) {
}
