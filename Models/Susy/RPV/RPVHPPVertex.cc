// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVHPPVertex class.
//

#include "RPVHPPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVHPPVertex::RPVHPPVertex() {
  orderInGs(0);
  orderInGem(3);
}

IBPtr RPVHPPVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVHPPVertex::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void RPVHPPVertex::persistentOutput(PersistentOStream & os) const {
}

void RPVHPPVertex::persistentInput(PersistentIStream & is, int) {
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<RPVHPPVertex,VVSLoopVertex>
describeHerwigRPVHPPVertex("Herwig::RPVHPPVertex", "HwRPV.so");

void RPVHPPVertex::Init() {

  static ClassDocumentation<RPVHPPVertex> documentation
    ("There is no documentation for the RPVHPPVertex class");

}

void RPVHPPVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			       tcPDPtr particle3) {


}
