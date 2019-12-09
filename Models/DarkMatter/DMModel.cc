// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMModel class.
//

#include "DMModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DMModel::DMModel() {}

IBPtr DMModel::clone() const {
  return new_ptr(*this);
}

IBPtr DMModel::fullclone() const {
  return new_ptr(*this);
}

void DMModel::persistentOutput(PersistentOStream & os) const {
}

void DMModel::persistentInput(PersistentIStream & is, int) {
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMModel,BSMModel>
describeHerwigDMModel("Herwig::DMModel", "HwDMModel.so");

void DMModel::Init() {

  static ClassDocumentation<DMModel> documentation
    ("The DMModel class is designed to implement a simple dark matter model"
     " with fermionic dark matter and a vector mediator, "
     "as described in  arXiv:1911.11147",
     "The DMModel class is designed to implement a simple dark matter model"
     " with fermionic dark matter and a vector mediator, "
     "as described in \\cite{Plehn:2019jeo}",
     "\\bibitem{Plehn:2019jeo}"
     "T.~Plehn, P.~Reimitz and P.~Richardson,"
     "%``Hadronic Footprint of GeV-Mass Dark Matter,''"
     "arXiv:1911.11147 [hep-ph]."
     "%%CITATION = ARXIV:1911.11147;%%");

}

