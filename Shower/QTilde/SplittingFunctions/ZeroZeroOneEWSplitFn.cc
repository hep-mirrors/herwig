// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZeroZeroOneEWSplitFn class.
//

#include "ZeroZeroOneEWSplitFn.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Models/StandardModel/SMFFHVertex.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

IBPtr ZeroZeroOneEWSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr ZeroZeroOneEWSplitFn::fullclone() const {
  return new_ptr(*this);
}

void ZeroZeroOneEWSplitFn::persistentOutput(PersistentOStream & os) const {
  os << _couplingValueIm << _couplingValueRe;
}

void ZeroZeroOneEWSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> _couplingValueIm >> _couplingValueRe;
}

// The following static variable is needed for the type description system in ThePEG.
DescribeClass<ZeroZeroOneEWSplitFn,Sudakov1to2FormFactor>
describeHerwigZeroZeroOneEWSplitFn("Herwig::ZeroZeroOneEWSplitFn", "HwShower.so");


void ZeroZeroOneEWSplitFn::Init() {

  static ClassDocumentation<ZeroZeroOneEWSplitFn> documentation
    ("The ZeroZeroOneEWSplitFn class implements purly beyond SM electroweak splittings H->H'V");

  static Parameter<ZeroZeroOneEWSplitFn, double> interfaceCouplingValueIm
    ("CouplingValue.Im",
     "The numerical value (imaginary part) of the splitting coupling to be imported for BSM splittings",
     &ZeroZeroOneEWSplitFn::_couplingValueIm, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

  static Parameter<ZeroZeroOneEWSplitFn, double> interfaceCouplingValueRe
    ("CouplingValue.Re",
     "The numerical value (real part) of the splitting coupling to be imported for BSM splittings",
     &ZeroZeroOneEWSplitFn::_couplingValueRe, 0.0, -1.0E6, +1.0E6,
     false, false, Interface::limited);

}
