// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwMEBase class.
//

#include "HwMEBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

AbstractNoPIOClassDescription<HwMEBase> HwMEBase::initHwMEBase;
// Definition of the static class description member.

void HwMEBase::Init() {

  static ClassDocumentation<HwMEBase> documentation
    ("The HwMEBase class provides the base class for the"
     " implementation of matrix elements in Herwig++");

}

