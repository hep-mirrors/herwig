// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeMECorrection class.
//

#include "QTildeMECorrection.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"

using namespace Herwig;

AbstractNoPIOClassDescription<QTildeMECorrection> QTildeMECorrection::initQTildeMECorrection;
// Definition of the static class description member.

void QTildeMECorrection::Init() {

  static ClassDocumentation<QTildeMECorrection> documentation
    ("There is no documentation for the QTildeMECorrection class");

}

