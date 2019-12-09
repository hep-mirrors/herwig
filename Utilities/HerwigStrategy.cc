// -*- C++ -*-
//
// HerwigStrategy.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2008-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HerwigStrategy class.
//

#include "HerwigStrategy.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/Repository.h"

using namespace Herwig;

const string HerwigStrategy::version = 
#include "hgstamp.inc"
"";

const std::string HerwigStrategy::versionstring() const {
      return HerwigStrategy::version + " / " + Repository::version();
}


IBPtr HerwigStrategy::clone() const {
  return new_ptr(*this);
}

IBPtr HerwigStrategy::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<HerwigStrategy,ThePEG::Strategy>
describeHerwigHerwigStrategy("Herwig::HerwigStrategy", "Herwig.so");

void HerwigStrategy::Init() {
  static ClassDocumentation<HerwigStrategy> interfaceDescription
    ("The default strategy for Herwig.",
     "Herwig~\\cite{Bahr:2008pv}", 
     "\\bibitem{Bahr:2008pv}\n"
     "  M.~Bahr {\\it et al.},\n"
     "  ``Herwig Physics and Manual,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 58} (2008) 639\n"
     "  [arXiv:0803.0883 [hep-ph]].\n"
     "  %%CITATION = EPHJA,C58,639;%%\n"
     );
}

