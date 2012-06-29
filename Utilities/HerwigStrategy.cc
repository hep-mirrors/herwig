// -*- C++ -*-
//
// HerwigStrategy.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2008-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HerwigStrategy class.
//

#include "HerwigStrategy.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ParticleData.h"

// sets the static version info string 'HerwigStrategy::version'
#include "versionstring.h"

using namespace Herwig;

IBPtr HerwigStrategy::clone() const {
  return new_ptr(*this);
}

IBPtr HerwigStrategy::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<HerwigStrategy> HerwigStrategy::initHerwigStrategy;

void HerwigStrategy::Init() {
  static ClassDocumentation<HerwigStrategy> interfaceDescription
    ("The default strategy for Herwig.",
     "Herwig++~\\cite{Bahr:2008pv}", 
     "%\\cite{Bahr:2008pv}\n"
     "\\bibitem{Bahr:2008pv}\n"
     "  M.~Bahr {\\it et al.},\n"
     "  ``Herwig++ Physics and Manual,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 58} (2008) 639\n"
     "  [arXiv:0803.0883 [hep-ph]].\n"
     "  %%CITATION = EPHJA,C58,639;%%\n"
     );
}

