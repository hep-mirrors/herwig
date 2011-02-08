// -*- C++ -*-
//
// StandardMatchers.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//  Ensures the StandardMatchers get created
//
#include "ThePEG/PDT/Matcher.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "StandardMatchers.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "Matcher.tcc"
#endif

using namespace Herwig;
using namespace ThePEG;

namespace {

void dummy() {

  static MatchPhoton m00;
  static MatchTop    m01;
  static MatchHadron m02;
  static MatchWBoson m03;
  static MatchZBoson m04;
  static MatchHiggsBoson m05;
}

}
