// -*- C++ -*-
//
// StandardMatchers.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//  Ensures the StandardMatchers get created
//
#include "ThePEG/PDT/Matcher.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "StandardMatchers.h"

using namespace Herwig;
using namespace ThePEG;

namespace {

void dummy() {

  static MatchPhoton m00;
  static MatchBottom m01;
  static MatchTop    m02;
  static MatchHadron m03;
  static MatchWBoson m04;
  static MatchZBoson m05;
  static MatchHiggsBoson m06;
  static MatchChargedLepton m07;
}

}
