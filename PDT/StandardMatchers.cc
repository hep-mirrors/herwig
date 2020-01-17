// -*- C++ -*-
//
// StandardMatchers.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
//  Ensures the StandardMatchers get created
//
#include "ThePEG/PDT/Matcher.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "StandardMatchers.h"

using namespace Herwig;

#define THEPEG_MATCH_DESC(T)                           \
/**                                                    \
 * This template specialization registers the Matcher  \
 */                                                    \
template <>                                            \
NoPIOClassDescription<T> T::initMatcher                \
 = NoPIOClassDescription<T>();                         \



namespace ThePEG {
  THEPEG_MATCH_DESC(MatchPhoton)
  THEPEG_MATCH_DESC(MatchBottom)
  THEPEG_MATCH_DESC(MatchTop)
  THEPEG_MATCH_DESC(MatchHadron)
  THEPEG_MATCH_DESC(MatchWBoson)
  THEPEG_MATCH_DESC(MatchZBoson)
  THEPEG_MATCH_DESC(MatchHiggsBoson)
  THEPEG_MATCH_DESC(MatchChargedLepton)
  THEPEG_MATCH_DESC(MatchLightParticle)
}

namespace {
  using namespace ThePEG;

  static MatchPhoton m00;
  static MatchBottom m01;
  static MatchTop    m02;
  static MatchHadron m03;
  static MatchWBoson m04;
  static MatchZBoson m05;
  static MatchHiggsBoson m06;
  static MatchChargedLepton m07;
  static MatchLightParticle m08;

}
