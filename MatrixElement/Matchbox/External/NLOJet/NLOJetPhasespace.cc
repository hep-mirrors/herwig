// -*- C++ -*-
//
// NLOJetPhasespace.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetPhasespace class.
//

#include "NLOJetPhasespace.h"

using namespace Herwig;

//DescribeClass<NLOJetPhasespace<N,I,F>,MatchboxPhasespace>
//describeNLOJetPhasespace("Herwig::NLOJetPhasespace<N,I,F>", "HwMatchboxNLOJet.so");

DescribeClass<NLOJetPhasespace<0,0,0>,MatchboxPhasespace>
describeNLOJetPhasespace_n0i0f0("Herwig::NLOJetPhasespace<0,0,0>", "HwMatchboxNLOJet.so");

DescribeClass<NLOJetPhasespace<0,2,0>,MatchboxPhasespace>
describeNLOJetPhasespace_n0i2f0("Herwig::NLOJetPhasespace<0,2,0>", "HwMatchboxNLOJet.so");

/*
DescribeClass<NLOJetPhasespace<1,1,0>,MatchboxPhasespace>
describeNLOJetPhasespace_n1i1f0("Herwig::NLOJetPhasespace<1,1,0>", "HwMatchboxNLOJet.so");
*/

DescribeClass<NLOJetPhasespace<2,2,0>,MatchboxPhasespace>
describeNLOJetPhasespace_n2i2f0("Herwig::NLOJetPhasespace<2,2,0>", "HwMatchboxNLOJet.so");

