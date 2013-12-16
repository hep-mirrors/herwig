// -*- C++ -*-
//
// NLOJetMEBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEBase class.
//

#include "NLOJetMEBase.h"

using namespace Herwig;

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<NLOJetMEBase<0,0,0>,Herwig::MatchboxMEBase>
describeNLOJetMEBase_n0i0f0("Herwig::NLOJetMEBase<0,0,0>", "HwMatchboxNLOJet.so");

DescribeAbstractClass<NLOJetMEBase<0,2,0>,Herwig::MatchboxMEBase>
describeNLOJetMEBase_n0i2f0("Herwig::NLOJetMEBase<0,2,0>", "HwMatchboxNLOJet.so");

DescribeAbstractClass<NLOJetMEBase<2,2,0>,Herwig::MatchboxMEBase>
describeNLOJetMEBase_n2i2f0("Herwig::NLOJetMEBase<2,2,0>", "HwMatchboxNLOJet.so");

