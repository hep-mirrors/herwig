// -*- C++ -*-
//
// PostClustering.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PostClustering class.
//

#include "PostClustering.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PostClustering.tcc"
#endif


using namespace Herwig;

PostClustering::~PostClustering() {}

AbstractNoPIOClassDescription<PostClustering> PostClustering::initPostClustering;
// Definition of the static class description member.

void PostClustering::Init() {

  static ClassDocumentation<PostClustering> documentation
    ("PostClustering defines an interface for things to be applied to"
     " all particles in a cascade reconstruction after a specific branching"
     " has been performed.");

}

