// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoGGGSplitFun class.
//

#include "GtoGGGSplitFun.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


GtoGGGSplitFun::~GtoGGGSplitFun() {}


AbstractClassDescription<GtoGGGSplitFun> GtoGGGSplitFun::initGtoGGGSplitFun;
// Definition of the static class description member.


void GtoGGGSplitFun::Init() {

  static ClassDocumentation<GtoGGGSplitFun> documentation
    ("This abstract class provides the common functionality for Initial State and ",
     "Final State  G->GGG  splitting functions.");

}

