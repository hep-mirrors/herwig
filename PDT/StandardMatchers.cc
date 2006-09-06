// -*- C++ -*-
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
}

}
