// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NonBShowerVeto class.
//

#include "NonBShowerVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

IBPtr NonBShowerVeto::clone() const {
  return new_ptr(*this);
}

IBPtr NonBShowerVeto::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeNoPIOClass<NonBShowerVeto,FullShowerVeto>
  describeHerwigNonBShowerVeto("Herwig::NonBShowerVeto", "HwShowerVeto.so");

void NonBShowerVeto::Init() {

  static ClassDocumentation<NonBShowerVeto> documentation
    ("The NonBShowerVeto class vetos the parton-shower when no b (anti)quarks have been produced");

}

bool NonBShowerVeto::vetoShower() {
  // loop over final-state
  for(vector<tPPtr>::const_iterator it=finalState().begin(); it!=finalState().end();++it) {
    // don't veto if find a b (anti)quark
    if(abs((**it).id())==5) {
      return false;
    }
  }
  // no b (anti)quarks veto shower
  return true;
}
