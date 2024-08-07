// -*- C++ -*-
//
// GluonMassGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GluonMassGenerator class.
//

#include "GluonMassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ClusterHadronizationHandler.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr GluonMassGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr GluonMassGenerator::fullclone() const {
  return new_ptr(*this);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void GluonMassGenerator::persistentOutput(PersistentOStream &) const {}

void GluonMassGenerator::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<GluonMassGenerator,HandlerBase>
  describeHerwigGluonMassGenerator("Herwig::GluonMassGenerator", "Herwig.so");

void GluonMassGenerator::Init() {

  static ClassDocumentation<GluonMassGenerator> documentation
    ("Dynamic gluon mass generation");

}

