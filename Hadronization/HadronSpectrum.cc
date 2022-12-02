// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HadronSpectrum class.
//

#include "HadronSpectrum.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Interface/RefVector.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HadronSpectrum::HadronSpectrum() 
  : Interfaced() {}

HadronSpectrum::~HadronSpectrum() {}

void HadronSpectrum::doinit() {
  Interfaced::doinit();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HadronSpectrum::persistentOutput(PersistentOStream & os) const {
  os << _table << _partons << _forbidden;
}

void HadronSpectrum::persistentInput(PersistentIStream & is, int) {
  is >> _table >> _partons >> _forbidden;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<HadronSpectrum,Interfaced>
  describeHerwigHadronSpectrum("Herwig::HadronSpectrum", "Herwig.so");

void HadronSpectrum::Init() {

  static ClassDocumentation<HadronSpectrum> documentation
    ("There is no documentation for the HadronSpectrum class");

  static RefVector<HadronSpectrum,ParticleData> interfacePartons
    ("Partons",
     "The partons which are to be considered as the consistuents of the hadrons.",
     &HadronSpectrum::_partons, -1, false, false, true, false, false);

  static RefVector<HadronSpectrum,ParticleData> interfaceForbidden
    ("Forbidden",
     "The PDG codes of the particles which cannot be produced in the hadronization.",
     &HadronSpectrum::_forbidden, -1, false, false, true, false, false);

}

