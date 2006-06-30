// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayerBase class.
//

#include "HwDecayerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HwDecayerBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

bool HwDecayerBase::accept(const DecayMode & dm) const {
  return false;
}

ParticleVector HwDecayerBase::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  return children;
}


void HwDecayerBase::persistentOutput(PersistentOStream & os) const {
  os << _initialize;
}

void HwDecayerBase::persistentInput(PersistentIStream & is, int) {
  is >> _initialize;
}

AbstractClassDescription<HwDecayerBase> HwDecayerBase::initHwDecayerBase;
// Definition of the static class description member.

void HwDecayerBase::Init() {

  static ClassDocumentation<HwDecayerBase> documentation
    ("The HwDecayerBase class is the base class for Decayers in Hw++.");

  static Switch<HwDecayerBase,bool> interfaceInitialize
    ("Initialize",
     "Initialization of the phase space calculation",
     &HwDecayerBase::_initialize, false, false, false);
  static SwitchOption interfaceInitializeon
    (interfaceInitialize,
     "on",
     "At initialisation find max weight and optimise the integration",
     true);
  static SwitchOption interfaceInitializeoff
    (interfaceInitialize,
     "off",
     "Use the maximum weight and channel weights supplied for the integration",
     false);

}

void HwDecayerBase::dataBaseOutput(ofstream & output,bool header) const
{
  // header for MySQL
  if(header){output << "update decayers set parameters=\"";}
  // photon generator if it exists
  // footer for MySQL
  if(header){output << " where ThePEGName=\" " << fullName() << "\";";}
}
