// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaHadronicMass class.
//

#include "BtoSGammaHadronicMass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BtoSGammaHadronicMass.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

BtoSGammaHadronicMass::~BtoSGammaHadronicMass() {}

void BtoSGammaHadronicMass::persistentOutput(PersistentOStream & os) const {
  os << _minMass << _maxMass;
}

void BtoSGammaHadronicMass::persistentInput(PersistentIStream & is, int) {
  is >> _minMass >> _maxMass;
}

AbstractClassDescription<BtoSGammaHadronicMass> BtoSGammaHadronicMass::initBtoSGammaHadronicMass;
// Definition of the static class description member.

void BtoSGammaHadronicMass::Init() {

  static ClassDocumentation<BtoSGammaHadronicMass> documentation
    ("The BtoSGammaHadronicMass class is the base class for the implementation"
     " of models of the hadronic spectrum in B to s gamma decays.");

  static Parameter<BtoSGammaHadronicMass,Energy> interfaceMinimumMass
    ("MinimumMass",
     "The minimum value of the hadronic mass",
     &BtoSGammaHadronicMass::_minMass, MeV, 825*MeV, 825*MeV, 5300*MeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaHadronicMass,Energy> interfaceMaximumMass
    ("MaximumMass",
     "The maximum value of the hadronic mass",
     &BtoSGammaHadronicMass::_maxMass, MeV, 5300*MeV, 825*MeV, 5300*MeV,
     false, false, Interface::limited);


}

void BtoSGammaHadronicMass::dataBaseOutput(ofstream & output,bool header,
					   bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::BtoSGammaHadronicMass " << fullName() << " \n";}
  output << "set " << fullName() << ":MinimumMass " << _minMass/MeV << " \n";
  output << "set " << fullName() << ":MaximumMass " << _maxMass/MeV << " \n";
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
