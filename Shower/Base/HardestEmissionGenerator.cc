// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardestEmissionGenerator class.
//

#include "HardestEmissionGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/HardTree.h"
#include "Herwig++/Shower/Base/Evolver.h"

using namespace Herwig;

AbstractClassDescription<HardestEmissionGenerator> HardestEmissionGenerator::initHardestEmissionGenerator;
// Definition of the static class description member.

void HardestEmissionGenerator::Init() {

  static ClassDocumentation<HardestEmissionGenerator> documentation
    ("The HardestEmissionGenerator class is the base class for the generation"
     "of the hardest emission in the POWHEG shower approach");

}

void HardestEmissionGenerator::setEvolver(tEvolverPtr in) {
  _evolver=in;
}

void HardestEmissionGenerator::persistentOutput(PersistentOStream & os) const {
  os << _evolver;
}

void HardestEmissionGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _evolver;
}

void HardestEmissionGenerator::rebind(const TranslationMap & trans) {
  _evolver = trans.translate(_evolver);
  Interfaced::rebind(trans);
}

IVector HardestEmissionGenerator::getReferences() {
  IVector ret = Interfaced::getReferences();
  ret.push_back(_evolver);
  return ret;
}
