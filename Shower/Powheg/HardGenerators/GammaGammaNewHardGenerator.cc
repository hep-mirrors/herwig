// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGammaNewHardGenerator class.
//

#include "GammaGammaNewHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Powheg/MEPP2GammaGammaPowhegNEW.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/ShowerTree.h"

using namespace Herwig;

GammaGammaNewHardGenerator::GammaGammaNewHardGenerator() {}

IBPtr GammaGammaNewHardGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr GammaGammaNewHardGenerator::fullclone() const {
  return new_ptr(*this);
}

void GammaGammaNewHardGenerator::persistentOutput(PersistentOStream & os) const {
}

void GammaGammaNewHardGenerator::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<GammaGammaNewHardGenerator> 
GammaGammaNewHardGenerator::initGammaGammaNewHardGenerator;
// Definition of the static class description member.

void GammaGammaNewHardGenerator::Init() {

  static ClassDocumentation<GammaGammaNewHardGenerator> documentation
    ("There is no documentation for the GammaGammaNewHardGenerator class");

}

HardTreePtr GammaGammaNewHardGenerator::generateHardest(ShowerTreePtr tree) {
  ThePEG::Ptr<MEPP2GammaGammaPowhegNEW>::transient_const_pointer meptr = 
    dynamic_ptr_cast<ThePEG::Ptr<MEPP2GammaGammaPowhegNEW>::
    transient_const_pointer>(ShowerHandler::currentHandler()->
			     currentSubProcess()->handler());
  assert(meptr);
  return meptr->generateHardest(tree);
}

bool GammaGammaNewHardGenerator::canHandle(ShowerTreePtr) {
  ThePEG::Ptr<MEPP2GammaGammaPowhegNEW>::transient_const_pointer meptr = 
    dynamic_ptr_cast<ThePEG::Ptr<MEPP2GammaGammaPowhegNEW>::
    transient_const_pointer>(ShowerHandler::currentHandler()->
			     currentSubProcess()->handler());
  return meptr;
}
