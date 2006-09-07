// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeModel class.
//

#include "QTildeModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/SudakovFormFactor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"

using namespace Herwig;

IBPtr QTildeModel::clone() const {
  return new_ptr(*this);
}

IBPtr QTildeModel::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<QTildeModel> QTildeModel::initQTildeModel;
// Definition of the static class description member.

void QTildeModel::Init() {

  static ClassDocumentation<QTildeModel> documentation
    ("The QTildeModel class is the ShowerModel object for the Herwig++ shower.");

}

void QTildeModel::checkConsistency() throw(InitException) {
}
