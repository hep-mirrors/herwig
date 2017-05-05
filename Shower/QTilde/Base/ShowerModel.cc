// -*- C++ -*-
//
// ShowerModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerModel class.
//

#include "ShowerModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "KinematicsReconstructor.h"
#include "PartnerFinder.h"
#include "Herwig/Shower/Core/Base/SudakovFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeAbstractClass<ShowerModel,Interfaced>
describeShowerModel ("Herwig::ShowerModel","HwShower.so");

void ShowerModel::persistentOutput(PersistentOStream & os) const {
  os << _reconstructor << _partnerfinder << _sudakovs;
}

void ShowerModel::persistentInput(PersistentIStream & is, int) {
  is >> _reconstructor >> _partnerfinder >> _sudakovs;
}

void ShowerModel::doinit() {
  Interfaced::doinit();
  checkConsistency();
}

void ShowerModel::Init() {

  static ClassDocumentation<ShowerModel> documentation
    ("The ShowerModel class contains the references for the classes which"
     " are specific to the shower evolution scheme.");

  static Reference<ShowerModel,KinematicsReconstructor> interfaceKinematicsReconstructor
    ("KinematicsReconstructor",
     "Reference to the KinematicsReconstructor object",
     &ShowerModel::_reconstructor, false, false, true, false, false);

  static Reference<ShowerModel,PartnerFinder> interfacePartnerFinder
    ("PartnerFinder",
     "Reference to the PartnerFinder object",
     &ShowerModel::_partnerfinder, false, false, true, false, false);

  static RefVector<ShowerModel,SudakovFormFactor> interfaceSudakovFormFactors
    ("SudakovFormFactors",
     "Vector of references to the SudakovFormFactor objects",
     &ShowerModel::_sudakovs, -1, false, false, true, false, false);

}

