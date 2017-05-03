// -*- C++ -*-
//
// QTildeModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeModel class.
//

#include "QTildeModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "QTildeReconstructor.h"
#include "QTildeFinder.h"
#include "QTildeSudakov.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<QTildeModel,Herwig::ShowerModel>
describeQTildeModel ("Herwig::QTildeModel","HwShower.so");

IBPtr QTildeModel::clone() const {
  return new_ptr(*this);
}

IBPtr QTildeModel::fullclone() const {
  return new_ptr(*this);
}

void QTildeModel::Init() {

  static ClassDocumentation<QTildeModel> documentation
    ("The QTildeModel class is the ShowerModel object for the Herwig shower.");

}

void QTildeModel::checkConsistency() {
  // check KinematicsReconstructor
  if(!dynamic_ptr_cast<Ptr<QTildeReconstructor>::pointer>(kinematicsReconstructor()))
    Throw<InitException>() << "KinematicsReconstructor must be either "
			 << "QTildeKinematicsReconstructor or a class inheriting from it"
			 << "in QTildeModel::checkConsistency()";
  // check PartnerFinder
  if(!dynamic_ptr_cast<Ptr<QTildeFinder>::pointer>(partnerFinder()))
    Throw<InitException>() << "PartnerFinder must be either "
			   << "QTildeFinder or a class inheriting from it"
			   << "in QTildeModel::checkConsistency()";
  // Sudakov form factors
  vector<SudakovPtr>::const_iterator sit;
  for(sit=sudakovFormFactors().begin();sit!=sudakovFormFactors().end();++sit) {
    if(!dynamic_ptr_cast<Ptr<QTildeSudakov>::pointer>(*sit))
      Throw<InitException>() << "SudakovFormFactors must be either "
			     << "QTildeSudakov or a class inheriting from it"
			     << "in QTildeModel::checkConsistency()"; 
  }
}
