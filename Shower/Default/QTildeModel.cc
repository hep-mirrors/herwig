// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeModel class.
//

#include "QTildeModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "QTildeReconstructor.h"
#include "QTildeFinder.h"
#include "QTildeSudakov.h"
#include "Herwig++/Shower/Default/MECorrections/QTildeMECorrection.h"
#include "ThePEG/Utilities/Throw.h"

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
    ("The QTildeModel class is the ShowerModel object for the Herwig++ shower.",
     "The Shower evolution was perform using the algorithm suggested in "
     "\\cite{Gieseke:2003rz}.",
     "\\bibitem{Gieseke:2003rz} S.~Gieseke, P.~Stephens and B.~Webber,"
     "JHEP {\\bf 0312} (2003) 045.");

}

void QTildeModel::checkConsistency() throw(InitException) {
  // check KinematicsReconstructor
  if(!dynamic_ptr_cast<QTildeReconstructorPtr>(kinematicsReconstructor()))
    Throw<InitException>() << "KinematicsReconstructor must be either "
			 << "QTildeKinematicsReconstructor or a class inheriting from it"
			 << "in QTildeModel::checkConsistency()";
  // check PartnerFinder
  if(!dynamic_ptr_cast<QTildeFinderPtr>(partnerFinder()))
    Throw<InitException>() << "PartnerFinder must be either "
			   << "QTildeFinder or a class inheriting from it"
			   << "in QTildeModel::checkConsistency()";
  // Sudakov form factors
  vector<SudakovPtr>::const_iterator sit;
  for(sit=sudakovFormFactors().begin();sit!=sudakovFormFactors().end();++sit) {
    if(!dynamic_ptr_cast<QTildeSudakovPtr>(*sit))
      Throw<InitException>() << "SudakovFormFactors must be either "
			     << "QTildeSudakov or a class inheriting from it"
			     << "in QTildeModel::checkConsistency()"; 
  }
  // Matrix element corrections
  // check KinematicsReconstructor
  vector<MECorrectionPtr>::const_iterator mit;
  for(mit=meCorrections().begin();mit!=meCorrections().end();++mit) {
  if(!dynamic_ptr_cast<QTildeMECorrectionPtr>(*mit)) {
    Throw<InitException>() << "meCorrections must be either "
			   << "TildeMECorrection or a class inheriting from it"
			   << "in QTildeModel::checkConsistency()"; 
  }
  }
}
