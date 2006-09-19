// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranModel class.
//

#include "FortranModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "FortranFinder.h"
#include "FortranReconstructor.h"
#include "FortranSudakov.h"
#include "ThePEG/Utilities/Throw.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"


using namespace Herwig;

NoPIOClassDescription<FortranModel> FortranModel::initFortranModel;
// Definition of the static class description member.

void FortranModel::Init() {

  static ClassDocumentation<FortranModel> documentation
    ("There is no documentation for the FortranModel class");

}

void FortranModel::checkConsistency() throw(InitException) {
  // check KinematicsReconstructor
  if(!dynamic_ptr_cast<FortranReconstructorPtr>(kinematicsReconstructor()))
    Throw<InitException>() << "KinematicsReconstructor must be either "
			 << "FortranReconstructor or a class inheriting from it"
			 << "in FortranModel::checkConsistency()";
  // check PartnerFinder
  if(!dynamic_ptr_cast<FortranFinderPtr>(partnerFinder()))
    Throw<InitException>() << "PartnerFinder must be either "
			   << "FortranFinder or a class inheriting from it"
			   << "in FortranModel::checkConsistency()";
  // Sudakov form factors
  vector<SudakovPtr>::const_iterator sit;
  for(sit=sudakovFormFactors().begin();sit!=sudakovFormFactors().end();++sit) {
    if(!dynamic_ptr_cast<FortranSudakovPtr>(*sit))
      Throw<InitException>() << "SudakovFormFactors must be either "
			     << "FortranSudakov or a class inheriting from it"
			     << "in FortranModel::checkConsistency()"; 
  }
  // no me corrections yet
  if(!meCorrections().empty())
    Throw<InitException>() << "Fortran ME corrections not yet implemented"
			   << Exception::runerror;


//   // Matrix element corrections
//   // check KinematicsReconstructor
//   vector<MECorrectionPtr>::const_iterator mit;
//   for(mit=meCorrections().begin();mit!=meCorrections().end();++mit) {
//   if(!dynamic_ptr_cast<FortranMECorrectionPtr>(*mit)) {
//     Throw<InitException>() << "meCorrections must be either "
// 			   << "TildeMECorrection or a class inheriting from it"
// 			   << "in FortranModel::checkConsistency()"; 
//   }
//  }
}
