// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanMECorrection class.
//

#include "DrellYanMECorrection.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


DrellYanMECorrection::~DrellYanMECorrection() {}


ClassDescription<DrellYanMECorrection> DrellYanMECorrection::initDrellYanMECorrection;
// Definition of the static class description member.


void DrellYanMECorrection::Init() {

  static ClassDocumentation<DrellYanMECorrection> documentation
    ("This class is responsible for Drell-Yan matrix element correction.");

}


//---------------------------------------------------------------------------------

bool DrellYanMECorrection::
softMECorrection( Lorentz5Momentum pEmissionHardestSoFar ) throw(Veto, Stop, Exception) {
  bool reject = false;
  
  //***LOOKHERE*** WRITE THE CODE
  
  return reject;
}





