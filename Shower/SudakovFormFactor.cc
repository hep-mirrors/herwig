// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


SudakovFormFactor::~SudakovFormFactor() {}


AbstractClassDescription<SudakovFormFactor> SudakovFormFactor::initSudakovFormFactor;
// Definition of the static class description member.


void SudakovFormFactor::Init() {

  static ClassDocumentation<SudakovFormFactor> documentation
    ("This abstract class is the base class for all the other Sudakov form factors classes.");

}


void SudakovFormFactor::reset() {
  _q = Energy();
  _z = 0.0;
  _phi = 0.0;
} 


void SudakovFormFactor::setupLookupTables() {}

