// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingFunction class.
//

#include "SplittingFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace Herwig;

AbstractNoPIOClassDescription<SplittingFunction> SplittingFunction::initSplittingFunction;
// Definition of the static class description member.

void SplittingFunction::Init() {

  static ClassDocumentation<SplittingFunction> documentation
    ("The SplittingFunction class is the based class for 1->2 splitting functions"
     " in Herwig++");

}

double SplittingFunction::generatePhi(ShowerParticle &,ShoKinPtr,
				      const double, const Energy,
				      const IdList &, const RhoDMatrix &,
				      const double) {
  cerr << "Using SplittingFunction::generatePhi()" << fullName() << "\n";
  exit(1);
  return Constants::twopi*UseRandom::rnd();
}

DecayMatrixElement SplittingFunction::matrixElement(ShowerParticle &,ShoKinPtr,
						    const double, const Energy, 
						    const IdList &,
						    const RhoDMatrix &, const double) {
  throw Exception() << "SplittingFunction::matrixElement called for " << fullName()
		    << Exception::runerror;
}
