// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MassCutOff class.
//

#include "MassCutOff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr MassCutOff::clone() const {
  return new_ptr(*this);
}

IBPtr MassCutOff::fullclone() const {
  return new_ptr(*this);
}

void MassCutOff::persistentOutput(PersistentOStream & os) const {
	os << ounit(vgcut_,GeV) << ounit(vqcut_,GeV);
}

void MassCutOff::persistentInput(PersistentIStream & is, int) {
	is >> iunit(vgcut_,GeV) >> iunit(vqcut_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MassCutOff,SudakovCutOff>
describeHerwigMassCutOff("Herwig::MassCutOff", "HwShower.so");

void MassCutOff::Init() {

  static ClassDocumentation<MassCutOff> documentation
    ("There is no documentation for the MassCutOff class");


  static Parameter<MassCutOff,Energy> interfaceGluonVirtualityCut
    ("GluonVirtualityCut",
     "For the FORTRAN cut-off option the minimum virtuality of the gluon",
     &MassCutOff::vgcut_, GeV, 0.85*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<MassCutOff,Energy> interfaceQuarkVirtualityCut
    ("QuarkVirtualityCut",
     "For the FORTRAN cut-off option the minimum virtuality added to"
     " the mass for particles other than the gluon",
     &MassCutOff::vqcut_, GeV, 0.85*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);
  

}

const vector<Energy> & MassCutOff::virtualMasses(const IdList & ids) {
  static vector<Energy> output;
  output.clear();
  for(auto id : ids) {
      output.push_back(id->mass());
      output.back() += id->id()==ParticleID::g ? vgcut_ : vqcut_;
  }
  return output;
}

