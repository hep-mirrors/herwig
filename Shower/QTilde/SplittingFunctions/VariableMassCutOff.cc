// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VariableMassCutOff class.
//

#include "VariableMassCutOff.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr VariableMassCutOff::clone() const {
  return new_ptr(*this);
}

IBPtr VariableMassCutOff::fullclone() const {
  return new_ptr(*this);
}

void VariableMassCutOff::persistentOutput(PersistentOStream & os) const {
	os << a_ << b_ << ounit(c_,GeV) << ounit(kinCutoffScale_, GeV);
}

void VariableMassCutOff::persistentInput(PersistentIStream & is, int) {
	is >> a_ >> b_ >> iunit(c_,GeV) >> iunit(kinCutoffScale_, GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VariableMassCutOff,SudakovCutOff>
describeHerwigVariableMassCutOff("Herwig::VariableMassCutOff", "HwShower.so");

void VariableMassCutOff::Init() {

  static ClassDocumentation<VariableMassCutOff> documentation
    ("There is no documentation for the VariableMassCutOff class");


  static Parameter<VariableMassCutOff,double> interfaceaParameter
    ("aParameter",
     "The a parameter for the kinematic cut-off",
     &VariableMassCutOff::a_, 0.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<VariableMassCutOff,double> interfacebParameter
    ("bParameter",
     "The b parameter for the kinematic cut-off",
     &VariableMassCutOff::b_, 2.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<VariableMassCutOff,Energy> interfacecParameter
    ("cParameter",
     "The c parameter for the kinematic cut-off",
     &VariableMassCutOff::c_, GeV, 0.3*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<VariableMassCutOff,Energy>
    interfaceKinScale ("cutoffKinScale",
		       "kinematic cutoff scale for the parton shower phase"
		       " space (unit [GeV])",
		       &VariableMassCutOff::kinCutoffScale_, GeV, 
		       2.3*GeV, 0.001*GeV, 10.0*GeV,false,false,false);
  

}


const vector<Energy> & VariableMassCutOff::virtualMasses(const IdList & ids) {
  static vector<Energy> output;
  output.clear();

  for(auto id : ids)
      output.push_back(id->mass());

  Energy kinCutoff = kinematicCutOff(
                kinCutoffScale_,
                *std::max_element(output.begin(),output.end())
        );

  for(auto & el : output)
      el = max(kinCutoff, el);
 
  return output;
}
