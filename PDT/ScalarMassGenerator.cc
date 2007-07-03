// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMassGenerator class.
//

#include "ScalarMassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarMassGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void ScalarMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,GeV) << ounit(_mass1,GeV) << ounit(_mass2,GeV)
     << ounit(_mplus,GeV) << ounit(_mminus,GeV) 
     << ounit(_m2plus,GeV2) << ounit(_m2minus,GeV2)
     << ounit(_term,GeV2);
}

void ScalarMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,GeV) >> iunit(_mass1,GeV) >> iunit(_mass2,GeV)
     >> iunit(_mplus,GeV) >> iunit(_mminus,GeV) 
     >> iunit(_m2plus,GeV2) >> iunit(_m2minus,GeV2)
     >> iunit(_term,GeV2);
}

ClassDescription<ScalarMassGenerator> ScalarMassGenerator::initScalarMassGenerator;
// Definition of the static class description member.

void ScalarMassGenerator::Init() {

  static ClassDocumentation<ScalarMassGenerator> documentation
    ("The ScalarMassGenerator class inherits from the "
     "GenericMassGenerator class and includes finite width effects "
     "for the scalar f_0 and a_0 mesons.");

  static ParVector<ScalarMassGenerator,Energy> interfacecoupling
    ("Coupling",
     "The coupling",
     &ScalarMassGenerator::_coupling,
     GeV, 0, 0*GeV, 0*GeV, Constants::MaxEnergy, false, false, true);

  static ParVector<ScalarMassGenerator,Energy> interfacemass1
    ("Mass1",
     "The mass for first particle",
     &ScalarMassGenerator::_mass1,
     GeV, 0, 0*GeV, 0*GeV, Constants::MaxEnergy, false, false, true);

  static ParVector<ScalarMassGenerator,Energy> interfacemass2
    ("Mass2",
     "The mass for second particle",
     &ScalarMassGenerator::_mass2,
     GeV, 0, 0*GeV, 0*GeV, Constants::MaxEnergy, false, false, true);

}


void ScalarMassGenerator::dataBaseOutput(ofstream & output)
{
  output << "update Mass_Generators set parameters=\"";
  output << "set " << fullName() << ":BreitWignerShape "   << _BWshape << "\n";
  output << "set " << fullName() << ":MaximumWeight " << _maxwgt    << "\n";
  output << "set " << fullName() << ":NGenerate "   << _ngenerate << "\n";
  for(unsigned int ix=0;ix<_coupling.size();++ix)
    {output << "insert " << fullName() 
	    << ":Couplings " 
	    << ix << " " << _coupling[ix]/GeV << "\n";}

  for(unsigned int ix=0;ix<_mass1.size();++ix)
    {output << "insert " << fullName() 
	    << ":Mass1 " 
	    << ix << " " << _mass1[ix]/GeV << "\n";}
  for(unsigned int ix=0;ix<_mass2.size();++ix)
    {output << "insert " << fullName() 
	    << ":Mass2 " 
	    << ix << " " << _mass2[ix]/GeV << "\n";}
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

}
