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

ScalarMassGenerator::~ScalarMassGenerator() {}

void ScalarMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _mass1 << _mass2 << _mplus << _mminus << _m2plus << _m2minus
     << _term;
}

void ScalarMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _mass1 >> _mass2 >> _mplus >> _mminus >> _m2plus >> _m2minus
     >> _term;
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
     0, 0, 0, 0, 1.E12, false, false, true);

  static ParVector<ScalarMassGenerator,Energy> interfacemass1
    ("Mass1",
     "The mass for first particle",
     &ScalarMassGenerator::_mass1,
     0, 0, 0, 0,  1.E12, false, false, true);

  static ParVector<ScalarMassGenerator,Energy> interfacemass2
    ("Mass2",
     "The mass for second particle",
     &ScalarMassGenerator::_mass2,
     0, 0, 0, 0,  1.E12, false, false, true);

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
	    << ix << " " << _coupling[ix] << "\n";}

  for(unsigned int ix=0;ix<_mass1.size();++ix)
    {output << "insert " << fullName() 
	    << ":Mass1 " 
	    << ix << " " << _mass1[ix] << "\n";}
  for(unsigned int ix=0;ix<_mass2.size();++ix)
    {output << "insert " << fullName() 
	    << ":Mass2 " 
	    << ix << " " << _mass2[ix] << "\n";}
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

}
