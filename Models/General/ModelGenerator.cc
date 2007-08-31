// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ModelGenerator class.
//

#include "ModelGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "BSMWidthGenerator.h"
#include "Herwig++/PDT/GenericMassGenerator.h"

using namespace Herwig;

void ModelGenerator::persistentOutput(PersistentOStream & os) const {
  os << _theHPConstructor << _theDecayConstructor << _theParticles 
     << _theRPConstructor << _theOffshell;
}

void ModelGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _theHPConstructor >> _theDecayConstructor >> _theParticles
     >> _theRPConstructor >> _theOffshell;
}

bool ModelGenerator::preInitialize() const {
  return true;
}

ClassDescription<ModelGenerator> ModelGenerator::initModelGenerator;
// Definition of the static class description member.

void ModelGenerator::Init() {

  static ClassDocumentation<ModelGenerator> documentation
    ("There is no documentation for the ModelGenerator class");
  
  static Reference<ModelGenerator,Herwig::HardProcessConstructor> 
    interfaceHardProcessConstructor
    ("HardProcessConstructor",
     "Pointer to the object that constructs the hard process",
     &ModelGenerator::_theHPConstructor, false, false, true, true);

  static Reference<ModelGenerator,Herwig::DecayConstructor> 
     interfaceDecayConstructor
     ("DecayConstructor",
      "Pointer to DecayConstructor helper class",
      &ModelGenerator::_theDecayConstructor, false, false, true, false);
  
  static RefVector<ModelGenerator,ThePEG::ParticleData> interfaceModelParticles
    ("ModelParticles",
     "Pointers to particles contained in model",
     &ModelGenerator::_theParticles, -1, false, false, true, false);

  static Reference<ModelGenerator,Herwig::ResonantProcessConstructor> 
    interfaceResonantProcessConstructor
    ("ResonantProcessConstructor",
     "Pointer to the object that constructs the resonant process(es)",
     &ModelGenerator::_theRPConstructor, false, false, true, true);
  
  
  static RefVector<ModelGenerator,ParticleData> interfaceOffshell
    ("Offshell",
     "The particles to treat as off-shell",
     &ModelGenerator::_theOffshell, -1, false, false, true, false);

}

void ModelGenerator::doinit() throw(InitException) {
  Interfaced::doinit();
  if(_theHPConstructor) {
    _theHPConstructor->init();
    _theHPConstructor->constructDiagrams();
  }
  if(_theRPConstructor) {
    _theRPConstructor->init();
    _theRPConstructor->constructResonances();
  }
  if( _theParticles.empty() ) return;
  //create mass and width generators for the requested particles
  PDVector::iterator pit = _theOffshell.begin();
  PDVector::iterator pend = _theOffshell.end();
  for(; pit != pend; ++pit) {
    PDPtr inpart = *pit;
    GenericMassGeneratorPtr mgen = dynamic_ptr_cast<GenericMassGeneratorPtr>
      (generator()->preinitCreate("Herwig::GenericMassGenerator",
				  inpart->fullName() + string("-MGen")));
    generator()->preinitInterface(mgen, "Particle", "set",
				  inpart->fullName());
    generator()->preinitInterface(inpart, "Mass_generator","set",
				  inpart->fullName() + string("-MGen"));
    generator()->preinitInterface(mgen, "Initialize", "set",
				  "Initialization");

    BSMWidthGeneratorPtr wgen = dynamic_ptr_cast<BSMWidthGeneratorPtr>
      (generator()->preinitCreate("Herwig::BSMWidthGenerator",
				  inpart->fullName() + string("-WGen")));
    generator()->preinitInterface(wgen, "Particle", "set",
				  inpart->fullName());
    generator()->preinitInterface(inpart, "Width_generator","set",
				  inpart->fullName() + string("-WGen"));
    generator()->preinitInterface(wgen, "Initialize", "set",
				  "Initialization");
  }
  //create decayers and decaymodes (if necessary)
  _theDecayConstructor->createDecayers(_theParticles);
  string filename = CurrentGenerator::current().filename() + 
    string("-BSMDecayModes.out");
  ofstream ofs(filename.c_str());
  ofs << "# BSM Model Decay Modes\n";
  ofs << "#\n#";
  pit = _theParticles.begin();
  pend = _theParticles.end();
  for( ; pit != pend; ++pit) {
    tPDPtr parent = *pit;
    if( parent->decaySelector().empty() ) {
      parent->stable(true);
      parent->width(0.0*MeV);
      parent->massGenerator(tGenericMassGeneratorPtr());
      parent->widthGenerator(tGenericWidthGeneratorPtr());
    }
    else
      writeDecayModes(ofs, parent);

    if(parent->massGenerator()) parent->massGenerator()->reset();
    if(parent->widthGenerator()) parent->widthGenerator()->reset();

//     cerr << "Testing mass and width gens " << parent->PDGName() << "\n"
// 	 << "Mass gen " << parent->massGenerator() << "  ";
//     if(parent->massGenerator()) cerr << parent->massGenerator()->state() << '\n';
//     cerr<< "Width gen " << parent->widthGenerator() << " ";
//     if( parent->widthGenerator() )
//       cerr << parent->widthGenerator()->state() << endl;
      
  }
  
}

void ModelGenerator::writeDecayModes(ofstream & ofs, tcPDPtr parent) const {
  ofs << " Parent: " << parent->PDGName() << "  Mass (GeV): " 
      << parent->mass()/GeV << "  Total Width (GeV): " 
      << parent->width()/GeV << endl;
  ofs << std::left << std::setw(48) << '#' << "BR" << endl; 
  Selector<tDMPtr>::const_iterator dit = parent->decaySelector().begin();
  Selector<tDMPtr>::const_iterator dend = parent->decaySelector().end();
  for(; dit != dend; ++dit)
    ofs << std::setw(48) << (*dit).second->tag() << (*dit).second->brat() 
	<< '\n';
  ofs << "#\n#";
  
}


