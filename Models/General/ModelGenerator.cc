// -*- C++ -*-
//
// ModelGenerator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ModelGenerator class.
//

#include "ModelGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "BSMWidthGenerator.h"
#include "Herwig++/PDT/GenericMassGenerator.h"
#include "ThePEG/Repository/BaseRepository.h"

using namespace Herwig;

IBPtr ModelGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr ModelGenerator::fullclone() const {
  return new_ptr(*this);
}

void ModelGenerator::persistentOutput(PersistentOStream & os) const {
  os << _theHPConstructor << _theDecayConstructor << _theParticles 
     << _theRPConstructor << _theOffshell << _theOffsel << _theBRnorm
     << _theNpoints << _theIorder << _theBWshape;
}

void ModelGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _theHPConstructor >> _theDecayConstructor >> _theParticles
     >> _theRPConstructor >> _theOffshell >> _theOffsel >> _theBRnorm
     >> _theNpoints >> _theIorder >> _theBWshape;
}

bool ModelGenerator::preInitialize() const {
  return true;
}

ClassDescription<ModelGenerator> ModelGenerator::initModelGenerator;
// Definition of the static class description member.

void ModelGenerator::Init() {

  static ClassDocumentation<ModelGenerator> documentation
    ("This class controls the the use of BSM physics.",
     "BSM physics was produced using the algorithm of "
     "\\cite{Gigg2007:cr}",
     "\\bibitem{Gigg:2007cr} M.~Gigg and P.~Richardson, \n"
     "Eur.\\ Phys.\\ J.\\  C {\\bf 51} (2007) 989.\n"
     "%%CITATION = EPHJA,C51,989;%%");
  
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
    ("DecayParticles",
     "ParticleData pointers to the particles requiring spin correlation "
     "decayers. If decay modes do not exist they will also be created.",
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

  static Switch<ModelGenerator,int> interfaceWhichOffshell
    ("WhichOffshell",
     "A switch to determine which particles to create mass and width "
     "generators for.",
     &ModelGenerator::_theOffsel, 0, false, false);
  static SwitchOption interfaceWhichOffshellSelected
    (interfaceWhichOffshell,
     "Selected",
     "Only create mass and width generators for the particles specified",
     0);
  static SwitchOption interfaceWhichOffshellAll
    (interfaceWhichOffshell,
     "All",
     "Treat all particles specified in the DecayParticles "
     "list as off-shell",
     1);
  
  static Switch<ModelGenerator,bool> interfaceBRNormalize
    ("BRNormalize",
     "Whether to normalize the partial widths to BR*total width for an "
     "on-shell particle",
     &ModelGenerator::_theBRnorm, true, false, false);
  static SwitchOption interfaceBRNormalizeNormalize
    (interfaceBRNormalize,
     "Yes",
     "Normalize the partial widths",
     true);
  static SwitchOption interfaceBRNormalizeNoNormalize
    (interfaceBRNormalize,
     "No",
     "Do not normalize the partial widths",
     false);

  static Parameter<ModelGenerator,int> interfacePoints
    ("InterpolationPoints",
     "Number of points to use for interpolation tables when needed",
     &ModelGenerator::_theNpoints, 50, 5, 1000,
     false, false, true);
  
  static Parameter<ModelGenerator,unsigned int> 
    interfaceInterpolationOrder
    ("InterpolationOrder", "The interpolation order for the tables",
     &ModelGenerator::_theIorder, 1, 1, 5,
     false, false, Interface::limited);

  static Switch<ModelGenerator,int> interfaceBreitWignerShape
    ("BreitWignerShape",
     "Controls the shape of the mass distribution generated",
     &ModelGenerator::_theBWshape, 0, false, false);
  static SwitchOption interfaceBreitWignerShapeDefault
    (interfaceBreitWignerShape,
     "Default",
     "Running width with q in numerator and denominator width factor",
     0);
  static SwitchOption interfaceBreitWignerShapeFixedWidth
    (interfaceBreitWignerShape,
     "FixedWidth",
     "Use a fixed width",
     1);
  static SwitchOption interfaceBreitWignerShapeNoq
    (interfaceBreitWignerShape,
     "Noq",
     "Use M rather than q in the numerator and denominator width factor",
     2);
  static SwitchOption interfaceBreitWignerShapeNoNumerator
    (interfaceBreitWignerShape,
     "NoNumerator",
     "Neglect the numerator factors",
     3);
}

void ModelGenerator::doinit() throw(InitException) {
  Interfaced::doinit();
  useMe();
  //create mass and width generators for the requested particles
  PDVector::iterator pit, pend;
  if( _theOffsel == 0 ) {
    pit = _theOffshell.begin();
    pend = _theOffshell.end();
  }
  else {
    pit = _theParticles.begin();
    pend = _theParticles.end();
  }
  for(; pit != pend; ++pit)
    createWidthGenerator(*pit);

  //create decayers and decaymodes (if necessary)
  if( _theDecayConstructor )
    _theDecayConstructor->createDecayers(_theParticles);

  // write out decays with spin correlations and set particles
  // that have no decay modes to stable.
  string filename = CurrentGenerator::current().filename() + 
    string("-BSMModelInfo.out");
  ofstream ofs(filename.c_str(), ios::out|ios::app);
  ofs << "# The decay modes listed below will have spin\n"
      << "# correlations included when they are generated.\n#\n#";
  pit = _theParticles.begin();
  pend = _theParticles.end();
  for( ; pit != pend; ++pit) {
    tPDPtr parent = *pit;
    // Check decays for ones where quarks cannot be put on constituent
    // mass-shell
    checkDecays(parent);
    parent->touch();
    parent->update();
    if( parent->CC() ) parent->CC()->synchronize();

    if( parent->decaySelector().empty() ) {
      parent->stable(true);
      parent->width(0.0*MeV);
      parent->massGenerator(tGenericMassGeneratorPtr());
      parent->widthGenerator(tGenericWidthGeneratorPtr());
    }
    else
      writeDecayModes(ofs, parent);

    if( parent->massGenerator() ) {
      parent->massGenerator()->reset();
      ofs << "# " <<parent->PDGName() << " will be considered off-shell.\n#";
    }
    if( parent->widthGenerator() ) parent->widthGenerator()->reset();
  }
  //Now construct hard processes given that we know which
  //objects have running widths
  if(_theHPConstructor) {
    _theHPConstructor->init();
    _theHPConstructor->constructDiagrams();
  }
  if(_theRPConstructor) {
    _theRPConstructor->init();
    _theRPConstructor->constructResonances();
  }

}

void ModelGenerator::checkDecays(PDPtr parent) {
  DecaySet::iterator dit = parent->decayModes().begin();
  DecaySet::iterator dend = parent->decayModes().end();
  Energy oldwidth(parent->width()), newwidth(parent->width());
  bool rescalebrat(false);
  for(; dit != dend; ++dit ) {
    if( !(**dit).on() ) continue;
    Energy release((**dit).parent()->mass());
    PDVector::const_iterator pit = (**dit).orderedProducts().begin();
    PDVector::const_iterator pend =(**dit).orderedProducts().end();
    for( ; pit != pend; ++pit ) {
      release -= (**pit).constituentMass();
    }
    if( release < 0.*MeV ) {
      rescalebrat = true;
      newwidth -= (**dit).brat()*oldwidth;
      generator()->preinitInterface(*dit, "OnOff", "set", "Off");
      generator()->preinitInterface(*dit, "BranchingRatio", 
				    "set", "0.0");
    }
  }

  if( rescalebrat ) {
    dit = parent->decayModes().begin();
    dend = parent->decayModes().end();
    for( ; dit != dend; ++dit ) {
      if( !(**dit).on() ) continue;
      double newbrat = ((**dit).brat())*oldwidth/newwidth;
      ostringstream brf;
      brf << newbrat;
      generator()->preinitInterface(*dit, "BranchingRatio",
				    "set", brf.str());
    }
    parent->width(newwidth);
    parent->widthCut(5*newwidth);
  }
  
}

void ModelGenerator::writeDecayModes(ofstream & ofs, tcPDPtr parent) const {
  ofs << " Parent: " << parent->PDGName() << "  Mass (GeV): " 
      << parent->mass()/GeV << "  Total Width (GeV): " 
      << parent->width()/GeV << endl;
  ofs << std::left << std::setw(40) << '#' 
      << std::left << std::setw(20) << "Partial Width/GeV"
      << "BR\n"; 
  Selector<tDMPtr>::const_iterator dit = parent->decaySelector().begin();
  Selector<tDMPtr>::const_iterator dend = parent->decaySelector().end();
  for(; dit != dend; ++dit)
    ofs << std::left << std::setw(40) << (*dit).second->tag() 
	<< std::left << std::setw(20) << (*dit).second->brat()*parent->width()/GeV 
	<< (*dit).second->brat() << '\n';
  ofs << "#\n#";
}

void ModelGenerator::createWidthGenerator(tPDPtr p) {
  string wn = p->fullName() + string("-WGen");
  string mn = p->fullName() + string("-MGen");
  GenericMassGeneratorPtr mgen = dynamic_ptr_cast<GenericMassGeneratorPtr>
    (generator()->preinitCreate("Herwig::GenericMassGenerator", mn));
  BSMWidthGeneratorPtr wgen = dynamic_ptr_cast<BSMWidthGeneratorPtr>
    (generator()->preinitCreate("Herwig::BSMWidthGenerator", wn));

  //set the particle interface
  generator()->preinitInterface(mgen, "Particle", "set", p->fullName());
  generator()->preinitInterface(wgen, "Particle", "set", p->fullName());

  //set the generator interfaces in the ParticleData object
  generator()->preinitInterface(p, "Mass_generator","set", mn);
  generator()->preinitInterface(p, "Width_generator","set", wn);
  //allow the branching fraction of this particle type to vary
  p->variableRatio(true);
  if( p->CC() ) p->CC()->variableRatio(true);
  
  //initialize the generators
  generator()->preinitInterface(mgen, "Initialize", "set", "Yes");
  generator()->preinitInterface(wgen, "Initialize", "set", "Yes");

  string norm = _theBRnorm ? "Yes" : "No";
  generator()->preinitInterface(wgen, "BRNormalize", "set", norm);
  ostringstream os;
  os << _theNpoints;
  generator()->preinitInterface(wgen, "InterpolationPoints", "set", 
				  os.str());
  os.str("");
  os << _theIorder;
  generator()->preinitInterface(wgen, "InterpolationOrder", "set",
				  os.str());
  os.str("");
  os << _theBWshape;
  generator()->preinitInterface(mgen, "BreitWignerShape", "set", 
				  os.str());
}
