// -*- C++ -*-
//
// GenericMassGenerator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericMassGenerator class.
//

#include "GenericMassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/Rebinder.h"

using namespace Herwig;
using namespace ThePEG;

void GenericMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << _particle << _antiparticle 
     << ounit(_lowermass,GeV) << ounit(_uppermass,GeV) << _maxwgt 
     << _BWshape << _ngenerate 
     << ounit(_mass,GeV) << ounit(_width,GeV) 
     << ounit(_mass2,GeV2) << ounit(_mwidth,GeV2) 
     << _ninitial << _initialize << _widthgen << _widthopt;
}

void GenericMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _particle >> _antiparticle 
     >> iunit(_lowermass,GeV) >> iunit(_uppermass,GeV) >> _maxwgt 
     >> _BWshape >> _ngenerate 
     >> iunit(_mass,GeV) >> iunit(_width ,GeV)
     >> iunit(_mass2,GeV2) >> iunit(_mwidth ,GeV2)
     >> _ninitial >> _initialize >> _widthgen >> _widthopt;
}

ClassDescription<GenericMassGenerator> GenericMassGenerator::initGenericMassGenerator;
// Definition of the static class description member.


void GenericMassGenerator::Init() {

  static ClassDocumentation<GenericMassGenerator> documentation
    ("The GenericMassGenerator class is the main class for"
     " mass generation in Herwig++.");

  static Parameter<GenericMassGenerator,string> interfaceParticle
    ("Particle",
     "The particle data object for this class",
     0, "", true, false, 
     &GenericMassGenerator::setParticle, 
     &GenericMassGenerator::getParticle);


  static Switch<GenericMassGenerator,bool> interfaceInitialize
    ("Initialize",
     "Initialize the calculation of the maximum weight etc",
     &GenericMassGenerator::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialization
    (interfaceInitialize,
     "Yes",
     "Do the initialization",
     true);
  static SwitchOption interfaceInitializeNoInitialization
    (interfaceInitialize,
     "No",
     "Don't do the initalization",
     false);

  static Switch<GenericMassGenerator,int> interfaceBreitWignerShape
    ("BreitWignerShape",
     "Controls the shape of the mass distribution generated",
     &GenericMassGenerator::_BWshape, 0, false, false);
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

  static Parameter<GenericMassGenerator,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the unweighting",
     &GenericMassGenerator::_maxwgt, 1.0, 0.0, 1000.0,
     false, false, true);

  static Parameter<GenericMassGenerator,int> interfaceNGenerate
    ("NGenerate",
     "The number of tries to generate the mass",
     &GenericMassGenerator::_ngenerate, 100, 0, 10000,
     false, false, true);
 
  static Parameter<GenericMassGenerator,int> interfaceNInitial
    ("NInitial",
     "Number of tries for the initialisation",
     &GenericMassGenerator::_ninitial, 1000, 0, 100000,
     false, false, true);

  static Switch<GenericMassGenerator,bool> interfaceWidthOption
    ("WidthOption",
     "Which width to use",
     &GenericMassGenerator::_widthopt, false, false, false);
  static SwitchOption interfaceWidthOptionNominalWidth
    (interfaceWidthOption,
     "NominalWidth",
     "Use the normal width from the particle data object",
     false);
  static SwitchOption interfaceWidthOptionPhysicalWidth
    (interfaceWidthOption,
     "PhysicalWidth",
     "Use the width calculated at the on-shell mass",
     true);

}

bool GenericMassGenerator::accept(const ParticleData & in) const {
  if(!_particle) return false;
  return in.id() == _particle->id() || 
    ( _particle->CC() && _particle->CC()->id() == in.id() );
}

void GenericMassGenerator::doinit() throw(InitException) {
  MassGenerator::doinit();
  // get the antiparticle
  _antiparticle=_particle->CC();
  // the width generator
  _particle->init();
  _widthgen=_particle->widthGenerator();
  if(_widthgen){_widthgen->init();}
  // local storage of particle properties for speed
  _mass=_particle->mass();
  _width= _widthopt ? _particle->generateWidth(_mass) : _particle->width();
  _mass2=_mass*_mass;
  _mwidth=_mass*_width;
  _lowermass = _mass-_particle->widthLoCut();
  _uppermass = _mass+_particle->widthUpCut();
  // print out messagw if doing the initialisation
  if(_initialize) {
    // zero the maximum weight
    _maxwgt=0.;
    // storage of variables for the loop
    double wgt=0.,swgt=0.,sqwgt=0.;
    Energy mdummy;
    // perform the initialisation
    for(int ix=0;ix<_ninitial;++ix) {
      mdummy=mass(wgt,*_particle,3);
      swgt  += wgt;
      sqwgt += sqr(wgt);
      if(wgt>_maxwgt) _maxwgt=wgt;
    }
    swgt=swgt/_ninitial;
    sqwgt=sqrt(max(0.,sqwgt/_ninitial-swgt*swgt)/_ninitial);
  }
}
  
void GenericMassGenerator::dataBaseOutput(ofstream & output, bool header) {
  if(header) output << "update Mass_Generators set parameters=\"";
  output << "set " << fullName() << ":BreitWignerShape "   << _BWshape << "\n";
  output << "set " << fullName() << ":MaximumWeight " << _maxwgt    << "\n";
  output << "set " << fullName() << ":NGenerate "   << _ngenerate << "\n";
  output << "set " << fullName() << ":WidthOption " << _widthopt << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void GenericMassGenerator::dofinish() {
  if(_initialize) {
    string fname = CurrentGenerator::current().filename() + 
      string("-") + name() + string(".output");
    ofstream output(fname.c_str());
    dataBaseOutput(output,true);
  }
  MassGenerator::dofinish();
}

void GenericMassGenerator::setParticle(string p) {
  _particle = Repository::findParticle(StringUtils::basename(p));
  if ( ! _particle ) 
    Throw<InterfaceException>() 
      << "Could not set Particle interface "
      << "for the object \"" << name()
      << "\". Particle \"" << StringUtils::basename(p) << "\" not found."
      << Exception::runerror;
}

string GenericMassGenerator::getParticle() const {
  return _particle ? _particle->fullName() : "";
}

void GenericMassGenerator::rebind(const TranslationMap & trans)
  throw(RebindException) {
  _particle = trans.translate(_particle);
  _antiparticle = trans.translate(_antiparticle);
  MassGenerator::rebind(trans);
}

IVector GenericMassGenerator::getReferences() {
  IVector ret = MassGenerator::getReferences();
  ret.push_back(_particle);
  ret.push_back(_antiparticle);
  return ret;
}
