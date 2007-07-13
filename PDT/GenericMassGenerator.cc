// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericMassGenerator class.
//

#include "GenericMassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GenericMassGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void GenericMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << _particle << _antiparticle 
     << ounit(_lowermass,GeV) << ounit(_uppermass,GeV) << _maxwgt 
     << _BWshape << _ngenerate 
     << ounit(_mass,GeV) << ounit(_width,GeV) 
     << ounit(_mass2,GeV2) << ounit(_mwidth,GeV2) 
     << _ninitial << _initialize << _widthgen;
}

void GenericMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _particle >> _antiparticle 
     >> iunit(_lowermass,GeV) >> iunit(_uppermass,GeV) >> _maxwgt 
     >> _BWshape >> _ngenerate 
     >> iunit(_mass,GeV) >> iunit(_width ,GeV)
     >> iunit(_mass2,GeV2) >> iunit(_mwidth ,GeV2)
     >> _ninitial >> _initialize >> _widthgen;
}

ClassDescription<GenericMassGenerator> GenericMassGenerator::initGenericMassGenerator;
// Definition of the static class description member.

void GenericMassGenerator::Init() {

  static ClassDocumentation<GenericMassGenerator> documentation
    ("The GenericMassGenerator class is the main class for"
     " mass generation in Herwig++.");

  static Reference<GenericMassGenerator,ParticleData> interfaceParticle
    ("Particle",
     "The particle data object for this class",
     &GenericMassGenerator::_particle, false, false, true, false, false);

  static Switch<GenericMassGenerator,bool> interfaceInitialize
    ("Initialize",
     "Initialize the calculation of the maximum weight etc",
     &GenericMassGenerator::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialization
    (interfaceInitialize,
     "Initialization",
     "Do the initialization",
     true);
  static SwitchOption interfaceInitializeNoInitialization
    (interfaceInitialize,
     "NoInitialization",
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

}


bool GenericMassGenerator::accept(const ParticleData & in) const
{
  if(!_particle){return false;}
  bool allowed=false;
  if(in.id()==_particle->id()){allowed=true;}
  else if(_particle->CC()&&(_particle->CC())->id()==in.id()){allowed=true;}
  return allowed;
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
  _width=_particle->width();
  _mass2=_mass*_mass;
  _mwidth=_mass*_width;
  _lowermass = _mass-_particle->widthLoCut();
  _uppermass = _mass+_particle->widthUpCut();
  // print out messagw if doing the initialisation
  if(_initialize)
    {
      // zero the maximum weight
      _maxwgt=0.;
      // storage of variables for the loop
      double wgt=0.,swgt=0.,sqwgt=0.;
      Energy mdummy;
      // perform the initialisation
      for(int ix=0;ix<_ninitial;++ix)
	{
	  mdummy=mass(*_particle,wgt,3);
	  swgt+=wgt;
	  sqwgt+=wgt*wgt;
	  if(wgt>_maxwgt){_maxwgt=wgt;}
	}
      swgt=swgt/_ninitial;
      sqwgt=sqrt(max(0.,sqwgt/_ninitial-swgt*swgt)/_ninitial);
    }
}

void GenericMassGenerator::dataBaseOutput(ofstream & output)
{
  output << "update Mass_Generators set parameters=\"";
  output << "set " << fullName() << ":BreitWignerShape "   << _BWshape << "\n";
  output << "set " << fullName() << ":MaximumWeight " << _maxwgt    << "\n";
  output << "set " << fullName() << ":NGenerate "   << _ngenerate << "\n";
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

}
