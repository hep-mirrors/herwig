// -*- C++ -*-
//
// GenericMassGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "GenericWidthGenerator.h"

using namespace Herwig;
using namespace ThePEG;

GenericMassGenerator::GenericMassGenerator()
  : maxWgt_(),BWShape_(0),nGenerate_(100),
    lowerMass_(),upperMass_(),
    mass_(),width_(),mass2_(),mWidth_(),
    nInitial_(1000),
    initialize_(false), output_(false), widthOpt_(false) {
}

GenericMassGenerator::~GenericMassGenerator() {
}

IBPtr GenericMassGenerator::clone() const {return new_ptr(*this);}
IBPtr GenericMassGenerator::fullclone() const {return new_ptr(*this);}

void GenericMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << particle_
     << ounit(lowerMass_,GeV) << ounit(upperMass_,GeV) << maxWgt_ 
     << BWShape_ << nGenerate_ 
     << ounit(mass_,GeV) << ounit(width_,GeV) 
     << ounit(mass2_,GeV2) << ounit(mWidth_,GeV2) 
     << nInitial_ << initialize_ << output_ << widthGen_ << widthGenB_
     << widthOpt_;
}

void GenericMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> particle_
     >> iunit(lowerMass_,GeV) >> iunit(upperMass_,GeV) >> maxWgt_ 
     >> BWShape_ >> nGenerate_ 
     >> iunit(mass_,GeV) >> iunit(width_ ,GeV)
     >> iunit(mass2_,GeV2) >> iunit(mWidth_ ,GeV2)
     >> nInitial_ >> initialize_ >> output_ >> widthGen_ >> widthGenB_
     >> widthOpt_;
}

ClassDescription<GenericMassGenerator> GenericMassGenerator::initGenericMassGenerator;
// Definition of the static class description member.


void GenericMassGenerator::Init() {

  static ClassDocumentation<GenericMassGenerator> documentation
    ("The GenericMassGenerator class is the main class for"
     " mass generation in Herwig.");

  static Parameter<GenericMassGenerator,string> interfaceParticle
    ("Particle",
     "The particle data object for this class",
     0, "", true, false, 
     &GenericMassGenerator::setParticle, 
     &GenericMassGenerator::getParticle);

  static Switch<GenericMassGenerator,bool> interfaceInitialize
    ("Initialize",
     "Initialize the calculation of the maximum weight etc",
     &GenericMassGenerator::initialize_, false, false, false);
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

  static Switch<GenericMassGenerator,bool> interfaceOutput
    ("Output",
     "Output the setup",
     &GenericMassGenerator::output_, false, false, false);
  static SwitchOption interfaceOutputYes
    (interfaceOutput,
     "Yes",
     "Output the data",
     true);
  static SwitchOption interfaceOutputNo
    (interfaceOutput,
     "No",
     "Don't output the data",
     false);

  static Switch<GenericMassGenerator,int> interfaceBreitWignerShape
    ("BreitWignerShape",
     "Controls the shape of the mass distribution generated",
     &GenericMassGenerator::BWShape_, 0, false, false);
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
     &GenericMassGenerator::maxWgt_, 1.0, 0.0, 1000.0,
     false, false, true);

  static Parameter<GenericMassGenerator,int> interfaceNGenerate
    ("NGenerate",
     "The number of tries to generate the mass",
     &GenericMassGenerator::nGenerate_, 100, 0, 10000,
     false, false, true);
 
  static Parameter<GenericMassGenerator,int> interfaceNInitial
    ("NInitial",
     "Number of tries for the initialisation",
     &GenericMassGenerator::nInitial_, 1000, 0, 100000,
     false, false, true);

  static Switch<GenericMassGenerator,bool> interfaceWidthOption
    ("WidthOption",
     "Which width to use",
     &GenericMassGenerator::widthOpt_, false, false, false);
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
  if(!particle_) return false;
  return in.id() == particle_->id() || 
    ( particle_->CC() && particle_->CC()->id() == in.id() );
}

void GenericMassGenerator::doinit() {
  MassGenerator::doinit();
  if(particle_->massGenerator()!=this) return;
  // the width generator
  particle_->init();
  widthGen_=particle_->widthGenerator();
  if(widthGen_) widthGen_->init();
  widthGenB_ = dynamic_ptr_cast<GenericWidthGeneratorPtr>(widthGen_);
  // local storage of particle properties for speed
  mass_=particle_->mass();
  width_= widthOpt_ ? particle_->generateWidth(mass_) : particle_->width();
  mass2_=mass_*mass_;
  mWidth_=mass_*width_;
  lowerMass_ = mass_-particle_->widthLoCut();
  upperMass_ = mass_+particle_->widthUpCut();
  // print out messagw if doing the initialisation
  if(initialize_) {
    // zero the maximum weight
    maxWgt_=0.;
    // storage of variable for the loop
    double wgt=0.;
    // perform the initialisation
    // (mass() contains a ThePEG::Random call)
    for(int ix=0; ix<nInitial_; ++ix) {
      mass(wgt, *particle_, 3);
      if ( wgt > maxWgt_ ) 
	maxWgt_ = wgt;
    }
  }
}
  
void GenericMassGenerator::dataBaseOutput(ofstream & output, bool header) {
  if(header) output << "update Mass_Generators set parameters=\"";
  output << "newdef " << fullName() << ":BreitWignerShape "   << BWShape_ << "\n";
  output << "newdef " << fullName() << ":MaximumWeight " << maxWgt_    << "\n";
  output << "newdef " << fullName() << ":NGenerate "   << nGenerate_ << "\n";
  output << "newdef " << fullName() << ":WidthOption " << widthOpt_ << "\n";
  if(header) output << "\n\" where BINARY ThePEGFullName=\"" 
		    << fullName() << "\";" << endl;
}

void GenericMassGenerator::dofinish() {
  if(output_) {
    string fname = CurrentGenerator::current().filename() + 
      string("-") + name() + string(".output");
    ofstream output(fname.c_str());
    dataBaseOutput(output,true);
  }
  MassGenerator::dofinish();
}

void GenericMassGenerator::setParticle(string p) {
  if ( (particle_ = Repository::GetPtr<tPDPtr>(p)) ) return;
  particle_ = Repository::findParticle(StringUtils::basename(p));
  if ( ! particle_ ) 
    Throw<InterfaceException>() 
      << "Could not set Particle interface "
      << "for the object \"" << name()
      << "\". Particle \"" << StringUtils::basename(p) << "\" not found."
      << Exception::runerror;
}

string GenericMassGenerator::getParticle() const {
  return particle_ ? particle_->fullName() : "";
}

void GenericMassGenerator::rebind(const TranslationMap & trans)
  {
  particle_ = trans.translate(particle_);
  MassGenerator::rebind(trans);
}

IVector GenericMassGenerator::getReferences() {
  IVector ret = MassGenerator::getReferences();
  ret.push_back(particle_);
  return ret;
}

pair<Energy,Energy> GenericMassGenerator::width(Energy q,int shape) const {
  Energy gam;
  if(shape==3 || BWShape_==1 || !widthGen_ ) {
    gam = width_;
  }
  else if(!widthGenB_) {
    gam = widthGen_->width(*particle_,q);
  }
  else {
    return widthGenB_->width(q,*particle_);
  }
  pair<Energy,Energy> output;
  output.first  = width_;
  output.second = width_;
  // take BR's into account
  double brSum=0.;
  for(DecaySet::const_iterator it = particle_->decayModes().begin();
      it!=particle_->decayModes().end();++it) {
    if((**it).on()) brSum += (**it).brat();
  }
  output.first *= brSum;
  return output;
}
