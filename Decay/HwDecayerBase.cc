// -*- C++ -*-
//
// HwDecayerBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayerBase class.
//

#include "HwDecayerBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/ShowerHandler.h"

using namespace Herwig;

bool HwDecayerBase::accept(const DecayMode & dm) const {
  // get the primary products
  tPDVector products=dm.orderedProducts();
  // add products for which the decay mode is all ready specified
  if(!dm.cascadeProducts().empty()) {
    for(ModeMSet::const_iterator mit=dm.cascadeProducts().begin();
	mit!=dm.cascadeProducts().end();++mit) {
      products.push_back(const_ptr_cast<PDPtr>((**mit).parent()));
    }
  }
  // can this mode be handled ?
  return accept(dm.parent(),products);
}

ParticleVector HwDecayerBase::decay(const DecayMode & dm,
				    const Particle & p) const {
  // handling of the decay including the special features of the
  // DecayMode  
  // get the primary products
  tPDVector products=dm.orderedProducts();
  // add products for which the decay mode is all ready specified
  if(!dm.cascadeProducts().empty()) {
    for(ModeMSet::const_iterator mit=dm.cascadeProducts().begin();
	mit!=dm.cascadeProducts().end();++mit) {
      products.push_back(const_ptr_cast<PDPtr>((**mit).parent()));
    }
  }
  // perform the primary decay
  ParticleVector output=decay(p,products);
  // perform the secondary decays
  if(!dm.cascadeProducts().empty()) {
    unsigned int iloc=dm.orderedProducts().size();
    for(ModeMSet::const_iterator mit=dm.cascadeProducts().begin();
	mit!=dm.cascadeProducts().end();++mit) {
      if(!(*mit)->decayer()) 
	throw Exception() << "Decay mode " << (**mit).tag() 
			  << "does not have a decayer, can't perform"
			  << "decay in  HwDecayerBase::decay()"
			  << Exception::eventerror;
      ParticleVector children=(*mit)->decayer()->decay(**mit,*output[iloc]);
      for(unsigned int ix=0;ix<children.size();++ix) {
	output[iloc]->addChild(children[ix]);
      }
      ++iloc;
    }
  }
  return output;
}

void HwDecayerBase::persistentOutput(PersistentOStream & os) const {
  os << _initialize << _dbOutput;
}

void HwDecayerBase::persistentInput(PersistentIStream & is, int) {
  is >> _initialize >> _dbOutput;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<HwDecayerBase,Decayer>
describeHerwigHwDecayerBase("Herwig::HwDecayerBase", "Herwig.so");

void HwDecayerBase::Init() {

  static ClassDocumentation<HwDecayerBase> documentation
    ("The HwDecayerBase class is the base class for Decayers in Hw++.");

  static Switch<HwDecayerBase,bool> interfaceInitialize
    ("Initialize",
     "Initialization of the phase space calculation",
     &HwDecayerBase::_initialize, false, false, false);
  static SwitchOption interfaceInitializeon
    (interfaceInitialize,
     "Yes",
     "At initialisation find max weight and optimise the integration",
     true);
  static SwitchOption interfaceInitializeoff
    (interfaceInitialize,
     "No",
     "Use the maximum weight and channel weights supplied for the integration",
     false);

  static Switch<HwDecayerBase,bool> interfaceDatabaseOutput
    ("DatabaseOutput",
     "Whether to print the database information",
     &HwDecayerBase::_dbOutput, false, false, false);
  static SwitchOption interfaceDatabaseOutputYes
    (interfaceDatabaseOutput,
     "Yes",
     "Output information on the decayer initialization",
     true);
  static SwitchOption interfaceDatabaseOutputNo
    (interfaceDatabaseOutput,
     "No",
     "Do not output information about the decayer initialization",
     false);
}

void HwDecayerBase::dofinish() {
  Decayer::dofinish();
  if(initialize() && databaseOutput()) {
    string fname = CurrentGenerator::current().filename() + 
      string("-") + name() + string(".output");
    ofstream output(fname.c_str());
    dataBaseOutput(output,true);
  }
}

bool HwDecayerBase::softMatrixElementVeto(PPtr , PPtr,
					  const bool &,
					  const Energy & ,
					  const vector<tcPDPtr> & ,
					  const double & ,
					  const Energy & ,
					  const Energy & ) {
  return false;
}

RealEmissionProcessPtr HwDecayerBase::generateHardest(RealEmissionProcessPtr) {
  return RealEmissionProcessPtr();
}

void HwDecayerBase::initializeMECorrection(RealEmissionProcessPtr , double & ,
			    double & ) {
  assert(false);
}

RealEmissionProcessPtr HwDecayerBase::applyHardMatrixElementCorrection(RealEmissionProcessPtr) {
  assert(false);
  return RealEmissionProcessPtr();
}

void HwDecayerBase::fixRho(RhoDMatrix & rho) const {
  if(ShowerHandler::currentHandlerIsSet() &&
     !ShowerHandler::currentHandler()->spinCorrelations())
    rho.reset();
}
