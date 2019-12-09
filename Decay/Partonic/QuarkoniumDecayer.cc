// -*- C++ -*-
//
// QuarkoniumDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyDecayer class.
//

#include "QuarkoniumDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/UseRandom.h>
#include <cassert>

using namespace Herwig;

QuarkoniumDecayer::QuarkoniumDecayer() : MECode(0) {} 

IBPtr QuarkoniumDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr QuarkoniumDecayer::fullclone() const {
  return new_ptr(*this);
}

void QuarkoniumDecayer::Init() {
  
  static ClassDocumentation<QuarkoniumDecayer> documentation
    ("The QuarkoniumDecayer performs partonic decays of quarkonium"
     " resonances");
  
  static Switch<QuarkoniumDecayer,int> interfaceMECode
    ("MECode",
     "The code for the ME type to use in the decay",
     &QuarkoniumDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Use a phase-space distribution",
     0);
  static SwitchOption interfaceMECodeOrePowell
    (interfaceMECode,
     "OrePowell",
     "Use the Ore-Powell matrix element",
     130);
  
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QuarkoniumDecayer,PartonicDecayerBase>
describeHerwigQuarkoniumDecayer("Herwig::QuarkoniumDecayer", "HwPartonicDecay.so");

bool QuarkoniumDecayer::accept(tcPDPtr, const tPDVector & children) const {
  return (children.size() == 3 || children.size() == 2);
}

ParticleVector QuarkoniumDecayer::decay(const Particle & p,
					const tPDVector & children) const {
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
  }
  assert(partons.size()==2 || partons.size()==3);
  Lorentz5Momentum products[3];
  Energy gluMass = getParticleData(ParticleID::g)->constituentMass();
  for(unsigned int i = 0; i<partons.size(); i++) {
    if(partons[i]->id() == ParticleID::g) products[i].setMass(gluMass);
    else                                  products[i].setMass(partons[i]->mass());
  }
  if(partons.size() == 3) {
    // 3-gluon or 2-gluon + photon decay
    // Ore & Powell orthopositronium matrix element
    if(MECode == 130) { 
      double x1, x2, x3, test;
      do {
	// if decay fails return empty vector.
	if (! Kinematics::threeBodyDecay(p.momentum(), products[0],
					 products[1],  products[2]))
	  return ParticleVector();
	x1 = 2.*(p.momentum()*products[0])/sqr(p.mass());
	x2 = 2.*(products[1]*products[2])/sqr(p.mass());
	x3 = 2. - x1 - x2;
	test = sqr(x1*(1.-x1)) + sqr(x2*(1.-x2)) + sqr(x3*(1.-x3));
	test /= sqr(x1*x2*x3);
      } 
      while(test < 2.*UseRandom::rnd());
    }
    else {
      if (! Kinematics::threeBodyDecay(p.momentum(), products[0], 
				       products[1],products[2]))
	return ParticleVector();
    }
    // test the momenta
    for(unsigned int i = 0; i<partons.size(); i++)
      partons[i]->set5Momentum(products[i]);
    // Now set colour connections
    if(partons[2]->id() == ParticleID::g) {
      partons[0]->colourNeighbour(partons[1]);
      partons[1]->colourNeighbour(partons[2]);
      partons[2]->colourNeighbour(partons[0]);
    } 
    else {
      partons[0]->colourNeighbour(partons[1]);
      partons[0]->antiColourNeighbour(partons[1]);
    }
  } 
  // two decay children
  // 2 gluon or q-qbar decay
  else {
    double Theta, Phi;
    Kinematics::generateAngles(Theta, Phi);
    Energy p1 = partons[0]->mass();
    Energy p2 = partons[1]->mass();
    if(p1 == ZERO) p1 = gluMass;
    if(p2 == ZERO) p2 = gluMass;
    if (! Kinematics::twoBodyDecay(p.momentum(), p1, p2, Theta, Phi,
				   products[0], products[1]))
      return ParticleVector();
    for(unsigned int i = 0; i<partons.size(); i++)
      partons[i]->set5Momentum(products[i]);
    int first(0), second(1); 
    if(partons[0]->id() < 0) swap(first,second);
    partons[first]->antiColourNeighbour(partons[second]);
    if(abs(partons[first]->id()) == ParticleID::g) 
      partons[first]->colourNeighbour(partons[second]);
  }
  return partons;
}
   
void QuarkoniumDecayer::persistentOutput(PersistentOStream &os) const { 
  os << MECode;
}

void QuarkoniumDecayer::persistentInput(PersistentIStream &is, int) { 
  is >> MECode;
}

void QuarkoniumDecayer::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
