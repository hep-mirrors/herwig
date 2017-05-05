// -*- C++ -*-
//
// HeavyDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyDecayer class.
//

#include "HeavyDecayer.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/UseRandom.h>
#include <ThePEG/Utilities/Maths.h>

using namespace Herwig;

HeavyDecayer::HeavyDecayer() : MECode(0) {} 

IBPtr HeavyDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr HeavyDecayer::fullclone() const {
  return new_ptr(*this);
}

void HeavyDecayer::Init() {

  static ClassDocumentation<HeavyDecayer> documentation
    ("Class to decay all particles in HERWIG by the algorithms used in HERWIG 6.4");
  
  static Switch<HeavyDecayer,int> interfaceMECode
    ("MECode",
     "The code for the ME type to use in the decay",
     &HeavyDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Use a phase-space distribution",
     0);
  static SwitchOption interfaceMECodeWeak
    (interfaceMECode,
     "Weak",
     "Use the weak V-A matrix element",
     100);
}

ClassDescription<HeavyDecayer> HeavyDecayer::initHeavyDecayer;

bool HeavyDecayer::accept(tcPDPtr parent, const tPDVector & children) const { 
  long id = parent->id();
  int flav1, flav2;
  if((id / 1000)%10) {
    flav1 = (id/1000)%10;
    flav2 = (id/10)%100;
  } 
  else {
    flav1 = id/100;  
    flav2 = (id/10)%10;
  }
  if(!flav1 || !flav2) return false;
  return children.size()==4;
}

ParticleVector HeavyDecayer::decay(const Particle & p,
				   const tPDVector & children) const {
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
  }
  Lorentz5Momentum products[4];
  for(int i = 0; i<4; i++) products[i].setMass(partons[i]->mass());
  // Fraction of momentum to spectator
  double xs = partons[3]->mass()/p.momentum().e();
  // Fraction of momentum to heavy quark (which is weakly decaying)
  double xb = 1.-xs;
  partons[3]->setMomentum(p.momentum()*xs);
  // Get the particle that is decaying
  long idSpec = partons[3]->id();
  long idQ = (p.id()/1000)%10;
  if(!idQ) idQ = (p.id()/100)%10;
  // Now the odd case of a B_c where the c decays, not the b
  if(idSpec == idQ) idQ = (p.id()/10)%10;
  PPtr inter = getParticleData(idQ)->produceParticle(p.momentum()*xb);
  // Three Body Decay
  // Free Massless (V-A)*(V-A) ME
  if(MECode == 100) {
    Energy Mw = getParticleData(ParticleID::Wplus)->mass();
    Energy GamW = getParticleData(ParticleID::Wplus)->width();
    Energy2 EMwSq = sqr(Mw);
    Energy4 GMwSq = sqr(Mw*GamW);
    Energy4 EmLim = GMwSq + sqr(EMwSq - sqr(inter->mass()-p.mass()));
    Energy4 EmTest;
    do {
      Kinematics::threeBodyDecay(inter->momentum(),products[0], products[1], 
				 products[2], &VAWt);
      Energy2 pw2 = (products[0]+products[1]).m2();
      EmTest = sqr(pw2-EMwSq);
    } 
    while((EmTest+GMwSq)*rnd() > EmLim);
  }
  else Kinematics::threeBodyDecay(inter->momentum(),products[0], products[1],
				  products[2]);
  // set the momenta of the products
  for(int i = 0; i<3; i++) partons[i]->setMomentum(products[i]);
  // Set up colour connections based on the diagram above and input order
  if(partons[0]->coloured()) {
    if(partons[0]->id() > 0) partons[0]->antiColourNeighbour(partons[1]);
    else                     partons[0]->colourNeighbour(    partons[1]);
  }
  if(partons[2]->coloured()) {
    if(partons[2]->id() > 0) partons[2]->antiColourNeighbour(partons[3]);
    else                     partons[2]->colourNeighbour(partons[3]);
  }
  return partons;
}
   
void HeavyDecayer::persistentOutput(PersistentOStream &os) const {
  os << MECode; 
}

void HeavyDecayer::persistentInput(PersistentIStream &is, int) {
  is >> MECode; 
}

double HeavyDecayer::VAWt(Energy2 t0, Energy2 t1, Energy2 t2, InvEnergy4 t3) {
  return (t1-t0)*(t0-t2)*t3; 
}

void HeavyDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
