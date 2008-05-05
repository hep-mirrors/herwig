// -*- C++ -*-
//
// WeakPartonicDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakPartonicDecayer class.
//

#include "WeakPartonicDecayer.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ConstituentParticleData.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

bool WeakPartonicDecayer::accept(tcPDPtr parent, const PDVector & prod) const {
  // check we can find the flavours of the quarks in the decaying meson
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
  // if two decay products one must be in triplet and one antitriplet
  if(prod.size()==2) {
    if((prod[0]->iColour()==PDT::Colour3&&prod[1]->iColour()==PDT::Colour3bar)||
       (prod[0]->iColour()==PDT::Colour3bar&&prod[1]->iColour()==PDT::Colour3))
      return true;
  }
  else if(prod.size()==3) {
    if(((prod[0]->iColour()==PDT::Colour3   &&prod[2]->iColour()==PDT::Colour3bar)||
	(prod[0]->iColour()==PDT::Colour3bar&&prod[2]->iColour()==PDT::Colour3))
       &&prod[1]->iColour()==PDT::Colour8) return true;
  }
  else if(prod.size()==4) {
    // first two particles should be leptons or q qbar
    if((prod[0]->id()>=11&&prod[0]->id()<=16&&prod[1]->id()<=-11&&prod[1]->id()>=-16)||
       (prod[1]->id()>=11&&prod[1]->id()<=16&&prod[0]->id()<=-11&&prod[0]->id()>=-16)||
       (prod[0]->iColour()==PDT::Colour3    &&prod[1]->iColour()==PDT::Colour3bar   )||
       (prod[1]->iColour()==PDT::Colour3    &&prod[0]->iColour()==PDT::Colour3bar   )) {
      // third particle quark and fourth colour anti-triplet or
      // thrid particle antiquark and fourth colour triplet
      if((prod[2]->iColour()==PDT::Colour3bar&&prod[3]->iColour()==PDT::Colour3   )||
	 (prod[2]->iColour()==PDT::Colour3   &&prod[3]->iColour()==PDT::Colour3bar))
	return true;
    }
    else return false;
  }
  return false;
}

ParticleVector WeakPartonicDecayer::decay(const Particle & parent,
					  const PDVector & children) const {
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
  }
  // these products have the mass but should have constituent mass
  for(unsigned int ix=0;ix<partons.size();++ix) {
    Energy cmass=partons[ix]->dataPtr()->constituentMass();
    partons[ix]->set5Momentum(Lorentz5Momentum(cmass));
  }
  // 2-body decays
  if(partons.size()==2) {
    double ctheta,phi;
    Lorentz5Momentum pout[2];
    for(unsigned int ix=0;ix<2;++ix) pout[ix].setMass(partons[ix]->mass());
    Kinematics::generateAngles(ctheta,phi);
    Kinematics::twoBodyDecay(parent.momentum(),pout[0].mass(),pout[1].mass(),
			     ctheta,phi,pout[0],pout[1]);
    for(unsigned int ix=0; ix<2;++ix) partons[ix]->setMomentum(pout[ix]);
    if(partons[0]->coloured()) {
      if(partons[0]->id() > 0) partons[0]->antiColourNeighbour(partons[1]);
      else                     partons[0]->    colourNeighbour(partons[1]);
    }
  }
  // 3-body decays
  else if(partons.size()==3) {
    // set masses of products
    Lorentz5Momentum pout[3],pin(parent.momentum());
    for(unsigned int ix=0;ix<3;++ix){pout[ix].setMass(partons[ix]->mass());}
    double xs(partons[2]->mass()/pin.e()),xb(1.-xs);
    pout[2]=xs*pin;
    // Get the particle quark that is decaying
    long idQ, idSpec;
    idSpec = partons[2]->id();
    idQ = (parent.id()/1000)%10;
    if(!idQ) idQ = (parent.id()/100)%10;
    // Now the odd case of a B_c where the c decays, not the b
    if(idSpec == idQ) idQ = (parent.id()/10)%10;
    // momentum of the decaying quark
    PPtr inter = getParticleData(idQ)->produceParticle(parent.momentum()*xb);
    // two body decay of heavy quark
    double ctheta,phi;
    Kinematics::generateAngles(ctheta,phi);
    Kinematics::twoBodyDecay(inter->momentum(),pout[0].mass(),pout[1].mass(),
			     ctheta,phi,pout[0],pout[1]);
    // set the momenta of the decay products
    for(unsigned int ix=0; ix<3;++ix) partons[ix]->setMomentum(pout[ix]);
    // make the colour connections
    // quark first
    if(partons[0]->data().iColour()==PDT::Colour3) {
      partons[0]->antiColourNeighbour(partons[1]);
      partons[1]->colourNeighbour(partons[0]);
      partons[1]->antiColourNeighbour(partons[2]);
      partons[2]->colourNeighbour(partons[1]);
    }
    // antiquark first
    else {
      partons[0]->colourNeighbour(partons[1]);
      partons[1]->antiColourNeighbour(partons[0]);
      partons[1]->colourNeighbour(partons[2]);
      partons[2]->antiColourNeighbour(partons[1]);
    }
  }
  // 4-body decays
  else if(partons.size()==4) {
    // set masses of products
    Lorentz5Momentum pout[4],pin(parent.momentum());
    for(unsigned int ix=0;ix<4;++ix) pout[ix].setMass(partons[ix]->mass());
    double xs(partons[3]->mass()/pin.e()),xb(1.-xs);
    pout[3]=xs*pin;
    // Get the particle quark that is decaying
    long idQ, idSpec;
    idSpec = partons[3]->id();
    idQ = (abs(parent.id())/1000)%10;
    if(!idQ) idQ = (abs(parent.id())/100)%10;
    // Now the odd case of a B_c where the c decays, not the b
    if(abs(idSpec) == idQ) idQ = (abs(parent.id())/10)%10;
    // change sign if antiquark
    if((idSpec>0&&idSpec<6)||idSpec<-6) idQ = -idQ;
    // momentum of the decaying quark
    PPtr inter = getParticleData(idQ)->produceParticle(parent.momentum()*xb);
    // include matrix element
    // only the phase space factor don't bother with the W propagator
    if(MECode==100) {
      // charges of the exchanged W boson
      int c1 = inter->data().iCharge()-partons[2]->data().iCharge();
      int c2 = partons[0]->data().iCharge()+partons[1]->data().iCharge();
      // normal decay
      if(c1==c2&&abs(c1)==3)
	Kinematics::threeBodyDecay(inter->momentum(),pout[1],pout[0],pout[2],&VAWt);
      // colour rearranged decay
      else {
	int c3 = inter->data().iCharge()-partons[1]->data().iCharge();
	int c4 = partons[0]->data().iCharge()+partons[2]->data().iCharge();
	if(c3==c4&&abs(c3)==3)
	  Kinematics::threeBodyDecay(inter->momentum(),pout[2],pout[0],pout[1],&VAWt);
	else {
	  generator()->log() << "Unknown order for colour rearranged decay"
			     << " in WeakPartonicDecayer::decay()\n";
	  generator()->log() << c1 << " " << c2 << " " << c3 << " " << c4 << "\n";
	  generator()->log() << parent << "\n" << *inter << "\n";
	  for(unsigned int ix=0;ix<4;++ix) generator()->log() << *partons[ix] << "\n";
	  throw Exception()  << "Unknown order for colour rearranged decay"
			     << " in WeakPartonicDecayer::decay() "
			     << Exception::runerror;
	}
      }
    }
    // flat phase space
    else
      Kinematics::threeBodyDecay(inter->momentum(),pout[1],pout[0],pout[2]);
    // set the momenta of the decay products
    for(unsigned int ix=0; ix<4;++ix) partons[ix]->setMomentum(pout[ix]);
    if(partons[0]->coloured()) {
      if(partons[0]->data().iColour()==PDT::Colour3)
	partons[0]->antiColourNeighbour(partons[1]);
      else
	partons[0]->    colourNeighbour(partons[1]);
    }
    if(partons[2]->data().iColour()==PDT::Colour3) {
      partons[2]->antiColourNeighbour(partons[3]);
    }
    else {
      partons[2]->    colourNeighbour(partons[3]);
    }
  }
  return partons;
}
  
void WeakPartonicDecayer::persistentOutput(PersistentOStream & os) const {
  os << MECode;
}

void WeakPartonicDecayer::persistentInput(PersistentIStream & is, int) {
  is >> MECode;
}

ClassDescription<WeakPartonicDecayer> WeakPartonicDecayer::initWeakPartonicDecayer;
// Definition of the static class description member.

void WeakPartonicDecayer::Init() {

  static ClassDocumentation<WeakPartonicDecayer> documentation
    ("The WeakPartonicDecayer class performs partonic decays of hadrons containing a "
     "heavy quark.");

  static Switch<WeakPartonicDecayer,int> interfaceMECode
    ("MECode",
     "The code for the type of matrix element to be used.",
     &WeakPartonicDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Phase space decays",
     0);
  static SwitchOption interfaceMECodeWeak
    (interfaceMECode,
     "Weak",
     "Weak matrix element",
     100);

}

double WeakPartonicDecayer::VAWt(Energy2 t0, Energy2 t1, Energy2 t2, InvEnergy4 t3) {
  return (t1-t0)*(t0-t2)*t3; 
}

void WeakPartonicDecayer::dataBaseOutput(ofstream & output,
					 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "set " << fullName() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

