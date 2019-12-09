// -*- C++ -*-
//
// BtoSGammaDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaDecayer class.
//

#include "BtoSGammaDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr BtoSGammaDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr BtoSGammaDecayer::fullclone() const {
  return new_ptr(*this);
}

bool BtoSGammaDecayer::accept(tcPDPtr , const tPDVector & children) const {
  // should be three decay products
  if(children.size()!=3) return false;
  // photon should be last
  if(children[2]->id()!=ParticleID::gamma) return false;
  // strange should be first
  if(abs(children[0]->id())!=ParticleID::s) return false;
  // first and second should form a colour singlet
  if((children[0]->iColour()==PDT::Colour3&&
      children[1]->iColour()==PDT::Colour3bar)||
     (children[1]->iColour()==PDT::Colour3&&
      children[0]->iColour()==PDT::Colour3bar)) return true;
  else return false;
}

ParticleVector BtoSGammaDecayer::decay(const Particle & parent,
				       const tPDVector & prod) const {
  ParticleVector children;
  for(unsigned int ix=0;ix<prod.size();++ix) {
    children.push_back(prod[ix]->produceParticle());
  }
  // momenta of the decay products
  Lorentz5Momentum pout[3],phad;
  pout[0].setMass(children[0]->dataPtr()->constituentMass());
  pout[1].setMass(children[1]->dataPtr()->constituentMass());
  pout[2].setMass(ZERO);
  // first calculate the hadronic mass spectrum
  phad.setMass(_hadronicmass->hadronicMass(parent.mass(),pout[0].mass()+pout[1].mass()));
  // two body decay to hadronic cluster and photon
  double ctheta,phi;
  Kinematics::generateAngles(ctheta,phi);
  Kinematics::twoBodyDecay(parent.momentum(),pout[2].mass(),phad.mass(),
			   ctheta,phi,pout[2],phad);
  // two body decay of the cluster
  Kinematics::generateAngles(ctheta,phi);
  Kinematics::twoBodyDecay(phad,pout[0].mass(),pout[1].mass(),
			   ctheta,phi,pout[0],pout[1]);
  // set momenta of decay products
  for(unsigned int ix=0;ix<3;++ix){children[ix]->setMomentum(pout[ix]);}
  // make the colour connections
  // quark first
  if(children[0]->data().iColour()==PDT::Colour3) {
    children[0]->antiColourNeighbour(children[1]);
    children[1]->colourNeighbour(children[0]);
  }
  // antiquark first
  else {
    children[0]->colourNeighbour(children[1]);
    children[1]->antiColourNeighbour(children[0]);
  }
  return children;
}


void BtoSGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _hadronicmass;
}

void BtoSGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hadronicmass;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BtoSGammaDecayer,PartonicDecayerBase>
describeHerwigBtoSGammaDecayer("Herwig::BtoSGammaDecayer", "HwPartonicDecay.so");

void BtoSGammaDecayer::Init() {

  static ClassDocumentation<BtoSGammaDecayer> documentation
    ("The BtoSGammaDecayer class performs to the exclusive decay B to s gamma");

  static Reference<BtoSGammaDecayer,BtoSGammaHadronicMass> interfaceHadronicMass
    ("HadronicMass",
     "Pointer to the object computing the hadronic mass spectrum.",
     &BtoSGammaDecayer::_hadronicmass, false, false, true, false, false);

}

void BtoSGammaDecayer::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  _hadronicmass->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":HadronicMass " 
	 << _hadronicmass->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
