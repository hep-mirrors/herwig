// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoCharmoniumDecayer class.
//

#include "BtoCharmoniumDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Kinematics.h"

using namespace Herwig;

inline IBPtr BtoCharmoniumDecayer::clone() const {
  return new_ptr(*this);
}

inline IBPtr BtoCharmoniumDecayer::fullclone() const {
  return new_ptr(*this);
}

BtoCharmoniumDecayer::BtoCharmoniumDecayer() {}

void BtoCharmoniumDecayer::persistentOutput(PersistentOStream & os) const {
}

void BtoCharmoniumDecayer::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<BtoCharmoniumDecayer> BtoCharmoniumDecayer::initBtoCharmoniumDecayer;
// Definition of the static class description member.

void BtoCharmoniumDecayer::Init() {

  static ClassDocumentation<BtoCharmoniumDecayer> documentation
    ("There is no documentation for the BtoCharmoniumDecayer class");

}

bool BtoCharmoniumDecayer::accept(tcPDPtr parent, const tPDVector & prod) const {
  // must be three decay products
  if(prod.size()!=3) return false;
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
  // first must be a singet
  if(prod[0]->iColour()!=PDT::Colour0) return false;
  // the second and third must be triplet/antitriplet
  if(!((prod[1]->iColour()==PDT::Colour3bar&&prod[2]->iColour()==PDT::Colour3   )||
       (prod[1]->iColour()==PDT::Colour3   &&prod[2]->iColour()==PDT::Colour3bar)))
    return false;
  // finally check enough phase space
  if(prod[0]->mass()+prod[1]->constituentMass()+prod[2]->constituentMass()
     >parent->mass()) return false;
  // if got here should be O.K.
  return true;
}

ParticleVector BtoCharmoniumDecayer::decay(const Particle & parent,
					  const tPDVector & children) const {
  // make the particles
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
    // these products have the mass but should have constituent mass
    partons[ix]->set5Momentum(Lorentz5Momentum(children[ix]->constituentMass()));
  }
  // get the momenta of the decaying quark and the spectator
  Lorentz5Momentum pin(parent.momentum());
  double xs(partons[2]->mass()/pin.e()),xb(1.-xs);
  Lorentz5Momentum pspect(xs*pin),pdec(xb*pin);
  pspect.setMass(partons[2]->mass());
  pdec.rescaleMass();
  // Get the particle quark that is decaying
  long idQ, idSpec;
  idSpec = partons[2]->id();
  idQ = (abs(parent.id())/1000)%10;
  if(!idQ) idQ = (abs(parent.id())/100)%10;
  // Now the odd case of a B_c where the c decays, not the b
  if(abs(idSpec) == idQ) idQ = (abs(parent.id())/10)%10;
  // change sign if spectator quark or antidiquark
  if((idSpec>0&&idSpec<6)||idSpec<-6) idQ = -idQ;
  vector<Lorentz5Momentum> pout(2,Lorentz5Momentum());
  double CosAngle, AzmAngle;
  Kinematics::generateAngles(CosAngle,AzmAngle);
  Kinematics::twoBodyDecay(pdec,partons[0]->mass(),partons[1]->mass(),
			   CosAngle,AzmAngle,pout[0],pout[1]);
  partons[0]->set5Momentum(pout[0]);
  partons[1]->set5Momentum(pout[1]);
  partons[2]->set5Momentum(pspect );
  if(partons[1]->dataPtr()->iColour()==PDT::Colour3)
    partons[1]->antiColourNeighbour(partons[2]);
  else
    partons[1]->    colourNeighbour(partons[2]);
  for(unsigned int ix=0;ix<partons.size();++ix) 
    cerr << *partons[ix] << "\n";
  return partons;
}

void BtoCharmoniumDecayer::dataBaseOutput(ofstream & output,
					 bool header) const {
}
