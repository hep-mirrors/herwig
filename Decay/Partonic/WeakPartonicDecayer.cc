// -*- C++ -*-
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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WeakPartonicDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

bool WeakPartonicDecayer::accept(const DecayMode & dm) const {
  // check we can find the flavours of the quarks in the decaying meson
  long id = dm.parent()->id();
  int flav1, flav2;
  if((id / 1000)%10) 
    {
      flav1 = (id/1000)%10;
      flav2 = (id/10)%100;
    } 
  else 
    {
      flav1 = id/100;  
      flav2 = (id/10)%10;
    }
  if(!flav1 || !flav2) return false;
  tcPDPtr prod[4];
  // if two decay products one must be in triplet and one antitriplet
  for(unsigned int ix=0;ix<4&&ix<dm.products().size();++ix)
    {prod[ix]=dm.orderedProducts()[ix];}
  if(dm.products().size()==2)
    {
      if((prod[0]->iColour()==PDT::Colour3&&prod[1]->iColour()==PDT::Colour3bar)||
	 (prod[0]->iColour()==PDT::Colour3bar&&prod[1]->iColour()==PDT::Colour3))
	{return true;}
    }
  else if(dm.products().size()==3)
    {
      if(((prod[0]->iColour()==PDT::Colour3   &&prod[2]->iColour()==PDT::Colour3bar)||
	  (prod[0]->iColour()==PDT::Colour3bar&&prod[2]->iColour()==PDT::Colour3))
	 &&prod[1]->iColour()==PDT::Colour8)
	{return true;}
    }
  else if(dm.products().size()==4)
    {
      // first two particles should be leptons or q qbar
      if((prod[0]->id()>=11&&prod[0]->id()<=16&&prod[1]->id()<=-11&&prod[1]->id()>=-16)||
	 (prod[1]->id()>=11&&prod[1]->id()<=16&&prod[0]->id()<=-11&&prod[0]->id()>=-16)||
	 (prod[0]->iColour()==PDT::Colour3    &&prod[1]->iColour()==PDT::Colour3bar   )||
	 (prod[1]->iColour()==PDT::Colour3    &&prod[0]->iColour()==PDT::Colour3bar   ))
	{
	  // third particle quark and fourth colour anti-triplet or
	  // thrid particle antiquark and fourth colour triplet
	  if((prod[2]->iColour()==PDT::Colour3bar&&prod[3]->iColour()==PDT::Colour3   )||
	     (prod[2]->iColour()==PDT::Colour3   &&prod[3]->iColour()==PDT::Colour3bar))
	    {return true;}
	}
      else{return false;}
    }
  return false;
}

ParticleVector WeakPartonicDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  // these products have the mass but should have constituent mass
  for(unsigned int ix=0;ix<children.size();++ix)
    {
      Energy cmass=(children[ix]->dataPtr())->constituentMass();
      children[ix]->set5Momentum(Lorentz5Momentum(0.,0.,0.,cmass,cmass));
    }
  // special for the gluon
  if(children[1]->id()==ParticleID::g)
    {
      Energy cmass(_globalParameters->effectiveGluonMass());
      children[1]->set5Momentum(Lorentz5Momentum(0.,0.,0.,cmass,cmass));
    }
  // 2-body decays
  if(children.size()==2)
    {
      double ctheta,phi;
      Lorentz5Momentum pout[2];
      for(unsigned int ix=0;ix<2;++ix){pout[ix].setMass(children[ix]->mass());}
      Kinematics::generateAngles(ctheta,phi);
      Kinematics::twoBodyDecay(parent.momentum(),pout[0].mass(),pout[1].mass(),
			       ctheta,phi,pout[0],pout[1]);
      for(unsigned int ix=0; ix<2;++ix) children[ix]->setMomentum(pout[ix]);
      if(children[0]->coloured()) 
	{
	  if(children[0]->id() > 0){children[0]->antiColourNeighbour(children[1]);}
	  else                     {children[0]->    colourNeighbour(children[1]);}
	}
    }
  // 3-body decays
  else if(children.size()==3)
    {
      // set masses of products
      Lorentz5Momentum pout[3],pin(parent.momentum());
      for(unsigned int ix=0;ix<3;++ix){pout[ix].setMass(children[ix]->mass());}
      double xs(children[2]->mass()/pin.e()),xb(1.-xs);
      pout[2]=xs*pin;
      // Get the particle quark that is decaying
      long idQ, idSpec;
      idSpec = children[2]->id();
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
      for(unsigned int ix=0; ix<3;++ix) children[ix]->setMomentum(pout[ix]);
      // make the colour connections
      // quark first
      if(children[0]->data().iColour()==PDT::Colour3)
	{
	  children[0]->antiColourNeighbour(children[1]);
	  children[1]->colourNeighbour(children[0]);
	  children[1]->antiColourNeighbour(children[2]);
	  children[2]->colourNeighbour(children[1]);
	}
      // antiquark first
      else
	{
	  children[0]->colourNeighbour(children[1]);
	  children[1]->antiColourNeighbour(children[0]);
	  children[1]->colourNeighbour(children[2]);
	  children[2]->antiColourNeighbour(children[1]);
	}
    }
  // 4-body decays
  else if(children.size()==4)
    {
      // set masses of products
      Lorentz5Momentum pout[4],pin(parent.momentum());
      for(unsigned int ix=0;ix<4;++ix){pout[ix].setMass(children[ix]->mass());}
      double xs(children[3]->mass()/pin.e()),xb(1.-xs);
      pout[3]=xs*pin;
      // Get the particle quark that is decaying
      long idQ, idSpec;
      idSpec = children[3]->id();
      idQ = (parent.id()/1000)%10;
      if(!idQ) idQ = (parent.id()/100)%10;
      // Now the odd case of a B_c where the c decays, not the b
      if(idSpec == idQ) idQ = (parent.id()/10)%10;
      // momentum of the decaying quark
      PPtr inter = getParticleData(idQ)->produceParticle(parent.momentum()*xb);
      // include matrix element
      if(MECode==100)
	{
	  // only the phase space factor don't bother with the W propagator
	  Kinematics::threeBodyDecay(inter->momentum(),pout[1],pout[0],pout[2],&VAWt);
	}
      // flat phase space
      else
	{Kinematics::threeBodyDecay(inter->momentum(),pout[1],pout[0],pout[2]);}
      // set the momenta of the decay products
      for(unsigned int ix=0; ix<4;++ix) children[ix]->setMomentum(pout[ix]);
      if(children[0]->coloured()) 
	{
	  if(children[0]->data().iColour()==PDT::Colour3)
	    {children[0]->antiColourNeighbour(children[1]);}
	  else
	    {children[0]->    colourNeighbour(children[1]);}
	}
      if(children[2]->data().iColour()==PDT::Colour3)
	{children[2]->antiColourNeighbour(children[3]);}
      else
	{children[2]->    colourNeighbour(children[3]);}
    }
  return children;
}


void WeakPartonicDecayer::persistentOutput(PersistentOStream & os) const {
  os << MECode << _globalParameters;
}

void WeakPartonicDecayer::persistentInput(PersistentIStream & is, int) {
  is >> MECode >> _globalParameters;
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

  static Reference<WeakPartonicDecayer,GlobalParameters> 
    interfaceGlobalParameters("GlobalParameters", 
		      "A reference to the GlobalParameters object", 
		      &Herwig::WeakPartonicDecayer::_globalParameters,
		      false, false, true, false);

}

double WeakPartonicDecayer::VAWt(double *temp) 
{return (temp[1]-temp[0])*(temp[0]-temp[2])*temp[3];}


}

