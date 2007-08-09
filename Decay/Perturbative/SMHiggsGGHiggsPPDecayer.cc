// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsGGHiggsPPDecayer class.
//

#include "SMHiggsGGHiggsPPDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/PDT/DecayMode.h"


using namespace Herwig;
using namespace ThePEG::Helicity;

bool SMHiggsGGHiggsPPDecayer::accept(tcPDPtr parent, const PDVector & children) const {
  bool acc(false);
  int idp = parent->id();
  int id0 = children[0]->id();
  int id1 = children[1]->id();
  if((idp == ParticleID::h0 && id0 == ParticleID::g && 
      id1 == ParticleID::g)||
     (idp == ParticleID::h0 && id0 == ParticleID::gamma 
      && id1 == ParticleID::gamma)) 
    {acc = true;}
  return acc;
}

ParticleVector SMHiggsGGHiggsPPDecayer::decay(const Particle & parent,
					      const PDVector & children) const {
  int imode(0);
  if(children[0]->id() == ParticleID::gamma && 
     children[1]->id() == ParticleID::gamma)
    {imode = 1;}
  ParticleVector out(generate(true,false,imode,parent));
  //colour flow
  if(children[0]->id() == ParticleID::g &&
     children[1]->id() == ParticleID::g) {
    out[0]->colourNeighbour(out[1]);
    out[0]->antiColourNeighbour(out[1]);
  }
  return out;
}

void SMHiggsGGHiggsPPDecayer::persistentOutput(PersistentOStream & os) const {
  os << _hggvertex  << _hppvertex << _h0wgt;
}

void SMHiggsGGHiggsPPDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hggvertex >> _hppvertex >> _h0wgt;
}

ClassDescription<SMHiggsGGHiggsPPDecayer>
SMHiggsGGHiggsPPDecayer::initSMHiggsGGHiggsPPDecayer;

void SMHiggsGGHiggsPPDecayer::Init() {

  static ClassDocumentation<SMHiggsGGHiggsPPDecayer> documentation
    ("This is an implentation of h0->gg or h0->gamma,gamma "
     "decayer using the SMHGGVertex.");
  
  static Reference<SMHiggsGGHiggsPPDecayer,SMHGGVertex> 
    interfaceSMHGGVertex
    ("SMHGGVertex",
     "Pointer to SMHGGVertex",
     &SMHiggsGGHiggsPPDecayer::_hggvertex, false, false, true, 
     false, false);
  
  static Reference<SMHiggsGGHiggsPPDecayer,SMHPPVertex> 
    interfaceSMHPPVertex
    ("SMHPPVertex",
     "Pointer to SMHPPVertex",
     &SMHiggsGGHiggsPPDecayer::_hppvertex, false, false, true, 
     false, false);
  
}

double SMHiggsGGHiggsPPDecayer::me2(bool vertex, const int, 
				    const Particle & part,
				    const ParticleVector & decay) const {
  RhoDMatrix rhoH(PDT::Spin0);
  rhoH.average();
  vector<VectorWaveFunction> V1,V2;
  ScalarWaveFunction hwave(const_ptr_cast<tPPtr>(&part),
			   rhoH,incoming,true,vertex);
  VectorWaveFunction(V1,decay[0],outgoing,true,true,vertex);
  VectorWaveFunction(V2,decay[1],outgoing,true,true,vertex);
  //Set up decay matrix
  DecayMatrixElement higgs(PDT::Spin0,PDT::Spin1,PDT::Spin1);
  Energy2 scale(part.mass()*part.mass());
  unsigned int v1hel,v2hel;
  for(v1hel = 0;v1hel < 3;v1hel+=2) {
    for(v2hel = 0;v2hel < 3;v2hel+=2) {
      
      if(decay[0]->id() == ParticleID::g &&
	 decay[1]->id() == ParticleID::g) {
	higgs(0,v1hel,v2hel) = _hggvertex->evaluate(scale,hwave,
						    V1[v1hel],
						    V2[v2hel]);
      }
      else {
	higgs(0,v1hel,v2hel) = _hppvertex->evaluate(scale,hwave,
						    V1[v1hel],
						    V2[v2hel]);
      }
    }
  }
  //store matrix element
  ME(higgs);
  double output = higgs.contract(rhoH).real()*UnitRemoval::E2/scale;
  //colour factor (N^2 - 1)/4
  if(decay[0]->id() == ParticleID::g &&
     decay[1]->id() == ParticleID::g) {
    output *= 2.;
  }
  //symmetric final states
  output /= 2.;
  return output;
}

