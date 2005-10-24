// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHDecayer class.
//

#include "SMHDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SMHDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/PDT/DecayMode.h"


namespace Herwig {
  using namespace ThePEG;
  using ThePEG::Helicity::RhoDMatrix;
  using Helicity::VectorWaveFunction;
  using Helicity::ScalarWaveFunction;
  using Helicity::Direction;
  using Helicity::incoming;
  using Helicity::outgoing;
  
  SMHDecayer::~SMHDecayer() {}

  bool SMHDecayer::accept(const DecayMode & dm) const {
    bool acc(false);
    int idp = dm.parent()->id();
    int id0 = dm.orderedProducts()[0]->id();
    int id1 = dm.orderedProducts()[1]->id();
    if((idp == ParticleID::h0 && id0 == ParticleID::g && 
	id1 == ParticleID::g)||
       (idp == ParticleID::h0 && id0 == ParticleID::gamma 
	&& id1 == ParticleID::gamma)) 
      {acc = true;}
    return acc;
  }
  
  ParticleVector SMHDecayer::decay(const DecayMode & dm,
				     const Particle & parent) const {
    int imode(0);
    if(dm.orderedProducts()[0]->id() == ParticleID::gamma &&
       dm.orderedProducts()[1]->id() == ParticleID::gamma)
      {imode = 1;}
    ParticleVector out(generate(true,false,imode,parent));
    //colour flow
    if(dm.orderedProducts()[0]->id() == ParticleID::g &&
       dm.orderedProducts()[1]->id() == ParticleID::g) {
    out[0]->colourNeighbour(out[1]);
    out[0]->antiColourNeighbour(out[1]);
    }
    return out;
  }
  
  void SMHDecayer::persistentOutput(PersistentOStream & os) const {
    os << _HGGVertex  << _HPPVertex << _Hwgt;
  }
  
  void SMHDecayer::persistentInput(PersistentIStream & is, int) {
    is >> _HGGVertex >> _HPPVertex >> _Hwgt;
  }
  
  ClassDescription<SMHDecayer>SMHDecayer::initSMHDecayer;

  void SMHDecayer::Init() {
       static ClassDocumentation<SMHDecayer> documentation
      ("This is an implentation of h0->gg or h0->gamma,gamma "
       "decayer using the SMHGGVertex.");

     static Reference<SMHDecayer,Helicity::SMHGGVertex> interfaceSMHGGVertex
       ("SMHGGVertex",
        "Pointer to SMHGGVertex",
        &SMHDecayer::_HGGVertex, false, false, true, true, false);
     

     static Reference<SMHDecayer,Helicity::SMHPPVertex> interfaceSMHPPVertex
       ("SMHPPVertex",
	"Pointer to SMHPPVertex",
	&SMHDecayer::_HPPVertex, false, false, true, true, false);

  }

  double SMHDecayer::me2(bool vertex, const int ichan, 
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
	  higgs(0,v1hel,v2hel) = _HGGVertex->evaluate(scale,hwave,
						      V1[v1hel],
						      V2[v2hel]);
	}
	else {
	  higgs(0,v1hel,v2hel) = _HPPVertex->evaluate(scale,hwave,
						      V1[v1hel],
						      V2[v2hel]);
	}
      }
    }
    //store matrix element
    ME(higgs);
    double output = higgs.contract(rhoH).real()/scale;
    //colour factor (N^2 - 1)/4
    if(decay[0]->id() == ParticleID::g &&
       decay[1]->id() == ParticleID::g) {
      output *= 2.;
    }
    //symmetric final states
    output /= 2.;
    return output;
  }
}
