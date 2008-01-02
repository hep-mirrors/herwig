// -*- C++ -*-
//
// SMHiggsGGHiggsPPDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
  os << _hggvertex  << _hppvertex << _h0wgt << _hwidth;
}

void SMHiggsGGHiggsPPDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hggvertex >> _hppvertex >> _h0wgt >> _hwidth;
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
	higgs(0,v1hel,v2hel) = _hggvertex->evaluate(scale,V1[v1hel],
						    V2[v2hel],hwave);
      }
      else {
	higgs(0,v1hel,v2hel) = _hppvertex->evaluate(scale,V1[v1hel],
						    V2[v2hel],hwave);
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
  // normalize if width generator
  if(_hwidth) {
    if(decay[0]->id() == ParticleID::g) 
      output *=UnitRemoval::E/_hwidth->partialWidth(part.mass(),13);
    else
      output *=UnitRemoval::E/_hwidth->partialWidth(part.mass(),12);
  }
  return output;
}

void SMHiggsGGHiggsPPDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
  if(higgs->widthGenerator()) {
    _hwidth=dynamic_ptr_cast<SMHiggsWidthGeneratorPtr>(higgs->widthGenerator());
  }
  if(_hggvertex) {
    _hggvertex->init();
  }
  else {
    throw InitException() << "SMHiggsGGHiggsPPDecayer::doinit() - " 
			  << "_hggvertex is null";
  }
  if(_hppvertex) {
    _hppvertex->init();
  }
  else {
    throw InitException() << "SMHiggsGGHiggsPPDecayer::doinit() - " 
			  << "_hppvertex is null";
  }
  //set up decay modes
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3);
  vector<double> wgt(0);
  //glu,glu mode
  extpart[0] = getParticleData(ParticleID::h0);
  extpart[1] = getParticleData(ParticleID::g);
  extpart[2] = getParticleData(ParticleID::g);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[0],wgt);
  //gamma,gamma mode
  extpart[1] = getParticleData(ParticleID::gamma);
  extpart[2] = getParticleData(ParticleID::gamma);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[1],wgt);
}

void SMHiggsGGHiggsPPDecayer::doinitrun() {
  _hggvertex->initrun();
  _hppvertex->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _h0wgt[ix] = mode(ix)->maxWeight();
    }
  }
}
