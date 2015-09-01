// -*- C++ -*-
//
// SMHiggsGGHiggsPPDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/Decay/GeneralDecayMatrixElement.h"


using namespace Herwig;
using namespace ThePEG::Helicity;

bool SMHiggsGGHiggsPPDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  int idp = parent->id();
  int id0 = children[0]->id();
  int id1 = children[1]->id();
  if((idp == ParticleID::h0 && 
      id0 == ParticleID::g     && id1 == ParticleID::g) ||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::gamma && id1 == ParticleID::gamma)) 
    return true;
  else
    return false;
}

ParticleVector SMHiggsGGHiggsPPDecayer::decay(const Particle & parent,
					      const tPDVector & children) const {
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

double SMHiggsGGHiggsPPDecayer::me2(const int, 
				    const Particle & part,
				    const ParticleVector & decay,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _swave = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					  incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::constructSpinInfo(_vwave[ix],decay[ix],
					    outgoing,true,true);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vwave[ix],decay[ix],outgoing,true);
  //Set up decay matrix
  Energy2 scale(sqr(part.mass()));
  unsigned int v1hel,v2hel;
  for(v1hel = 0;v1hel < 3;v1hel+=2) {
    for(v2hel = 0;v2hel < 3;v2hel+=2) {
      if(decay[0]->id() == ParticleID::g &&
	 decay[1]->id() == ParticleID::g) {
	(*ME())(0,v1hel,v2hel) = _hggvertex->evaluate(scale,_vwave[0][v1hel],
						   _vwave[1][v2hel],_swave);
      }
      else {
	(*ME())(0,v1hel,v2hel) = _hppvertex->evaluate(scale,_vwave[0][v1hel],
						   _vwave[1][v2hel],_swave);
      }
    }
  }
  //store matrix element
  double output = ME()->contract(_rho).real()*UnitRemoval::E2/scale;
  //colour factor (N^2 - 1)/4
  if(decay[0]->id() == ParticleID::g &&
     decay[1]->id() == ParticleID::g) {
    output *= 8.;
  }
  //symmetric final states
  output /= 2.;

  return output;
}

void SMHiggsGGHiggsPPDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
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
  tPDVector extpart(3);
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
