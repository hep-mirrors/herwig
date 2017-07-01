// -*- C++ -*-
//
// SMHiggsGGHiggsPPDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsGGHiggsPPDecayer class.
//

#include "SMHiggsGGHiggsPPDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<SMHiggsGGHiggsPPDecayer,DecayIntegrator>
describeHerwigSMHiggsGGHiggsPPDecayer("Herwig::SMHiggsGGHiggsPPDecayer",
				      "HwPerturbativeHiggsDecay.so");

bool SMHiggsGGHiggsPPDecayer::accept(tcPDPtr parent,
				       const tPDVector & children) const {
  int idp = parent->id();
  int id0 = children[0]->id();
  int id1 = children[1]->id();
  if((idp == ParticleID::h0 && 
      id0 == ParticleID::g     && id1 == ParticleID::g) ||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::gamma && id1 == ParticleID::gamma)||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::Z0 && id1 == ParticleID::gamma)||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::gamma && id1 == ParticleID::Z0)) 
    return true;
  else
    return false;
}

ParticleVector SMHiggsGGHiggsPPDecayer::decay(const Particle & parent,
					      const tPDVector & children) const {
  int imode(2);
  if(children[0]->id() == ParticleID::gamma && 
     children[1]->id() == ParticleID::gamma)
    imode = 1;
  else if(children[0]->id() ==ParticleID::g)
    imode = 0;
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
  os << _hggvertex << _hppvertex << _hzpvertex  << _h0wgt;
}

void SMHiggsGGHiggsPPDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hggvertex >> _hppvertex >> _hzpvertex >> _h0wgt;
}

ClassDescription<SMHiggsGGHiggsPPDecayer>
SMHiggsGGHiggsPPDecayer::initSMHiggsGGHiggsPPDecayer;

void SMHiggsGGHiggsPPDecayer::Init() {

  static ClassDocumentation<SMHiggsGGHiggsPPDecayer> documentation
    ("This is an implentation of h0->gg or h0->gamma,gamma "
     "decayer using the SMHGGVertex.");
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHGGVertex
    ("SMHGGVertex",
     "Pointer to SMHGGVertex",
     &SMHiggsGGHiggsPPDecayer::_hggvertex, false, false, true, 
     false, false);
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHPPVertex
    ("SMHPPVertex",
     "Pointer to SMHPPVertex",
     &SMHiggsGGHiggsPPDecayer::_hppvertex, false, false, true, 
     false, false);
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHZPVertex
    ("SMHZPVertex",
     "Pointer to SMHZPVertex",
     &SMHiggsGGHiggsPPDecayer::_hzpvertex, false, false, true, 
     false, false);
  
  static ParVector<SMHiggsGGHiggsPPDecayer,double> interfaceMaxWeights
    ("MaxWeights",
     "Maximum weights for the various decays",
     &SMHiggsGGHiggsPPDecayer::_h0wgt, 3, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

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
					    outgoing,true,
					    decay[ix]->id()!=ParticleID::Z0);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vwave[ix],decay[ix],outgoing,decay[ix]->id()!=ParticleID::Z0);
  //Set up decay matrix
  Energy2 scale(sqr(part.mass()));
  unsigned int v1hel,v2hel;
  AbstractVVSVertexPtr vertex;
  unsigned int vstep1(2),vstep2(2);
  double sym(1.);
  if(decay[0]->id() == ParticleID::g &&
     decay[1]->id() == ParticleID::g) {
    vertex = _hggvertex;
    sym = 2.;
  }
  else if(decay[0]->id() == ParticleID::gamma &&
	  decay[1]->id() == ParticleID::gamma) {
    vertex = _hppvertex;
    sym = 2.;
  }
  else if(decay[0]->id() == ParticleID::Z0 &&
	  decay[1]->id() == ParticleID::gamma) {
    vertex = _hzpvertex;
    vstep1 = 1;
  }
  else if(decay[1]->id() == ParticleID::Z0 &&
	  decay[0]->id() == ParticleID::gamma) {
    vertex = _hzpvertex;
    vstep2 = 1;
  }
  else
    assert(false);
  // loop over the helicities of the outgoing bosons
  for(v1hel = 0;v1hel < 3;v1hel+=vstep1) {
    for(v2hel = 0;v2hel < 3;v2hel+=vstep2) {
      (*ME())(0,v1hel,v2hel) = vertex->evaluate(scale,_vwave[0][v1hel],
						_vwave[1][v2hel],_swave);
    }
  }
  //store matrix element
  double output = ME()->contract(_rho).real()*UnitRemoval::E2/scale;
  //colour factor (N^2 - 1)/4
  if(decay[0]->id() == ParticleID::g) output *= 8.;
  //symmetric final states
  output /= sym;
  // return the answer
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
  if(_hzpvertex) {
    _hzpvertex->init();
  }
  else {
    throw InitException() << "SMHiggsGGHiggsZPDecayer::doinit() - " 
			  << "_hzpvertex is null";
  }
  //set up decay modes
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  vector<double> wgt(0);
  // glu,glu mode
  extpart[0] = getParticleData(ParticleID::h0);
  extpart[1] = getParticleData(ParticleID::g);
  extpart[2] = getParticleData(ParticleID::g);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[0],wgt);
  // gamma,gamma mode
  extpart[1] = getParticleData(ParticleID::gamma);
  extpart[2] = getParticleData(ParticleID::gamma);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[1],wgt);
  // Z0,gamma mode
  extpart[1] = getParticleData(ParticleID::Z0);
  extpart[2] = getParticleData(ParticleID::gamma);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,_h0wgt[2],wgt);
}

void SMHiggsGGHiggsPPDecayer::doinitrun() {
  _hggvertex->initrun();
  _hppvertex->initrun();
  _hzpvertex->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _h0wgt[ix] = mode(ix)->maxWeight();
    }
  }
}

void SMHiggsGGHiggsPPDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  for(unsigned int ix=0;ix<_h0wgt.size();++ix) {
    os << "newdef " << name() << ":MaxWeights " << ix << " "
	   << _h0wgt[ix] << "\n";
  }
  os << "newdef " << name() << ":SMHGGVertex " << _hggvertex->fullName() << "\n";
  os << "newdef " << name() << ":SMHPPVertex " << _hppvertex->fullName() << "\n";
  os << "newdef " << name() << ":SMHZPVertex " << _hzpvertex->fullName() << "\n";
  DecayIntegrator::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}
