// -*- C++ -*-
//
// SMHiggsFermionsDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsFermionsDecayer class.
//

#include "SMHiggsFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Decay/DecayVertex.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMHiggsFermionsDecayer::SMHiggsFermionsDecayer() {
  _maxwgt.resize(9);
  _maxwgt[0]=0.;
  _maxwgt[1]=0;		
  _maxwgt[2]=0;		
  _maxwgt[3]=0.0194397;	
  _maxwgt[4]=0.463542;	
  _maxwgt[5]=0.;		
  _maxwgt[6]=6.7048e-09; 
  _maxwgt[7]=0.00028665; 
  _maxwgt[8]=0.0809643;  
}

void SMHiggsFermionsDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "SMHiggsFermionsDecayer needs the StandardModel class"
			  << " to be either the Herwig++ one or a class inheriting"
			  << " from it";
  _hvertex = hwsm->vertexFFH();
  // make sure they are initialized
  _hvertex->init();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
  if(higgs->widthGenerator()) {
    _hwidth=dynamic_ptr_cast<SMHiggsWidthGeneratorPtr>(higgs->widthGenerator());
  }
  // set up the decay modes
  vector<double> wgt(0);
  unsigned int imode=0;
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  int iy;
  extpart[0]=higgs;
  for(unsigned int istep=0;istep<11;istep+=10) {
    for(unsigned ix=1;ix<7;++ix) {
      if(istep<10||ix%2!=0) {
	iy = ix+istep;
	extpart[1]=getParticleData( iy);
	extpart[2]=getParticleData(-iy);
	mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	addMode(mode,_maxwgt[imode],wgt);
	++imode;
      }
    }
  }
}
  
bool SMHiggsFermionsDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  if(parent->id()!=ParticleID::h0||children.size()!=2) return false;
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if(id1==-id2&&(abs(id1)<=6||(abs(id1)>=11&&abs(id1)<=16)))
    return true;
  else
    return false;
}

ParticleVector SMHiggsFermionsDecayer::decay(const Particle & parent,
					     const tPDVector & children) const {
  // id's of the decaying particles
  tPDVector::const_iterator pit(children.begin());
  int id1((**pit).id());
  int imode=-1;
  if(abs(id1)<=6)                     imode = abs(id1)-1;
  else if(abs(id1)>=11&&abs(id1)<=16) imode = (abs(id1)-11)/2+6;
  ParticleVector output(generate(false,false,imode,parent));
  // set up the colour flow
  if(output[0]->hasColour())      output[0]->antiColourNeighbour(output[1]);
  else if(output[1]->hasColour()) output[1]->antiColourNeighbour(output[0]);
  return output;
}


void SMHiggsFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _maxwgt << _hvertex << _hwidth;
}

void SMHiggsFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _maxwgt >> _hvertex >> _hwidth;
}

ClassDescription<SMHiggsFermionsDecayer> SMHiggsFermionsDecayer::initSMHiggsFermionsDecayer;
// Definition of the static class description member.

void SMHiggsFermionsDecayer::Init() {

  static ClassDocumentation<SMHiggsFermionsDecayer> documentation
    ("The SMHiggsFermionsDecayer class implements the decat of the Standard Model"
     " Higgs boson to the Standard Model fermions");

  static ParVector<SMHiggsFermionsDecayer,double> interfaceMaxWeights
    ("MaxWeights",
     "Maximum weights for the various decays",
     &SMHiggsFermionsDecayer::_maxwgt, 9, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

}

// return the matrix element squared
double SMHiggsFermionsDecayer::me2(bool vertex, const int, const Particle & inpart,
				   const ParticleVector & decay) const {
  RhoDMatrix rhoin;
  // check if the incoming particle has a spin info and if not create it
  ScalarWaveFunction inwave = ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),
						 rhoin,incoming,true,vertex);
  // construct the spinors for the outgoing particles
  int iferm,ianti;
  if(decay[0]->id()<0){iferm=1;ianti=0;}
  else{iferm=0;ianti=1;}
  vector<SpinorWaveFunction   > awave;
  vector<SpinorBarWaveFunction> fwave;
  SpinorWaveFunction   (awave,decay[ianti],outgoing,true,vertex);
  SpinorBarWaveFunction(fwave,decay[iferm],outgoing,true,vertex);
  // compute the matrix element
  DecayMatrixElement newme(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half);
  Energy2 scale(inpart.mass()*inpart.mass());
  unsigned int ifm,ia;
  for(ifm=0;ifm<2;++ifm)
    {
      for(ia=0;ia<2;++ia)
	{
	  if(iferm>ianti)
	    {newme(0,ia,ifm)=_hvertex->evaluate(scale,awave[ia],fwave[ifm],inwave);}
	  else
	    {newme(0,ifm,ia)=_hvertex->evaluate(scale,awave[ia],fwave[ifm],inwave);}
	}
    }
  ME(newme);
  int id = abs(decay[0]->id());
  double output=(newme.contract(rhoin)).real()*UnitRemoval::E2/scale;
  if(id <=6) output*=3.;
  // normalize if width generator
  if(_hwidth) {
    if(id<=6) 
      output *= inpart.data().width()/_hwidth->partialWidth(inpart.mass(),id);
    else if(id>=11&&id<=15&&(id-9)%2==0) 
      output *= inpart.data().width()/_hwidth->partialWidth(inpart.mass(),(id+3)/2);
  }
  // test of the partial width
  /*
  Ptr<Herwig::StandardModel>::transient_const_pointer 
    hwsm=dynamic_ptr_cast<Ptr<Herwig::StandardModel>::transient_const_pointer>(standardModel());
  double g2(hwsm->alphaEM(scale)*4.*pi/hwsm->sin2ThetaW());
  Energy mass(hwsm->mass(scale,decay[0]->dataPtr())),
    mw(getParticleData(ParticleID::Wplus)->mass());
  double beta(sqrt(1.-4.*decay[0]->mass()*decay[0]->mass()/scale));
  Energy test(g2*mass*mass*beta*beta*beta*inpart.mass()/32./pi/mw/mw);
  if(abs(decay[0]->id())<=6){test *=3.;}
  cout << "testing the answer " << output << "     " 
       << test
       << endl;
  */
  return output;
}

void SMHiggsFermionsDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    os << "set " << fullName() << ":MaxWeights " << ix << " "
	   << _maxwgt[ix] << "\n";
  }
  DecayIntegrator::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

void SMHiggsFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _maxwgt[ix] = mode(ix)->maxWeight();
    }
  }
}
