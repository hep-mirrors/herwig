// -*- C++ -*-
//
// SMHiggsWWDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsWWDecayer class.
//

#include "SMHiggsWWDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
typedef Selector<tDMPtr> DecaySelector;

ClassDescription<SMHiggsWWDecayer> SMHiggsWWDecayer::initSMHiggsWWDecayer;
// Definition of the static class description member.

void SMHiggsWWDecayer::Init() {

  static ClassDocumentation<SMHiggsWWDecayer> documentation
    ("The SMHiggsWWDecayer class performs the decay of the Standard Model Higgs"
     " boson to W+W- and Z0Z0");

  static ParVector<SMHiggsWWDecayer,double> interfaceWMaximum
    ("WMaximum",
     "The maximum weight for H-> W+W- decays",
     &SMHiggsWWDecayer::_wmax, 2, 1.0, 0.0, 10000.0,
     false, false, Interface::limited);

  static ParVector<SMHiggsWWDecayer,double> interfaceZMaximum
    ("ZMaximum",
     "The maximum weight for H-> Z0Z0 decays",
     &SMHiggsWWDecayer::_zmax, 2, 1.0, 0.0, 10000.0,
     false, false, Interface::limited);
}

SMHiggsWWDecayer::SMHiggsWWDecayer() : _wmax(2,1.00), _zmax(2,1.00)
{}

void SMHiggsWWDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) 
    throw InitException() << "SMHiggsWWDecayer needs the StandardModel class"
			  << " to be either the Herwig one or a class inheriting"
			  << " from it";
  _theFFWVertex = hwsm->vertexFFW();
  _theFFZVertex = hwsm->vertexFFZ();
  _theHVVVertex = hwsm->vertexWWH();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
  // the W+W- decays
  for(unsigned int ix=0;ix<2;++ix) {
    tPDPtr h0     = getParticleData(ParticleID::h0);
    tPDPtr wplus  = getParticleData(ParticleID::Wplus);
    tPDPtr wminus = getParticleData(ParticleID::Wminus);
    DecaySelector wpDecay =  wplus->decaySelector();
    DecaySelector wmDecay = wminus->decaySelector();
    tPDVector extpart(5);
    extpart[0]=h0;
    DecayPhaseSpaceModePtr mode;
    DecayPhaseSpaceChannelPtr newchannel;
    vector<double> wgt(1,1.);
    unsigned int imode=0;
    for(DecaySelector::const_iterator wp=wpDecay.begin();wp!=wpDecay.end();++wp) {
      // extract the decay products of W+
      tPDVector prod=(*wp).second->orderedProducts();
      if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
      extpart[1]=prod[0];
      extpart[2]=prod[1];
      for(DecaySelector::const_iterator wm=wmDecay.begin();wm!=wmDecay.end();++wm) {
	// extract the decay products of W-
	tPDVector prod=(*wm).second->orderedProducts();
	if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
	extpart[3]=prod[0];
	extpart[4]=prod[1];
	// create the new mode
	mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	// create the phase space channel
	newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	newchannel->addIntermediate(extpart[0],0,0.0,-1,-2);
	if(ix==0) {
	  newchannel->addIntermediate(wplus     ,1,0., 1, 2);
	  newchannel->addIntermediate(wminus    ,1,0., 3, 4);
	}
	else {
	  newchannel->addIntermediate(wplus     ,0,0., 1, 2);
	  newchannel->addIntermediate(wminus    ,0,0., 3, 4);
	}
	mode->addChannel(newchannel);
	addMode(mode,_wmax[ix],wgt);
	// insert mode into selector
	_ratio.push_back(wp->second->brat()*wm->second->brat());
	if(ix==0) _wdecays.insert (_ratio.back(),imode);
	++imode;
      }
    }
    // the Z0Z0 decays
    tPDPtr Z0=getParticleData(ParticleID::Z0);
    DecaySelector Z0Decay = Z0->decaySelector();
    for(DecaySelector::const_iterator z1=Z0Decay.begin();z1!=Z0Decay.end();++z1) {
      // extract the decay products of Z0
      tPDVector prod=(*z1).second->orderedProducts();
      if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
      extpart[1]=prod[0];
      extpart[2]=prod[1];
      for(DecaySelector::const_iterator z2=Z0Decay.begin();z2!=Z0Decay.end();++z2) {
	// extract the decay products of Z0
	tPDVector prod=(*z2).second->orderedProducts();
	if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
	extpart[3]=prod[0];
	extpart[4]=prod[1];
	// create the new mode
	mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	// create the phase space channel
	newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	newchannel->addIntermediate(extpart[0],0,0.0,-1,-2);
	if(ix==0) {
	  newchannel->addIntermediate(Z0    ,1,0., 1, 2);
	  newchannel->addIntermediate(Z0    ,1,0., 3, 4);
	}
	else {
	  newchannel->addIntermediate(Z0    ,0,0., 1, 2);
	  newchannel->addIntermediate(Z0    ,0,0., 3, 4);
	}
	mode->addChannel(newchannel);
	addMode(mode,_zmax[ix],wgt);
	// insert mode into selector
	_ratio.push_back(z1->second->brat()*z2->second->brat());
	if(ix==0) _zdecays.insert (_ratio.back(),imode);
	++imode;
      }
    }
  }
}

bool SMHiggsWWDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  // if not two decay products return false
  if(children.size()!=2) return false;
  // if not decaying higgs return false
  if(parent->id()!=ParticleID::h0) return false;
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if((id1==-id2&&abs(id1)==ParticleID::Wplus)||
     (id1== id2&&    id1 ==ParticleID::Z0))
    return true;
  else
    return false;
}

void SMHiggsWWDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _theFFZVertex << _theHVVVertex 
     << _wdecays << _zdecays << _ratio << _wmax << _zmax;
}

void SMHiggsWWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFWVertex >> _theFFZVertex >> _theHVVVertex 
     >> _wdecays >> _zdecays >> _ratio >> _wmax >> _zmax;
}

ParticleVector SMHiggsWWDecayer::decay(const Particle & parent,
				       const tPDVector & children) const {
  // select the decay modes of the gauge bosons
  unsigned int imode;
  if(abs(children[0]->id())==ParticleID::Wplus)
    imode=_wdecays.select(UseRandom::rnd());
  else
    imode=_zdecays.select(UseRandom::rnd());
  // use different phase space for low/high mass higgs
  if(parent.mass()>1.8*children[0]->mass()) 
    imode+=_wdecays.size()+_zdecays.size();
  // generate the kinematics
  return generate(true,false,imode,parent);
}

double SMHiggsWWDecayer::me2(const int, const Particle & inpart,
			     const ParticleVector & decay,
			     MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
					 PDT::Spin1Half,PDT::Spin1Half)));
  // check if Z or W decay
  bool Z0=decay[0]->id()==-decay[1]->id();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
    _swave = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_fwave1,decay[0],outgoing,true);
    SpinorWaveFunction   ::
      constructSpinInfo(_awave1,decay[1],outgoing,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_fwave2,decay[2],outgoing,true);
    SpinorWaveFunction   ::
      constructSpinInfo(_awave2,decay[3],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_fwave1,decay[0],outgoing);
  SpinorWaveFunction   ::
    calculateWaveFunctions(_awave1,decay[1],outgoing);
  SpinorBarWaveFunction::
    calculateWaveFunctions(_fwave2,decay[2],outgoing);
  SpinorWaveFunction   ::
    calculateWaveFunctions(_awave2,decay[3],outgoing);
  // get the intermediates and vertex
  tcPDPtr inter[2];
  AbstractFFVVertexPtr vert;
  if(Z0) {
    inter[0]=getParticleData(ParticleID::Z0);
    inter[1]=inter[0];
    vert=_theFFZVertex;
  }
  else {
    inter[0]=getParticleData(ParticleID::Wplus);
    inter[1]=getParticleData(ParticleID::Wminus);
    vert=_theFFWVertex;
  }


  // construct the spinors for the outgoing particles
  Energy2 scale0(sqr(inpart.mass()));
  Energy2 scale1((decay[0]->momentum()+decay[1]->momentum()).m2());
  Energy2 scale2((decay[2]->momentum()+decay[3]->momentum()).m2());
  // for decays to quarks ensure boson is massive enough to
  // put quarks on constituent mass-shell
  if(scale1<sqr(decay[0]->dataPtr()->constituentMass()+
		decay[1]->dataPtr()->constituentMass())) return 0.;
  if(scale2<sqr(decay[2]->dataPtr()->constituentMass()+
		decay[3]->dataPtr()->constituentMass())) return 0.;
  // compute the boson currents
  VectorWaveFunction curr1[2][2],curr2[2][2];
  unsigned int ohel1,ohel2,ohel3,ohel4;
  for(ohel1=0;ohel1<2;++ohel1) {
    for(ohel2=0;ohel2<2;++ohel2) {
      curr1[ohel1][ohel2]=vert->evaluate(scale1,1,inter[0],
					 _awave1[ohel2],_fwave1[ohel1]);
      curr2[ohel1][ohel2]=vert->evaluate(scale2,1,inter[1],
					 _awave2[ohel2],_fwave2[ohel1]);
    }
  }
  // compute the matrix element
  for(ohel1=0;ohel1<2;++ohel1) {
    for(ohel2=0;ohel2<2;++ohel2) {
      for(ohel3=0;ohel3<2;++ohel3) {
	for(ohel4=0;ohel4<2;++ohel4) {
	  (*ME())(0,ohel1,ohel2,ohel3,ohel4)=
	    _theHVVVertex->evaluate(scale0,curr1[ohel1][ohel2],
				    curr2[ohel3][ohel4],_swave);
	}
      }
    }
  }
  double output=(ME()->contract(_rho)).real()*scale0*UnitRemoval::InvE2;
  // set up the colour flows
  if(decay[0]->coloured()) {
    output*=3.;
    decay[0]->antiColourNeighbour(decay[1]);
  }
  if(decay[2]->coloured()) {
    output*=3.;
    decay[2]->antiColourNeighbour(decay[3]);
  }
  // divide out the gauge boson branching ratios
  output/=_ratio[imode()];
  // if Z0 decays identical particle factor
  if(Z0) output*=0.5;
  // return the answer
  return output;
}

void SMHiggsWWDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<2;++ix) {
    os << "newdef " << name() << ":WMaximum "    << ix << " " << _wmax[ix]  << "\n";
    os << "newdef " << name() << ":ZMaximum "    << ix << " " << _zmax[ix]  << "\n";
  }
  DecayIntegrator::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

void SMHiggsWWDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  for(unsigned int ix=0;ix<2;++ix) {
    _zmax[ix]=0.;
    _wmax[ix]=0.;
  }
  unsigned int ntest=_wdecays.size()+_zdecays.size();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      unsigned int iloc = ix<ntest ? 0 : 1;
      if(mode(ix)->externalParticles(1)->iCharge()==
	 -mode(ix)->externalParticles(2)->iCharge()) {
	_zmax[iloc]=max(mode(ix)->maxWeight(),_zmax[iloc]);
      }
      else {
	_wmax[iloc]=max(mode(ix)->maxWeight(),_wmax[iloc]);
      }
    }
  }
}
