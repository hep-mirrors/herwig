// -*- C++ -*-
//
// SMHiggsWWDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsWWDecayer class.
//

#include "SMHiggsWWDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
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

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHiggsWWDecayer,PerturbativeDecayer>
describeHerwigSMHiggsWWDecayer("Herwig::SMHiggsWWDecayer", "HwPerturbativeHiggsDecay.so");

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
  PerturbativeDecayer::doinit();
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
    tPDPtr wplus  = getParticleData(ParticleID::Wplus);
    tPDPtr wminus = getParticleData(ParticleID::Wminus);
    DecaySelector wpDecay =  wplus->decaySelector();
    DecaySelector wmDecay = wminus->decaySelector();
    unsigned int imode=0;
    for(DecaySelector::const_iterator wp=wpDecay.begin();wp!=wpDecay.end();++wp) {
      // extract the decay products of W+
      tPDVector prod=(*wp).second->orderedProducts();
      if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
      tPDVector out={prod[0],prod[1],tPDPtr(),tPDPtr()};
      for(DecaySelector::const_iterator wm=wmDecay.begin();wm!=wmDecay.end();++wm) {
	// extract the decay products of W-
	tPDVector prod=(*wm).second->orderedProducts();
	if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
	out[2] = prod[0];
	out[3] = prod[1];
	// create the new mode
	PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(higgs,out,_wmax[ix]));
	// create the phase space channel
	PhaseSpaceChannel phase((PhaseSpaceChannel(mode),0,wplus,0,wminus,1,1,1,2,2,3,2,4));
	mode->addChannel(phase);
	if(ix==0) {
	  phase.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	  phase.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	}
	addMode(mode);
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
      tPDVector out = {prod[0],prod[1],tPDPtr(),tPDPtr()};
      for(DecaySelector::const_iterator z2=Z0Decay.begin();z2!=Z0Decay.end();++z2) {
	// extract the decay products of Z0
	tPDVector prod=(*z2).second->orderedProducts();
	if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
	out[2] = prod[0];
	out[3] = prod[1];
	// create the new mode
	PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(higgs,out,_zmax[ix]));
	// create the phase space channel
	PhaseSpaceChannel phase((PhaseSpaceChannel(mode),0,Z0,0,Z0,1,1,1,2,2,3,2,4));
	mode->addChannel(phase);
	if(ix==0) {
	  phase.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	  phase.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	}
	addMode(mode);
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
  ParticleVector decay = generate(true,false,imode,parent);
  // set up the colour flows
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->children()[0]->coloured()) {
      decay[ix]->children()[0]->antiColourNeighbour(decay[ix]->children()[1]);
    }
  }
  return decay;
}

void SMHiggsWWDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
  					  incoming,true);
  SpinorBarWaveFunction::
    constructSpinInfo(_fwave1,decay[0],outgoing,true);
  SpinorWaveFunction   ::
    constructSpinInfo(_awave1,decay[1],outgoing,true);
  SpinorBarWaveFunction::
    constructSpinInfo(_fwave2,decay[2],outgoing,true);
  SpinorWaveFunction   ::
    constructSpinInfo(_awave2,decay[3],outgoing,true);
}

double SMHiggsWWDecayer::me2(const int,const Particle & part,
				   const tPDVector & outgoing,
				   const vector<Lorentz5Momentum> & momenta,
				   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
  					 PDT::Spin1Half,PDT::Spin1Half)));
  // check if Z or W decay
  bool Z0=outgoing[0]->id()==-outgoing[1]->id();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _swave = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(_rho);
  }

  SpinorBarWaveFunction fw1(momenta[0],outgoing[0],Helicity::outgoing);
  SpinorBarWaveFunction fw2(momenta[2],outgoing[2],Helicity::outgoing);
  SpinorWaveFunction    aw1(momenta[1],outgoing[1],Helicity::outgoing);
  SpinorWaveFunction    aw2(momenta[3],outgoing[3],Helicity::outgoing);
  _fwave1.resize(2);
  _fwave2.resize(2);
  _awave1.resize(2);
  _awave2.resize(2);  
  for(unsigned int ihel=0;ihel<2;++ihel) {
    fw1.reset(ihel);   _fwave1[ihel] = fw1;
    fw2.reset(ihel);   _fwave2[ihel] = fw2;
    aw1.reset(ihel);   _awave1[ihel] = aw1;
    aw2.reset(ihel);   _awave2[ihel] = aw2;
  }  
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
  Energy2 scale0(sqr(part.mass()));
  Energy2 scale1((momenta[0]+momenta[1]).m2());
  Energy2 scale2((momenta[2]+momenta[3]).m2());
  // for decays to quarks ensure boson is massive enough to
  // put quarks on constituent mass-shell
  if(scale1<sqr(outgoing[0]->constituentMass()+
  		outgoing[1]->constituentMass())) return 0.;
  if(scale2<sqr(outgoing[2]->constituentMass()+
  		outgoing[3]->constituentMass())) return 0.;
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
  // colour factors
  if(outgoing[0]->coloured()) output *= 3.;
  if(outgoing[2]->coloured()) output *= 3.;
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
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

void SMHiggsWWDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  for(unsigned int ix=0;ix<2;++ix) {
    _zmax[ix]=0.;
    _wmax[ix]=0.;
  }
  unsigned int ntest=_wdecays.size()+_zdecays.size();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      unsigned int iloc = ix<ntest ? 0 : 1;
      if(mode(ix)->outgoing()[0]->iCharge()==
	 -mode(ix)->outgoing()[1]->iCharge()) {
	_zmax[iloc]=max(mode(ix)->maxWeight(),_zmax[iloc]);
      }
      else {
	_wmax[iloc]=max(mode(ix)->maxWeight(),_wmax[iloc]);
      }
    }
  }
}
