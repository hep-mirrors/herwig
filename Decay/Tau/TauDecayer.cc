// -*- C++ -*-
//
// TauDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauDecayer class.
//
//  Author: Peter Richardson
//

#include "TauDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Decay/DecayVertex.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TauDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // set up the phase-space channels
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr channel;
  tPDVector extpart,ptemp;
  extpart.push_back(getParticleData(ParticleID::tauminus));
  extpart.push_back(getParticleData(ParticleID::nu_tau));
  Energy mtau(extpart[0]->mass());
  double maxweight;
  vector<double> channelwgts;
  int iq(0),ia(0);
  _modemap.clear();
  unsigned int ix,iy;
  bool done;
  vector<double>::iterator start,end;
  for(ix=0;ix<_current->numberOfModes();++ix) {
    // get the external particles for this mode
    extpart.resize(2);
    ptemp=_current->particles(-3,ix,iq,ia);
    for(iy=0;iy<ptemp.size();++iy) extpart.push_back(ptemp[iy]);
    // create the mode
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    // create the first piece of the channel
    channel = new_ptr(DecayPhaseSpaceChannel(mode));
    channel->addIntermediate(extpart[0],0,0.0,-1,1);
    done=_current->createMode(-3,ix,mode,2,1,channel,mtau);
    if(done) {
      // the maximum weight and the channel weights
      // the maximum
      maxweight = _wgtmax.size()>numberModes() ? _wgtmax[numberModes()] : 0;
      // the weights for the channel
      if(_wgtloc.size()>numberModes()&&
	 _wgtloc[numberModes()]+mode->numberChannels()<=_weights.size()) {
	start=_weights.begin()+_wgtloc[numberModes()];
	end  = start+mode->numberChannels();
	channelwgts=vector<double>(start,end);
      }
      else {
	channelwgts.resize(mode->numberChannels(),1./(mode->numberChannels()));
      }
      _modemap.push_back(ix);
      // special for the two body modes
      if(extpart.size()==3) {
	channelwgts.clear();
	mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      }
      addMode(mode,maxweight,channelwgts);
    }
  }
  _current->reset();
  _current->touch();
  _current->update();
}

void TauDecayer::doinitrun() {
  _current->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _weights.clear();_wgtloc.clear();_wgtmax.clear();
    unsigned int ix,iy;
    for(ix=0;ix<numberModes();++ix) {
      _wgtmax.push_back(mode(ix)->maxWeight());
      _wgtloc.push_back(_weights.size());
      for(iy=0;iy<mode(ix)->numberChannels();++iy) {
	_weights.push_back(mode(ix)->channelWeight(iy));
      }
    }
  }
}

bool TauDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  bool allowed(false);
  // find the neutrino 
  int idnu(0),idtemp,idin(parent->id());
  vector<int> idother;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)==16) idnu=idtemp; 
    else                idother.push_back(idtemp);
  }
  if((idnu==ParticleID::nu_tau    && idin==ParticleID::tauminus)||
     (idnu==ParticleID::nu_taubar && idin==ParticleID::tauplus )) {
    allowed=_current->accept(idother);
  }
  return allowed;
}


int TauDecayer::modeNumber(bool & cc,tcPDPtr parent, const tPDVector & children) const {
  int imode(-1);
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp;vector<int> idother;
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)!=16) idother.push_back(idtemp);
  }
  unsigned int itemp=_current->decayMode(idother);
  for(unsigned int ix=0;ix<_modemap.size();++ix) {
    if(_modemap[ix]==itemp) imode=ix;
  }
  // perform the decay
  cc=parent->id()==ParticleID::tauplus;
  return imode;
}


void TauDecayer::persistentOutput(PersistentOStream & os) const {
  os << _modemap << _current << _wgtloc 
     << _wgtmax << _weights << _polOpt << _tauMpol << _tauPpol;
}

void TauDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _modemap >> _current >> _wgtloc 
     >> _wgtmax >> _weights >> _polOpt >> _tauMpol >> _tauPpol;
}

ClassDescription<TauDecayer> TauDecayer::initTauDecayer;
// Definition of the static class description member.

void TauDecayer::Init() {

  static ClassDocumentation<TauDecayer> documentation
    ("The TauDecayer class is designed to use a weak current"
     " to perform the decay of the tau.");

  static Reference<TauDecayer,WeakDecayCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &TauDecayer::_current, false, false, true, false, false);

  static ParVector<TauDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &TauDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<TauDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &TauDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<TauDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &TauDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);

  static Switch<TauDecayer,bool> interfacePolarizationOption
    ("PolarizationOption",
     "Option of forcing the polarization of the tau leptons, N.B. you"
     " should only use this option for making distributions for"
     " comparision if you really know what you are doing.",
     &TauDecayer::_polOpt, false, false, false);
  static SwitchOption interfacePolarizationOptionDefault
    (interfacePolarizationOption,
     "Default",
     "Don't force the polarization use the full spin density matrices"
     " to get the right answer",
     false);
  static SwitchOption interfacePolarizationOptionForce
    (interfacePolarizationOption,
     "Force",
     "Force the polarizations",
     true);

  static Parameter<TauDecayer,double> interfaceTauMinusPolarization
    ("TauMinusPolarization",
     "The polarization of the tau-, left=-1, right=+1 if this is forced.",
     &TauDecayer::_tauMpol, 0.0, -1.0, 1.0,
     false, false, Interface::limited);


  static Parameter<TauDecayer,double> interfaceTauPlusPolarization
    ("TauPlusPolarization",
     "The polarization of the tau+, left=-1, right=+1 if this is forced.",
     &TauDecayer::_tauPpol, 0.0, -1.0, 1.0,
     false, false, Interface::limited);

}

// combine the currents to give the matrix element
double TauDecayer::me2(const int ichan,const Particle & inpart,
		       const ParticleVector & decay,
		       MEOption meopt) const {
  // map the mode to those in the current
  int mode(_modemap[imode()]);
  // get the particles for the hadronic current
  ParticleVector hadpart(decay.begin()+1,decay.end());
  Energy q;
  // extract info on the decaying particle
  if(meopt==Initialize) {
    // spin density matrix for the decaying particle
    _rho = RhoDMatrix(PDT::Spin1Half);
    if(inpart.id()==ParticleID::tauminus)
      SpinorWaveFunction   ::calculateWaveFunctions(_inspin,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inbar ,_rho,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
    if(_polOpt) {
      _rho(0,1) = _rho(1,0) = 0.;
      if(inpart.id()==ParticleID::tauminus) {
	_rho(0,0) = 0.5*(1.-_tauMpol);
	_rho(1,1) = 0.5*(1.+_tauMpol);
      }
      else {
	_rho(0,0) = 0.5*(1.+_tauPpol);
	_rho(1,1) = 0.5*(1.-_tauPpol);
      }
    }
    // work out the mapping for the hadron vector
    _constants = vector<unsigned int>(decay.size()+1);
    _ispin     = vector<PDT::Spin   >(decay.size());
    int itemp(1);
    unsigned int ix(decay.size());
    do {
      --ix;
      _ispin[ix]     = decay[ix]->data().iSpin();
      itemp         *= _ispin[ix];
      _constants[ix] = itemp;
    }
    while(ix>0);
    _constants[decay.size()] = 1;
    _constants[0           ] = _constants[1];
  }
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,_ispin)));  
  // connect the spininfo up if needed
  if(meopt==Terminate) {
    if(inpart.id()==ParticleID::tauminus) {
      SpinorWaveFunction   ::
	constructSpinInfo(_inspin,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::
	constructSpinInfo(_inbar,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(_inbar ,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::
	constructSpinInfo(_inspin,decay[0],outgoing,true);
    }
    _current->current(mode,ichan,q,hadpart,meopt);
    return 0.;
  }
  // calculate the spinors for the decay products
  if(inpart.id()==ParticleID::tauminus)
    SpinorBarWaveFunction::calculateWaveFunctions(_inbar ,decay[0],outgoing);
  else   
    SpinorWaveFunction   ::calculateWaveFunctions(_inspin,decay[0],outgoing);
  // calculate the hadron current
  vector<LorentzPolarizationVectorE> 
    hadron(_current->current(mode,ichan,q,hadpart,meopt));
  // prefactor
  double pre = sqr(pow(inpart.mass()/q,int(hadpart.size()-2)));
  // calculate the lepton current
  LorentzPolarizationVectorE lepton[2][2];
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      if(inpart.id()==15) 
	lepton[ix][iy]=2.*_inspin[ix].leftCurrent(_inbar[iy]); 
      else                
	lepton[iy][ix]=2.*_inspin[ix].leftCurrent(_inbar[iy]); 
    }
  }
  // compute the matrix element
  vector<unsigned int> ihel(decay.size()+1);
  for(unsigned int hhel=0;hhel<hadron.size();++hhel) {
    // map the index for the hadrons to a helicity state
    for(unsigned int ix=decay.size();ix>1;--ix) {
      ihel[ix]=(hhel%_constants[ix-1])/_constants[ix];
    }
    // loop over the helicities of the tau and neutrino and set up the matrix 
    // element
    for(ihel[1]=0;ihel[1]<2;++ihel[1]){
      for(ihel[0]=0;ihel[0]<2;++ihel[0]) {
	(*ME())(ihel)= lepton[ihel[0]][ihel[1]].dot(hadron[hhel])*
	  SM().fermiConstant();
      }
    }
  }
  // multiply by the CKM element
  int iq,ia;
  _current->decayModeInfo(mode,iq,ia);
  double ckm(1.);
  if(iq<=6) {
    if(iq%2==0) ckm = SM().CKM(iq/2-1,(abs(ia)-1)/2);
    else        ckm = SM().CKM(abs(ia)/2-1,(iq-1)/2);
  }
  return 0.5*pre*ckm*(ME()->contract(_rho)).real();
}
  
// output the setup information for the particle database
void TauDecayer::dataBaseOutput(ofstream & output,bool header) const {
  unsigned int ix;
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  for(ix=0;ix<_wgtloc.size();++ix) {
    output << "insert " << name() << ":WeightLocation " << ix << " " 
	   << _wgtloc[ix] << "\n";
  }
  for(ix=0;ix<_wgtmax.size();++ix) {
    output << "insert " << name() << ":MaximumWeight "  << ix << " " 
	   << _wgtmax[ix] << "\n";
  }
  for(ix=0;ix<_weights.size();++ix) {
    output << "insert " << name() << ":Weights "        << ix << " " 
	   << _weights[ix] << "\n";
  }
  _current->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":WeakCurrent " << _current->name() << " \n";
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
}
