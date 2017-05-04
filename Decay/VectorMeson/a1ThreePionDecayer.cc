// -*- C++ -*-
//
// a1ThreePionDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1ThreePionDecayer class.
//

#include "a1ThreePionDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

a1ThreePionDecayer::a1ThreePionDecayer() 
  : _rhomass(1,0.7761*GeV), _rhowidth(1,0.1445*GeV), _sigmamass(0.8*GeV),
    _sigmawidth(0.8*GeV), _psigma(ZERO), _mpi(ZERO), _mpi2(ZERO),
    _lambda2(1.2*GeV2), _a1mass2(1.23*1.23*GeV2),
    _zsigma(0.), _zmag(1.3998721), _zphase(0.43585036),
    _rhomag(1,1.), _rhophase(1,0.), _coupling(90.44), 
    _localparameters(true), _zerowgts(3), _onewgts(7), _twowgts(7),
    _threewgts(8) ,_zeromax(19.144), _onemax(7.83592), 
    _twomax(6.64804), _threemax(6.66296) {
  // weights for the different channels
  _threewgts[0] = 0.153071;
  _threewgts[1] = 0.165741;
  _threewgts[2] = 0.107509;
  _threewgts[3] = 0.10275 ;
  _threewgts[4] = 0.109738;
  _threewgts[5] = 0.11254 ;
  _threewgts[6] = 0.125344;
  _threewgts[7] = 0.123307;
  _onewgts[0] = 0.19616 ;
  _onewgts[1] = 0.191408;
  _onewgts[2] = 0.12137 ;
  _onewgts[3] = 0.115498;
  _onewgts[4] = 0.12729 ;
  _onewgts[5] = 0.127183;
  _onewgts[6] = 0.12109 ;
  _zerowgts[0] = 0.339108;
  _zerowgts[1] = 0.335601;
  _zerowgts[2] = 0.325291;
  _twowgts[0] = 0.188163;
  _twowgts[1] = 0.192479;
  _twowgts[2] = 0.121658;
  _twowgts[3] = 0.12135 ;
  _twowgts[4] = 0.127298;
  _twowgts[5] = 0.124835;
  _twowgts[6] = 0.124217;
  // generation of intermediates
  generateIntermediates(true);
}

void a1ThreePionDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    // get the weights for the different channels
    for(unsigned int ix=0;ix<_zerowgts.size();++ix)
      _zerowgts[ix]=mode(0)->channelWeight(ix);
    for(unsigned int ix=0;ix<_onewgts.size();++ix)
      _onewgts[ix]=mode(1)->channelWeight(ix);
    for(unsigned int ix=0;ix<_twowgts.size();++ix)
      _twowgts[ix]=mode(2)->channelWeight(ix);
    for(unsigned int ix=0;ix<_threewgts.size();++ix)
      _threewgts[ix]=mode(3)->channelWeight(ix);
    // get the maximum weight
    _zeromax  = mode(0)->maxWeight();
    _onemax   = mode(1)->maxWeight();
    _twomax   = mode(2)->maxWeight();
    _threemax = mode(3)->maxWeight();
  }
}

void a1ThreePionDecayer::doinit() {
  DecayIntegrator::doinit();
  // particles we need for the external state
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  // possible intermediate particles
  // the different rho resonances
  tPDPtr rhop[3] = {getParticleData(213),getParticleData(100213),
		    getParticleData(30213)};
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),
		    getParticleData(30113)};
  tPDPtr rhom[3] = {getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // the sigma
  tPDPtr sigma = getParticleData(9000221);
  // set up the phase space integration
  tPDVector extpart(4);
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr newchannel;
  // decay mode a_0 -> pi0 pi0 pi0
  extpart[0]=a10;
  extpart[1]=pi0;
  extpart[2]=pi0;
  extpart[3]=pi0;
  if(sigma) {
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,1);
    newchannel->addIntermediate(sigma,0,0.0,2,3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,2);
    newchannel->addIntermediate(sigma,0,0.0,1,3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,3);
    newchannel->addIntermediate(sigma,0,0.0,1,2);
    mode->addChannel(newchannel);
  }
  if(_zerowgts.size()!=mode->numberChannels()) 
    _zerowgts=vector<double>(mode->numberChannels(),1./mode->numberChannels());
  addMode(mode,_zeromax,_zerowgts);
  // decay mode a_1+ -> pi+ pi0 pi0
  extpart[0]=a1p;
  extpart[1]=pi0;
  extpart[2]=pi0;
  extpart[3]=pip;
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho+ channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,1);
    newchannel->addIntermediate(rhop[ix],0,0.0,2,3);
    mode->addChannel(newchannel);
    // second rho+ channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,2);
    newchannel->addIntermediate(rhop[ix],0,0.0,1,3);
    mode->addChannel(newchannel);
  }
  // the sigma channel
  if(sigma) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,3);
    newchannel->addIntermediate(sigma,0,0.0,1,2);
    mode->addChannel(newchannel);
  }
  if(_onewgts.size()!=mode->numberChannels()) 
    _onewgts=vector<double>(mode->numberChannels(),1./mode->numberChannels());
  addMode(mode,_onemax,_onewgts);
  // decay mode a_1 -> pi+ pi- pi0
  extpart[0]=a10;
  extpart[1]=pip;
  extpart[2]=pim;
  extpart[3]=pi0;
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,1);
    newchannel->addIntermediate(rhom[ix],0,0.0,2,3);
    mode->addChannel(newchannel);
    // second channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,2);
    newchannel->addIntermediate(rhop[ix],0,0.0,1,3);
    mode->addChannel(newchannel);
  }
  // sigma channel
  if(sigma) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,3);
    newchannel->addIntermediate(sigma,0,0.0,1,2);
    mode->addChannel(newchannel);
  }
  if(_twowgts.size()!=mode->numberChannels()) 
    _twowgts=vector<double>(mode->numberChannels(),1./mode->numberChannels());
  addMode(mode,_twomax,_twowgts);
  // decay mode a_1+ -> pi+ pi+ pi-
  extpart[0]=a1p;
  extpart[1]=pip;
  extpart[2]=pip;
  extpart[3]=pim;
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  for(unsigned int ix=0;ix<3;++ix) {
    // the neutral rho channels
    if(!rho0[ix]) continue;
    // first channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,1);
    newchannel->addIntermediate(rho0[ix],0,0.0,2,3);
    mode->addChannel(newchannel);
    // interchanged channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,2);
    newchannel->addIntermediate(rho0[ix],0,0.0,1,3);
    mode->addChannel(newchannel);
  }
  // the sigma channels
  if(sigma) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,1);
    newchannel->addIntermediate(sigma,0,0.0,2,3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,2);
    newchannel->addIntermediate(sigma,0,0.0,1,3);
    mode->addChannel(newchannel);
  }
  if(_threewgts.size()!=mode->numberChannels()) 
    _threewgts=vector<double>(mode->numberChannels(),1./mode->numberChannels());
  addMode(mode,_threemax,_threewgts);
  // set up the parameters 
  _mpi=getParticleData(ParticleID::piplus)->mass();
  _mpi2=sqr(_mpi);
  if(_localparameters) {
    if(_rhomass.size()<_rhocoupling.size()) {
      unsigned int itemp=_rhomass.size();
      _rhomass.resize(_rhocoupling.size());
      _rhowidth.resize(_rhocoupling.size());
      for(unsigned int ix=itemp;ix<_rhocoupling.size();++ix) {
	_rhomass[ix]=rhop[ix]->mass();
	_rhowidth[ix]=rhop[ix]->width();
      }
      // reset the intermediates in the phase space integration if needed
      resetIntermediate(sigma,_sigmamass,_sigmawidth);
      for(unsigned int iy=0;iy<_rhocoupling.size();++iy) {
	resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
	resetIntermediate(rhop[iy],_rhomass[iy],_rhowidth[iy]);
	resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
      }
    }
  }
  else {
    _a1mass2=sqr(getParticleData(ParticleID::a_1plus)->mass());
    if(sigma) {
      _sigmamass=sigma->mass();
      _sigmawidth=sigma->width();
    }
    _rhomass.resize(_rhocoupling.size());
    _rhowidth.resize(_rhocoupling.size());
    for(unsigned int ix=0;ix<_rhocoupling.size();++ix) {
      _rhomass[ix]=rhop[ix]->mass();
      _rhowidth[ix]=rhop[ix]->width();
    }
  }
  // parameters for the resonances
  // for the sigma
  _psigma=Kinematics::pstarTwoBodyDecay(_sigmamass,_mpi,_mpi);
  // for the rho
  _prho.resize(_rhomass.size());_hm2.resize(_rhomass.size());
  _dhdq2m2.resize(_rhomass.size());_rhoD.resize(_rhomass.size());
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    _prho[ix]    = Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpi,_mpi);
    _hm2[ix]     = hFunction(_rhomass[ix]);
    _dhdq2m2[ix] = dhdq2Parameter(ix);
    _rhoD[ix]    = DParameter(ix);
  }
  // convert the magnitude and phase of z into a phase
  _zsigma = _zmag*Complex(cos(_zphase),sin(_zphase));
  // convert rho couplings
  for(unsigned int ix=0;ix<_rhomag.size();++ix) {
    _rhocoupling.push_back(_rhomag[ix]*Complex(cos(_rhophase[ix]),sin(_rhophase[ix])));
  }
}
  
int a1ThreePionDecayer::modeNumber(bool & cc,tcPDPtr parent,
				       const tPDVector & children) const {
  if(children.size()!=3) return -1;
  int id(parent->id());
  // check the pions
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,npi0(0),npiplus(0),npiminus(0);
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
  }
  int imode(-1);
  // a_1+ decay modes
  if(id==ParticleID::a_1plus) {
    cc=false;
    if(npiplus==1&&npi0==2)          imode=1;
    else if(npiplus==2&&npiminus==1) imode=3;
  }
  // a_1- modes
  else if(id==ParticleID::a_1minus) {
    cc=true;
    if(npiminus==1&&npi0==2)         imode=1;
    else if(npiminus==2&&npiplus==1) imode=3;
  }
  // a_0 modes
  else if(id==ParticleID::a_10) {
    cc=false;
    if(npiminus==1&&npiplus==1&&npi0==1) imode=2;
    else if(npi0==3)                     imode=0;
  }
  return imode;
}
  
void a1ThreePionDecayer::persistentOutput(PersistentOStream & os) const {
   os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV) << ounit(_prho,GeV) 
      << ounit(_hm2,GeV2) << ounit(_rhoD,GeV2) << _dhdq2m2 <<  ounit(_sigmamass,GeV)
      << ounit(_sigmawidth,GeV) << ounit(_psigma,GeV) << ounit(_mpi,GeV)
      << ounit(_mpi2,GeV2) << ounit(_lambda2,GeV2) << ounit(_a1mass2,GeV2) << _zsigma  
      << _rhocoupling << _coupling << _localparameters << _zerowgts << _onewgts 
      << _twowgts << _threewgts << _zeromax << _zmag << _zphase
      << _onemax << _twomax << _threemax << _coupling << _rhomag << _rhophase;
}
  
void a1ThreePionDecayer::persistentInput(PersistentIStream & is, int) {
   is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) >> iunit(_prho,GeV) 
      >> iunit(_hm2,GeV2) >> iunit(_rhoD,GeV2) >> _dhdq2m2 >>  iunit(_sigmamass,GeV)
      >> iunit(_sigmawidth,GeV) >> iunit(_psigma,GeV) >> iunit(_mpi,GeV) 
      >> iunit(_mpi2,GeV2) >> iunit(_lambda2,GeV2) >> iunit(_a1mass2,GeV2) >> _zsigma
      >> _rhocoupling >> _coupling >> _localparameters >> _zerowgts >> _onewgts 
      >> _twowgts >> _threewgts >> _zeromax >> _zmag >> _zphase
      >> _onemax >> _twomax >> _threemax >> _coupling >> _rhomag >> _rhophase;
}

ClassDescription<a1ThreePionDecayer> a1ThreePionDecayer::inita1ThreePionDecayer;
// Definition of the static class description member.
  
void a1ThreePionDecayer::Init() {
    
  static ClassDocumentation<a1ThreePionDecayer> documentation
    ("The a1ThreePionDecayer class is designed to decay the a_1 "
     "resonance to three pions using a model based on that used in the modelling "
     "of tau->4 pions.",
     "The decay of the $a_1$ resonance to three pions uses a model based on"
     "tau to four pions, \\cite{Bondar:2002mw}.",
     "%\\cite{Bondar:2002mw}\n"
     "\\bibitem{Bondar:2002mw}\n"
     "  A.~E.~Bondar, S.~I.~Eidelman, A.~I.~Milstein, T.~Pierzchala, N.~I.~Root, Z.~Was and M.~Worek,\n"
     "   ``Novosibirsk hadronic currents for tau --> 4pi channels of tau decay\n"
     "  %library TAUOLA,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 146}, 139 (2002)\n"
     "  [arXiv:hep-ph/0201149].\n"
     "  %%CITATION = CPHCB,146,139;%%\n"
     );

  static Switch<a1ThreePionDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1ThreePionDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use the local values",
     true);
  static SwitchOption interfaceLocalParametersDefault
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particleData objects",
     false);

  static Parameter<a1ThreePionDecayer,double> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1ThreePionDecayer::_coupling, 90.44, 0.0, 1000.0,
     false, false, true);

  static ParVector<a1ThreePionDecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1ThreePionDecayer::_rhomass,
     GeV, 0, ZERO, ZERO, 10000*GeV, false, false, true);

  static ParVector<a1ThreePionDecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1ThreePionDecayer::_rhowidth,
     GeV, 0, ZERO, ZERO, 10000*GeV, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "The magnitude of the rho couplings",
     &a1ThreePionDecayer::_rhomag, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static ParVector<a1ThreePionDecayer,double> interfaceRhoPhase
    ("RhoPhase",
     "The phase of the rho coupling",
     &a1ThreePionDecayer::_rhophase, -1, 0., 0.0, 2.*Constants::pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionDecayer,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of the mass scale squared to use in the form-factor",
     &a1ThreePionDecayer::_lambda2, GeV2, 1.2*GeV2, 0.0001*GeV2, 10.0*GeV2,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy2> interfacea1mass2
    ("a1mass2",
     "The local value of the square of the a_1 mass",
     &a1ThreePionDecayer::_a1mass2, GeV2, 1.5129*GeV2, 0.5*GeV2, 10.0*GeV2,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy> interfaceSigmaMass
    ("SigmaMass",
     "The local value of the sigma mass",
     &a1ThreePionDecayer::_sigmamass, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy> interfaceSigmaWidth
    ("SigmaWidth",
     "The local value of the sigma width",
     &a1ThreePionDecayer::_sigmawidth, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacezerowgts
    ("AllNeutralWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi0pi0pi0",
     &a1ThreePionDecayer::_zerowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1ThreePionDecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1ThreePionDecayer::_twowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1ThreePionDecayer::_threewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceZeroMax
    ("ZeroMax",
     "The maximum weight for the integration fo the channel a_1^0->pi0pi0pi0",
     &a1ThreePionDecayer::_zeromax, 107.793, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1ThreePionDecayer::_onemax, 1088.96, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1ThreePionDecayer::_twomax, 1750.73, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1ThreePionDecayer::_threemax, 739.334, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceSigmaMagnitude
    ("SigmaMagnitude",
     "magnitude of the relative sigma coupling",
     &a1ThreePionDecayer::_zmag, 1.3998721, 0.0, 10.0e20,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceSigmaPhase
    ("SigmaPhase",
     "phase of the relative sigma coupling",
     &a1ThreePionDecayer::_zphase, 0.43585036, 0.0, Constants::twopi,
     false, false, true);
}

double a1ThreePionDecayer::me2(const int ichan,
			       const Particle & inpart,
			       const ParticleVector & decay,
			       MEOption meopt) const {
  useMe();
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
						const_ptr_cast<tPPtr>(&inpart),
						incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    // set up the spin information for the decay products
    for(unsigned int ix=0;ix<3;++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return 0.;
  }
  // momentum of the incoming particle
  Lorentz5Momentum Q=inpart.momentum();
  // momenta of the intermediates
  Energy2 s1=(decay[1]->momentum()+decay[2]->momentum()).m2();
  Energy2 s2=(decay[0]->momentum()+decay[2]->momentum()).m2();
  Energy2 s3=(decay[0]->momentum()+decay[1]->momentum()).m2();
  Energy2 dot01=Q*decay[0]->momentum();
  Energy2 dot02=Q*decay[1]->momentum();
  Energy2 dot03=Q*decay[2]->momentum();
  // vector for the output
  LorentzVector<complex<Energy3> > output;
  // a_10 -> pi0pi0pi0
  if(imode()==0) {
    //the breit-wigners
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    Complex sig3=sigmaBreitWigner(s3);
    // compute the vector
    LorentzPolarizationVectorE tmpoutput;
    if(ichan<0) {
      tmpoutput= sig1*(decay[0]->momentum())+sig2*(decay[1]->momentum())
	+sig3*(decay[2]->momentum());
    }
    else if(ichan==0) tmpoutput=sig1*(decay[0]->momentum());
    else if(ichan==1) tmpoutput=sig2*(decay[1]->momentum());
    else if(ichan==2) tmpoutput=sig3*(decay[2]->momentum());
    // the coupling z and identical particle factor
    output = tmpoutput * _zsigma* 1./sqrt(6.) *Q.mass2();
  }
  // a_1+ -> pi0pi0pi+
  else if(imode()==1) {
    // scalar propagator
    Complex sig1 = sigmaBreitWigner(s3);
    // sigma terms
    if(ichan<0||ichan==6) 
      output = _zsigma*Q.mass2()*sig1*decay[2]->momentum();
    // the rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1=_rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2=_rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output +=rho1*(dot03*(decay[1]->momentum())-
		       dot02*(decay[2]->momentum()));
      }
      if(ichan<0||ichan==2*ix+1){
	output +=rho2*(dot03*(decay[0]->momentum())-
		       dot01*(decay[2]->momentum()));
      }
    }
    // the identical particle factor
    output *= 1./sqrt(2.);
  }
  // a_10->pi+pi-pi0
  else if(imode()==2) {
    // the sigma terms
    Complex sig1=sigmaBreitWigner(s3);
    if(ichan<0||ichan==6)
      output = _zsigma*Q.mass2()*sig1*decay[2]->momentum();
    // rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1=_rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2=_rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output+=rho1*(dot03*(decay[1]->momentum())
		      -dot02*(decay[2]->momentum()));
      }
      if(ichan<0||ichan==2*ix+1) {
	output+=rho2*(dot03*(decay[0]->momentum())
		      -dot01*(decay[2]->momentum()));
      }
    }
  }
  // a1+ -> pi+pi+pi-
  else if(imode()==3) {
    // the scalar propagators 
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    // sigma terms
    LorentzPolarizationVectorE tmpoutput;
    if(ichan<0||ichan==6) tmpoutput+=sig1*(decay[0]->momentum());
    if(ichan<0||ichan==7) tmpoutput+=sig2*(decay[1]->momentum());
    output = tmpoutput * _zsigma * Q.mass2();
    // rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1 = _rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2 = _rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output-=rho1*( dot03*(decay[1]->momentum())-
		       dot02*(decay[2]->momentum()));
      }
      if(ichan<0||ichan==2*ix+1) {
	output-=rho2*( dot03*(decay[0]->momentum())-
		       dot01*(decay[2]->momentum()));
      }
    }
    // the identical particle factor
    output *= 1./sqrt(2.);
  }
  // form-factor
  LorentzPolarizationVector outputFinal 
    = output * a1FormFactor(Q.mass2())*_coupling/(Q.mass()*sqr(_rhomass[0]));
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix)
    (*ME())(ix,0,0,0)=outputFinal.dot(_vectors[ix]);
  // return the answer
  double out = ME()->contract(_rho).real();
  // test of the answer
//   double test = threeBodyMatrixElement(imode(),sqr(inpart.mass()),s3,s2,s1,
// 				       decay[0]->mass(),decay[1]->mass(), 
// 				       decay[2]->mass());
//   if(ichan<0) cerr << "testing matrix element " << inpart.PDGName() << " -> "
// 		   << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
// 		   << decay[2]->PDGName() << " " << out << " " << test << " " 
// 		   << (out-test)/(out+test) << "\n"; 
  return out;
}

// matrix element for the running a_1 width
double a1ThreePionDecayer::
threeBodyMatrixElement(const int iopt,const Energy2 q2, const Energy2 s3,
		       const Energy2 s2,const Energy2 s1, const Energy m1, 
		       const Energy m2 ,const Energy m3) const {
  Energy6 meout(0.*pow<3,1>(GeV2));
  Energy2 m12(sqr(m1)),m22(sqr(m2)),m32(sqr(m3));
  Energy2 dot01(q2-s1+m12),dot02(q2-s2+m22),dot03(q2-s3+m32),
    dot12(s3-m12-m22),dot13(s2-m12-m32),dot23(s1-m22-m32);
  if(iopt==0) {
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    Complex sig3=sigmaBreitWigner(s3);
    Energy2 metemp = 
      real(0.25*sig1*conj(sig1)*lambda(q2,s1,m12)/q2+
	   0.25*sig2*conj(sig2)*lambda(q2,s2,m22)/q2+
	   0.25*sig3*conj(sig3)*lambda(q2,s3,m32)/q2+
 	   sig1*conj(sig2)*(-dot12+0.5*dot01*dot02/q2)+
 	   sig1*conj(sig3)*(-dot13+0.5*dot01*dot03/q2)+
 	   sig2*conj(sig3)*(-dot23+0.5*dot02*dot03/q2));
    meout = metemp*real(_zsigma*conj(_zsigma))/6.*sqr(q2);
  }
  else if(iopt==1||iopt==2) {
    // the sigma terms
    Complex sig=sigmaBreitWigner(s3);
    Complex rho1,rho2;
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      rho1 += _rhocoupling[ix]*rhoBreitWigner(s1,ix);
      rho2 += _rhocoupling[ix]*rhoBreitWigner(s2,ix);
    }
    meout =
      0.25*lambda(q2,m32,s3)*q2*norm(_zsigma*sig)+
      0.25*norm(rho1)*(dot23*dot02*dot03-m32*sqr(dot02)-m22*sqr(dot03))+
      0.25*norm(rho2)*(dot13*dot01*dot03-m32*sqr(dot01)-m12*sqr(dot03))-
      0.5*real(_zsigma*sig*conj(rho1))*q2*(dot03*dot23-2.*m32*dot02)-
      0.5*real(_zsigma*sig*conj(rho2))*q2*(dot03*dot13-2.*m32*dot01)-
      0.25*real(rho1*conj(rho2))*(sqr(dot03)*dot12-dot03*dot02*dot13
				  -dot03*dot01*dot23+2.*m32*dot02*dot01);
    if(iopt==1) meout *= 0.5;
  }
  else if(iopt==3) {
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    Complex rho1,rho2;
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      rho1 += _rhocoupling[ix]*rhoBreitWigner(s1,ix);
      rho2 += _rhocoupling[ix]*rhoBreitWigner(s2,ix);
    }
    meout =
      0.25*lambda(q2,m12,s1)*q2*norm(_zsigma*sig1)+
      0.25*lambda(q2,m22,s2)*q2*norm(_zsigma*sig2)+
      0.25*norm(rho1)*(dot23*dot02*dot03-m32*sqr(dot02)-m22*sqr(dot03))+
      0.25*norm(rho2)*(dot13*dot01*dot03-m32*sqr(dot01)-m12*sqr(dot03))-
      0.25*real(rho1*conj(rho2))*(sqr(dot03)*dot12-dot03*dot02*dot13
				  -dot03*dot01*dot23+2.*m32*dot02*dot01)-
      real(_zsigma*sig1*conj(_zsigma*sig2))*q2*(q2*dot12-0.5*dot02*dot01)+
      0.5*real(_zsigma*sig1*conj(rho1))*q2*(dot03*dot12-dot02*dot13)+
      0.5*real(_zsigma*sig2*conj(rho1))*q2*(2.*dot03*m22-dot02*dot23)+
      0.5*real(_zsigma*sig1*conj(rho2))*q2*(2.*dot03*m12-dot01*dot13)+
      0.5*real(_zsigma*sig2*conj(rho2))*q2*(dot03*dot12-dot01*dot23);
    meout *= 0.5;
  }
  return meout*a1FormFactor(q2)*sqr(_coupling/sqr(_rhomass[0]))/q2/3.;
}

WidthCalculatorBasePtr 
a1ThreePionDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit) {
    if(abs((**pit).id())==ParticleID::piplus) ++ncharged;
  }
  // integrator to perform the integral
  vector<double> inweights(2,0.5);
  vector<int> intype;intype.push_back(2);intype.push_back(3);
  vector<Energy> inmass(2,_rhomass[0]),inwidth(2,_rhowidth[0]);
  vector<double> inpow(2,0.0);
  Energy m[3];
  Energy mpi0=getParticleData(ParticleID::pi0)->mass();
  Energy mpic=getParticleData(ParticleID::piplus)->mass();
  m[0] = ncharged<2 ? mpi0 : mpic;
  m[1] = m[0];
  m[2] = (ncharged==0||ncharged==2) ? mpi0 : mpic;
  return new_ptr(ThreeBodyAllOnCalculator<a1ThreePionDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,ncharged,m[0],m[1],m[2]));
}

// output the setup information for the particle database
void a1ThreePionDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Coupling " << _coupling     << "\n";
  output << "newdef " << name() << ":Lambda2 "  << _lambda2/GeV2 << "\n";
  output << "newdef " << name() << ":a1mass2 "  << _a1mass2/GeV2 << "\n";
  output << "newdef " << name() << ":SigmaMass "  << _sigmamass/GeV  << "\n";
  output << "newdef " << name() << ":SigmaWidth " << _sigmawidth/GeV << "\n";
  output << "newdef " << name() << ":SigmaMagnitude " << _zmag << "\n";
  output << "newdef " << name() << ":SigmaPhase " << _zphase << "\n";
  for(unsigned int ix=0;ix<_rhomag.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoMagnitude " << ix << " " 
		    << _rhomag[ix] << "\n";
    else     output << "insert " << name() << ":RhoMagnitude " << ix << " " 
		    << _rhomag[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhophase.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoPhase " << ix << " " 
		    << _rhophase[ix] << "\n";
    else     output << "insert " << name() << ":RhoPhase " << ix << " " 
		    << _rhophase[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
  }
  // integration weights for the different channels
  for(unsigned int ix=0;ix<_zerowgts.size();++ix) {
    output << "newdef " << name() << ":AllNeutralWeights " 
	   << ix << " " << _zerowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_onewgts.size();++ix) {
    output << "newdef " << name() << ":OneChargedWeights " 
	   << ix << " " << _onewgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_twowgts.size();++ix) {
    output << "newdef " << name() << ":TwoChargedWeights " 
	   << ix << " " << _twowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_threewgts.size();++ix) {
    output << "newdef " << name() << ":ThreeChargedWeights " 
	   << ix << " " << _threewgts[ix] << "\n";
  }
  output << "newdef " << name() << ":ZeroMax "  << _zeromax  << "\n";
  output << "newdef " << name() << ":OneMax "   << _onemax   << "\n";
  output << "newdef " << name() << ":TwoMax "   << _twomax   << "\n";
  output << "newdef " << name() << ":ThreeMax " << _threemax << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
