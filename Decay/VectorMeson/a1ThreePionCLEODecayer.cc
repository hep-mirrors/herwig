// -*- C++ -*-
//
// a1ThreePionCLEODecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1ThreePionCLEODecayer class.
//

#include "a1ThreePionCLEODecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using Constants::pi;
 
a1ThreePionCLEODecayer::a1ThreePionCLEODecayer() 
  : _rhomass(2), _rhowidth(2), _f2mass(1.275*GeV), _f2width(0.185*GeV), 
    _pf2cc(ZERO), _pf200(ZERO), _f0mass(1.186*GeV), _f0width(0.350*GeV),
    _pf0cc(ZERO), _pf000(ZERO), _sigmamass(0.860*GeV), _sigmawidth(0.880*GeV), 
    _psigmacc(ZERO), _psigma00(ZERO), _mpi0(ZERO), _mpic(ZERO),
    _coupling(45.57/GeV), _rhomagP(2), _rhophaseP(2), _rhomagD(2),
    _rhophaseD(2), _f2mag(0.71/GeV2), _f2phase(0.56*pi), _f2coup(ZERO,ZERO),
    _f0mag(0.77), _f0phase(-0.54*pi), _f0coup(0.,0.), _sigmamag(2.10),
    _sigmaphase(0.23*pi), _sigmacoup(0.,0.), _localparameters(true),
    _zerowgts(9), _onewgts(9), _twowgts(9), _threewgts(12), _zeromax(13.0704),
    _onemax(6.91104), _twomax(6.94654), _threemax(6.40086) {
  // rho masses and widths
  _rhomass[0] = 0.7743*GeV; _rhowidth[0] = 0.1491*GeV;
  _rhomass[1] = 1.370 *GeV; _rhowidth[1] = 0.386 *GeV;
  // p-wave rho and rho prime
  _rhomagP[0] = 1.  ; _rhophaseP[0] = 0.;
  _rhomagP[1] = 0.12; _rhophaseP[1] = 0.99*pi;
  // d-wave rho and rho prime
  _rhomagD[0] = 0.37/GeV2; _rhophaseD[0] = -0.15*pi;
  _rhomagD[1] = 0.87/GeV2; _rhophaseD[1] =  0.53*pi;
  // set up the integration channels
  _zerowgts[0]  = 0.132162;_zerowgts[1]  = 0.116638;_zerowgts[2]  = 0.121088;
  _zerowgts[3]  = 0.10656 ;_zerowgts[4]  = 0.102577;_zerowgts[5]  = 0.101169;
  _zerowgts[6]  = 0.104587;_zerowgts[7]  = 0.104663;_zerowgts[8]  = 0.110557;
  _onewgts[0]   = 0.177017;_onewgts[1]   = 0.176011;_onewgts[2]   = 0.110129;
  _onewgts[3]   = 0.108023;_onewgts[4]   = 0.110553;_onewgts[5]   = 0.109976;
  _onewgts[6]   = 0.088634;_onewgts[7]   = 0.059104;_onewgts[8]   = 0.060553;
  _twowgts[0]   = 0.173357;_twowgts[1]   = 0.172283;_twowgts[2]   = 0.116031;
  _twowgts[3]   = 0.114642;_twowgts[4]   = 0.109058;_twowgts[5]   = 0.114073;
  _twowgts[6]   = 0.080946;_twowgts[7]   = 0.060135;_twowgts[8]   = 0.059477;
  _threewgts[0] = 0.125022;_threewgts[1] = 0.129911;_threewgts[2] = 0.074165;
  _threewgts[3] = 0.075813;_threewgts[4 ]= 0.071154;_threewgts[5 ]= 0.077730;
  _threewgts[6] = 0.082255;_threewgts[7 ]= 0.086761;_threewgts[8 ]= 0.067106;
  _threewgts[9] = 0.070171;_threewgts[10]= 0.070146;_threewgts[11]= 0.069767;
  // generation of intermediates
  generateIntermediates(true);
}
  
void a1ThreePionCLEODecayer::doinit() {
  DecayIntegrator::doinit();
  // pointers to the particles we need as external particles
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
  // the f_2
  tPDPtr f2=getParticleData(225);
  // the f_0
  tPDPtr f0=getParticleData(10221);
  // set up the integration channels
  tPDVector extpart(4);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  // decay mode a_10 -> pi0 pi0 pi0
  extpart[0]=a10;
  extpart[1]=pi0;
  extpart[2]=pi0;
  extpart[3]=pi0;
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  // there are six sigma channels
  tPDPtr temp;
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==0)      temp = sigma;
    else if(ix==1) temp = f2;
    else if(ix==2) temp = f0;
    if(temp) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(a10,0,0.0,-1,1);
      newchannel->addIntermediate(temp,0,0.0,2,3);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(a10,0,0.0,-1,2);
      newchannel->addIntermediate(temp,0,0.0,1,3);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(a10,0,0.0,-1,3);
      newchannel->addIntermediate(temp,0,0.0,1,2);
      mode->addChannel(newchannel);
    }
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
  //  the f_2  channel
  if(f2) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,3);
    newchannel->addIntermediate(f2,0,0.0,1,2);
    mode->addChannel(newchannel);
  }
  // the f_0 channel
  if(f0) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,3);
    newchannel->addIntermediate(f0,0,0.0,1,2);
    mode->addChannel(newchannel);
  }
  if(_onewgts.size()!=mode->numberChannels()) 
    _onewgts=vector<double>(mode->numberChannels(),1./mode->numberChannels());
  addMode(mode,_onemax,_onewgts);
  // decay mode a_10 -> pi+ pi- pi0
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
  // f_2 channel
  if(f2) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,3);
    newchannel->addIntermediate(f2,0,0.0,1,2);
    mode->addChannel(newchannel);
  }
  // f_0 channel
  if(f0) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,3);
    newchannel->addIntermediate(f0,0,0.0,1,2);
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
  // the f_2 channels
  if(f2) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,1);
    newchannel->addIntermediate(f2,0,0.0,2,3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,2);
    newchannel->addIntermediate(f2,0,0.0,1,3);
    mode->addChannel(newchannel);
  }
  // the f_0 channel
  if(f0) {
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,1);
    newchannel->addIntermediate(f0,0,0.0,2,3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,2);
    newchannel->addIntermediate(f0,0,0.0,1,3);
    mode->addChannel(newchannel);
  }
  if(_threewgts.size()!=mode->numberChannels()) 
    _threewgts=vector<double>(mode->numberChannels(),1./mode->numberChannels());
  addMode(mode,_threemax,_threewgts);
  // if using local parameters set the values in the phase space channels
  if(_localparameters) {
    for(unsigned int iy=0;iy<_rhomass.size();++iy) {
      resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
      resetIntermediate(rhop[iy],_rhomass[iy],_rhowidth[iy]);
      resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
    }
    resetIntermediate(sigma,_sigmamass,_sigmawidth);
    resetIntermediate(f2,_f2mass,_f2width);
    resetIntermediate(f0,_f0mass,_f0width);
    // make sure the rho array has enough masses
    if(_rhomass.size()<3) {
      for(unsigned int ix=_rhomass.size();ix<3;++ix) {
	_rhomass.push_back(rhop[ix]->mass());
	_rhowidth.push_back(rhop[ix]->width());
      }
    }
  }
  // set the local variables if needed
  else {
    // masses and widths for the particles
    _rhomass.resize(3);_rhowidth.resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      _rhomass[ix]=rhop[ix]->mass();
      _rhowidth[ix]=rhop[ix]->width();
    }
    if(f2) {
      _f2mass=f2->mass();
      _f2width=f2->width();
    }
    if(f0) {
      _f0mass=f0->mass();
      _f0width=f0->width();
    }
    if(sigma) {
      _sigmamass=sigma->mass();
      _sigmawidth=sigma->width();
    }
  }
  // parameters for the breit-wigners
  _mpic=pip->mass();
  _mpi0=pi0->mass();
  // momenta of the decay products for on-shell particles
  _psigmacc = Kinematics::pstarTwoBodyDecay(_sigmamass,_mpic,_mpic);
  _psigma00 = Kinematics::pstarTwoBodyDecay(_sigmamass,_mpi0,_mpi0);
  _pf2cc    = Kinematics::pstarTwoBodyDecay(_f2mass   ,_mpic,_mpic);
  _pf200    = Kinematics::pstarTwoBodyDecay(_f2mass   ,_mpi0,_mpi0);
  _pf0cc    = Kinematics::pstarTwoBodyDecay(_f0mass   ,_mpic,_mpic);
  _pf000    = Kinematics::pstarTwoBodyDecay(_f0mass   ,_mpi0,_mpi0); 
  _prhocc.resize(3);_prhoc0.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    _prhocc[ix] = Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpic,_mpic);
    _prhoc0[ix] = Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpic,_mpi0);
  }
  // couplings for the different modes
  Complex ii(0.,1.);
  _rhocoupP.resize(_rhomagP.size());
  for(unsigned int ix=0;ix<_rhomagP.size();++ix)
    _rhocoupP[ix]=_rhomagP[ix]*(cos(_rhophaseP[ix])+ii*sin(_rhophaseP[ix]));
  _rhocoupD.resize(_rhomagD.size());
  for(unsigned int ix=0;ix<_rhomagD.size();++ix)
    _rhocoupD[ix]=_rhomagD[ix]*(cos(_rhophaseD[ix])+ii*sin(_rhophaseD[ix]));
  _f0coup=_f0mag*(cos(_f0phase)+ii*sin(_f0phase));
  _f2coup=_f2mag*(cos(_f2phase)+ii*sin(_f2phase));
  _sigmacoup=_sigmamag*(cos(_sigmaphase)+ii*sin(_sigmaphase));
}


inline void a1ThreePionCLEODecayer::doinitrun() {
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

int a1ThreePionCLEODecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void a1ThreePionCLEODecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV) 
     << ounit(_prhocc,GeV) << ounit(_prhoc0,GeV) 
     << ounit(_f2mass,GeV) << ounit(_f2width,GeV) 
     << ounit(_pf2cc,GeV) 
     << ounit(_pf200,GeV) << ounit(_f0mass,GeV) 
     << ounit(_f0width,GeV) << ounit(_pf0cc,GeV) << ounit(_pf000,GeV) 
     << ounit(_sigmamass,GeV) << ounit(_sigmawidth,GeV)
     << ounit(_psigmacc,GeV) << ounit(_psigma00,GeV) 
     << ounit(_mpi0,GeV) << ounit(_mpic,GeV) 
     << ounit(_coupling,1/GeV) << _rhomagP << _rhophaseP 
     << _rhocoupP << ounit(_rhomagD ,1/GeV2)<< _rhophaseD 
     << ounit(_rhocoupD,1/GeV2) << ounit(_f2mag,1/GeV2)
     << _f2phase << ounit(_f2coup,1/GeV2)
     << _f0mag << _f0phase << _f0coup 
     << _sigmamag << _sigmaphase << _sigmacoup 
     << _localparameters 
     << _zerowgts << _onewgts << _twowgts << _threewgts 
     << _zeromax << _onemax << _twomax << _threemax;
}
  
void a1ThreePionCLEODecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) 
     >> iunit(_prhocc,GeV) >> iunit(_prhoc0,GeV) 
     >> iunit(_f2mass,GeV) >> iunit(_f2width,GeV) 
     >> iunit(_pf2cc,GeV)
     >> iunit(_pf200,GeV) >> iunit(_f0mass,GeV) 
     >> iunit(_f0width,GeV) >> iunit(_pf0cc,GeV) >> iunit(_pf000,GeV) 
     >> iunit(_sigmamass,GeV) >> iunit(_sigmawidth,GeV) 
     >> iunit(_psigmacc,GeV) >> iunit(_psigma00,GeV) 
     >> iunit(_mpi0,GeV) >> iunit(_mpic,GeV) 
     >> iunit(_coupling,1/GeV) >> _rhomagP >> _rhophaseP
     >> _rhocoupP >> iunit(_rhomagD,1/GeV2) >> _rhophaseD 
     >> iunit(_rhocoupD,1/GeV2)>>iunit(_f2mag,1/GeV2) 
     >> _f2phase >> iunit(_f2coup,1/GeV2)
     >> _f0mag >> _f0phase >> _f0coup
     >> _sigmamag >> _sigmaphase >> _sigmacoup
     >> _localparameters 
     >> _zerowgts >> _onewgts >> _twowgts >> _threewgts
     >> _zeromax >> _onemax >> _twomax >> _threemax;
}
  
ClassDescription<a1ThreePionCLEODecayer> 
a1ThreePionCLEODecayer::inita1ThreePionCLEODecayer;
// Definition of the static class description member.
  
void a1ThreePionCLEODecayer::Init() {
  
  static ClassDocumentation<a1ThreePionCLEODecayer> documentation
    ("The a1ThreePionCLEODecayer class performs the decay of the "
     "a_1 to three pions using the model of CLEO",
     "The decay of a_1 to three pions was modelled after \\cite{Asner:1999kj}.",
     "%\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );

  static ParVector<a1ThreePionCLEODecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1ThreePionCLEODecayer::_rhomass,
     GeV, 0, ZERO, -10000*GeV, 10000*GeV, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1ThreePionCLEODecayer::_rhowidth,
     GeV, 0, ZERO, -10000*GeV, 10000*GeV, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &a1ThreePionCLEODecayer::_f2mass, GeV, 1.275*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &a1ThreePionCLEODecayer::_f2width, GeV, 0.185*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_0Mass
    ("f_0Mass",
     "The mass of the f_0 meson",
     &a1ThreePionCLEODecayer::_f0mass, GeV, 1.186*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_0Width
    ("f_0Width",
     "The width of the f_0 meson",
     &a1ThreePionCLEODecayer::_f0width, GeV, 0.350*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma meson",
     &a1ThreePionCLEODecayer::_sigmamass, GeV, 0.860*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma meson",
     &a1ThreePionCLEODecayer::_sigmawidth, GeV, 0.880*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1ThreePionCLEODecayer::_coupling, 1./GeV, 45.57/GeV, -0./GeV, 1000./GeV,
     false, false, false);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhomagP
    ("RhoPWaveMagnitude",
     "The magnitude of the couplings for the p-wave rho currents",
     &a1ThreePionCLEODecayer::_rhomagP,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhophaseP
    ("RhoPWavePhase",
     "The phase of the couplings for the p-wave rho currents",
     &a1ThreePionCLEODecayer::_rhophaseP,
     0, 0, 0, -pi, pi, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhomagD,
     1/MeV2, 0, ZERO, ZERO, 10000/MeV2, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhophaseD,
     0, 0, 0, -pi, pi, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &a1ThreePionCLEODecayer::_f0phase, 0.54*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &a1ThreePionCLEODecayer::_f2phase, 0.56*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &a1ThreePionCLEODecayer::_sigmaphase, 0.23*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef0Magnitude
    ("f0Magnitude",
     "The magnitude of the f_0 scalar current",
     &a1ThreePionCLEODecayer::_f0mag, 0.77, 0.0, 10,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,InvEnergy2> interfacef2Magnitude
    ("f2Magnitude",
     "The magnitude of the f_2 tensor current",
     &a1ThreePionCLEODecayer::_f2mag, 1./GeV2, 0.71/GeV2, ZERO, 10./GeV2,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacesigmaMagnitude
    ("sigmaMagnitude",
     "The magnitude of the sigma scalar current",
     &a1ThreePionCLEODecayer::_sigmamag, 2.1, 0.0, 10,
     false, false, true);

  static Switch<a1ThreePionCLEODecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1ThreePionCLEODecayer::_localparameters, true, false, false);
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

  static ParVector<a1ThreePionCLEODecayer,double> interfacezerowgts
    ("AllNeutralWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zerowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1ThreePionCLEODecayer::_threewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceZeroMax
    ("ZeroMax",
     "The maximum weight for the integration fo the channel a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zeromax, 0.0716349E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onemax,1.23756E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twomax,2.43819E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1ThreePionCLEODecayer::_threemax, 1.38754E3, 0.0, 10000.0,
     false, false, true);
}

double a1ThreePionCLEODecayer::me2(const int ichan,
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
  Energy2 q2=Q.mass2();
  // identify the mesons
  // calculate the invariants and form factors
  Lorentz5Momentum ps1=decay[1]->momentum()+decay[2]->momentum();
  Lorentz5Momentum ps2=decay[0]->momentum()+decay[2]->momentum();
  Lorentz5Momentum ps3=decay[0]->momentum()+decay[1]->momentum();
  ps1.rescaleMass();
  ps2.rescaleMass();
  ps3.rescaleMass();
  Energy2 s1=ps1.mass2(),s2=ps2.mass2(),s3=ps3.mass2();
  complex<InvEnergy> F1,F2,F3;
  formFactors(imode(),ichan,q2,s1,s2,s3,F1,F2,F3);
  // use the form-factors to compute the current
  LorentzPolarizationVector output=
    LorentzPolarizationVector(F1*decay[1]->momentum())-LorentzPolarizationVector(F1*decay[2]->momentum())
    -LorentzPolarizationVector(F2*decay[2]->momentum())+LorentzPolarizationVector(F2*decay[0]->momentum())
    +LorentzPolarizationVector(F3*decay[0]->momentum())-LorentzPolarizationVector(F3*decay[1]->momentum());
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix) 
    (*ME())(ix,0,0,0)=output.dot(_vectors[ix]);
  // answer
  double out = ME()->contract(_rho).real();
  // test of the answer
//   double test = threeBodyMatrixElement(imode(),sqr(inpart.mass()),
// 				       s3,s2,s1,decay[0]->mass(),decay[1]->mass(), 
// 				       decay[2]->mass());
//   if(ichan<0) cerr << "testing matrix element " << inpart.PDGName() << " -> "
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << decay[2]->PDGName() << out << " " << test << " " 
//        << (out-test)/(out+test) << "\n";  
  // return the answer
  return out;
}

// matrix element for the running a_1 width
double a1ThreePionCLEODecayer::
threeBodyMatrixElement(const int iopt,const Energy2 q2, const Energy2 s3,
		       const Energy2 s2,const Energy2 s1, const Energy m1, 
		       const Energy m2 ,const Energy m3) const {
  Energy2 m12=m1*m1,m22=m2*m2,m32=m3*m3;
  // calculate the form factors
  complex<InvEnergy> F1,F2,F3;
  formFactors(iopt,-1,q2,s1,s2,s3,F1,F2,F3);
  // analytic calculation of the matrix element
  double dot1=( F1*conj(F1)*(2.*m22+2.*m32-s1)+F2*conj(F2)*(2.*m12+2.*m32-s2)
		+F3*conj(F3)*(2.*m12+2.*m22-s3)-F1*conj(F2)*( s1+s2-s3-4.*m32)
		+F1*conj(F3)*( s1-s2+s3-4.*m22)-F2*conj(F3)*(-s1+s2+s3-4.*m12)).real();
  complex<Energy> dot2 = 0.5*(F1*(s3-m32-s2+m22)-F2*(s1-m12-s3+m32)+F3*(s2-m22-s1+m12));
  return (-dot1+(dot2*conj(dot2)).real()/q2)/3.;
}

WidthCalculatorBasePtr 
a1ThreePionCLEODecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit) {
    if(abs((**pit).id())==ParticleID::piplus) ++ncharged;
  }
  // integrator to perform the integral
  vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
  vector<int> intype;intype.push_back(2);intype.push_back(3);
  vector<Energy> inmass(2,_rhomass[0]),inwidth(2,_rhowidth[0]);
  vector<double> inpow(2,0.0);
  Energy m[3];
  m[0] = ncharged<2 ? _mpi0 : _mpic;
  m[1] = m[0];
  m[2] = (ncharged==0||ncharged==2) ? _mpi0 : _mpic;
  return new_ptr(ThreeBodyAllOnCalculator<a1ThreePionCLEODecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,ncharged,m[0],m[1],m[2]));
}

// calculate the form factos
void a1ThreePionCLEODecayer::formFactors(int iopt,int ichan,
					 Energy2 q2,Energy2 s1,Energy2 s2,
					 Energy2 s3,
					 complex<InvEnergy> & FF1,
					 complex<InvEnergy> & FF2,
					 complex<InvEnergy> & FF3) const {
  Complex F1, F2, F3;
  InvEnergy fact = _coupling;
  // a_1^0 pi0 pi0 pi0 mode
  if(iopt==0) {
    fact*=1./sqrt(6.);
    // compute the breit wigners we need
    Complex sigbws1 = sigmaBreitWigner(s1,1);
    Complex sigbws2 = sigmaBreitWigner(s2,1);
    Complex sigbws3 = sigmaBreitWigner(s3,1);
    Complex f0bws1  = f0BreitWigner(s1,1);
    Complex f0bws2  = f0BreitWigner(s2,1);
    Complex f0bws3  = f0BreitWigner(s3,1);
    Complex f2bws1  = f2BreitWigner(s1,1);
    Complex f2bws2  = f2BreitWigner(s2,1);
    Complex f2bws3  = f2BreitWigner(s3,1);
    if(ichan<0) {
      // the scalar terms
      F1=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	-2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      F2=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
      F3=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	+2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      // the tensor terms
      complex<Energy2> Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
      complex<Energy2> Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
      complex<Energy2> Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
      F1 += Complex(_f2coup*( 0.5*(s3-s2)*f2bws1-Dfact2+Dfact3));
      F2 += Complex(_f2coup*( 0.5*(s3-s1)*f2bws2-Dfact1+Dfact3));
      F3 += Complex(_f2coup*(-0.5*(s1-s2)*f2bws3-Dfact1+Dfact2));
    }
    else if(ichan==0) {
      F2=-2./3.*_sigmacoup*sigbws1;
      F3=-2./3.*_sigmacoup*sigbws1;
    }
    else if(ichan==1) {
      F1=-2./3.*_sigmacoup*sigbws2;
      F3=+2./3.*_sigmacoup*sigbws2;
    }
    else if(ichan==2) {
      F1= 2./3.*_sigmacoup*sigbws3;
      F2= 2./3.*_sigmacoup*sigbws3;
    }
    else if(ichan==3) {
      complex<Energy2> Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
      F1+=Complex(_f2coup*0.5*(s3-s2)*f2bws1);
      F2-=Complex(_f2coup*Dfact1); 
      F3-=Complex(_f2coup*Dfact1);
    }
    else if(ichan==4) {
      complex<Energy2> Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
      F2+=Complex(_f2coup*0.5*(s3-s1)*f2bws2);
      F1-=Complex(_f2coup*Dfact2);
      F3+=Complex(_f2coup*Dfact2);
    }
    else if(ichan==5) {
      complex<Energy2> Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
      F3+=Complex(-_f2coup*0.5*(s1-s2)*f2bws3);
      F1+=Complex(_f2coup*Dfact3);
      F2+=Complex(_f2coup*Dfact3);
    }
    else if(ichan==6) {
      F2=-2./3.*_f0coup*f0bws1;
      F3=-2./3.*_f0coup*f0bws1;
    }
    else if(ichan==7) {
      F1=-2./3.*_f0coup*f0bws2;
      F3=+2./3.*_f0coup*f0bws2;
    }
    else if(ichan==8) {
      F1= 2./3.*_f0coup*f0bws3;
      F2= 2./3.*_f0coup*f0bws3;
    }
  }
  // a_1^+ -> pi0 pi0 pi+
  else if(iopt==1) {
    fact *= 1./sqrt(2.);
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bw,sigbw,f2bw;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,1);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,1);
    }
    f0bw  = f0BreitWigner(s3,1);
    sigbw = sigmaBreitWigner(s3,1);
    f2bw  = f2BreitWigner(s3,1);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1+=_rhocoupP[ix]*rhos1bw[ix];
	F2+=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1+=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2+=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3+=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
      F1+=scalar;
      F2+=scalar;
      // the tensor terms
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      if(ires<_rhocoupP.size()) F1+=_rhocoupP[ires]*rhos1bw[ires];
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F2+=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      if(ires<_rhocoupP.size()) F2+=_rhocoupP[ires]*rhos2bw[ires];
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F1+=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3-=Complex(_rhocoupD[ires]*Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F1+=2./3.*_sigmacoup*sigbw;
      F2+=2./3.*_sigmacoup*sigbw;
    }
    else if(ichan==7) {
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan==8) {
      F1+=2./3.*_f0coup*f0bw;
      F2+=2./3.*_f0coup*f0bw;
    }
  }
  // a_1^0 ->pi+pi-pi0
  else if(iopt==2) {
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bw,sigbw,f2bw;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,1);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,1);
    }
    f0bw  =f0BreitWigner(s3,0);
    sigbw =sigmaBreitWigner(s3,0);
    f2bw  =f2BreitWigner(s3,0);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1+=_rhocoupP[ix]*rhos1bw[ix];
	F2+=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
      Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1+=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2+=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3+=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
      F1+=scalar;
      F2+=scalar;
      // the tensor terms
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      if(ires<_rhocoupP.size()) F1+=_rhocoupP[ires]*rhos1bw[ires];
      Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
      if(ires<_rhocoupD.size()) {
	F2+=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      if(ires<_rhocoupP.size()) F2+=_rhocoupP[ires]*rhos2bw[ires];
      Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
      if(ires<_rhocoupD.size()) {
	F1+=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3-=Complex(_rhocoupD[ires]*-Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F1+=2./3.*_sigmacoup*sigbw;
      F2+=2./3.*_sigmacoup*sigbw;
    }
    else if(ichan==7) {
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan==8) {
      F1+=2./3.*_f0coup*f0bw;
      F2+=2./3.*_f0coup*f0bw;
    }
  }
  // a_1^+ -> pi+ pi+ pi- mode
  else {
    fact *= 1./sqrt(2.);
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bws1,sigbws1,f2bws1,f0bws2,sigbws2,f2bws2;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,0);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,0);
    }
    f0bws1  =f0BreitWigner(s1,0);
    sigbws1 =sigmaBreitWigner(s1,0);
    f2bws1  =f2BreitWigner(s1,0);
    f0bws2  =f0BreitWigner(s2,0);
    sigbws2 =sigmaBreitWigner(s2,0);
    f2bws2  =f2BreitWigner(s2,0);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1-=_rhocoupP[ix]*rhos1bw[ix];
	F2-=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=1./3.*(s1-s3);
      Energy2 Dfact2=1./3.*(s2-s3);
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1-=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2-=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3-=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      F1-=2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      F2-=2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
      F3+=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	+2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      // the tensor terms
      complex<Energy2> sfact1 
	= 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      complex<Energy2> sfact2 
	= 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1+=Complex(_f2coup*(0.5*(s3-s2)*f2bws1-sfact2));
      F2+=Complex(_f2coup*(0.5*(s3-s1)*f2bws2-sfact1));
      F3+=Complex(_f2coup*(-sfact1+sfact2));
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      Energy2 Dfact2=1./3.*(s2-s3);
      if(ires<_rhocoupP.size()) F1-=_rhocoupP[ires]*rhos1bw[ires];
      if(ires<_rhocoupD.size()){
	F2-=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3-=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      Energy2 Dfact1=1./3.*(s1-s3);
      if(ires<_rhocoupP.size()) F2-=_rhocoupP[ires]*rhos2bw[ires];
      if(ires<_rhocoupD.size()) {
	F1-=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F2-=2./3.*_sigmacoup*sigbws1;
      F3-=2./3.*_sigmacoup*sigbws1;
    }
    else if(ichan==7) {
      F1-=2./3.*_sigmacoup*sigbws2;
      F3+=2./3.*_sigmacoup*sigbws2;
    }
    else if(ichan==8) {
      complex<Energy2> sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      F1+=Complex(_f2coup*0.5*(s3-s2)*f2bws1);
      F2-=Complex(_f2coup*sfact1);
      F3-=Complex(_f2coup*sfact1);
    }
    else if(ichan==9) {
      complex<Energy2> sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1-=Complex(_f2coup*sfact2);
      F2+=Complex(_f2coup*0.5*(s3-s1)*f2bws2);
      F3+=Complex(_f2coup*sfact2);
    }
    else if(ichan==10) {
      F2-=2./3.*_f0coup*f0bws1;
      F3-=2./3.*_f0coup*f0bws1;
    }
    else if(ichan==11) {
      F1-=2./3.*_f0coup*f0bws2;
      F3+=2./3.*_f0coup*f0bws2;
    }
  }
  FF1 = F1 * fact;
  FF2 = F2 * fact;
  FF3 = F3 * fact;
} 

// output the setup information for the particle database
void a1ThreePionCLEODecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // masses and widths of the intermediate particles
  output << "newdef " << name() << ":f_2Mass "    << _f2mass/GeV     << "\n";
  output << "newdef " << name() << ":f_2Width "   << _f2width/GeV    << "\n";
  output << "newdef " << name() << ":f_0Mass "    << _f0mass/GeV     << "\n";
  output << "newdef " << name() << ":f_0Width "   << _f0width/GeV    << "\n";
  output << "newdef " << name() << ":sigmaMass "  << _sigmamass/GeV  << "\n";
  output << "newdef " << name() << ":sigmaWidth " << _sigmawidth/GeV << "\n";
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
  }
  // couplings and phases for different channels
  output << "newdef " << name() << ":f0Phase " << _f0phase << "\n";
  output << "newdef " << name() << ":f2Phase " << _f2phase<< "\n";
  output << "newdef " << name() << ":sigmaPhase " << _sigmaphase<< "\n";
  output << "newdef " << name() << ":f0Magnitude " << _f0mag<< "\n";
  output << "newdef " << name() << ":f2Magnitude " << _f2mag*GeV2 << "\n";
  output << "newdef " << name() << ":sigmaMagnitude " << _sigmamag << "\n";
  output << "newdef " << name() << ":Coupling " << _coupling*GeV << "\n";
  for(unsigned int ix=0;ix<_rhomagP.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoPWaveMagnitude " << ix << " " 
		    << _rhomagP[ix] << "\n";
    else     output << "insert " << name() << ":RhoPWaveMagnitude " << ix << " " 
		    << _rhomagP[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhophaseP.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoPWavePhase " << ix << " " 
		    << _rhophaseP[ix] << "\n";
    else     output << "insert " << name() << ":RhoPWavePhase " << ix << " " 
		    << _rhophaseP[ix] << "\n";
  }  
  for(unsigned int ix=0;ix<_rhomagD.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoDWaveMagnitude " << ix << " " 
		    << _rhomagD[ix]*MeV2 << "\n";
    else     output << "insert " << name() << ":RhoDWaveMagnitude " << ix << " " 
		    << _rhomagD[ix]*MeV2 << "\n";
  }
  for(unsigned int ix=0;ix<_rhophaseD.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoDWavePhase " << ix << " " 
		    << _rhophaseD[ix] << "\n";
    else     output << "insert " << name() << ":RhoDWavePhase " << ix << " " 
		    << _rhophaseD[ix] << "\n";
  }
  // use local values of the masses etc.
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
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
  // maximum weights for the different  channels
  output << "newdef " << name() << ":ZeroMax "  << _zeromax  << "\n";
  output << "newdef " << name() << ":OneMax "   << _onemax   << "\n";
  output << "newdef " << name() << ":TwoMax "   << _twomax   << "\n";
  output << "newdef " << name() << ":ThreeMax " << _threemax << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
