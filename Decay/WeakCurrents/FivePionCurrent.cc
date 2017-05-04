// -*- C++ -*-
//
// FivePionCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FivePionCurrent class.
//

#include "FivePionCurrent.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

FivePionCurrent::FivePionCurrent() {
  // set the number of modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(3);
  // masses of the intermediates
  _rhomass    = 776*MeV;
  _a1mass     = 1260*MeV;
  _omegamass  = 782*MeV;
  _sigmamass  = 800*MeV;
  // widths of the intermediates
  _rhowidth   = 150*MeV;
  _a1width    = 400*MeV;
  _omegawidth = 8.5*MeV;
  _sigmawidth = 600*MeV;
  // use local values of the resonance parameters
  _localparameters=true;
  // include the rho Breit-Wigners in omega decay
  _rhoomega = true;
  // Normalisation parameters for the different currents
  _c =4.*GeV2;
  _c0=3.;
  // various meson coupling constants
  _fomegarhopi=0.07/MeV;
  _grhopipi=6.0;
  _garhopi=6.*GeV;
  _faaf=4.*GeV;
  _ffpipi=5.*GeV;
  _presigma = ZERO;
  _preomega = ZERO;
}

inline void FivePionCurrent::doinit() {
  WeakDecayCurrent::doinit();
  if(!_localparameters) {
    _rhomass    = getParticleData(ParticleID::rhominus)->mass();
    _rhowidth   = getParticleData(ParticleID::rhominus)->width();
    _omegamass  = getParticleData(ParticleID::omega)->mass();
    _omegawidth = getParticleData(ParticleID::omega)->width();
    _sigmamass  = getParticleData(9000221)->mass();
    _sigmawidth = getParticleData(9000221)->width();
    _a1mass    = getParticleData(ParticleID::a_1minus)->mass();
    _a1width   = getParticleData(ParticleID::a_1minus)->width();
  }
  // prefactors
  _presigma =  _c/sqr(sqr(_a1mass)*_sigmamass*_rhomass)*_faaf*_ffpipi*
    _garhopi*_grhopipi;      
  _preomega = _c0*_fomegarhopi*sqr(_grhopipi/(sqr(_rhomass)*_omegamass));
}

void FivePionCurrent::persistentOutput(PersistentOStream & os) const {
  static const InvEnergy7 InvGeV7 = pow<-7,1>(GeV);
  static const InvEnergy3 InvGeV3 = pow<-3,1>(GeV);
  os << ounit(_rhomass,GeV)  << ounit(_a1mass,GeV) << ounit(_omegamass,GeV) 
     << ounit(_sigmamass,GeV) << ounit(_rhowidth,GeV)  
     << ounit(_a1width,GeV) << ounit(_omegawidth,GeV) << ounit(_sigmawidth,GeV) 
     << _localparameters << ounit(_c,GeV2) << _c0 
     << ounit(_fomegarhopi,1/GeV) << _grhopipi << ounit(_garhopi,GeV) 
     << ounit(_faaf,GeV) << ounit(_ffpipi,GeV)
     << ounit(_preomega,InvGeV7) << ounit(_presigma,InvGeV3) << _rhoomega;
}

void FivePionCurrent::persistentInput(PersistentIStream & is, int) {
  static const InvEnergy7 InvGeV7 = pow<-7,1>(GeV);
  static const InvEnergy3 InvGeV3 = pow<-3,1>(GeV);
  is >> iunit(_rhomass,GeV)  >> iunit(_a1mass,GeV) >> iunit(_omegamass,GeV) 
     >> iunit(_sigmamass,GeV) >> iunit(_rhowidth,GeV)  
     >> iunit(_a1width,GeV) >> iunit(_omegawidth,GeV) >> iunit(_sigmawidth,GeV) 
     >> _localparameters >> iunit(_c,GeV2) >> _c0 
     >> iunit(_fomegarhopi,1/GeV) >> _grhopipi >> iunit(_garhopi,GeV) 
     >> iunit(_faaf,GeV) >> iunit(_ffpipi,GeV) 
     >> iunit(_preomega,InvGeV7) >> iunit(_presigma,InvGeV3) >> _rhoomega;
}

ClassDescription<FivePionCurrent> FivePionCurrent::initFivePionCurrent;
// Definition of the static class description member.

void FivePionCurrent::Init() {

  static ClassDocumentation<FivePionCurrent> documentation
    ("The FivePionCurrent class implements the model of hep-ph/0602162",
     "The model of \\cite{Kuhn:2006nw} was used for the hadronic five pion current.",
     "\\bibitem{Kuhn:2006nw} J.~H.~Kuhn and Z.~Was, hep-ph/0602162, (2006).");

  static Parameter<FivePionCurrent,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho meson",
     &FivePionCurrent::_rhomass, MeV, 776*MeV, 500*MeV, 1000*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The mass of the a_1 meson",
     &FivePionCurrent::_a1mass, MeV, 1260*MeV, 1000*MeV, 1500*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &FivePionCurrent::_omegamass, MeV, 782*MeV, 600*MeV, 900*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceSigmaMass
    ("SigmaMass",
     "The mass of the sigma meson",
     &FivePionCurrent::_sigmamass, MeV, 800*MeV, 400*MeV, 1200*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho meson",
     &FivePionCurrent::_rhowidth, MeV, 150*MeV, 100*MeV, 300*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The width of the a_1 meson",
     &FivePionCurrent::_a1width, MeV, 400*MeV, 100*MeV, 800*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &FivePionCurrent::_omegawidth, MeV, 8.5*MeV, 1.0*MeV, 20.0*MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceSigmaWidth
    ("SigmaWidth",
     "The width of the sigma meson",
     &FivePionCurrent::_sigmawidth, MeV, 600*MeV, 100*MeV, 1200*MeV,
     false, false, Interface::limited);

  static Switch<FivePionCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the meson masses and widths or those from the"
     " ParticleData objects.",
     &FivePionCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static Switch<FivePionCurrent,bool> interfaceRhoOmega
    ("RhoOmega",
     "Option for the treatment of the rho Breit-Wigners in the omega decay",
     &FivePionCurrent::_rhoomega, true, false, false);
  static SwitchOption interfaceRhoOmegaInclude
    (interfaceRhoOmega,
     "Yes",
     "Include the rho Breit-Wigners",
     true);
  static SwitchOption interfaceRhoOmegaOmit
    (interfaceRhoOmega,
     "No",
     "Don't include the rho Breit-Wigners",
     false);

  static Parameter<FivePionCurrent,Energy2> interfaceC
    ("C",
     "The normalisation parameter for the a_1 sigma current",
     &FivePionCurrent::_c, GeV2, 4.0*GeV2, 0.1*GeV2, 20.0*GeV2,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,double> interfaceC0
    ("C0",
     "The normalisation constant for the omega-rho current",
     &FivePionCurrent::_c0, 3., 0.1, 10.0,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,InvEnergy> interfacegomegarhopi
    ("fomegarhopi",
     "The coupling of omega-rho-pi",
     &FivePionCurrent::_fomegarhopi, 1./MeV, 0.07/MeV, 0.01/MeV, 0.2/MeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,double> interfacegrhopipi
    ("grhopipi",
     "The coupling for rho-pi-pi",
     &FivePionCurrent::_grhopipi, 6.0, 1.0, 20.0,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfacegarhopi
    ("garhopi",
     "The coupling of a-rho-pi",
     &FivePionCurrent::_garhopi, GeV, 6.0*GeV, 0.1*GeV, 20.0*GeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfacefaaf
    ("faaf",
     "The coupling of a-a-f",
     &FivePionCurrent::_faaf, GeV, 4.0*GeV, 0.1*GeV, 20.0*GeV,
     false, false, Interface::limited);

  static Parameter<FivePionCurrent,Energy> interfaceffpipi
    ("ffpipi",
     "The coupling of f-pi-pi",
     &FivePionCurrent::_ffpipi, GeV, 5.0*GeV, 0.1*GeV, 20.0*GeV,
     false, false, Interface::limited);
}

bool FivePionCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check five products
  if(id.size()!=5){return false;}
  int npiminus=0,npiplus=0,npi0=0;
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID:: piplus){++npiplus;}
    else if(id[ix]==ParticleID::piminus){++npiminus;}
    else if(id[ix]==ParticleID::pi0){++npi0;}
  }
  if(npiplus>npiminus) swap(npiplus,npiminus);
  if(     npiminus==3&&npiplus==2&&npi0==0) allowed=true;
  else if(npiminus==2&&npiplus==1&&npi0==2) allowed=true;
  else if(npiminus==1&&npiplus==0&&npi0==4) allowed=true;
  return allowed;
}

bool FivePionCurrent::createMode(int icharge, unsigned int imode,
				 DecayPhaseSpaceModePtr mode,
				 unsigned int iloc,unsigned int ires,
				 DecayPhaseSpaceChannelPtr phase,Energy upp) {
  // check the charge
  if(abs(icharge)!=3) return false;
  // check that the modes are kinematical allowed
  Energy min(ZERO);
  // 3 pi- 2pi+
  if(imode==0) {
    min=5.*getParticleData(ParticleID::piplus)->mass();
  }
  // 2pi- pi+ 2pi0
  else if(imode==1) {
    min=3.*getParticleData(ParticleID::piplus)->mass()
      +2.*getParticleData(ParticleID::pi0)->mass();
  }
  // pi- 4pi0
  else {
    min=   getParticleData(ParticleID::piplus)->mass()
      +4.*getParticleData(ParticleID::pi0)->mass();
  }
  if(min>upp) return false;
  // intermediates for the channels
  tPDPtr omega(getParticleData(ParticleID::omega)),rhop,rhom,
    rho0(getParticleData(ParticleID::rho0)),a1m,a10(getParticleData(ParticleID::a_10)),
    sigma(getParticleData(9000221));
  if(icharge==3) {
    rhop = getParticleData(ParticleID::rhominus);
    rhom = getParticleData(ParticleID::rhoplus);
    a1m  = getParticleData(ParticleID::a_1plus);
  }
  else {
    rhop = getParticleData(ParticleID::rhoplus);
    rhom = getParticleData(ParticleID::rhominus);
    a1m  = getParticleData(ParticleID::a_1minus);
  }
  DecayPhaseSpaceChannelPtr newchannel;
  // all charged mode
  if(imode==0) {
    if(sigma) {
      // 1st two a_1 sigma channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc  );
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+1, iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+2);
      mode->addChannel(newchannel);
      // 2nd two a_1 sigma channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc  );
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+1, iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc  );
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+4, iloc+2);
      mode->addChannel(newchannel);
      // 3rd two a_1 sigma channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc  );
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+4, iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+2);
      mode->addChannel(newchannel);
      // 4th two a_1 sigma channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc  );
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+1, iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+3);
      mode->addChannel(newchannel);
      // 5th two a_1 sigma channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc  );
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+1, iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc  );
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+4, iloc+3);
      mode->addChannel(newchannel);
      // 6th two a_1 sigma channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc  );
      newchannel->addIntermediate(rho0   ,0,0.0, iloc+4, iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+3);
      mode->addChannel(newchannel);
    }
  }
  // 2 neutral mode
  else if(imode==1) {
    // first three omega channels
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+4);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+2);
    newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+1);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+4);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+1);
    newchannel->addIntermediate(rhop   ,0,0.0, iloc  , iloc+2);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+4);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc  );
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
    mode->addChannel(newchannel);
    // second three omega channels
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+4);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+2);
    newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+3);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+4);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+3);
    newchannel->addIntermediate(rhop   ,0,0.0, iloc  , iloc+2);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+4);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc  );
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
    mode->addChannel(newchannel);
    // third  three omega channels
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+4);
    newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+1);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+1);
    newchannel->addIntermediate(rhop   ,0,0.0, iloc  , iloc+4);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc  );
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+4);
    mode->addChannel(newchannel);
    // fourth three omega channels
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+4);
    newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+3);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc+3);
    newchannel->addIntermediate(rhop   ,0,0.0, iloc  , iloc+4);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
    newchannel->addIntermediate(omega  ,0,0.0,-ires-3, iloc  );
    newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+4);
    mode->addChannel(newchannel);
    // first two sigma a_1 channels
    if(sigma) {
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+2);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+4);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
      mode->addChannel(newchannel);
      // second two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+3);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+2);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+4);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+3);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
      mode->addChannel(newchannel);
      // third two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+2, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+3);
      newchannel->addIntermediate(rho0   ,0,0.0, iloc  , iloc+1);
      mode->addChannel(newchannel);
    }
  }
  // 4 neutral mode
  else {
    if(sigma) {
      // 1st two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc  , iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc  );
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
      mode->addChannel(newchannel);
      // 2nd two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc  );
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+4, iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc  );
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
      mode->addChannel(newchannel);
      // 3rd two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+1);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+3);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+1, iloc+2);
      mode->addChannel(newchannel);
      // 4th two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+1, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+3);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc  , iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+1, iloc+4);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc  );
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
      mode->addChannel(newchannel);
      // 5th two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc  , iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc+3, iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc  );
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+4, iloc+2);
      mode->addChannel(newchannel);
      // 6th two sigma a_1 channels
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+3);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+4, iloc+2);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-1,-ires-2);
      newchannel->addIntermediate(sigma  ,0,0.0, iloc  , iloc+1);
      newchannel->addIntermediate(a1m    ,0,0.0,-ires-3, iloc+4);
      newchannel->addIntermediate(rhom   ,0,0.0, iloc+3, iloc+2);
      mode->addChannel(newchannel);
    }
  }
  // reset the parameters of the resonances if using local values
  if(_localparameters) {
    mode->resetIntermediate(rhom,_rhomass,_rhowidth);
    mode->resetIntermediate(rhop,_rhomass,_rhowidth);
    mode->resetIntermediate(rho0,_rhomass,_rhowidth);
    mode->resetIntermediate(omega,_omegamass,_omegawidth);
    mode->resetIntermediate(a1m,_a1mass,_a1width);
    mode->resetIntermediate(a10,_a1mass,_a1width);
    if(sigma) mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector FivePionCurrent::particles(int icharge, unsigned int imode,int,int) {
  // particle data objects for the pions
  tPDPtr piplus (getParticleData(ParticleID::piplus ));
  tPDPtr pi0    (getParticleData(ParticleID::pi0    ));
  tPDPtr piminus(getParticleData(ParticleID::piminus));
  if(icharge==3) swap(piplus,piminus);
  tPDVector output(5);
  // all charged
  if(imode==0) {
    output[0]=piminus;
    output[1]=piminus;
    output[2]=piplus;;
    output[3]=piplus;
    output[4]=piminus;
  }
  // two neutral
  else if(imode==1) {
    output[0]=piplus;
    output[1]=piminus;
    output[2]=pi0;;
    output[3]=piminus;
    output[4]=pi0;
  }
  // four neutral
  else {
    output[0]=pi0;
    output[1]=pi0;
    output[2]=piminus;;
    output[3]=pi0;
    output[4]=pi0;
  }
  return output;
}

// the decay mode
unsigned int FivePionCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::pi0) ++npi;
  }
  return npi/2;
}

// output the information for the database
void FivePionCurrent::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::FivePionCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":RhoMass "    << _rhomass/MeV << "\n";
  output << "newdef " << name() << ":A1Mass  "    << _a1mass/MeV  << "\n";
  output << "newdef " << name() << ":SigmaMass  " << _sigmamass/MeV  << "\n";
  output << "newdef " << name() << ":OmegaMass  " << _omegamass/MeV  << "\n";
  output << "newdef " << name() << ":RhoWidth "    << _rhowidth/MeV << "\n";
  output << "newdef " << name() << ":A1Width  "    << _a1width/MeV  << "\n";
  output << "newdef " << name() << ":SigmaWidth  " << _sigmawidth/MeV  << "\n";
  output << "newdef " << name() << ":OmegaWidth  " << _omegawidth/MeV  << "\n";
  output << "newdef " << name() << ":LocalParameters " <<  _localparameters << "\n";
  output << "newdef " << name() << ":RhoOmega " << _rhoomega << "\n";
  output << "newdef " << name() << ":C " << _c/GeV2 << "\n";
  output << "newdef " << name() << ":C0 " << _c0 << "\n";
  output << "newdef " << name() << ":fomegarhopi " <<_fomegarhopi*MeV << "\n";
  output << "newdef " << name() << ":grhopipi " <<_grhopipi << "\n";
  output << "newdef " << name() << ":garhopi " << _garhopi/GeV << "\n";
  output << "newdef " << name() << ":faaf " <<_faaf/GeV << "\n";
  output << "newdef " << name() << ":ffpipi " << _ffpipi/GeV << "\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
}

vector<LorentzPolarizationVectorE> 
FivePionCurrent::current(const int imode,const int ichan,
			 Energy & scale,const ParticleVector & decay,
			 DecayIntegrator::MEOption meopt) const {
  useMe();
  if(meopt==DecayIntegrator::Terminate) {
    for(unsigned int ix=0;ix<5;++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  LorentzVector<complex<InvEnergy2> > output;
  Lorentz5Momentum q1(decay[0]->momentum());
  Lorentz5Momentum q2(decay[1]->momentum());
  Lorentz5Momentum q3(decay[2]->momentum());
  Lorentz5Momentum q4(decay[3]->momentum());
  Lorentz5Momentum q5(decay[4]->momentum());
  // total momentum
  Lorentz5Momentum Q(q1+q2+q3+q4+q5);
  Q.rescaleMass();
  scale=Q.mass();
  // decide which decay mode
  if(imode==0) {
    if(ichan<0) {
      output=
	a1SigmaCurrent(0,Q,q1,q2,q3,q4,q5)+
	a1SigmaCurrent(0,Q,q5,q2,q3,q4,q1)+
	a1SigmaCurrent(0,Q,q1,q5,q3,q4,q2)+
	a1SigmaCurrent(0,Q,q1,q2,q4,q3,q5)+
	a1SigmaCurrent(0,Q,q5,q2,q4,q3,q1)+
	a1SigmaCurrent(0,Q,q1,q5,q4,q3,q2);
    }
    else if(ichan==0 ) output=a1SigmaCurrent(2,Q,q1,q2,q3,q4,q5);
    else if(ichan==1 ) output=a1SigmaCurrent(1,Q,q1,q2,q3,q4,q5);
    else if(ichan==2 ) output=a1SigmaCurrent(2,Q,q5,q2,q3,q4,q1);
    else if(ichan==3 ) output=a1SigmaCurrent(1,Q,q5,q2,q3,q4,q1);
    else if(ichan==4 ) output=a1SigmaCurrent(2,Q,q1,q5,q3,q4,q2);
    else if(ichan==5 ) output=a1SigmaCurrent(1,Q,q1,q5,q3,q4,q2);
    else if(ichan==6 ) output=a1SigmaCurrent(2,Q,q1,q2,q4,q3,q5);
    else if(ichan==7 ) output=a1SigmaCurrent(1,Q,q1,q2,q4,q3,q5);
    else if(ichan==8 ) output=a1SigmaCurrent(2,Q,q5,q2,q4,q3,q1);
    else if(ichan==9 ) output=a1SigmaCurrent(1,Q,q5,q2,q4,q3,q1);
    else if(ichan==10) output=a1SigmaCurrent(2,Q,q1,q5,q4,q3,q2);
    else if(ichan==11) output=a1SigmaCurrent(1,Q,q1,q5,q4,q3,q2);
    // identical particle symmetry factor
    output/=sqrt(12.);
  }
  else if(imode==1) {
    if(ichan<0) {
      output=
	rhoOmegaCurrent(0,Q,q1,q2,q3,q4,q5)
	+rhoOmegaCurrent(0,Q,q1,q4,q3,q2,q5)
	+rhoOmegaCurrent(0,Q,q1,q2,q5,q4,q3)
	+rhoOmegaCurrent(0,Q,q1,q4,q5,q2,q3)
	+a1SigmaCurrent(0,Q,q2,q4,q1,q3,q5)
	+a1SigmaCurrent(0,Q,q3,q5,q2,q1,q4)
	+a1SigmaCurrent(0,Q,q3,q5,q4,q1,q2);
    }
    else if(ichan==0 ) output=rhoOmegaCurrent(3,Q,q1,q2,q3,q4,q5);
    else if(ichan==1 ) output=rhoOmegaCurrent(2,Q,q1,q2,q3,q4,q5);
    else if(ichan==2 ) output=rhoOmegaCurrent(1,Q,q1,q2,q3,q4,q5);
    else if(ichan==3 ) output=rhoOmegaCurrent(3,Q,q1,q4,q3,q2,q5);
    else if(ichan==4 ) output=rhoOmegaCurrent(2,Q,q1,q4,q3,q2,q5);
    else if(ichan==5 ) output=rhoOmegaCurrent(1,Q,q1,q4,q3,q2,q5);
    else if(ichan==6 ) output=rhoOmegaCurrent(3,Q,q1,q2,q5,q4,q3);
    else if(ichan==7 ) output=rhoOmegaCurrent(2,Q,q1,q2,q5,q4,q3);
    else if(ichan==8 ) output=rhoOmegaCurrent(1,Q,q1,q2,q5,q4,q3);
    else if(ichan==9 ) output=rhoOmegaCurrent(3,Q,q1,q4,q5,q2,q3);
    else if(ichan==10) output=rhoOmegaCurrent(2,Q,q1,q4,q5,q2,q3);
    else if(ichan==11) output=rhoOmegaCurrent(1,Q,q1,q4,q5,q2,q3);
    else if(ichan==12) output=a1SigmaCurrent(2,Q,q3,q5,q4,q1,q2);
    else if(ichan==13) output=a1SigmaCurrent(1,Q,q3,q5,q4,q1,q2);
    else if(ichan==14) output=a1SigmaCurrent(2,Q,q3,q5,q2,q1,q4);
    else if(ichan==15) output=a1SigmaCurrent(1,Q,q3,q5,q2,q1,q4);
    else if(ichan==16) output=a1SigmaCurrent(2,Q,q2,q4,q1,q3,q5);
    else if(ichan==17) output=a1SigmaCurrent(1,Q,q2,q4,q1,q3,q5);
    // identical particle symmetry factor
    output/=2.;
  }
  else if(imode==2) {
    if(ichan<0) {
      output=
	a1SigmaCurrent(0,Q,q1,q2,q3,q4,q5)+
	a1SigmaCurrent(0,Q,q5,q2,q3,q4,q1)+
	a1SigmaCurrent(0,Q,q2,q4,q3,q1,q5)+
	a1SigmaCurrent(0,Q,q1,q4,q3,q2,q5)+
	a1SigmaCurrent(0,Q,q1,q5,q3,q4,q2)+
	a1SigmaCurrent(0,Q,q4,q5,q3,q1,q2);
    }
    else if(ichan==0 ) output=a1SigmaCurrent(1,Q,q1,q2,q3,q4,q5);
    else if(ichan==1 ) output=a1SigmaCurrent(2,Q,q1,q2,q3,q4,q5);
    else if(ichan==2 ) output=a1SigmaCurrent(1,Q,q5,q2,q3,q4,q1);
    else if(ichan==3 ) output=a1SigmaCurrent(2,Q,q5,q2,q3,q4,q1);
    else if(ichan==4 ) output=a1SigmaCurrent(2,Q,q2,q4,q3,q1,q5);
    else if(ichan==5 ) output=a1SigmaCurrent(1,Q,q2,q4,q3,q1,q5);
    else if(ichan==6 ) output=a1SigmaCurrent(1,Q,q1,q4,q3,q2,q5);
    else if(ichan==7 ) output=a1SigmaCurrent(2,Q,q1,q4,q3,q2,q5);
    else if(ichan==8 ) output=a1SigmaCurrent(1,Q,q1,q5,q3,q4,q2);
    else if(ichan==9 ) output=a1SigmaCurrent(2,Q,q1,q5,q3,q4,q2);
    else if(ichan==10) output=a1SigmaCurrent(2,Q,q4,q5,q3,q1,q2);
    else if(ichan==11) output=a1SigmaCurrent(1,Q,q4,q5,q3,q1,q2);
    // identical particle symmetry factor
    output/=sqrt(24.);
  }
  else {
    throw DecayIntegratorError() << "Unknown decay mode in the " 
				 << "FivePionCurrent::"
				 << "hadronCurrent()" << Exception::abortnow;
  }
  // normalise and return the current
  return vector<LorentzPolarizationVectorE>(1, output * pow<3,1>(scale));
}
