// -*- C++ -*-
//
// DtoKPiPiE691.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiE691 class.
//

#include "DtoKPiPiE691.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

DtoKPiPiE691::DtoKPiPiE691() {
  // amplitudes and phases for D+ -> K-pi+pi+
  _a1NR    = 1.  ; _phi1NR    =    0.;
  _a1K892  = 0.78; _phi1K892  = - 60.;
  _a1K1430 = 0.53; _phi1K1430 =  132.;
  _a1K1680 = 0.47; _phi1K1680 = - 51.;
  // amplitudes and phases for D0 -> K-pi+pi0
  _a2NR    = 1.00; _phi2NR    =    0.;
  _a2K8920 = 3.19; _phi2K8920 =  167.;
  _a2K892m = 2.96; _phi2K892m = -112.;
  _a2rho   = 8.56; _phi2rho   =   40.;
  // amplitudes and phases for D0 -> Kbar0 pi+pi-
  _a3NR    = 1.00; _phi3NR    =    0.;
  _a3K892  = 2.31; _phi3K892  =  109.;
  _a3rho   = 1.59; _phi3rho   = -123.;
  // masses and widths
  _localparameters=true;
  _mK8920 = 0.8961 *GeV; _wK8920 = 0.0505*GeV;
  _mK892m = 0.89159*GeV; _wK892m = 0.0498*GeV;
  _mK1680 = 1.714  *GeV; _wK1680 = 0.323 *GeV;
  _mK1430 = 1.429  *GeV; _wK1430 = 0.287 *GeV;
  _mrho0  = 0.7681 *GeV; _wrho0  = 0.1515*GeV;
  _mrhop  = 0.7681 *GeV; _wrhop  = 0.1515*GeV;
  // generate the intermediate particles
  generateIntermediates(true);
}

void DtoKPiPiE691::doinit() {
  DecayIntegrator::doinit();
  // complex amplitudes calculated from magnitudes and phases
  double fact = Constants::pi/180.;
  // D+ -> K-pi+pi+
  _c1NR   = _a1NR   *Complex(cos(_phi1NR   *fact),sin(_phi1NR   *fact)); 
  _c1K892 = _a1K892 *Complex(cos(_phi1K892 *fact),sin(_phi1K892 *fact));
  _c1K1430= _a1K1430*Complex(cos(_phi1K1430*fact),sin(_phi1K1430*fact));
  _c1K1680= _a1K1680*Complex(cos(_phi1K1680*fact),sin(_phi1K1680*fact));
  // D0 -> K-pi+pi0
  _c2NR   = _a2NR   *Complex(cos(_phi2NR   *fact),sin(_phi2NR   *fact));
  _c2K8920= _a2K8920*Complex(cos(_phi2K8920*fact),sin(_phi2K8920*fact));
  _c2K892m= _a2K892m*Complex(cos(_phi2K892m*fact),sin(_phi2K892m*fact));
  _c2rho  = _a2rho  *Complex(cos(_phi2rho  *fact),sin(_phi2rho  *fact));
  // D0 -> Kbar0pi+pi-
  _c3NR   = _a3NR   *Complex(cos(_phi3NR   *fact),sin(_phi3NR   *fact));
  _c3K892 = _a3K892 *Complex(cos(_phi3K892 *fact),sin(_phi3K892 *fact));
  _c3rho  = _a3rho  *Complex(cos(_phi3rho  *fact),sin(_phi3rho  *fact));
  // intermediate resonances
  tPDPtr k8920 = getParticleData(ParticleID::Kstarbar0 );
  tPDPtr k892m = getParticleData(ParticleID::Kstarminus);
  tPDPtr k1430 = getParticleData(ParticleID::Kstar_0bar0);
  tPDPtr k1680 = getParticleData(-30313);
  tPDPtr rho0  = getParticleData(ParticleID::rho0);
  tPDPtr rhop  = getParticleData(ParticleID::rhoplus);
  // D+ -> K-pi+pi+
  tPDVector extpart(4);
  extpart[0]=getParticleData(ParticleID::Dplus);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piplus);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  int ix=0;
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k8920,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k8920,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k1430) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1430,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1430,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k1680) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1680,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1680,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // add the mode
  vector<double> wtemp;
  if(ix<=int(_weights.size())) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit,wit+ix);
  }
  else {
    wtemp=vector<double>(ix,1./double(ix));
  }
  if(_maxwgt.empty()) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[0],wtemp);
  // D0 -> K-pi+pi0
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::pi0);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  int iy=ix;
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k8920,0,0., 1,2);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(rhop) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rhop,0,0., 2,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  // add the mode
  if(iy<=int(_weights.size())) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit+ix,wit+iy);
  }
  else {
    wtemp=vector<double>(iy-ix,1./double(iy-ix));
  }
  if(_maxwgt.size()<2) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[1],wtemp);
  // D0 -> Kbar0 pi+ pi-
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kbar0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piminus);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  ix=iy;
  if(rho0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho0,0,0., 2,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  // add the mode
  if(iy<=int(_weights.size())) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit+ix,wit+iy);
  }
  else {
    wtemp=vector<double>(iy-ix,1./double(iy-ix));
  }
  if(_maxwgt.size()<3) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[2],wtemp);
  // reset the resonance parameters in the integration if needed
  if(_localparameters) {
    resetIntermediate(k8920,_mK8920,_wK8920);
    resetIntermediate(k892m,_mK892m,_wK892m);
    resetIntermediate(k1430,_mK1680,_wK1680);
    resetIntermediate(k1680,_mK1430,_wK1430);
    resetIntermediate(rho0 ,_mrho0 ,_wrho0 );
    resetIntermediate(rhop ,_mrhop ,_wrhop );
  }
  // get values from the ParticleData objects if needed
  else {
    _mK8920 = k8920->mass();
    _mK892m = k892m->mass();
    _mK1680 = k1430->mass();
    _mK1430 = k1680->mass();
    _mrho0  = rho0 ->mass();
    _mrhop  = rhop ->mass();
    _wK8920 = k8920->width();
    _wK892m = k892m->width();
    _wK1680 = k1430->width();
    _wK1430 = k1680->width();
    _wrho0  = rho0 ->width();
    _wrhop  = rhop ->width();
  }
}

void DtoKPiPiE691::persistentOutput(PersistentOStream & os) const {
  os << _a1NR << _phi1NR << _a1K892 << _phi1K892 << _a1K1430 << _phi1K1430 
     << _a1K1680 << _phi1K1680 << _a2NR << _phi2NR << _a2K8920 << _phi2K8920 
     << _a2K892m << _phi2K892m << _a2rho << _phi2rho << _a3NR << _phi3NR 
     << _a3K892 << _phi3K892 << _a3rho << _phi3rho << _c1NR << _c1K892 << _c1K1430 
     << _c1K1680 << _c2NR << _c2K8920 << _c2K892m << _c2rho << _c3NR << _c3K892 
     << _c3rho << _localparameters << ounit(_mK8920,GeV) << ounit(_wK8920,GeV) 
     << ounit(_mK892m,GeV) << ounit(_wK892m,GeV) << ounit(_mK1680,GeV) 
     << ounit(_wK1680,GeV) << ounit(_mK1430,GeV) << ounit(_wK1430,GeV) 
     << ounit(_mrho0 ,GeV) << ounit(_wrho0 ,GeV) << ounit(_mrhop ,GeV) 
     << ounit(_wrhop ,GeV) << _maxwgt << _weights;
}

void DtoKPiPiE691::persistentInput(PersistentIStream & is, int) {
  is >> _a1NR >> _phi1NR >> _a1K892 >> _phi1K892 >> _a1K1430 >> _phi1K1430 
     >> _a1K1680 >> _phi1K1680 >> _a2NR >> _phi2NR >> _a2K8920 >> _phi2K8920 
     >> _a2K892m >> _phi2K892m >> _a2rho >> _phi2rho >> _a3NR >> _phi3NR 
     >> _a3K892 >> _phi3K892 >> _a3rho >> _phi3rho >> _c1NR >> _c1K892 >> _c1K1430 
     >> _c1K1680 >> _c2NR >> _c2K8920 >> _c2K892m >> _c2rho >> _c3NR >> _c3K892 
     >> _c3rho >> _localparameters >> iunit(_mK8920,GeV) >> iunit(_wK8920,GeV) 
     >> iunit(_mK892m,GeV) >> iunit(_wK892m,GeV) >> iunit(_mK1680,GeV) 
     >> iunit(_wK1680,GeV) >> iunit(_mK1430,GeV) >> iunit(_wK1430,GeV) 
     >> iunit(_mrho0 ,GeV) >> iunit(_wrho0 ,GeV) >> iunit(_mrhop ,GeV) 
     >> iunit(_wrhop ,GeV) >> _maxwgt >> _weights;
}

ClassDescription<DtoKPiPiE691> DtoKPiPiE691::initDtoKPiPiE691;
// Definition of the static class description member.

void DtoKPiPiE691::Init() {

  static ClassDocumentation<DtoKPiPiE691> documentation
    ("The DtoKPiPiE691 class implements the model of E691 for D -> K pi pi"
     "Dalitz decays",
     "The model of \\cite{Anjos:1992kb} for $D\\to K\\pi\\pi$ decays was used",
     "\\bibitem{Anjos:1992kb} J.~C.~Anjos {\\it et al.}  [E691 Collaboration],"
     "Phys.\\ Rev.\\  D {\\bf 48} (1993) 56.");

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipNonResonantMagnitude
    ("KmPipPipNonResonantMagnitude",
     "The magnitude of the non-resonant component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_a1NR, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipNonResonantPhase
    ("KmPipPipNonResonantPhase",
     "The phase of the non-resonant component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_phi1NR, 0., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipK892Magnitude
    ("KmPipPipK892Magnitude",
     "The magnitude of the K*(892) component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_a1K892, 0.78, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipK892Phase
    ("KmPipPipK892Phase",
     "The phase of the K*(892) component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_phi1K892, -60., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipK1430Magnitude
    ("KmPipPipK1430Magnitude",
     "The magnitude of the K*0(1430) component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_a1K1430, 0.53, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipK1430Phase
    ("KmPipPipK1430Phase",
     "The phase of the K*0(1430) component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_phi1K1430, 132., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipK1680Magnitude
    ("KmPipPipK1680Magnitude",
     "The magnitude of the K*(1680) component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_a1K1680, 0.47, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPipK1680Phase
    ("KmPipPipK1680Phase",
     "The phase of the K*(1680) component for D+ -> K- pi+ pi+",
     &DtoKPiPiE691::_phi1K1680, - 51., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0NonResonantMagnitude
    ("KmPipPi0NonResonantMagnitude",
     "The magnitude of the non-resonant component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_a2NR, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0NonResonantPhase
    ("KmPipPi0NonResonantPhase",
     "The phase of the non-resonant component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_phi2NR, 0., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0K8920Magnitude
    ("KmPipPi0K8920Magnitude",
     "The magnitude of the K*(892)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_a2K8920, 3.19, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0K8920Phase
    ("KmPipPi0K8920Phase",
     "The phase of the K*(892)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_phi2K8920, 167., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0K892mMagnitude
    ("KmPipPi0K892mMagnitude",
     "The magnitude of the K*(892)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_a2K892m, 2.96, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0K892mPhase
    ("KmPipPi0K892mPhase",
     "The phase of the K*(892)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_phi2K892m, -112., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0RhoMagnitude
    ("KmPipPi0RhoMagnitude",
     "The magnitude of the rho component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_a2rho, 8.56, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceKmPipPi0RhoPhase
    ("KmPipPi0RhoPhase",
     "The phase of the rho component for D0 -> K- pi+ pi0",
     &DtoKPiPiE691::_phi2rho, 40., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceK0PipPimNonResonantMagnitude
    ("K0PipPimNonResonantMagnitude",
     "The magnitude of the non-resonant component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiE691::_a3NR, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceK0PipPimNonResonantPhase
    ("K0PipPimNonResonantPhase",
     "The phase of the non-resonant component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiE691::_phi3NR, 0., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceK0PipPimK892Magnitude
    ("K0PipPimK892Magnitude",
     "The magnitude of the K*(892) component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiE691::_a3K892, 2.31, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceK0PipPimK892Phase
    ("K0PipPimK892Phase",
     "The phase of the K*(892)0 component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiE691::_phi3K892, 109., -180., 180.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceK0PipPimRhoMagnitude
    ("K0PipPimRhoMagnitude",
     "The magnitude of the rho component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiE691::_a3rho, 1.59, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,double> interfaceK0PipPimRhoPhase
    ("K0PipPimRhoPhase",
     "The phase of the rho component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiE691::_phi3rho, -123., -180., 180.,
     false, false, Interface::limited);

  static Switch<DtoKPiPiE691,bool> interfaceLocalParameters
    ("LocalParameters",
     "Whether to use local values for the masses and widths or"
     " those from the ParticleData objects",
     &DtoKPiPiE691::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the ParticleData objects",
     false);

  static Parameter<DtoKPiPiE691,Energy> interfaceK8920Mass
    ("K8920Mass",
     "The mass of the K*(892)0",
     &DtoKPiPiE691::_mK8920, GeV, 0.8961 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceK8920Width
    ("K8920Width",
     "The width of the K*(892)0",
     &DtoKPiPiE691::_wK8920, GeV, 0.0505*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceK892MinusMass
    ("K892MinusMass",
     "The mass of the K*(892)-",
     &DtoKPiPiE691::_mK892m, GeV, 0.89159*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceK892MinusWidth
    ("K892MinusWidth",
     "The width of the K*(892)-",
     &DtoKPiPiE691::_wK892m, GeV, 0.0498*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceK1680Mass
    ("K1680Mass",
     "The mass of the K*(1680)",
     &DtoKPiPiE691::_mK1680, GeV, 1.714  *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceK1680Width
    ("K1680Width",
     "The width of the K*(1680)",
     &DtoKPiPiE691::_wK1680, GeV, 0.323 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceK1430Mass
    ("K1430Mass",
     "The mass of the K*(1430)",
     &DtoKPiPiE691::_mK1430, GeV, 1.429  *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceK1430Width
    ("K1430Width",
     "The width of the K*(1430)",
     &DtoKPiPiE691::_wK1430, GeV, 0.287 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceRho0Mass
    ("Rho0Mass",
     "The mass of the rho0",
     &DtoKPiPiE691::_mrho0, GeV, 0.7681 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceRho0Width
    ("Rho0Width",
     "The width of the rho0",
     &DtoKPiPiE691::_wrho0, GeV, 0.1515*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceRhoPlusMass
    ("RhoPlusMass",
     "The mass of the rho+",
     &DtoKPiPiE691::_mrhop, GeV, 0.7681 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE691,Energy> interfaceRhoPlusWidth
    ("RhoPlusWidth",
     "The width of the rho+",
     &DtoKPiPiE691::_wrhop, GeV, 0.1515*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiE691,double> interfaceMaximumWeights
    ("MaximumWeights",
     "The maximum weights for the unweighting of the decays",
     &DtoKPiPiE691::_maxwgt, -1, 1.0, 0.0, 1.0e11,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiE691,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &DtoKPiPiE691::_weights, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
}

int DtoKPiPiE691::modeNumber(bool & cc,tcPDPtr parent,
			     const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be D0 or D+
  if(abs(id0)!=ParticleID::D0&&abs(id0)!=ParticleID::Dplus) return -1;
  cc = id0<0;
  // must be three decay products
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int npip(0),npim(0),nkm(0),nk0(0),npi0(0);
  int id;
  for( ;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id          ==ParticleID::piplus)  ++npip;
    else if(id     ==ParticleID::pi0)     ++npi0;
    else if(id     ==ParticleID::piminus) ++npim;
    else if(abs(id)==ParticleID::K0)      ++nk0;
    else if(id     ==ParticleID::K_L0)    ++nk0;
    else if(id     ==ParticleID::K_S0)    ++nk0;
    else if(abs(id)==ParticleID::Kplus)   ++nkm;
  }
  if(abs(id0)==ParticleID::Dplus) {
    if((id0==ParticleID::Dplus &&(nkm==1&&npip==2))||
       (id0==ParticleID::Dminus&&(nkm==1&&npim==2))) return 0;
    else return -1;
  }
  else {
    if(npim==1&&npip==1&&nk0==1) return  2;
    else if(nkm==1&&(npip+npim)==1&&npi0==1) return 1;
    else                        return -1;
  }
}

double DtoKPiPiE691::me2(const int ichan,
			 const Particle & inpart,
			 const ParticleVector & decay,
			 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return 0.;
  }
  Complex amp;
  // D+ -> K-pi+pi+
  if(imode()==0) {
    Lorentz5Momentum pres1=decay[0]->momentum()+decay[1]->momentum();
    pres1.rescaleMass();
    double ct1 =-decayAngle(inpart.momentum(),pres1,decay[0]->momentum());
    Lorentz5Momentum pres2=decay[0]->momentum()+decay[2]->momentum();
    pres2.rescaleMass();
    double ct2 =-decayAngle(inpart.momentum(),pres2,decay[0]->momentum());
    if(ichan<0) {
      amp = _c1NR*sqrt(2.)
	+_c1K892 *amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920)
	+_c1K892 *amplitude(1,ct2,pres2.mass(),_wK8920,_mK8920)
	+_c1K1430*amplitude(0,ct1,pres1.mass(),_wK1430,_mK1430)
	+_c1K1430*amplitude(0,ct2,pres2.mass(),_wK1430,_mK1430)
	+_c1K1680*amplitude(1,ct1,pres1.mass(),_wK1680,_mK1680)
	+_c1K1680*amplitude(1,ct2,pres2.mass(),_wK1680,_mK1680);
    }
    else if(ichan==0) {
      amp=_c1K892 *amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920);
    }
    else if(ichan==1) {
      amp=_c1K892 *amplitude(1,ct2,pres2.mass(),_wK8920,_mK8920);
    }
    else if(ichan==2) {
      amp=_c1K1430*amplitude(1,ct1,pres1.mass(),_wK1430,_mK1430);
    }
    else if(ichan==3) {
      amp=_c1K1430*amplitude(1,ct2,pres2.mass(),_wK1430,_mK1430);
    }
    else if(ichan==4) {
      amp=_c1K1680*amplitude(1,ct1,pres1.mass(),_wK1680,_mK1680);
    }
    else if(ichan==5) {
      amp=_c1K1680*amplitude(1,ct2,pres2.mass(),_wK1680,_mK1680);
    }
  }
  // D0 -> K-pi+pi0
  else if(imode()==1) {
    Lorentz5Momentum pres1=decay[0]->momentum()+decay[1]->momentum();
    pres1.rescaleMass();
    double ct1 = decayAngle(inpart.momentum(),pres1,decay[0]->momentum());
    Lorentz5Momentum pres2=decay[0]->momentum()+decay[2]->momentum();
    pres2.rescaleMass();
    double ct2 = decayAngle(inpart.momentum(),pres2,decay[2]->momentum());
    Lorentz5Momentum pres3=decay[1]->momentum()+decay[2]->momentum();
    pres3.rescaleMass();
    double ct3 = decayAngle(inpart.momentum(),pres3,decay[2]->momentum());
    if(ichan<0) {
      amp = _c2NR
	+_c2K8920*amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920)
	+_c2K892m*amplitude(1,ct2,pres2.mass(),_wK892m,_mK892m)
	+_c2rho  *amplitude(1,ct3,pres3.mass(),_wrhop ,_mrhop );
    }
    else if(ichan==0) {
      amp = _c2K8920*amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920);
    }
    else if(ichan==1) {
      amp = _c2K892m*amplitude(1,ct2,pres2.mass(),_wK892m,_mK892m);
    }
    else if(ichan==2) {
      amp = _c2rho  *amplitude(1,ct3,pres3.mass(),_wrhop ,_mrhop);
    }
  }
  // D0 -> Kbar0pi+pi-
  else if(imode()==2) {
    Lorentz5Momentum pres1=decay[0]->momentum()+decay[2]->momentum();
    pres1.rescaleMass();
    double ct1 = decayAngle(inpart.momentum(),pres1,decay[0]->momentum());
    Lorentz5Momentum pres2=decay[1]->momentum()+decay[2]->momentum();
    pres2.rescaleMass();
    double ct2 = decayAngle(inpart.momentum(),pres2,decay[1]->momentum());
    if(ichan<0) {
      amp = _c3NR
	+_c3K892*amplitude(1,ct1,pres1.mass(),_wK892m,_mK892m)
	+_c3rho *amplitude(1,ct2,pres2.mass(),_wrho0 ,_mrho0 );
    }
    else if(ichan==0) {
      amp = _c3rho *amplitude(1,ct2,pres2.mass(),_wrho0 ,_mrho0 );
    }
    else if(ichan==1) {
      amp = _c3K892*amplitude(1,ct1,pres1.mass(),_wK892m,_mK892m);
    }
  }
  // now compute the matrix element
  (*ME())(0,0,0,0)=amp;
  return norm(amp);
}

void DtoKPiPiE691::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // parameters
  output << "newdef " << name() << ":KmPipPipNonResonantMagnitude " 
	 << _a1NR      << "\n";
  output << "newdef " << name() << ":KmPipPipNonResonantPhase     " 
	 << _phi1NR    << "\n";
  output << "newdef " << name() << ":KmPipPipK892Magnitude        " 
	 << _a1K892    << "\n";
  output << "newdef " << name() << ":KmPipPipK892Phase            " 
	 << _phi1K892  << "\n";
  output << "newdef " << name() << ":KmPipPipK1430Magnitude       " 
	 << _a1K1430   << "\n";
  output << "newdef " << name() << ":KmPipPipK1430Phase           " 
	 << _phi1K1430 << "\n";
  output << "newdef " << name() << ":KmPipPipK1680Magnitude       " 
	 << _a1K1680   << "\n";
  output << "newdef " << name() << ":KmPipPipK1680Phase           " 
	 << _phi1K1680 << "\n";
  output << "newdef " << name() << ":KmPipPi0NonResonantMagnitude " 
	 << _a2NR      << "\n";
  output << "newdef " << name() << ":KmPipPi0NonResonantPhase     " 
	 << _phi2NR    << "\n";
  output << "newdef " << name() << ":KmPipPi0K8920Magnitude       " 
	 << _a2K8920   << "\n";
  output << "newdef " << name() << ":KmPipPi0K8920Phase           " 
	 << _phi2K8920 << "\n";
  output << "newdef " << name() << ":KmPipPi0K892mMagnitude       " 
	 << _a2K892m   << "\n";
  output << "newdef " << name() << ":KmPipPi0K892mPhase           " 
	 << _phi2K892m << "\n";
  output << "newdef " << name() << ":KmPipPi0RhoMagnitude         " 
	 << _a2rho     << "\n";
  output << "newdef " << name() << ":KmPipPi0RhoPhase             " 
	 << _phi2rho   << "\n";
  output << "newdef " << name() << ":K0PipPimNonResonantMagnitude " 
	 << _a3NR      << "\n";
  output << "newdef " << name() << ":K0PipPimNonResonantPhase     " 
	 << _phi3NR    << "\n";
  output << "newdef " << name() << ":K0PipPimK892Magnitude        " 
	 << _a3K892    << "\n";
  output << "newdef " << name() << ":K0PipPimK892Phase            " 
	 << _phi3K892  << "\n";
  output << "newdef " << name() << ":K0PipPimRhoMagnitude         " 
	 << _a3rho     << "\n";
  output << "newdef " << name() << ":K0PipPimRhoPhase             " 
	 << _phi3rho   << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":K8920Mass      " << _mK8920/GeV << "\n";
  output << "newdef " << name() << ":K8920Width     " << _wK8920/GeV << "\n";
  output << "newdef " << name() << ":K892MinusMass  " << _mK892m/GeV << "\n";
  output << "newdef " << name() << ":K892MinusWidth " << _wK892m/GeV << "\n";
  output << "newdef " << name() << ":K1680Mass      " << _mK1680/GeV << "\n";
  output << "newdef " << name() << ":K1680Width     " << _wK1680/GeV << "\n";
  output << "newdef " << name() << ":K1430Mass      " << _mK1430/GeV << "\n";
  output << "newdef " << name() << ":K1430Width     " << _wK1430/GeV << "\n";
  output << "newdef " << name() << ":Rho0Mass       " << _mrho0 /GeV << "\n";
  output << "newdef " << name() << ":Rho0Width      " << _wrho0 /GeV << "\n";
  output << "newdef " << name() << ":RhoPlusMass    " << _mrhop /GeV << "\n";
  output << "newdef " << name() << ":RhoPlusWidth   " << _wrhop /GeV << "\n";
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    output << "insert " << name() << ":MaximumWeights " 
	   << ix << " " << _maxwgt[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_weights.size();++ix) {
    output << "insert " << name() << ":Weights " 
	   << ix << " " << _weights[ix] << "\n";
  }
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}

void DtoKPiPiE691::doinitrun() {
  DecayIntegrator::doinitrun();
  _weights.resize(mode(0)->numberChannels()+mode(1)->numberChannels()+
		  mode(2)->numberChannels());
  _maxwgt.resize(3);
  unsigned int iy=0;
  for(unsigned int ix=0;ix<3;++ix) {
    _maxwgt[ix]=mode(ix)->maxWeight();
    for(unsigned int iz=0;iz<mode(ix)->numberChannels();++iz) {
      _weights[iy]=mode(ix)->channelWeight(iz);
      ++iy;
    }
  }
}
