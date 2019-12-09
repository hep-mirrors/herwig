// -*- C++ -*-
//
// DtoKPiPiMarkIII.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiMarkIII class.
//

#include "DtoKPiPiMarkIII.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

Complex DtoKPiPiMarkIII::amplitude(bool rho, Energy mD, 
					  Energy mA , Energy mB , Energy mC ,
					  Energy mAB, Energy mAC, Energy mBC,
					  Energy mres, Energy wres) const {
  InvEnergy radius = rho ? _rrho : _rKstar;
  Energy  pAB  = Kinematics::pstarTwoBodyDecay(mAB ,mA,mB);
  Energy  pR   = Kinematics::pstarTwoBodyDecay(mres,mA,mB);
  Energy2 mgam = wres*sqr(mres)/mAB*Math::powi(pAB/pR,3)*
    (1.+sqr(radius*pR))/(1.+sqr(radius*pAB));
  Energy2 s = (sqr(mAC)-sqr(mBC)-(sqr(mD)-sqr(mC))*(sqr(mA)-sqr(mB))/sqr(mres))*
    sqrt((1.+sqr(radius*pR))/(1.+sqr(radius*pAB))); 
  Complex output=s/((sqr(mres)- sqr(mAB)-complex<Energy2>(ZERO,mgam)));
  return output;
}

DtoKPiPiMarkIII::DtoKPiPiMarkIII() {
  // Amplitudes and phases for D0 -> K- pi+ pi0
  _a1rho    = 1.0000; _phi1rho    =   0.;
  _a1Kstarm = 0.4018; _phi1Kstarm = 154.;
  _a1Kstar0 = 0.4244; _phi1Kstar0 =   7.;
  _a1NR     = 2.0693; _phi1NR     =  52.;
  // Amplitudes and phases for D0 -> Kbar0 pi+ pi-
  _a2rho    = 0.0975; _phi2rho    =  93.;
  _a2Kstar  = 0.2225; _phi2Kstar  =   0.;
  _a2NR     = 1.0000; _phi2NR     =   0.;
  // Amplitudes and phases for D+ -> Kbar0 pi+ pi0
  _a3rho    = 1.0000; _phi3rho    =   0.;
  _a3Kstar  = 0.5617; _phi3Kstar  =  43.;
  _a3NR     = 2.7250; _phi3NR     = 250.;
  // Amplitudes and phases for D+ -> K- pi+ pi+
  _a4Kstar  = 0.04749; _phi4Kstar  = 105.;
  _a4NR     = 1.00000; _phi4NR     =   0.;
  // masses and widths of the resonances
  _localparameters = true;
  _mrhop   = 0.770*GeV; _wrhop   = 0.1533*GeV;
  _mrho0   = 0.770*GeV; _wrho0   = 0.1533*GeV;
  _mKstarm = 0.8921*GeV; _wKstarm = 0.0511*GeV;
  _mKstar0 = 0.8695*GeV; _wKstar0 = 0.0502*GeV;
  // radii of the mesons
  _rrho   = 5.*mm*1e-12/hbarc; 
  _rKstar = 2.*mm*1e-12/hbarc;
  // intermediates
  generateIntermediates(true);
}


void DtoKPiPiMarkIII::persistentOutput(PersistentOStream & os) const {
  os << _a1rho << _phi1rho << _a1Kstarm << _phi1Kstarm << _a1Kstar0 << _phi1Kstar0 
     << _a1NR << _phi1NR << _c1rho << _c1Kstarm << _c1Kstar0 << _c1NR << _a2rho 
     << _phi2rho << _a2Kstar << _phi2Kstar << _a2NR << _phi2NR << _c2rho << _c2Kstar 
     << _c2NR << _a3rho << _phi3rho << _a3Kstar << _phi3Kstar << _a3NR << _phi3NR 
     << _c3rho << _c3Kstar << _c3NR << _a4Kstar << _phi4Kstar << _a4NR << _phi4NR 
     << _c4Kstar << _c4NR << _localparameters << ounit(_mrhop,GeV) << ounit(_wrhop,GeV) 
     << ounit(_mrho0,GeV) << ounit(_wrho0,GeV) << ounit(_mKstarm,GeV) 
     << ounit(_wKstarm,GeV) << ounit(_mKstar0,GeV) << ounit(_wKstar0,GeV) 
     << ounit(_rrho,1./GeV) << ounit(_rKstar,1./GeV) << _maxwgt << _weights;
}

void DtoKPiPiMarkIII::persistentInput(PersistentIStream & is, int) {
  is >> _a1rho >> _phi1rho >> _a1Kstarm >> _phi1Kstarm >> _a1Kstar0 >> _phi1Kstar0 
     >> _a1NR >> _phi1NR >> _c1rho >> _c1Kstarm >> _c1Kstar0 >> _c1NR >> _a2rho 
     >> _phi2rho >> _a2Kstar >> _phi2Kstar >> _a2NR >> _phi2NR >> _c2rho >> _c2Kstar 
     >> _c2NR >> _a3rho >> _phi3rho >> _a3Kstar >> _phi3Kstar >> _a3NR >> _phi3NR 
     >> _c3rho >> _c3Kstar >> _c3NR >> _a4Kstar >> _phi4Kstar >> _a4NR >> _phi4NR 
     >> _c4Kstar >> _c4NR >> _localparameters >> iunit(_mrhop,GeV) >> iunit(_wrhop,GeV) 
     >> iunit(_mrho0,GeV) >> iunit(_wrho0,GeV) >> iunit(_mKstarm,GeV) 
     >> iunit(_wKstarm,GeV) >> iunit(_mKstar0,GeV) >> iunit(_wKstar0,GeV) 
     >> iunit(_rrho,1./GeV) >> iunit(_rKstar,1./GeV) >> _maxwgt >> _weights;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DtoKPiPiMarkIII,DecayIntegrator>
describeHerwigDtoKPiPiMarkIII("Herwig::DtoKPiPiMarkIII", "HwSMDecay.so");

void DtoKPiPiMarkIII::Init() {

  static ClassDocumentation<DtoKPiPiMarkIII> documentation
    ("The DtoKPiPiMarkIII class performs the D -> K pi pi decays using the fit"
     "of the MarkIII collaboration",
     "The fit of the Mark III collaboration \\cite{Adler:1987sd} was used for "
     "$D\\to K\\pi\\pi$ decays",
     "\\bibitem{Adler:1987sd} J.~Adler {\\it et al.} [MARK-III Collaboration], "
     "Phys.\\ Lett.\\  B {\\bf 196} (1987) 107.");

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0RhoMagnitude
    ("KmPipPi0RhoMagnitude",
     "The magnitude of the rho component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_a1rho, 1.00, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0RhoPhase
    ("KmPipPi0RhoPhase",
     "The phase of the rho component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_phi1rho,  0.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0KstarmMagnitude
    ("KmPipPi0KstarmMagnitude",
     "The magnitude of the K*(892)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_a1Kstarm, 0.371, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0KstarmPhase
    ("KmPipPi0KstarmPhase",
     "The phase of the K*(892)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_phi1Kstarm, 154.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0Kstar0Magnitude
    ("KmPipPi0Kstar0Magnitude",
     "The magnitude of the K*(892)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_a1Kstar0, 0.391, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0Kstar0Phase
    ("KmPipPi0Kstar0Phase",
     "The phase of the K*(892)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_phi1Kstar0, 7.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0NonResonantMagnitude
    ("KmPipPi0NonResonantMagnitude",
     "The magnitude of the non-resonant component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_a1NR, 1.889, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPi0NonResonantPhase
    ("KmPipPi0NonResonantPhase",
     "The phase of the non-resonant component for D0 -> K- pi+ pi0",
     &DtoKPiPiMarkIII::_phi1NR, 52.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPimRhoMagnitude
    ("K0PipPimRhoMagnitude",
     "The magnitude of the rho component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_a2rho, 1.000, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPimRhoPhase
    ("K0PipPimRhoPhase",
     "The phase of the rho component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_phi2rho,  93.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPimKstarMagnitude
    ("K0PipPimKstarMagnitude",
     "The magnitude of the K*(892) component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_a2Kstar, 2.106, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPimKstarPhase
    ("K0PipPimKstarPhase",
     "The phase of the K*(892)0 component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_phi2Kstar, 0.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPimNonResonantMagnitude
    ("K0PipPimNonResonantMagnitude",
     "The magnitude of the non-resonant component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_a2NR, 9.379, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPimNonResonantPhase
    ("K0PipPimNonResonantPhase",
     "The phase of the non-resonant component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_phi2NR, 0.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPi0RhoMagnitude
    ("K0PipPi0RhoMagnitude",
     "The magnitude of the rho component for D+ -> Kbar0 pi+ pi0",
     &DtoKPiPiMarkIII::_a3rho, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPi0RhoPhase
    ("K0PipPi0RhoPhase",
     "The phase of the rho component for D+ -> Kbar0 pi+ pi0",
     &DtoKPiPiMarkIII::_phi3rho, 0.0, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPi0KstarMagnitude
    ("K0PipPi0KstarMagnitude",
     "The magnitude of the K* component for D+ -> Kbar0 pi+ pi0",
     &DtoKPiPiMarkIII::_a3Kstar, 0.517, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPi0KstarPhase
    ("K0PipPi0KstarPhase",
     "The phase of the K* component for D+ -> Kbar0 pi+ pi0",
     &DtoKPiPiMarkIII::_phi3Kstar, 43.0, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPi0NonResonantMagnitude
    ("K0PipPi0NonResonantMagnitude",
     "The magnitude of the non-resonant component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_a3NR, 2.490, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceK0PipPi0NonResonantPhase
    ("K0PipPi0NonResonantPhase",
     "The phase of the non-resonant component for D0 -> Kbar0 pi+pi-",
     &DtoKPiPiMarkIII::_phi3NR, 250.,0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPipNonResonantMagnitude
    ("KmPipPipNonResonantMagnitude",
     "The magnitude of the non-resonant component for D+ -> K- pi+ pi+",
     &DtoKPiPiMarkIII::_a4NR, 1.00, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPipNonResonantPhase
    ("KmPipPipNonResonantPhase",
     "The phase of the non-resonant component for D+ -> K- pi+ pi+",
     &DtoKPiPiMarkIII::_phi4NR, 0., 0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPipKstarMagnitude
    ("KmPipPipKstarMagnitude",
     "The magnitude of the K*(892) component for D+ -> K- pi+ pi+",
     &DtoKPiPiMarkIII::_a4Kstar, 0.049, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,double> interfaceKmPipPipKstarPhase
    ("KmPipPipKstarPhase",
     "The phase of the K*(892) component for D+ -> K- pi+ pi+",
     &DtoKPiPiMarkIII::_phi4Kstar, 105.,0., 360.,
     false, false, Interface::limited);

  static Switch<DtoKPiPiMarkIII,bool> interfaceLocalParameters
    ("LocalParameters",
     "Whether to use local values for the masses and widths or"
     " those from the ParticleData objects",
     &DtoKPiPiMarkIII::_localparameters, true, false, false);
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

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceKstar0Mass
    ("Kstar0Mass",
     "The mass of the K*(892)0",
     &DtoKPiPiMarkIII::_mKstar0, GeV, 0.8965 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceKstar0Width
    ("Kstar0Width",
     "The width of the K*(892)0",
     &DtoKPiPiMarkIII::_wKstar0, GeV, 0.0502*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceKstarMinusMass
    ("KstarMinusMass",
     "The mass of the K*(892)-",
     &DtoKPiPiMarkIII::_mKstarm, GeV, 0.8921*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceKstarMinusWidth
    ("KstarMinusWidth",
     "The width of the K*(892)-",
     &DtoKPiPiMarkIII::_wKstarm, GeV, 0.0511*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceRho0Mass
    ("Rho0Mass",
     "The mass of the rho0",
     &DtoKPiPiMarkIII::_mrho0, GeV, 0.770 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceRho0Width
    ("Rho0Width",
     "The width of the rho0",
     &DtoKPiPiMarkIII::_wrho0, GeV, 0.1533*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceRhoPlusMass
    ("RhoPlusMass",
     "The mass of the rho+",
     &DtoKPiPiMarkIII::_mrhop, GeV, 0.770 *GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,Energy> interfaceRhoPlusWidth
    ("RhoPlusWidth",
     "The width of the rho+",
     &DtoKPiPiMarkIII::_wrhop, GeV, 0.1533*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,InvEnergy> interfaceRhoRadius
    ("RhoRadius",
     "The radius of the rho for the Blatt-Weisskopf factor",
     &DtoKPiPiMarkIII::_rrho, mm*1e-12/hbarc, 5.0*mm*1e-12/hbarc, 0.0*mm*1e-12/hbarc, 10.0*mm*1e-12/hbarc,
     true, false, Interface::limited);

  static Parameter<DtoKPiPiMarkIII,InvEnergy> interfaceKstarRadius
    ("KstarRadius",
     "The radius of the K* for the Blatt-Weisskopf factor",
     &DtoKPiPiMarkIII::_rKstar, mm*1e-12/hbarc, 2.0*mm*1e-12/hbarc, 0.0*mm*1e-12/hbarc, 10.0*mm*1e-12/hbarc,
     true, false, Interface::limited);

  static ParVector<DtoKPiPiMarkIII,double> interfaceMaximumWeights
    ("MaximumWeights",
     "The maximum weights for the unweighting of the decays",
     &DtoKPiPiMarkIII::_maxwgt, -1, 1.0, 0.0, 1.0e11,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiMarkIII,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &DtoKPiPiMarkIII::_weights, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
}

void DtoKPiPiMarkIII::doinit() {
  DecayIntegrator::doinit();
  double fact = Constants::pi/180.;
  // amplitudes for D0 -> K- pi+ pi0
  _c1rho    = _a1rho   *Complex(cos(_phi1rho   *fact),sin(_phi1rho   *fact));
  _c1Kstarm = _a1Kstarm*Complex(cos(_phi1Kstarm*fact),sin(_phi1Kstarm*fact));
  _c1Kstar0 = _a1Kstar0*Complex(cos(_phi1Kstar0*fact),sin(_phi1Kstar0*fact));
  _c1NR     = _a1NR    *Complex(cos(_phi1NR    *fact),sin(_phi1NR    *fact));
  // amplitudes for D0 -> Kbar0 pi+ pi-
  _c2rho    = _a2rho   *Complex(cos(_phi2rho   *fact),sin(_phi2rho   *fact));
  _c2Kstar  = _a2Kstar *Complex(cos(_phi2Kstar *fact),sin(_phi2Kstar *fact));
  _c2NR     = _a2NR    *Complex(cos(_phi2NR    *fact),sin(_phi2NR    *fact));
  // amplitudes for D+ -> Kbar0 pi+ pi0
  _c3rho    = _a3rho   *Complex(cos(_phi3rho   *fact),sin(_phi3rho   *fact));
  _c3Kstar  = _a3Kstar *Complex(cos(_phi3Kstar *fact),sin(_phi3Kstar *fact));
  _c3NR     = _a3NR    *Complex(cos(_phi3NR    *fact),sin(_phi3NR    *fact));
  // amplitudes for D+ -> K- pi+ pi+
  _c4Kstar  = _a4Kstar *Complex(cos(_phi4Kstar *fact),sin(_phi4Kstar *fact));
  _c4NR     = _a4NR    *Complex(cos(_phi4NR    *fact),sin(_phi4NR    *fact));
  // intermediate resonances
  tPDPtr k8920 = getParticleData(ParticleID::Kstarbar0 );
  tPDPtr k892m = getParticleData(ParticleID::Kstarminus);
  tPDPtr rho0  = getParticleData(ParticleID::rho0      );
  tPDPtr rhop  = getParticleData(ParticleID::rhoplus   );
  // D0 -> K- pi+ pi0
  tPDVector extpart(4);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::pi0);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  unsigned int ix=0;
  if(rhop) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rhop,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k8920,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
  }
  // add the mode
  vector<double> wtemp;
  if(ix<=_weights.size()) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit,wit+ix);
  }
  else {
    wtemp=vector<double>(ix,1./double(ix));
  }
  if(_maxwgt.empty()) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[0],wtemp);
  // D0 -> Kbar0 pi+ pi-
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kbar0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piminus);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  unsigned int iy=ix;
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
  if(iy<=_weights.size()) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit+ix,wit+iy);
  }
  else {
    wtemp=vector<double>(iy-ix,1./double(iy-ix));
  }
  if(_maxwgt.size()<2) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[1],wtemp);
  // D+ -> Kbar0pi+pi0
  ix=iy;
  extpart[0]=getParticleData(ParticleID::Dplus);
  extpart[1]=getParticleData(ParticleID::Kbar0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::pi0);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  if(rhop) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rhop,0,0., 2,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k8920,0,0., 1,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  // add the mode
  if(iy<=_weights.size()) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit+ix,wit+iy);
  }
  else {
    wtemp=vector<double>(iy-ix,1./double(iy-ix));
  }
  if(_maxwgt.size()<3) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[2],wtemp);
  // D+ -> K-pi+pi+
  ix=iy;
  extpart[0]=getParticleData(ParticleID::Dplus);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piplus);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k8920,0,0., 1,2);
    mode->addChannel(newchannel);
    ++iy;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k8920,0,0., 1,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  // add the mode
  if(iy<=_weights.size()) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit+ix,wit+iy);
  }
  else {
    wtemp=vector<double>(iy-ix,1./double(iy-ix));
  }
  if(_maxwgt.size()<4) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[3],wtemp);
  // reset the resonance parameters in the integration if needed
  if(_localparameters) {
    resetIntermediate(k8920,_mKstar0,_wKstar0);
    resetIntermediate(k892m,_mKstarm,_wKstarm);
    resetIntermediate(rho0 ,_mrho0 ,_wrho0 );
    resetIntermediate(rhop ,_mrhop ,_wrhop );
  }
  else {
    _mKstar0 = k8920->mass();
    _mKstarm = k892m->mass();
    _mrho0   = rho0 ->mass();
    _mrhop   = rhop ->mass();
    _wKstar0 = k8920->width();
    _wKstarm = k892m->width();
    _wrho0   = rho0 ->width();
    _wrhop   = rhop ->width();
  }
}

int DtoKPiPiMarkIII::modeNumber(bool & cc, tcPDPtr parent, 
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
       (id0==ParticleID::Dminus&&(nkm==1&&npim==2))) 
      return 3;
    else if((id0==ParticleID::Dplus &&nk0==1&&npip==1&&npi0==1)||
	    (id0==ParticleID::Dminus&&nk0==1&&npim==1&&npi0==1))
      return 2;
    else 
      return -1;
  }
  else {
    if(npim==1&&npip==1&&nk0==1) return  1;
    else if(nkm==1&&(npip+npim)==1&&npi0==1) return 0;
    else                        return -1;
  }
}

double DtoKPiPiMarkIII::me2( const int ichan,const Particle & inpart,
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
  // compute the invariant masses needed to calulate the amplitudes
  Energy mA  = decay[0]->mass();
  Energy mB  = decay[1]->mass();
  Energy mC  = decay[2]->mass();
  Energy mD  = inpart.mass();
  Energy mAB = (decay[0]->momentum()+decay[1]->momentum()).m();
  Energy mAC = (decay[0]->momentum()+decay[2]->momentum()).m();
  Energy mBC = (decay[1]->momentum()+decay[2]->momentum()).m();
  // compute the amplitudes for the resonaces present in both models
  Complex amp;
  // D0 -> K- pi+ pi0
  if(imode()==0) {
    if(ichan<0) {
      amp = _c1NR
	-_c1rho   *amplitude(true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrhop  ,_wrhop  )
 	-_c1Kstarm*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstarm,_wKstarm)
 	-_c1Kstar0*amplitude(false,mD,mA,mB,mC,mAB,mAC,mBC,_mKstar0,_wKstar0);
    }
    else if(ichan==0) {
      amp = _c1rho   *amplitude(true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrhop  ,_wrhop  );
    }
    else if(ichan==1) {
      amp = _c1Kstarm*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstarm,_wKstarm);
    }
    else {
      amp = _c1Kstar0*amplitude(false,mD,mA,mB,mC,mAB,mAC,mBC,_mKstar0,_wKstar0);
    }
  }
  // D0 -> Kbar0 pi+pi-
  else if(imode()==1) {
    if(ichan<0) {
      amp = _c2NR*(-1.+2.*UseRandom::rndbool())
	-_c2rho  *amplitude(true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrho0  ,_wrho0  )
 	+_c2Kstar*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstarm,_wKstarm);
    }
    else if(ichan==0) {
      amp =_c2rho  *amplitude(true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrho0  ,_wrho0  );
    }
    else {
      amp =_c2Kstar*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstarm,_wKstarm);
    }
  }
  // D+ -> Kbar0 pi+ pi0
  else if(imode()==2) {
    if(ichan<0) {
      amp = _c3NR
 	-_c3rho  *amplitude(true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrhop  ,_wrhop  )
  	-_c3Kstar*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstar0,_wKstar0);
    }
    else if(ichan==0) {
      amp = _c3rho  *amplitude(true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrhop,_wrhop);
    }
    else {
      amp = _c3Kstar*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstar0,_wKstar0);
    }
  }
  // D+ -> K- pi+ pi+
  else {
    if(ichan<0) {
      amp = _c4NR
	+_c4Kstar*amplitude(false,mD,mA,mB,mC,mAB,mAC,mBC,_mKstar0,_wKstar0)
 	+_c4Kstar*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstar0,_wKstar0);
    }
    else if(ichan==0) {
      amp =_c4Kstar*amplitude(false,mD,mA,mB,mC,mAB,mAC,mBC,_mKstar0,_wKstar0);
    }
    else {
      amp =_c4Kstar*amplitude(false,mD,mA,mC,mB,mAC,mAB,mBC,_mKstar0,_wKstar0);
    }
  }
  // now compute the matrix element
  (*ME())(0,0,0,0)=amp;
  return norm(amp);
}

void DtoKPiPiMarkIII::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // parameters
  output << "newdef " << name() << ":KmPipPi0RhoMagnitude " << _a1rho << "\n";
  output << "newdef " << name() << ":KmPipPi0RhoPhase " << _phi1rho << "\n";
  output << "newdef " << name() << ":KmPipPi0KstarmMagnitude " << _a1Kstarm << "\n";
  output << "newdef " << name() << ":KmPipPi0KstarmPhase " << _phi1Kstarm << "\n";
  output << "newdef " << name() << ":KmPipPi0Kstar0Magnitude " << _a1Kstar0 << "\n";
  output << "newdef " << name() << ":KmPipPi0Kstar0Phase " << _phi1Kstar0 << "\n";
  output << "newdef " << name() << ":KmPipPi0NonResonantMagnitude " << _a1NR << "\n";
  output << "newdef " << name() << ":KmPipPi0NonResonantPhase " << _phi1NR << "\n";
  output << "newdef " << name() << ":K0PipPimRhoMagnitude " << _a2rho << "\n";
  output << "newdef " << name() << ":K0PipPimRhoPhase " << _phi2rho << "\n";
  output << "newdef " << name() << ":K0PipPimKstarMagnitude " << _a2Kstar << "\n";
  output << "newdef " << name() << ":K0PipPimKstarPhase " << _phi2Kstar << "\n";
  output << "newdef " << name() << ":K0PipPimNonResonantMagnitude " << _a2NR << "\n";
  output << "newdef " << name() << ":K0PipPimNonResonantPhase " << _phi2NR << "\n";
  output << "newdef " << name() << ":K0PipPi0RhoMagnitude " << _a3rho << "\n";
  output << "newdef " << name() << ":K0PipPi0RhoPhase " << _phi3rho << "\n";
  output << "newdef " << name() << ":K0PipPi0KstarMagnitude " << _a3Kstar << "\n";
  output << "newdef " << name() << ":K0PipPi0KstarPhase " << _phi3Kstar << "\n";
  output << "newdef " << name() << ":K0PipPi0NonResonantMagnitude " << _a3NR << "\n";
  output << "newdef " << name() << ":K0PipPi0NonResonantPhase " << _phi3NR << "\n";
  output << "newdef " << name() << ":KmPipPipNonResonantMagnitude " << _a4NR << "\n";
  output << "newdef " << name() << ":KmPipPipNonResonantPhase " << _phi4NR << "\n";
  output << "newdef " << name() << ":KmPipPipKstarMagnitude " << _a4Kstar << "\n";
  output << "newdef " << name() << ":KmPipPipKstarPhase " << _phi4Kstar << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Kstar0Mass " << _mKstar0/GeV << "\n";
  output << "newdef " << name() << ":Kstar0Width " << _wKstar0/GeV << "\n";
  output << "newdef " << name() << ":KstarMinusMass " << _mKstarm/GeV << "\n";
  output << "newdef " << name() << ":KstarMinusWidth " << _wKstarm/GeV << "\n";
  output << "newdef " << name() << ":Rho0Mass " << _mrho0/GeV << "\n";
  output << "newdef " << name() << ":Rho0Width " << _wrho0/GeV << "\n";
  output << "newdef " << name() << ":RhoPlusMass " << _mrhop/GeV << "\n";
  output << "newdef " << name() << ":RhoPlusWidth " << _wrhop/GeV << "\n";
  output << "newdef " << name() << ":RhoRadius " << _rrho/(mm*1e-12/hbarc) << "\n";
  output << "newdef " << name() << ":KstarRadius " << _rKstar/(mm*1e-12/hbarc) << "\n";
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

void DtoKPiPiMarkIII::doinitrun() {
  DecayIntegrator::doinitrun();
  _weights.resize(mode(0)->numberChannels()+mode(1)->numberChannels()+
		  mode(2)->numberChannels()+mode(3)->numberChannels());
  _maxwgt.resize(4);
  unsigned int iy=0;
  for(unsigned int ix=0;ix<4;++ix) {
    _maxwgt[ix]=mode(ix)->maxWeight();
    for(unsigned int iz=0;iz<mode(ix)->numberChannels();++iz) {
      _weights[iy]=mode(ix)->channelWeight(iz);
      ++iy;
    }
  }
}

