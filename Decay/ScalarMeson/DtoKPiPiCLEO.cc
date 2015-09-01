// -*- C++ -*-
//
// DtoKPiPiCLEO.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiCLEO class.
//

#include "DtoKPiPiCLEO.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

DtoKPiPiCLEO::DtoKPiPiCLEO() : _c1NR(), _c1rho(), _c1Kstarm(), _c1Kstar0(), 
			       _c1K1430m(), _c1K14300(), _c1rho1700(), _c1K1680(), 
			       _c2Kstarp(), _c2rho(), _c2omega(), _c2Kstarm(),
			       _c2f980(), _c2f2(), _c2f1370(), _c2K14300(),
			       _c2K14302(), _c2K1680(), _c2NR(), _rD0(), _rres() {
  // use local values for masses and widths
  _localparameters=true;
  // masses and widths
  _momega   =  782.57*MeV; _womega   =   8.44*MeV;
  _mf2      = 1275.4 *MeV; _wf2      = 185.1 *MeV;
  _mf1370   = 1310   *MeV; _wf1370   = 272.0 *MeV;
  _mK14300  = 1412   *MeV; _wK14300  = 294   *MeV;
  _mK14302  = 1425.6 *MeV; _wK14302  =  98.5 *MeV;
  _mK1680   = 1717   *MeV; _wK1680   = 322   *MeV;
  _mrho1700 = 1700   *MeV; _wrho1700 = 240   *MeV;
  _mK8920   =  896.1 *MeV; _wK8920   =  50.5 *MeV;
  _mK892A   =  891.5 *MeV; _wK892A   =  50   *MeV;
  _mK892B   =  891.66*MeV; _wK892B   =  50.8 *MeV;
  _mrhoA    =  770   *MeV; _wrhoA    = 150.7 *MeV;
  _mrhoB    =  769.3 *MeV; _wrhoB    = 150.2 *MeV;
  _mf980    =  977.00*MeV; _wf980    =  50.  *MeV; 
  _f0opt=false; _gpi=0.09; _gK=0.02; 
  // amplitudes and phases for D0 -> K-pi+pi0
  _a1NR      = 1.75     ; _phi1NR      =  31.2;
  _a1rho     = 1.00     ; _phi1rho     =   0. ;
  _a1Kstarm  = 0.44     ; _phi1Kstarm  = 163  ;
  _a1Kstar0  = 0.39     ; _phi1Kstar0  =  -0.2;
  _a1K1430m  = 0.77*GeV2; _phi1K1430m  =  55.5;
  _a1K14300  = 0.85*GeV2; _phi1K14300  = 166  ;
  _a1rho1700 = 2.50     ; _phi1rho1700 = 171  ;
  _a1K1680   = 2.50     ; _phi1K1680   = 103  ;
  // amplitudes and phases for D0 -> K0pi+pi-
  _a2Kstarp   = 0.11     ; _phi2Kstarp = 321;
  _a2rho      = 1.00     ; _phi2rho    =   0;
  _a2omega    = 0.037    ; _phi2omega  = 114;
  _a2Kstarm   = 1.56     ; _phi2Kstarm = 150;
  _a2f980     = 0.34*GeV2; _phi2f980   = 188;
  _a2f2       = 0.7/GeV2 ; _phi2f2     = 308;
  _a2f1370    = 1.8*GeV2 ; _phi2f1370  =  85;
  _a2K14300   = 2.0*GeV2 ; _phi2K14300 =   3;
  _a2K14302   = 1.0/GeV2 ; _phi2K14302 = 335;
  _a2K1680    = 5.6      ; _phi2K1680  = 174;
  _a2NR       = 1.1      ; _phi2NR     = 340;
  // radial sizes
  _rD0  = 5.0/GeV;
  _rres = 1.5/GeV;
  // zero masses
  _mpi=ZERO;
  _mkp=ZERO;
  _mk0=ZERO;
  // intermediates
  generateIntermediates(true);
}

void DtoKPiPiCLEO::doinit() {
  DecayIntegrator::doinit();
  // complex amplitudes for K-pi+pi0
  double fact = Constants::pi/180.;
  _c1NR      = _a1NR     *Complex(cos(_phi1NR     *fact),sin(_phi1NR     *fact));
  _c1rho     = _a1rho    *Complex(cos(_phi1rho    *fact),sin(_phi1rho    *fact));
  _c1Kstarm  = _a1Kstarm *Complex(cos(_phi1Kstarm *fact),sin(_phi1Kstarm *fact));
  _c1Kstar0  = _a1Kstar0 *Complex(cos(_phi1Kstar0 *fact),sin(_phi1Kstar0 *fact));
  _c1K1430m  = _a1K1430m *Complex(cos(_phi1K1430m *fact),sin(_phi1K1430m *fact));
  _c1K14300  = _a1K14300 *Complex(cos(_phi1K14300 *fact),sin(_phi1K14300 *fact));
  _c1rho1700 = _a1rho1700*Complex(cos(_phi1rho1700*fact),sin(_phi1rho1700*fact));
  _c1K1680   = _a1K1680  *Complex(cos(_phi1K1680  *fact),sin(_phi1K1680  *fact));
  // complex amplitudes for D0 -> K0pi+pi-
  _c2Kstarp  = _a2Kstarp*Complex(cos(_phi2Kstarp*fact),sin(_phi2Kstarp*fact));
  _c2rho     = _a2rho   *Complex(cos(_phi2rho   *fact),sin(_phi2rho   *fact));
  _c2omega   = _a2omega *Complex(cos(_phi2omega *fact),sin(_phi2omega *fact));
  _c2Kstarm  = _a2Kstarm*Complex(cos(_phi2Kstarm*fact),sin(_phi2Kstarm*fact));
  _c2f980    = _a2f980  *Complex(cos(_phi2f980  *fact),sin(_phi2f980  *fact));
  _c2f2      = _a2f2    *Complex(cos(_phi2f2    *fact),sin(_phi2f2    *fact));
  _c2f1370   = _a2f1370 *Complex(cos(_phi2f1370 *fact),sin(_phi2f1370 *fact));
  _c2K14300  = _a2K14300*Complex(cos(_phi2K14300*fact),sin(_phi2K14300*fact));
  _c2K14302  = _a2K14302*Complex(cos(_phi2K14302*fact),sin(_phi2K14302*fact));
  _c2K1680   = _a2K1680 *Complex(cos(_phi2K1680 *fact),sin(_phi2K1680 *fact));
  _c2NR      = _a2NR    *Complex(cos(_phi2NR    *fact),sin(_phi2NR    *fact));
  // pion and kaon masses
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mkp = getParticleData(ParticleID::Kplus )->mass();
  _mk0 = getParticleData(ParticleID::K0    )->mass();
  // resonances for the channels 
  tPDPtr k892m   = getParticleData(ParticleID::Kstarminus);
  tPDPtr k892p   = getParticleData(ParticleID::Kstarplus);
  tPDPtr k8920   = getParticleData(ParticleID::Kstarbar0);
  tPDPtr rho770  = getParticleData(ParticleID::rhoplus);
  tPDPtr rho0    = getParticleData(ParticleID::rho0);
  tPDPtr rho1700 = getParticleData(30213);
  tPDPtr k1680m  = getParticleData(-30323);
  tPDPtr k1430m0 = getParticleData(ParticleID::Kstar_0minus);
  tPDPtr k1430m2 = getParticleData(ParticleID::Kstar_2minus);
  tPDPtr k143000 = getParticleData(ParticleID::Kstar_0bar0 );
  tPDPtr omega   = getParticleData(ParticleID::omega);
  tPDPtr f980    = getParticleData(9010221);
  tPDPtr f1370   = getParticleData(10221);
  tPDPtr f2      = getParticleData(ParticleID::f_2);
  DecayPhaseSpaceChannelPtr newchannel;
  // D0 -> K- pi+ pi0
  tPDVector extpart(4);
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::pi0);
  DecayPhaseSpaceModePtr mode1 = new_ptr(DecayPhaseSpaceMode(extpart,this));
  int ix=0;
  if(rho770) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode1));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho770,0,0., 2,3);
    mode1->addChannel(newchannel);
    ++ix;
  }
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode1));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode1->addChannel(newchannel);
    ++ix;
  }
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode1));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k8920,0,0., 1,2);
    mode1->addChannel(newchannel);
    ++ix;
  }
  if(k1430m0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode1));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1430m0,0,0., 1,3);
    mode1->addChannel(newchannel);
    ++ix;
  }
  if(k143000) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode1));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k143000,0,0., 1,2);
    mode1->addChannel(newchannel);
    ++ix;
  }
  if(rho1700) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode1));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho1700,0,0., 2,3);
    mode1->addChannel(newchannel);
    ++ix;
  }
  if(k1680m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode1));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1680m,0,0., 1,3);
    mode1->addChannel(newchannel);
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
  addMode(mode1,_maxwgt[0],wtemp);
  // D0 -> Kbar0 pi+ pi-
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kbar0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piminus);
  DecayPhaseSpaceModePtr mode2 = new_ptr(DecayPhaseSpaceMode(extpart,this));
  int iy=ix;
  if(k892p) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k892p,0,0., 1,2);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(rho0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho0,0,0., 2,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(omega) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(omega,0,0., 2,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(f980) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(f980,0,0., 2,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(f2) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(f2,0,0., 2,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(f1370) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(f1370,0,0., 2,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(k1430m0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1430m0,0,0., 1,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(k1430m2) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1430m2,0,0., 1,3);
    mode2->addChannel(newchannel);
    ++iy;
  }
  if(k1680m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode2));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1680m,0,0., 1,3);
    mode2->addChannel(newchannel);
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
  addMode(mode2,_maxwgt[1],wtemp);
  if(!_localparameters) {
    _momega   = omega  ->mass();
    _mf980    = f980   ->mass();
    _mf2      = f2     ->mass();
    _mf1370   = f1370  ->mass();
    _mK14300  = k1430m0->mass();
    _mK14302  = k1430m2->mass();
    _mK1680   = k1680m ->mass();
    _mrho1700 = rho1700->mass();
    _mK8920   = k8920  ->mass();
    _mK892A   = k892p  ->mass();
    _mK892B   = k892p  ->mass();
    _mrhoA    = rho770 ->mass();
    _mrhoB    = rho0   ->mass();
    _womega   = omega  ->width();
    _wf980    = f980   ->width();
    _wf2      = f2     ->width();
    _wf1370   = f1370  ->width();
    _wK14300  = k1430m0->width();
    _wK14302  = k1430m2->width();
    _wK1680   = k1680m ->width();
    _wrho1700 = rho1700->width();
    _wK8920   = k8920  ->width();
    _wK892A   = k892p  ->width();
    _wK892B   = k892p  ->width();
    _wrhoA    = rho770 ->width();
    _wrhoB    = rho0   ->width();
  }
  else {
    mode1->resetIntermediate(rho770 ,_mrhoA   ,_wrhoA   );
    mode1->resetIntermediate(k892m  ,_mK892A  ,_wK892A  );
    mode1->resetIntermediate(k8920  ,_mK8920  ,_wK8920  );
    mode1->resetIntermediate(k1430m0,_mK14300 ,_wK14300 );
    mode1->resetIntermediate(k143000,_mK14300 ,_wK14300 );
    mode1->resetIntermediate(rho1700,_mrho1700,_wrho1700);
    mode1->resetIntermediate(k1680m ,_mK1680  ,_wK1680  );
    mode2->resetIntermediate(k892p  ,_mK892B  ,_wK892B  );
    mode2->resetIntermediate(rho0   ,_mrhoB   ,_wrhoB   );
    mode2->resetIntermediate(omega  ,_momega  ,_womega  );
    mode2->resetIntermediate(k892m  ,_mK892B  ,_wK892B  );
    mode2->resetIntermediate(f980   ,_mf980   ,_wf980   );
    mode2->resetIntermediate(f2     ,_mf2     ,_wf2     );
    mode2->resetIntermediate(f1370  ,_mf1370  ,_wf1370  );
    mode2->resetIntermediate(k1430m0,_mK14300 ,_wK14300 );
    mode2->resetIntermediate(k1430m2,_mK14302 ,_wK14302 );
    mode2->resetIntermediate(k1680m ,_mK1680  ,_wK1680  );
  }
}

void DtoKPiPiCLEO::persistentOutput(PersistentOStream & os) const {
  os << ounit(_momega,GeV) << ounit(_womega,GeV) << ounit(_mf980,GeV) 
     << _gpi << _gK  << ounit(_mf2,GeV) << ounit(_wf2,GeV) << ounit(_mf1370,GeV) 
     << ounit(_wf1370,GeV) << ounit(_mK14300,GeV) << ounit(_wK14300,GeV) 
     << ounit(_mK14302,GeV) << ounit(_wK14302,GeV) << ounit(_mK1680,GeV) 
     << ounit(_wK1680,GeV) << ounit(_mrho1700,GeV) << ounit(_wrho1700,GeV) 
     << ounit(_mK8920,GeV) << ounit(_wK8920,GeV) << ounit(_mK892A,GeV) 
     << ounit(_wK892A,GeV) << ounit(_mK892B,GeV) << ounit(_wK892B,GeV) 
     << ounit(_mrhoA,GeV) << ounit(_wrhoA,GeV) << ounit(_mrhoB,GeV) 
     << ounit(_wrhoB,GeV) << _a1NR << _phi1NR << _a1rho << _phi1rho 
     << _a1Kstarm << _phi1Kstarm << _a1Kstar0 << _phi1Kstar0 << ounit(_a1K1430m,GeV2) 
     << _phi1K1430m << ounit(_a1K14300,GeV2) << _phi1K14300 << _a1rho1700 
     << _phi1rho1700 << _a1K1680 << _phi1K1680 << _c1NR << _c1rho << _c1Kstarm 
     << _c1Kstar0 << ounit(_c1K1430m,GeV2) << ounit(_c1K14300,GeV2) 
     << _c1rho1700 << _c1K1680 << _a2Kstarp << _phi2Kstarp << _a2rho 
     << _phi2rho << _a2omega << _phi2omega << _a2Kstarm << _phi2Kstarm 
     << ounit(_a2f980,GeV2) << _phi2f980 << ounit(_a2f2,1./GeV2) << _phi2f2 
     << ounit(_a2f1370,GeV2) << _phi2f1370  << ounit(_a2K14300,GeV2) << _phi2K14300 
     << ounit(_a2K14302,1./GeV2) << _phi2K14302 << _a2K1680 << _phi2K1680 << _a2NR 
     << _phi2NR << _c2Kstarp << _c2rho << _c2omega << _c2Kstarm << ounit(_c2f980,GeV2)
     << ounit(_c2f2,1./GeV2) << ounit(_c2f1370,GeV2) << ounit(_c2K14300,GeV2) 
     << ounit(_c2K14302,1./GeV2) << _c2K1680 << _c2NR << _maxwgt << _weights 
     << ounit(_rD0,1./GeV) << ounit(_rres,1./GeV) << ounit(_mpi,GeV) 
     << ounit(_mkp,GeV) << ounit(_mk0,GeV) << ounit(_wf980,GeV) << _f0opt 
     << _localparameters;
}

void DtoKPiPiCLEO::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_momega,GeV) >> iunit(_womega,GeV) >> iunit(_mf980,GeV) 
     >> _gpi >> _gK  >> iunit(_mf2,GeV) >> iunit(_wf2,GeV) >> iunit(_mf1370,GeV) 
     >> iunit(_wf1370,GeV) >> iunit(_mK14300,GeV) >> iunit(_wK14300,GeV) 
     >> iunit(_mK14302,GeV) >> iunit(_wK14302,GeV) >> iunit(_mK1680,GeV) 
     >> iunit(_wK1680,GeV) >> iunit(_mrho1700,GeV) >> iunit(_wrho1700,GeV) 
     >> iunit(_mK8920,GeV) >> iunit(_wK8920,GeV) >> iunit(_mK892A,GeV) 
     >> iunit(_wK892A,GeV) >> iunit(_mK892B,GeV) >> iunit(_wK892B,GeV) 
     >> iunit(_mrhoA,GeV) >> iunit(_wrhoA,GeV) >> iunit(_mrhoB,GeV) 
     >> iunit(_wrhoB,GeV) >> _a1NR >> _phi1NR >> _a1rho >> _phi1rho 
     >> _a1Kstarm >> _phi1Kstarm >> _a1Kstar0 >> _phi1Kstar0 >> iunit(_a1K1430m,GeV2) 
     >> _phi1K1430m >> iunit(_a1K14300,GeV2) >> _phi1K14300 >> _a1rho1700 
     >> _phi1rho1700 >> _a1K1680 >> _phi1K1680 >> _c1NR >> _c1rho >> _c1Kstarm 
     >> _c1Kstar0 >> iunit(_c1K1430m,GeV2) >> iunit(_c1K14300,GeV2) 
     >> _c1rho1700 >> _c1K1680 >> _a2Kstarp >> _phi2Kstarp >> _a2rho 
     >> _phi2rho >> _a2omega >> _phi2omega >> _a2Kstarm >> _phi2Kstarm 
     >> iunit(_a2f980,GeV2) >> _phi2f980 >> iunit(_a2f2,1./GeV2) >> _phi2f2 
     >> iunit(_a2f1370,GeV2) >> _phi2f1370  >> iunit(_a2K14300,GeV2) >> _phi2K14300 
     >> iunit(_a2K14302,1./GeV2) >> _phi2K14302 >> _a2K1680 >> _phi2K1680 >> _a2NR 
     >> _phi2NR >> _c2Kstarp >> _c2rho >> _c2omega >> _c2Kstarm >> iunit(_c2f980,GeV2)
     >> iunit(_c2f2,1./GeV2) >> iunit(_c2f1370,GeV2) >> iunit(_c2K14300,GeV2) 
     >> iunit(_c2K14302,1./GeV2) >> _c2K1680 >> _c2NR >> _maxwgt >> _weights 
     >> iunit(_rD0,1./GeV) >> iunit(_rres,1./GeV) >> iunit(_mpi,GeV) 
     >> iunit(_mkp,GeV) >> iunit(_mk0,GeV) >> iunit(_wf980,GeV) >> _f0opt 
     >> _localparameters;
}

ClassDescription<DtoKPiPiCLEO> DtoKPiPiCLEO::initDtoKPiPiCLEO;
// Definition of the static class description member.

void DtoKPiPiCLEO::Init() {

  static ClassDocumentation<DtoKPiPiCLEO> documentation
    ("The DtoKPiPiCLEO class implements the models of CLEO for"
     " D0 -> Kbar0 pi+pi- and D0 -> K- pi+ pi0",
     "The CLEO fits of \\cite{Muramatsu:2002jp} and \\cite{Kopp:2000gv} were"
     " used for the decays $D^0\\to\\bar{K}^0\\pi^+\\pi^-$ and"
     " $D^0\\to K^-\\pi^+\\pi^0$.",
     "\\bibitem{Muramatsu:2002jp} H.~Muramatsu {\\it et al.}  "
     "[CLEO Collaboration],Phys.\\ Rev.\\ Lett.\\  {\\bf 89} (2002) 251802"
     "[Erratum-ibid.\\  {\\bf 90} (2003) 059901] [arXiv:hep-ex/0207067].\n"
     "\\bibitem{Kopp:2000gv} S.~Kopp {\\it et al.}  [CLEO Collaboration], "
     "Phys.\\ Rev.\\  D {\\bf 63} (2001) 092001 [arXiv:hep-ex/0011065]."
     );

  static Switch<DtoKPiPiCLEO,bool> interfaceLocalParameters
    ("LocalParameters",
     "Whether to use local values for the masses and widths or"
     " those from the ParticleData objects",
     &DtoKPiPiCLEO::_localparameters, true, false, false);
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

  static Parameter<DtoKPiPiCLEO,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &DtoKPiPiCLEO::_momega, MeV, 782.57*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacef980Mass
    ("f980Mass",
     "The mass of the f_0(980) meson",
     &DtoKPiPiCLEO::_mf980, MeV, 977.00*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &DtoKPiPiCLEO::_mf2, MeV, 1275.4 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacef1370Mass
    ("f1370Mass",
     "The mass of the f_0(1370) meson",
     &DtoKPiPiCLEO::_mf1370, MeV, 1310   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceK_01430Mass
    ("K_01430Mass",
     "The mass of the K_0(1430) meson",
     &DtoKPiPiCLEO::_mK14300, MeV, 1412   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceK_21430Mass
    ("K_21430Mass",
     "The mass of the K_2(1430) meson",
     &DtoKPiPiCLEO::_mK14302, MeV, 1425.6 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstar1680Mass
    ("Kstar1680Mass",
     "The mass of the K*(1680) meson",
     &DtoKPiPiCLEO::_mK1680, MeV, 1717   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacerho1700Mass
    ("rho1700Mass",
     "The mass of the rho(1700) meson",
     &DtoKPiPiCLEO::_mrho1700, MeV, 1700   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstar0892Mass
    ("Kstar0892Mass",
     "The mass of the K*0(892) meson",
     &DtoKPiPiCLEO::_mK8920, MeV, 896.1 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstarPlus892AMass
    ("KstarPlus892AMass",
     "The mass of the K*+(892) meson in D0 -> K-pi+pi0",
     &DtoKPiPiCLEO::_mK892A, MeV, 891.5 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstarPlus892BMass
    ("KstarPlus892BMass",
     "The mass of the K*+(892) meson in D0 -> K0pi+pi-",
     &DtoKPiPiCLEO::_mK892B, MeV, 891.66*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacerhoPlusMass
    ("RhoPlusMass",
     "The mass of the rho+ meson in D0 -> K-pi+pi0",
     &DtoKPiPiCLEO::_mrhoA, MeV, 770   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacerho0Mass
    ("Rho0Mass",
     "The mass of the rho+ meson in D0 -> K0pi+pi-",
     &DtoKPiPiCLEO::_mrhoB, MeV, 769.3 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &DtoKPiPiCLEO::_womega, MeV, 8.44*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacef980Width
    ("f980Width",
     "The width of the f_0(980) meson",
     &DtoKPiPiCLEO::_wf980, MeV,  50.  *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &DtoKPiPiCLEO::_wf2, MeV, 185.1 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacef1370Width
    ("f1370Width",
     "The width of the f_0(1370) meson",
     &DtoKPiPiCLEO::_wf1370, MeV, 272.0 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceK_01430Width
    ("K_01430Width",
     "The width of the K_0(1430) meson",
     &DtoKPiPiCLEO::_wK14300, MeV, 294   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceK_21430Width
    ("K_21430Width",
     "The width of the K_2(1430) meson",
     &DtoKPiPiCLEO::_wK14302, MeV,  98.5 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstar1680Width
    ("Kstar1680Width",
     "The width of the K*(1680) meson",
     &DtoKPiPiCLEO::_wK1680, MeV, 322   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacerho1700Width
    ("rho1700Width",
     "The width of the rho(1700) meson",
     &DtoKPiPiCLEO::_wrho1700, MeV, 240   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstar0892Width
    ("Kstar0892Width",
     "The width of the K*0(892) meson",
     &DtoKPiPiCLEO::_wK8920, MeV, 50.5 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstarPlus892AWidth
    ("KstarPlus892AWidth",
     "The width of the K*+(892) meson in D0 -> K-pi+pi0",
     &DtoKPiPiCLEO::_wK892A, MeV,  50   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfaceKstarPlus892BWidth
    ("KstarPlus892BWidth",
     "The width of the K*+(892) meson in D0 -> K0pi+pi-",
     &DtoKPiPiCLEO::_wK892B, MeV, 50.8 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacerhoPlusWidth
    ("RhoPlusWidth",
     "The width of the rho+ meson in D0 -> K-pi+pi0",
     &DtoKPiPiCLEO::_wrhoA, MeV, 150.7 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy> interfacerho0Width
    ("Rho0Width",
     "The width of the rho+ meson in D0 -> K0pi+pi-",
     &DtoKPiPiCLEO::_wrhoB, MeV, 150.2 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfacegPi
    ("gPi",
     "The g_pi coupling for the f_0(980) width",
     &DtoKPiPiCLEO::_gpi, 0.09, 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfacegK
    ("gK",
     "The g_K coupling for the f_0(980) width",
     &DtoKPiPiCLEO::_gK, 0.02, 0.0, 1.,
     false, false, Interface::limited);

  static Switch<DtoKPiPiCLEO,bool> interfacef0Option
    ("f0Option",
     "Option for the treatment of the f_0(980) width",
     &DtoKPiPiCLEO::_f0opt, true, false, false);
  static SwitchOption interfacef0OptionCoupled
    (interfacef0Option,
     "Coupled",
     "Use the coupling pion and kaon channels",
     true);
  static SwitchOption interfacef0OptionSWave
    (interfacef0Option,
     "SWave",
     "Use a simple s-wave running width",
     false);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedNonResonantAmplitude
    ("ChargedNonResonantAmplitude",
     "Amplitude for the non-resonant component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1NR, 1.75, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedNonResonantPhase
    ("ChargedNonResonantPhase",
     "Phase for the non-resonant component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1NR, 31.2, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedRhoAmplitude
    ("ChargedRhoAmplitude",
     "Amplitude for the rho+ component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1rho, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedRhoPhase
    ("ChargedRhoPhase",
     "Phase for the rho+ component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1rho, 0.0, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedKStarMinusAmplitude
    ("ChargedKStarMinusAmplitude",
     "Amplitude for the K*(892)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1Kstarm, 0.44, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedKStarMinusPhase
    ("ChargedKStarMinusPhase",
     "Phase for the K*(892)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1Kstarm, 163, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedKStar0Amplitude
    ("ChargedKStar0Amplitude",
     "Amplitude for the K*(892)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1Kstar0, 0.39, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedKStar0Phase
    ("ChargedKStar0Phase",
     "Phase for the K*(892)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1Kstar0, -0.2, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy2> interfaceChargedK_0MinusAmplitude
    ("ChargedK_0MinusAmplitude",
     "Amplitude for the K_0(1430)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1K1430m, GeV2, 0.77*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedK_0MinusPhase
    ("ChargedK_0MinusPhase",
     "Phase for the K_0(1430)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1K1430m, 55.5, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy2> interfaceChargedK_00Amplitude
    ("ChargedK_00Amplitude",
     "Amplitude for the K_0(1430)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1K14300, GeV2, 0.85*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedK_00Phase
    ("ChargedK_00Phase",
     "Phase for the K_0(1430)0 component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1K14300, 166, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedRho1700Amplitude
    ("ChargedRho1700Amplitude",
     "Amplitude for the rho1700+ component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1rho1700, 2.5, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedRho1700Phase
    ("ChargedRho1700Phase",
     "Phase for the rho1700+ component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1rho1700, 171., -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedK1680MinusAmplitude
    ("ChargedK1680MinusAmplitude",
     "Amplitude for the K*(1680)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_a1K1680, 2.5, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceChargedK1680MinusPhase
    ("ChargedK1680MinusPhase",
     "Phase for the K*(1680)- component for D0 -> K- pi+ pi0",
     &DtoKPiPiCLEO::_phi1K1680, 103, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralKStarPlusAmplitude
    ("NeutralKStarPlusAmplitude",
     "Amplitude for the K*(892)+ component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2Kstarp, 0.11, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralKStarPlusPhase
    ("NeutralKStarPlusPhase",
     "Phase for the K*(892)+ component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2Kstarp, 321., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralRhoAmplitude
    ("NeutralRhoAmplitude",
     "Amplitude for the rho0 component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2rho, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralRhoPhase
    ("NeutralRhoPhase",
     "Phase for the rho0 component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2rho, 0.0, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralOmegaAmplitude
    ("NeutralOmegaAmplitude",
     "Amplitude for the omega component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2omega, 0.037, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralOmegaPhase
    ("NeutralOmegaPhase",
     "Phase for the omega component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2omega, 114., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralKStarMinusAmplitude
    ("NeutralKStarMinusAmplitude",
     "Amplitude for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2Kstarm, 1.56, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralKStarMinusPhase
    ("NeutralKStarMinusPhase",
     "Phase for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2Kstarm, 150., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy2> interfaceNeutralf980Amplitude
    ("Neutralf980Amplitude",
     "Amplitude for the f_0(980) component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2f980, GeV2, 0.34*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralf980Phase
    ("Neutralf980Phase",
     "Phase for the f_0(980) component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2f980, 188., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,InvEnergy2> interfaceNeutralf2Amplitude
    ("Neutralf2Amplitude",
     "Amplitude for the f_2 component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2f2, 1./GeV2, 0.7/GeV2, ZERO, 10.0/GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralf2Phase
    ("Neutralf2Phase",
     "Phase for the f_0(2) component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2f2, 308., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy2> interfaceNeutralf1370Amplitude
    ("Neutralf1370Amplitude",
     "Amplitude for the f_0(1370) component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2f1370, GeV2, 1.8*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralf1370Phase
    ("Neutralf1370Phase",
     "Phase for the f_0(1370) component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2f1370, 85., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,Energy2> interfaceNeutralKK_0MinusAmplitude
    ("NeutralKK_0MinusAmplitude",
     "Amplitude for the K_0(1430)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2K14300, GeV2, 2.0*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralKK_0MinusPhase
    ("NeutralKK_0MinusPhase",
     "Phase for the K_0(1430)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2K14300, 3, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,InvEnergy2> interfaceNeutralKK_2MinusAmplitude
    ("NeutralKK_2MinusAmplitude",
     "Amplitude for the K_2(1430)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2K14302, 1./GeV2, 1.0/GeV2, ZERO, 10.0/GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralKK_2MinusPhase
    ("NeutralKK_2MinusPhase",
     "Phase for the K_2(1430)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2K14302, 335, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralK1680MinusAmplitude
    ("NeutralK1680MinusAmplitude",
     "Amplitude for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2K1680, 5.6, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralK1680MinusPhase
    ("NeutralK1680MinusPhase",
     "Phase for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2K1680, 174, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralNonResonantAmplitude
    ("NeutralNonResonantAmplitude",
     "Amplitude for the non-resonant component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_a2NR, 1.1, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,double> interfaceNeutralNonResonantPhase
    ("NeutralNonResonantPhase",
     "Phase for the non-resonant component for D0 -> Kbar0 pi+ pi-",
     &DtoKPiPiCLEO::_phi2NR, 340, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,InvEnergy> interfaceDRadius
    ("DRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the D",
     &DtoKPiPiCLEO::_rD0, 1./GeV, 5./GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiCLEO,InvEnergy> interfaceResonanceRadius
    ("ResonanceRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the"
     "intermediate resonances",
     &DtoKPiPiCLEO::_rres, 1./GeV, 1.5/GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiCLEO,double> interfaceMaximumWeights
    ("MaximumWeights",
     "The maximum weights for the unweighting of the decays",
     &DtoKPiPiCLEO::_maxwgt, -1, 1.0, 0.0, 1.0e11,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiCLEO,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &DtoKPiPiCLEO::_weights, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
}

int DtoKPiPiCLEO::modeNumber(bool & cc,tcPDPtr parent,
			     const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be D0
  if(abs(id0)!=ParticleID::D0) return -1;
  cc = id0==ParticleID::Dbar0;
  // must be three decay products
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int npip(0),npim(0),nkm(0),nk0(0),npi0(0);
  for( ;pit!=children.end();++pit) {
    id0=(**pit).id();
    if(id0          ==ParticleID::piplus)  ++npip;
    else if(id0     ==ParticleID::pi0)     ++npi0;
    else if(id0     ==ParticleID::piminus) ++npim;
    else if(abs(id0)==ParticleID::K0)      ++nk0;
    else if(id0     ==ParticleID::K_L0)    ++nk0;
    else if(id0     ==ParticleID::K_S0)    ++nk0;
    else if(abs(id0)==ParticleID::Kplus)   ++nkm;
  }
  if(npim==1&&npip==1&&nk0==1) return  1;
  else if(nkm==1&&(npip+npim)==1&&npi0==1) return 0;
  else                        return -1;
}

double DtoKPiPiCLEO::me2(const int ichan,
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
  // compute the invariant masses needed to calulate the amplitudes
  Energy mD  = inpart.mass();
  Energy mA  = decay[0]->mass();
  Energy mB  = decay[1]->mass();
  Energy mC  = decay[2]->mass();
  Energy mAB = (decay[0]->momentum()+decay[1]->momentum()).m();
  Energy mAC = (decay[0]->momentum()+decay[2]->momentum()).m();
  Energy mBC = (decay[1]->momentum()+decay[2]->momentum()).m();
  // compute the amplitudes for the resonaces present in both models
  Complex amp(0);
  // calculate the matrix element
  if(imode()==0) {
    if(ichan<0) {
      amp = _c1NR
	+_c1rho        *amplitude(1,false,mD,mB,mC,mA,mBC,mAB,mAC,_mrhoA   ,_wrhoA   )
	+_c1rho1700    *amplitude(1,false,mD,mB,mC,mA,mBC,mAB,mAC,_mrho1700,_wrho1700)
	+_c1Kstarm     *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK892A  ,_wK892A  )
	+_c1K1430m/GeV2*amplitude(0,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14300 ,_wK14300 )
	+_c1K1680      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1680  ,_wK1680)
	+_c1Kstar0     *amplitude(1,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK8920  ,_wK8920  )
	+_c1K14300/GeV2*amplitude(0,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14300 ,_wK14300 );
    }
    else if(ichan==0) {
      amp= _c1rho       *amplitude(1,false,mD,mB,mC,mA,mBC,mAB,mAC,_mrhoA   ,_wrhoA   );
    }
    else if(ichan==1) {
      amp=_c1Kstarm     *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK892A  ,_wK892A  ) ;
    }
    else if(ichan==2) {
      amp=_c1Kstar0     *amplitude(1,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK8920  ,_wK8920  ) ;
    }
    else if(ichan==3) {
      amp=_c1K1430m/GeV2*amplitude(0,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14300 ,_wK14300 );
    }
    else if(ichan==4) {
      amp=_c1K14300/GeV2*amplitude(0,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14300 ,_wK14300 );
    }
    else if(ichan==5) {
      amp=_c1rho1700    *amplitude(1,false,mD,mB,mC,mA,mBC,mAB,mAC,_mrho1700,_wrho1700);
    }
    else if(ichan==6) {
      amp=_c1K1680      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1680  ,_wK1680);
    }
  }
  else {
    if(ichan<0) {
      amp = _c2NR
	+_c2Kstarm     *amplitude(1,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK892B ,_wK892B )
	+Complex(_c2K14300/GeV2)*amplitude(0,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK14300,_wK14300)
	+Complex(_c2K14302*GeV2)*amplitude(2,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK14302,_wK14302)
	+_c2K1680      *amplitude(1,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK1680 ,_wK1680 )
	+_c2Kstarp     *amplitude(1,false ,mD,mA,mB,mC,mAB,mAC,mBC,_mK892B ,_wK892B )
	+_c2rho        *amplitude(1,false ,mD,mB,mC,mA,mBC,mAB,mAC,_mrhoB  ,_wrhoB  )
	+_c2omega      *amplitude(1,false ,mD,mB,mC,mA,mBC,mAB,mAC,_momega ,_womega )
	+Complex(_c2f980/GeV2  )*amplitude(0,_f0opt,mD,mB,mC,mA,mBC,mAB,mAC,_mf980  ,_wf980  )
	+Complex(_c2f1370/GeV2 )*amplitude(0,false ,mD,mB,mC,mA,mBC,mAB,mAC,_mf1370 ,_wf1370 )
	+Complex(_c2f2*GeV2    )*amplitude(2,false ,mD,mB,mC,mA,mBC,mAB,mAC,_mf2    ,_wf2    );
    }
    else if(ichan==0) {
      amp=_c2Kstarp     *amplitude(1,false ,mD,mA,mB,mC,mAB,mAC,mBC,_mK892B ,_wK892B );
    }
    else if(ichan==1) {
      amp=_c2rho        *amplitude(1,false ,mD,mB,mC,mA,mBC,mAB,mAC,_mrhoB  ,_wrhoB  );
    }
    else if(ichan==2) {
      amp=_c2omega      *amplitude(1,false ,mD,mB,mC,mA,mBC,mAB,mAC,_momega ,_womega );
    }
    else if(ichan==3) {
      amp=_c2Kstarm     *amplitude(1,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK892B ,_wK892B );
    }
    else if(ichan==4) {
      amp=_c2f980/GeV2  *amplitude(0,_f0opt,mD,mB,mC,mA,mBC,mAB,mAC,_mf980  ,_wf980  );
    }
    else if(ichan==5) {
      amp=_c2f2*GeV2    *amplitude(2,false ,mD,mB,mC,mA,mBC,mAB,mAC,_mf2    ,_wf2    );
    }
    else if(ichan==6) {
      amp=_c2f1370/GeV2 *amplitude(0,false ,mD,mB,mC,mA,mBC,mAB,mAC,_mf1370 ,_wf1370 );
    }
    else if(ichan==7) {
      amp=_c2K14300/GeV2*amplitude(0,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK14300,_wK14300);
    }
    else if(ichan==8) {
      amp=_c2K14302*GeV2*amplitude(2,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK14302,_wK14302);
    }
    else if(ichan==9) {
      amp=_c2K1680      *amplitude(1,false ,mD,mA,mC,mB,mAC,mAB,mBC,_mK1680 ,_wK1680 );
    }
  }
  // now compute the matrix element
  (*ME())(0,0,0,0)=amp;
  return norm(amp);
}

void DtoKPiPiCLEO::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // parameters
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":OmegaMass "          << _momega/MeV   << "\n";
  output << "newdef " << name() << ":f980Mass "           << _mf980/MeV    << "\n";
  output << "newdef " << name() << ":f_2Mass "            << _mf2/MeV      << "\n";
  output << "newdef " << name() << ":f1370Mass "          << _mf1370/MeV   << "\n";
  output << "newdef " << name() << ":K_01430Mass "        << _mK14300/MeV  << "\n";
  output << "newdef " << name() << ":K_21430Mass "        << _mK14302/MeV  << "\n";
  output << "newdef " << name() << ":Kstar1680Mass "      << _mK1680/MeV   << "\n";
  output << "newdef " << name() << ":rho1700Mass "        << _mrho1700/MeV << "\n";
  output << "newdef " << name() << ":Kstar0892Mass "      << _mK8920/MeV   << "\n";
  output << "newdef " << name() << ":KstarPlus892AMass "  << _mK892A/MeV   << "\n";
  output << "newdef " << name() << ":KstarPlus892BMass "  << _mK892B/MeV   << "\n";
  output << "newdef " << name() << ":RhoPlusMass "        << _mrhoA/MeV    << "\n";
  output << "newdef " << name() << ":Rho0Mass "           << _mrhoB/MeV    << "\n";
  output << "newdef " << name() << ":OmegaWidth "         << _womega/MeV   << "\n";
  output << "newdef " << name() << ":f980Width "          << _wf980/MeV    << "\n";
  output << "newdef " << name() << ":f_2Width "           << _wf2/MeV      << "\n";
  output << "newdef " << name() << ":f1370Width "         << _wf1370/MeV   << "\n";
  output << "newdef " << name() << ":K_01430Width "       << _wK14300/MeV  << "\n";
  output << "newdef " << name() << ":K_21430Width "       << _wK14302/MeV  << "\n";
  output << "newdef " << name() << ":Kstar1680Width "     << _wK1680/MeV   << "\n";
  output << "newdef " << name() << ":rho1700Width "       << _wrho1700/MeV << "\n";
  output << "newdef " << name() << ":Kstar0892Width "     << _wK8920/MeV   << "\n";
  output << "newdef " << name() << ":KstarPlus892AWidth " << _wK892A/MeV   << "\n";
  output << "newdef " << name() << ":KstarPlus892BWidth " << _wK892B/MeV   << "\n";
  output << "newdef " << name() << ":RhoPlusWidth "       << _wrhoA/MeV    << "\n";
  output << "newdef " << name() << ":Rho0Width "          << _wrhoB/MeV    << "\n";
  output << "newdef " << name() << ":gPi " << _gpi << "\n";
  output << "newdef " << name() << ":gK " << _gK << "\n";
  output << "newdef " << name() << ":f0Option " << _f0opt << "\n";
  output << "newdef " << name() << ":ChargedNonResonantAmplitude " << _a1NR << "\n";
  output << "newdef " << name() << ":ChargedNonResonantPhase " << _phi1NR<< "\n";
  output << "newdef " << name() << ":ChargedRhoAmplitude " << _a1rho<< "\n";
  output << "newdef " << name() << ":ChargedRhoPhase " << _phi1rho<< "\n";
  output << "newdef " << name() << ":ChargedKStarMinusAmplitude " << _a1Kstarm<< "\n";
  output << "newdef " << name() << ":ChargedKStarMinusPhase " << _phi1Kstarm<< "\n";
  output << "newdef " << name() << ":ChargedKStar0Amplitude " << _a1Kstar0<< "\n";
  output << "newdef " << name() << ":ChargedKStar0Phase " << _phi1Kstar0<< "\n";
  output << "newdef " << name() << ":ChargedK_0MinusAmplitude " 
	 << _a1K1430m/GeV2 << "\n";
  output << "newdef " << name() << ":ChargedK_0MinusPhase " << _phi1K1430m<< "\n";
  output << "newdef " << name() << ":ChargedK_00Amplitude " 
	 << _a1K14300/GeV2 << "\n";
  output << "newdef " << name() << ":ChargedK_00Phase " << _phi1K14300<< "\n";
  output << "newdef " << name() << ":ChargedRho1700Amplitude " << _a1rho1700<< "\n";
  output << "newdef " << name() << ":ChargedRho1700Phase " << _phi1rho1700<< "\n";
  output << "newdef " << name() << ":ChargedK1680MinusAmplitude " << _a1K1680<< "\n";
  output << "newdef " << name() << ":ChargedK1680MinusPhase " << _phi1K1680<< "\n";
  output << "newdef " << name() << ":NeutralKStarPlusAmplitude " << _a2Kstarp<< "\n";
  output << "newdef " << name() << ":NeutralKStarPlusPhase " << _phi2Kstarp<< "\n";
  output << "newdef " << name() << ":NeutralRhoAmplitude " << _a2rho<< "\n";
  output << "newdef " << name() << ":NeutralRhoPhase " << _phi2rho<< "\n";
  output << "newdef " << name() << ":NeutralOmegaAmplitude " << _a2omega<< "\n";
  output << "newdef " << name() << ":NeutralOmegaPhase " <<_phi2omega << "\n";
  output << "newdef " << name() << ":NeutralKStarMinusAmplitude " << _a2Kstarm<< "\n";
  output << "newdef " << name() << ":NeutralKStarMinusPhase " << _phi2Kstarm<< "\n";
  output << "newdef " << name() << ":Neutralf980Amplitude " << _a2f980/GeV2<< "\n";
  output << "newdef " << name() << ":Neutralf980Phase " << _phi2f980<< "\n";
  output << "newdef " << name() << ":Neutralf2Amplitude " << _a2f2*GeV2<< "\n";
  output << "newdef " << name() << ":Neutralf2Phase " << _phi2f2<< "\n";
  output << "newdef " << name() << ":Neutralf1370Amplitude " << _a2f1370/GeV2<< "\n";
  output << "newdef " << name() << ":Neutralf1370Phase " << _phi2f1370<< "\n";
  output << "newdef " << name() << ":NeutralKK_0MinusAmplitude " 
	 << _a2K14300/GeV2 << "\n";
  output << "newdef " << name() << ":NeutralKK_0MinusPhase " << _phi2K14300 << "\n";
  output << "newdef " << name() << ":NeutralKK_2MinusAmplitude " 
	 << _a2K14302*GeV2<< "\n";
  output << "newdef " << name() << ":NeutralKK_2MinusPhase " << _phi2K14302 << "\n";
  output << "newdef " << name() << ":NeutralK1680MinusAmplitude " << _a2K1680<< "\n";
  output << "newdef " << name() << ":NeutralK1680MinusPhase " << _phi2K1680<< "\n";
  output << "newdef " << name() << ":NeutralNonResonantAmplitude " << _a2NR<< "\n";
  output << "newdef " << name() << ":NeutralNonResonantPhase " << _phi2NR << "\n";
  output << "newdef " << name() << ":DRadius " << _rD0*GeV << "\n";
  output << "newdef " << name() << ":ResonanceRadius " << _rres*GeV << "\n";
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

void DtoKPiPiCLEO::doinitrun() {
  DecayIntegrator::doinitrun();
  _weights.resize(mode(0)->numberChannels()+mode(1)->numberChannels());
  _maxwgt.resize(2);
  unsigned int iy=0;
  for(unsigned int ix=0;ix<2;++ix) {
    _maxwgt[ix]=mode(ix)->maxWeight();
    for(unsigned int iz=0;iz<mode(ix)->numberChannels();++iz) {
      _weights[iy]=mode(ix)->channelWeight(iz);
      ++iy;
    }
  }
}

Complex DtoKPiPiCLEO::amplitude(int ispin,bool f0, Energy mD,
					Energy mA , Energy mB , Energy mC ,
					Energy mAB, Energy mAC, Energy mBC,
					Energy mres, Energy wres) const{
  // compute the production momenta
  Energy pDR  = Kinematics::pstarTwoBodyDecay(mD,mres,mC);
  Energy pDAB = Kinematics::pstarTwoBodyDecay(mD,mAB ,mC);
  // and the decay momenta
  Energy pAB = Kinematics::pstarTwoBodyDecay(mAB ,mA,mB);
  Energy pR  = Kinematics::pstarTwoBodyDecay(mres,mA,mB);
  double Fd(1.),Fr(1.),s(1.);
  switch(ispin) {
  case 0:
    // default values of parameters are correct
    break;
  case 1:
    Fr = sqrt((1.+sqr(_rres*pR ))/(1.+sqr(_rres*pAB )));
    Fd = sqrt((1.+sqr(_rD0 *pDR))/(1.+sqr(_rD0 *pDAB)));
    s = ((mAC-mBC)*(mAC+mBC)+(mD-mC)*(mD+mC)*(mB-mA)*(mB+mA)/sqr(mres))/GeV2;
    break;
  case 2:
    Fr = sqrt((9.+3.*sqr(_rres*pR  )+Math::powi(_rres*pR  ,4))/
	      (9.+3.*sqr(_rres*pAB )+Math::powi(_rres*pAB ,4)));
    Fd = sqrt((9.+3.*sqr(_rD0 *pDR )+Math::powi(_rD0 *pDR ,4))/
	      (9.+3.*sqr(_rD0 *pDAB)+Math::powi(_rD0 *pDAB,4)));
    s = sqr(((mBC-mAC)*(mBC+mAC)+(mD-mC)*(mD+mC)*(mA-mB)*(mA+mB)/sqr(mres))/GeV2)
      -(mAB*mAB-2.*mD*mD-2.*mC*mC+sqr((mD-mC)*(mD+mC))/sqr(mres))*
       (mAB*mAB-2.*mA*mA-2.*mB*mB+sqr((mA-mB)*(mA+mB))/sqr(mres))/3./GeV2/GeV2;
    break;
  default:
    throw Exception() << "D0toKPiPiCLEO::amplitude spin is too high ispin = " 
		      << ispin << Exception::runerror;
  }
  // calculate the width term
  complex<Energy2> bw;
  if(!f0) {
    Energy2 mwid=wres*Math::powi(pAB/pR,2*ispin+1)*sqr(Fr*mres)/mAB;
    bw = sqr(mres)-sqr(mAB)-complex<Energy2>(ZERO,mwid);
  }
  else {
    Energy Gamma_pi = _gpi*sqrt(0.25*sqr(mAB)-sqr(_mpi));
    Energy Gamma_K  = 0.5*_gK *(sqrt(0.25*sqr(mAB)-sqr(_mkp))+
				sqrt(0.25*sqr(mAB)-sqr(_mk0)));
    bw = sqr(mres)-sqr(mAB)-complex<Energy2>(ZERO,mres*(Gamma_pi+Gamma_K));
  }
  return s*Fr*Fd*GeV2/bw;
}
