// -*- C++ -*-
//
// a1ThreePionDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1ThreePionDecayer class.
//

#include "a1ThreePionDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
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
    _localparameters(true),
    _zerowgts ({0.339108,0.335601,0.325291}),
    _onewgts  ({0.19616 ,0.191408,0.12137 ,0.115498,0.12729 ,0.127183,0.12109 }),
    _twowgts  ({0.188163,0.192479,0.121658,0.12135 ,0.127298,0.124835,0.124217}),
    _threewgts({0.153071,0.165741,0.107509,0.10275 ,0.109738,0.11254 ,0.125344,0.123307}) ,
    _zeromax(19.144), _onemax(7.83592), 
    _twomax(6.64804), _threemax(6.66296) {
  // generation of intermediates
  generateIntermediates(true);
}

void a1ThreePionDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    // get the weights for the different channels
    for(unsigned int ix=0;ix<_zerowgts.size();++ix)
      _zerowgts[ix]=mode(0)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_onewgts.size();++ix)
      _onewgts[ix]=mode(1)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_twowgts.size();++ix)
      _twowgts[ix]=mode(2)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_threewgts.size();++ix)
      _threewgts[ix]=mode(3)->channels()[ix].weight();
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
  // decay mode a_0 -> pi0 pi0 pi0
  tPDPtr in = a10;
  tPDVector out={pi0,pi0,pi0};
  if(sigma) {
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,_zeromax));
    if(_zerowgts.size()!=3) _zerowgts=vector<double>(3,1./3.);
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,sigma,0,1,1,2,1,3));
    c1.weight(_zerowgts[0]);
    mode->addChannel(c1);
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,sigma,0,2,1,1,1,3));
    c2.weight(_zerowgts[1]);
    mode->addChannel(c2);
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
    c3.weight(_zerowgts[2]);
    mode->addChannel(c3);
    addMode(mode);
  }
  // decay mode a_1+ -> pi+ pi0 pi0
  in = a1p;
  out = {pi0,pi0,pip};
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,_onemax));
  unsigned int nrho(0);
  for(unsigned int ix=0;ix<3;++ix) if(rhop[ix]) ++nrho;
  if(_onewgts.size()!=2*nrho+1) _onewgts=vector<double>(2*nrho+1,1./(2.*nrho+1.));
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho+ channel
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rhop[ix],0,1,1,2,1,3));
    c1.weight(_onewgts[2*ix]);
    mode->addChannel(c1);
    // second rho+ channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rhop[ix],0,2,1,1,1,3));
    c2.weight(_onewgts[2*ix+1]);
    mode->addChannel(c2);
  }
  // the sigma channel
  if(sigma) {
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
    c3.weight(_onewgts.back());
    mode->addChannel(c3);
  }
  addMode(mode);
  // decay mode a_1 -> pi+ pi- pi0
  in = a10;
  out = {pip,pim,pi0};
  mode = new_ptr(PhaseSpaceMode(in,out,_twomax));
  if(_twowgts.size()!=2*nrho+1) _twowgts=vector<double>(2*nrho+1,1./(2.*nrho+1.));
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho channel
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rhop[ix],0,1,1,2,1,3));
    c1.weight(_twowgts[2*ix]);
    mode->addChannel(c1);
    // second channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rhom[ix],0,2,1,1,1,3));
    c2.weight(_twowgts[2*ix+1]);
    mode->addChannel(c2);
  }
  // sigma channel
  if(sigma) {
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
    c3.weight(_twowgts.back());
    mode->addChannel(c3);
  }
  addMode(mode);
  // decay mode a_1+ -> pi+ pi+ pi-
  in = a1p;
  out = {pip,pip,pim};
  mode = new_ptr(PhaseSpaceMode(in,out,_threemax));
  nrho = 0;
  for(unsigned int ix=0;ix<3;++ix) if(rho0[ix]) ++nrho;
  if(_threewgts.size()!=2*nrho+2) _threewgts=vector<double>(2*nrho+2,1./(2.*nrho+2.));
  for(unsigned int ix=0;ix<3;++ix) {
    // the neutral rho channels
    if(!rho0[ix]) continue;
    // the neutral rho channels
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rho0[ix],0,1,1,2,1,3));
    c1.weight(_threewgts[2*ix]);
    mode->addChannel(c1);
    // interchanged channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rho0[ix],0,2,1,1,1,3));
    c2.weight(_threewgts[2*ix+1]);
    mode->addChannel(c2);
  }
  // the sigma channels
  if(sigma) {
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,1,1,2,1,3));
    c3.weight(_threewgts[6]);
    mode->addChannel(c3);
    PhaseSpaceChannel c4((PhaseSpaceChannel(mode),0,sigma,0,2,1,1,1,3));
    c4.weight(_threewgts[7]);
    mode->addChannel(c4);
  }
  addMode(mode);
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

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<a1ThreePionDecayer,DecayIntegrator>
describeHerwiga1ThreePionDecayer("Herwig::a1ThreePionDecayer", "HwVMDecay.so");
  
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

void   a1ThreePionDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double a1ThreePionDecayer::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  useMe();
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  // momentum of the incoming particle
  Lorentz5Momentum Q=part.momentum();
  // momenta of the intermediates
  Energy2 s1=(momenta[1]+momenta[2]).m2();
  Energy2 s2=(momenta[0]+momenta[2]).m2();
  Energy2 s3=(momenta[0]+momenta[1]).m2();
  Energy2 dot01=Q*momenta[0];
  Energy2 dot02=Q*momenta[1];
  Energy2 dot03=Q*momenta[2];
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
      tmpoutput= sig1*(momenta[0])+sig2*(momenta[1])
	+sig3*(momenta[2]);
    }
    else if(ichan==0) tmpoutput=sig1*(momenta[0]);
    else if(ichan==1) tmpoutput=sig2*(momenta[1]);
    else if(ichan==2) tmpoutput=sig3*(momenta[2]);
    // the coupling z and identical particle factor
    output = tmpoutput * _zsigma* 1./sqrt(6.) *Q.mass2();
  }
  // a_1+ -> pi0pi0pi+
  else if(imode()==1) {
    // scalar propagator
    Complex sig1 = sigmaBreitWigner(s3);
    // sigma terms
    if(ichan<0||ichan==6) 
      output = _zsigma*Q.mass2()*sig1*momenta[2];
    // the rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1=_rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2=_rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output +=rho1*(dot03*(momenta[1])-
		       dot02*(momenta[2]));
      }
      if(ichan<0||ichan==2*ix+1){
	output +=rho2*(dot03*(momenta[0])-
		       dot01*(momenta[2]));
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
      output = _zsigma*Q.mass2()*sig1*momenta[2];
    // rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1=_rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2=_rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output+=rho1*(dot03*(momenta[1])
		      -dot02*(momenta[2]));
      }
      if(ichan<0||ichan==2*ix+1) {
	output+=rho2*(dot03*(momenta[0])
		      -dot01*(momenta[2]));
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
    if(ichan<0||ichan==6) tmpoutput+=sig1*(momenta[0]);
    if(ichan<0||ichan==7) tmpoutput+=sig2*(momenta[1]);
    output = tmpoutput * _zsigma * Q.mass2();
    // rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1 = _rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2 = _rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output-=rho1*( dot03*(momenta[1])-
		       dot02*(momenta[2]));
      }
      if(ichan<0||ichan==2*ix+1) {
	output-=rho2*( dot03*(momenta[0])-
		       dot01*(momenta[2]));
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
  // double test = threeBodyMatrixElement(imode(),sqr(part.mass()),s3,s2,s1,
  // 				       momenta[0].mass(),momenta[1].mass(), 
  // 				       momenta[2].mass());
  // if(ichan<0) cerr << "testing matrix element " << part.PDGName() << " -> "
  // 		   << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  // 		   << outgoing[2]->PDGName() << " " << out << " " << test << " " 
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
