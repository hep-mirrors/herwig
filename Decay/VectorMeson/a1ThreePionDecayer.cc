// -*- C++ -*-
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

using namespace Herwig;
using namespace ThePEG::Helicity;

a1ThreePionDecayer::a1ThreePionDecayer() 
  : _rhomass(1,0.7761*GeV), _rhowidth(1,0.1445*GeV), _sigmamass(0.8*GeV),
    _sigmawidth(0.8*GeV), _lambda2(1.2*GeV2), _a1mass2(1.23*1.23*GeV2),
    _zsigma(1.269,0.591), _rhocoupling(1,1.), _coupling(80.76), 
    _localparameters(true), _zerowgts(3), _onewgts(7), _twowgts(7),
    _threewgts(8) ,_zeromax(107.793), _onemax(1088.96), 
    _twomax(1750.73), _threemax(739.334) {
  // weights for the different channels
  _threewgts[0] = 0.140498;
  _threewgts[1] = 0.140619;
  _threewgts[2] = 0.113731;
  _threewgts[3] = 0.113435;
  _threewgts[4] = 0.121519;
  _threewgts[5] = 0.121633;
  _threewgts[6] = 0.124574;
  _threewgts[7] = 0.123991;
  _onewgts[0]  = 0.170911;
  _onewgts[1]  = 0.170866;
  _onewgts[2] = 0.128304;
  _onewgts[3] = 0.127611;
  _onewgts[4] = 0.134997;
  _onewgts[5] = 0.134879;
  _onewgts[6] = 0.132432;
  _zerowgts[0] = 0.333382;
  _zerowgts[1] = 0.332914;
  _zerowgts[2] = 0.333704;
  _twowgts[0] = 0.170642;
  _twowgts[1] = 0.170989;
  _twowgts[2] = 0.127737;
  _twowgts[3] = 0.127527;
  _twowgts[4] = 0.135416;
  _twowgts[5] = 0.135579;
  _twowgts[6] = 0.132111;
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

void a1ThreePionDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // set up the integration channels
  PDVector extpart;
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
  extpart.resize(4);
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
}
  
int a1ThreePionDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  if(dm.products().size()!=3) return -1;
  int id(dm.parent()->id());
  // check the pions
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
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
     << _twowgts << _threewgts << _zeromax
     << _onemax << _twomax << _threemax << _coupling;
}
  
void a1ThreePionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) >> iunit(_prho,GeV) 
     >> iunit(_hm2,GeV2) >> iunit(_rhoD,GeV2) >> _dhdq2m2 >>  iunit(_sigmamass,GeV)
     >> iunit(_sigmawidth,GeV) >> iunit(_psigma,GeV) >> iunit(_mpi,GeV) 
     >> iunit(_mpi2,GeV2) >> iunit(_lambda2,GeV2) >> iunit(_a1mass2,GeV2) >> _zsigma 
     >> _rhocoupling >> _coupling >> _localparameters >> _zerowgts >> _onewgts 
     >> _twowgts >> _threewgts >> _zeromax
     >> _onemax >> _twomax >> _threemax >> _coupling;
}

ClassDescription<a1ThreePionDecayer> a1ThreePionDecayer::inita1ThreePionDecayer;
// Definition of the static class description member.
  
void a1ThreePionDecayer::Init() {
    
  static ClassDocumentation<a1ThreePionDecayer> documentation
    ("The a1ThreePionDecayer class is designed to decay the a_1 "
     "resonance to three pions using a model based on that used in the modelling "
     "of tau->4 pions.");

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
     "Default",
     "Use the values from the particleData objects",
     false);

  static Parameter<a1ThreePionDecayer,double> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1ThreePionDecayer::_coupling, 80.76, 0.0, 1000.0,
     false, false, true);

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
     &a1ThreePionDecayer::_sigmamass, GeV, 0.8*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy> interfaceSigmaWidth
    ("SigmaWidth",
     "The local value of the sigma width",
     &a1ThreePionDecayer::_sigmawidth, GeV, 0.8*GeV, 0.0*GeV, 10.0*GeV,
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

  static ParVector<a1ThreePionDecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1ThreePionDecayer::_rhomass,
     MeV, 0, 0*MeV, -10000*MeV, 10000*MeV, false, false, true);

  static ParVector<a1ThreePionDecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1ThreePionDecayer::_rhowidth,
     MeV, 0, 0*MeV, -10000*MeV, 10000*MeV, false, false, true);

}

double a1ThreePionDecayer::me2(bool vertex, const int ichan,
			       const Particle & inpart,
			       const ParticleVector & decay) const {
  // wavefunctions for the decaying particles
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // create the spin information for the decay products if needed
  unsigned int ix;
  if(vertex) {
    for(ix=0;ix<decay.size();++ix) {
      // workaround for gcc 3.2.3 bug
      //ALB {ScalarWaveFunction(decay[ix],outgoing,true,vertex);}}
	PPtr mytemp = decay[ix] ; ScalarWaveFunction(mytemp,outgoing,true,vertex);
    }
  }
  // momentum of the incoming particle
  Lorentz5Momentum Q=inpart.momentum();
  // vector for the output
  LorentzVector<complex<Energy3> > output;
  // a_10 -> pi0pi0pi0
  if(imode()==0) {
    // the momenta
    Lorentz5Momentum pa=decay[1]->momentum()+decay[2]->momentum();
    pa.rescaleMass();
    Lorentz5Momentum pb=decay[0]->momentum()+decay[2]->momentum();
    pb.rescaleMass();
    Lorentz5Momentum pc=decay[0]->momentum()+decay[1]->momentum();
    pc.rescaleMass();
    //the breit-wigners
    Complex sig1=sigmaBreitWigner(pa.mass2());
    Complex sig2=sigmaBreitWigner(pb.mass2());
    Complex sig3=sigmaBreitWigner(pc.mass2());
    // compute the vector
    LorentzPolarizationVectorE tmpoutput;
    if(ichan<0) {
      tmpoutput= sig1*(decay[0]->momentum())+sig2*(decay[1]->momentum())
	+sig3*(decay[2]->momentum());
    }
    else if(ichan==15) tmpoutput=sig1*(decay[0]->momentum());
    else if(ichan==16) tmpoutput=sig2*(decay[1]->momentum());
    else if(ichan==17) tmpoutput=sig3*(decay[2]->momentum());
    // the coupling z and identical particle factor
    output = tmpoutput * _zsigma* 1./sqrt(6.) *Q.mass2();
  }
  // a_1+ -> pi0pi0pi+
  else if(imode()==1) {
    Lorentz5Momentum pa=decay[2]->momentum()+decay[0]->momentum();
    pa.rescaleMass();
    Lorentz5Momentum pb=decay[2]->momentum()+decay[1]->momentum();
    pb.rescaleMass();
    Lorentz5Momentum pc=decay[0]->momentum()+decay[1]->momentum();
    pc.rescaleMass();
    // scalar propagator
    Complex sig1 = sigmaBreitWigner(pc.mass2());
    // sigma terms
    if(ichan<0||ichan==14) 
      output = _zsigma*Q.mass2()*sig1*decay[2]->momentum();
    // the rho terms
    complex<Energy2> dot01=Q*decay[0]->momentum();
    complex<Energy2> dot02=Q*decay[1]->momentum();
    complex<Energy2> dot03=Q*decay[2]->momentum();
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1=_rhocoupling[ix]*rhoBreitWigner(pa.mass2(),ix);
      Complex rho2=_rhocoupling[ix]*rhoBreitWigner(pb.mass2(),ix);
      if(ichan<0||ichan==8+2*ix) {
	output +=rho1*(dot03*(decay[1]->momentum())-
		       dot02*(decay[2]->momentum()));
      }
      if(ichan<0||ichan==9+2*ix){
	output +=rho2*(dot03*(decay[0]->momentum())-
		       dot01*(decay[2]->momentum()));
      }
    }
    // the identical particle factor
    output *= 1./sqrt(2.);
  }
  // a_10->pi+pi-pi0
  else if(imode()==2) {
      // the momenta
      Lorentz5Momentum pa=decay[2]->momentum()+decay[1]->momentum();
      pa.rescaleMass();
      Lorentz5Momentum pb=decay[2]->momentum()+decay[0]->momentum();
      pb.rescaleMass();
      Lorentz5Momentum pc=decay[0]->momentum()+decay[1]->momentum();
      pc.rescaleMass();
      // the sigma terms
      Complex sig1=sigmaBreitWigner(pc.mass2());
      if(ichan<0||ichan==24)
	{output = _zsigma*Q.mass2()*sig1*decay[2]->momentum();}
      // rho terms
      complex<Energy2> dot01=Q*decay[0]->momentum();
      complex<Energy2> dot02=Q*decay[1]->momentum();
      complex<Energy2> dot03=Q*decay[2]->momentum();
      for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
	Complex rho1=_rhocoupling[ix]*rhoBreitWigner(pa.mass2(),ix);
	Complex rho2=_rhocoupling[ix]*rhoBreitWigner(pb.mass2(),ix);
	if(ichan<0||ichan==18+2*ix) {
	  output+=rho1*(dot03*(decay[1]->momentum())
			-dot02*(decay[2]->momentum()));
	}
	if(ichan<0||ichan==19+2*ix) {
	  output+=rho2*(dot03*(decay[0]->momentum())
			-dot01*(decay[2]->momentum()));
	}
      }
  }
  // a1+ -> pi+pi+pi-
  else if(imode()==3) {
    // momenta of the intermediates
    Lorentz5Momentum pa=decay[1]->momentum()+decay[2]->momentum();
    pa.rescaleMass();
    Lorentz5Momentum pb=decay[0]->momentum()+decay[2]->momentum();
    pb.rescaleMass();
    Lorentz5Momentum pc=decay[0]->momentum()+decay[1]->momentum();
    pc.rescaleMass();
    // the scalar propagators 
    Complex sig1=sigmaBreitWigner(pa.mass2());
    Complex sig2=sigmaBreitWigner(pb.mass2());
    // sigma terms
    LorentzPolarizationVectorE tmpoutput;
    if(ichan<0||ichan==6) tmpoutput+=sig1*(decay[0]->momentum());
    if(ichan<0||ichan==7) tmpoutput+=sig2*(decay[1]->momentum());
    output = tmpoutput * _zsigma * Q.mass2();
    // rho terms
    complex<Energy2> dot01=Q*decay[0]->momentum();
    complex<Energy2> dot02=Q*decay[1]->momentum();
    complex<Energy2> dot03=Q*decay[2]->momentum();
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1 = _rhocoupling[ix]*rhoBreitWigner(pa.mass2(),ix);
      Complex rho2 = _rhocoupling[ix]*rhoBreitWigner(pb.mass2(),ix);
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
    = output * a1FormFactor(Q.mass2())*_coupling/(Q.mass()*_rhomass[0]*_rhomass[0]);
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  for(unsigned int ix=0;ix<3;++ix){newME(ix,0,0,0)=outputFinal.dot(invec[ix]);}
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}
