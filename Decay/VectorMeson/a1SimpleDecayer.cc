// -*- C++ -*-
//
// a1SimpleDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1SimpleDecayer class.
//

#include "a1SimpleDecayer.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/WidthCalculatorBase.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

a1SimpleDecayer::a1SimpleDecayer() 
  : _rhomass(3), _rhowidth(3), _rhowgts(3),_localparameters(true), 
    _coupling(47.95/GeV), _onemax(5.4474), _twomax(5.47784), 
    _threemax(5.40185), _onewgts(6), 
    _twowgts(6), _threewgts(6), _mpi(ZERO) {
  // rho masses, widths and weights
  _rhomass[0] = 0.773*GeV; _rhowidth[0] = 0.145*GeV; _rhowgts[0] =  1.0;  
  _rhomass[1] = 1.370*GeV; _rhowidth[1] = 0.510*GeV; _rhowgts[1] = -0.145;
  _rhomass[2] = 1.750*GeV; _rhowidth[2] = 0.120*GeV; _rhowgts[2] =  0.;
  // integration weights
  _onewgts[0] = 0.235562; _twowgts[0] = 0.236208; _threewgts[0] = 0.234259;
  _onewgts[1] = 0.231098; _twowgts[1] = 0.229481; _threewgts[1] = 0.233634;
  _onewgts[2] = 0.131071; _twowgts[2] = 0.131169; _threewgts[2] = 0.135922;
  _onewgts[3] = 0.131135; _twowgts[3] = 0.133604; _threewgts[3] = 0.129231;
  _onewgts[4] = 0.135841; _twowgts[4] = 0.132685; _threewgts[4] = 0.133949;
  _onewgts[5] = 0.135294; _twowgts[5] = 0.136854; _threewgts[5] = 0.133005;
  generateIntermediates(true);
}

void a1SimpleDecayer::doinit() {
  DecayIntegrator::doinit();
  // pointers to the particles we need as external particles
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  // the different rho resonances
  tPDPtr rhop[3] = {getParticleData(213),getParticleData(100213),
		    getParticleData(30213)};
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),
		    getParticleData(30113)};
  tPDPtr rhom[3] = {getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  tPDVector extpart(4);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
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
    newchannel->addIntermediate(a1p,0,0.0,-1,2);
    newchannel->addIntermediate(rhop[ix],0,0.0,1,3);
    mode->addChannel(newchannel);
    // second rho+ channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,1);
    newchannel->addIntermediate(rhop[ix],0,0.0,2,3);
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
    newchannel->addIntermediate(a10,0,0.0,-1,2);
    newchannel->addIntermediate(rhop[ix],0,0.0,1,3);
    mode->addChannel(newchannel);
    // second channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a10,0,0.0,-1,1);
    newchannel->addIntermediate(rhom[ix],0,0.0,2,3);
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
    if(!rho0[ix]) continue;
    // the neutral rho channels
    // first channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,2);
    newchannel->addIntermediate(rho0[ix],0,0.0,1,3);
    mode->addChannel(newchannel);
    // interchanged channel
    newchannel = new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(a1p,0,0.0,-1,1);
    newchannel->addIntermediate(rho0[ix],0,0.0,2,3);
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
      if(!rhop[ix]) continue;
      _rhomass[ix]=rhop[ix]->mass();
      _rhowidth[ix]=rhop[ix]->width();
    }
  }
  _mpi = pip->mass();
}

void a1SimpleDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    // get the weights for the different channels
    for(unsigned int ix=0;ix<_onewgts.size();++ix)
      _onewgts[ix]=mode(0)->channelWeight(ix);
    for(unsigned int ix=0;ix<_twowgts.size();++ix)
      _twowgts[ix]=mode(1)->channelWeight(ix);
    for(unsigned int ix=0;ix<_threewgts.size();++ix)
      _threewgts[ix]=mode(2)->channelWeight(ix);
    // get the maximum weight
    _onemax   = mode(0)->maxWeight();
    _twomax   = mode(1)->maxWeight();
    _threemax = mode(2)->maxWeight();
  }
}

void a1SimpleDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV) << _rhowgts 
     << _localparameters << ounit(_coupling,1./GeV) << _onemax
     << _twomax << _threemax << _onewgts << _twowgts << _threewgts
     << ounit(_mpi,GeV);
}

void a1SimpleDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) >> _rhowgts 
     >> _localparameters >> iunit(_coupling,1./GeV) >> _onemax
     >> _twomax >> _threemax >> _onewgts >> _twowgts >> _threewgts
     >> iunit(_mpi,GeV);
}

ClassDescription<a1SimpleDecayer> a1SimpleDecayer::inita1SimpleDecayer;
// Definition of the static class description member.

void a1SimpleDecayer::Init() {

  static ClassDocumentation<a1SimpleDecayer> documentation
    ("The a1SimpleDecayer class implements a simple model for the decay of"
     " the a_1 to three pions based on the approach of Kuhn and Santanmaria,"
     " Z.Phys. C48, 445 (1990)",
     "The decays of the $a_1$ were modelled using the approach of "
     "\\cite{Kuhn:1990ad}.\n",
     "\\bibitem{Kuhn:1990ad} J.~H.~Kuhn and A.~Santamaria,\n"
     "Z.\\ Phys.\\  C {\\bf 48} (1990) 445.\n"
     "%%CITATION = ZEPYA,C48,445;%%\n");

  static Switch<a1SimpleDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1SimpleDecayer::_localparameters, true, false, false);
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

  static Parameter<a1SimpleDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1SimpleDecayer::_coupling, 1./GeV, 47.95/GeV, ZERO, 100./GeV,
     false, false, Interface::limited);

  static ParVector<a1SimpleDecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1SimpleDecayer::_rhomass,
     MeV, 0, ZERO, ZERO, 10000*MeV, false, false, true);

  static ParVector<a1SimpleDecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1SimpleDecayer::_rhowidth,
     MeV, 0, ZERO, ZERO, 10000*MeV, false, false, true);

  static ParVector<a1SimpleDecayer,double> interfaceRhoWeights
    ("RhoWeights",
     "Weight for the different rho resonances",
     &a1SimpleDecayer::_rhowgts, -1, 0.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<a1SimpleDecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1SimpleDecayer::_onemax, 5.57613, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1SimpleDecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1SimpleDecayer::_twomax, 5.61218, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1SimpleDecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1SimpleDecayer::_threemax, 5.5384, 0.0, 10000.0,
     false, false, true);
  
  static ParVector<a1SimpleDecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1SimpleDecayer::_onewgts,
     0, 0, 0, 0., 1., false, false, true);

  static ParVector<a1SimpleDecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1SimpleDecayer::_twowgts,
     0, 0, 0, 0., 1., false, false, true);

  static ParVector<a1SimpleDecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1SimpleDecayer::_threewgts,
     0, 0, 0, 0., 1., false, false, true);

}

int a1SimpleDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(npiplus==1&&npi0==2)          imode=0;
    else if(npiplus==2&&npiminus==1) imode=2;
  }
  // a_1- modes
  else if(id==ParticleID::a_1minus) {
    cc=true;
    if(npiminus==1&&npi0==2)         imode=0;
    else if(npiminus==2&&npiplus==1) imode=2;
  }
  // a_0 modes
  else if(id==ParticleID::a_10) {
    cc=false;
    if(npiminus==1&&npiplus==1&&npi0==1) imode=1;
  }
  return imode;
}

double a1SimpleDecayer::me2(const int ichan,const Particle & inpart,
			    const ParticleVector & decay,MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
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
  Lorentz5Vector<complex<Energy> > current;
  Energy2 s1 = (decay[1]->momentum()+decay[2]->momentum()).m2();
  Energy2 s2 = (decay[0]->momentum()+decay[2]->momentum()).m2();
  if(ichan<0) {
    current = rhoFormFactor(s2,-1)*(decay[0]->momentum()-decay[2]->momentum())
    +rhoFormFactor(s1,-1)*(decay[1]->momentum()-decay[2]->momentum());
  }
  else if(ichan<3) {
    current = 
      rhoFormFactor(s2,ichan)*(decay[0]->momentum()-decay[2]->momentum());
  }
  else if(ichan<6) {
    current = 
      rhoFormFactor(s1,-1)*(decay[1]->momentum()-decay[2]->momentum());
  }
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix)
    (*ME())(ix,0,0,0)=_coupling*current.dot(_vectors[ix]);
  // matrix element and identical particle factor
  double output=ME()->contract(_rho).real();
  if(imode()!=1) output*=0.5;
  // test the output
//   double test = threeBodyMatrixElement(imode(),sqr(inpart.mass()),
// 				       s3,s2,s1,decay[0]->mass(),decay[1]->mass(), 
// 				       decay[2]->mass());
//   if(ichan<0) cerr << "testing matrix element " << inpart.PDGName() << " -> "
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << decay[2]->PDGName() << output << " " << test << " " 
//        << (output-test)/(output+test) << "\n";  
  // return the answer
  return output;
}

double a1SimpleDecayer::
threeBodyMatrixElement(const int iopt,const Energy2 q2, const Energy2 s3,
		       const Energy2 s2,const Energy2 s1, const Energy m1, 
		       const Energy m2 ,const Energy m3) const {
  Energy2 v12  = (s2-2.*sqr(m1)-2.*sqr(m3))+0.25*sqr(s1-s3-sqr(m1)+sqr(m3))/q2;
  Energy2 v22  = (s1-2.*sqr(m2)-2.*sqr(m3))+0.25*sqr(s2-s3-sqr(m2)+sqr(m3))/q2;
  Energy2 v1v2 = (0.5*q2-s3-0.5*(3*sqr(m3)-sqr(m1)-sqr(m2)))
    +0.25*(s1-s3-sqr(m1)+sqr(m3))*(s2-s3-sqr(m2)+sqr(m3))/q2;
  Complex rho1=rhoFormFactor(s2,-1);
  Complex rho2=rhoFormFactor(s1,-1);
  double me = sqr(_coupling)*real(v12*rho1*conj(rho1)+v22*rho2*conj(rho2)
				  +2.*v1v2*rho1*conj(rho2))/3.;
  if(iopt!=1) me *= 0.5;
  return me;
}

WidthCalculatorBasePtr
a1SimpleDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit) {
    if(abs((**pit).id())==ParticleID::piplus) ++ncharged;
  }
  --ncharged;
  // integrator to perform the integral
  vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
  vector<int> intype;intype.push_back(2);intype.push_back(3);
  vector<Energy> inmass(2,_rhomass[0]),inwidth(2,_rhowidth[0]);
  vector<double> inpow(2,0.0);
  Energy mpi0=getParticleData(ParticleID::pi0)->mass();
  Energy mpic=getParticleData(ParticleID::piplus)->mass();
  Energy m[3];
  m[0] = ncharged<2 ? mpi0 : mpic;
  m[1] = m[0];
  m[2] = (ncharged==0||ncharged==2) ? mpi0 : mpic;
  return new_ptr(ThreeBodyAllOnCalculator<a1SimpleDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,ncharged,m[0],m[1],m[2]));
}


void a1SimpleDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Coupling " << _coupling*GeV << "\n";
  output << "newdef " << name() << ":OneMax   " <<   _onemax << "\n";
  output << "newdef " << name() << ":TwoMax   " <<   _twomax << "\n";
  output << "newdef " << name() << ":ThreeMax " << _threemax << "\n";
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    if(ix<3) output << "newdef    " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/MeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/MeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    if(ix<3) output << "newdef    " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/MeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/MeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowgts.size();++ix) {
    if(ix<3) output << "newdef    " << name() << ":RhoWeights " << ix << " " 
		    << _rhowgts[ix] << "\n";
    else     output << "insert " << name() << ":RhoWeights " << ix << " " 
		    << _rhowgts[ix] << "\n";
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
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

// functions to return the Breit-Wigners
Complex a1SimpleDecayer::rhoFormFactor(Energy2 q2,int ires) const {
  Complex output(0.),norm(0.);
  for(unsigned int ix=0,N=min(3,int(_rhowgts.size()));ix<N;++ix)
    norm+=_rhowgts[ix];
  if(ires<0) {
    for(unsigned int ix=0,N=min(3,int(_rhowgts.size()));ix<N;++ix) {
      output+=_rhowgts[ix]*rhoBreitWigner(q2,ix);
    }
  }
  else {
    unsigned int temp(ires);
    if(temp<_rhowgts.size()&&temp<3)
      output=_rhowgts[temp]*rhoBreitWigner(q2,temp);
    else
      output=0.;
  }
  return output/norm;
}
