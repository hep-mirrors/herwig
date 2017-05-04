// -*- C++ -*-
//
// VectorMeson3PionDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson3PionDecayer class.
//

#include "VectorMeson3PionDecayer.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson3PionDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    double temp;
    unsigned int iy;
    for(unsigned int ix=0;ix<_incoming.size();++ix) {
      _maxwgt[ix]=mode(ix)->maxWeight();
      for(iy=0;iy<3;++iy) {
	if(mode(ix)->numberChannels()>3*iy+1) {
	  temp=mode(ix)->channelWeight(3*iy)+mode(ix)->channelWeight(3*iy+1)+
	    mode(ix)->channelWeight(3*iy+2);
	  temp/=3.;
	  switch(iy) {
	  case 0: _rho1wgt[ix]=temp; break;
	  case 1: _rho2wgt[ix]=temp; break;
	  case 2: _rho3wgt[ix]=temp; break;
	  }
	}
      }
    }
  }
}

VectorMeson3PionDecayer::VectorMeson3PionDecayer() 
  : _incoming(2), _coupling(2), _directcoupling(2), _directphase(2),
    _rho2coupling(2), _rho2phase(2), _rho3coupling(2), _rho3phase(2),
    _maxwgt(2), _rho1wgt(2), _rho2wgt(2), _rho3wgt(2), _rho1mass(2),
    _rho2mass(2), _rho3mass(2), _rho1width(2), _rho2width(2), 
    _rho3width(2), _defaultmass(2), _mpic(ZERO), _mpi0(ZERO) {
  // matrix element storage
  // omega decay
  _incoming[0] = 223;
  _coupling[0] = 178.71/GeV;
  _directcoupling[0] = 0.;_directphase[0] = 0.;
  _rho2coupling[0] = 0.;_rho2phase[0] = 0.;
  _rho3coupling[0] = 0.;_rho3phase[0] = 0.;
  _maxwgt[0] = 6.64168;
  _rho1wgt[0] =  1.0;
  _rho2wgt[0] = -1.0; 
  _rho3wgt[0] = -1.0;
  _rho1mass[0] = 0.7758*GeV; 
  _rho2mass[0] = 1.4650*GeV;
  _rho3mass[0] = 1.7000*GeV; 
  _rho1width[0] = 0.1503*GeV; 
  _rho2width[0] = 0.3100*GeV;
  _rho3width[0] = 0.2400*GeV;
  _defaultmass[0] = true;
  // phi decay
  _incoming[1] = 333;
  _coupling[1] = 8.788/GeV;
  _directcoupling[1] = 0.78; _directphase[1] = -2.47;
  _rho2coupling[1] = 0.;_rho2phase[1] = 0.;
  _rho3coupling[1] = 0.;_rho3phase[1] = 0.;
  _maxwgt[1] = 5.62103;
  _rho1wgt[1] =  1.0;
  _rho2wgt[1] = -1.0; 
  _rho3wgt[1] = -1.0;
  _rho1mass[1] = 0.7758*GeV; 
  _rho2mass[1] = 1.4500*GeV;
  _rho3mass[1] = 1.7000*GeV; 
  _rho1width[1] = 0.1439*GeV; 
  _rho2width[1] = 0.3100*GeV;
  _rho3width[1] = 0.2400*GeV;
  _defaultmass[1] = false;
  // initial size of the arrays
  _initsize=_coupling.size();
  // generation of intermediates
  generateIntermediates(true);
}

void VectorMeson3PionDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the consistence of the decay modes
  unsigned int isize=_incoming.size();
  if(isize!=_coupling.size()       || 
     isize!=_directcoupling.size() || isize!=_directphase.size()  ||
     isize!=_rho2coupling.size()   || isize!=_rho2phase.size()    ||
     isize!=_rho3coupling.size()   || isize!=_rho3phase.size()    ||
     isize!=_maxwgt.size()         || isize!=_rho1wgt.size()      ||
     isize!=_rho2wgt.size()        || isize!=_rho3wgt.size()      ||
     isize!=_rho1mass.size()       || isize!=_rho2mass.size()     ||
     isize!=_rho3mass.size()       || isize!=_rho1width.size()    ||
     isize!=_rho2width.size()      || isize!=_rho3width.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "VectorMeson3PionDecayer::doinit()" 
			  << Exception::abortnow;
  // calculate the parameters 
  // set the external particles
  tPDVector extpart(4);
  extpart[1]=getParticleData(ParticleID::pi0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piminus);
  // pointer to the different rho resonances
  // the rho0 resonances
  tPDPtr rho0[3]={getParticleData(113),getParticleData(100113),getParticleData(30113)};
  // the charged rho resonance
  tPDPtr rhom[3]={getParticleData(-213),getParticleData(-100213),
		  getParticleData(-30213)};
  tPDPtr rhop[3]={getParticleData(213),getParticleData(100213),getParticleData(30213)};
  // create the integration channels for the decay
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr newchannel;
  unsigned int iy,iz;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(int(_incoming[ix]));
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    // decide which rho resonances to add
    double  temp[3] = {_rho1wgt[ix]  ,_rho2wgt[ix]  ,_rho3wgt[ix]  };
    Energy  mass[3] = {_rho1mass[ix] ,_rho2mass[ix] ,_rho3mass[ix] };
    Energy width[3] = {_rho1width[ix],_rho2width[ix],_rho3width[ix]};
    vector<double> wgt;
    // set the mass parameters to the default if needed
    if(_defaultmass[ix]) {
      if(rhom[0]) {
	_rho1mass[ix] = rhom[0]->mass(); 
	_rho1width[ix] = rhom[0]->width();
      }
      if(rhom[1]) {
	_rho2mass[ix] = rhom[1]->mass(); 
	_rho2width[ix] = rhom[1]->width();
      }
      if(rhom[2]) {
	_rho3mass[ix] = rhom[2]->mass(); 
	_rho3width[ix] = rhom[2]->width();
      }
    }
    double sumwgt(0);
    for(iy=0;iy<3;++iy) {
      if(temp[iy]>0) sumwgt+=temp[iy];
    }
    for(iy=0;iy<3;++iy) {
      if(temp[iy]>0) {
	// set the weights for the channels
	for(iz=0;iz<3;++iz) wgt.push_back(temp[iy]/3./sumwgt);
	// rho0 channel
	newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	newchannel->addIntermediate(extpart[0],0,0.0,-1,1);
	newchannel->addIntermediate(rho0[iy]  ,0,0.0, 2,3);
	mode->addChannel(newchannel);
	if(!_defaultmass[ix]) resetIntermediate(rho0[iy],mass[iy],width[iy]);
	// rho+ channel
	newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	newchannel->addIntermediate(extpart[0],0,0.0,-1,3);
	newchannel->addIntermediate(rhop[iy]  ,0,0.0, 1,2);
	mode->addChannel(newchannel);
	if(!_defaultmass[ix]) resetIntermediate(rhop[iy],mass[iy],width[iy]);
	// rho- channel
	newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	newchannel->addIntermediate(extpart[0],0,0.0,-1,2);
	newchannel->addIntermediate(rhom[iy]  ,0,0.0, 1,3);
	mode->addChannel(newchannel);
	if(!_defaultmass[ix]) mode->resetIntermediate(rhom[iy],mass[iy],width[iy]);
	addMode(mode,_maxwgt[ix],wgt);
      }
    }
  }
  // work out the masses and constants for the running widths
  Energy pcm;
  _mpi0=getParticleData(ParticleID::pi0)->mass();
  _mpic=getParticleData(ParticleID::piplus)->mass();
  Complex ii(0.,1.);
  _rhomass.resize(_incoming.size()); _rhomass2.resize(_incoming.size());
  _rho0const.resize(_incoming.size()); _rhocconst.resize( _incoming.size());
  _ccoupling.resize(_incoming.size());
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    // set the masses of the rho resonances
    _rhomass[ix].push_back(_rho1mass[ix]);
    _rhomass[ix].push_back(_rho2mass[ix]);
    _rhomass[ix].push_back(_rho3mass[ix]);
    _rhomass2[ix].push_back(_rho1mass[ix]*_rho1mass[ix]);
    _rhomass2[ix].push_back(_rho2mass[ix]*_rho2mass[ix]);
    _rhomass2[ix].push_back(_rho3mass[ix]*_rho3mass[ix]);
    // set up the constants for the running width
    pcm=Kinematics::pstarTwoBodyDecay(_rho1mass[ix],_mpic,_mpic);
    _rho0const[ix].push_back(_rho1mass[ix]*_rho1mass[ix]*_rho1width[ix]/(pcm*pcm*pcm));
    pcm=Kinematics::pstarTwoBodyDecay(_rho2mass[ix],_mpic,_mpic);
    _rho0const[ix].push_back(_rho2mass[ix]*_rho2mass[ix]*_rho2width[ix]/(pcm*pcm*pcm));
    pcm=Kinematics::pstarTwoBodyDecay(_rho3mass[ix],_mpic,_mpic);
    _rho0const[ix].push_back(_rho3mass[ix]*_rho3mass[ix]*_rho3width[ix]/(pcm*pcm*pcm));
    pcm=Kinematics::pstarTwoBodyDecay(_rho1mass[ix],_mpi0,_mpic);
    _rhocconst[ix].push_back(_rho1mass[ix]*_rho1mass[ix]*_rho1width[ix]/(pcm*pcm*pcm));
    pcm=Kinematics::pstarTwoBodyDecay(_rho2mass[ix],_mpi0,_mpic);
    _rhocconst[ix].push_back(_rho2mass[ix]*_rho2mass[ix]*_rho2width[ix]/(pcm*pcm*pcm));
    pcm=Kinematics::pstarTwoBodyDecay(_rho3mass[ix],_mpi0,_mpic);
    _rhocconst[ix].push_back(_rho3mass[ix]*_rho3mass[ix]*_rho3width[ix]/(pcm*pcm*pcm));
    // set the complex coupling constants
    _ccoupling[ix].push_back(1./_rhomass2[ix][0]);
    _ccoupling[ix].push_back( _rho2coupling[ix]/_rhomass2[ix][0]*
			      (cos(_rho2phase[ix])+ii*sin(_rho2phase[ix])));
    _ccoupling[ix].push_back( _rho3coupling[ix]/_rhomass2[ix][0]*
			      (cos(_rho3phase[ix])+ii*sin(_rho3phase[ix])));
    _ccoupling[ix].push_back(_directcoupling[ix]/_rhomass2[ix][0]*
			     (cos(_directphase[ix])+ii*sin(_directphase[ix])));
  }
}

int VectorMeson3PionDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  cc=false;
  // must be three outgoing particles
  if(children.size()!=3) return -1;
  // check the id's of the outgoing particles
  int id;
  unsigned int npi0(0),npip(0),npim(0);
  tPDVector::const_iterator pit  = children.begin();
  for(;pit!=children.end();++pit) {
    id = (*pit)->id();
    if(id==ParticleID::pi0)          ++npi0;
    else if(id==ParticleID::piplus)  ++npip;
    else if(id==ParticleID::piminus) ++npim;
  }
  if(!(npi0==1&&npip==1&&npim==1)) return -1;
  unsigned int ix(0);
  id=parent->id();
  int imode(-1);
  do {
    if(_incoming[ix]==id) imode=ix;
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  return imode;
}

void VectorMeson3PionDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << ounit(_coupling,1/MeV) << _directcoupling << _rho2coupling 
     << _rho3coupling << _directphase << _rho2phase << _rho3phase 
     << _maxwgt << _rho1wgt << _rho2wgt << _rho3wgt << ounit(_rho1mass,GeV) 
     << ounit(_rho2mass,GeV) << ounit(_rho3mass,GeV) << ounit(_rho1width,GeV) 
     << ounit(_rho2width,GeV) << ounit(_rho3width,GeV) << _defaultmass 
     << _rho0const << _rhocconst << ounit(_rhomass,GeV) << ounit(_rhomass2,GeV2) 
     << ounit(_ccoupling,1/MeV2) << ounit(_mpi0,MeV) << ounit(_mpic,MeV);
}

void VectorMeson3PionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> iunit(_coupling,1/MeV) >> _directcoupling >> _rho2coupling 
     >> _rho3coupling >> _directphase >> _rho2phase >> _rho3phase 
     >> _maxwgt >> _rho1wgt >> _rho2wgt >> _rho3wgt >> iunit(_rho1mass,GeV) 
     >> iunit(_rho2mass,GeV) >> iunit(_rho3mass,GeV) >> iunit(_rho1width,GeV) 
     >> iunit(_rho2width,GeV) >> iunit(_rho3width,GeV) >> _defaultmass
     >> _rho0const >> _rhocconst >> iunit(_rhomass,GeV) >> iunit(_rhomass2,GeV2) 
     >> iunit(_ccoupling,1/MeV2) >> iunit(_mpi0,MeV) >> iunit(_mpic,MeV);
}

ClassDescription<VectorMeson3PionDecayer> 
VectorMeson3PionDecayer::initVectorMeson3PionDecayer;
// Definition of the static class description member.

void VectorMeson3PionDecayer::Init() {

  static ClassDocumentation<VectorMeson3PionDecayer> documentation
    ("The VectorMeson3PionDecayer class is designed for the decay "
     "of I=0 vector mesons to three pions via a current taking into account the "
     "rho and a possible direct term",
     "The decay of I=0 vector mesons to three pions via a current taking into account the "
     "rho and a possible direct term is taken from \\cite{Aloisio:2003ur}.",
     "%\\cite{Aloisio:2003ur}\n"
     "\\bibitem{Aloisio:2003ur}\n"
     "  A.~Aloisio {\\it et al.}  [KLOE Collaboration],\n"
     "  %``Study of the decay Phi --> pi+ pi- pi0 with the KLOE detector,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 561}, 55 (2003)\n"
     "  [Erratum-ibid.\\  B {\\bf 609}, 449 (2005)]\n"
     "  [arXiv:hep-ex/0303016].\n"
     "  %%CITATION = PHLTA,B561,55;%%\n"
     );
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson3PionDecayer::_incoming,
     0, 0, 0, 0, 1000000, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay, this is the coupling of the decaying "
     "particle to the lowest lying rho multiplet.",
     &VectorMeson3PionDecayer::_coupling,
     1./GeV, -1, 1./GeV,-1000./GeV, 1000./GeV, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceDirectCoupling
    ("DirectCoupling",
     "The magnitude of the coupling of the direct term with respect to the "
     "coupling of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_directcoupling,
     0, 0, 0, -100000., 100000., false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Coupling
    ("Rho2Coupling",
     "The magnitude of the coupling of the second rho multiplet with respect to the "
     "lowest lying  multiplet",
     &VectorMeson3PionDecayer::_rho2coupling,
     0, 0, 0, -10., 10., false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Coupling
    ("Rho3Coupling",
     "The magntiude of the coupling of the third rho multiplet with respect to the "
     "lowest lying multiplet",
     &VectorMeson3PionDecayer::_rho3coupling,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceDirectPhase
    ("DirectPhase",
     "The phase of the coupling of the direct term with respect to the "
     "coupling of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_directphase,
     0, 0, 0, -Constants::twopi, Constants::twopi, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Phase
    ("Rho2Phase",
     "The phasee of the coupling of the second rho multiplet with respect to the "
     "lowest lying  multiplet",
     &VectorMeson3PionDecayer::_rho2phase,
     0, 0, 0, -Constants::twopi, Constants::twopi, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Phase
    ("Rho3Phase",
     "The phase of the coupling of the third rho multiplet with respect to the "
     "lowest lying multiplet",
     &VectorMeson3PionDecayer::_rho3phase,
     0, 0, 0, -Constants::twopi, Constants::twopi, false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceMaxWgt
    ("MaxWeight",
     "The maximum weight for the integration of the channel",
     &VectorMeson3PionDecayer::_maxwgt,
     0, 0, 0, 0., 100., false, false, true);
  
  static ParVector<VectorMeson3PionDecayer,double> interfaceRho1Wgt
    ("Rho1Weight",
     "The weight for the lowest lying rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho1wgt,
     0, 0, 0, 0, 1., false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho2Wgt
    ("Rho2Weight",
     "The weight for the second rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho2wgt,
     0, 0, 0, -2., 1., false, false, true);

  static ParVector<VectorMeson3PionDecayer,double> interfaceRho3Wgt
    ("Rho3Weight",
     "The weight for the third rho multiplet's in the integration",
     &VectorMeson3PionDecayer::_rho3wgt,
     0, 0, 0, -2., 1., false, false, true);

  static ParVector<VectorMeson3PionDecayer,Energy> interfaceRho1Mass
    ("Rho1Mass",
     "The mass of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_rho1mass,
     GeV, -1, 0.77*GeV, ZERO, 5.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,Energy> interfaceRho2Mass
    ("Rho2Mass",
     "The mass of the second rho multiplet",
     &VectorMeson3PionDecayer::_rho2mass,
     GeV, -1, 0.77*GeV, ZERO, 5.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,Energy> interfaceRho3Mass
    ("Rho3Mass",
     "The mass of the third rho multiplet",
     &VectorMeson3PionDecayer::_rho3mass,
     GeV, -1, 0.77*GeV, ZERO, 5.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,Energy> interfaceRho1Width
    ("Rho1Width",
     "The width of the lowest lying rho multiplet",
     &VectorMeson3PionDecayer::_rho1width,
     GeV, -1, 0.15*GeV, ZERO, 1.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,Energy> interfaceRho2Width
    ("Rho2Width",
     "The width of the second rho multiplet",
     &VectorMeson3PionDecayer::_rho2width,
     GeV, -1, 0.15*GeV, ZERO, 1.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,Energy> interfaceRho3Width
    ("Rho3Width",
     "The width of the third rho multiplet",
     &VectorMeson3PionDecayer::_rho3width,
     GeV, -1, 0.15*GeV, ZERO, 1.*GeV, false, false, true);

  static ParVector<VectorMeson3PionDecayer,bool> interfaceDefaultParam
    ("DefaultParameters",
     "If true the default rho masses and widths are used, otherwise "
     "the values specified in the rhomass and width arrays are used,",
     &VectorMeson3PionDecayer::_defaultmass,
     0, 0, 1, 0, 1, false, false, true);

}

double VectorMeson3PionDecayer::me2(const int ichan,
				    const Particle & inpart,
				    const ParticleVector & decay,
				    MEOption meopt) const {
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
  // compute the matrix element
  // work out the prefactor
  complex<InvEnergy2> pre(ZERO);
  Complex resfact,ii(0.,1.);
  if(ichan<0){pre=_ccoupling[imode()][3];}
  Energy pcm;
  // work out the direct invariant masses needed
  Energy mrho0(sqrt(decay[1]->momentum().m2(decay[2]->momentum())));
  Energy mrhop(sqrt(decay[1]->momentum().m2(decay[0]->momentum())));
  Energy mrhom(sqrt(decay[2]->momentum().m2(decay[0]->momentum())));
  // contribution of the resonances
  int ichannow(-3);
  for(unsigned int ix=0;ix<3;++ix) {
    ichannow+=3;
    if((ix==0 && _rho1wgt[imode()]>0.) || (ix==1 && _rho2wgt[imode()]>0.) ||
       (ix==2 && _rho3wgt[imode()]>0.)) {
      if(ichan<0) {
	// rho0 contribution
	pcm = Kinematics::pstarTwoBodyDecay(mrho0,_mpic,_mpic);
	resfact = _rhomass2[imode()][ix]/
	  (mrho0*mrho0-_rhomass2[imode()][ix]
	   +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
	// rho+ contribution
	pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
	resfact+= _rhomass2[imode()][ix]/
	  (mrhop*mrhop-_rhomass2[imode()][ix]
	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
	// rho- contribution
	pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
	resfact+= _rhomass2[imode()][ix]/
	  (mrhom*mrhom-_rhomass2[imode()][ix]
	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhom);
	// add the contribution
      }
      else if(ichan==ichannow) {
	pcm = Kinematics::pstarTwoBodyDecay(mrho0,_mpic,_mpic);
	resfact = _rhomass2[imode()][ix]/
	  (mrho0*mrho0-_rhomass2[imode()][ix]
	   +ii*pcm*pcm*pcm*_rho0const[imode()][ix]/mrho0);
      }
      else if(ichan==ichannow+1) {
	pcm = Kinematics::pstarTwoBodyDecay(mrhop,_mpic,_mpi0);
	resfact+= _rhomass2[imode()][ix]/
	  (mrhop*mrhop-_rhomass2[imode()][ix]
	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhop);
      }
      else if(ichan==ichannow+2) {
	pcm = Kinematics::pstarTwoBodyDecay(mrhom,_mpic,_mpi0);
	resfact+= _rhomass2[imode()][ix]/
	  (mrhom*mrhom-_rhomass2[imode()][ix]
	   +ii*pcm*pcm*pcm*_rhocconst[imode()][ix]/mrhom);
      }
      pre += resfact * _ccoupling[imode()][ix];
      ichannow+=3;
    }
  }
  // polarization vector piece
  LorentzPolarizationVector 
    scalar=_coupling[imode()]*pre*epsilon(decay[0]->momentum(),
					  decay[1]->momentum(),
					  decay[2]->momentum());
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix) 
    (*ME())(ix,0,0,0)=scalar.dot(_vectors[ix]);
  // return the answer
  return ME()->contract(_rho).real();
}

double VectorMeson3PionDecayer::
threeBodyMatrixElement(const int imode, const Energy2 q2,
		       const  Energy2 s3, const Energy2 s2, const Energy2 s1, const 
		       Energy , const Energy , const Energy ) const {
  Lorentz5Momentum p1,p2,p3; Energy2 ee1,ee2,ee3;Energy pp1,pp2,pp3;
  Energy q(sqrt(q2));
  Energy2 mpi2c(_mpic*_mpic),mpi20(_mpi0*_mpi0);
  p1.setE(0.5*(q2+mpi20-s1)/q); ee1=p1.e()*p1.e(); pp1=sqrt(ee1-mpi20);
  p2.setE(0.5*(q2+mpi2c-s2)/q); ee2=p2.e()*p2.e(); pp2=sqrt(ee2-mpi2c);
  p3.setE(0.5*(q2+mpi2c-s3)/q); ee3=p3.e()*p3.e(); pp3=sqrt(ee3-mpi2c);
  // take momentum of 1 parallel to z axis
  p1.setX(ZERO);p1.setY(ZERO);p1.setZ(pp1);
  // construct 2 
  double cos2(0.5*(ee1+ee2-ee3-mpi20)/pp1/pp2);
  p2.setX(pp2*sqrt(1.-cos2*cos2)); p2.setY(ZERO); p2.setZ(-pp2*cos2);
  // construct 3
  double cos3(0.5*(ee1-ee2+ee3-mpi20)/pp1/pp3);
  p3.setX(-pp3*sqrt(1.-cos3*cos3)); p3.setY(ZERO); p3.setZ(-pp3*cos3); 
  // compute the prefactor
  complex<InvEnergy2> pre(_ccoupling[imode][3]);
  Complex resfact,ii(0.,1.);
  // rho0 contribution
  Energy pcm,mrho1(sqrt(s1)),mrho2(sqrt(s2)),mrho3(sqrt(s3));
  for(unsigned int ix=0;ix<3;++ix) {
    // rho0 contribution
    pcm = Kinematics::pstarTwoBodyDecay(mrho1,_mpic,_mpic);
    resfact = _rhomass2[imode][ix]/(mrho1*mrho1-_rhomass2[imode][ix]
				    +ii*pcm*pcm*pcm*_rho0const[imode][ix]/mrho1);
    // rho+ contribution
    pcm = Kinematics::pstarTwoBodyDecay(mrho2,_mpic,_mpi0);
    resfact+= _rhomass2[imode][ix]/(mrho2*mrho3-_rhomass2[imode][ix]
				    +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho2);
    // rho- contribution
    pcm = Kinematics::pstarTwoBodyDecay(mrho3,_mpic,_mpi0);
    resfact+= _rhomass2[imode][ix]/(mrho3*mrho3-_rhomass2[imode][ix]
				    +ii*pcm*pcm*pcm*_rhocconst[imode][ix]/mrho3);
    // add the contribution
    pre+=resfact *_ccoupling[imode][ix];
  }
  LorentzPolarizationVector current =
    _coupling[imode]*(pre*epsilon(p1,p2,p3));
  Complex temp(current.dot(current.conjugate()));
  return -temp.real()/3.;
} 

WidthCalculatorBasePtr 
VectorMeson3PionDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int imode(-1),id(dm.parent()->id());
  unsigned int ix=0;
  do {
    if(_incoming[ix]==id) imode=ix; 
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  // construct the integrator
  vector<double> inweights(3,1./3.);
  vector<int> intype;intype.push_back(1);intype.push_back(2);intype.push_back(3);
  Energy mrho(getParticleData(ParticleID::rhoplus)->mass());
  Energy wrho(getParticleData(ParticleID::rhoplus)->width());
  vector<Energy> inmass(3,mrho);
  vector<Energy> inwidth(3,wrho);
  vector<double> inpow(2,0.0);
  //tcDecayIntegratorPtr decayer(this);
  WidthCalculatorBasePtr output(
    new_ptr(ThreeBodyAllOnCalculator<VectorMeson3PionDecayer>
	    (inweights,intype,inmass,inwidth,inpow,
	     *this,imode,_mpi0,_mpic,_mpic)));
  return output;
}

void VectorMeson3PionDecayer::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " 
	     << ix << " " << _incoming[ix] << endl;
      output << "newdef " << name() << ":Coupling " 
	     << ix << " " << _coupling[ix]*GeV << endl;
      output << "newdef " << name() << ":DirectCoupling " 
	     << ix << " " << _directcoupling[ix] << endl;
      output << "newdef " << name() << ":Rho2Coupling " 
	     << ix << " " << _rho2coupling[ix] << endl;
      output << "newdef " << name() << ":Rho3Coupling " 
	     << ix << " " << _rho3coupling[ix] << endl;
      output << "newdef " << name() << ":DirectPhase " 
	     << ix << " " << _directphase[ix] << endl;
      output << "newdef " << name() << ":Rho2Phase " 
	     << ix << " " << _rho2phase[ix] << endl;
      output << "newdef " << name() << ":Rho3Phase " 
	     << ix << " " << _rho3phase[ix] << endl;
      output << "newdef " << name() << ":MaxWeight " 
	     << ix << " " << _maxwgt[ix] << endl;
      output << "newdef " << name() << ":Rho1Weight " 
	     << ix << " " << _rho1wgt[ix] << endl;
      output << "newdef " << name() << ":Rho2Weight " 
	     << ix << " " << _rho2wgt[ix] << endl;
      output << "newdef " << name() << ":Rho3Weight " 
	     << ix << " " << _rho3wgt[ix] << endl;
      output << "newdef " << name() << ":Rho1Mass " 
	     << ix << " " << _rho1mass[ix]/GeV << endl;
      output << "newdef " << name() << ":Rho2Mass " 
	     << ix << " " << _rho2mass[ix]/GeV<< endl;
      output << "newdef " << name() << ":Rho3Mass " 
	     << ix << " " << _rho3mass[ix]/GeV<< endl;
      output << "newdef " << name() << ":Rho1Width " 
	     << ix << " " << _rho1width[ix]/GeV << endl;
      output << "newdef " << name() << ":Rho2Width " 
	     << ix << " " << _rho2width[ix]/GeV << endl;
      output << "newdef " << name() << ":Rho3Width " 
	     << ix << " " << _rho3width[ix]/GeV << endl;
      output << "newdef " << name() << ":DefaultParameters " 
	     << ix << " " << _defaultmass[ix] << endl;
    }
    else {
      output << "insert " << name() << ":Incoming " 
	     << ix << " " << _incoming[ix] << endl;
      output << "insert " << name() << ":Coupling " 
	     << ix << " " << _coupling[ix]*GeV << endl;
      output << "insert " << name() << ":DirectCoupling " 
	     << ix << " " << _directcoupling[ix] << endl;
      output << "insert " << name() << ":Rho2Coupling " 
	     << ix << " " << _rho2coupling[ix] << endl;
      output << "insert " << name() << ":Rho3Coupling " 
	     << ix << " " << _rho3coupling[ix] << endl;
      output << "insert " << name() << ":DirectPhase " 
	     << ix << " " << _directphase[ix] << endl;
      output << "insert " << name() << ":Rho2Phase " 
	     << ix << " " << _rho2phase[ix] << endl;
      output << "insert " << name() << ":Rho3Phase " 
	     << ix << " " << _rho3phase[ix] << endl;
      output << "insert " << name() << ":MaxWeight " 
	     << ix << " " << _maxwgt[ix] << endl;
      output << "insert " << name() << ":Rho1Weight " 
	     << ix << " " << _rho1wgt[ix] << endl;
      output << "insert " << name() << ":Rho2Weight " 
	     << ix << " " << _rho2wgt[ix] << endl;
      output << "insert " << name() << ":Rho3Weight " 
	     << ix << " " << _rho3wgt[ix] << endl;
      output << "insert " << name() << ":Rho1Mass " 
	     << ix << " " << _rho1mass[ix]/GeV << endl;
      output << "insert " << name() << ":Rho2Mass " 
	     << ix << " " << _rho2mass[ix]/GeV<< endl;
      output << "insert " << name() << ":Rho3Mass " 
	     << ix << " " << _rho3mass[ix]/GeV<< endl;
      output << "insert " << name() << ":Rho1Width " 
	     << ix << " " << _rho1width[ix]/GeV << endl;
      output << "insert " << name() << ":Rho2Width " 
	     << ix << " " << _rho2width[ix]/GeV << endl;
      output << "insert " << name() << ":Rho3Width " 
	     << ix << " " << _rho3width[ix]/GeV << endl;
      output << "insert " << name() << ":DefaultParameters " 
	     << ix << " " << _defaultmass[ix] << endl;
    }
  }
  if(header){output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;}
}
