// -*- C++ -*-
//
// TwoMesonRhoKStarCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoMesonRhoKStarCurrent class.
//
//  Author: Peter Richardson
//

#include "TwoMesonRhoKStarCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

TwoMesonRhoKStarCurrent::TwoMesonRhoKStarCurrent() {
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-1);
  addDecayMode(2,-3);
  setInitialModes(5);
  // the weights of the different resonances in the matrix elements
  _pimag.push_back(  1.0);_pimag.push_back(  0.167);_pimag.push_back(  0.05);
  _piphase.push_back(0.0);_piphase.push_back(  180);_piphase.push_back(0.0);
  _kmag.push_back(  1.0);_kmag.push_back(  0.038);_kmag.push_back(  0.00);
  _kphase.push_back(0.0);_kphase.push_back(  180);_kphase.push_back(0.0);
  // models to use
  _pimodel = 0;_kmodel=0;
  // parameter for the masses (use the parameters freom the CLEO fit 
  // rather than the PDG masses etc)
  _rhoparameters=true;
  _rhomasses.push_back(774.6*MeV);_rhomasses.push_back(1408*MeV);
  _rhomasses.push_back(1700*MeV);
  _rhowidths.push_back(149*MeV);_rhowidths.push_back(502*MeV);
  _rhowidths.push_back(235*MeV);
  _kstarparameters=true;
  _kstarmasses.push_back(0.8921*GeV);_kstarmasses.push_back(1.700*GeV);
  _kstarwidths.push_back(0.0513*GeV);_kstarwidths.push_back(0.235*GeV);
}

void TwoMesonRhoKStarCurrent::doinit() {
  WeakDecayCurrent::doinit();
  // check consistency of parametrers
  if(_rhomasses.size()!=_rhowidths.size()||
     _kstarmasses.size()!=_kstarwidths.size()) {
    throw InitException() << "Inconsistent parameters in TwoMesonRhoKStarCurrent"
			  << "::doinit()" << Exception::abortnow;
  }
  // the resonances
  tPDPtr res[6]={getParticleData(-213   ),getParticleData(-100213),
		 getParticleData(-30213 ),getParticleData(-323   ),
		 getParticleData(-100323),getParticleData(-30323 )};
  // reset the masses in the form-factors if needed
  if(_rhoparameters&&_rhomasses.size()<3) {
    for(unsigned int ix=_rhomasses.size();ix<3;++ix) {
      if(res[ix]) _rhomasses.push_back(res[ix]->mass() );
      if(res[ix]) _rhowidths.push_back(res[ix]->width());
    }
  }
  else if(!_rhoparameters) {
    _rhomasses.clear();_rhowidths.clear();
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]) _rhomasses.push_back(res[ix]->mass() );
      if(res[ix]) _rhowidths.push_back(res[ix]->width());
    }
  }
  // then the Kstar resonances
  if(_kstarparameters&&_kstarmasses.size()<3) {
    for(unsigned int ix=_kstarmasses.size();ix<3;++ix) {
      if(res[ix+3]) _kstarmasses.push_back(res[ix+3]->mass());
      if(res[ix+3]) _kstarwidths.push_back(res[ix+3]->width());
    }
  }
  else if(!_kstarparameters) {
    _kstarmasses.clear();_kstarwidths.clear();
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix+3]) _kstarmasses.push_back(res[ix+3]->mass());
      if(res[ix+3]) _kstarwidths.push_back(res[ix+3]->width());
    }
  }
  // set up for the Breit Wigners
  Energy mpi0(   getParticleData(ParticleID::pi0   )->mass());
  Energy mpiplus(getParticleData(ParticleID::piplus)->mass());
  Energy mk0(    getParticleData(ParticleID::K0    )->mass());
  // rho resonances
  for(unsigned int ix=0;ix<3;++ix) {
    _mass.push_back(_rhomasses[ix]);
    _width.push_back(_rhowidths[ix]);
    _mass2.push_back(_mass[ix]*_mass[ix]);
    _massw.push_back(_mass[ix]*_width[ix]);
    _massa.push_back(mpi0);
    _massb.push_back(mpiplus);
    _mom.push_back(pcm(ix,_mass[ix]));
    _hm2.push_back(GSModelhFunction(ix,_mass[ix]));
    _dparam.push_back(GSModelDParameter(ix));
    _dhdq2.push_back(GSModeldhdq2Parameter(ix));
  }
  // Kstar resonances
  for(unsigned int ix=0;ix<3;++ix) {
    _mass.push_back(_kstarmasses[ix]);
    _width.push_back(_kstarwidths[ix]);
    _mass2.push_back(_mass[ix+3]*_mass[ix+3]);
    _massw.push_back(_mass[ix+3]*_width[ix+3]);
    _massa.push_back(mk0);
    _massb.push_back(mpiplus);
    _mom.push_back(pcm(ix+3,_mass[ix+3]));
    _hm2.push_back(GSModelhFunction(ix+3,_mass[ix+3]));
    _dparam.push_back(GSModelDParameter(ix+3));
    _dhdq2.push_back(GSModeldhdq2Parameter(ix+3));
  }
  // weights for the rho channels
  if(_pimag.size()!=_piphase.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
			  << " rho channel must be the same size in "
			  << "TwoMesonRhoKStarCurrent::doinit()" << Exception::runerror;
  _piwgt.resize(_pimag.size());
  for(unsigned int ix=0;ix<_pimag.size();++ix) {
    double angle = _piphase[ix]/180.*Constants::pi;
    _piwgt[ix] = _pimag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
  // weights for the K* channels
  if(_kmag.size()!=_kphase.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
			  << " K* channel must be the same size in "
			  << "TwoMesonRhoKStarCurrent::doinit()" << Exception::runerror;
  _kwgt.resize(_kmag.size());
  for(unsigned int ix=0;ix<_kmag.size();++ix) {
    double angle = _kphase[ix]/180.*Constants::pi;
    _kwgt[ix] = _kmag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
}

void TwoMesonRhoKStarCurrent::persistentOutput(PersistentOStream & os) const {
  os << _pimodel << _kmodel << _piwgt << _pimag << _piphase << _kwgt << _kmag 
     << _kphase << _rhoparameters << _kstarparameters << ounit(_rhomasses,GeV) << ounit(_rhowidths,GeV) 
     << ounit(_kstarmasses,GeV) << ounit(_kstarwidths,GeV) 
     << ounit(_mass,GeV) << ounit(_width,GeV) << ounit(_mass2,GeV2) << ounit(_massw,GeV2) 
     << ounit(_massa,GeV) <<ounit(_massb,GeV) << ounit(_mom,GeV) << ounit(_dhdq2,1/GeV2) 
     << _hm2 << _dparam;
}

void TwoMesonRhoKStarCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _pimodel >> _kmodel >> _piwgt >> _pimag >> _piphase >> _kwgt >> _kmag 
     >> _kphase >> _rhoparameters >> _kstarparameters >> iunit(_rhomasses,GeV) >> iunit(_rhowidths,GeV) 
     >> iunit(_kstarmasses,GeV) >> iunit(_kstarwidths,GeV) 
     >> iunit(_mass,GeV) >> iunit(_width,GeV) >> iunit(_mass2,GeV2) >> iunit(_massw,GeV2) 
     >> iunit(_massa,GeV) >> iunit(_massb,GeV) >> iunit(_mom,GeV) >> iunit(_dhdq2,1/GeV2) 
     >> _hm2 >> _dparam;
}

ClassDescription<TwoMesonRhoKStarCurrent> TwoMesonRhoKStarCurrent::initTwoMesonRhoKStarCurrent;
// Definition of the static class description member.

void TwoMesonRhoKStarCurrent::Init() {

  static ParVector<TwoMesonRhoKStarCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_rhomasses, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoMesonRhoKStarCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_rhowidths, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoMesonRhoKStarCurrent,Energy> interfaceKstarMasses
    ("KstarMasses",
     "The masses of the different K* resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_kstarmasses, MeV, -1, 891.66*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoMesonRhoKStarCurrent,Energy> interfaceKstarWidths
    ("KstarWidths",
     "The widths of the different K* resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_kstarwidths, MeV, -1, 50.8*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static Switch<TwoMesonRhoKStarCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values for the rho meson masses and widths",
     &TwoMesonRhoKStarCurrent::_rhoparameters, true, false, false);
  static SwitchOption interfaceRhoParameterstrue
    (interfaceRhoParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceRhoParametersParticleData
    (interfaceRhoParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static Switch<TwoMesonRhoKStarCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values for the Kstar meson masses and widths",
     &TwoMesonRhoKStarCurrent::_kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);
  
  static ParVector<TwoMesonRhoKStarCurrent,double> interfacePiMagnitude
    ("PiMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_pimag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoMesonRhoKStarCurrent,double> interfacePiPhase
    ("PiPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoMesonRhoKStarCurrent::_piphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoMesonRhoKStarCurrent,double> interfaceKMagnitude
    ("KMagnitude",
     "Magnitude of the weight of the different resonances for the K pi channel",
     &TwoMesonRhoKStarCurrent::_kmag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoMesonRhoKStarCurrent,double> interfaceKPhase
    ("KPhase",
     "Phase of the weight of the different resonances for the K pi channel",
     &TwoMesonRhoKStarCurrent::_kphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static Switch<TwoMesonRhoKStarCurrent,int> interfacePiModel
    ("PiModel",
     "The model to use for the propagator for the pion modes.",
     &TwoMesonRhoKStarCurrent::_pimodel, 0, false, false);
  static SwitchOption interfacePiModelKuhn
    (interfacePiModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfacePiModelGounaris
    (interfacePiModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);
  
  static Switch<TwoMesonRhoKStarCurrent,int> interfaceKModel
    ("KModel",
     "The model to use for the propagator for the kaon modes.",
     &TwoMesonRhoKStarCurrent::_kmodel, 0, false, false);
  static SwitchOption interfaceKModelKuhn
    (interfaceKModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfaceKModelGounaris
    (interfaceKModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);

  static ClassDocumentation<TwoMesonRhoKStarCurrent> documentation
    ("The TwoMesonRhoKStarCurrent class is designed to implement weak"
     "decay to two scalar mesons using the models of either Kuhn and "
     "Santamaria (Z. Phys. C48, 445 (1990)) or Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968).  The mixing parameters are taken from "
     "Phys. Rev. D61:112002,2000 (CLEO), although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.",
     "The weak "
     "decay current to two scalar mesons is implemented "
     "using the models of either Kuhn and "
     "Santamaria \\cite{Kuhn:1990ad} or Gounaris and Sakurai \\cite{Gounaris:1968mw}. "
     "The mixing parameters are taken from "
     "\\cite{Asner:1999kj}, although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.",
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Gounaris:1968mw}\n"
     "\\bibitem{Gounaris:1968mw}\n"
     "  G.~J.~Gounaris and J.~J.~Sakurai,\n"
     "   ``Finite width corrections to the vector meson dominance prediction for rho\n"
     "  %$\\to$ e+ e-,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 21}, 244 (1968).\n"
     "  %%CITATION = PRLTA,21,244;%%\n"
     "%\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );



}

// complete the construction of the decay mode for integration
bool TwoMesonRhoKStarCurrent::createMode(int icharge, unsigned int imode,
					 DecayPhaseSpaceModePtr mode,
					 unsigned int iloc,unsigned int,
					 DecayPhaseSpaceChannelPtr phase,Energy upp) {
  if(abs(icharge)!=3) return false; 
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==2) {
    part[0]=getParticleData(ParticleID::K0);
    part[1]=getParticleData(ParticleID::piplus);
  }
  else if(imode==3) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::K0);
  }
  else if(imode==4) {
    part[0]=getParticleData(ParticleID::eta);
    part[1]=getParticleData(ParticleID::Kplus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  DecayPhaseSpaceChannelPtr newchannel;
  // set the resonances
  // two pion or  K+ K0 decay
  tPDPtr res[3];
  if(imode==0||imode==3) {
    res[0]=getParticleData(213);
    res[1]=getParticleData(100213);
    res[2]=getParticleData(30213);
  }
  // K+ pi0 or K0 pi+ or K eta decay
  else if(imode==1||imode==2||imode==4) {
    res[0]=getParticleData(323);
    res[1]=getParticleData(100323);
    res[2]=getParticleData(30323);
  }
  else {
    throw Exception() << "Failure of initialisation in TwoMesonRhoKStarCurrent" 
		      << Exception::abortnow;
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(res[ix]) {
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(res[ix],0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
    }
  }
  // reset the masses in the intergrators if needed
  // for the rho 
  if(_rhoparameters&&(imode==0||imode==3)) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_rhomasses.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_rhomasses[ix],_rhowidths[ix]);
      }
    }
  }
  // for the K*
  else if(_kstarparameters&&imode!=0&&imode!=3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_kstarmasses.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_kstarmasses[ix],_kstarwidths[ix]);
      }
    }
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector TwoMesonRhoKStarCurrent::particles(int icharge, unsigned int imode,
					    int,int) {
  if(abs(icharge)!=3) return tPDVector();
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==2) {
    output[0]=getParticleData(ParticleID::K0);
    output[1]=getParticleData(ParticleID::piplus);
  }
  else if(imode==3) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::Kbar0);
  }
  else {
    output[0]=getParticleData(ParticleID::eta);
    output[1]=getParticleData(ParticleID::Kplus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
TwoMesonRhoKStarCurrent::current(const int imode, const int ichan,
				 Energy & scale,const ParticleVector & outpart,
				 DecayIntegrator::MEOption meopt) const {
  useMe();
  if(meopt==DecayIntegrator::Terminate) {
    for(unsigned int ix=0;ix<2;++ix)
      ScalarWaveFunction::constructSpinInfo(outpart[ix],outgoing,true);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(outpart[0]->momentum()-outpart[1]->momentum());
  Lorentz5Momentum psum (outpart[0]->momentum()+outpart[1]->momentum());
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  LorentzPolarizationVector vect;
  // calculate the current
  Complex FPI(0.),denom(0.);
  // rho
  if(imode==0||imode==3) {
    if(ichan<0) {
      for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix) {
	FPI+=_piwgt[ix]*BreitWigner(q2,_pimodel,0,ix);
	denom+=_piwgt[ix];
      }
    }
    else if(ichan<int(_piwgt.size())&&ichan<3) {
      FPI=_piwgt[ichan]*BreitWigner(q2,_pimodel,0,ichan);
      for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix) denom+=_piwgt[ix];
    }
  }
  // K*
  else {
    if(ichan<0) {
      for(unsigned int ix=0;ix<_kwgt.size()&&ix<3;++ix) {
	FPI+=_kwgt[ix]*BreitWigner(q2,_kmodel,1,ix);
	denom+=_kwgt[ix];
      }
    }
    else if (ichan<int(_kwgt.size())&&ichan<3) {
      FPI=_kwgt[ichan]*BreitWigner(q2,_kmodel,1,ichan);
      for(unsigned int ix=0;ix<_kwgt.size()&&ix<3;++ix) denom+=_kwgt[ix];
    }
  }
  // additional prefactors
  FPI/=denom;
  // pion mode
  if(imode==0)           FPI *= sqrt(2.0);
  // single kaon modes
  else if(imode==1)      FPI *= sqrt(0.5);
  else if(imode==2)      FPI *= 1.       ;
    // two kaon modes
  else if(imode==3)      FPI *= 1.       ;
  // the kaon eta mode
  else if(imode==4)      FPI *=sqrt(1.5);
  // compute the current
  pdiff-=psum;
  return vector<LorentzPolarizationVectorE>(1,FPI*pdiff);
}
   
bool TwoMesonRhoKStarCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check there are only two particles
  if(id.size()!=2){return false;}
  // pion modes
  if((id[0]==ParticleID::piminus && id[1]==ParticleID::pi0)     ||
     (id[0]==ParticleID::pi0     && id[1]==ParticleID::piminus) ||
     (id[0]==ParticleID::piplus  && id[1]==ParticleID::pi0)     ||
     (id[0]==ParticleID::pi0     && id[1]==ParticleID::piplus)) allowed=true;
  // single charged kaon
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kplus)) allowed=true;
  // single neutral kaon
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::Kbar0)   ||
	  (id[0]==ParticleID::Kbar0   && id[1]==ParticleID::piminus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::K0)      ||
	  (id[0]==ParticleID::K0      && id[1]==ParticleID::piplus)) allowed=true;
  // two kaons
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::K0)     ||
	  (id[0]==ParticleID::K0     && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kbar0)  ||
	  (id[0]==ParticleID::Kbar0  && id[1]==ParticleID::Kplus)) allowed=true;
  // charged kaon and eta
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kplus)) allowed=true;
  return allowed;
}

// the decay mode
unsigned int TwoMesonRhoKStarCurrent::decayMode(vector<int> idout) {
  unsigned int imode(0),nkaon(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::K0) {
      imode=2;
      ++nkaon;
    }
    else if (abs(idout[ix])==ParticleID::Kplus) {
      imode=1;
      ++nkaon;
    }
    else if (idout[ix]==ParticleID::eta) {
      imode=4;
      break;
    }
  }
  if(nkaon==2) imode=3;
  return imode;
}

// output the information for the database
void TwoMesonRhoKStarCurrent::dataBaseOutput(ofstream & output,bool header,
					     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoMesonRhoKStarCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<_rhomasses.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << _rhomasses[ix]/MeV << "\n";
  }
  for(ix=0;ix<_rhowidths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << _rhowidths[ix]/MeV << "\n";
  }
  for(ix=0;ix<_kstarmasses.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarMasses " << ix << " " << _kstarmasses[ix]/MeV << "\n";
  }
  for(ix=0;ix<_kstarwidths.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarWidths " << ix << " " << _kstarwidths[ix]/MeV << "\n";
  }
  output << "newdef " << name() << ":RhoParameters " << _rhoparameters << "\n";
  output << "newdef " << name() << ":KstarParameters " << _kstarparameters << "\n";
  for(ix=0;ix<_piwgt.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PiMagnitude " << ix << " " << _pimag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PiPhase "     << ix << " " << _piphase[ix] << "\n";
  }
  for(ix=0;ix<_kwgt.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KMagnitude " << ix << " " << _kmag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KPhase "     << ix << " " << _kphase[ix] << "\n";
  }
  output << "newdef " << name() << ":PiModel " << _pimodel << "\n";
  output << "newdef " << name() << ":KModel  " << _kmodel  << "\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
