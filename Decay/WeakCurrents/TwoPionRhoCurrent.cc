// -*- C++ -*-
//
// TwoPionRhoCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionRhoCurrent class.
//
//  Author: Peter Richardson
//

#include "TwoPionRhoCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
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

TwoPionRhoCurrent::TwoPionRhoCurrent() {
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(4);
  // the weights of the different resonances in the matrix elements
  _pimag   = {1.0,0.167,0.05};
  _piphase = {0.0,180  ,0.0};
  // models to use
  _pimodel = 0;
  // parameter for the masses (use the parameters freom the CLEO fit 
  // rather than the PDG masses etc)
  _rhoparameters=true;
  _rhomasses = {774.6*MeV,1408*MeV,1700*MeV};
  _rhowidths = {149*MeV,502*MeV,235*MeV};
}

void TwoPionRhoCurrent::doinit() {
  WeakDecayCurrent::doinit();
  // check consistency of parametrers
  if(_rhomasses.size()!=_rhowidths.size()) {
    throw InitException() << "Inconsistent parameters in TwoPionRhoCurrent"
			  << "::doinit()" << Exception::abortnow;
  }
  // the resonances
  tPDPtr res[3]={getParticleData(-213   ),
		 getParticleData(-100213),
		 getParticleData(-30213 )};
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
  // set up for the Breit Wigners
  Energy mpi0(   getParticleData(ParticleID::pi0   )->mass());
  Energy mpiplus(getParticleData(ParticleID::piplus)->mass());
  // rho resonances
  for(unsigned int ix=0;ix<3;++ix) {
    _mass.push_back(_rhomasses[ix]);
    _width.push_back(_rhowidths[ix]);
    _massa.push_back(mpi0);
    _massb.push_back(mpiplus);
    _hres.push_back(Resonance::Hhat(sqr(_mass.back()),_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _dh.push_back(Resonance::dHhatds(_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _h0.push_back(Resonance::H(ZERO,_mass.back(),_width.back(),_massa.back(),_massb.back(),_dh.back(),_hres.back()));
  }
  // weights for the rho channels
  if(_pimag.size()!=_piphase.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
			  << " rho channel must be the same size in "
			  << "TwoPionRhoCurrent::doinit()" << Exception::runerror;
  _piwgt.resize(_pimag.size());
  for(unsigned int ix=0;ix<_pimag.size();++ix) {
    double angle = _piphase[ix]/180.*Constants::pi;
    _piwgt[ix] = _pimag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
}

void TwoPionRhoCurrent::persistentOutput(PersistentOStream & os) const {
  os << _pimodel << _piwgt << _pimag << _piphase
     << _rhoparameters << ounit(_rhomasses,GeV) << ounit(_rhowidths,GeV)
     << ounit(_mass,GeV) << ounit(_width,GeV)
     << ounit(_massa,GeV) <<ounit(_massb,GeV)
     << _dh << ounit(_hres,GeV2) << ounit(_h0,GeV2);
}

void TwoPionRhoCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _pimodel >> _piwgt >> _pimag >> _piphase
     >> _rhoparameters >> iunit(_rhomasses,GeV) >> iunit(_rhowidths,GeV) 
     >> iunit(_mass,GeV) >> iunit(_width,GeV)
     >> iunit(_massa,GeV) >> iunit(_massb,GeV)
     >> _dh >> iunit(_hres,GeV2) >> iunit(_h0,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionRhoCurrent,WeakDecayCurrent>
describeHerwigTwoPionRhoCurrent("Herwig::TwoPionRhoCurrent", "HwWeakCurrents.so");

void TwoPionRhoCurrent::Init() {

  static ParVector<TwoPionRhoCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoPionRhoCurrent::_rhomasses, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoPionRhoCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoPionRhoCurrent::_rhowidths, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static Switch<TwoPionRhoCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values for the rho meson masses and widths",
     &TwoPionRhoCurrent::_rhoparameters, true, false, false);
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
  
  static ParVector<TwoPionRhoCurrent,double> interfacePiMagnitude
    ("PiMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoPionRhoCurrent::_pimag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoPionRhoCurrent,double> interfacePiPhase
    ("PiPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoPionRhoCurrent::_piphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static Switch<TwoPionRhoCurrent,int> interfacePiModel
    ("PiModel",
     "The model to use for the propagator for the pion modes.",
     &TwoPionRhoCurrent::_pimodel, 0, false, false);
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

  static ClassDocumentation<TwoPionRhoCurrent> documentation
    ("The TwoPionRhoCurrent class is designed to implement weak"
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
bool TwoPionRhoCurrent::createMode(int icharge, unsigned int imode,
					 DecayPhaseSpaceModePtr mode,
					 unsigned int iloc,unsigned int,
					 DecayPhaseSpaceChannelPtr phase,Energy upp) {
  if((imode<=1&&abs(icharge)!=3) ||
     (imode>1 && icharge !=0)) return false; 
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::K0);
  }
  else if(imode==2 || imode==3 ) {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::piminus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  DecayPhaseSpaceChannelPtr newchannel;
  // set the resonances
  // two pion or  K+ K0 decay
  tPDPtr res[3]={getParticleData(213),getParticleData(100213),getParticleData(30213)};
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
  if(_rhoparameters) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_rhomasses.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_rhomasses[ix],_rhowidths[ix]);
      }
    }
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector TwoPionRhoCurrent::particles(int icharge, unsigned int imode,
					     int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==2) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::Kbar0);
  }
  else {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::piminus);
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
TwoPionRhoCurrent::current(const int imode, const int ichan,
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
  if(ichan<0) {
    for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix) {
      FPI+=_piwgt[ix]*BreitWigner(q2,_pimodel,ix);
      denom+=_piwgt[ix];
    }
  }
  else if(ichan<int(_piwgt.size())&&ichan<3) {
    FPI=_piwgt[ichan]*BreitWigner(q2,_pimodel,ichan);
    for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix) denom+=_piwgt[ix];
  }
  // additional prefactors
  FPI/=denom;
  // pion mode
  if(imode==0)           FPI *= sqrt(2.0);
    // two kaon modes
  else if(imode==1)      FPI *= 1.       ;
  // compute the current
  pdiff-=psum;
  return vector<LorentzPolarizationVectorE>(1,FPI*pdiff);
}
   
bool TwoPionRhoCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // pion modes
  if((abs(id[0])==ParticleID::piplus  &&     id[1] ==ParticleID::pi0   ) ||
     (    id[0] ==ParticleID::pi0     && abs(id[1])==ParticleID::piplus))
    return true;
  // single charged kaon
  else if((abs(id[0])==ParticleID::Kplus  &&     id[1] ==ParticleID::pi0  ) ||
	  (    id[0] ==ParticleID::pi0    && abs(id[1])==ParticleID::Kplus))
    return true;
  // single neutral kaon
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::Kbar0)   ||
	  (id[0]==ParticleID::Kbar0   && id[1]==ParticleID::piminus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::K0)      ||
	  (id[0]==ParticleID::K0      && id[1]==ParticleID::piplus))
    return true;
  // two kaons
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::K0)     ||
	  (id[0]==ParticleID::K0     && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::Kbar0)  ||
	  (id[0]==ParticleID::Kbar0  && id[1]==ParticleID::Kplus))
    return true;
  // charged kaon and eta
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kplus))
    return true;
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::piplus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::piminus))
    return true;
  else
    return false;
}

// the decay mode
unsigned int TwoPionRhoCurrent::decayMode(vector<int> idout) {
  unsigned int imode(0),nkaon(0),npi(0);
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
    else if(abs(idout[ix])==ParticleID::piplus) {
      ++npi;
    }
  }
  if(nkaon==2)    return 3;
  else if(npi==2) return 5;
  else            return imode;
}

// output the information for the database
void TwoPionRhoCurrent::dataBaseOutput(ofstream & output,bool header,
					     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionRhoCurrent " 
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
  output << "newdef " << name() << ":RhoParameters " << _rhoparameters << "\n";
  for(ix=0;ix<_piwgt.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PiMagnitude " << ix << " " << _pimag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":PiPhase "     << ix << " " << _piphase[ix] << "\n";
  }
  output << "newdef " << name() << ":PiModel " << _pimodel << "\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
