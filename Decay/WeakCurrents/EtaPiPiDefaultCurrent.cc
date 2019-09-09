// -*- C++ -*-
//
// EtaPiPiDefaultCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiDefaultCurrent class.
//

#include "EtaPiPiDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"
using namespace Herwig;
using namespace ThePEG;

DescribeClass<EtaPiPiDefaultCurrent,WeakCurrent>
describeHerwigEtaPiPiDefaultCurrent("Herwig::EtaPiPiDefaultCurrent",
				    "HwWeakCurrents.so");

EtaPiPiDefaultCurrent::EtaPiPiDefaultCurrent() {
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi=ZERO;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts = { 1.0,-0.145,0.};
  _rhoF5wgts   = {-26.,   6.5,1.};
  // local values of the rho parameters
  _rhoF123masses = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rhoF123widths = {0.145*GeV,0.510*GeV,0.120*GeV};
  _rhoF5masses   = {0.773*GeV,1.500*GeV,1.750*GeV};
  _rhoF5widths   = {0.145*GeV,0.220*GeV,0.120*GeV};
}

void EtaPiPiDefaultCurrent::doinit() {
  WeakCurrent::doinit();
  // the particles we will use a lot
  tPDPtr piplus(getParticleData(ParticleID::piplus));
  // masses for the running widths
  _mpi=piplus->mass();
}

void EtaPiPiDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts << _rhoF5wgts << ounit(_fpi,GeV) << ounit(_mpi,GeV)
     << ounit(_rhoF123masses,GeV) << ounit(_rhoF5masses,GeV) 
     << ounit(_rhoF123widths,GeV) << ounit(_rhoF5widths,GeV);
}

void EtaPiPiDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >> _rhoF5wgts >> iunit(_fpi,GeV) >> iunit(_mpi,GeV)
     >> iunit(_rhoF123masses,GeV) >> iunit(_rhoF5masses,GeV) 
     >> iunit(_rhoF123widths,GeV) >> iunit(_rhoF5widths,GeV);
}

void EtaPiPiDefaultCurrent::Init() {
        
  static ClassDocumentation<EtaPiPiDefaultCurrent> documentation
    ("The EtaPiPiDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.",
     "The three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, "
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "and pi- Kbar0 pi0, pi- pi0 eta "
     "use the same currents as \\cite{Jadach:1993hs,Kuhn:1990ad,Decker:1992kj}.",
     "%\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Decker:1992kj}\n"
     "\\bibitem{Decker:1992kj}\n"
     "  R.~Decker, E.~Mirkes, R.~Sauer and Z.~Was,\n"
     "  %``Tau decays into three pseudoscalar mesons,''\n"
     "  Z.\\ Phys.\\  C {\\bf 58}, 445 (1993).\n"
     "  %%CITATION = ZEPYA,C58,445;%%\n"
     );
  
  static ParVector<EtaPiPiDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &EtaPiPiDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,double> interfaceF5RhoWgt
    ("F5RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &EtaPiPiDefaultCurrent::_rhoF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF5masses
    ("rhoF5masses",
     "The masses for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF5masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<EtaPiPiDefaultCurrent,Energy> interfacerhoF5widths
    ("rhoF5widths",
     "The widths for the rho resonances if used local values",
     &EtaPiPiDefaultCurrent::_rhoF5widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<EtaPiPiDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &EtaPiPiDefaultCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool EtaPiPiDefaultCurrent::createMode(int icharge, tcPDPtr resonance,
				       FlavourInfo flavour,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return false;
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero ) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       ) return false;
  // make sure that the decays are kinematically allowed
  int iq(0),ia(0);
  tPDVector part = particles(icharge,imode,iq,ia);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // set up the resonances
  tPDPtr res[3];
  if(icharge==0) {
    res[0] =getParticleData(113);
    res[1] =getParticleData(100113);
    res[2] =getParticleData(30113);
  }
  else {
    res[0] =getParticleData(213);
    res[1] =getParticleData(100213);
    res[2] =getParticleData(30213);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<3;++ix) {
  	if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
      }
    }
  }
  // channels for pi- pi0 eta
  for(unsigned int ix=0;ix<3;++ix) {
    if(resonance && resonance != res[ix]) continue;
    for(unsigned int iy=0;iy<3;++iy) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,iloc+3,ires+1,res[iy],
			ires+2,iloc+1,ires+2,iloc+2));
    }
  }
  // reset the rho masses
  for(unsigned int ix=0;ix<_rhoF5masses.size();++ix)
    mode->resetIntermediate(res[ix],_rhoF5masses[ix],_rhoF5widths[ix]);
  return true;
}

void EtaPiPiDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPiPiDefaultCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":FPi "     << _fpi/MeV     << "\n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123RhoWeight " << ix << " " << _rhoF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F5RhoWeight " << ix << " " << _rhoF5wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123masses " << ix 
	   << " " << _rhoF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123widths " << ix << " " 
	   << _rhoF123widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF5masses " << ix << " " 
	   << _rhoF5masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF5widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF5widths " << ix << " " 
	   << _rhoF5widths[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
EtaPiPiDefaultCurrent::current(tcPDPtr resonance,
			       FlavourInfo flavour,
			       const int imode, const int ichan, Energy & scale, 
			       const tPDVector & outgoing,
			       const vector<Lorentz5Momentum> & momenta,
			       DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode==0) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode==1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode==1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return vector<LorentzPolarizationVectorE>();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q = momenta[0] + momenta[1] + momenta[2];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s3 = (momenta[0]+momenta[1]).m2();
  // the form factor
  Complex F5(0.);
  int ires1(-1),ires2(-1);
  if(ichan>=0) {
    ires1 = ichan/3;
    ires2 = ichan%3;
  }
  else {
    if(resonance) {
      switch(resonance->id()/1000) {
      case 0:
	ires1 = 0;
	break;
      case 100:
	ires1 = 1;
	break;
      case 30 :
	ires1 = 2;
	break;
      default:
	assert(false);
      }
    }
  }
  F5 = BrhoF5(q2,ires1)*BrhoF123(s3,ires2)*sqrt(2./3.);
  // constant the current
  LorentzPolarizationVector vect = -Complex(0.,1.)*F5/sqr(Constants::twopi)/pow<3,1>(_fpi)*
    Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()*vect);
}

bool EtaPiPiDefaultCurrent::accept(vector<int> id) {
  // check there are only three particles
  if(id.size()!=3) return false;
  unsigned int npip(0),npim(0),npi0(0),neta(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::eta)     ++neta;
  }
  if( (npip==1&&npim==1&&neta==1) ||
      (npi0==1&&npim+npip==1&&neta==1))
    return true;
  else
    return false;
}

unsigned int EtaPiPiDefaultCurrent::decayMode(vector<int> id) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==2) return 1;
  else       return 0;
}

tPDVector EtaPiPiDefaultCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector output(3);
  output[0]=getParticleData(ParticleID::piplus);
  output[2]=getParticleData(ParticleID::eta);
  if(imode==0) {
    output[1]=getParticleData(ParticleID::pi0);
    if(icharge==-3) {
      for(unsigned int ix=0;ix<output.size();++ix) {
	if(output[ix]->CC()) output[ix]=output[ix]->CC();
      }
    }
  }
  else {
    output[1]=getParticleData(ParticleID::piminus);
  }
  return output;
}
