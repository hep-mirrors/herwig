// -*- C++ -*-
//
// TwoPionPhotonCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionPhotonCurrent class.
//
//  Author: Peter Richardson
//

#include "TwoPionPhotonCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

TwoPionPhotonCurrent::TwoPionPhotonCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // weight of the resonances in the current
  _resweights = {1.0,-0.1,0.0};
  // parameters of the rho resonaces
  _rhomasses = {0.773*GeV,1.70*GeV};
  _rhowidths = {0.145*GeV,0.26*GeV};
  // parameters fo the omega resonance
  _omegamass=782*MeV;_omegawidth=8.5*MeV;
  // couplings
  _grho   = 0.11238947*GeV2;
  _grhoomegapi = 12.924/GeV;
  // parameters for the resonance used in the integration
  _intmass  = 1.2*GeV;
  _intwidth = 0.35*GeV;
}

void TwoPionPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(_grho,GeV2) << ounit(_grhoomegapi,1/GeV) << _resweights 
     << ounit(_rhomasses,GeV) << ounit(_rhowidths,GeV) << ounit(_omegamass,GeV) 
     << ounit(_omegawidth,GeV) << ounit(_intmass,GeV) 
     << ounit(_intwidth,GeV) ;
}

void TwoPionPhotonCurrent::persistentInput(PersistentIStream & is, int) { 
  is >> iunit(_grho,GeV2) >> iunit(_grhoomegapi,1/GeV) >> _resweights 
     >> iunit(_rhomasses,GeV) >> iunit(_rhowidths,GeV) >> iunit(_omegamass,GeV) 
     >> iunit(_omegawidth,GeV) >> iunit(_intmass,GeV)
     >> iunit(_intwidth,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionPhotonCurrent,WeakCurrent>
describeHerwigTwoPionPhotonCurrent("Herwig::TwoPionPhotonCurrent", "HwWeakCurrents.so");

void TwoPionPhotonCurrent::Init() {

  static ParVector<TwoPionPhotonCurrent,double> interfacereswgt
    ("Weights",
     "The weights of the different resonances for the decay tau -> nu pi pi gamma",
     &TwoPionPhotonCurrent::_resweights,
     0, 0, 0, -1000, 1000, false, false, true);

  static ParVector<TwoPionPhotonCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the decay tau ->  pi pi photon",
     &TwoPionPhotonCurrent::_rhomasses, MeV, -1, 773.*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoPionPhotonCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the decay tau -> nu pi pi photon",
     &TwoPionPhotonCurrent::_rhowidths, MeV, -1, 145.*MeV, ZERO, 1000.*MeV,
     false, false, true);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceomegamass
    ("omegamass",
     "The mass of the omega",
     &TwoPionPhotonCurrent::_omegamass, GeV, 0.782*GeV, ZERO, 1.0*GeV,
     false, false, true);
  
  static Parameter<TwoPionPhotonCurrent,Energy> interfaceomegawidth
    ("omegawidth",
     "The width of the omega for the decay tau- -> pi pi photon",
     &TwoPionPhotonCurrent::_omegawidth, GeV, 0.0085*GeV, ZERO, 1.*GeV,
     false, false, false);
  
  static ClassDocumentation<TwoPionPhotonCurrent> documentation
    ("The TwoPionPhotonCurrent class implements the decay "
     "tau+/- -> pi+/- pi0 gamma via an omega.",
     "The decay $\\tau^\\pm \\to \\omega \\to \\pi^\\pm \\pi^0 \\gamma$ "
     "is modelled after \\cite{Jadach:1993hs}.",
     "  %\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     );

  static Parameter<TwoPionPhotonCurrent,Energy2> interfacegrho
    ("grho",
     "The rho meson decay constant.",
     &TwoPionPhotonCurrent::_grho, GeV2, 0.11238947*GeV2, -1.*GeV2, 1.*GeV2,
     false, false, false);

  static Parameter<TwoPionPhotonCurrent,InvEnergy> interfacegrhoomegapi
    ("grhoomegapi",
     "The rho-omega-pi coupling",
     &TwoPionPhotonCurrent::_grhoomegapi, 1./GeV, 12.924/GeV,
     -100./GeV, 100./GeV,
     false, false, false);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceIntegrationMass
    ("IntegrationMass",
     "Mass of the pseudoresonance used to improve integration effciency",
     &TwoPionPhotonCurrent::_intmass, GeV, 1.4*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<TwoPionPhotonCurrent,Energy> interfaceIntegrationWidth
    ("IntegrationWidth",
     "Width of the pseudoresonance used to improve integration effciency",
     &TwoPionPhotonCurrent::_intwidth, GeV, 0.5*GeV, ZERO, 10.0*GeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool TwoPionPhotonCurrent::createMode(int icharge, tcPDPtr resonance,
				      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
				      unsigned int imode,PhaseSpaceModePtr mode,
				      unsigned int iloc,int ires,
				      PhaseSpaceChannel phase, Energy upp ) {
  assert(!resonance);
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    switch(i3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
    default:
      return false;
    }
  }
  // check that the mode is are kinematical allowed
  Energy min(getParticleData(ParticleID::piplus)->mass()+
  	     getParticleData(ParticleID::pi0   )->mass());
  if(min>upp) return false;
  // set up the integration channels;
  tPDPtr omega(getParticleData(ParticleID::omega));
  tPDPtr rho;
  if(icharge==-3)     rho = getParticleData(-213);
  else if(icharge==0) rho = getParticleData( 113);
  else if(icharge==3) rho = getParticleData( 213);
  mode->addChannel((PhaseSpaceChannel(phase),ires,rho,
		    ires+1,omega,ires+1,iloc+1,
		    ires+2,iloc+2,ires+2,iloc+3));
  // reset the masses and widths of the resonances if needed
  mode->resetIntermediate(rho,_intmass,_intwidth);
  // set up the omega masses and widths
  mode->resetIntermediate(omega,_omegamass,_omegawidth);
  return true;
}

// the particles produced by the current
tPDVector TwoPionPhotonCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
		       getParticleData(ParticleID::pi0),
		       getParticleData(ParticleID::gamma)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void TwoPionPhotonCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-decay[2]->momentum()
						     ,ix,Helicity::outgoing);
  }
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[2],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
TwoPionPhotonCurrent::current(tcPDPtr resonance,
			      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			      const int imode, const int ,Energy & scale, 
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  assert(!resonance);
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge()+outgoing[2]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    switch(i3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  useMe();
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  }
  // locate the particles
  Lorentz5Momentum pout(momenta[1]+momenta[2]+momenta[0]);
  // overall hadronic mass
  pout.rescaleMass();
  scale=pout.mass();
  Energy2 q2(pout.m2());
  // mass of the omega
  pout = momenta[1]+momenta[2];
  pout.rescaleMass();
  Energy2 s2(pout.m2());
  // compute the prefactor
  complex<InvEnergy3> prefactor(-FFunction(ZERO)*FFunction(q2)*scale*
				sqrt(Constants::twopi*generator()->standardModel()->alphaEM())*
				BreitWigner(s2,10));
  // dot products which don't depend on the polarization vector
  Energy2 dot12(momenta[2]*momenta[1]);
  Energy2 dot13(momenta[2]*momenta[0]);
  Energy2 dot23(momenta[1]*momenta[0]);
  Energy2 mpi2 = sqr(momenta[0].mass());
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix!=1) {
      // obtain the dot products we need
      complex<Energy> dote2 = temp[ix]*momenta[1];
      complex<Energy> dote3 = temp[ix]*momenta[0];
      // now compute the coefficients
      complex<Energy4> coeffa = mpi2*dot13-dot12*(dot23-dot13);
      complex<Energy3> coeffb = dote2*dot13-dote3*dot12;
      complex<Energy3> coeffc = dote2*dot23-dote3*(mpi2+dot12);
      // finally compute the current
      ret[ix]= prefactor*(coeffa*temp[ix]
			  -coeffb*momenta[1]
			  +coeffc*momenta[2]);
    }
    else
      ret[ix]=LorentzPolarizationVectorE();
  }
  return ret;
}

bool TwoPionPhotonCurrent::accept(vector<int> id) {
  if(id.size()!=3){return false;}
  unsigned int npiplus(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::gamma)  ++ngamma;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return (npiplus==1&&ngamma==1&&npi0==1) ||
    (npi0==2&&ngamma==1);
}

unsigned int TwoPionPhotonCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::gamma)   ++ngamma;
  }
  if((npip==1 || npim == 1) && npi0==1 && ngamma==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void TwoPionPhotonCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionPhotonCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":omegamass "    << _omegamass/GeV << "\n";
  output << "newdef " << name() << ":omegawidth "    << _omegawidth/GeV << "\n";
  output << "newdef " << name() << ":grho "    << _grho/GeV2 << "\n";
  output << "newdef " << name() << ":grhoomegapi "    << _grhoomegapi*GeV << "\n";
  output << "newdef " << name() << ":IntegrationMass "  << _intmass/GeV  << "\n";
  output << "newdef " << name() << ":IntegrationWidth " << _intwidth/GeV  << "\n";
  unsigned int ix;
  for(ix=0;ix<_resweights.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":Weights " << ix 
		    << " " << _resweights[ix] << "\n";
    else     output << "insert " << name() << ":Weights " << ix 
		    << " " << _resweights[ix] << "\n";
  }
  for(ix=0;ix<_rhomasses.size();++ix) {
    if(ix<2) output << "newdef " << name() << ":RhoMasses " << ix 
		    << " " << _rhomasses[ix]/MeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix 
		    << " " << _rhomasses[ix]/MeV << "\n";
  }
  for(ix=0;ix<_rhowidths.size();++ix) {
    if(ix<2) output << "newdef " << name() << ":RhoWidths " << ix 
		    << " " << _rhowidths[ix]/MeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix 
		    << " " << _rhowidths[ix]/MeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
