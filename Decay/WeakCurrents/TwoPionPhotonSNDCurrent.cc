// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionPhotonSNDCurrent class.
//

#include "TwoPionPhotonSNDCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

TwoPionPhotonSNDCurrent::TwoPionPhotonSNDCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // amplitudes for the weights in the current
  amp_   = {1.,0.175,0.014};
  phase_ = {0.,124.,-63.};
  // rho masses and widths
  rhoMasses_ = {0.77526*GeV,1.510*GeV,1.720*GeV};
  rhoWidths_ = {0.1491 *GeV,0.44 *GeV,0.25 *GeV};
  // coupling
  gRhoOmegaPi_   = 15.9/GeV;
  fRho_        = 4.9583;
  gGammaOmegaPi_ = 0.695821538653/GeV;
  fRho_        = 4.9583;
  // omega parameters
  omegaMass_  = 782.65*MeV;
  omegaWidth_ = 8.49 *MeV;
}

IBPtr TwoPionPhotonSNDCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoPionPhotonSNDCurrent::fullclone() const {
  return new_ptr(*this);
}

void TwoPionPhotonSNDCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  wgts_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    wgts_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

void TwoPionPhotonSNDCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << amp_ << phase_ << wgts_ << fRho_
     << ounit(gRhoOmegaPi_,1./GeV) << ounit(gGammaOmegaPi_,1./GeV)
     << ounit(omegaMass_,GeV) << ounit(omegaWidth_,GeV) << ounit(mpi_,GeV);
}

void TwoPionPhotonSNDCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> amp_ >> phase_ >> wgts_ >> fRho_
     >> iunit(gRhoOmegaPi_,1./GeV) >> iunit(gGammaOmegaPi_,1./GeV)
     >> iunit(omegaMass_,GeV) >> iunit(omegaWidth_,GeV) >> iunit(mpi_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionPhotonSNDCurrent,WeakCurrent>
describeHerwigTwoPionPhotonSNDCurrent("Herwig::TwoPionPhotonSNDCurrent",
				      "HwWeakCurrents.so");

void TwoPionPhotonSNDCurrent::Init() {

  static ClassDocumentation<TwoPionPhotonSNDCurrent> documentation
    ("The TwoPionPhotonSNDCurrent class provides the weka current for"
     "pi pi gamma using the model of SND",
     "The current based on \\cite{Achasov:2016zvn} for $\\pi\\pi^0\\gamma$ was used.\n",
     "\\bibitem{Achasov:2016zvn}"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Updated measurement of the $e^+e^- \\to \\omega \\pi^0 \\to \\pi^0\\pi^0\\gamma$ cross section with the SND detector,''\n"
     "Phys.\\ Rev.\\ D {\\bf 94} (2016) no.11,  112001\n"
     "doi:10.1103/PhysRevD.94.112001\n"
     "[arXiv:1610.00235 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.94.112001;%%\n"
     "%12 citations counted in INSPIRE as of 22 Aug 2018\n");

  static ParVector<TwoPionPhotonSNDCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons",
     &TwoPionPhotonSNDCurrent::rhoMasses_, GeV, -1, 775.26*MeV,
     0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<TwoPionPhotonSNDCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons",
     &TwoPionPhotonSNDCurrent::rhoWidths_, GeV, -1, 0.1491*GeV,
     0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<TwoPionPhotonSNDCurrent,double> interfaceAmplitudes
    ("Amplitudes",
     "THe amplitudes for the different rho resonances",
     &TwoPionPhotonSNDCurrent::amp_, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static ParVector<TwoPionPhotonSNDCurrent,double> interfacePhase
    ("Phase",
     "The phases for the different rho resonances in degrees",
     &TwoPionPhotonSNDCurrent::phase_, -1, 0., 0.0, 360.,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,double> interfacefRho
    ("fRho",
     "The coupling of the photon and the rho meson",
     &TwoPionPhotonSNDCurrent::fRho_, 4.9583, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,InvEnergy> interfacegRhoOmegaPi
    ("gRhoOmegaPi",
     "The coupling rho-omega-pi",
     &TwoPionPhotonSNDCurrent::gRhoOmegaPi_, 1./GeV,
     15.9/GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,InvEnergy> interfacegGammaOmegaPi
    ("gGammaOmegaPi",
     "The coupling gamma-omega-pi",
     &TwoPionPhotonSNDCurrent::gGammaOmegaPi_, 1./GeV,
     0.695821538653/GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &TwoPionPhotonSNDCurrent::omegaMass_, GeV, 0.78265*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoPionPhotonSNDCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &TwoPionPhotonSNDCurrent::omegaWidth_, GeV, 8.49*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool TwoPionPhotonSNDCurrent::createMode(int icharge, tcPDPtr resonance,
					 FlavourInfo flavour,
					 unsigned int imode,PhaseSpaceModePtr mode,
					 unsigned int iloc,int ires,
					 PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return false;
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return false;
      break;
    default:
      return false;
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return false;
  // check that the mode is are kinematical allowed
  Energy min(getParticleData(ParticleID::piplus)->mass()+
  	     getParticleData(ParticleID::pi0   )->mass());
  if(min>upp) return false;
  // set up the integration channels;
  tPDPtr omega(getParticleData(ParticleID::omega));
  vector<tPDPtr> rho;
  if(icharge==-3)
    rho = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  else if(icharge==0)
    rho = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  else if(icharge==3)
    rho = {getParticleData( 213),getParticleData( 100213),getParticleData( 30213)};
  for(unsigned int ix=0;ix<3;++ix) {
    if(resonance && resonance!=rho[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,rho[ix],
		      ires+1,omega,ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
    // channel with the pions exchanged
    if(icharge==0) 
      mode->addChannel((PhaseSpaceChannel(phase),ires,rho[ix],
			ires+1,omega,ires+1,iloc+2,
			ires+2,iloc+1,ires+2,iloc+3));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<3;++ix) {
    mode->resetIntermediate(rho[ix],rhoMasses_[ix],rhoWidths_[ix]);
  }
  // set up the omega masses and widths
  mode->resetIntermediate(omega,omegaMass_,omegaWidth_);
  return true;
}

// the particles produced by the current
tPDVector TwoPionPhotonSNDCurrent::particles(int icharge, unsigned int imode,int,int) {
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

void TwoPionPhotonSNDCurrent::constructSpinInfo(ParticleVector decay) const {
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
TwoPionPhotonSNDCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge()+outgoing[2]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode!=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>1 || icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>1 || icharge ==3) return vector<LorentzPolarizationVectorE>();
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
  useMe();
  // polarization vectors of the photon
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]+momenta[2]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  unsigned int imin=0, imax = wgts_.size();
  if(ichan>0) {
    if(outgoing[0]!=outgoing[1])
      imin = ichan;
    else
      imin = ichan/2;
    imax = imin+1;
  }
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imin = 0;
      break;
    case 100:
      imin = 1;
      break;
    case 30 :
      imin = 2;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  vector<LorentzPolarizationVectorE> ret(3);
  // need to include exchange of identical particles for the I_3=0 case
  for(int iorder=0;iorder<2;++iorder) {
    Lorentz5Momentum pout(momenta[2]);
    if(outgoing[0]==outgoing[1]) {
      if(ichan>=0&& ichan%2!=iorder) continue;
    }
    else if(iorder==1) continue;
    // add pion momentum
    if(iorder==0) pout += momenta[1];
    else          pout += momenta[0];
    // mass of the omega
    pout.rescaleMass();
    Energy2 s2(pout.m2());
    // compute the rho width
    Energy2 mr2(sqr(rhoMasses_[0]));
    Energy grho = rhoWidths_[0]*mr2/q2*pow(max(double((q2-4.*sqr(mpi_))/(mr2-4.*sqr(mpi_))),0.),1.5);
    Energy qw = Kinematics::pstarTwoBodyDecay(q.mass(),pout.mass(),mpi_);
    grho += pow<3,1>(qw)*sqr(gRhoOmegaPi_)/12./Constants::pi;
    // compute the prefactor
    complex<InvEnergy4> pre = gRhoOmegaPi_*gGammaOmegaPi_/fRho_*
      Resonance::BreitWignerFW(s2,omegaMass_,omegaWidth_)/sqr(omegaMass_);
    if(imode==0) pre *=sqrt(2.);
    Complex bw(0.);
    for(unsigned int ix=imin;ix<imax;++ix) {
      Energy wid = ix==0 ? grho : rhoWidths_[ix];
      Energy2 mR2 = sqr(rhoMasses_[ix]);
      bw += mR2*wgts_[ix]/(mR2-q2-Complex(0.,1.)*q.mass()*wid);
    }
    pre *=bw;
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix==1) continue;
      LorentzVector<complex<Energy2> > v2 = Helicity::epsilon(pout,temp[ix],momenta[2]);
      ret[ix] += pre*scale*Helicity::epsilon(q,v2,pout);
    }
  }
  return ret;
}

bool TwoPionPhotonSNDCurrent::accept(vector<int> id) {
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

unsigned int TwoPionPhotonSNDCurrent::decayMode(vector<int> id) {
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
void TwoPionPhotonSNDCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionPhotonSNDCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":RhoMasses " << ix 
		    << " " << rhoMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix 
		    << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":RhoWidths " << ix 
		    << " " << rhoWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix 
		    << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":Amplitudes " << ix 
		    << " " << amp_[ix] << "\n";
    else     output << "insert " << name() << ":Amplitudes " << ix 
		    << " " << amp_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<3) output << "newdef " << name() << ":Phases " << ix 
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phases " << ix 
		    << " " << phase_[ix] << "\n";
  }
  output << "newdef " << name() << ":fRho "    << fRho_ << "\n";
  output << "newdef " << name() << ":gRhoOmegaPi "    << gRhoOmegaPi_*GeV << "\n";
  output << "newdef " << name() << ":gGammaOmegaPi "    << gGammaOmegaPi_*GeV << "\n";
  output << "newdef " << name() << ":OmegaMass "    << omegaMass_/GeV << "\n";
  output << "newdef " << name() << ":OmegaWidth "    << omegaWidth_/GeV << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
