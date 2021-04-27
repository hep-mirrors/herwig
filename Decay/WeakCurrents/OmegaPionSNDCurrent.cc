// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPionSNDCurrent class.
//

#include "OmegaPionSNDCurrent.h"
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

OmegaPionSNDCurrent::OmegaPionSNDCurrent() {
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
  //fRho_        = 4.9583;
  fRho_ = 5.06325;//evaluated with alphaEM in considered energy range
}

IBPtr OmegaPionSNDCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr OmegaPionSNDCurrent::fullclone() const {
  return new_ptr(*this);
}

void OmegaPionSNDCurrent::doinit() {
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

void OmegaPionSNDCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << amp_ << phase_ << wgts_ << fRho_
     << ounit(gRhoOmegaPi_,1./GeV) << ounit(mpi_,GeV);
}

void OmegaPionSNDCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> amp_ >> phase_ >> wgts_ >> fRho_
     >> iunit(gRhoOmegaPi_,1./GeV) >> iunit(mpi_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<OmegaPionSNDCurrent,WeakCurrent>
describeHerwigOmegaPionSNDCurrent("Herwig::OmegaPionSNDCurrent",
				  "HwWeakCurrents.so");

void OmegaPionSNDCurrent::Init() {

  static ClassDocumentation<OmegaPionSNDCurrent> documentation
    ("The OmegaPionSNDCurrent class provides a current for omega pi"
     " using the model of SND",
     "The current based on \\cite{Achasov:2016zvn} for $\\omega\\pi$ was used.\n",
     "\\bibitem{Achasov:2016zvn}"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Updated measurement of the $e^+e^- \\to \\omega \\pi^0 \\to \\pi^0\\pi^0\\gamma$ cross section with the SND detector,''\n"
     "Phys.\\ Rev.\\ D {\\bf 94} (2016) no.11,  112001\n"
     "doi:10.1103/PhysRevD.94.112001\n"
     "[arXiv:1610.00235 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.94.112001;%%\n"
     "%12 citations counted in INSPIRE as of 22 Aug 2018\n");

  static ParVector<OmegaPionSNDCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons",
     &OmegaPionSNDCurrent::rhoMasses_, GeV, -1, 775.26*MeV,
     0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OmegaPionSNDCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons",
     &OmegaPionSNDCurrent::rhoWidths_, GeV, -1, 0.1491*GeV,
     0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OmegaPionSNDCurrent,double> interfaceAmplitudes
    ("Amplitudes",
     "THe amplitudes for the different rho resonances",
     &OmegaPionSNDCurrent::amp_, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static ParVector<OmegaPionSNDCurrent,double> interfacePhase
    ("Phase",
     "The phases for the different rho resonances in degrees",
     &OmegaPionSNDCurrent::phase_, -1, 0., -360., 360.,
     false, false, Interface::limited);

  static Parameter<OmegaPionSNDCurrent,double> interfacefRho
    ("fRho",
     "The coupling of the photon and the rho meson",
     &OmegaPionSNDCurrent::fRho_, 4.9583, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<OmegaPionSNDCurrent,InvEnergy> interfacegRhoOmegaPi
    ("gRhoOmegaPi",
     "The coupling rho-omega-pi",
     &OmegaPionSNDCurrent::gRhoOmegaPi_, 1./GeV,
     15.9/GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

}

// complete the construction of the decay mode for integration
bool OmegaPionSNDCurrent::createMode(int icharge, tcPDPtr resonance,
					 FlavourInfo flavour,
					 unsigned int imode,PhaseSpaceModePtr mode,
					 unsigned int iloc,int ires,
					 PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode >= 1))
    return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown || flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
    // and check I_3
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
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return false;
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::omega)->massMin();
  if(imode==0)
    min += getParticleData(ParticleID::piplus)->mass();
  else
    min += getParticleData(ParticleID::pi0   )->mass();
  if(min>upp) return false;
  // set up the integration channels;
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
		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<3;++ix) {
    mode->resetIntermediate(rho[ix],rhoMasses_[ix],rhoWidths_[ix]);
  }
  return true;
}

// the particles produced by the current
tPDVector OmegaPionSNDCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
		       getParticleData(ParticleID::omega)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void OmegaPionSNDCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum()
						     ,ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents
vector<LorentzPolarizationVectorE>
OmegaPionSNDCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale,
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check the charge
  if((abs(icharge)!=3 && imode == 0) ||
     (   icharge!=0   && imode == 1))
    return vector<LorentzPolarizationVectorE>();
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown || flavour.I3!=IsoSpin::I3Unknown) {
    if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
    // and check I_3
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
  // other flavour stuff
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero     ) return vector<LorentzPolarizationVectorE>();
  useMe();
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix)
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  // locate the particles
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  // compute the rho width
  Energy2 mr2(sqr(rhoMasses_[0]));
  Energy grho = rhoWidths_[0]*mr2/q2*pow((q2-4.*sqr(mpi_))/(mr2-4.*sqr(mpi_)),1.5);
  Energy qw = Kinematics::pstarTwoBodyDecay(q.mass(),momenta[0].mass(),momenta[1].mass());
  grho += pow<3,1>(qw)*sqr(gRhoOmegaPi_)/12./Constants::pi;
  unsigned int imin=0, imax = wgts_.size();
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
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
  // compute the prefactor
  complex<InvEnergy> pre = gRhoOmegaPi_/fRho_;
  Complex bw(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy wid = ix==0 ? grho : rhoWidths_[ix];
    Energy2 mR2 = sqr(rhoMasses_[ix]);
    bw += mR2*wgts_[ix]/(mR2-q2-Complex(0.,1.)*q.mass()*wid);
  }
  pre = pre * bw;
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] = pre*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  if(imode==0) pre *=sqrt(2.);
  return ret;
}

bool OmegaPionSNDCurrent::accept(vector<int> id) {
  if(id.size()!=2){return false;}
  unsigned int npiplus(0),npi0(0),nomega(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::omega)  ++nomega;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return nomega==1 && (npiplus==1||npi0==1);
}

unsigned int OmegaPionSNDCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),nomega(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::omega)   ++nomega;
  }
  if((npip==1 || npim == 1) && nomega==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void OmegaPionSNDCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::OmegaPionSNDCurrent " << name()
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
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}
