// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhiPiCurrent class.
//

#include "PhiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PhiPiCurrent::PhiPiCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // amplitudes for the weights in the current
  amp_   = {0.045/GeV,0.0315/GeV,0./GeV};
  phase_ = {180.,0.,180.};
  br4pi_     = {0.,0.33,0.};
  // rho masses and widths
  rhoMasses_ = {0.77526*GeV,1.593*GeV,1.909*GeV};
  rhoWidths_ = {0.1491 *GeV,0.203*GeV,0.048*GeV};
}

IBPtr PhiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr PhiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void PhiPiCurrent::doinit() {
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

void PhiPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(wgts_,1./GeV)
     << ounit(mpi_,GeV) << br4pi_;
}

void PhiPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(wgts_,1./GeV)
     >> iunit(mpi_,GeV) >> br4pi_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PhiPiCurrent,WeakCurrent>
describeHerwigPhiPiCurrent("Herwig::PhiPiCurrent",
			   "HwWeakCurrents.so");

void PhiPiCurrent::Init() {

  static ClassDocumentation<PhiPiCurrent> documentation
    ("The PhiPiCurrent class implements a model of the current for phi pi"
     "based on the model of Phys.Rev. D77 (2008) 092002, 2008.",
     "The current for $\\phi\\pi$ based on \\cite{Aubert:2007ym} was used.",
     "\\bibitem{Aubert:2007ym}\n"
     "B.~Aubert {\\it et al.} [BaBar Collaboration],\n"
     "%``Measurements of $e^{+} e^{-} \\to K^{+} K^{-} \\eta$,"
     " $K^{+} K^{-} \\pi^0$ and $K^0_{s} K^\\pm \\pi^\\mp$ "
     "cross sections using initial state radiation events,''\n"
     "Phys.\\ Rev.\\ D {\\bf 77} (2008) 092002\n"
     "doi:10.1103/PhysRevD.77.092002\n"
     "[arXiv:0710.4451 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.77.092002;%%\n"
     "%153 citations counted in INSPIRE as of 27 Aug 2018\n");

  static ParVector<PhiPiCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons",
     &PhiPiCurrent::rhoMasses_, GeV, 3, 775.26*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PhiPiCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons",
     &PhiPiCurrent::rhoWidths_, GeV, 3, 149.1*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PhiPiCurrent,InvEnergy> interfaceAmplitudes
    ("Amplitudes",
     "The amplitudes for the different resonances",
     &PhiPiCurrent::amp_, 1./GeV, 3, 0./GeV, 0./GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<PhiPiCurrent,double> interfacePhase
    ("Phase",
     "The phases for the different rho resonances in degrees",
     &PhiPiCurrent::phase_, 3, 0., 0.0, 360.,
     false, false, Interface::limited);
  
  static ParVector<PhiPiCurrent,double> interfaceBR4Pi
    ("BR4Pi",
     "The branching ratios to 4 pi for the various resonances",
     &PhiPiCurrent::br4pi_, 3, 0., 0.0, 1.0,
     false, false, Interface::limited);


}

// complete the construction of the decay mode for integration
bool PhiPiCurrent::createMode(int icharge, tcPDPtr resonance,
			      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			      unsigned int imode,PhaseSpaceModePtr mode,
			      unsigned int iloc,int ires,
			      PhaseSpaceChannel phase, Energy upp ) {
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
  Energy min = getParticleData(ParticleID::phi)->massMin();
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
tPDVector PhiPiCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
		       getParticleData(ParticleID::phi)};
  if(imode==0) {
    if(icharge==3)       extpart[0] = getParticleData(ParticleID::piplus );
    else if(icharge==-3) extpart[0] = getParticleData(ParticleID::piminus);
  }
  else {
    extpart[0] = getParticleData(ParticleID::pi0);
  }
  return extpart;
}

void PhiPiCurrent::constructSpinInfo(ParticleVector decay) const {
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
PhiPiCurrent::current(tcPDPtr resonance,
		      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
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
  for(unsigned int ix=0;ix<3;++ix)
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  // locate the particles
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  // work out the channel
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
  complex<InvEnergy> pre(ZERO);
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy2 mR2 = sqr(rhoMasses_[ix]);
    Energy wid = rhoWidths_[ix]*
      (1.-br4pi_[ix]+ br4pi_[ix]*mR2/q2*pow((q2-16.*sqr(mpi_))/(mR2-16.*sqr(mpi_)),1.5));
    pre += wgts_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*wid);
  }
  vector<LorentzPolarizationVectorE> ret(3);
  if(imode==0) pre *= sqrt(2.);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] = pre*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool PhiPiCurrent::accept(vector<int> id) {
  if(id.size()!=2){return false;}
  unsigned int npiplus(0),npi0(0),nphi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::phi)  ++nphi;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return nphi==1 && (npiplus==1||npi0==1);
}

unsigned int PhiPiCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),nphi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::phi)   ++nphi;
  }
  if((npip==1 || npim == 1) && nphi==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void PhiPiCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::PhiPiCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
    output << "newdef " << name() << ":RhoMasses " << ix 
	   << " " << rhoMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
    output << "newdef " << name() << ":RhoWidths " << ix 
	   << " " << rhoWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    output << "newdef " << name() << ":Amplitudes " << ix 
	   << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    output << "newdef " << name() << ":Phases " << ix 
	   << " " << phase_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    output << "newdef " << name() << ":BR4Pi " << ix 
	   << " " << br4pi_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
