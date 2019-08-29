// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PionPhotonCurrent class.
//

#include "PionPhotonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

PionPhotonCurrent::PionPhotonCurrent() {
  // modes handled
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {0.77526*GeV,0.78284*GeV,1.01952*GeV,1.45*GeV,1.70*GeV};
  // widths for the resonances
  resWidths_ = {0.1491 *GeV,0.00868*GeV,0.00421*GeV,0.40*GeV,0.30*GeV};
  // amplitudes
  amp_   = {0.0426/GeV,0.0434/GeV,0.00303/GeV,0.00523/GeV,ZERO};
  // phases
  phase_ = {-12.7,0.,158.,180.,0.};
}

IBPtr PionPhotonCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr PionPhotonCurrent::fullclone() const {
  return new_ptr(*this);
}

void PionPhotonCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  couplings_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    couplings_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

void PionPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV)
     << ounit(mpi_,GeV);
}

void PionPhotonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV)
     >> iunit(mpi_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PionPhotonCurrent,WeakCurrent>
describeHerwigPionPhotonCurrent("Herwig::PionPhotonCurrent",
				"HwWeakCurrents.so");

void PionPhotonCurrent::Init() {

  static ClassDocumentation<PionPhotonCurrent> documentation
    ("The PionPhotonCurrent class implements a current based"
     " on the model of SND for pion+photon",
     "The current based on the model of \\cite{Achasov:2016bfr}"
     " for pion and photon was used.",
     "\\bibitem{Achasov:2016bfr}\n"
     "M.~N.~Achasov {\\it et al.} [SND Collaboration],\n"
     "%``Study of the reaction $e^+e^- \\to \\pi^0\\gamma$ with the SND detector at the VEPP-2M collider,''\n"
     "Phys.\\ Rev.\\ D {\\bf 93} (2016) no.9,  092001\n"
     "doi:10.1103/PhysRevD.93.092001\n"
     "[arXiv:1601.08061 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.93.092001;%%\n"
     "%20 citations counted in INSPIRE as of 23 Aug 2018\n");

  static ParVector<PionPhotonCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &PionPhotonCurrent::resMasses_, GeV, 5, 775.26*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PionPhotonCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &PionPhotonCurrent::resWidths_, GeV, 5, 149.1*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<PionPhotonCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &PionPhotonCurrent::amp_, 1./GeV, 5, 1./GeV, 0.0/GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<PionPhotonCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in degrees",
     &PionPhotonCurrent::phase_, 5, 0., -360.0, 360.0,
     false, false, Interface::limited);

}

// complete the construction of the decay mode for integration
bool PionPhotonCurrent::createMode(int icharge, tcPDPtr resonance,
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
    if(imode==0) {
      if(flavour.I!=IsoSpin::IOne) return false;
    }
    else {
      if(flavour.I!=IsoSpin::IOne &&
	 flavour.I!=IsoSpin::IZero) return false;
    }
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
  // check that the mode is are kinematical allowed
  Energy min = imode==0 ?
    getParticleData(ParticleID::piplus)->mass() :
    getParticleData(ParticleID::pi0   )->mass();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res;
  if(imode==0) {
    if(icharge==-3) res.push_back(getParticleData(-213));
    else            res.push_back(getParticleData( 213));
  }
  else {
    if(flavour.I==IsoSpin::IUnknown||flavour.I==IsoSpin::IOne)
       res.push_back(getParticleData(113));
    if(flavour.I==IsoSpin::IUnknown||flavour.I==IsoSpin::IZero) {
      res.push_back(getParticleData(   223));
      res.push_back(getParticleData(   333));
      res.push_back(getParticleData(100223));
      res.push_back(getParticleData( 30223));
    }
  }
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
 		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<res.size();++ix) {
    int ires(0);
    if(res[ix]->id()==213)         ires=1;
    else if(res[ix]->id()==   223) ires=2;
    else if(res[ix]->id()==100223) ires=3;
    else if(res[ix]->id()== 30223) ires=4;
    mode->resetIntermediate(res[ix],resMasses_[ires],resWidths_[ires]);
  }
  return true;
}

// the particles produced by the current
tPDVector PionPhotonCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart = {tPDPtr(),
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

void PionPhotonCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum(),
						     ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
PionPhotonCurrent::current(tcPDPtr resonance,
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
  if(flavour.I!=IsoSpin::IUnknown) {
    if(imode==0) {
      if(flavour.I!=IsoSpin::IOne) return vector<LorentzPolarizationVectorE>();
    }
    else {
      if(flavour.I!=IsoSpin::IOne &&
	 flavour.I!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
    }
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
  useMe();
  // polarization vectors of the photon
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) ++ix;
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  unsigned int imin = 0;
  unsigned int imax = imode==0 ? 1 : 5;
  if(flavour.I==IsoSpin::IOne)
    imax = 1;
  else if(flavour.I==IsoSpin::IZero) {
    imin = 1;
  }
  if(ichan>0) {
    if(flavour.I==IsoSpin::IZero)
      imin = ichan+1;
    else
      imin = ichan;
    imax=imin+1;
  }
  if(resonance) {
    switch(abs(resonance->id())) {
    case 113: case 213 :
      imin=0;
      break;
    case 223:
      imin=1;
      break;
    case 333:
      imin=2;
      break;
    case 100223:
      imin = 3;
      break;
    case 30223 :
      imin = 4;
      break;
    default:
      assert(false);
    }
    imax=imin+1;
  }
  // compute the form factor
  complex<InvEnergy> formFactor(ZERO);
  // loop over the resonances
  for(unsigned int ix=imin;ix<imax;++ix) {
    Energy2 mR2(sqr(resMasses_[ix]));
    // compute the width
    Energy width(ZERO);
    // rho
    if(ix==0) {
      width = resWidths_[0]*mR2/q2*pow(max(double((q2-4.*sqr(mpi_))/(mR2-4.*sqr(mpi_))),0.),1.5);
    }
    else {
      width = resWidths_[ix];
    }
    formFactor += couplings_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*width);
  }
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) continue;
    ret[ix] += formFactor*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool PionPhotonCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int npiplus(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::gamma)  ++ngamma;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return ngamma == 1 && (npiplus==1 || npi0==1);
}

unsigned int PionPhotonCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),npi0(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)         ++npip;
    else if(id[ix]==ParticleID::piminus)   ++npim;
    else if(id[ix]==ParticleID::pi0)       ++npi0;
    else if(id[ix]==ParticleID::gamma)   ++ngamma;
  }
  if((npip==1 || npim == 1) && ngamma==1)
    return 0;
  else
    return 1;
}

// output the information for the database
void PionPhotonCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::PionPhotonCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<resMasses_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<resWidths_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
    else     output << "insert " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<5) output << "newdef " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
