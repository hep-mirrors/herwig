// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPhotonCurrent class.
//

#include "EtaPhotonCurrent.h"
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

EtaPhotonCurrent::EtaPhotonCurrent() {
  // modes handled
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {0.77526*GeV,0.78284*GeV,1.01952*GeV,1.465*GeV,1.70*GeV};
  // widths for the resonances
  resWidths_ = {0.1491 *GeV,0.00868*GeV,0.00421*GeV,0.40*GeV,0.30*GeV};
  // amplitudes
  amp_   = {0.0861/GeV,0.00824/GeV,0.0158/GeV,0.0147/GeV,ZERO};
  // phases
  phase_ = {0.,11.3,170.,61.,0.};
}

IBPtr EtaPhotonCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaPhotonCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaPhotonCurrent::doinit() {
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

void EtaPhotonCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV)
     << ounit(mpi_,GeV);
}

void EtaPhotonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV)
     >> iunit(mpi_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPhotonCurrent,WeakCurrent>
describeHerwigEtaPhotonCurrent("Herwig::EtaPhotonCurrent",
				"HwWeakCurrents.so");

void EtaPhotonCurrent::Init() {

  static ClassDocumentation<EtaPhotonCurrent> documentation
    ("The EtaPhotonCurrent class implements a current based"
     " on the model of SND for pion+photon",
     "The current based on the model of \\cite{Achasov:2006dv"
     " for eta and photon was used.",
     "\\bibitem{Achasov:2006dv}\n"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Study of the e+ e- ---> eta gamma process with SND detector at the VEPP-2M e+ e- collider,''\n"
     "Phys.\\ Rev.\\ D {\\bf 74} (2006) 014016\n"
     "doi:10.1103/PhysRevD.74.014016\n"
     "[hep-ex/0605109].\n"
     "%%CITATION = doi:10.1103/PhysRevD.74.014016;%%\n"
     "%25 citations counted in INSPIRE as of 23 Aug 2018\n");


  static ParVector<EtaPhotonCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &EtaPhotonCurrent::resMasses_, GeV, 5, 775.26*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhotonCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &EtaPhotonCurrent::resWidths_, GeV, 5, 149.1*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhotonCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &EtaPhotonCurrent::amp_, 1./GeV, 5, 1./GeV, 0.0/GeV, 100./GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhotonCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in degrees",
     &EtaPhotonCurrent::phase_, 5, 0., 0.0, 360.0,
     false, false, Interface::limited);

}

// complete the construction of the decay mode for integration
bool EtaPhotonCurrent::createMode(int icharge, tcPDPtr resonance,
				   IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
				   unsigned int, PhaseSpaceModePtr mode,
				   unsigned int iloc,int ires,
				   PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IZero) return false;
  }
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    if(i3!=IsoSpin::I3Zero) return false;
  }
  // check that the mode is are kinematical allowed
  Energy min = getParticleData(ParticleID::eta)->mass();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(113),
		   getParticleData(   223),
		   getParticleData(   333),
		   getParticleData(100223),
		   getParticleData( 30223)};
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
 		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses and widths of the resonances if needed
  for(unsigned int ix=0;ix<res.size();++ix) {
    mode->resetIntermediate(res[ix],resMasses_[ix],resWidths_[ix]);
  }
  return true;
}

// the particles produced by the current
tPDVector EtaPhotonCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  return {getParticleData(ParticleID::eta),getParticleData(ParticleID::gamma)};
}

void EtaPhotonCurrent::constructSpinInfo(ParticleVector decay) const {
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
EtaPhotonCurrent::current(tcPDPtr resonance,
			   IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			   const int, const int ichan, Energy & scale, 
			   const tPDVector & ,
			   const vector<Lorentz5Momentum> & momenta,
			   DecayIntegrator::MEOption) const {
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
  }
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    if(i3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
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
  unsigned int imax = couplings_.size();
  if(ichan>0) {
    imin = ichan;
    imax = imin+1;
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

bool EtaPhotonCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int neta(0),ngamma(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::eta)   ++neta;
    else if(id[ix]==ParticleID::gamma) ++ngamma;
  }
  return ngamma == 1 && neta==1;
}

unsigned int EtaPhotonCurrent::decayMode(vector<int>) {
  return 0;
}

// output the information for the database
void EtaPhotonCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPhotonCurrent " << name() 
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
