// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaOmegaCurrent class.
//

#include "EtaOmegaCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EtaOmegaCurrent::EtaOmegaCurrent() {
  addDecayMode(3,-3);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {1.425*GeV,1.67*GeV};
  // widths for the resonances
  resWidths_ = {215*MeV  , 114*MeV};
  // amplitudes
  amp_   = {0.0863/GeV,0.0655/GeV};
  // phases
  phase_ = {0.,180.};
}

IBPtr EtaOmegaCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaOmegaCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaOmegaCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  couplings_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    double phi = phase_[ix]/180.*Constants::pi;
    couplings_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  }
}

void EtaOmegaCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV);
}

void EtaOmegaCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaOmegaCurrent,WeakCurrent>
describeHerwigEtaOmegaCurrent("Herwig::EtaOmegaCurrent",
			    "HwWeakCurrents.so");

void EtaOmegaCurrent::Init() {

  static ClassDocumentation<EtaOmegaCurrent> documentation
    ("The EtaOmegaCurrent class implements a current based"
     " on the model of SND for eta + omega "
     "The current based on the model of \\cite{Achasov:2016qvd}"
     " for eta and omega was used.",
     "\\bibitem{Achasov:2016qvd}\n"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Measurement of the $e^+e^- \\to \\omega\\eta$ cross section below $\\sqrt{s}=2$ GeV,''\n"
     "Phys.\\ Rev.\\ D {\\bf 94} (2016) no.9,  092002\n"
     "doi:10.1103/PhysRevD.94.092002\n"
     "[arXiv:1607.00371 [hep-ex]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.94.092002;%%\n"
     "%18 citations counted in INSPIRE as of 12 Oct 2018\n");

  static ParVector<EtaOmegaCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &EtaOmegaCurrent::resMasses_, GeV, 1, 1680*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaOmegaCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &EtaOmegaCurrent::resWidths_, GeV, 1, 150*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaOmegaCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &EtaOmegaCurrent::amp_, 0.00115/GeV, 1, 1./GeV, 0.0/GeV, 10/GeV,
     false, false, Interface::limited);

  static ParVector<EtaOmegaCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in degrees",
     &EtaOmegaCurrent::phase_, 1, 0., 0.0, 360.0,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool EtaOmegaCurrent::createMode(int icharge, tcPDPtr resonance,
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
  Energy min = getParticleData(ParticleID::eta)->mass()+
               getParticleData(ParticleID::omega)->massMin();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(100223),getParticleData(30223)};
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
tPDVector EtaOmegaCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  return {getParticleData(ParticleID::eta),getParticleData(ParticleID::omega)};
}

void EtaOmegaCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum(),
						     ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
EtaOmegaCurrent::current(tcPDPtr resonance,
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
    case 100223:
      imin=0;
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
    Energy width = resWidths_[ix];
    formFactor += couplings_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*width);
  }
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] += formFactor*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  return ret;
}

bool EtaOmegaCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int neta(0),nomega(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::eta)   ++neta;
    else if(id[ix]==ParticleID::omega) ++nomega;
  }
  return nomega == 1 && neta==1;
}

unsigned int EtaOmegaCurrent::decayMode(vector<int>) {
  return 0;
}

// output the information for the database
void EtaOmegaCurrent::dataBaseOutput(ofstream & output,bool header,
				   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaOmegaCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<resMasses_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceMasses " << ix 
		    << " " << resMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<resWidths_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
    else     output << "insert " << name() << ":ResonanceWidths " << ix 
		    << " " << resWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
    else     output << "insert " << name() << ":Amplitude " << ix 
		    << " " << amp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<phase_.size();++ix) {
    if(ix<1) output << "newdef " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
    else     output << "insert " << name() << ":Phase " << ix 
		    << " " << phase_[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
