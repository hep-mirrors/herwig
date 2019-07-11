// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPhiCurrent class.
//

#include "EtaPhiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EtaPhiCurrent::EtaPhiCurrent() {
  addDecayMode(3,-3);
  setInitialModes(3);
  // Masses for the resonances
  resMasses_ = {1.670*GeV,2.14*GeV};
  // widths for the resonances
  resWidths_ = {122*MeV,43.5*MeV};
  // amplitudes
  amp_   = {0.175/GeV,0.00409/GeV};
  // phases
  phase_ = {0.,2.19};
}

IBPtr EtaPhiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaPhiCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaPhiCurrent::doinit() {
  WeakCurrent::doinit();
  assert(phase_.size()==amp_.size());
  couplings_.clear();
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    couplings_.push_back(amp_[ix]*(cos(phase_[ix])+ii*sin(phase_[ix])));
  }
}

void EtaPhiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(resMasses_,GeV) << ounit(resWidths_,GeV)
     << ounit(amp_,1./GeV) << phase_ << ounit(couplings_,1./GeV);
}

void EtaPhiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(resMasses_,GeV) >> iunit(resWidths_,GeV)
     >> iunit(amp_,1./GeV) >> phase_ >> iunit(couplings_,1./GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPhiCurrent,WeakCurrent>
describeHerwigEtaPhiCurrent("Herwig::EtaPhiCurrent",
			    "HwWeakCurrents.so");

void EtaPhiCurrent::Init() {

  static ClassDocumentation<EtaPhiCurrent> documentation
    ("The EtaPhiCurrent class implements a current based"
     " on the model of SND for eta + phi "
     "The current based on the model of \\cite{Achasov:2018ygm}"
     " for eta and phi was used.",
     "\\bibitem{Achasov:2018ygm}\n"
     "M.~N.~Achasov {\\it et al.},\n"
     "%``Measurement of the $e^+e^−\\to\\eta K^+K^−$ Cross Section by Means of the SND Detector,''\n"
     "Phys.\\ Atom.\\ Nucl.\\  {\\bf 81} (2018) no.2,  205\n"
     " [Yad.\\ Fiz.\\  {\\bf 81} (2018) no.2,  195].\n"
     "doi:10.1134/S1063778818020023\n"
     "%%CITATION = doi:10.1134/S1063778818020023;%%\n");

  static ParVector<EtaPhiCurrent,Energy> interfaceResonanceMasses
    ("ResonanceMasses",
     "The masses of the resonances for the form factor",
     &EtaPhiCurrent::resMasses_, GeV, 1, 1680*MeV, 0.5*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhiCurrent,Energy> interfaceResonanceWidths
    ("ResonanceWidths",
     "The widths of the resonances for the form factor",
     &EtaPhiCurrent::resWidths_, GeV, 1, 150*MeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhiCurrent,InvEnergy> interfaceAmplitude
    ("Amplitude",
     "The amplitudes of the couplings",
     &EtaPhiCurrent::amp_, 0.00115/GeV, 1, 1./GeV, 0.0/GeV, 10/GeV,
     false, false, Interface::limited);

  static ParVector<EtaPhiCurrent,double> interfacePhase
    ("Phase",
     "The phases of the couplings in radians",
     &EtaPhiCurrent::phase_, 1, 0., 0.0, 2.*Constants::pi,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool EtaPhiCurrent::createMode(int icharge, tcPDPtr resonance,
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
               getParticleData(ParticleID::phi)->massMin();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(100333)};
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
tPDVector EtaPhiCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  return {getParticleData(ParticleID::eta),getParticleData(ParticleID::phi)};
}

void EtaPhiCurrent::constructSpinInfo(ParticleVector decay) const {
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
EtaPhiCurrent::current(tcPDPtr resonance,
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
    case 100333:
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

bool EtaPhiCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  unsigned int neta(0),nphi(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::eta)   ++neta;
    else if(id[ix]==ParticleID::phi) ++nphi;
  }
  return nphi == 1 && neta==1;
}

unsigned int EtaPhiCurrent::decayMode(vector<int>) {
  return 0;
}

// output the information for the database
void EtaPhiCurrent::dataBaseOutput(ofstream & output,bool header,
				   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPhiCurrent " << name() 
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
