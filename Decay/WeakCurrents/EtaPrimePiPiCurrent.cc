// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPrimePiPiCurrent class.
//

#include "EtaPrimePiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EtaPrimePiPiCurrent::EtaPrimePiPiCurrent() : fpi_(93.3*MeV) {
  rhoMasses_ = {0.77549*GeV,1.54*GeV ,1.76*GeV,2.11*GeV};
  rhoWidths_ = {0.1494 *GeV,0.356*GeV,.113*GeV,.176*GeV};
  amp_       = {1.,0.,0.,0.02};
  phase_     = {0.,Constants::pi,Constants::pi,Constants::pi};
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
}

IBPtr EtaPrimePiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr EtaPrimePiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void EtaPrimePiPiCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(rhoMasses_.size()   != rhoWidths_.size())
    throw InitException() << "Inconsistent parameters in EtaPrimePiPiCurrent"
			  << "::doinit()" << Exception::abortnow;
  // weights for the rho channels
  if(amp_.size()!=phase_.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
  			  << " rho channel must be the same size in "
  			  << "EtaPrimePiPiCurrent::doinit()" << Exception::runerror;
  // combine mags and phase
  weights_.clear();
  for(unsigned int ix=0;ix<amp_.size();++ix) {
    weights_.push_back(amp_[ix]*(cos(phase_[ix])+Complex(0.,1.)*sin(phase_[ix])));
  }
  Complex denom = std::accumulate(weights_.begin(),weights_.end(),Complex(0.));
  for(unsigned int ix=0;ix<weights_.size();++ix)
    weights_[ix] /=denom;
}

void EtaPrimePiPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << weights_ << amp_ << phase_ << ounit(fpi_,MeV)
     << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV);
}

void EtaPrimePiPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> weights_ >> amp_ >> phase_ >> iunit(fpi_,MeV)
     >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPrimePiPiCurrent,WeakCurrent>
describeHerwigEtaPrimePiPiCurrent("Herwig::EtaPrimePiPiCurrent", "HwWeakCurrents.so");

void EtaPrimePiPiCurrent::Init() {

  static ClassDocumentation<EtaPrimePiPiCurrent> documentation
    ("There is no documentation for the EtaPrimePiPiCurrent class");

  static ParVector<EtaPrimePiPiCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::rhoMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<EtaPrimePiPiCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::rhoWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<EtaPrimePiPiCurrent,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::amp_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<EtaPrimePiPiCurrent,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &EtaPrimePiPiCurrent::phase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Parameter<EtaPrimePiPiCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &EtaPrimePiPiCurrent::fpi_, MeV, 93.3*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool EtaPrimePiPiCurrent::createMode(int icharge, tcPDPtr resonance,
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
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero)
    return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       )
    return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero       )
    return false;
  // make sure that the decays are kinematically allowed
  int iq(0),ia(0);
  tPDVector part = particles(icharge,imode,iq,ia);
  Energy min=ZERO;
  for(tPDPtr p : part) min += p->massMin();
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
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],ires+1,res[0],ires+1,iloc+3,
		      ires+2,iloc+1,ires+2,iloc+2));
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()&&res[ix]) {
      mode->resetIntermediate(res[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  return true;
}

// the particles produced by the current
tPDVector EtaPrimePiPiCurrent::particles(int icharge, unsigned int imode,
				    int,int) {
  tPDVector output(3);
  output[0]=getParticleData(ParticleID::piplus);
  output[2]=getParticleData(ParticleID::etaprime);
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

// hadronic current
vector<LorentzPolarizationVectorE> 
EtaPrimePiPiCurrent::current(tcPDPtr resonance,
			FlavourInfo flavour,
			const int imode, const int ichan,Energy & scale, 
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
  Lorentz5Momentum q=momenta[0]+momenta[1]+momenta[2];
  q.rescaleMass();
  Energy2 s1 = (momenta[0]+momenta[1]).m2();
  Energy2 Q2 = q.mass2();
  Energy  Q  = q.mass();
  Complex BW1 = Resonance::BreitWignerPWave(s1,rhoMasses_[0],rhoWidths_[0],
					    momenta[0].mass(),momenta[1].mass());
  vector<Complex> BWs = {Resonance::BreitWignerPWave(Q2,rhoMasses_[0],rhoWidths_[0],
						     momenta[0].mass(),momenta[1].mass()),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[1],
						  rhoWidths_[1]*pow(Q/rhoMasses_[1],3)),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[2],
						  rhoWidths_[2]*pow(Q/rhoMasses_[2],3)),
			 Resonance::BreitWignerFW(Q2,rhoMasses_[3],
						  rhoWidths_[3]*pow(Q/rhoMasses_[3],3))};
  unsigned int imin=0,imax=4;
  if(resonance) {
    switch(resonance->id()/1000) {
    case 0:
      imax = 1;
      break;
    case 100:
      imin = 1;
      imax = 2;
      break;
    case 30 :
      imin = 2;
      imax = 3;
      break;
    default:
      assert(false);
    }
  }
  if(ichan>0&&ichan!=3) {
    imin = ichan;
    imax = ichan+1;
  }
  else if(ichan==3) {
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  // form factor
  Complex fact(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    fact += weights_[ix]*BWs[ix];
  }
  fact *= -0.25*Complex(0.,1.)/sqr(Constants::pi)*sqrt(2./3.)*BW1;
  if(imode==0) fact *=sqrt(2.);
  scale=Q;
  LorentzPolarizationVectorE output = fact/pow<3,1>(fpi_)*Q*
    Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  return vector<LorentzPolarizationVectorE>(1,output);
}
   
bool EtaPrimePiPiCurrent::accept(vector<int> id) {
  // check there are only three particles
  if(id.size()!=3) return false;
  unsigned int npip(0),npim(0),npi0(0),neta(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::etaprime)     ++neta;
  }
  if( (npip==1&&npim==1&&neta==1) ||
      (npi0==1&&npim+npip==1&&neta==1))
    return true;
  else
    return false;
}

// the decay mode
unsigned int EtaPrimePiPiCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==2) return 1;
  else       return 0;
}

// output the information for the database
void EtaPrimePiPiCurrent::dataBaseOutput(ofstream & output,bool header,
				    bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::EtaPrimePiPiCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<3)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<weights_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMagnitude " << ix << " " << amp_[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoPhase "     << ix << " " << phase_[ix] << "\n";
  }
  output << "newdef " << name() << ":FPi " << fpi_/MeV << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

