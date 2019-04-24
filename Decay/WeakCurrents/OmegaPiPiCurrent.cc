// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPiPiCurrent class.
//

#include "OmegaPiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

OmegaPiPiCurrent::OmegaPiPiCurrent() : mRes_(1.62*GeV) {
  wRes_ = 0.288*GeV;
  gRes_ = 2.83;
  mSigma_ = 0.6*GeV;
  wSigma_ = 1.0*GeV;
  mf0_ = 0.98*GeV;
  gPiPi_ = 0.331;
  gKK_ = 0.144;
  gf0_ = 0.85;
}

IBPtr OmegaPiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr OmegaPiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void OmegaPiPiCurrent::persistentOutput(PersistentOStream & os) const {
}

void OmegaPiPiCurrent::persistentInput(PersistentIStream & is, int) {
}

void OmegaPiPiCurrent::doinit() {
  WeakCurrent::doinit();
  // assert(phase_.size()==amp_.size());
  // couplings_.clear();
  // Complex ii(0.,1.);
  // for(unsigned int ix=0;ix<amp_.size();++ix) {
  //   double phi = phase_[ix]/180.*Constants::pi;
  //   couplings_.push_back(amp_[ix]*(cos(phi)+ii*sin(phi)));
  // }
}


// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<OmegaPiPiCurrent,WeakCurrent>
describeHerwigOmegaPiPiCurrent("Herwig::OmegaPiPiCurrent", "HwWeakCurrents.so");

void OmegaPiPiCurrent::Init() {

  static ClassDocumentation<OmegaPiPiCurrent> documentation
    ("There is no documentation for the OmegaPiPiCurrent class");

}



// complete the construction of the decay mode for integration
bool OmegaPiPiCurrent::createMode(int icharge, tcPDPtr resonance,
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
  Energy min = getParticleData(ParticleID::omega)->massMin()+
    2.*getParticleData(ParticleID::pi0)->mass();
  if(min>upp) return false;
  // resonances for the intermediate channels
  tPDVector res = {getParticleData(30223)};
  tPDVector res2 = {getParticleData(9000221),getParticleData(9010221)};
  // set up the integration channels;
  for(unsigned int ix=0;ix<res.size();++ix) {
    if(resonance && resonance!=res[ix]) continue;
    for(unsigned int iy=0;iy<res2.size();++iy) {
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
   		      ires+1,iloc+1,ires+1,res2[iy],ires+2,iloc+2,ires+2,iloc+3));
    }
  }
  return true;
}

// the particles produced by the current
tPDVector OmegaPiPiCurrent::particles(int icharge, unsigned int imode,int,int) {
  assert(icharge==0 && imode<=1);
  if(imode==0) 
    return {getParticleData(ParticleID::omega),
	    getParticleData(ParticleID::piplus),
	    getParticleData(ParticleID::piminus)};
  else if(imode==1) 
    return {getParticleData(ParticleID::omega),
	getParticleData(ParticleID::pi0),
	getParticleData(ParticleID::pi0)};
}

void OmegaPiPiCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[0]->momentum(),
  						     ix,Helicity::outgoing);
  }
  VectorWaveFunction::constructSpinInfo(temp,decay[0],
  					outgoing,true,true);
  for(unsigned int ix=1;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
OmegaPiPiCurrent::current(tcPDPtr resonance,
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
  // polarization vectors of the omega
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  // total momentum of the system
  Lorentz5Momentum q(momenta[0]+momenta[1]+momenta[2]);
  // overall hadronic mass
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2(q.m2());
  // resonance factor for s channel resonance
  Energy2 mR2=sqr(mRes_);
  Complex pre= mR2*gRes_/(q2-mR2 + Complex(0.,1.)*scale*wRes_);
  // compute the form factor
  complex<Energy> formFactor(ZERO);

  formFactor = pre*scale*gSigma_;

  // needs to be multiplied by thing inside || in 1.9
  
  //cm energy for intermediate f0 channel
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  
  //sigma meson
  Energy2 mSigma2 = sqr(mSigma_);
  Complex Sigma_form = mSigma2/(s1-mSigma2 + Complex(0.,1.)*scale*wSigma_);
  
  //f0(980) following Phys. Lett. B 63, 224 (1976)
  Energy2 mf02 = sqr(mf0_);
  Energy2 mPi2 = momenta[1].m2();
  Energy pcm_pipi = 0.5*sqrt(s1-4*mPi2);
  Energy pcm_f0 = 0.5*sqrt(mf0_*mf0_-4*mPi2);
  Energy GPiPi = gPiPi_*pcm_pipi;
  Energy G0 = gPiPi_*pcm_f0;
  Energy m_kaon = getParticleData(311)->mass();
  Energy2 m_kaon2 = sqr(m_kaon);
  Energy GKK;
  Complex f0_form;
  if(0.25*s1>m_kaon2){
    GKK = gKK_*sqrt(0.25*s1-m_kaon2);
    f0_form = gf0_*mf0_*sqrt(G0*GPiPi)/(s1-mf02+Complex(0.,1.)*mf0_*(GPiPi+GKK));
  }
  else{
    GKK = gKK_*sqrt(m_kaon2-0.25*s1);
    f0_form = gf0_*mf0_*sqrt(G0*GPiPi)/(s1-mf02+Complex(0.,1.)*mf0_*(GPiPi+Complex(0.,1.)*GKK));
  }
  formFactor *=(Sigma_form+f0_form);

//   // loop over the resonances
//   for(unsigned int ix=imin;ix<imax;++ix) {
//     Energy2 mR2(sqr(resMasses_[ix]));
//     // compute the width
//     Energy width = resWidths_[ix];
//     formFactor += couplings_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*width);
//   }
  
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] = formFactor*temp[ix];
  }
  return ret;
}

bool OmegaPiPiCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
   unsigned int nomega(0),npip(0),npim(0),npi0(0);
   for(unsigned int ix=0;ix<id.size();++ix) {
     if(abs(id[ix])==ParticleID::piminus) ++npim;
     if(abs(id[ix])==ParticleID::pi0    ) ++npi0;
     if(abs(id[ix])==ParticleID::piplus ) ++npip;
     else if(id[ix]==ParticleID::omega)  ++nomega;
   }
   return nomega == 1 && (npi0==2 || (npip==1&&npim==1));
}

unsigned int OmegaPiPiCurrent::decayMode(vector<int> id) {
   unsigned int npi0(0);
   for(unsigned int ix=0;ix<id.size();++ix) {
     if(abs(id[ix])==ParticleID::pi0    ) ++npi0;
   }
   return npi0==2 ? 1 : 0;
}

// output the information for the database
void OmegaPiPiCurrent::dataBaseOutput(ofstream & output,bool header,
				   bool create) const {
  // if(header) output << "update decayers set parameters=\"";
  // if(create) output << "create Herwig::OmegaPiPiCurrent " << name() 
  // 		    << " HwWeakCurrents.so\n";
  // for(unsigned int ix=0;ix<resMasses_.size();++ix) {
  //   if(ix<1) output << "newdef " << name() << ":ResonanceMasses " << ix 
  // 		    << " " << resMasses_[ix]/GeV << "\n";
  //   else     output << "insert " << name() << ":ResonanceMasses " << ix 
  // 		    << " " << resMasses_[ix]/GeV << "\n";
  // }
  // for(unsigned int ix=0;ix<resWidths_.size();++ix) {
  //   if(ix<1) output << "newdef " << name() << ":ResonanceWidths " << ix 
  // 		    << " " << resWidths_[ix]/GeV << "\n";
  //   else     output << "insert " << name() << ":ResonanceWidths " << ix 
  // 		    << " " << resWidths_[ix]/GeV << "\n";
  // }
  // for(unsigned int ix=0;ix<amp_.size();++ix) {
  //   if(ix<1) output << "newdef " << name() << ":Amplitude " << ix 
  // 		    << " " << amp_[ix]*GeV << "\n";
  //   else     output << "insert " << name() << ":Amplitude " << ix 
  // 		    << " " << amp_[ix]*GeV << "\n";
  // }
  // for(unsigned int ix=0;ix<phase_.size();++ix) {
  //   if(ix<1) output << "newdef " << name() << ":Phase " << ix 
  // 		    << " " << phase_[ix] << "\n";
  //   else     output << "insert " << name() << ":Phase " << ix 
  // 		    << " " << phase_[ix] << "\n";
  // }
  // WeakCurrent::dataBaseOutput(output,false,false);
  // if(header) output << "\n\" where BINARY ThePEGName=\"" 
  // 		    << fullName() << "\";" << endl;
}
