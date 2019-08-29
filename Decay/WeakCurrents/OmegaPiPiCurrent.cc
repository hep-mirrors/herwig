// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaPiPiCurrent class.
//

#include "OmegaPiPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

OmegaPiPiCurrent::OmegaPiPiCurrent() :mRes_(1.62*GeV) { 
  wRes_ = 0.288*GeV;
  gRes_ = 2.83;
  mSigma_ = 0.6*GeV;
  wSigma_ = 1.0*GeV;
  gSigma_ = 1.0;
  mf0_ = 0.98*GeV;
  gPiPi_ = 0.331;
  gKK_ = 0.144;
  gf0_ = 0.85;
  addDecayMode(1,-1);
  addDecayMode(1,-1);
  setInitialModes(2);
}

IBPtr OmegaPiPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr OmegaPiPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void OmegaPiPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(mRes_,GeV) << ounit(wRes_,GeV) << gRes_
     << ounit(mSigma_,GeV) << ounit(wSigma_,GeV) << ounit(mf0_,GeV)
     << gPiPi_ << gKK_ << gSigma_ << gf0_;
}

void OmegaPiPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mRes_,GeV) >> iunit(wRes_,GeV) >> gRes_
     >> iunit(mSigma_,GeV) >> iunit(wSigma_,GeV) >> iunit(mf0_,GeV)
     >> gPiPi_ >> gKK_ >> gSigma_ >> gf0_;
}

void OmegaPiPiCurrent::doinit() {
  WeakCurrent::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<OmegaPiPiCurrent,WeakCurrent>
describeHerwigOmegaPiPiCurrent("Herwig::OmegaPiPiCurrent", "HwWeakCurrents.so");

void OmegaPiPiCurrent::Init() {

  static ClassDocumentation<OmegaPiPiCurrent> documentation
    ("The OmegaPiPiCurrent class provides the current for I=0 omega pi pi");
  
  static Parameter<OmegaPiPiCurrent,Energy> interfacemRes
    ("mRes",
     "The mass of the s-channel resonance",
     &OmegaPiPiCurrent::mRes_, GeV, 1.62*GeV, 0.*GeV, 10.*GeV,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy> interfacewRes
    ("wRes",
     "The width of the s-channel resonance",
     &OmegaPiPiCurrent::wRes_, GeV, 0.288*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,double> interfacegRes
    ("gRes",
     "The coupling of the s-channel resonance",
     &OmegaPiPiCurrent::gRes_, 2.83, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,Energy> interfacemSigma
    ("mSigma",
     "The mass of the Sigma",
     &OmegaPiPiCurrent::mSigma_, GeV, 0.6*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,Energy> interfacewSigma
    ("wSigma",
     "The width of the Sigma",
     &OmegaPiPiCurrent::wSigma_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,double> interfacegSigma
    ("gSigma",
     "The coupling of the Sigma resonance",
     &OmegaPiPiCurrent::gSigma_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,Energy> interfacemf0
    ("mf0",
     "The mass of the f_0(980)",
     &OmegaPiPiCurrent::mf0_, GeV, 0.98*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<OmegaPiPiCurrent,double> interfacegf0
    ("gf0",
     "The coupling of the f_0(980) meson",
     &OmegaPiPiCurrent::gf0_, 0.85, 0.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,double> interfacegPiPi
    ("gPiPi",
     "The coupling of the f_0(980) to pipi",
     &OmegaPiPiCurrent::gPiPi_, .331, 0.0, 10.0,
     false, false, Interface::limited);
  
  static Parameter<OmegaPiPiCurrent,double> interfacegKK
    ("gKK",
     "The coupling of the f_0(980) to KK",
     &OmegaPiPiCurrent::gKK_, .144, 0.0, 10.0,
     false, false, Interface::limited);
  
}



// complete the construction of the decay mode for integration
bool OmegaPiPiCurrent::createMode(int icharge, tcPDPtr resonance,
			       FlavourInfo flavour,
			       unsigned int, PhaseSpaceModePtr mode,
			       unsigned int iloc,int ires,
			       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IZero) return false;
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown && flavour.I3!=IsoSpin::I3Zero) return false;
  // and other flavour
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
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
  else
    assert(false);
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
			  FlavourInfo flavour,
			  const int, const int ichan, Energy & scale, 
			  const tPDVector & ,
			  const vector<Lorentz5Momentum> & momenta,
			  DecayIntegrator::MEOption) const {
  // no isospin/flavour here
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IZero) return vector<LorentzPolarizationVectorE>();
  if(flavour.I3!=IsoSpin::I3Unknown && flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
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

  //formFactor = pre*scale*gSigma_;

  // needs to be multiplied by thing inside || in 1.9
  
  //cm energy for intermediate f0 channel
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy sqrs1 = sqrt(s1);
  
  //sigma meson
  Energy2 mSigma2 = sqr(mSigma_);
  Complex Sigma_form = mSigma2/(s1-mSigma2 + Complex(0.,1.)*sqrs1*wSigma_);
  
  //f0(980) following Phys. Lett. B 63, 224 (1976)
  Energy2 mf02 = sqr(mf0_);
  Energy mPi = getParticleData(211)->mass();
  Energy2 mPi2 = sqr(mPi);
  Energy pcm_pipi = 0.5*sqrt(s1-4*mPi2);
  Energy pcm_f0 = 0.5*sqrt(mf0_*mf0_-4*mPi2);
  Energy GPiPi = gPiPi_*pcm_pipi;
  Energy G0 = gPiPi_*pcm_f0;
  Energy m_kaon0 = getParticleData(311)->mass();
  Energy2 m_kaon02 = sqr(m_kaon0);
  Energy m_kaonP = getParticleData(321)->mass();
  Energy2 m_kaonP2 = sqr(m_kaonP);
  Energy GKK;
  Complex f0_form;
  if(0.25*s1>m_kaon02){
    GKK = 0.5*gKK_*(sqrt(0.25*s1-m_kaon02)+sqrt(0.25*s1-m_kaonP2));
    f0_form = gf0_*mf0_*sqrt(G0*GPiPi)/(s1-mf02+Complex(0.,1.)*mf0_*(GPiPi+GKK));
  }
  else{
    if(0.25*s1>m_kaonP2){
      f0_form = gf0_*mf0_*sqrt(G0*GPiPi)/(s1-mf02+Complex(0.,1.)*mf0_*(GPiPi+0.5*gKK_*(Complex(0.,1.)*sqrt(m_kaon02-0.25*s1)+sqrt(0.25*s1-m_kaonP2))));
    }
    else{
      GKK = 0.5*gKK_*(sqrt(m_kaon02-0.25*s1)+sqrt(m_kaonP2-0.25*s1));
      f0_form = gf0_*mf0_*sqrt(G0*GPiPi)/(s1-mf02+Complex(0.,1.)*mf0_*(GPiPi+Complex(0.,1.)*GKK));
    }
  }
  formFactor =(Sigma_form+f0_form)*pre*scale*gSigma_;

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
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::OmegaPiPiCurrent " << name() 
  		    << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":mRes "   << " " << mRes_/GeV   << "\n";
  output << "newdef " << name() << ":wRes "   << " " << wRes_/GeV   << "\n";
  output << "newdef " << name() << ":mSigma " << " " << mSigma_/GeV << "\n";
  output << "newdef " << name() << ":wSigma " << " " << wSigma_/GeV << "\n";
  output << "newdef " << name() << ":mf0 "    << " " << mf0_/GeV    << "\n";
  output << "newdef " << name() << ":gRes "   << " " << gRes_       << "\n";
  output << "newdef " << name() << ":gSigma " << " " << gSigma_     << "\n";
  output << "newdef " << name() << ":gf0 "    << " " << gf0_        << "\n";
  output << "newdef " << name() << ":gPiPi "  << " " << gPiPi_      << "\n";
  output << "newdef " << name() << ":gKK "    << " " << gKK_        << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
   		    << fullName() << "\";" << endl;
}
