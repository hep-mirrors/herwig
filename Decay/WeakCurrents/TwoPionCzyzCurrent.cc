// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoPionCzyzCurrent class.
//

#include "TwoPionCzyzCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

TwoPionCzyzCurrent::TwoPionCzyzCurrent()
  : omegaMag_(18.7e-4), omegaPhase_(0.106),
    omegaMass_(782.4*MeV),omegaWidth_(8.33*MeV), beta_(2.148),
    nMax_(1000) {
  // various parameters
  rhoMag_  =  {1.,1.,0.59,0.048,0.40,0.43};
  rhoPhase_ = {0.,0.,-2.20,-2.0,-2.9,1.19}; 
  rhoMasses_ = {773.37*MeV,1490*MeV, 1870*MeV,2120*MeV,2321*MeV,2567*MeV};
  rhoWidths_ = { 147.1*MeV, 429*MeV,  357*MeV, 300*MeV, 444*MeV, 491*MeV};
  // set up for the modes in the base class
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(3);
}

IBPtr TwoPionCzyzCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoPionCzyzCurrent::fullclone() const {
  return new_ptr(*this);
}

void TwoPionCzyzCurrent::persistentOutput(PersistentOStream & os) const {
  os << beta_ << omegaWgt_ << omegaMag_ << omegaPhase_
     << ounit(omegaMass_,GeV) << ounit(omegaWidth_,GeV)
     << rhoWgt_ << rhoMag_ << rhoPhase_
     << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << ounit(mass_,GeV) << ounit(width_,GeV) << coup_
     << dh_ << ounit(hres_,GeV2) << ounit(h0_,GeV2) << nMax_;
}

void TwoPionCzyzCurrent::persistentInput(PersistentIStream & is, int) {
  is >> beta_ >> omegaWgt_ >> omegaMag_ >> omegaPhase_
     >> iunit(omegaMass_,GeV) >> iunit(omegaWidth_,GeV)
     >> rhoWgt_ >> rhoMag_ >> rhoPhase_
     >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> iunit(mass_,GeV) >> iunit(width_,GeV) >> coup_
     >> dh_ >> iunit(hres_,GeV2) >> iunit(h0_,GeV2) >> nMax_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoPionCzyzCurrent,WeakDecayCurrent>
  describeHerwigTwoPionCzyzCurrent("Herwig::TwoPionCzyzCurrent", "HwWeakCurrents.so");

void TwoPionCzyzCurrent::Init() {

  static ClassDocumentation<TwoPionCzyzCurrent> documentation
    ("The TwoPionCzyzCurrent class uses the currents from "
     "PRD 81 094014 for the weak current with two pions",
     "The current for two pions from \\cite{Czyz:2010hj} was used.",
     "%\\cite{Czyz:2010hj}\n"
     "\\bibitem{Czyz:2010hj}\n"
     "H.~Czyz, A.~Grzelinska and J.~H.~Kuhn,\n"
     "%``Narrow resonances studies with the radiative return method,''\n"
     "Phys.\\ Rev.\\ D {\\bf 81} (2010) 094014\n"
     "doi:10.1103/PhysRevD.81.094014\n"
     "[arXiv:1002.0279 [hep-ph]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.81.094014;%%\n"
     "%28 citations counted in INSPIRE as of 30 Jul 2018\n");

  static ParVector<TwoPionCzyzCurrent,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoMasses_, MeV, -1, 775.8*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<TwoPionCzyzCurrent,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoWidths_, MeV, -1, 150.3*MeV, ZERO, 1000.*MeV,
     false, false, true);
  
  static ParVector<TwoPionCzyzCurrent,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "Magnitude of the weight of the different resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoMag_, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<TwoPionCzyzCurrent,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase of the weight of the different resonances for the pi pi channel",
     &TwoPionCzyzCurrent::rhoPhase_, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static Parameter<TwoPionCzyzCurrent,unsigned int> interfacenMax
    ("nMax",
     "The maximum number of resonances to include in the sum,"
     " should be approx infinity",
     &TwoPionCzyzCurrent::nMax_, 1000, 10, 10000,
     false, false, Interface::limited);
  
  static Parameter<TwoPionCzyzCurrent,double> interfacebeta
    ("beta",
     "The beta parameter for the couplings",
     &TwoPionCzyzCurrent::beta_, 2.148, 0.0, 100.,
     false, false, Interface::limited);
  
  static Parameter<TwoPionCzyzCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &TwoPionCzyzCurrent::omegaMass_, MeV,782.4*MeV, 0.0*MeV, 100.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<TwoPionCzyzCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The mass of the omega meson",
     &TwoPionCzyzCurrent::omegaWidth_, MeV, 8.33*MeV, 0.0*MeV, 100.0*MeV,
     false, false, Interface::limited);

  static Parameter<TwoPionCzyzCurrent,double> interfaceOmegaMagnitude
    ("OmegaMagnitude",
     "The magnitude of the omega couplings",
     &TwoPionCzyzCurrent::omegaMag_, 18.7e-4, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TwoPionCzyzCurrent,double> interfaceOmegaPhase
    ("OmegaPhase",
     "The magnitude of the omega couplings",
     &TwoPionCzyzCurrent::omegaPhase_, 0.106, 0.0, 2.*Constants::pi,
     false, false, Interface::limited);

}

void TwoPionCzyzCurrent::doinit() {
  WeakDecayCurrent::doinit();
  // check consistency of parametrers
  if(rhoMasses_.size()!=rhoWidths_.size())
    throw InitException() << "Inconsistent parameters in TwoPionCzyzCurrent"
			  << "::doinit()" << Exception::abortnow;
  // weights for the rho channels
  if(rhoMag_.size()!=rhoPhase_.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
  			  << " rho channel must be the same size in "
  			  << "TwoPionCzyzCurrent::doinit()" << Exception::runerror;
  Complex rhoSum(0.);
  for(unsigned int ix=0;ix<rhoMag_.size();++ix) {
    rhoWgt_.push_back(rhoMag_[ix]*(cos(rhoPhase_[ix])+Complex(0.,1.)*sin(rhoPhase_[ix])));
    if(ix>0) rhoSum +=rhoWgt_.back();
  }
  omegaWgt_ = omegaMag_*(cos(omegaPhase_)+Complex(0.,1.)*sin(omegaPhase_));
  // set up the masses and widths of the rho resonances
  double gamB(tgamma(2.-beta_));
  Complex cwgt(0.);
  Energy mpi(getParticleData(ParticleID::piplus)->mass());
  for(unsigned int ix=0;ix<nMax_;++ix) {
    // this is gam(2-beta+n)/gam(n+1)
    if(ix>0) {
      gamB *= ((1.-beta_+double(ix)))/double(ix);
    }
    Complex c_n = tgamma(beta_-0.5) /(0.5+double(ix)) / sqrt(Constants::pi) *
      sin(Constants::pi*(beta_-1.-double(ix)))/Constants::pi*gamB;
    if(ix%2!=0) c_n *= -1.;
    // set the masses and widths
    // calc for higher resonances
    if(ix>=rhoMasses_.size()) {
      mass_ .push_back(rhoMasses_[0]*sqrt(1.+2.*double(ix)));
      width_.push_back(rhoWidths_[0]/rhoMasses_[0]*mass_.back());
    }
    // input for lower ones
    else {
      mass_ .push_back(rhoMasses_[ix]);
      width_.push_back(rhoWidths_[ix]);
      if(ix>0) cwgt += c_n;
    }
    // parameters for the gs propagators
    hres_.push_back(Resonance::Hhat(sqr(mass_.back()),mass_.back(),width_.back(),mpi,mpi));
    dh_  .push_back(Resonance::dHhatds(mass_.back(),width_.back(),mpi,mpi));
    h0_.push_back(Resonance::H(ZERO,mass_.back(),width_.back(),mpi,mpi,dh_.back(),hres_.back()));
    coup_.push_back(c_n);
  }
  // fix up the early weights
  for(unsigned int ix=1;ix<rhoMasses_.size();++ix) {
    coup_[ix] = rhoWgt_[ix]*cwgt/rhoSum;
  }
}

// complete the construction of the decay mode for integration
bool TwoPionCzyzCurrent::createMode(int icharge, unsigned int imode,
					 DecayPhaseSpaceModePtr mode,
					 unsigned int iloc,unsigned int,
					 DecayPhaseSpaceChannelPtr phase,Energy upp) {
  if((imode==0 && abs(icharge)!=3) ||
     (imode>0  && icharge !=0)) return false; 
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else  {
    part[0]=getParticleData(ParticleID::piplus);
    part[1]=getParticleData(ParticleID::piminus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  DecayPhaseSpaceChannelPtr newchannel;
  // set the resonances
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
    if(res[ix]) {
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(res[ix],0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
    }
  }
  // reset the masses in the intergrators
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix<rhoMasses_.size()&&res[ix]) {
      mode->resetIntermediate(res[ix],rhoMasses_[ix],rhoWidths_[ix]);
    }
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector TwoPionCzyzCurrent::particles(int icharge, unsigned int imode,
					     int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::piminus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
TwoPionCzyzCurrent::current(const int imode, const int ichan,
			    Energy & scale,const ParticleVector & outpart,
			    DecayIntegrator::MEOption meopt) const {
  useMe();
  if(meopt==DecayIntegrator::Terminate) {
    for(unsigned int ix=0;ix<2;++ix)
      ScalarWaveFunction::constructSpinInfo(outpart[ix],outgoing,true);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(outpart[0]->momentum()-outpart[1]->momentum());
  Lorentz5Momentum psum (outpart[0]->momentum()+outpart[1]->momentum());
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  LorentzPolarizationVector vect;
  // calculate the current
  // compute the current
  pdiff-=psum;
  Energy ma = outpart[0]->mass();
  Energy mb = outpart[1]->mass();
  Complex FPI=Fpi(q2,imode,ichan,ma,mb);
  return vector<LorentzPolarizationVectorE>(1,FPI*pdiff);
}
   
bool TwoPionCzyzCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // pion modes
  if((abs(id[0])==ParticleID::piplus  &&     id[1] ==ParticleID::pi0   ) ||
     (    id[0] ==ParticleID::pi0     && abs(id[1])==ParticleID::piplus))
    return true;
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::piplus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::piminus))
    return true;
  else
    return false;
}

// the decay mode
unsigned int TwoPionCzyzCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==2) return 1;
  else       return 0;
}

// output the information for the database
void TwoPionCzyzCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::TwoPionCzyzCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<rhoMasses_.size();++ix) {
    if(ix<6)  output << "newdef ";
    else      output << "insert ";
    output << name() << ":RhoMasses " << ix << " " << rhoMasses_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWidths_.size();++ix) {
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoWidths " << ix << " " << rhoWidths_[ix]/MeV << "\n";
  }
  for(ix=0;ix<rhoWgt_.size();++ix) {
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoMagnitude " << ix << " " << rhoMag_[ix]   << "\n";
    if(ix<6) output << "newdef ";
    else     output << "insert ";
    output << name() << ":RhoPhase "     << ix << " " << rhoPhase_[ix] << "\n";
  }
  output << "newdef " << name() << ":OmegaMass " << omegaMass_/MeV << "\n";
  output << "newdef " << name() << ":OmegaWidth " << omegaWidth_/MeV << "\n";
  output << "newdef " << name() << ":OmegaMagnitude " << omegaMag_ << "\n";
  output << "newdef " << name() << ":OmegaPhase " << omegaPhase_ << "\n";
  output << "newdef " << name() << ":nMax " << nMax_ << "\n";
  output << "newdef " << name() << ":beta " << beta_ << "\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

Complex TwoPionCzyzCurrent::Fpi(Energy2 q2,const int imode, const int ichan,
				Energy ma, Energy mb) const {
  Complex FPI(0.);
  for(unsigned int ix=0;ix<coup_.size();++ix) {
    if(ichan>=0&&ix!=abs(ichan)) continue;
    Complex term = coup_[ix]*Resonance::BreitWignerGS(q2,mass_[ix],width_[ix],
						      ma,mb,h0_[ix],dh_[ix],hres_[ix]);
    // include rho-omega if needed
    if(ix==0&&imode!=0)
      term *= 1./(1.+omegaWgt_)*(1.+omegaWgt_*Resonance::BreitWignerFW(q2,omegaMass_,omegaWidth_));
    FPI += term;
  }
  // factor for cc mode
  if(imode==0)           FPI *= sqrt(2.0);
  return FPI;
}