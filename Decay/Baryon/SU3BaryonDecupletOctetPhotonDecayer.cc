// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonDecupletOctetPhotonDecayer class.
// 

#include "SU3BaryonDecupletOctetPhotonDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SU3BaryonDecupletOctetPhotonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  double wgtmax;
  vector<double> wgt(0);
  for(unsigned int ix=0;ix<incomingB_.size();++ix) {
    extpart[0]=getParticleData(incomingB_[ix]);
    extpart[1]=getParticleData(outgoingB_[ix]);
    extpart[2]=getParticleData(ParticleID::gamma);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    wgtmax= maxweight_.size()>numberModes() ? 
      maxweight_[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
  }
}

void SU3BaryonDecupletOctetPhotonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

SU3BaryonDecupletOctetPhotonDecayer::SU3BaryonDecupletOctetPhotonDecayer() {
  // the coupling
  C_=1.0/GeV;
  // the relative parities of the two baryon multiplets
  parity_=true;
  // PDG codes for the various octet baryons
  proton_   = 2212;
  neutron_  = 2112;
  sigma0_   = 3212;
  sigmap_   = 3222;
  sigmam_   = 3112;
  lambda_   = 3122;
  xi0_      = 3322;
  xim_      = 3312;
  // PDG codes for the various decuplet baryons
  deltapp_  = 2224;
  deltap_   = 2214;
  delta0_   = 2114;
  deltam_   = 1114;
  sigmasp_  = 3224;
  sigmas0_  = 3214;
  sigmasm_  = 3114;
  omega_    = 3334;
  xism_     = 3314;
  xis0_     = 3324;
  // intermediates
  generateIntermediates(false);
}

int SU3BaryonDecupletOctetPhotonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  if(incomingB_.size()==0) setupModes(0); 
  // must be two outgoing particles
  if(children.size()!=2||(children[0]->id()!=ParticleID::gamma&&
			  children[1]->id()!=ParticleID::gamma)) return -1;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  int iout = id1==ParticleID::gamma ? id2 : id1;
  unsigned int ix(0);
  cc=false;
  int imode(-1);
  do {
    if(id0==incomingB_[ix]) {
      if(iout==outgoingB_[ix]) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-incomingB_[ix]) {
      if(iout==-outgoingB_[ix]) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incomingB_.size()&&imode<0);
  return imode;
}

void SU3BaryonDecupletOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(C_,1./GeV) << parity_ << proton_ << neutron_ << sigma0_ << sigmap_ 
     << sigmam_ << lambda_ << xi0_ << xim_ << deltapp_ << deltap_ << delta0_ << deltam_
     << sigmasp_ << sigmas0_ << sigmasm_ << omega_ << xism_ << xis0_ << incomingB_ 
     << outgoingB_ << maxweight_ << ounit(prefactor_,1./GeV);
}

void SU3BaryonDecupletOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(C_,1./GeV) >> parity_ >> proton_ >> neutron_ >> sigma0_ >> sigmap_ 
     >> sigmam_ >> lambda_ >> xi0_ >> xim_ >> deltapp_ >> deltap_ >> delta0_ >> deltam_
     >> sigmasp_ >> sigmas0_ >> sigmasm_ >> omega_ >> xism_ >> xis0_ >> incomingB_ 
     >> outgoingB_ >> maxweight_ >> iunit(prefactor_,1./GeV);
}

ClassDescription<SU3BaryonDecupletOctetPhotonDecayer> 
SU3BaryonDecupletOctetPhotonDecayer::initSU3BaryonDecupletOctetPhotonDecayer;
// Definition of the static class description member.

void SU3BaryonDecupletOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetPhotonDecayer> documentation
    ("The SU3BaryonDecupletOctetPhotonDecayer class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a photon.");

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,InvEnergy> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetPhotonDecayer::C_, 1.0/GeV, 1.0/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetPhotonDecayer::parity_, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "The multiplets have the same parity.",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "The multiplets have different parities.",
     false);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::proton_, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::neutron_, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmap_, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigma0_, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmam_, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::lambda_, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xi0_, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xim_, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::deltapp_, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::deltap_, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::delta0_, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::deltam_, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmasp_, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmas0_, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmasm_, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::omega_, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xis0_, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xism_, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetPhotonDecayer::maxweight_,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonDecupletOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const {
  A3=0.;
  B3=0.;
  if(parity_) {
    A1 = 0.;
    B1 = -prefactor_[imode]*(m0+m1);
    A2 = 0.;
    B2 = prefactor_[imode]*(m0+m1);
  }
  else {
    A1= prefactor_[imode]*(m0-m1);
    B1=0.;
    A2= prefactor_[imode]*(m0+m1);
    B2=0.;
  }
}

// set up the decay modes
void SU3BaryonDecupletOctetPhotonDecayer::setupModes(unsigned int iopt) const {
  if(incomingB_.size()!=0&&iopt==0) return;
  if(iopt==1) {
    outgoingB_.clear();
    incomingB_.clear();
  }
  vector<InvEnergy> factor;
  vector<int> intemp,outtemp;
  double ortw(1./sqrt(12.)),orr(1./sqrt(3.));
  // decays of the delta+
  intemp.push_back(deltap_);outtemp.push_back(proton_);
  factor.push_back(C_*orr);
  // decays of the delta0
  intemp.push_back(delta0_);outtemp.push_back(neutron_);
  factor.push_back(C_*orr);
  // sigma*+
  intemp.push_back(sigmasp_);outtemp.push_back(sigmap_);
  factor.push_back(-C_*orr);
  // sigma*0
  intemp.push_back(sigmas0_);outtemp.push_back(lambda_);
  factor.push_back(-C_*.5);
  intemp.push_back(sigmas0_);outtemp.push_back(sigma0_);
  factor.push_back(C_*ortw);
  // xi*0
  intemp.push_back(xis0_);outtemp.push_back(xi0_);
  factor.push_back(-C_*orr);
  // set up the modes
  tPDVector extpart(2);
  for(unsigned int ix=0;ix<intemp.size();++ix) {
    if(intemp[ix]!=0&&outtemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()) {
	incomingB_.push_back(intemp[ix]);
	outgoingB_.push_back(outtemp[ix]);
	if(iopt==1) prefactor_.push_back(factor[ix]);
      }
    }
  }
}
void SU3BaryonDecupletOctetPhotonDecayer::dataBaseOutput(ofstream & output,
							 bool header) const {
  if(header) output << "update decayers set parameters=\""; 
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Ccoupling " << C_*GeV<< "\n";
  output << "newdef " << name() << ":Parity " << parity_<< "\n";
  output << "newdef " << name() << ":Proton " << proton_ << "\n";
  output << "newdef " << name() << ":Neutron " << neutron_ << "\n";
  output << "newdef " << name() << ":Sigma+ " << sigmap_ << "\n";
  output << "newdef " << name() << ":Sigma0 " << sigma0_ << "\n";
  output << "newdef " << name() << ":Sigma- " << sigmam_ << "\n";
  output << "newdef " << name() << ":Lambda " << lambda_ << "\n";
  output << "newdef " << name() << ":Xi0 " << xi0_ << "\n";
  output << "newdef " << name() << ":Xi- " << xim_ << "\n";
  output << "newdef " << name() << ":Delta++ " << deltapp_ << "\n";
  output << "newdef " << name() << ":Delta+ " << deltap_ << "\n";
  output << "newdef " << name() << ":Delta0 " << delta0_ << "\n";
  output << "newdef " << name() << ":Delta- " << deltam_ << "\n";
  output << "newdef " << name() << ":Sigma*+ " << sigmasp_ << "\n";
  output << "newdef " << name() << ":Sigma*0 " << sigmas0_ << "\n";
  output << "newdef " << name() << ":Sigma*- " << sigmasm_ << "\n";
  output << "newdef " << name() << ":Omega " << omega_ << "\n";
  output << "newdef " << name() << ":Xi*0 " << xis0_ << "\n";
  output << "newdef " << name() << ":Xi*- " << xism_ << "\n";
  for(unsigned int ix=0;ix<maxweight_.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
