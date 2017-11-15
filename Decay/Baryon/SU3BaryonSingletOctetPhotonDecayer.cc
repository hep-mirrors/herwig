// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonSingletOctetPhotonDecayer class.
//

#include "SU3BaryonSingletOctetPhotonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonSingletOctetPhotonDecayer::SU3BaryonSingletOctetPhotonDecayer() {
  // the coupling
  _c=0.252/GeV;
  // the relative parities of the two baryon multiplets
  _parity=false;
  // PDG codes for the various ground state baryons
  _sigma0   = 3212;
  _lambda   = 3122;
  // PDG codes for the excited baryon
  _elambda  = 3124;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonSingletOctetPhotonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  double wgtmax;
  vector<double> wgt(0);
  for(unsigned int ix=0;ix<_outgoingB.size();++ix) {
    extpart[0]=getParticleData(_elambda);
    extpart[1]=getParticleData(_outgoingB[ix]);
    extpart[2]=getParticleData(ParticleID::gamma);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    if(_maxweight.size()>numberModes()){wgtmax=_maxweight[numberModes()];}
    else{wgtmax=1.;}
    addMode(mode,wgtmax,wgt);
  }
}

void SU3BaryonSingletOctetPhotonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonSingletOctetPhotonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_outgoingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id()),iout;
  if(id1==ParticleID::gamma){iout=id2;}
  else if(id2==ParticleID::gamma){iout=id1;}
  else{return imode;}
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_elambda){
      if(iout==_outgoingB[ix]) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_elambda) {
      if(iout==-_outgoingB[ix]) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_outgoingB.size()&&imode<0);
  return imode;
}

void SU3BaryonSingletOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_c,1./GeV) << _parity << _sigma0 << _lambda << _elambda << _outgoingB 
     << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonSingletOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_c,1./GeV) >> _parity >> _sigma0 >> _lambda >> _elambda >> _outgoingB 
     >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonSingletOctetPhotonDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonSingletOctetPhotonDecayer("Herwig::SU3BaryonSingletOctetPhotonDecayer", "HwBaryonDecay.so");

void SU3BaryonSingletOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonSingletOctetPhotonDecayer> documentation
    ("The SU3BaryonSingletOctetPhotonDecayer class performs the decay"
     " of a singlet baryon to an octet baryon and a photon.");

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The C coupling of the baryon resonances.",
     &SU3BaryonSingletOctetPhotonDecayer::_c, 1./GeV, 0.252/GeV, -10./GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonSingletOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonSingletOctetPhotonDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_elambda, 3124, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonSingletOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonSingletOctetPhotonDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-1
void SU3BaryonSingletOctetPhotonDecayer::
halfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
		       Complex&A1,Complex&A2,Complex&B1,Complex&B2) const {
  if(_parity) {
    A1=    _prefactor[imode]*(m0+m1);
    B1=0.;
    A2=-2.*_prefactor[imode]*(m0+m1);
    B2=0.;
  }
  else {
    A1=0.;
    B1=    _prefactor[imode]*(m1-m0);
    A2=0.;
    B2=-2.*_prefactor[imode]*(m0+m1);
  }
}

// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonSingletOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const {
  A3=0.;B3=0.;
  if(_parity) {
    A1=0.;
    B1=-_prefactor[imode]*(m0+m1);
    A2=0.;
    B2= _prefactor[imode]*(m0+m1);
  }
  else {
    A1=_prefactor[imode]*(m0-m1);B1=0.;
    A2=_prefactor[imode]*(m0+m1);B2=0.;
  }
}

// set up the decay modes
void SU3BaryonSingletOctetPhotonDecayer::setupModes(unsigned int iopt) const {
  if(_outgoingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.clear();}
  // set up for the various different decay modes
  vector<int> outtemp;
  vector<InvEnergy> factor;
  if(_elambda==0)
    throw DecayIntegratorError() << "Invalid incoming particle in "
				 << "SU3BaryonSingletOctetScalarDecayer::" 
				 << "setupModes()" << Exception::abortnow;
  // decays of the excited lambda
  outtemp.push_back(_sigma0);factor.push_back(_c/sqrt(2.));
  outtemp.push_back(_lambda);factor.push_back(_c/sqrt(6.));
  tPDVector extpart(2);extpart[0]=getParticleData(_elambda);
  int inspin(extpart[0]->iSpin()),outspin;
  for(unsigned int ix=0;ix<outtemp.size();++ix) {
    if(outtemp[ix]!=0) {
      extpart[1]=getParticleData(outtemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()) {
	_outgoingB.push_back(outtemp[ix]);
	if(iopt==1) {
	  outspin = extpart[1]->iSpin();
	  factor[ix] *=2.;
	  if(inspin==2&&outspin==2)      _prefactor.push_back(2.*factor[ix]);
	  else if(inspin==4&&outspin==2) _prefactor.push_back(   factor[ix]);
	  else
	    throw DecayIntegratorError() 
	      << "Invalid combination of spins in "
	      << "SU3BaryonSingletOctetScalarDecayer::" 
	      << "setupModes()" << Exception::abortnow;
	}
      }
    }
  }
}

void SU3BaryonSingletOctetPhotonDecayer::dataBaseOutput(ofstream & output,
							bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Coupling " << _c*GeV << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":ExcitedLambda " << _elambda << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
