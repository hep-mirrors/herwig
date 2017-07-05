// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3OctetOctetPhotonDecayer class.
//

#include "SU3BaryonOctetOctetPhotonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonOctetOctetPhotonDecayer::SU3BaryonOctetOctetPhotonDecayer() {
  // default values of the parameters
  // these are the values for first excited multiplet
  // the couplings of the anticommutator and communtator terms
  _lf=-0.009/GeV;
  _ld=-0.024/GeV;
  // the relative parities of the two baryon multiplets
  _parity=true;
  // PDG codes for the various ground state baryons
  _proton   = 2212;
  _neutron  = 2112;
  _sigma0   = 3212;
  _sigmap   = 3222;
  _sigmam   = 3112;
  _lambda   = 3122;
  _xi0      = 3322;
  _xim      = 3312;
  // PDG codes for the various excited baryons
  _eproton  = 12212;
  _eneutron = 12112;
  _esigma0  = 13212;
  _esigmap  = 13222;
  _esigmam  = 13112;
  _elambda  = 23122;
  _exi0     = 13322;
  _exim     = 13312;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonOctetOctetPhotonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  double wgtmax;
  vector<double> wgt(0);
  for(unsigned int ix=0;ix<_incomingB.size();++ix) {
    extpart[0]=getParticleData(_incomingB[ix]);
    extpart[1]=getParticleData(_outgoingB[ix]);
    extpart[2]=getParticleData(ParticleID::gamma);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    wgtmax = _maxweight.size()>numberModes() ? 
      _maxweight[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
    // testing code
//     Energy MR = extpart[0]->mass();
//     Energy MB = extpart[1]->mass();
//     Energy kg = 0.5*(sqr(MR)-sqr(MB))/MR;
//     Energy2 Tji = 128.*sqr((sqr(MR)-sqr(MB))/2.*_prefactor[ix]/4.);
//     Energy width = 1./8./Constants::pi/sqr(MR)*kg*Tji;
//     generator()->log() << "testing " << extpart[0]->PDGName() << "->"
// 		       << extpart[1]->PDGName() << " " 
// 		       << extpart[2]->PDGName() << " " << width/MeV
// 		       << " MeV\n";
  }
}

void SU3BaryonOctetOctetPhotonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonOctetOctetPhotonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id()),id1(children[0]->id()),id2(children[1]->id()),iout;
  if(id1==ParticleID::gamma){iout=id2;}
  else if(id2==ParticleID::gamma){iout=id1;}
  else{return imode;}
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_incomingB[ix]) {
      if(iout==_outgoingB[ix]) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_incomingB[ix]) {
      if(iout==-_outgoingB[ix]) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incomingB.size()&&imode<0);
  return imode;
}

void SU3BaryonOctetOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_lf,1./GeV) << ounit(_ld,1./GeV) <<  _parity << _proton << _neutron 
     << _sigma0 << _sigmap << _sigmam << _lambda << _xi0 << _xim << _eproton 
     << _eneutron << _esigma0 << _esigmap << _esigmam << _elambda << _exi0 << _exim 
     << _incomingB << _outgoingB << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonOctetOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_lf,1./GeV) >> iunit(_ld,1./GeV) >>  _parity >> _proton >> _neutron 
     >> _sigma0 >> _sigmap >> _sigmam >> _lambda >> _xi0 >> _xim >> _eproton 
     >> _eneutron >> _esigma0 >> _esigmap >> _esigmam >> _elambda >> _exi0 >> _exim 
     >> _incomingB >> _outgoingB >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonOctetOctetPhotonDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonOctetOctetPhotonDecayer("Herwig::SU3BaryonOctetOctetPhotonDecayer", "HwBaryonDecay.so");

void SU3BaryonOctetOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetOctetPhotonDecayer> documentation
    ("The SU3BaryonOctetOctetPhotonDecayer class is designed for the "
     "radiative decay of an octet baryon to another octet baryon.");

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,InvEnergy> interfaceFcoupling
    ("Fcoupling",
     "The F coupling of the baryon resonances",
     &SU3BaryonOctetOctetPhotonDecayer::_lf, 1./GeV, -0.009/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,InvEnergy> interfaceDcoupling
    ("Dcoupling",
     "The D coupling of the baryon resonances",
     &SU3BaryonOctetOctetPhotonDecayer::_ld, 1./GeV, -0.024/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonOctetOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetOctetPhotonDecayer::_parity, true, false, false);
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

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedProton
    ("ExcitedProton",
     "The PDG code for the heavier proton-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_eproton, 12212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedNeutron
    ("ExcitedNeutron",
     "The PDG code for the heavier neutron-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_eneutron, 12112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigmap
    ("ExcitedSigma+",
     "The PDG code for the heavier Sigma+-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigmap, 13222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigma0
    ("ExcitedSigma0",
     "The PDG code for the heavier Sigma0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigma0, 13212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigmam
    ("ExcitedSigma-",
     "The PDG code for the heavier Sigma--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigmam, 13112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_elambda, 23122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedXi0
    ("ExcitedXi0",
     "The PDG code for the heavier Xi0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_exi0, 13322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedXim
    ("ExcitedXi-",
     "The PDG code for the heavier Xi--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_exim, 13312, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonOctetOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetOctetPhotonDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-1
void SU3BaryonOctetOctetPhotonDecayer::halfHalfVectorCoupling(int imode,Energy m0,
							      Energy m1, Energy,
							      Complex&A1,
							      Complex&A2,Complex&B1,
							      Complex&B2) const {
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
void SU3BaryonOctetOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1, Energy,
			    Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const {
  if(_parity) {
    A1=0.;
    B1=-_prefactor[imode]*(m0+m1);
    A2=0.;
    B2= _prefactor[imode]*(m0+m1);
  }
  else {
    A1= _prefactor[imode]*(m0-m1);
    B1=0.;
    A2= _prefactor[imode]*(m0+m1);
    B2=0.;
  }
  A3=0.;
  B3=0.;
}

// set up the decay modes
void SU3BaryonOctetOctetPhotonDecayer::setupModes(unsigned int iopt) const {
  if(_incomingB.size()!=0&&iopt==0) return;
  if(iopt==1) {
    _outgoingB.clear();
    _incomingB.clear();
  }
  // set up for the various different decay modes
  vector<InvEnergy> factor;
  vector<int> intemp,outtemp;
  // decays of the excited proton
  intemp.push_back(_eproton);outtemp.push_back(_proton);
  factor.push_back(_lf+_ld/3.);
  // decays of the excited neutron
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);
  factor.push_back(-2.*_ld/3.);
  // decays of the excited lambda
  intemp.push_back(_elambda);outtemp.push_back(_sigma0);
  factor.push_back(2.*_ld/sqrt(3.));
  intemp.push_back(_elambda);outtemp.push_back(_lambda);
  factor.push_back(-_ld/3.);
  // decays of the excited sigma+
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);
  factor.push_back(_lf+_ld/3.);
  // decays of the excited sigma0
  intemp.push_back(_esigma0);outtemp.push_back(_sigma0);
  factor.push_back(_lf/3.);
  intemp.push_back(_esigma0);outtemp.push_back(_lambda);
  factor.push_back(2.*_ld/sqrt(3.));
  // decays of the excited simga-
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);
  factor.push_back(-_ld-_lf/3.);
  // decays of the excited xi-
  intemp.push_back(_exim);outtemp.push_back(_xim);
  factor.push_back(_ld/3.-_lf);
  // decays of the excited xi0
  intemp.push_back(_exi0);outtemp.push_back(_xi0);
  factor.push_back(-2.*_ld/3.);
  int inspin,outspin;
  tPDVector extpart(2);
  for(unsigned int ix=0;ix<intemp.size();++ix) {
    if(intemp[ix]!=0&&outtemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()) {
	_incomingB.push_back(intemp[ix]);
	_outgoingB.push_back(outtemp[ix]);
	if(iopt==1) {
	  inspin  = extpart[0]->iSpin();
	  outspin = extpart[1]->iSpin();
	  factor[ix] *=2.;
	  if(inspin==2&&outspin==2)      _prefactor.push_back(2.*factor[ix]);
	  else if(inspin==4&&outspin==2) _prefactor.push_back(factor[ix]);
	  else throw DecayIntegratorError() 
	    << "Invalid combination of spins in "
	    << "SU3BaryonOctetOctetPhotonDecayer::" 
	    << "setupModes()" << Exception::abortnow;}
      }
    }
  }
}

void SU3BaryonOctetOctetPhotonDecayer::dataBaseOutput(ofstream & output,
						      bool header) const {
  if(header) output << "update decayers set parameters=\""; 
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Fcoupling " << _lf*GeV << "\n";
  output << "newdef " << name() << ":Dcoupling " << _ld*GeV << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Proton " << _proton << "\n";
  output << "newdef " << name() << ":Neutron " << _neutron << "\n";
  output << "newdef " << name() << ":Sigma+ " << _sigmap << "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Sigma- " << _sigmam << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":Xi0 " << _xi0 << "\n";
  output << "newdef " << name() << ":Xi- " << _xim << "\n"; 
  output << "newdef " << name() << ":ExcitedProton " << _eproton << "\n";
  output << "newdef " << name() << ":ExcitedNeutron " << _eneutron << "\n";
  output << "newdef " << name() << ":ExcitedSigma+ " << _esigmap << "\n";
  output << "newdef " << name() << ":ExcitedSigma0 " << _esigma0 << "\n";
  output << "newdef " << name() << ":ExcitedSigma- " << _esigmam << "\n";
  output << "newdef " << name() << ":ExcitedLambda " << _elambda << "\n";
  output << "newdef " << name() << ":ExcitedXi0 " << _exi0 << "\n";
  output << "newdef " << name() << ":ExcitedXi- " << _exim << "\n"; 
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
