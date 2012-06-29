// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonOctetOctetScalarDecayer class.
//

#include "SU3BaryonOctetOctetScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SU3BaryonOctetOctetScalarDecayer::SU3BaryonOctetOctetScalarDecayer() {
  // default values of the parameters
  // these are the values for first excited multiplet
  // the couplings of the anticommutator and communtator terms
  _sf= 0.11;
  _sd= 0.60;
  // the relative parities of the two baryon multiplets
  _parity=true;
  // the pion decay constant
  _fpi=130.7*MeV;
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

void SU3BaryonOctetOctetScalarDecayer::doinit() {
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
    extpart[2]=getParticleData(_outgoingM[ix]);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    wgtmax = 
      _maxweight.size()>numberModes() ? _maxweight[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
    // testing code
//     Energy MR=extpart[0]->mass();
//     Energy MB=extpart[1]->mass();
//     Energy Mp=extpart[2]->mass();
//     Energy kp = 0.5/MR*sqrt((sqr(MR)-sqr(MB+Mp))*(sqr(MR)-sqr(MB-Mp)));
//     Energy width;
//     if(_parity) {
//       width = sqr((MR+MB)/MR)*kp/(8.*Constants::pi)*
// 	sqr(_prefactor[ix])*(sqr(MR-MB)-sqr(Mp));
//     }
//     else {
//     width  = sqr((MR-MB)/MR)*kp/(8.*Constants::pi)*
// 	sqr(_prefactor[ix])*(sqr(MR+MB)-sqr(Mp));
//     } 
//     generator()->log() << "Partial Width for "
// 		       << extpart[0]->PDGName() << "->"
// 		       << extpart[1]->PDGName() << " "
// 		       << extpart[2]->PDGName() << " "
// 		       << width/GeV << "\n";
  }
}
void SU3BaryonOctetOctetScalarDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) 
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonOctetOctetScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_incomingB[ix]) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_incomingB[ix]) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) {
	imode=ix;
	cc=true;
      }
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incomingB.size()&&imode<0);
  return imode;
}

void SU3BaryonOctetOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _sf << _sd << _parity << ounit(_fpi,GeV) << _proton << _neutron 
     << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _eproton << _eneutron << _esigma0 
     << _esigmap << _esigmam << _elambda << _exi0 << _exim << _incomingB << _outgoingB 
     << _outgoingM << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonOctetOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _sf >> _sd >> _parity >> iunit(_fpi,GeV) >> _proton >> _neutron 
     >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _eproton >> _eneutron >> _esigma0 
     >> _esigmap >> _esigmam >> _elambda >> _exi0 >> _exim >> _incomingB >> _outgoingB 
     >> _outgoingM >> _maxweight >> iunit(_prefactor,1./GeV);
}

ClassDescription<SU3BaryonOctetOctetScalarDecayer> 
SU3BaryonOctetOctetScalarDecayer::initSU3BaryonOctetOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonOctetOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetOctetScalarDecayer> documentation
    ("The SU3BaryonOctetOctetScalarDecayer class is designed for the"
     " decay of excited baryon resonances assuming SU(3) symmetry");

  static Parameter<SU3BaryonOctetOctetScalarDecayer,double> interfaceFcoupling
    ("Fcoupling",
     "The F coupling of the baryon resonances",
     &SU3BaryonOctetOctetScalarDecayer::_sf, 0.60, -20.0, 20.0,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,double> interfaceDcoupling
    ("Dcoupling",
     "The D coupling of the baryon resonances",
     &SU3BaryonOctetOctetScalarDecayer::_sd, 0.11, -20.0, 20.0,
     false, false, true);

  static Switch<SU3BaryonOctetOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetOctetScalarDecayer::_parity, true, false, false);
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

  static Parameter<SU3BaryonOctetOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonOctetOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedProton
    ("ExcitedProton",
     "The PDG code for the heavier proton-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_eproton, 12212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedNeutron
    ("ExcitedNeutron",
     "The PDG code for the heavier neutron-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_eneutron, 12112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigmap
    ("ExcitedSigma+",
     "The PDG code for the heavier Sigma+-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigmap, 13222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigma0
    ("ExcitedSigma0",
     "The PDG code for the heavier Sigma0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigma0, 13212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigmam
    ("ExcitedSigma-",
     "The PDG code for the heavier Sigma--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigmam, 13112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_elambda, 23112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedXi0
    ("ExcitedXi0",
     "The PDG code for the heavier Xi0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_exi0, 13322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedXim
    ("ExcitedXi-",
     "The PDG code for the heavier Xi--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_exim, 13312, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonOctetOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 10000., false, false, true);
}

// couplings for spin-1/2 to spin-1/2 spin-0
void SU3BaryonOctetOctetScalarDecayer::halfHalfScalarCoupling(int imode,Energy m0,
							      Energy m1,Energy,
							      Complex& A,
							      Complex& B) const {
  if(_parity) {
    A=0.;
    B=_prefactor[imode]*(m0+m1);
  }
  else {
    A=_prefactor[imode]*(m0-m1);
    B=0.;
  }
}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonOctetOctetScalarDecayer::threeHalfHalfScalarCoupling(int imode,Energy m0,
								   Energy m1,Energy,
								   Complex& A,
								   Complex& B) const {
  if(_parity) {
    A=_prefactor[imode]*(m0+m1);
    B=0.;
  }
  else {
    A=0.;
    B=_prefactor[imode]*(m0+m1);
  }
}

// set up the decay modes
void SU3BaryonOctetOctetScalarDecayer::setupModes(unsigned int iopt) const {
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1) {
    _outgoingB.clear();
    _incomingB.clear();
    _outgoingM.clear();
  }
  // set up for the various different decay modes
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.));
  // decays of the excited proton
  intemp.push_back(_eproton);outtemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back((_sd+_sf));
  intemp.push_back(_eproton);outtemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(-ort*(_sd+_sf));
  intemp.push_back(_eproton);outtemp.push_back(_sigmap);mestemp.push_back(311);
  factor.push_back(_sd-_sf);
  intemp.push_back(_eproton);outtemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_eproton);outtemp.push_back(_proton);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_eproton);outtemp.push_back(_lambda);mestemp.push_back(321);
  factor.push_back(-ors*(_sd+3.*_sf));
  // decays of the excited neutron
  intemp.push_back(_eneutron);outtemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back((_sd+_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_sigmam);mestemp.push_back(321);
  factor.push_back(_sd-_sf);
  intemp.push_back(_eneutron);outtemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_eneutron);outtemp.push_back(_lambda);mestemp.push_back(311);
  factor.push_back(-ors*(_sd+3.*_sf));
  // decays of the excited lambda
  intemp.push_back(_elambda);outtemp.push_back(_sigma0);mestemp.push_back(111);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_lambda);mestemp.push_back(221);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_elambda);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_elambda);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(-ors*(3.*_sf+_sd));
  intemp.push_back(_elambda);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(-ors*(3.*_sf+_sd));
  // decays of the excited sigma+
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_xi0);mestemp.push_back(321);
  factor.push_back(_sd+_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_proton);mestemp.push_back(-311);
  factor.push_back(_sd-_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(2.*ors*_sd);
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(2.*ors*_sd);
  // decays of the excited sigma0
  intemp.push_back(_esigma0);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigma0);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigma0);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_sigma0);mestemp.push_back(221);
  factor.push_back(2.*ors*_sd);
  intemp.push_back(_esigma0);outtemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(2.*ors*_sd);
  // decays of the excited simga-
  intemp.push_back(_esigmam);outtemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_xim);mestemp.push_back(311);
  factor.push_back(_sd+_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_neutron);mestemp.push_back(-321);
  factor.push_back(_sd-_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(2.*_sd*ors);
  // decays of the excited xi-
  intemp.push_back(_exim);outtemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_sd+_sf);
  intemp.push_back(_exim);outtemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_exim);outtemp.push_back(_xi0);mestemp.push_back(-211);
  factor.push_back(_sd-_sf);
  intemp.push_back(_exim);outtemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_exim);outtemp.push_back(_lambda);mestemp.push_back(-321);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_exim);outtemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf+_sd));
  // decays of the excited xi0
  intemp.push_back(_exi0);outtemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(_sd+_sf);
  intemp.push_back(_exi0);outtemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_exi0);outtemp.push_back(_xim);mestemp.push_back(211);
  factor.push_back(_sd-_sf);
  intemp.push_back(_exi0);outtemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_exi0);outtemp.push_back(_lambda);mestemp.push_back(-311);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_exi0);outtemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf+_sd));
  int inspin,outspin;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<intemp.size();++ix) {
    if(intemp[ix]!=0&&outtemp[ix]!=0&&mestemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      extpart[2]=getParticleData(mestemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin()) {
	_incomingB.push_back(intemp[ix]);
	_outgoingB.push_back(outtemp[ix]);
	_outgoingM.push_back(mestemp[ix]);
	if(iopt==1) {
	  inspin  = extpart[0]->iSpin();
	  outspin = extpart[1]->iSpin();
	  if(inspin==2&&outspin==2)
	    _prefactor.push_back(ort*factor[ix]/_fpi);
	  else if(inspin==4&&outspin==2)
	    _prefactor.push_back(ort*factor[ix]/_fpi);
	  else
	    throw DecayIntegratorError()<< "Invalid combination of spins in "
					<< "SU3BaryonOctetOctetScalarDecayer::" 
					<< "setupModes()" 
					<< Exception::abortnow;
	}
      }
    }
  }
}
 
void SU3BaryonOctetOctetScalarDecayer::dataBaseOutput(ofstream & output,
						      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Fcoupling " << _sf << "\n";
  output << "newdef " << name() << ":Dcoupling " << _sd << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
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
