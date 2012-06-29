// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonSingletOctetScalarDecayer class.
//

#include "SU3BaryonSingletOctetScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonSingletOctetScalarDecayer::SU3BaryonSingletOctetScalarDecayer() {
  // the coupling
  _c=0.39;
  // the relative parities of the two baryon multiplets
  _parity=false;
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
  // PDG codes for the excited baryon
  _elambda  = 13122;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonSingletOctetScalarDecayer::doinit() {
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
    extpart[2]=getParticleData(_outgoingM[ix]);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    wgtmax = _maxweight.size()>numberModes() 
      ? _maxweight[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
  }
}

void SU3BaryonSingletOctetScalarDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonSingletOctetScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
						   const tPDVector & children) const {
  int imode(-1);
  if(_outgoingB.size()==0) setupModes(0);
  // must be two outgoing particles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_elambda) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_elambda) {
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
  while(ix<_outgoingB.size()&&imode<0);
  return imode;
}

void SU3BaryonSingletOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _c << _parity << ounit(_fpi,GeV) << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _elambda << _outgoingB 
     << _outgoingM << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonSingletOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _c >> _parity >> iunit(_fpi,GeV) >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _elambda >> _outgoingB 
     >> _outgoingM >> _maxweight >> iunit(_prefactor,1./GeV);
}

ClassDescription<SU3BaryonSingletOctetScalarDecayer> SU3BaryonSingletOctetScalarDecayer::initSU3BaryonSingletOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonSingletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonSingletOctetScalarDecayer> documentation
    ("The SU3BaryonSingletOctetScalarDecayer class is designed for"
     "the decay of an excited SU(3) singlet baryon");

  static Parameter<SU3BaryonSingletOctetScalarDecayer,double> interfaceCcoupling
    ("Coupling",
     "The C coupling of the baryon resonances",
     &SU3BaryonSingletOctetScalarDecayer::_c, 0.39, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonSingletOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonSingletOctetScalarDecayer::_parity, true, false, false);
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

  static Parameter<SU3BaryonSingletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonSingletOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_elambda, 13122, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonSingletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonSingletOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void SU3BaryonSingletOctetScalarDecayer::
halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy,
		       Complex& A, Complex& B) const {
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
void SU3BaryonSingletOctetScalarDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0, Energy m1,Energy,
			    Complex& A, Complex& B) const {
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
void SU3BaryonSingletOctetScalarDecayer::setupModes(unsigned int iopt) const {
  if(_outgoingB.size()!=0&&iopt==0) return;
  if(iopt==1) {
    _outgoingB.clear();
    _outgoingM.clear();
  }
  // set up for the various different decay modes
  vector<int> outtemp,mestemp;
  double rt(sqrt(2.));
  if(_elambda==0)
    throw DecayIntegratorError() << "Invalid incoming particle in "
				 << "SU3BaryonSingletOctetScalarDecayer::" 
				 << "setupModes()" << Exception::abortnow;
  // decays of the excited lambda
  outtemp.push_back(_sigma0);mestemp.push_back(111);
  outtemp.push_back(_sigmap);mestemp.push_back(-211);
  outtemp.push_back(_sigmam);mestemp.push_back(211);
  outtemp.push_back(_lambda);mestemp.push_back(221);
  outtemp.push_back(_xim);mestemp.push_back(321);
  outtemp.push_back(_xi0);mestemp.push_back(311);
  outtemp.push_back(_proton);mestemp.push_back(-321);
  outtemp.push_back(_neutron);mestemp.push_back(-311);
  tPDVector extpart(3);
  extpart[0]=getParticleData(_elambda);
  int inspin(extpart[0]->iSpin()),outspin;
  for(unsigned int ix=0;ix<outtemp.size();++ix) {
    if(outtemp[ix]!=0&&mestemp[ix]!=0) {
      extpart[1]=getParticleData(outtemp[ix]);
      extpart[2]=getParticleData(mestemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin()) {
	_outgoingB.push_back(outtemp[ix]);
	_outgoingM.push_back(mestemp[ix]);
	if(iopt==1) {
	  outspin = extpart[1]->iSpin();
	  if(inspin==2&&outspin==2)      _prefactor.push_back(_c*rt/_fpi);
	  else if(inspin==4&&outspin==2) _prefactor.push_back(_c*rt/_fpi);
	  else throw DecayIntegratorError() 
	    << "Invalid combination of spins in "
	    << "SU3BaryonSingletOctetScalarDecayer::" 
	    << "setupModes()" << Exception::abortnow;
	}
      }
    }
  }
}

void SU3BaryonSingletOctetScalarDecayer::dataBaseOutput(ofstream & output,
							bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Coupling " << _c << "\n";
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
  output << "newdef " << name() << ":ExcitedLambda " << _elambda << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix) 
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  if(header) 
    output << "\n\" where BINARY ThePEGName=\"" 
	   << fullName() << "\";" << endl;
}
