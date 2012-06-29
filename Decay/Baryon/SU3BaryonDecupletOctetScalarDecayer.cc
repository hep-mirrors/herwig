// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonDecupletOctetScalarDecayer class.
//

#include "SU3BaryonDecupletOctetScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonDecupletOctetScalarDecayer::SU3BaryonDecupletOctetScalarDecayer() {
  // couplings and off-shell parameter
  _c=1.5;
  // the relative parities of the two baryon multiplets
  _parity=true;
  // the pion decay constant
  _fpi=92.4*MeV;
  // PDG codes for the various octet baryons
  _proton   = 2212;
  _neutron  = 2112;
  _sigma0   = 3212;
  _sigmap   = 3222;
  _sigmam   = 3112;
  _lambda   = 3122;
  _xi0      = 3322;
  _xim      = 3312;
  // PDG codes for the various decuplet baryons
  _deltapp  = 2224;
  _deltap   = 2214;
  _delta0   = 2114;
  _deltam   = 1114;
  _sigmasp  = 3224;
  _sigmas0  = 3214;
  _sigmasm  = 3114;
  _omega    = 3334;
  _xism     = 3314;
  _xis0     = 3324;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonDecupletOctetScalarDecayer::doinit() {
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
    wgtmax = _maxweight.size()>numberModes()? 
      _maxweight[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
  }
}

void SU3BaryonDecupletOctetScalarDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonDecupletOctetScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
  cc=false;
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

void SU3BaryonDecupletOctetScalarDecayer::
persistentOutput(PersistentOStream & os) const {
  os << _c << _parity << ounit(_fpi,GeV) << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _outgoingM << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonDecupletOctetScalarDecayer::
persistentInput(PersistentIStream & is, int) {
  is >> _c >> _parity >> iunit(_fpi,GeV) >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _outgoingM >> _maxweight >> iunit(_prefactor,1./GeV);
}

ClassDescription<SU3BaryonDecupletOctetScalarDecayer> 
SU3BaryonDecupletOctetScalarDecayer::initSU3BaryonDecupletOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonDecupletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetScalarDecayer> documentation
    ("The SU3BaryonDecupletOctetScalarDecayer class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a pseudoscalar"
     " meson from the lightest multiplet.");

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,double> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetScalarDecayer::_c, 1.5, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetScalarDecayer::_parity, true, false, false);
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

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonDecupletOctetScalarDecayer::_fpi, MeV, 92.4*MeV, 50.0*MeV, 150.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltapp, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltap, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_delta0, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltam, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasp, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmas0, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasm, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_omega, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xis0, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xism, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonDecupletOctetScalarDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex& A, Complex& B) const {
  if(_parity) {
    A = _prefactor[imode]*(m0+m1);
    B = 0.;
  }
  else {
    A = 0.;
    B = _prefactor[imode]*(m0+m1);
  }
}

// set up the decay modes
void SU3BaryonDecupletOctetScalarDecayer::setupModes(unsigned int iopt) const {
  if(_incomingB.size()!=0&&iopt==0) return;
  if(iopt==1) {
    _outgoingB.clear();
    _incomingB.clear();
    _outgoingM.clear();
  }
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.)),orr(1./sqrt(3.));
  // decays of the delta++
  intemp.push_back(_deltapp);outtemp.push_back(_proton);mestemp.push_back(211);
  factor.push_back(_c);
  intemp.push_back(_deltapp);outtemp.push_back(_sigmap);mestemp.push_back(321);
  factor.push_back(-_c);
  // decays of the delta+
  intemp.push_back(_deltap);outtemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back(_c*orr);
  intemp.push_back(_deltap);outtemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigmap);mestemp.push_back(311);
  factor.push_back(_c*orr);
  // decays of the delta0
  intemp.push_back(_delta0);outtemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back(_c*orr);
  intemp.push_back(_delta0);outtemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigmam);mestemp.push_back(321);
  factor.push_back(_c*orr);
  // decays of the delta-
  intemp.push_back(_deltam);outtemp.push_back(_neutron);mestemp.push_back(-211);
  factor.push_back(-_c);
  intemp.push_back(_deltam);outtemp.push_back(_sigmam);mestemp.push_back(311);
  factor.push_back(_c);
  // sigma*+
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_proton);mestemp.push_back(-311);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_xi0);mestemp.push_back(321);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(_c*ort);
  intemp.push_back(_sigmasp);outtemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(_c*ort);
  // sigma*0
  intemp.push_back(_sigmas0);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(_c*ort);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigma0);mestemp.push_back(221);
  factor.push_back(_c*ort);
  // sigma*-
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_neutron);mestemp.push_back(-321);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_xim);mestemp.push_back(311);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(_c*ort);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(_c*ort);
  // xi*0
  intemp.push_back(_xis0);outtemp.push_back(_xim);mestemp.push_back(211);
  factor.push_back(_c*orr);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_xis0);outtemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(_c*orr);
  intemp.push_back(_xis0);outtemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(_c*ors);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(_c*ort);
  intemp.push_back(_xis0);outtemp.push_back(_lambda);mestemp.push_back(-311);
  factor.push_back(_c*ort);
  // xi*-
  intemp.push_back(_xism);outtemp.push_back(_xi0);mestemp.push_back(-211);
  factor.push_back(_c*orr);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_xism);outtemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_c*orr);
  intemp.push_back(_xism);outtemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(_c*ors);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(_c*ort);
  intemp.push_back(_xism);outtemp.push_back(_lambda);mestemp.push_back(-321);
  factor.push_back(_c*ort);
  // set up the modes
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
	if(iopt==1) _prefactor.push_back(factor[ix]/_fpi);
      }
    }
  }
}

void SU3BaryonDecupletOctetScalarDecayer::dataBaseOutput(ofstream & output, 
							 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Ccoupling " << _c << "\n";
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
  output << "newdef " << name() << ":Delta++ " << _deltapp << "\n";
  output << "newdef " << name() << ":Delta+ " << _deltap << "\n";
  output << "newdef " << name() << ":Delta0 " << _delta0 << "\n";
  output << "newdef " << name() << ":Delta- " << _deltam << "\n";
  output << "newdef " << name() << ":Sigma*+ " << _sigmasp << "\n";
  output << "newdef " << name() << ":Sigma*0 " << _sigmas0 << "\n";
  output << "newdef " << name() << ":Sigma*- " << _sigmasm << "\n";
  output << "newdef " << name() << ":Omega " << _omega << "\n";
  output << "newdef " << name() << ":Xi*0 " << _xis0 << "\n";
  output << "newdef " << name() << ":Xi*- " << _xism << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix)
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
