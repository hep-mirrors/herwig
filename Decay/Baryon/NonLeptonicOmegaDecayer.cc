// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NonLeptonicOmegaDecayer class.
//

#include "NonLeptonicOmegaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

NonLeptonicOmegaDecayer::NonLeptonicOmegaDecayer() {
  // couplings for the decays
  _fstar  =  0.5e-7;
  _dstar  = -3.98e-7;
  _omegad = -1.16e-8/1.8;
  _omegaf =  1.53e-8/1.8;
  _cbstar =  1.35;
  _sc     = -0.85;
  _c      =  1.50;
  _fpi    =  92.4*MeV;
  _hc     =  0.39e-7*GeV;
  _hpi    =  3.2e-7;
  _d      =  0.44e-7*GeV;
  _f      = -0.50e-7*GeV;
  // massses of the particles
  _mlambda = 1115.683*MeV;
  _mxi     = 1314.830*MeV;
  _momega  = 1672.450*MeV;
  _mxistar = 1531.800*MeV;
  _mpip    =  139.570*MeV;
  _mpi0    =  134.977*MeV;
  _mkp     =  493.667*MeV;
  _mk0     =  497.648*MeV;
  _mbstar  =  1620   *MeV;
  _mr      =  1500   *MeV;
  // use local values for the masses for the couplings
  _localmasses=true;
  // the PDG codes for the modes
  _incomingB = 3334;
  _outgoingB.resize(3);_outgoingM.resize(3);_maxweight.resize(3);
  _outgoingB[0] = 3122;_outgoingM[0] =-321;_maxweight[0]=1.5;
  _outgoingB[1] = 3322;_outgoingM[1] =-211;_maxweight[1]=0.4;
  _outgoingB[2] = 3312;_outgoingM[2] = 111;_maxweight[2]=0.2;
  // intermediates
  generateIntermediates(false);
}

void NonLeptonicOmegaDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

void NonLeptonicOmegaDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // reset the masses if needed
  if(!_localmasses) {
    _mlambda = getParticleData(3122)->mass();
    _mxi     = getParticleData(3322)->mass();
    _momega  = getParticleData(3334)->mass();
    _mxistar = getParticleData(3324)->mass();
    _mpip    = getParticleData(211)->mass();
    _mpi0    = getParticleData(111)->mass();
    _mkp     = getParticleData(321)->mass();
    _mk0     = getParticleData(311)->mass();
  }
  // calculate the couplings
  _a.resize(3);_b.resize(3);
  // couplings for lambda K (N.B. sign of B due to gamma_5 defn)
  _a[0] =  0.5*_c/sqrt(3.)/_fpi*((_d-3.*_f)/(_mlambda-_mxi)
			       +_hc/(_momega-_mxistar))
    +0.5*_cbstar/sqrt(3.)/_fpi*(_dstar-3.*_fstar)/(_mlambda/_mbstar-1.);
  _b[0] =  0.5*_sc/sqrt(3.)/_fpi*(_omegad-3.*_omegaf)/(_mlambda/_mr-1.);
  // couplings for xi0 pi-
  _a[1] = _c/   sqrt(2.)/_fpi*(_hc/3./(_momega-_mxistar)
			       +_hpi*_mpip*_mpip/2./(_mkp*_mkp-_mpip*_mpip));
  _b[1] = ZERO;
  // couplings for xi- pi0
  _a[2] = _c/2./         _fpi*(_hc/3./(_momega-_mxistar)
			       +_hpi*_mpi0*_mpi0/2./(_mk0*_mk0-_mpi0*_mpi0));
  _b[2] = ZERO;
  // set up the decay modes
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  double wgtmax;
  vector<double> wgt(0);
  for(unsigned int ix=0;ix<_outgoingB.size();++ix) {
    extpart[0]=getParticleData(_incomingB);
    extpart[1]=getParticleData(_outgoingB[ix]);
    extpart[2]=getParticleData(_outgoingM[ix]);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    wgtmax = _maxweight.size()>numberModes() ? _maxweight[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
  }
}

int NonLeptonicOmegaDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  do {
    if(id0==_incomingB) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) imode=ix;
    }
    else if(id0==-_incomingB) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) imode=ix;
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) imode=ix;
    }
    ++ix;
  }
  while(ix<_outgoingB.size()&&imode<0);
  // charge conjugation
  cc=id0<0;
  // return the answer
  return imode;
}

void NonLeptonicOmegaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _dstar << _fstar << _omegad << _omegaf << _cbstar << _sc << _c 
     << ounit(_fpi,GeV) 
     << ounit(_hc,GeV) << _hpi << ounit(_d,GeV) << ounit(_f,GeV) << ounit(_mlambda,GeV) 
     << ounit(_mxi,GeV) << ounit(_momega,GeV) << ounit(_mxistar,GeV) << ounit(_mpip,GeV) 
     << ounit(_mpi0,GeV) << ounit(_mkp,GeV) << ounit(_mk0,GeV) << ounit(_mbstar,GeV) 
     << ounit(_mr,GeV) << _localmasses << _incomingB << _outgoingB << _outgoingM 
     << ounit(_a,1./GeV) << ounit(_b,1./GeV) << _maxweight;
}

void NonLeptonicOmegaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _dstar >> _fstar >> _omegad >> _omegaf >> _cbstar >> _sc >> _c 
     >> iunit(_fpi,GeV) 
     >> iunit(_hc,GeV) >> _hpi >> iunit(_d,GeV) >> iunit(_f,GeV) >> iunit(_mlambda,GeV) 
     >> iunit(_mxi,GeV) >> iunit(_momega,GeV) >> iunit(_mxistar,GeV) >> iunit(_mpip,GeV) 
     >> iunit(_mpi0,GeV) >> iunit(_mkp,GeV) >> iunit(_mk0,GeV) >> iunit(_mbstar,GeV) 
     >> iunit(_mr,GeV) >> _localmasses >> _incomingB >> _outgoingB >> _outgoingM 
     >> iunit(_a,1./GeV) >> iunit(_b,1./GeV) >> _maxweight;
}

ClassDescription<NonLeptonicOmegaDecayer> 
NonLeptonicOmegaDecayer::initNonLeptonicOmegaDecayer;
// Definition of the static class description member.

void NonLeptonicOmegaDecayer::Init() {

  static ClassDocumentation<NonLeptonicOmegaDecayer> documentation
    ("The NonLeptonicOmegaDecayer class performs the non-leptonic decays"
     " of the omega based on the results of hep-ph/9905398.",
     "The  non-leptonic decays of the Omega baryon were simulated "
     "using the NonLeptonicOmegaDecayer class based on the results of"
     "\\cite{Borasoy:1999ip}.",
     "\\bibitem{Borasoy:1999ip}\n"
     "B.~Borasoy and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 60} (1999) 054021 [arXiv:hep-ph/9905398].\n"
     "%%CITATION = PHRVA,D60,054021;%%\n");

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceDStar
    ("DStar",
     "The d* coupling from hep-ph/9905398 multiplied by MB*.",
     &NonLeptonicOmegaDecayer::_dstar, -3.98e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceFStar
    ("FStar",
     "The f* coupling from hep-ph/9905398 multiplied by MB*.",
     &NonLeptonicOmegaDecayer::_fstar, 0.5e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceomegad
    ("omegad",
     "The omega_d coupling from hep-ph/9905398 multiplied by MR.",
     &NonLeptonicOmegaDecayer::_omegad, -1.16e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceomegaf
    ("omegaf",
     "The omega_f coupling from hep-ph/9905398 multiplied by MR.",
     &NonLeptonicOmegaDecayer::_omegaf, 1.53e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceCBstar
    ("CBstar",
     "The C_b* coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_cbstar, 1.35, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfacesc
    ("sc",
     "The sc coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_sc,-0.85, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceC
    ("C",
     "The C coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_c, 1.5, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &NonLeptonicOmegaDecayer::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfacehc
    ("hc",
     "The h_c coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_hc, GeV, 0.39e-7*GeV, -10.0e-7*GeV, 10.0e-7*GeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfacehpi
    ("hpi",
     "The h_pi coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_hpi, 3.2e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaced
    ("d",
     "The d coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_d, GeV, 0.44e-7*GeV, -1e-6*GeV, 1e-6*GeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfacef
    ("f",
     "The f coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_f, GeV, -0.50e-7*GeV, -1e-6*GeV, 1e-6*GeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMLambda
    ("MLambda",
     "The mass of the Lambda baryon",
     &NonLeptonicOmegaDecayer::_mlambda, MeV, 1115.683*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMXi
    ("MXi",
     "The mass of the Xi baryon",
     &NonLeptonicOmegaDecayer::_mxi, MeV, 1314.830*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMOmega
    ("MOmega",
     "The mass of the Omega baryon",
     &NonLeptonicOmegaDecayer::_momega, MeV, 1672.450*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMXiStar
    ("MXiStar",
     "The mass of the XiStar baryon",
     &NonLeptonicOmegaDecayer::_mxistar, MeV, 1531.800*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMpiplus
    ("Mpiplus",
     "The mass of the charged pion",
     &NonLeptonicOmegaDecayer::_mpip, MeV, 139.57*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMKplus
    ("MKplus",
     "The mass of the charged kaon",
     &NonLeptonicOmegaDecayer::_mkp, MeV, 493.667*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMpi0
    ("Mpi0",
     "The mass of the neutral pion",
     &NonLeptonicOmegaDecayer::_mpi0, MeV, 134.977*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMK0
    ("MK0",
     "The mass of the neutral kaon",
     &NonLeptonicOmegaDecayer::_mk0, MeV, 497.648*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMBstar
    ("MBstar",
     "The mass of the excited B* resonnaces",
     &NonLeptonicOmegaDecayer::_mbstar, MeV, 1620.0*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMR
    ("MR",
     "The mass of the excited R resonnaces",
     &NonLeptonicOmegaDecayer::_mr, MeV, 1620.0*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Switch<NonLeptonicOmegaDecayer,bool> interfaceLocalMasses
    ("LocalMasses",
     "Use local values of all the masses for the couplings.",
     &NonLeptonicOmegaDecayer::_localmasses, true, false, false);
  static SwitchOption interfaceLocalMassesLocal
    (interfaceLocalMasses,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalMassesNonLocal
    (interfaceLocalMasses,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static ParVector<NonLeptonicOmegaDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &NonLeptonicOmegaDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void NonLeptonicOmegaDecayer::threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,
							  Energy,Complex& A,
							  Complex& B) const {
  useMe();
  A = _a[imode]*(m0+m1);
  B = _b[imode]*(m0+m1);
}

void NonLeptonicOmegaDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  output << "newdef " << name() << ":DStar " << _dstar<< "\n";
  output << "newdef " << name() << ":FStar " << _fstar << "\n";
  output << "newdef " << name() << ":omegad " << _omegad<< "\n";
  output << "newdef " << name() << ":omegaf " << _omegaf<< "\n";
  output << "newdef " << name() << ":CBstar " << _cbstar<< "\n";
  output << "newdef " << name() << ":sc " << _sc << "\n";
  output << "newdef " << name() << ":C " << _c << "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":hc " << _hc/GeV << "\n";
  output << "newdef " << name() << ":hpi " << _hpi<< "\n";
  output << "newdef " << name() << ":d " << _d/GeV << "\n";
  output << "newdef " << name() << ":f " << _f/GeV << "\n";
  output << "newdef " << name() << ":MLambda " << _mlambda/MeV << "\n";
  output << "newdef " << name() << ":MXi " << _mxi/MeV << "\n";
  output << "newdef " << name() << ":MOmega " << _momega/MeV << "\n";
  output << "newdef " << name() << ":MXiStar " << _mxistar/MeV << "\n";
  output << "newdef " << name() << ":Mpiplus " << _mpip/MeV << "\n";
  output << "newdef " << name() << ":MKplus " << _mkp/MeV << "\n";
  output << "newdef " << name() << ":Mpi0 " << _mpi0/MeV << "\n";
  output << "newdef " << name() << ":MK0 " << _mk0/MeV << "\n";
  output << "newdef " << name() << ":MBstar " << _mbstar/MeV << "\n";
  output << "newdef " << name() << ":MR " << _mr/MeV << "\n";
  output << "newdef " << name() << ":LocalMasses " << _localmasses << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}
