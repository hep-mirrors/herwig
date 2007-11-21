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

namespace Herwig {
using namespace ThePEG;

NonLeptonicOmegaDecayer::~NonLeptonicOmegaDecayer() {}

void NonLeptonicOmegaDecayer::doinit() throw(InitException) {
  Baryon1MesonDecayerBase::doinit();
  // reset the masses if needed
  if(!_localmasses)
    {
      _Mlambda = getParticleData(3122)->mass();
      _Mxi     = getParticleData(3322)->mass();
      _Momega  = getParticleData(3334)->mass();
      _MXistar = getParticleData(3324)->mass();
      _mpip    = getParticleData(211)->mass();
      _mpi0    = getParticleData(111)->mass();
      _MKp     = getParticleData(321)->mass();
      _MK0     = getParticleData(311)->mass();
    }
  // calculate the couplings
  _A.resize(3);_B.resize(3);
  // couplings for lambda K (N.B. sign of B due to gamma_5 defn)
  _A[0] = 0.5*_C/sqrt(3.)/_fpi*((_d-3.*_f)/(_Mlambda-_Mxi)
			       +_hc/(_Momega-_MXistar))
    +_CBstar/2./sqrt(3.)/_fpi*(_dstar-3.*_fstar)/(_Mlambda/_MBstar-1.);
  _B[0] =0./MeV;
    //-0.5*_sc/sqrt(3.)/_fpi*(_omegad-3.*_omegaf)/(_Mlambda/_MR-1.);
  // couplings for xi0 pi-
  _A[1] = _C/   sqrt(2.)/_fpi*(_hc/3./(_Momega-_MXistar)
			       +_hpi*_mpip*_mpip/2./(_MKp*_MKp-_mpip*_mpip));
  _B[1] = 0./MeV;
  // couplings for xi- pi0
  _A[2] = _C/2./         _fpi*(_hc/3./(_Momega-_MXistar)
			       +_hpi*_mpi0*_mpi0/2./(_MK0*_MK0-_mpi0*_mpi0));
  _B[2] = 0./MeV;
  // set up the decay modes
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  double wgtmax;
  vector<double> wgt(0);
  for(unsigned int ix=0;ix<_outgoingB.size();++ix)
    {
      extpart[0]=getParticleData(_incomingB);
      extpart[1]=getParticleData(_outgoingB[ix]);
      extpart[2]=getParticleData(_outgoingM[ix]);
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
      if(_maxweight.size()>numberModes()){wgtmax=_maxweight[numberModes()];}
      else{wgtmax=1.;}
      addMode(mode,wgtmax,wgt);
    }
}

int NonLeptonicOmegaDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const PDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  do
    {
      if(id0==_incomingB)
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;}
	}
      else if(id0==-_incomingB)
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){imode=ix;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){imode=ix;}
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
  os << _dstar << _fstar << _omegad << _omegaf << _CBstar << _sc << _C << ounit(_fpi,GeV) 
     << ounit(_hc,GeV) << _hpi << ounit(_d,GeV) << ounit(_f,GeV) << ounit(_Mlambda,GeV) 
     << ounit(_Mxi,GeV) << ounit(_Momega,GeV) << ounit(_MXistar,GeV) << ounit(_mpip,GeV) 
     << ounit(_mpi0,GeV) << ounit(_MKp,GeV) << ounit(_MK0,GeV) << ounit(_MBstar,GeV) 
     << ounit(_MR,GeV) << _localmasses << _incomingB << _outgoingB << _outgoingM 
     << ounit(_A,1./GeV) << ounit(_B,1./GeV) << _maxweight;
}

void NonLeptonicOmegaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _dstar >> _fstar >> _omegad >> _omegaf >> _CBstar >> _sc >> _C >> iunit(_fpi,GeV) 
     >> iunit(_hc,GeV) >> _hpi >> iunit(_d,GeV) >> iunit(_f,GeV) >> iunit(_Mlambda,GeV) 
     >> iunit(_Mxi,GeV) >> iunit(_Momega,GeV) >> iunit(_MXistar,GeV) >> iunit(_mpip,GeV) 
     >> iunit(_mpi0,GeV) >> iunit(_MKp,GeV) >> iunit(_MK0,GeV) >> iunit(_MBstar,GeV) 
     >> iunit(_MR,GeV) >> _localmasses >> _incomingB >> _outgoingB >> _outgoingM 
     >> iunit(_A,1./GeV) >> iunit(_B,1./GeV) >> _maxweight;
}

ClassDescription<NonLeptonicOmegaDecayer> NonLeptonicOmegaDecayer::initNonLeptonicOmegaDecayer;
// Definition of the static class description member.

void NonLeptonicOmegaDecayer::Init() {

  static ClassDocumentation<NonLeptonicOmegaDecayer> documentation
    ("The NonLeptonicOmegaDecayer class performs the non-leptonic decays"
     " of the omega.");

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
     "The C_B* coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_CBstar, 1.35, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfacesc
    ("sc",
     "The sc coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_sc,-0.85, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceC
    ("C",
     "The C coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_C, 1.5, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &NonLeptonicOmegaDecayer::_fpi, MeV, 92.4*MeV, 0.0*MeV, 200.0*MeV,
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
     &NonLeptonicOmegaDecayer::_Mlambda, MeV, 1115.683*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMXi
    ("MXi",
     "The mass of the Xi baryon",
     &NonLeptonicOmegaDecayer::_Mxi, MeV, 1314.830*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMOmega
    ("MOmega",
     "The mass of the Omega baryon",
     &NonLeptonicOmegaDecayer::_Momega, MeV, 1672.450*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMXiStar
    ("MXiStar",
     "The mass of the XiStar baryon",
     &NonLeptonicOmegaDecayer::_MXistar, MeV, 1531.800*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMpiplus
    ("Mpiplus",
     "The mass of the charged pion",
     &NonLeptonicOmegaDecayer::_mpip, MeV, 139.57*MeV, 0.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMKplus
    ("MKplus",
     "The mass of the charged kaon",
     &NonLeptonicOmegaDecayer::_MKp, MeV, 493.667*MeV, 0.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMpi0
    ("Mpi0",
     "The mass of the neutral pion",
     &NonLeptonicOmegaDecayer::_mpi0, MeV, 134.977*MeV, 0.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMK0
    ("MK0",
     "The mass of the neutral kaon",
     &NonLeptonicOmegaDecayer::_MK0, MeV, 497.648*MeV, 0.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMBstar
    ("MBstar",
     "The mass of the excited B* resonnaces",
     &NonLeptonicOmegaDecayer::_MBstar, MeV, 1670.0*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMR
    ("MR",
     "The mass of the excited R resonnaces",
     &NonLeptonicOmegaDecayer::_MR, MeV, 1620.0*MeV, 0.0*MeV, 10000.0*MeV,
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
							  Complex& B) const
 {A=_A[imode]*(m0+m1);B=_B[imode]*(m0+m1);}

  void NonLeptonicOmegaDecayer::dataBaseOutput(ofstream & output,bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "set " << fullName() << ":DStar " << _dstar<< "\n";
  output << "set " << fullName() << ":FStar " << _fstar << "\n";
  output << "set " << fullName() << ":omegad " << _omegad<< "\n";
  output << "set " << fullName() << ":omegaf " << _omegaf<< "\n";
  output << "set " << fullName() << ":CBstar " << _CBstar<< "\n";
  output << "set " << fullName() << ":sc " << _sc << "\n";
  output << "set " << fullName() << ":C " << _C << "\n";
  output << "set " << fullName() << ":Fpi " << _fpi/MeV << "\n";
  output << "set " << fullName() << ":hc " << _hc/GeV << "\n";
  output << "set " << fullName() << ":hpi " << _hpi<< "\n";
  output << "set " << fullName() << ":d " << _d/GeV << "\n";
  output << "set " << fullName() << ":f " << _f/GeV << "\n";
  output << "set " << fullName() << ":MLambda " << _Mlambda/MeV << "\n";
  output << "set " << fullName() << ":MXi " << _Mxi/MeV << "\n";
  output << "set " << fullName() << ":MOmega " << _Momega/MeV << "\n";
  output << "set " << fullName() << ":MXiStar " << _MXistar/MeV << "\n";
  output << "set " << fullName() << ":Mpiplus " << _mpip/MeV << "\n";
  output << "set " << fullName() << ":MKplus " << _MKp/MeV << "\n";
  output << "set " << fullName() << ":Mpi0 " << _mpi0/MeV << "\n";
  output << "set " << fullName() << ":MK0 " << _MK0/MeV << "\n";
  output << "set " << fullName() << ":MBstar " << _MBstar/MeV << "\n";
  output << "set " << fullName() << ":MR " << _MR/MeV << "\n";
  output << "set " << fullName() << ":LocalMasses " << _localmasses << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix)
    {output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	    << _maxweight[ix] << "\n";}
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
