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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "NonLeptonicOmegaDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

NonLeptonicOmegaDecayer::~NonLeptonicOmegaDecayer() {}

bool NonLeptonicOmegaDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(id0==_incomingB)
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){allowed=true;}
	}
      else if(id0==-_incomingB)
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){allowed=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){allowed=true;}
	}
      ++ix;
    }
  while(ix<_outgoingB.size()&&!allowed);
  return allowed;
}

ParticleVector NonLeptonicOmegaDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode=-1;
  int id=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;bool cc;
  do 
    {
      if(id==_incomingB)
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id==-_incomingB)
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){imode=ix;cc=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){imode=ix;cc=true;}
	}
      ++ix;
    }
  while(ix<_outgoingB.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void NonLeptonicOmegaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _dstar << _fstar << _omegad << _omegaf << _CBstar << _sc << _C << _fpi << _hc 
     << _hpi <<_d << _f << _Mlambda << _Mxi << _Momega << _MXistar << _mpip << _mpi0 
     << _MKp << _MK0 << _MBstar << _MR << _localmasses << _incomingB << _outgoingB 
     << _outgoingM << _A << _B << _maxweight;
}

void NonLeptonicOmegaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _dstar >> _fstar >> _omegad >> _omegaf >> _CBstar >> _sc >> _C >> _fpi >> _hc 
     >> _hpi >> _d >> _f >> _Mlambda >> _Mxi >> _Momega >> _MXistar >> _mpip >> _mpi0 
     >> _MKp >> _MK0 >> _MBstar >> _MR >> _localmasses >> _incomingB >> _outgoingB 
     >> _outgoingM >> _A >> _B >> _maxweight;
}

ClassDescription<NonLeptonicOmegaDecayer> NonLeptonicOmegaDecayer::initNonLeptonicOmegaDecayer;
// Definition of the static class description member.

void NonLeptonicOmegaDecayer::Init() {

  static ClassDocumentation<NonLeptonicOmegaDecayer> documentation
    ("The \\classname{NonLeptonicOmegaDecayer} class performs the non-leptonic decays"
     " of the omega.");

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceDStar
    ("DStar",
     "The d* coupling from hep-ph/9905398 multiplied by MB*.",
     &NonLeptonicOmegaDecayer::_dstar, 0.6, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceFStar
    ("FStar",
     "The f* coupling from hep-ph/9905398 multiplied by MB*.",
     &NonLeptonicOmegaDecayer::_fstar, -4., -10.0, 10.0,
     false, false, true);


  static Parameter<NonLeptonicOmegaDecayer,double> interfaceomegad
    ("omegad",
     "The omega_d coupling from hep-ph/9905398 multiplied by MR.",
     &NonLeptonicOmegaDecayer::_omegad, 0.82, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceomegaf
    ("omegaf",
     "The omega_f coupling from hep-ph/9905398 multiplied by MR.",
     &NonLeptonicOmegaDecayer::_omegaf,-2.17, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceCBstar
    ("CBstar",
     "The C_B* coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_CBstar, 1.35, 0.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfacesc
    ("sc",
     "The sc coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_sc,-0.85, 0.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceC
    ("C",
     "The C coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_CBstar, 1.5, 0.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &NonLeptonicOmegaDecayer::_fpi, MeV, 92.4*MeV, 0.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfacehc
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
     &NonLeptonicOmegaDecayer::_MBstar, MeV, 1660.0*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMR
    ("MR",
     "The mass of the excited R resonnaces",
     &NonLeptonicOmegaDecayer::_MR, MeV, 1700.0*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, true);

  static Switch<NonLeptonicOmegaDecayer,bool> interfaceLocalMasses
    ("LocalMasses",
     "Use local values of all the masses for the couplings.",
     &NonLeptonicOmegaDecayer::_localmasses, 0, false, false);
  static SwitchOption interfaceLocalMassesLocal
    (interfaceLocalMasses,
     "Local",
     "Use local values",
     1);
  static SwitchOption interfaceLocalMassesNonLocal
    (interfaceLocalMasses,
     "NonLocal",
     "Use values from the particle data objects",
     0);

  static ParVector<NonLeptonicOmegaDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &NonLeptonicOmegaDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void NonLeptonicOmegaDecayer::halfThreeHalfScalarCoupling(int imode,Complex& A,Complex& B) const
{A=_A[imode];B=_B[imode];}

}
