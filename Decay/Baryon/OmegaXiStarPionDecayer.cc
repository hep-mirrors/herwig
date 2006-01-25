// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaXiStarPionDecayer class.
//

#include "OmegaXiStarPionDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "OmegaXiStarPionDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

OmegaXiStarPionDecayer::~OmegaXiStarPionDecayer() {}

bool OmegaXiStarPionDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  if(id0==_idin)
    {
      if((id1==_idout&&id2==-211)||
	 (id2==_idout&&id1==-211)){allowed=true;}
    }
  else if(id0==-_idin)
    {
      if((id1==-_idout&&id2==211)||
	 (id2==-_idout&&id1==211)){allowed=true;}
    }
  return allowed;
}

ParticleVector OmegaXiStarPionDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode(0);
  bool cc(parent.id()!=_idin);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void OmegaXiStarPionDecayer::persistentOutput(PersistentOStream & os) const {
  os << _Acomm << _AP << _AS << _BP << _BS << _idin << _idout <<  _wgtmax;
}

void OmegaXiStarPionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _Acomm >> _AP >> _AS >> _BP >> _BS >> _idin >> _idout >>  _wgtmax;
}

ClassDescription<OmegaXiStarPionDecayer> OmegaXiStarPionDecayer::initOmegaXiStarPionDecayer;
// Definition of the static class description member.

void OmegaXiStarPionDecayer::Init() {

  static ClassDocumentation<OmegaXiStarPionDecayer> documentation
    ("The OmegaXiStarPionDecayer class performs the weak decay"
     " of the Omega to Xi*0 and pi-");

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAcomm
    ("Acomm",
     "The Acomm coupling for the decay",
     &OmegaXiStarPionDecayer::_Acomm, 20.91e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAP
    ("AP",
     "The A_P coupling for the decay",
     &OmegaXiStarPionDecayer::_AP, -9.20e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAS
    ("AS",
     "The A_S coupling for the decay",
     &OmegaXiStarPionDecayer::_AS, -6.32e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceBP
    ("BP",
     "The B_P coupling for the decay",
     &OmegaXiStarPionDecayer::_BP, 230.1e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceBS
    ("BS",
     "The B_S coupling for the decay",
     &OmegaXiStarPionDecayer::_BS, -100.8e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the decay",
     &OmegaXiStarPionDecayer::_wgtmax, 0.0032, 0., 100.,
     false, false, false);

  static Parameter<OmegaXiStarPionDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDF code for the incoming baryon",
     &OmegaXiStarPionDecayer::_idin, 3334, 0, 1000000,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,int> interfaceOutgoing
    ("Outgoing",
     "The PDF code for the outgoing baryon",
     &OmegaXiStarPionDecayer::_idout, 3324, 0, 1000000,
     false, false, true);
}

// couplings for spin-3/2 to spin-3/2 spin-0
void OmegaXiStarPionDecayer::
threeHalfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
				 Complex&A1,Complex&A2,Complex&B1,Complex&B2) const
{
  A2=0.;
  B2=0.;
  A1=_Acomm+_AP+_AS;
  B1=_BP+_BS;
}

  void OmegaXiStarPionDecayer::dataBaseOutput(ofstream & output,bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "set " << fullName() << ":Acomm " << _Acomm << "\n";
  output << "set " << fullName() << ":AP " << _AP << "\n";
  output << "set " << fullName() << ":AS " << _AS << "\n";
  output << "set " << fullName() << ":BP " << _BP << "\n";
  output << "set " << fullName() << ":BS " << _BS << "\n";
  output << "set " << fullName() << ":MaximumWeight " << _wgtmax << "\n";
  output << "set " << fullName() << ":Incoming " << _idin << "\n";
  output << "set " << fullName() << ":Outgoing " << _idout << "\n";
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
