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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonDecupletOctetPhotonDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

SU3BaryonDecupletOctetPhotonDecayer::~SU3BaryonDecupletOctetPhotonDecayer() {}

bool SU3BaryonDecupletOctetPhotonDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  if(_incomingB.size()==0){setupModes(1);}
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
      if(id0==_incomingB[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==ParticleID::gamma)||
	     (id2==_outgoingB[ix]&&id1==ParticleID::gamma)){allowed=true;}
	}
      else if(id0==-_incomingB[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==ParticleID::gamma)||
	     (id2==-_outgoingB[ix]&&id1==ParticleID::gamma)){allowed=true;}
	}
      ++ix;
    }
  while(ix<_incomingB.size()&&!allowed);
  return allowed;
}

ParticleVector SU3BaryonDecupletOctetPhotonDecayer::decay(const DecayMode & dm,
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
      if(id==_incomingB[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==ParticleID::gamma)||
	     (id2==_outgoingB[ix]&&id1==ParticleID::gamma)){imode=ix;cc=false;}
	}
      else if(id==-_incomingB[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==ParticleID::gamma)||
	     (id2==-_outgoingB[ix]&&id1==ParticleID::gamma)){imode=ix;cc=true;}
	}
      ++ix;
    }
  while(ix<_incomingB.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void SU3BaryonDecupletOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _maxweight << _A1 << _A2 << _A3 << _B1 << _B2 << _B3;
}

void SU3BaryonDecupletOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _maxweight >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3;
}

ClassDescription<SU3BaryonDecupletOctetPhotonDecayer> SU3BaryonDecupletOctetPhotonDecayer::initSU3BaryonDecupletOctetPhotonDecayer;
// Definition of the static class description member.

void SU3BaryonDecupletOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetPhotonDecayer> documentation
    ("The \\classname{SU3BaryonDecupletOctetPhotonDecayer} class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a photon.");

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,InvEnergy> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetPhotonDecayer::_C, 1./GeV, 1.0/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetPhotonDecayer::_parity, 0, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "The multiplets have the same parity.",
     0);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "The multiplets have different parities.",
     1);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_proton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_neutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_lambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xim, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_deltapp, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_deltap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_delta0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_deltam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmasp, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmas0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmasm, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_omega, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xis0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xism, 0, -100000, 100000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetPhotonDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}



// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonDecupletOctetPhotonDecayer::
halfThreeHalfVectorCoupling(int imode,Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const
{
  A1=_A1[imode]; A2=_A2[imode]; A3=_A3[imode];
  B1=_B1[imode]; B2=_B2[imode]; B3=_B3[imode];
}


// set up the decay modes
void SU3BaryonDecupletOctetPhotonDecayer::setupModes(unsigned int iopt) const
{
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_incomingB.resize(0);}
  vector<double> factor;
  vector<int> intemp,outtemp;
  double ortw(1./sqrt(12.)),orr(1./sqrt(3.));
  // decays of the delta+
  intemp.push_back(_deltap);outtemp.push_back(_proton);
  factor.push_back(_C*orr);
  // decays of the delta0
  intemp.push_back(_delta0);outtemp.push_back(_neutron);
  factor.push_back(_C*orr);
  // sigma*+
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);
  factor.push_back(-_C*orr);
  // sigma*0
  intemp.push_back(_sigmas0);outtemp.push_back(_lambda);
  factor.push_back(-_C*.5);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigma0);
  factor.push_back(_C*ortw);
  // xi*0
  intemp.push_back(_xis0);outtemp.push_back(_xi0);
  factor.push_back(-_C*orr);
  // set up the modes
  Energy m0,m1;
  for(unsigned int ix=0;ix<intemp.size();++ix)
    {
      if(intemp[ix]!=0&&outtemp[ix]!=0)
	{
	  _incomingB.push_back(intemp[ix]);
	  _outgoingB.push_back(outtemp[ix]);
	  if(iopt==1)
	    {
	      m0 = getParticleData(_incomingB.back())->mass();
	      m1 = getParticleData(_outgoingB.back())->mass();
	      if(_parity==0)
		{
		  _A1.push_back(0.);
		  _A2.push_back(0.);
		  _B1.push_back(-factor[ix]*(m0+m1));
		  _B2.push_back(factor[ix]*(m0+m1));
		}
	      else
		{
		  _A1.push_back(factor[ix]*(m0-m1));
		  _A2.push_back(factor[ix]*(m0+m1));
		  _B1.push_back(0.);
		  _B2.push_back(0.);
		}
	      _A3.push_back(0.);_B3.push_back(0.);
	    }
	}
    }
}
}
