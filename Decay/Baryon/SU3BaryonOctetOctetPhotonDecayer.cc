// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3OctetOctetPhotonDecayer class.
//

#include "SU3BaryonOctetOctetPhotonDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonOctetOctetPhotonDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

SU3BaryonOctetOctetPhotonDecayer::~SU3BaryonOctetOctetPhotonDecayer() {}

bool SU3BaryonOctetOctetPhotonDecayer::accept(const DecayMode & dm) const {
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

ParticleVector SU3BaryonOctetOctetPhotonDecayer::decay(const DecayMode & dm,
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


void SU3BaryonOctetOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _lf << _ld <<  _parity << _proton << _neutron << _sigma0 << _sigmap << _sigmam 
     << _lambda << _xi0 << _xim << _eproton << _eneutron << _esigma0 << _esigmap 
     << _esigmam << _elambda << _exi0 << _exim << _incomingB << _outgoingB << _maxweight
     << _A1 << _A2 << _A3 << _B1 << _B2 << _B3;
}

void SU3BaryonOctetOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _lf >> _ld >>  _parity >> _proton >> _neutron >> _sigma0 >> _sigmap >> _sigmam 
     >> _lambda >> _xi0 >> _xim >> _eproton >> _eneutron >> _esigma0 >> _esigmap 
     >> _esigmam >> _elambda >> _exi0 >> _exim >> _incomingB >> _outgoingB >> _maxweight
     >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3;
}

ClassDescription<SU3BaryonOctetOctetPhotonDecayer> SU3BaryonOctetOctetPhotonDecayer::initSU3BaryonOctetOctetPhotonDecayer;
// Definition of the static class description member.

void SU3BaryonOctetOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetOctetPhotonDecayer> documentation
    ("The \\classname{SU3BaryonOctetOctetPhotonDecayer} class is designed for the "
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

  static Switch<SU3BaryonOctetOctetPhotonDecayer,int> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetOctetPhotonDecayer::_parity, 0, false, false);
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

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_proton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_neutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_lambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_xi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_xim, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedProton
    ("ExcitedProton",
     "The PDG code for the heavier proton-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_eproton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedNeutron
    ("ExcitedNeutron",
     "The PDG code for the heavier neutron-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_eneutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigmap
    ("ExcitedSigma+",
     "The PDG code for the heavier Sigma+-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigma0
    ("ExcitedSigma0",
     "The PDG code for the heavier Sigma0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigmam
    ("ExcitedSigma-",
     "The PDG code for the heavier Sigma--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_elambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedXi0
    ("ExcitedXi0",
     "The PDG code for the heavier Xi0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_exi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedXim
    ("ExcitedXi-",
     "The PDG code for the heavier Xi--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_exim, 0, -100000, 100000,
     false, false, true);

  static ParVector<SU3BaryonOctetOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetOctetPhotonDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-1
void SU3BaryonOctetOctetPhotonDecayer::halfHalfVectorCoupling(int imode,Complex&A1,
							      Complex&A2,Complex&B1,
							      Complex&B2) const
 {A1=_A1[imode];A2=_A2[imode];B1=_B1[imode];B2=_B2[imode];}


// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonOctetOctetPhotonDecayer::
halfThreeHalfVectorCoupling(int imode,Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const
{
  A1=_A1[imode]; A2=_A2[imode]; A3=_A3[imode];
  B1=_B1[imode]; B2=_B2[imode]; B3=_B3[imode];
}

// set up the decay modes
void SU3BaryonOctetOctetPhotonDecayer::setupModes(unsigned int iopt) const
{
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_incomingB.resize(0);}
  // set up for the various different decay modes
  vector<double> factor;
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
  Energy m0,m1;
  int inspin,outspin;
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
	      inspin  = getParticleData(_incomingB.back())->iSpin();
	      outspin = getParticleData(_outgoingB.back())->iSpin();
	      /*
		if(m0>m1)
		{
		cout << "testing the width " 
		<< _incomingB.back() << "   " << _outgoingB.back() << "   " 
		<< m0 << "  " << m1 << "    " 
		<< 1./8./pi/m0/m0/2./m0*(m0*m0-m1*m1)*
		128.*factor[ix]*factor[ix]/4.*(m0*m0-m1*m1)*(m0*m0-m1*m1) << endl;
		}
	      */
	      factor[ix] *=4.;
	      if(inspin==2&&outspin==2)
		{
		  if(_parity==0)
		    {
		      _A1.push_back(0.);
		      _B1.push_back(factor[ix]*(m1-m0));
		      _A2.push_back(0.);
		      _B2.push_back(-2.*factor[ix]*(m0+m1));
		    }
		  else
		    {
		      _A1.push_back(factor[ix]*(m0+m1));
		      _B1.push_back(0.);
		      _A2.push_back(-2.*factor[ix]*(m0+m1));
		      _B2.push_back(0.);
		    }
		  _A3.push_back(0.);_B3.push_back(0.);
		}
	      else if(inspin==4&&outspin==2)
		{
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
	      else
		{throw DecayIntegratorError() << "Invalid combination of spins in "
					      << "SU3BaryonOctetOctetPhotonDecayer::" 
					      << "setupModes()" << Exception::abortnow;}
	    }
	}
    }
}
}
