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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonSingletOctetScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

SU3BaryonSingletOctetScalarDecayer::~SU3BaryonSingletOctetScalarDecayer() {}

bool SU3BaryonSingletOctetScalarDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  if(_outgoingB.size()==0){setupModes(1);}
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
      if(id0==_elambda)
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){allowed=true;}
	}
      else if(id0==-_elambda)
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
  return false;
}

ParticleVector SU3BaryonSingletOctetScalarDecayer::decay(const DecayMode & dm,
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
      if(id==_elambda)
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id==-_elambda)
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


void SU3BaryonSingletOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _fpi << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _elambda << _outgoingB 
     << _outgoingM << _maxweight << _A1 << _A2 << _A3 << _B1 << _B2 << _B3;
}

void SU3BaryonSingletOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _fpi >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _elambda >> _outgoingB 
     >> _outgoingM >> _maxweight >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3;
}

ClassDescription<SU3BaryonSingletOctetScalarDecayer> SU3BaryonSingletOctetScalarDecayer::initSU3BaryonSingletOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonSingletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonSingletOctetScalarDecayer> documentation
    ("The \\classname{SU3BaryonSingletOctetScalarDecayer} class is designed for the decay"
     " of an excited SU(3) singlet baryon");

  static Parameter<SU3BaryonSingletOctetScalarDecayer,double> interfaceCcoupling
    ("Coupling",
     "The C coupling of the baryon resonances",
     &SU3BaryonSingletOctetScalarDecayer::_C, 0.0, -10.0, 10.0,
     false, false, true);
  static Switch<SU3BaryonSingletOctetScalarDecayer,int> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonSingletOctetScalarDecayer::_parity, 0, false, false);
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

  static Parameter<SU3BaryonSingletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonSingletOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_proton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_neutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_lambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xim, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_elambda, 0, -100000, 100000,
     false, false, true);

  static ParVector<SU3BaryonSingletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonSingletOctetScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void SU3BaryonSingletOctetScalarDecayer::halfHalfScalarCoupling(int imode,Complex& A,Complex& B) const
 {A=_A1[imode];B=_B1[imode];}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonSingletOctetScalarDecayer::halfThreeHalfScalarCoupling(int imode,Complex& A,Complex& B) const
{A=_A1[imode];B=_B1[imode];}

// set up the decay modes
void SU3BaryonSingletOctetScalarDecayer::setupModes(unsigned int iopt) const
{
  if(_outgoingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_outgoingM.resize(0);}
  // set up for the various different decay modes
  vector<int> outtemp,mestemp;
  double rt(sqrt(2.));
  if(_elambda==0){throw DecayIntegratorError() << "Invalid incoming particle in "
					       << "SU3BaryonSingletOctetScalarDecayer::" 
					       << "setupModes()" << Exception::abortnow;}
  // decays of the excited lambda
  outtemp.push_back(_sigma0);mestemp.push_back(111);
  outtemp.push_back(_sigmap);mestemp.push_back(-211);
  outtemp.push_back(_sigmam);mestemp.push_back(211);
  outtemp.push_back(_lambda);mestemp.push_back(221);
  outtemp.push_back(_xim);mestemp.push_back(321);
  outtemp.push_back(_xi0);mestemp.push_back(311);
  outtemp.push_back(_proton);mestemp.push_back(-321);
  outtemp.push_back(_neutron);mestemp.push_back(-311);
  Energy m0(getParticleData(_elambda)->mass()),m1;
  int inspin(getParticleData(_elambda)->iSpin()),outspin;
  for(unsigned int ix=0;ix<outtemp.size();++ix)
    {
      if(outtemp[ix]!=0&&mestemp[ix]!=0)
	{
	  _outgoingB.push_back(outtemp[ix]);
	  _outgoingM.push_back(mestemp[ix]);
	  if(iopt==1)
	    {
	      m1 = getParticleData(_outgoingB.back())->mass();
	      outspin = getParticleData(_outgoingB.back())->iSpin();
	      if(inspin==2&&outspin==2)
		{
		  if(_parity==0)
		    {
		      _A1.push_back(0.);
		      _B1.push_back(_C*rt/_fpi*(m0+m1));
		    }
		  else
		    {
		      _A1.push_back(_C*rt/_fpi*(m0-m1));
		      _B1.push_back(0.);
		    }
		  _A2.push_back(0.);_B2.push_back(0.);
		  _A3.push_back(0.);_B3.push_back(0.);
		}
	      else if(inspin==4&&outspin==2)
		{
		  if(_parity==0)
		    {
		      _A1.push_back(_C*rt/_fpi*(m0+m1));
		      _B1.push_back(0.);
		    }
		  else
		    {
		      _A1.push_back(0.);
		      _B1.push_back(_C*rt/_fpi*(m0+m1));
		    }
		  _A2.push_back(0.);_B2.push_back(0.);
		  _A3.push_back(0.);_B3.push_back(0.);
		}
	    }
	  else
	    {throw DecayIntegratorError() << "Invalid combination of spins in "
					   << "SU3BaryonSingletOctetScalarDecayer::" 
					   << "setupModes()" << Exception::abortnow;}
	}
    }
}

}
