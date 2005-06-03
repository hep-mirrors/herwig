// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonSingletOctetPhotonDecayer class.
//

#include "SU3BaryonSingletOctetPhotonDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonSingletOctetPhotonDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

SU3BaryonSingletOctetPhotonDecayer::~SU3BaryonSingletOctetPhotonDecayer() 
{}

bool SU3BaryonSingletOctetPhotonDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  if(_outgoingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id()),iout;
  if(id1==ParticleID::gamma){iout=id2;}
  else if(id2==ParticleID::gamma){iout=id1;}
  else{return false;}
  unsigned int ix(0);
  do
    {
      if(id0==_elambda){if(iout==_outgoingB[ix]){allowed=true;}}
      else if(id0==-_elambda){if(iout==-_outgoingB[ix]){allowed=true;}}
      ++ix;
    }
  while(ix<_outgoingB.size()&&!allowed);
  return allowed;
}

ParticleVector SU3BaryonSingletOctetPhotonDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode(-1),id(parent.id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id()),iout;
  unsigned int ix(0);bool cc;
  if(id1==ParticleID::gamma){iout=id2;}
  else{iout=id1;}
  do 
    {
      if(id==_elambda){if(iout==_outgoingB[ix]){imode=ix;cc=false;}}
      else if(id==-_elambda){if(iout==-_outgoingB[ix]){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_outgoingB.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void SU3BaryonSingletOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _sigma0 << _lambda << _elambda << _outgoingB 
     << _maxweight << _prefactor;
}

void SU3BaryonSingletOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _sigma0 >> _lambda >> _elambda >> _outgoingB 
     >> _maxweight >> _prefactor;
}

ClassDescription<SU3BaryonSingletOctetPhotonDecayer> SU3BaryonSingletOctetPhotonDecayer::initSU3BaryonSingletOctetPhotonDecayer;
// Definition of the static class description member.

void SU3BaryonSingletOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonSingletOctetPhotonDecayer> documentation
    ("The \\classname{SU3BaryonSingletOctetPhotonDecayer} class performs the decay"
     " of a singlet baryon to an octet baryon and a photon.");

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The C coupling of the baryon resonances.",
     &SU3BaryonSingletOctetPhotonDecayer::_C, 1./GeV, 0.252/GeV, -10./GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonSingletOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonSingletOctetPhotonDecayer::_parity, true, false, false);
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

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_elambda, 3124, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonSingletOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonSingletOctetPhotonDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-1
void SU3BaryonSingletOctetPhotonDecayer::
halfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
		       Complex&A1,Complex&A2,Complex&B1,Complex&B2) const
{
  if(_parity)
    {
      A1=    _prefactor[imode]*(m0+m1);B1=0.;
      A2=-2.*_prefactor[imode]*(m0+m1);B2=0.;
    }
  else
    {
      A1=0.;B1=    _prefactor[imode]*(m1-m0);
      A2=0.;B2=-2.*_prefactor[imode]*(m0+m1);
    }
}

// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonSingletOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy m2,
			    Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const
{
  A3=0.;B3=0.;
  if(_parity)
    {
      A1=0.;B1=-_prefactor[imode]*(m0+m1);
      A2=0.;B2= _prefactor[imode]*(m0+m1);
    }
  else
    {
      A1=_prefactor[imode]*(m0-m1);B1=0.;
      A2=_prefactor[imode]*(m0+m1);B2=0.;
    }
}

// set up the decay modes
void SU3BaryonSingletOctetPhotonDecayer::setupModes(unsigned int iopt) const
{
  if(_outgoingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);}
  // set up for the various different decay modes
  vector<int> outtemp;
  vector<double> factor;
  if(_elambda==0){throw DecayIntegratorError() << "Invalid incoming particle in "
					       << "SU3BaryonSingletOctetScalarDecayer::" 
					       << "setupModes()" << Exception::abortnow;}
  // decays of the excited lambda
  outtemp.push_back(_sigma0);factor.push_back(_C/sqrt(2.));
  outtemp.push_back(_lambda);factor.push_back(_C/sqrt(6.));
  PDVector extpart(2);extpart[0]=getParticleData(_elambda);
  int inspin(extpart[0]->iSpin()),outspin;
  for(unsigned int ix=0;ix<outtemp.size();++ix)
    {
      if(outtemp[ix]!=0)
	{
	  extpart[1]=getParticleData(outtemp[ix]);
	  if(extpart[0]->massMax()>extpart[1]->massMin())
	    {
	      _outgoingB.push_back(outtemp[ix]);
	      if(iopt==1)
		{
		  outspin = extpart[1]->iSpin();
		  factor[ix] *=2.;
		  if(inspin==2&&outspin==2)
		    {_prefactor.push_back(2.*factor[ix]);}
		  else if(inspin==4&&outspin==2)
		    {_prefactor.push_back(factor[ix]);}
		  else
		    {throw DecayIntegratorError() 
			<< "Invalid combination of spins in "
			<< "SU3BaryonSingletOctetScalarDecayer::" 
			<< "setupModes()" << Exception::abortnow;}
		}
	    }
	}
    }
}

void SU3BaryonSingletOctetPhotonDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":Points " << _npoint << "\n";
  output << "set " << fullName() << ":Coupling " << _C*GeV << "\n";
  output << "set " << fullName() << ":Parity " << _parity<< "\n";
  output << "set " << fullName() << ":Sigma0 " << _sigma0 << "\n";
  output << "set " << fullName() << ":Lambda " << _lambda << "\n";
  output << "set " << fullName() << ":ExcitedLambda " << _elambda << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix)
    {output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	    << _maxweight[ix] << "\n";}
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
