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

int SU3BaryonDecupletOctetPhotonDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(dm.products().size()!=2){return imode;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id()),iout;
  if(id1==ParticleID::gamma){iout=id2;}
  else if(id2==ParticleID::gamma){iout=id1;}
  else{return imode;}
  unsigned int ix(0);
  cc=false;
  do
    {
      if(id0==_incomingB[ix]){if(iout==_outgoingB[ix]){imode=ix;cc=false;}}
      else if(id0==-_incomingB[ix]){if(iout==-_outgoingB[ix]){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incomingB.size()&&imode<0);
  return imode;
}

void SU3BaryonDecupletOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_C,1./GeV) << _parity << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonDecupletOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_C,1./GeV) >> _parity >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _maxweight >> iunit(_prefactor,1./GeV);
}

ClassDescription<SU3BaryonDecupletOctetPhotonDecayer> SU3BaryonDecupletOctetPhotonDecayer::initSU3BaryonDecupletOctetPhotonDecayer;
// Definition of the static class description member.

void SU3BaryonDecupletOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetPhotonDecayer> documentation
    ("The SU3BaryonDecupletOctetPhotonDecayer class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a photon.");

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,InvEnergy> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetPhotonDecayer::_C, 1.1/GeV, 1.0/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetPhotonDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "The multiplets have the same parity.",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "The multiplets have different parities.",
     false);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_deltapp, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_deltap, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_delta0, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_deltam, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmasp, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmas0, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_sigmasm, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_omega, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xis0, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::_xism, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetPhotonDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonDecupletOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
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
      A1= _prefactor[imode]*(m0-m1);B1=0.;
      A2= _prefactor[imode]*(m0+m1);B2=0.;
    }
}


// set up the decay modes
void SU3BaryonDecupletOctetPhotonDecayer::setupModes(unsigned int iopt) const
{
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_incomingB.resize(0);}
  vector<InvEnergy> factor;
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
  PDVector extpart(2);
  for(unsigned int ix=0;ix<intemp.size();++ix)
    {
      if(intemp[ix]!=0&&outtemp[ix]!=0)
	{
	  extpart[0]=getParticleData(intemp[ix]);
	  extpart[1]=getParticleData(outtemp[ix]);
	  if(extpart[0]->massMax()>extpart[1]->massMin())
	    {
	      _incomingB.push_back(intemp[ix]);
	      _outgoingB.push_back(outtemp[ix]);
	      if(iopt==1)
		{_prefactor.push_back(factor[ix]);}
	    }
	}
    }
}
void SU3BaryonDecupletOctetPhotonDecayer::dataBaseOutput(ofstream & output,
							 bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "set " << fullName() << ":Ccoupling " << _C*GeV<< "\n";
  output << "set " << fullName() << ":Parity " << _parity<< "\n";
  output << "set " << fullName() << ":Proton " << _proton << "\n";
  output << "set " << fullName() << ":Neutron " << _neutron << "\n";
  output << "set " << fullName() << ":Sigma+ " << _sigmap << "\n";
  output << "set " << fullName() << ":Sigma0 " << _sigma0 << "\n";
  output << "set " << fullName() << ":Sigma- " << _sigmam << "\n";
  output << "set " << fullName() << ":Lambda " << _lambda << "\n";
  output << "set " << fullName() << ":Xi0 " << _xi0 << "\n";
  output << "set " << fullName() << ":Xi- " << _xim << "\n";
  output << "set " << fullName() << ":Delta++ " << _deltapp << "\n";
  output << "set " << fullName() << ":Delta+ " << _deltap << "\n";
  output << "set " << fullName() << ":Delta0 " << _delta0 << "\n";
  output << "set " << fullName() << ":Delta- " << _deltam << "\n";
  output << "set " << fullName() << ":Sigma*+ " << _sigmasp << "\n";
  output << "set " << fullName() << ":Sigma*0 " << _sigmas0 << "\n";
  output << "set " << fullName() << ":Sigma*- " << _sigmasm << "\n";
  output << "set " << fullName() << ":Omega " << _omega << "\n";
  output << "set " << fullName() << ":Xi*0 " << _xis0 << "\n";
  output << "set " << fullName() << ":Xi*- " << _xism << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix)
    {output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	    << _maxweight[ix] << "\n";}
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
