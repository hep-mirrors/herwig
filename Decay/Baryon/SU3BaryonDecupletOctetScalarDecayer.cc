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

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonDecupletOctetScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig{
using namespace ThePEG;

SU3BaryonDecupletOctetScalarDecayer::~SU3BaryonDecupletOctetScalarDecayer() {}

bool SU3BaryonDecupletOctetScalarDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  if(_incomingB.size()==0){setupModes(0);}
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
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){allowed=true;}
	}
      else if(id0==-_incomingB[ix])
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
  while(ix<_incomingB.size()&&!allowed);
  return allowed;
}

ParticleVector SU3BaryonDecupletOctetScalarDecayer::decay(const DecayMode & dm,
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
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id==-_incomingB[ix])
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
  while(ix<_incomingB.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void SU3BaryonDecupletOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _fpi << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _outgoingM << _maxweight << _A1 << _A2 << _A3 << _B1 << _B2 << _B3;
}

void SU3BaryonDecupletOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _fpi >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _outgoingM >> _maxweight >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3;
}

ClassDescription<SU3BaryonDecupletOctetScalarDecayer> SU3BaryonDecupletOctetScalarDecayer::initSU3BaryonDecupletOctetScalarDecayer;
// Definition of the static class description member.

void SU3BaryonDecupletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetScalarDecayer> documentation
    ("The \\classname{SU3BaryonDecupletOctetScalarDecayer} class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a pseudoscalar"
     " meson from the lightest multiplet.");

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,double> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetScalarDecayer::_C, 1.3, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetScalarDecayer,int> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetScalarDecayer::_parity, 0, false, false);
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

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonDecupletOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_proton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_neutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_lambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xim, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltapp, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_delta0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasp, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmas0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasm, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_omega, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xis0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xism, 0, -100000, 100000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonDecupletOctetScalarDecayer::halfThreeHalfScalarCoupling(int imode,
								      Complex& A,
								      Complex& B) const
 {A=_A1[imode];B=_B1[imode];}

// set up the decay modes
void SU3BaryonDecupletOctetScalarDecayer::setupModes(unsigned int iopt) const
{
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_incomingB.resize(0);_outgoingM.resize(0);}
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.)),orr(1./sqrt(3.));
  // decays of the delta++
  intemp.push_back(_deltapp);outtemp.push_back(_proton);mestemp.push_back(211);
  factor.push_back(_C);
  intemp.push_back(_deltapp);outtemp.push_back(_sigmap);mestemp.push_back(321);
  factor.push_back(-_C);
  // decays of the delta+
  intemp.push_back(_deltap);outtemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back(_C*orr);
  intemp.push_back(_deltap);outtemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigmap);mestemp.push_back(311);
  factor.push_back(_C*orr);
  // decays of the delta0
  intemp.push_back(_delta0);outtemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back(_C*orr);
  intemp.push_back(_delta0);outtemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(_C*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigmam);mestemp.push_back(321);
  factor.push_back(_C*orr);
  // decays of the delta-
  intemp.push_back(_deltam);outtemp.push_back(_neutron);mestemp.push_back(-211);
  factor.push_back(-_C);
  intemp.push_back(_deltam);outtemp.push_back(_sigmam);mestemp.push_back(311);
  factor.push_back(_C);
  // sigma*+
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_proton);mestemp.push_back(-311);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_xi0);mestemp.push_back(321);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(_C*ort);
  intemp.push_back(_sigmasp);outtemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(_C*ort);
  // sigma*0
  intemp.push_back(_sigmas0);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(_C*ort);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_C*ort);
  // sigma*-
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_C*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_neutron);mestemp.push_back(-321);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_xim);mestemp.push_back(311);
  factor.push_back(_C*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(_C*ort);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(_C*ort);
  // xi*0
  intemp.push_back(_xis0);outtemp.push_back(_xim);mestemp.push_back(211);
  factor.push_back(_C*orr);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_xis0);outtemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(_C*orr);
  intemp.push_back(_xis0);outtemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(_C*ors);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(_C*ort);
  intemp.push_back(_xis0);outtemp.push_back(_lambda);mestemp.push_back(-311);
  factor.push_back(_C*ort);
  // xi*-
  intemp.push_back(_xism);outtemp.push_back(_xi0);mestemp.push_back(-211);
  factor.push_back(_C*orr);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(_C*ors);
  intemp.push_back(_xism);outtemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_C*orr);
  intemp.push_back(_xism);outtemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(_C*ors);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(_C*ort);
  intemp.push_back(_xism);outtemp.push_back(_lambda);mestemp.push_back(-321);
  factor.push_back(_C*ort);
  // set up the modes
  Energy msum(0.);
  for(unsigned int ix=0;ix<intemp.size();++ix)
    {
      if(intemp[ix]!=0&&outtemp[ix]!=0&&mestemp[ix]!=0)
	{
	  _incomingB.push_back(intemp[ix]);
	  _outgoingB.push_back(outtemp[ix]);
	  _outgoingM.push_back(mestemp[ix]);
	  if(iopt==1)
	    {
	      msum = getParticleData(_incomingB.back())->mass()+
		getParticleData(_outgoingB.back())->mass();
	      if(_parity==0)
		{_A1.push_back(factor[ix]*msum/_fpi);_B1.push_back(0.);}
	      else
		{_A1.push_back(0.);_B1.push_back(factor[ix]*msum/_fpi);}
	      _A2.push_back(0.);_B2.push_back(0.);
	      _A3.push_back(0.);_B3.push_back(0.);
	    }
	}
    }
}
}
