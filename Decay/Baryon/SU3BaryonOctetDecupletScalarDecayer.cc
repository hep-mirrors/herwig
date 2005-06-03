// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonOctetDecupletScalarDecayer class.
//

#include "SU3BaryonOctetDecupletScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SU3BaryonOctetDecupletScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

SU3BaryonOctetDecupletScalarDecayer::~SU3BaryonOctetDecupletScalarDecayer() {}

bool SU3BaryonOctetDecupletScalarDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
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

ParticleVector SU3BaryonOctetDecupletScalarDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode(-1),id(parent.id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);bool cc;
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


void SU3BaryonOctetDecupletScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _fpi << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _outgoingM << _maxweight << _prefactor;
}

void SU3BaryonOctetDecupletScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _fpi >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _outgoingM >> _maxweight >> _prefactor;
}

ClassDescription<SU3BaryonOctetDecupletScalarDecayer> SU3BaryonOctetDecupletScalarDecayer::initSU3BaryonOctetDecupletScalarDecayer;
// Definition of the static class description member.

void SU3BaryonOctetDecupletScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetDecupletScalarDecayer> documentation
    ("The \\classname{SU3BaryonOctetDecupletScalarDecayer} class is designed for the"
     " decay of an SU(3) octet baryon to an SU(3) decuplet baryon and a pseudoscalar"
     " meson from the lightest multiplet.");

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decuplet decays.",
     &SU3BaryonOctetDecupletScalarDecayer::_C, 1.35, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonOctetDecupletScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetDecupletScalarDecayer::_parity, true, false, false);
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

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonOctetDecupletScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_proton, 12212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_neutron, 12112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmap, 13222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigma0, 13212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmam, 13112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_lambda, 23122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xi0, 13322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xim, 13312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltapp, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltap, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_delta0, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltam, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmasp, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmas0, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmasm, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_omega, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xis0, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xism, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonOctetDecupletScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetDecupletScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonOctetDecupletScalarDecayer::halfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
								      Complex& A,
								      Complex& B) const
{
  if(_parity)
    {
      A=_prefactor[imode]*(m0+m1);
      B=0.;
    }
  else
    {
      A=0.;
      B=_prefactor[imode]*(m0+m1);
    }
}

// couplings for spin-3/2 to spin-3/2 spin-0
void SU3BaryonOctetDecupletScalarDecayer::
threeHalfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,Complex& A1,Complex& A2,
				 Complex& B1,Complex& B2) const
{
  A2=0.;B2=0.;
  if(_parity)
    {
      A1=0.;
      B1=_prefactor[imode]*(m0+m1);
    }
  else
    {
      A1=_prefactor[imode]*(m0-m1);
      B1=0.;
    }
}

// set up the decay modes
void SU3BaryonOctetDecupletScalarDecayer::setupModes(unsigned int iopt) const
{
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.resize(0);_incomingB.resize(0);_outgoingM.resize(0);}
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.)),orr(1./sqrt(3.));
  // decays to the delta++
  outtemp.push_back(_deltapp);intemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back(_C);
  outtemp.push_back(_deltapp);intemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(-_C);
  // decays to the delta+
  outtemp.push_back(_deltap);intemp.push_back(_neutron);mestemp.push_back(-211);
  factor.push_back(_C*orr);
  outtemp.push_back(_deltap);intemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(_C*rt*orr);
  outtemp.push_back(_deltap);intemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(_C*rt*orr);
  outtemp.push_back(_deltap);intemp.push_back(_sigmap);mestemp.push_back(-311);
  factor.push_back(_C*orr);
  // decays to the delta0
  outtemp.push_back(_delta0);intemp.push_back(_proton);mestemp.push_back(211);
  factor.push_back(_C*orr);
  outtemp.push_back(_delta0);intemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(_C*rt*orr);
  outtemp.push_back(_delta0);intemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(_C*rt*orr);
  outtemp.push_back(_delta0);intemp.push_back(_sigmam);mestemp.push_back(-321);
  factor.push_back(_C*orr);
  // decays to the delta-
  outtemp.push_back(_deltam);intemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back(-_C);
  outtemp.push_back(_deltam);intemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_C);
  // decays to the sigma*+
  outtemp.push_back(_sigmasp);intemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmasp);intemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmasp);intemp.push_back(_proton);mestemp.push_back(311);
  factor.push_back(_C*orr);
  outtemp.push_back(_sigmasp);intemp.push_back(_xi0);mestemp.push_back(-321);
  factor.push_back(_C*orr);
  outtemp.push_back(_sigmasp);intemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(_C*ort);
  outtemp.push_back(_sigmasp);intemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(_C*ort);
  // decays to the sigma*0
  outtemp.push_back(_sigmas0);intemp.push_back(_neutron);mestemp.push_back(311);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_proton);mestemp.push_back(321);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_xim);mestemp.push_back(-321);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_xi0);mestemp.push_back(-311);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_sigmam);mestemp.push_back(-211);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_sigmap);mestemp.push_back(211);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(_C*ort);
  outtemp.push_back(_sigmas0);intemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_C*ort);
  // decays to the sigma*-
  outtemp.push_back(_sigmasm);intemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmasm);intemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_C*ors);
  outtemp.push_back(_sigmasm);intemp.push_back(_neutron);mestemp.push_back(321);
  factor.push_back(_C*orr);
  outtemp.push_back(_sigmasm);intemp.push_back(_xim);mestemp.push_back(-311);
  factor.push_back(_C*orr);
  outtemp.push_back(_sigmasm);intemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(_C*ort);
  outtemp.push_back(_sigmasm);intemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(_C*ort);
  // decays to the xi*0
  outtemp.push_back(_xis0);intemp.push_back(_xim);mestemp.push_back(-211);
  factor.push_back(_C*orr);
  outtemp.push_back(_xis0);intemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(_C*ors);
  outtemp.push_back(_xis0);intemp.push_back(_sigmap);mestemp.push_back(321);
  factor.push_back(_C*orr);
  outtemp.push_back(_xis0);intemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(_C*ors);
  outtemp.push_back(_xis0);intemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(_C*ort);
  outtemp.push_back(_xis0);intemp.push_back(_lambda);mestemp.push_back(311);
  factor.push_back(_C*ort);
  // decays to the xi*-
  outtemp.push_back(_xism);intemp.push_back(_xi0);mestemp.push_back(211);
  factor.push_back(_C*orr);
  outtemp.push_back(_xism);intemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(_C*ors);
  outtemp.push_back(_xism);intemp.push_back(_sigmam);mestemp.push_back(311);
  factor.push_back(_C*orr);
  outtemp.push_back(_xism);intemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(_C*ors);
  outtemp.push_back(_xism);intemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(_C*ort);
  outtemp.push_back(_xism);intemp.push_back(_lambda);mestemp.push_back(321);
  factor.push_back(_C*ort);
  // set up the modes
  Energy m0,m1;
  int inspin,outspin;
  PDVector extpart(3);
  for(unsigned int ix=0;ix<outtemp.size();++ix)
    {
      if(outtemp[ix]!=0&&intemp[ix]!=0&&mestemp[ix]!=0)
	{
	  extpart[0]=getParticleData(intemp[ix]);
	  extpart[1]=getParticleData(outtemp[ix]);
	  extpart[2]=getParticleData(mestemp[ix]);
	  if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin())
	    {
	      _incomingB.push_back(intemp[ix]);
	      _outgoingB.push_back(outtemp[ix]);
	      _outgoingM.push_back(mestemp[ix]);
	      if(iopt==1)
		{
		  m0 = extpart[0]->mass();inspin  = extpart[0]->iSpin();
		  m1 = extpart[1]->mass();outspin = extpart[1]->iSpin();
		  if(inspin==2&&outspin==4)
		    {_prefactor.push_back(factor[ix]/_fpi);}
		  else if(inspin==4&&outspin==4)
		    {_prefactor.push_back(factor[ix]/_fpi);}
		}
	    }
	}
    }
}

void SU3BaryonOctetDecupletScalarDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":Points " << _npoint << "\n";
  output << "set " << fullName() << ":Coupling " << _C<< "\n";
  output << "set " << fullName() << ":Parity " << _parity<< "\n";
  output << "set " << fullName() << ":Fpi " << _fpi << "\n";
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
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
