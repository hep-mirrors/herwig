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
      //      cout << "testing the decays " << ix << "  " 
      //	   << getParticleData(_incomingB[ix])->PDGName() << " -> " 
      //	   << getParticleData(_outgoingB[ix])->PDGName() << "  "
      //   << getParticleData(_outgoingM[ix])->PDGName() << endl;
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


void SU3BaryonOctetDecupletScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _C << _parity << _fpi << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _outgoingM << _maxweight << _A1 << _A2 << _A3 << _B1 << _B2 << _B3;
}

void SU3BaryonOctetDecupletScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _C >> _parity >> _fpi >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _outgoingM >> _maxweight >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3;
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
     &SU3BaryonOctetDecupletScalarDecayer::_C, 1.3, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonOctetDecupletScalarDecayer,int> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetDecupletScalarDecayer::_parity, 0, false, false);
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

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonOctetDecupletScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_proton, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_neutron, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigma0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_lambda, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xi0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xim, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltapp, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltap, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_delta0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltam, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmasp, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmas0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmasm, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_omega, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xis0, 0, -100000, 100000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xism, 0, -100000, 100000,
     false, false, true);

  static ParVector<SU3BaryonOctetDecupletScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetDecupletScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonOctetDecupletScalarDecayer::halfThreeHalfScalarCoupling(int imode,
								      Complex& A,
								      Complex& B) const
 {A=_A1[imode];B=_B1[imode];}

// couplings for spin-3/2 to spin-3/2 spin-0
void SU3BaryonOctetDecupletScalarDecayer::
threeHalfThreeHalfScalarCoupling(int imode,Complex& A1,Complex& A2,
				 Complex& B1,Complex& B2) const
{
  A1 = _A1[imode];
  A2 = _A2[imode];
  B1 = _B1[imode];
  B2 = _B2[imode];
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
  for(unsigned int ix=0;ix<outtemp.size();++ix)
    {
      if(outtemp[ix]!=0&&intemp[ix]!=0&&mestemp[ix]!=0)
	{
	  _incomingB.push_back(intemp[ix]);
	  _outgoingB.push_back(outtemp[ix]);
	  _outgoingM.push_back(mestemp[ix]);
	  if(iopt==1)
	    {
	      m0 = getParticleData(_incomingB.back())->mass();
	      m1 = getParticleData(_outgoingB.back())->mass();
	      inspin  = getParticleData(_incomingB.back())->iSpin();
	      outspin = getParticleData(_outgoingB.back())->iSpin();
	      if(inspin==2&&outspin==4)
		{
		  if(_parity==0)
		    {_A1.push_back(factor[ix]*(m0+m1)/_fpi);_B1.push_back(0.);}
		  else
		    {_A1.push_back(0.);_B1.push_back(factor[ix]*(m0+m1)/_fpi);}
		  
		  _A2.push_back(0.);_B2.push_back(0.);
		  _A3.push_back(0.);_B3.push_back(0.);
		}
	      else if(inspin==4&&outspin==4)
		{
		  if(_parity==0)
		    {_A1.push_back(factor[ix]*(m0-m1)/_fpi);_B1.push_back(0.);}
		  else
		    {_A1.push_back(0.);_B1.push_back(factor[ix]*(m0+m1)/_fpi);}
		  _A2.push_back(0.);_B2.push_back(0.);
		  _A3.push_back(0.);_B3.push_back(0.);
		}
	    }
	}
    }
}
}
