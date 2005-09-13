// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeMesonCurrentBase class.
//

#include "ThreeMesonCurrentBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeMesonCurrentBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"


namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::ScalarSpinInfo;

ThreeMesonCurrentBase::~ThreeMesonCurrentBase() {}

void ThreeMesonCurrentBase::persistentOutput(PersistentOStream & os) const {
}

void ThreeMesonCurrentBase::persistentInput(PersistentIStream & is, int) {
}

AbstractClassDescription<ThreeMesonCurrentBase> ThreeMesonCurrentBase::initThreeMesonCurrentBase;
// Definition of the static class description member.

void ThreeMesonCurrentBase::Init() {
    
  static ClassDocumentation< ThreeMesonCurrentBase> documentation
    ("The \\classname{ThreeMesonCurrentBase} class is designed to be the "
     "base class for "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta.");

}

// the hadronic currents    
vector<LorentzPolarizationVector> 
ThreeMesonCurrentBase::current(bool vertex, const int imode, const int ichan, 
			       Energy & scale,const ParticleVector & decay) const
{
  // storage for the currents
  vector<LorentzPolarizationVector> temp;
  // spininfo for the particles
  if(vertex)
    {for(unsigned int ix=0;ix<3;++ix)
	{decay[ix]->spinInfo(new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true)));}}
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q=0;
  for(unsigned int ix=0;ix<decay.size();++ix){q+=decay[ix]->momentum();}
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.m2();
  Energy2 s1 = decay[1]->momentum().m2()+decay[2]->momentum().m2()+
    2.*decay[1]->momentum().dot(decay[2]->momentum());
  Energy2 s2 = decay[0]->momentum().m2()+decay[2]->momentum().m2()+
    2.*decay[0]->momentum().dot(decay[2]->momentum());
  Energy2 s3 = decay[0]->momentum().m2()+decay[1]->momentum().m2()+
    2.*decay[0]->momentum().dot(decay[1]->momentum());
  complex<double> F1,F2,F3,F4,F5;
  calculateFormFactors(ichan,imode,q2,s1,s2,s3,F1,F2,F3,F4,F5);
  //if(inpart.id()==ParticleID::tauplus){F5=conj(F5);}
  // the first three form-factors
  LorentzPolarizationVector vect;
  vect = (F2-F1)*decay[2]->momentum()
        +(F1-F3)*decay[1]->momentum()
        +(F3-F2)*decay[0]->momentum();
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect += (F4-dot)*q
    +F5*Helicity::EpsFunction::product(decay[0]->momentum(),decay[1]->momentum(),
				       decay[2]->momentum());
  // factor to get dimensions correct
  temp.push_back(q.mass()*vect);
  return temp;
}

bool ThreeMesonCurrentBase::accept(vector<int> id)
{
  bool allowed=false;
  int npiplus(0),npiminus(0),nkplus(0),nkminus(0),npi0(0),nk0(0),nk0bar(0),neta(0);
  for(unsigned int ix=0;ix<id.size();++ix)
    {
      if(id[ix]==ParticleID::piplus){++npiplus;}
      else if(id[ix]==ParticleID::piminus){++npiminus;}
      else if(id[ix]==ParticleID::Kplus){++nkplus;}
      else if(id[ix]==ParticleID::Kminus){++nkminus;}
      else if(id[ix]==ParticleID::pi0){++npi0;}
      else if(id[ix]==ParticleID::K0){++nk0;}
      else if(id[ix]==ParticleID::Kbar0){++nk0bar;}
      else if(id[ix]==ParticleID::eta){++neta;}
    }
  int imode(-1);
  if((npiplus==2&&npiminus==1)||(npiminus==2&&npiplus==1)){imode=0;}
  else if((npiplus==1&&npi0==2)||(npiminus==1&npi0==2)){imode=1;}
  else if((nkplus==1&&nkminus==1&&npiplus==1)||
	  (nkplus==1&nkminus==1&npiminus==1)){imode=2;}
  else if((nk0==1&&nk0bar==1&&npiplus==1)||(nk0==1&&nk0bar==1&npiminus==1)){imode=3;}
  else if((nkplus==1&&nk0bar==1&&npi0==1)||(nkminus==1&npi0==1&&nk0==1)){imode=4;}
  else if((nkplus==1&&npi0==2)||(npi0==2&nkminus==1)){imode=5;}
  else if((npiplus==1&&npiminus==1&nkplus==1)||
	  (nkminus==1&npiminus==1&npiplus==1)){imode=6;}
  else if((nk0==1&&npiplus==1&&npi0==1)||(npiminus==1&nk0bar==1&npi0==1)){imode=7;}
  else if((npiplus==1&&npi0==1&neta==1)||(npiminus==1&npi0==1&neta==1)){imode=8;}
  if(imode==-1){allowed=false;}
  else{allowed = acceptMode(imode);}
  return allowed;
}

unsigned int ThreeMesonCurrentBase::decayMode(vector<int> id)
{
  int npiplus(0),npiminus(0),nkplus(0),nkminus(0),npi0(0),nk0(0),nk0bar(0),neta(0);
  for(unsigned int ix=0;ix<id.size();++ix)
    {
      if(id[ix]==ParticleID::piplus){++npiplus;}
      else if(id[ix]==ParticleID::piminus){++npiminus;}
      else if(id[ix]==ParticleID::Kplus){++nkplus;}
      else if(id[ix]==ParticleID::Kminus){++nkminus;}
      else if(id[ix]==ParticleID::pi0){++npi0;}
      else if(id[ix]==ParticleID::K0){++nk0;}
      else if(id[ix]==ParticleID::Kbar0){++nk0bar;}
      else if(id[ix]==ParticleID::eta){++neta;}
    }
  int imode(-1);
  if((npiplus==2&&npiminus==1)||(npiminus==2&&npiplus==1)){imode=0;}
  else if((npiplus==1&&npi0==2)||(npiminus==1&npi0==2)){imode=1;}
  else if((nkplus==1&&nkminus==1&&npiplus==1)||
	  (nkplus==1&nkminus==1&npiminus==1)){imode=2;}
  else if((nk0==1&&nk0bar==1&&npiplus==1)||(nk0==1&&nk0bar==1&npiminus==1)){imode=3;}
  else if((nkplus==1&&nk0bar==1&&npi0==1)||(nkminus==1&npi0==1&&nk0==1)){imode=4;}
  else if((nkplus==1&&npi0==2)||(npi0==2&nkminus==1)){imode=5;}
  else if((npiplus==1&&npiminus==1&nkplus==1)||
	  (nkminus==1&npiminus==1&npiplus==1)){imode=6;}
  else if((nk0==1&&npiplus==1&&npi0==1)||(npiminus==1&nk0bar==1&npi0==1)){imode=7;}
  else if((npiplus==1&&npi0==1&neta==1)||(npiminus==1&npi0==1&neta==1)){imode=8;}
  return imode;
}

void ThreeMesonCurrentBase::dataBaseOutput(ofstream & output,bool header,
					   bool create) const
{WeakDecayCurrent::dataBaseOutput(output,header,create);}

}
