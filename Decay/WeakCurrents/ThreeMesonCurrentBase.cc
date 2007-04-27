// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeMesonCurrentBase class.
//

#include "ThreeMesonCurrentBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"


using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::ScalarSpinInfo;

ThreeMesonCurrentBase::ThreeMesonCurrentBase() {
  // the quarks for the different modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(12);
}

AbstractNoPIOClassDescription<ThreeMesonCurrentBase> 
ThreeMesonCurrentBase::initThreeMesonCurrentBase;
// Definition of the static class description member.

void ThreeMesonCurrentBase::Init() {
    
  static ClassDocumentation< ThreeMesonCurrentBase> documentation
    ("The ThreeMesonCurrentBase class is designed to be the "
     "base class for "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta, K0S pi- K0S, K0L pi- K0L, K0S pi- K0L");

}

// the hadronic currents    
vector<LorentzPolarizationVector> 
ThreeMesonCurrentBase::current(bool vertex, const int imode, const int ichan, 
			       Energy & scale,const ParticleVector & decay) const {
  // spininfo for the particles
  if(vertex) {
    for(unsigned int ix=0;ix<3;++ix) {
      decay[ix]->spinInfo(new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true)));
    }
  }
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q=0;
  for(unsigned int ix=0;ix<decay.size();++ix){q+=decay[ix]->momentum();}
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (decay[1]->momentum()+decay[2]->momentum()).m2();
  Energy2 s2 = (decay[0]->momentum()+decay[2]->momentum()).m2();
  Energy2 s3 = (decay[0]->momentum()+decay[1]->momentum()).m2();
  Complex F1(0.),F2(0.),F3(0.),F4(0.),F5(0.);
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
  vect += (F4-dot)*q;
  if(F5!=0.) vect += Complex(0.,1.)*F5*Helicity::EpsFunction::product(decay[0]->momentum(),
								      decay[1]->momentum(),
								      decay[2]->momentum());
  // factor to get dimensions correct
  return vector<LorentzPolarizationVector>(1,q.mass()*vect);
}

bool ThreeMesonCurrentBase::accept(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nk0(0),nk0bar(0),neta(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
    else if(id[ix]==ParticleID::eta)     ++neta;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  int imode(-1);
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) imode= 0;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) imode= 1;
  else if( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))                 imode= 2;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))              imode= 3;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))                 imode= 4;
  else if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   imode= 5;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )               imode= 6;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))             imode= 7;
  else if( (npip==1&&npi0==1&&neta==1) ||
	   (npim==1&&npi0==1&&neta==1))               imode= 8;
  else if( nks==2 && (npip==1||npim==1) )             imode= 9;
  else if( nkl==2 && (npip==1||npim==1) )             imode=10;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) )     imode=11;
  return imode==-1 ? false : acceptMode(imode);
}

unsigned int ThreeMesonCurrentBase::decayMode(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nk0(0),nk0bar(0),neta(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
    else if(id[ix]==ParticleID::eta)     ++neta;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  int imode(-1);
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) imode= 0;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) imode= 1;
  else if( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))                 imode= 2;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))              imode= 3;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))                 imode= 4;
  else if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   imode= 5;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )               imode= 6;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))             imode= 7;
  else if( (npip==1&&npi0==1&&neta==1) ||
	   (npim==1&&npi0==1&&neta==1))               imode= 8;
  else if( nks==2 && (npip==1||npim==1) )             imode= 9;
  else if( nkl==2 && (npip==1||npim==1) )             imode=10;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) )     imode=11;
  return imode;
}

void ThreeMesonCurrentBase::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  WeakDecayCurrent::dataBaseOutput(output,header,create);
}

PDVector ThreeMesonCurrentBase::particles(int icharge, unsigned int imode,int,int) {
  PDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::piminus);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kplus);
  }
  else if(imode==3) {
    extpart[0]=getParticleData(ParticleID::K0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==4) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::K0);
  }
  else if(imode==5) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::Kminus);
  }
  else if(imode==6) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==7) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::Kbar0);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else if(imode==8) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::eta);
  }
  else if(imode==9) {
    extpart[0]=getParticleData(ParticleID::K_S0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_S0);
  }
  else if(imode==10) {
    extpart[0]=getParticleData(ParticleID::K_L0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_L0);
  }
  else if(imode==11) {
    extpart[0]=getParticleData(ParticleID::K_S0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_L0);
  }
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}

