// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeMesonDefaultCurrent class.
//

#include "ThreeMesonDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeMesonDefaultCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

ThreeMesonDefaultCurrent::~ThreeMesonDefaultCurrent() {}

void ThreeMesonDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts << _KstarF123wgts << _rhoF5wgts << _KstarF5wgts
     << _rhoKstarwgt <<  _a1runwidth << _a1runq2 <<  _initializea1
     << _a1mass << _a1width << _K1mass << _K1width << _fpi << _mpi << _mK
     <<_rhoparameters << _rhoF123masses << _rhoF5masses << _rhoF123widths 
     << _rhoF5widths << _Kstarparameters << _KstarF123masses <<_KstarF5masses
     << _KstarF123widths << _KstarF5widths << _a1parameters << _K1parameters;
}

void ThreeMesonDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >> _KstarF123wgts >> _rhoF5wgts >> _KstarF5wgts
     >> _rhoKstarwgt >>  _a1runwidth >> _a1runq2 >>  _initializea1
     >> _a1mass >> _a1width >> _K1mass >> _K1width >> _fpi >> _mpi >> _mK
     >>_rhoparameters >> _rhoF123masses >> _rhoF5masses >> _rhoF123widths 
     >> _rhoF5widths >> _Kstarparameters >> _KstarF123masses >>_KstarF5masses
     >> _KstarF123widths >> _KstarF5widths >> _a1parameters >> _K1parameters;
}

ClassDescription<ThreeMesonDefaultCurrent> ThreeMesonDefaultCurrent::initThreeMesonDefaultCurrent;
// Definition of the static class description member.

void ThreeMesonDefaultCurrent::Init() {
        
  static ClassDocumentation<ThreeMesonDefaultCurrent> documentation
    ("The \\classname{ThreeMesonDefaultCurrent} class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.");
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF123KstarWgt
    ("F123KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_KstarF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF5RhoWgt
    ("F5RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_rhoF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF5KstarWgt
    ("F5KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_KstarF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Parameter<ThreeMesonDefaultCurrent,double> interfaceRhoKstarWgt
    ("RhoKstarWgt",
     "The relative weights of the rho and K* in the F5 form factor",
     &ThreeMesonDefaultCurrent::_rhoKstarwgt, -0.2, -1.0e12, 1.0e12,
     false, false, false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &ThreeMesonDefaultCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Initialization",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "NoInitialization",
     "Use the default values",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_rhoparameters, true, false, false);
  static SwitchOption interfaceRhoParameterstrue
    (interfaceRhoParameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceRhoParametersParticleData
    (interfaceRhoParameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_Kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
       "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfacea1Parameters
    ("a1Parameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_a1parameters, true, false, false);
  static SwitchOption interfacea1Parameterstrue
    (interfacea1Parameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfacea1ParametersParticleData
    (interfacea1Parameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceK1Parameters
    ("K1Parameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_K1parameters, true, false, false);
  static SwitchOption interfaceK1Parameterstrue
    (interfaceK1Parameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceK1ParametersParticleData
    (interfaceK1Parameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &ThreeMesonDefaultCurrent::_a1runwidth,
     0, 0, 0, 0, 10000000, false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &ThreeMesonDefaultCurrent::_a1runq2,
     0, 0, 0, 0, 10000000, false, false, true);
  
  
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &ThreeMesonDefaultCurrent::_a1width, MeV, 559.0*MeV, 0*MeV, 10000*MeV,
     false, false, false);
  
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &ThreeMesonDefaultCurrent::_a1mass, MeV, 1251.0*MeV, 0*MeV, 10000*MeV,
     false, false, false);
  
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceK1Width
    ("K1Width",
     "The K_1 width if using local values.",
     &ThreeMesonDefaultCurrent::_K1width, MeV, 559.0*MeV, 0*MeV, 10000*MeV,
     false, false, false);
  
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceK1Mass
    ("K1Mass",
     "The K_1 mass if using local values.",
     &ThreeMesonDefaultCurrent::_K1mass, MeV, 1251.0*MeV, 0*MeV, 10000*MeV,
     false, false, false);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF123masses,
     0, 0, 0, 0, 10000., false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF123widths,
     0, 0, 0, 0, 10000., false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfacerhoF5masses
    ("rhoF5masses",
     "The masses for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF5masses,
     0, 0, 0, 0, 10000., false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfacerhoF5widths
    ("rhoF5widths",
     "The widths for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF5widths,
     0, 0, 0, 0, 10000., false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceKstarF123masses
    ("KstarF123masses",
     "The masses for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF123masses,
     0, 0, 0, 0, 10000., false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceKstarF123widths
    ("KstarF123widths",
     "The widths for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF123widths,
     0, 0, 0, 0, 10000., false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceKstarF5masses
    ("KstarF5masses",
     "The masses for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF5masses,
     0, 0, 0, 0, 10000., false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceKstarF5widths
    ("KstarF5widths",
     "The widths for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF5widths,
     0, 0, 0, 0, 10000., false, false, true);

}
  
// modes handled by this class
bool ThreeMesonDefaultCurrent::acceptMode(int imode) const{return imode>=0&&imode<=8;}

// calculate the form-factors
void ThreeMesonDefaultCurrent::
calculateFormFactors(const int ichan, const int imode,
		     Energy2 q2, Energy2 s1, Energy2 s2, Energy2 s3,
		     Complex & F1, Complex & F2, Complex & F3, Complex & F4,
		     Complex & F5) const
{
  F1=0.;F2=0.;F3=0.;F4=0.;F5=0.;
  // calculate the pi- pi- pi+ factor
  if(imode==0)
    {
      Complex a1fact=a1BreitWigner(q2)*2./3.;
      if(ichan<0){F1= a1fact*BrhoF123(s1,-1);F2 =-a1fact*BrhoF123(s2,-1);}
      else if(ichan%2==0){F1 = a1fact*BrhoF123(s1,ichan/2);}
      else if(ichan%2==1){F2 =-a1fact*BrhoF123(s2,(ichan-1)/2);}
    }
  // calculate the pi0 pi0 pi- factor
  else if(imode==1)
    {
      Complex a1fact=a1BreitWigner(q2)*2./3.;
      if(ichan<0){F1 = a1fact*BrhoF123(s1,-1);F2 =-a1fact*BrhoF123(s2,-1);}
      else if(ichan%2==0){F1 = a1fact*BrhoF123(s1,ichan/2);}
      else if(ichan%2==1){F2 =-a1fact*BrhoF123(s2,(ichan-1)/2);}
    }
  // calculate the K- pi - K+ factor
  else if(imode==2)
    {
      Complex a1fact=a1BreitWigner(q2)*sqrt(2.)/3.;
      if(ichan<0)
	{
	  F1 =-a1fact*BKstarF123(s1,-1); F2 = a1fact*BrhoF123(s2,-1);
	  F5 = BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
	}
      else if(ichan%8==0){F1 =-a1fact*BKstarF123(s1,ichan/8);}
      else if(ichan%8==1){F2 = a1fact*BrhoF123(s2,(ichan-1)/8);}
      else if(ichan%8>=2){F5 = BrhoF5(q2,ichan/8)*FKrho(s1,s2,(ichan-2)%8)*sqrt(2.);}
    }
  // calculate the K0 pi- k0bar
  else if(imode==3)
    {
      Complex a1fact=a1BreitWigner(q2)*sqrt(2.)/3.;
      if(ichan<0)
	{
	  F1 =-a1fact*BKstarF123(s1,-1);F2 = a1fact*BrhoF123(s2,-1);
	  F5 =-BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
	}
      else if(ichan%8==0){F1 = -a1fact*BKstarF123(s1,ichan/8);}
      else if(ichan%8==1){F2 = a1fact*BrhoF123(s2,(ichan-1)/8);}
      else if(ichan%8>=2){F5 = -BrhoF5(q2,ichan/8)*FKrho(s1,s2,(ichan-2)%8)*sqrt(2.);}
    }
  // calculate the K- pi0 k0
  else if(imode==4)
    {
      Complex a1fact=a1BreitWigner(q2);
      if(ichan<0){F2 =-a1fact*BrhoF123(s2,-1);}
      else{F2 =-a1fact*BrhoF123(s2,ichan);}
      int iq,ia;
      decayModeInfo(imode,iq,ia);
    }
  // calculate the pi0 pi0 K-
  else if(imode==5)
    {
      Complex K1fact=K1BreitWigner(q2)/6.;
      if(ichan<0){F1 = K1fact*BKstarF123(s1,-1);F2 =-K1fact*BKstarF123(s2,-1);}
      else if(ichan%2==0){F1 = K1fact*BKstarF123(s1,ichan/2);}
      else{F2 =-K1fact*BKstarF123(s2,(ichan-1)/2);}
    }
  // calculate the K- pi- pi+
  else if(imode==6)
    {
      Complex K1fact=K1BreitWigner(q2)*sqrt(2.)/3.;
      if(ichan<0)
	{
	  F1 =-K1fact*BrhoF123(s1,-1);F2 = K1fact*BKstarF123(s2,-1);
	  F5 =-BKstarF123(q2,-1)*FKrho(s2,s1,-1)*sqrt(2.);
	}
      else if(ichan%8==0){F1 =-K1fact*BrhoF123(s1,ichan/8);}
      else if(ichan%8==1){F2 = K1fact*BKstarF123(s2,(ichan-1)/8);}
      else{F5 = -BKstarF123(q2,ichan/8)*FKrho(s2,s1,(ichan-2)%8)*sqrt(2.);}
    }
  // calculate the pi- K0bar pi0
  else if(imode==7)
    {
      Complex K1fact=K1BreitWigner(q2);
      if(ichan<0){F2 =-K1fact*BrhoF123(s2,-1);F5 =-2.*BKstarF123(q2,-1)*FKrho(s1,s2,-1);}
      else if(ichan%7==0){F2 =-K1fact*BrhoF123(s2,ichan/7);}
      else {F5 =-2.*BKstarF123(q2,ichan/7)*FKrho(s1,s2,(ichan-1)%7);}
    }
  // calculate the pi- pi0 eta
  else if(imode==8)
    {
      if(ichan<0){F5 = BrhoF5(q2,-1)*BrhoF123(s3,-1)*sqrt(2./3.);}
      else{F5 = BrhoF5(q2,ichan/3)*BrhoF123(s3,ichan%3)*sqrt(2./3.);}
    }
  // multiply by the prefactors
  F1/=_fpi;F2/=_fpi;F3/=_fpi;F4/=_fpi;F5/=_fpi;
  F5 =-F5*Complex(0.,1.)/4./pi/pi/_fpi/_fpi;
}

// complete the construction of the decay mode for integration
bool ThreeMesonDefaultCurrent::createMode(int icharge, unsigned int imode,
					  DecayPhaseSpaceModePtr mode,
					  unsigned int iloc,unsigned int ires,
					  DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  bool kineallowed=true;
  if(!acceptMode(imode)){return false;}
  int iq,ia;
  PDVector extpart=particles(1,imode,iq,ia);
  Energy min(0.);
  for(unsigned int ix=0;ix<extpart.size();++ix)
    {min+=extpart[0]->mass()-extpart[0]->widthLoCut();}
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  // the particles we will use a lot
  tPDPtr a1,k1;
  if(icharge==-3)
    {
      a1=getParticleData(ParticleID::a_1minus);
      k1=getParticleData(ParticleID::K_1minus);
    }
  else if(icharge==3)
    {
      a1=getParticleData(ParticleID::a_1plus);
      k1=getParticleData(ParticleID::K_1plus);
    }
  else
    {return false;}
  // the rho0 resonances
  tPDPtr rho0[3];
  rho0[0] = getParticleData(113);
  rho0[1] = getParticleData(100113);
  rho0[2] = getParticleData(30113);
  tPDPtr rhoc[3],Kstar0[3],Kstarc[3];
  if(icharge==-3)
    {
      // the charged rho resonances
      rhoc[0] = getParticleData(-213);
      rhoc[1] = getParticleData(-100213);
      rhoc[2] = getParticleData(-30213);
      // the K*0 resonances
      Kstar0[0] = getParticleData(313);
      Kstar0[1] = getParticleData(100313);
      Kstar0[2] = getParticleData(30313);
      // the charged K* resonances
      Kstarc[0] = getParticleData(-323);
      Kstarc[1] = getParticleData(-100323);
      Kstarc[2] = getParticleData(-30323);
    }
  else
    {
      // the charged rho resonances
      rhoc[0] = getParticleData(213);
      rhoc[1] = getParticleData(100213);
      rhoc[2] = getParticleData(30213);
      // the K*0 resonances
      Kstar0[0] = getParticleData(-313);
      Kstar0[1] = getParticleData(-100313);
      Kstar0[2] = getParticleData(-30313);
      // the charged K* resonances
      tPDPtr Kstarc[3];
      Kstarc[0] = getParticleData(323);
      Kstarc[1] = getParticleData(100323);
      Kstarc[2] = getParticleData(30323);
    }
  DecayPhaseSpaceChannelPtr newchannel;
  if(imode==0)
    {
      // channels for pi- pi- pi+
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc+1,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==1)
    {
      // channels for pi0 pi0 pi-
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc+1,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==2)
    {
      // channels for K- pi- K+
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstar0[ix],0,0.0, iloc+1,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rho0[ix],0,0.0,iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc+1,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(rho0[iy],0,0.0,iloc,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==3)
    {
      // channels for K0 pi- K0bar
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc+1,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(rho0[iy],0,0.0, iloc,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==4)
    {
      // channels for K- pi0 K0
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==5)
    {  
      // channels for pi0 pi0 K-
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1        ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1        ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==6)
    {
      // channels for K- pi- pi+
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1      ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc+1,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1        ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(Kstar0[ix],0,0.0, iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(rho0[iy]  ,0,0.0, iloc+1,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==7)
    {
      // channels for pi- kbar0 pi0
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	  newchannel->init();
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc+1,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(rhoc[iy]  ,0,0.0, iloc,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==8)
    {
      // channels for pi- pi0 eta
      for(unsigned int ix=0;ix<3;++ix)
	{
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix],0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(rho0[iy],0,0.0, iloc+1,iloc+2);
	      newchannel->init();
	      mode->addChannel(newchannel);
	    }
	}
    }
  if(_rhoparameters)
    {
      for(unsigned int ix=0;ix<_rhoF123masses.size();++ix)
	{
	  mode->resetIntermediate(rhoc[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
	  mode->resetIntermediate(rho0[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
	}
    }
  // K star parameters in the base class
  if(_Kstarparameters)
    {
      for(unsigned int ix=0;ix<_KstarF123masses.size();++ix)
	{
	  mode->resetIntermediate(Kstarc[ix],_KstarF123masses[ix],_KstarF123widths[ix]);
	  mode->resetIntermediate(Kstar0[ix],_KstarF123masses[ix],_KstarF123widths[ix]);
	}
    }
  return kineallowed;
}


  PDVector ThreeMesonDefaultCurrent::particles(int icharge, unsigned int imode,int iq,
					       int ia)
{
  PDVector extpart(3);
  if(imode==0)
    {
      extpart[0]=getParticleData(ParticleID::piminus);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::piplus);
    }
  else if(imode==1)
    {
      extpart[0]=getParticleData(ParticleID::pi0);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::piminus);
    }
  else if(imode==2)
    {
      extpart[0]=getParticleData(ParticleID::Kminus);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::Kplus);
    }
  else if(imode==3)
    {
      extpart[0]=getParticleData(ParticleID::K0);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::Kbar0);
    }
  else if(imode==4)
    {
      extpart[0]=getParticleData(ParticleID::Kminus);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::K0);
    }
  else if(imode==5)
    {
      extpart[0]=getParticleData(ParticleID::pi0);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::Kminus);
    }
  else if(imode==6)
    {
      extpart[0]=getParticleData(ParticleID::Kminus);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::piplus);
    }
  else if(imode==7)
    {
      extpart[0]=getParticleData(ParticleID::piminus);
      extpart[1]=getParticleData(ParticleID::Kbar0);
      extpart[2]=getParticleData(ParticleID::pi0);
    }
  else if(imode==8)
    {
      extpart[0]=getParticleData(ParticleID::piminus);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::eta);
    }
  // conjugate the particles if needed
  if(icharge==3)
    {
      for(unsigned int ix=0;ix<3;++ix)
	{if(extpart[0]->CC()){extpart[0]=extpart[0]->CC();}}
    }
  // return the answer
  return extpart;
}

void ThreeMesonDefaultCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/ThreeMesonDefaultCurrent " << fullName() << " \n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix)
    {output << "insert " << fullName() << ":F123RhoWeight " << ix 
	    << " " << _rhoF123wgts[ix] << "\n";}
  for(unsigned int ix=0;ix<_KstarF123wgts.size();++ix)
    {output << "insert " << fullName() << ":F123KstarWeight " << ix 
	    << " " << _KstarF123wgts[ix] << "\n";}
  for(unsigned int ix=0;ix<_rhoF5wgts.size();++ix)
    {output << "insert " << fullName() << ":F5RhoWeight " << ix 
	    << " " << _rhoF5wgts[ix] << "\n";}
  for(unsigned int ix=0;ix<_KstarF5wgts.size();++ix)
    {output << "insert " << fullName() << ":F5KstarWeight " << ix 
	    << " " << _KstarF5wgts[ix] << "\n";}
  output << "set " << fullName() << ":RhoKstarWgt " << _rhoKstarwgt << "\n";
  output << "set " << fullName() << ":Initializea1 " << _initializea1 << "\n";
  output << "set " << fullName() << ":RhoParameters " << _rhoparameters << "\n";
  output << "set " << fullName() << ":KstarParameters " << _Kstarparameters << "\n";
  output << "set " << fullName() << ":a1Parameters " << _a1parameters << "\n";
  output << "set " << fullName() << ":K1Parameters " << _K1parameters << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix)
    {output << "insert " << fullName() << ":a1RunningWidth " << ix 
	    << " " << _a1runwidth[ix] << "\n";}
  for(unsigned int ix=0;ix<_a1runq2.size();++ix)
    {output << "insert " << fullName() << ":a1RunningQ2 " << ix 
	    << " " << _a1runq2[ix] << "\n";}
  output << "set " << fullName() << ":A1Width " << _a1width/MeV << "\n";
  output << "set " << fullName() << ":A1Mass " << _a1mass/MeV << "\n";
  output << "set " << fullName() << ":K1Width " << _K1width/MeV << "\n";
  output << "set " << fullName() << ":K1Mass " << _K1mass/MeV << "\n";
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix)
    {output << "insert " << fullName() << ":rhoF123masses " << ix 
	    << " " << _rhoF123masses[ix] << "\n";}
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix)
    {output << "insert " << fullName() << ":rhoF123widths " << ix 
	    << " " << _rhoF123widths[ix] << "\n";}
  for(unsigned int ix=0;ix<_rhoF5masses.size();++ix)
    {output << "insert " << fullName() << ":rhoF5masses " << ix 
	    << " " << _rhoF5masses[ix] << "\n";}
  for(unsigned int ix=0;ix<_rhoF5widths.size();++ix)
    {output << "insert " << fullName() << ":rhoF5widths " << ix 
	    << " " << _rhoF5widths[ix] << "\n";}
  for(unsigned int ix=0;ix<_KstarF123masses.size();++ix)
    {output << "insert " << fullName() << ":KstarF123masses " << ix 
	    << " " << _KstarF123masses[ix] << "\n";}
  for(unsigned int ix=0;ix<_KstarF123widths.size();++ix)
    {output << "insert " << fullName() << ":KstarF123widths " << ix 
	    << " " << _KstarF123widths[ix] << "\n";}
  for(unsigned int ix=0;ix<_KstarF5masses.size();++ix)
    {output << "insert " << fullName() << ":KstarF5masses " << ix 
	    << " " << _KstarF5masses[ix] << "\n";}
  for(unsigned int ix=0;ix<_KstarF5widths.size();++ix)
    {output << "insert " << fullName() << ":KstarF5widths " << ix 
	    << " " << _KstarF5widths[ix] << "\n";}
}
  
}

// the functions for the integrands of the a_1 width
namespace Herwig {
using namespace Genfun;

FUNCTION_OBJECT_IMP(Defaulta1MatrixElement)

Defaulta1MatrixElement::Defaulta1MatrixElement(Ptr<Herwig::ThreeMesonDefaultCurrent>::pointer in)
  {_decayer=in;}

  
Defaulta1MatrixElement::~Defaulta1MatrixElement() { }
  
Defaulta1MatrixElement::Defaulta1MatrixElement(const Defaulta1MatrixElement & right)  {  }
  
unsigned int Defaulta1MatrixElement::dimensionality() const {return 7;}

double Defaulta1MatrixElement::operator ()(const Argument & a) const 
{return _decayer->a1MatrixElement(a[0],a[1],a[2],a[3],a[4],a[5],a[6]);}

}
