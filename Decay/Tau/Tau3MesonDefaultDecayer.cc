// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau3MesonDefaultDecayer class.
//

#include "Tau3MesonDefaultDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
  using namespace ThePEG;
  
  Tau3MesonDefaultDecayer::~Tau3MesonDefaultDecayer() {}
  
  void Tau3MesonDefaultDecayer::persistentOutput(PersistentOStream & os) const {
    os << _rhoF123wgts << _KstarF123wgts << _rhoF5wgts << _KstarF5wgts
       << _rhoKstarwgt <<  _a1runwidth << _a1runq2 <<  _initializea1
       << _a1mass << _a1width << _K1mass << _K1width << _sinfact << _cosfact
       << _fpi << _mpi << _mK  
       <<_rhoparameters << _rhoF123masses << _rhoF5masses << _rhoF123widths 
       << _rhoF5widths << _Kstarparameters << _KstarF123masses <<_KstarF5masses 
       << _KstarF123widths << _KstarF5widths << _a1parameters << _K1parameters
       << _pimpimpipchan<<_pi0pi0pimchan<<_KmpimKpchan<<_K0pimK0chan
       <<_Kmpi0K0chan<<_pi0pi0Kmchan<<_Kmpimpipchan<<_pimK0pi0chan
       << _pimpi0etachan << _pimpimpipwgts<<_pi0pi0pimwgts<<_KmpimKpwgts
       << _K0pimK0wgts<<_Kmpi0K0wgts<<_pi0pi0Kmwgts<<_Kmpimpipwgts<<_pimK0pi0wgts
       << _pimpi0etawgts  << _pimpimpipmax<<_pi0pi0pimmax<<_KmpimKpmax<<_K0pimK0max
       <<_Kmpi0K0max<<_pi0pi0Kmmax<<_Kmpimpipmax<<_pimK0pi0max<< _pimpi0etamax;
  }

  void Tau3MesonDefaultDecayer::persistentInput(PersistentIStream & is, int) {
    is >> _rhoF123wgts >> _KstarF123wgts >> _rhoF5wgts >> _KstarF5wgts
       >> _rhoKstarwgt >>  _a1runwidth >> _a1runq2 >>  _initializea1
       >> _a1mass >> _a1width >> _K1mass >> _K1width >> _sinfact >> _cosfact
       >> _fpi >> _mpi >> _mK
       >>_rhoparameters >> _rhoF123masses >> _rhoF5masses >> _rhoF123widths 
       >> _rhoF5widths >> _Kstarparameters >> _KstarF123masses >>_KstarF5masses 
       >> _KstarF123widths >> _KstarF5widths >> _a1parameters >> _K1parameters
       >> _pimpimpipchan>>_pi0pi0pimchan>>_KmpimKpchan>>_K0pimK0chan
       >>_Kmpi0K0chan>>_pi0pi0Kmchan>>_Kmpimpipchan>>_pimK0pi0chan
       >> _pimpi0etachan >> _pimpimpipwgts>>_pi0pi0pimwgts>>_KmpimKpwgts
       >> _K0pimK0wgts>>_Kmpi0K0wgts>>_pi0pi0Kmwgts>>_Kmpimpipwgts>>_pimK0pi0wgts
       >> _pimpi0etawgts  >> _pimpimpipmax>>_pi0pi0pimmax>>_KmpimKpmax>>_K0pimK0max
       >>_Kmpi0K0max>>_pi0pi0Kmmax>>_Kmpimpipmax>>_pimK0pi0max>> _pimpi0etamax;
  }
  
  ClassDescription<Tau3MesonDefaultDecayer> Tau3MesonDefaultDecayer::initTau3MesonDefaultDecayer;
  // Definition of the static class description member.
  
  void Tau3MesonDefaultDecayer::Init() {
        
    static ClassDocumentation<Tau3MesonDefaultDecayer> documentation
      ("The \\classname{Tau3MesonDefaultDecayer} class is designed to implement "
       "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
       "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
       "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.");

    static ParVector<Tau3MesonDefaultDecayer,double> interfaceF123RhoWgt
      ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &Tau3MesonDefaultDecayer::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);

    static ParVector<Tau3MesonDefaultDecayer,double> interfaceF123KstarWgt
      ("F123KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &Tau3MesonDefaultDecayer::_KstarF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);

    static ParVector<Tau3MesonDefaultDecayer,double> interfaceF5RhoWgt
      ("F5RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &Tau3MesonDefaultDecayer::_rhoF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);

    static ParVector<Tau3MesonDefaultDecayer,double> interfaceF5KstarWgt
      ("F5KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &Tau3MesonDefaultDecayer::_KstarF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);

    static Parameter<Tau3MesonDefaultDecayer,double> interfaceRhoKstarWgt
      ("RhoKstarWgt",
       "The relative weights of the rho and K* in the F5 form factor",
       &Tau3MesonDefaultDecayer::_rhoKstarwgt, -0.2, -1.0e12, 1.0e12,
       false, false, false);
    
    
    static Switch<Tau3MesonDefaultDecayer,bool> interfaceInitializea1
      ("Initializea1",
       "Initialise the calculation of the a_1 running width",
       &Tau3MesonDefaultDecayer::_initializea1, false, false, false);
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



    static Switch<Tau3MesonDefaultDecayer,bool> interfaceRhoParameters
      ("RhoParameters",
       "Use local values of the rho meson masses and widths",
       &Tau3MesonDefaultDecayer::_rhoparameters, true, false, false);
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

    static Switch<Tau3MesonDefaultDecayer,bool> interfaceKstarParameters
      ("KstarParameters",
       "Use local values of the rho meson masses and widths",
       &Tau3MesonDefaultDecayer::_Kstarparameters, true, false, false);
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

    static Switch<Tau3MesonDefaultDecayer,bool> interfacea1Parameters
      ("a1Parameters",
       "Use local values of the rho meson masses and widths",
       &Tau3MesonDefaultDecayer::_a1parameters, true, false, false);
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


    static Switch<Tau3MesonDefaultDecayer,bool> interfaceK1Parameters
      ("K1Parameters",
       "Use local values of the rho meson masses and widths",
       &Tau3MesonDefaultDecayer::_K1parameters, true, false, false);
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


    static ParVector<Tau3MesonDefaultDecayer,double> interfacea1RunningWidth
      ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &Tau3MesonDefaultDecayer::_a1runwidth,
     0, 0, 0, 0, 100000, false, false, true);

    static ParVector<Tau3MesonDefaultDecayer,double> interfacea1RunningQ2
      ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &Tau3MesonDefaultDecayer::_a1runq2,
     0, 0, 0, 0, 100000, false, false, true);


    static Parameter<Tau3MesonDefaultDecayer,Energy> interfaceA1Width
      ("A1Width",
       "The a_1 width if using local values.",
       &Tau3MesonDefaultDecayer::_a1width, MeV, 559.0*MeV, 0*MeV, 10000*MeV,
       false, false, false);
    
    static Parameter<Tau3MesonDefaultDecayer,Energy> interfaceA1Mass
      ("A1Mass",
       "The a_1 mass if using local values.",
       &Tau3MesonDefaultDecayer::_a1mass, MeV, 1251.0*MeV, 0*MeV, 10000*MeV,
       false, false, false);

    static Parameter<Tau3MesonDefaultDecayer,Energy> interfaceK1Width
      ("K1Width",
       "The K_1 width if using local values.",
       &Tau3MesonDefaultDecayer::_K1width, MeV, 559.0*MeV, 0*MeV, 10000*MeV,
       false, false, false);
    
    static Parameter<Tau3MesonDefaultDecayer,Energy> interfaceK1Mass
      ("K1Mass",
       "The K_1 mass if using local values.",
       &Tau3MesonDefaultDecayer::_K1mass, MeV, 1251.0*MeV, 0*MeV, 10000*MeV,
       false, false, false);
    
    
    static ParVector<Tau3MesonDefaultDecayer,double> interfacepimpimpipWeights
      ("pi-pi-pi+Weights",
       "The weights of the channels for the integration of the pi-pi-pi+ channel",
       &Tau3MesonDefaultDecayer::_pimpimpipwgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfacepi0pi0pimWeights
      ("pi0pi0pi+Weights",
       "The weights of the channels for the integration of the pi0pi0pi- channel",
       &Tau3MesonDefaultDecayer::_pi0pi0pimwgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfaceKmpi0K0Weights
      ("K-pi0K0Weights",
       "The weights of the channels for the integration of the K-pi0K0 channel",
       &Tau3MesonDefaultDecayer::_Kmpi0K0wgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfaceKmpimKpWeights
      ("K-pi-K+Weights",
       "The weights of the channels for the integration of the K-pi-K+ channel",
       &Tau3MesonDefaultDecayer::_KmpimKpwgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfacepi0pi0KmWeights
      ("pi0pi0K-Weights",
       "The weights of the channels for the integration of the pi0pi0K- channel",
       &Tau3MesonDefaultDecayer::_pi0pi0Kmwgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfaceK0pimK0Weights
      ("K0pi-K0Weights",
       "The weights of the channels for the integration of the K0pi-K0 channel",
       &Tau3MesonDefaultDecayer::_K0pimK0wgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfaceKmpimpipWeights
      ("K-pi-pi+Weights",
       "The weights of the channels for the integration of the K-pi-pi+ channel",
       &Tau3MesonDefaultDecayer::_Kmpimpipwgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfacepimK0pi0Weights
      ("pi-K0pi0Weights",
       "The weights of the channels for the integration of the pi-K0pi0 channel",
       &Tau3MesonDefaultDecayer::_pimK0pi0wgts,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfacepimpi0etaWeights
      ("pi-pi0etaWeights",
       "The weights of the channels for the integration of the pi-pi0eta channel",
       &Tau3MesonDefaultDecayer::_pimpi0etawgts,
       0, 0, 0, 0, 1.0, false, false, true);

    static Parameter<Tau3MesonDefaultDecayer,double> interfacepimpimpipMax
      ("pi-pi-pi+Max",
       "Maximum weight for the integration of the pi-pi-pi+ channel",
       &Tau3MesonDefaultDecayer::_pimpimpipmax, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfacepi0pi0pimMax
      ("pi0pi0pi-Max",
       "Maximum weight for the integration of the pi0pi0pi- channel",
       &Tau3MesonDefaultDecayer::_pi0pi0pimmax, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfaceKmpimKpMax
      ("K-pi-K+Max",
       "Maximum weight for the integration of the K-pi-K+ channel",
       &Tau3MesonDefaultDecayer::_KmpimKpmax, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfacepK0pimK0Max
      ("K0pi-K0Max",
       "Maximum weight for the integration of the K0pi-K0 channel",
       &Tau3MesonDefaultDecayer::_K0pimK0max, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfaceKmpi0K0Max
      ("K-pi0K0Max",
       "Maximum weight for the integration of the K-pi0K0 channel",
       &Tau3MesonDefaultDecayer::_Kmpi0K0max, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfacepi0pi0KmMax
      ("pi0pi0K-Max",
       "Maximum weight for the integration of the pi0pi0K- channel",
       &Tau3MesonDefaultDecayer::_pi0pi0Kmmax, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfaceKmpimpipMax
      ("K-pi-pi+Max",
       "Maximum weight for the integration of the K-pi-pi+ channel",
       &Tau3MesonDefaultDecayer::_Kmpimpipmax, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfacepimK0pi0Max
      ("pi-K0pi0Max",
       "Maximum weight for the integration of the pi-K0pi0 channel",
       &Tau3MesonDefaultDecayer::_pimK0pi0max, 1.0, 0, 1.0e12,
       false, false, false);
    static Parameter<Tau3MesonDefaultDecayer,double> interfacepimpi0etaMax
      ("pi-pi0etaMax",
       "Maximum weight for the integration of the pi-pi0peta channel",
       &Tau3MesonDefaultDecayer::_pimpi0etamax, 1.0, 0, 1.0e12,
       false, false, false);


    static ParVector<Tau3MesonDefaultDecayer,double> interfacerhoF123masses
      ("rhoF123masses",
       "The masses for the rho resonances if used local values",
       &Tau3MesonDefaultDecayer::_rhoF123masses,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfacerhoF123widths
      ("rhoF123widths",
       "The widths for the rho resonances if used local values",
       &Tau3MesonDefaultDecayer::_rhoF123widths,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfacerhoF5masses
      ("rhoF5masses",
       "The masses for the rho resonances if used local values",
       &Tau3MesonDefaultDecayer::_rhoF5masses,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfacerhoF5widths
      ("rhoF5widths",
       "The widths for the rho resonances if used local values",
       &Tau3MesonDefaultDecayer::_rhoF5widths,
       0, 0, 0, 0, 1.0, false, false, true);

    static ParVector<Tau3MesonDefaultDecayer,double> interfaceKstarF123masses
      ("KstarF123masses",
       "The masses for the Kstar resonances if used local values",
       &Tau3MesonDefaultDecayer::_KstarF123masses,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfaceKstarF123widths
      ("KstarF123widths",
       "The widths for the Kstar resonances if used local values",
       &Tau3MesonDefaultDecayer::_KstarF123widths,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfaceKstarF5masses
      ("KstarF5masses",
       "The masses for the Kstar resonances if used local values",
       &Tau3MesonDefaultDecayer::_KstarF5masses,
       0, 0, 0, 0, 1.0, false, false, true);
    static ParVector<Tau3MesonDefaultDecayer,double> interfaceKstarF5widths
      ("KstarF5widths",
       "The widths for the Kstar resonances if used local values",
       &Tau3MesonDefaultDecayer::_KstarF5widths,
       0, 0, 0, 0, 1.0, false, false, true);
  }

  // modes handled by this class
  bool Tau3MesonDefaultDecayer::acceptMode(int imode) const{return imode>=1&&imode<=9;}

  // calculate the form-factors
  void Tau3MesonDefaultDecayer::calculateFormFactors(const int imode, const int ichan,
						     Energy2 q2, Energy2 s1,
						     Energy2 s2, Energy2 s3,
						     Complex & F1,
						     Complex & F2,
						     Complex & F3,
						     Complex & F4,
						     Complex & F5) const
  {
    // calculate the pi- pi- pi+ factor
    if(imode==1)
      {
	Complex a1fact=a1BreitWigner(q2)*2./3.;
	F3 = 0.;
	F4 = 0.;
	F5 = 0.;
	if(ichan<0)
	  {
	    F1 = a1fact*BrhoF123(s1,-1);
	    F2 =-a1fact*BrhoF123(s2,-1);
	  }
	else if(ichan<3)
	  {
	    F1 = a1fact*BrhoF123(s1,ichan);
	    F2=0.;
	  }
	else if(ichan<6)
	  {
	    F1=0.;
	    F2=-a1fact*BrhoF123(s2,ichan-3);
	  }

      }
    // calculate the pi0 pi0 pi- factor
    else if(imode==2)
      {
	Complex a1fact=a1BreitWigner(q2)*2./3.;
	F3 = 0.;
	F4 = 0.;
	F5 = 0.;
	if(ichan<0)
	  {
	    F1 = a1fact*BrhoF123(s1,-1);
	    F2 =-a1fact*BrhoF123(s2,-1);
	  }
	else if(ichan<9)
	  {
	    F1 = 0.;
	    F2 =-a1fact*BrhoF123(s2,-1);
	  }
	else if(ichan<12)
	  {
	    F1 = 0.;
	    F2 =-a1fact*BrhoF123(s2,ichan-9);
	  }
      }
    // calculate the K- pi - K+ factor
    else if(imode==3)
      {
	Complex a1fact=a1BreitWigner(q2)*sqrt(2.)/3.;
	F3 = 0.;
	F4 = 0.;
	if(ichan<0)
	  {
	    F1 =-a1fact*BKstarF123(s1,-1);
	    F2 = a1fact*BrhoF123(s2,-1);
	    F5 = BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
	  }
	else if((ichan-3)%8==1)
	  {
	    F1 = -a1fact*BKstarF123(s1,(ichan-12)/8);
	    F2 = 0.;
	    F5 = 0.;
	  }
	else if((ichan-3)%8==2)
	  {
	    F1 = 0.;
	    F2 = a1fact*BrhoF123(s2,(ichan-12)/8-1);
	    F5 = 0.;
	  }
	else if((ichan-14)%8>=0)
	  {
	    F1 = 0.;
	    F2 = 0.;
	    F5 = BrhoF5(q2,(ichan-12)/8)*FKrho(s1,s2,(ichan-14)%8)*sqrt(2.);
	  }
      }
    // calculate the K0 pi- k0bar
    else if(imode==4)
      {
	Complex a1fact=a1BreitWigner(q2)*sqrt(2.)/3.;
	F3 = 0.;
	F4 = 0.;
	if(ichan<0)
	  {
	    F1 =-a1fact*BKstarF123(s1,-1);
	    F2 = a1fact*BrhoF123(s2,-1);
	    F5 =-BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
	  }
	else if((ichan-11)%8==1)
	  {
	    F1 = -a1fact*BKstarF123(s1,(ichan-20)/8);
	    F2 = 0.;
	    F5 = 0.;
	  }
	else if((ichan-11)%8==2)
	  {
	    F1 = 0.;
	    F2 = a1fact*BrhoF123(s2,(ichan-20)/8-1);
	    F5 = 0.;
	  }
	else if((ichan-22)%8>=0)
	  {
	    F1 = 0.;
	    F2 = 0.;
	    F5 = -BrhoF5(q2,(ichan-20)/8)*FKrho(s1,s2,(ichan-22)%8)*sqrt(2.);
	  }
      }
    // calculate the K- pi0 k0
    else if(imode==5)
      {
	Complex a1fact=a1BreitWigner(q2);
	F1 = 0.;
	F3 = 0.;
	F4 = 0.;
	F5 = 0.;
	if(ichan<0){F2 =-a1fact*BrhoF123(s2,-1);}
	else{F2 =-a1fact*BrhoF123(s2,ichan-60);}
      }
    // calculate the pi0 pi0 K-
    else if(imode==6)
      {
	Complex K1fact=K1BreitWigner(q2)/6.;
	F3 = 0.;
	F4 = 0.;
	F5 = 0.;
	if(ichan<0)
	  {
	    F1 = K1fact*BKstarF123(s1,-1);
	    F2 =-K1fact*BKstarF123(s2,-1);
	  }
	else if(ichan<66)
	  {
	    F1 = K1fact*BKstarF123(s1,ichan-63);
	    F2 = 0.;
	  }
	else
	  {
	    F1 = 0.;
	    F2 =-K1fact*BKstarF123(s2,ichan-66);
	  }
      }
    // calculate the K- pi- pi+
    else if(imode==7)
      {
	Complex K1fact=K1BreitWigner(q2)*sqrt(2.)/3.;
	F3 = 0.;
	F4 = 0.;
	if(ichan<0)
	  {
	    F1 =-K1fact*BrhoF123(s1,-1);
	    F2 = K1fact*BKstarF123(s2,-1);
	    F5 =-BKstarF123(q2,-1)*FKrho(s2,s1,-1)*sqrt(2.);
	  }
	else if((ichan-68)%8==1)
	  {
	    F1 = -K1fact*BrhoF123(s1,(ichan-69)/8);
	    F2 = 0.;
	    F5 = 0.;
	  }
	else if((ichan-68)%8==2)
	  {
	    F1 = 0.;
	    F2 = K1fact*BKstarF123(s2,(ichan-69)/8);
	    F5 = 0.;
	  }
	else
	  {
	    F1 = 0;
	    F2 = 0;
	    F5 = -BKstarF123(q2,(ichan-69)/8)*FKrho(s2,s1,(ichan-71)%8)*sqrt(2.);
	  }
      }
    // calculate the pi- K0bar pi0
    else if(imode==8)
      {
	Complex K1fact=K1BreitWigner(q2);
	F1 = 0.;
	F3 = 0.;
	F4 = 0.;
	if(ichan<0)
	  {
	    F2 =-K1fact*BrhoF123(s2,-1);
	    F5 =-2.*BKstarF123(q2,-1)*FKrho(s1,s2,-1);
	  }
	else if((ichan-92)%7)
	  {
	    F2 =-K1fact*BrhoF123(s2,(ichan-93)/7);
	    F5 =0.;
	    
	  }
	else 
	  {
	    F2 =0.;
	    F5 =-2.*BKstarF123(q2,(ichan-93)/7)*FKrho(s1,s2,(ichan-94)%7);
	  }
      }
    // calculate the pi- pi0 eta
    else if(imode==9)
      {
	F1 = 0.;
	F2 = 0.;
	F3 = 0.;
	F4 = 0.;
	if(ichan<0)
	  {F5 = BrhoF5(q2,-1)*BrhoF123(s3,-1)*sqrt(2./3.);}
	else
	  {F5 = BrhoF5(q2,(ichan-114)/3)*BrhoF123(s3,(ichan-114)%3)*sqrt(2./3.);}
      }
    // multiply by the prefactors
    if(imode<=5||imode==9)
      {F1*=_cosfact;F2*=_cosfact;F3*=_cosfact;F4*=_cosfact;F5*=_cosfact;}
    else
      {F1*=_sinfact;F2*=_sinfact;F3*=_sinfact;F4*=_sinfact;F5*=_sinfact;}
    F5 =-F5*Complex(0.,1.)/4./pi/pi/_fpi/_fpi;
  }
  // mapping of mode to integration channel
  int Tau3MesonDefaultDecayer::phaseSpaceMode(int in) const {return in-1;}
}

// the functions for the integrands of the a_1 width
namespace Herwig {
  using namespace Genfun;

  FUNCTION_OBJECT_IMP(tau3MesonDefaulta1MatrixElement)

    tau3MesonDefaulta1MatrixElement::tau3MesonDefaulta1MatrixElement(Ptr<Herwig::Tau3MesonDefaultDecayer>::pointer in)
  {_decayer=in;}

  
  tau3MesonDefaulta1MatrixElement::~tau3MesonDefaulta1MatrixElement() {
  }
  
  tau3MesonDefaulta1MatrixElement::tau3MesonDefaulta1MatrixElement(const tau3MesonDefaulta1MatrixElement & right) 
  {  }
  
  unsigned int tau3MesonDefaulta1MatrixElement::dimensionality() const {return 4;}

  double tau3MesonDefaulta1MatrixElement::operator ()(const Argument & a) const 
  {
    return _decayer->a1MatrixElement(a[0],a[1],a[2],a[3]);
  }

}
