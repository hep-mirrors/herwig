// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourPionNovosibirskCurrent class.
//

#include "FourPionNovosibirskCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"


#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FourPionNovosibirskCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;  

FourPionNovosibirskCurrent::~FourPionNovosibirskCurrent() {}

void FourPionNovosibirskCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhomass << _a1mass << _omegamass << _sigmamass << _rhowidth << _a1width
     << _omegawidth << _sigmawidth << _zsigma << _lambda2
     << _initializea1 << _localparameters << _a1runwidth << _a1runq2 << _onedlam2 
     << _a1massolam2 << _psigma << _mpi
     << _aomega << _athreec << _aonec << _bomega << _bthreec << _bonec 
     << _comega << _cthreec <<_conec << _omegaparam << _intwidth << _intmass
     << _mpi2 << _hm2 << _dhdq2m2 << _prho << _rhoD;
}

void FourPionNovosibirskCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhomass >> _a1mass >> _omegamass >> _sigmamass >> _rhowidth >> _a1width
     >> _omegawidth >> _sigmawidth >> _zsigma >> _lambda2
     >> _initializea1 >> _localparameters >> _a1runwidth >> _a1runq2 >> _onedlam2
     >> _a1massolam2 >>_psigma >> _mpi 
     >> _aomega >> _athreec >> _aonec >> _bomega >> _bthreec >> _bonec 
     >> _comega >> _cthreec >>_conec >> _omegaparam >> _intwidth >> _intmass
     >> _mpi2 >> _hm2 >> _dhdq2m2 >> _prho >> _rhoD;
}

ClassDescription<FourPionNovosibirskCurrent> FourPionNovosibirskCurrent::initFourPionNovosibirskCurrent;
// Definition of the static class description member.

void FourPionNovosibirskCurrent::Init() {

  static ClassDocumentation<FourPionNovosibirskCurrent> documentation
    ("The \\classname{FourPionNovosibirskCurrent} class performs the decay"
     " of the tau to four pions using currents based on the the"
     " Novosibirsk e+e- data");

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacerhoMass
    ("rhoMass",
     "The local value of the rho mass",
     &FourPionNovosibirskCurrent::_rhomass, GeV,0.7761*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacea1mass
    ("a1Mass",
     "The local value of the square of the a_1 mass",
     &FourPionNovosibirskCurrent::_a1mass, GeV, 1.2300*GeV, 0.5*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceSigmaMass
    ("sigmaMass",
     "The local value of the sigma mass",
     &FourPionNovosibirskCurrent::_sigmamass, GeV, 0.8*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceOmegaMass
    ("omegaMass",
     "The local value of the omega mass",
     &FourPionNovosibirskCurrent::_omegamass, GeV, 0.7820*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacerhoWidth
    ("rhoWidth",
     "The local value of the rho width",
     &FourPionNovosibirskCurrent::_rhowidth, GeV,0.1445*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacea1width
    ("a1Width",
     "The local value of the square of the a_1 width",
     &FourPionNovosibirskCurrent::_a1width, GeV, 0.45*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceSigmaWidth
    ("sigmaWidth",
     "The local value of the sigma width",
     &FourPionNovosibirskCurrent::_sigmawidth, GeV, 0.8*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceOmegaWidth
    ("omegaWidth",
     "The local value of the omega width",
     &FourPionNovosibirskCurrent::_omegawidth, GeV, 0.00841*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceIntegrationMass
    ("IntegrationMass",
     "Mass of the pseudoresonance used to improve integration effciency",
     &FourPionNovosibirskCurrent::_intmass, GeV, 1.4*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceIntegrationWidth
    ("IntegrationWidth",
     "Width of the pseudoresonance used to improve integration effciency",
     &FourPionNovosibirskCurrent::_intwidth, GeV, 0.5*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,double> interfaceSigmaMagnitude
    ("SigmaMagnitude",
     "magnitude of the relative sigma coupling",
     &FourPionNovosibirskCurrent::_zmag, 1.3998721, 0.0, 10.0,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,double> interfaceSigmaPhase
    ("SigmaPhase",
     "phase of the relative sigma coupling",
     &FourPionNovosibirskCurrent::_zphase, 0.43585036, 0.0, 2.*pi,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of the mass scale squared to use in the form-factor",
     &FourPionNovosibirskCurrent::_lambda2, GeV2, 1.2*GeV2, 0.0001*GeV2, 10.0*GeV2,
     false, false, true);

  static Switch<FourPionNovosibirskCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &FourPionNovosibirskCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use the local values",
     true);
  static SwitchOption interfaceLocalParametersDefault
    (interfaceLocalParameters,
     "Default",
     "Use the values from the particleData objects",
     false);

  static Switch<FourPionNovosibirskCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &FourPionNovosibirskCurrent::_initializea1, false, false, false);
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
  
  static ParVector<FourPionNovosibirskCurrent,double> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &FourPionNovosibirskCurrent::_a1runwidth,
     0, 0, 0, 0, 100000, false, false, true);
  
  static ParVector<FourPionNovosibirskCurrent,double> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &FourPionNovosibirskCurrent::_a1runq2,
     0, 0, 0, 0, 10000000, false, false, true);

}

// complete the construction of the decay mode for integration
bool FourPionNovosibirskCurrent::createMode(int icharge, unsigned int imode,
					    DecayPhaseSpaceModePtr mode,
					    unsigned int iloc,unsigned int ires,
					    DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  // check the charge
  if(abs(icharge)!=3){return false;}
  bool kineallowed=true;
  // check that the modes are kinematical allowed
  Energy min(0.);
  if(imode==0)
    {min=   getParticleData(ParticleID::piplus)->mass()
	+3.*getParticleData(ParticleID::pi0)->mass();}
  else
    {min=3.*getParticleData(ParticleID::piplus)->mass()
	+getParticleData(ParticleID::pi0)->mass();}
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  // intermediates for the channels
  tPDPtr omega,rhop,rhom,rho0,a1m,a10,sigma,rhot;
  omega= getParticleData(ParticleID::omega);
  rho0 = getParticleData(ParticleID::rho0);
  sigma= getParticleData(9000221); 
  a10  = getParticleData(ParticleID::a_10);
  if(icharge==3)
    {
      rhop = getParticleData(ParticleID::rhominus);
      rhom = getParticleData(ParticleID::rhoplus);
      a1m  = getParticleData(ParticleID::a_1plus);
      rhot = getParticleData(30213);
    }
  else if(icharge==-3)
    {
      rhop = getParticleData(ParticleID::rhoplus);
      rhom = getParticleData(ParticleID::rhominus);
      a1m  = getParticleData(ParticleID::a_1minus);
      rhot = getParticleData(-30213);
    }
  else{return false;}
  DecayPhaseSpaceChannelPtr newchannel;
  // the omega channels for the three charged pion mode
  // first  channel two channels with rho0
  if(imode==1)
    {
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(rho0    ,0,0.0, iloc+1,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(rho0    ,0,0.0, iloc,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      // second two channels with rho -
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      // third two channels with rho +
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      //  a_1 channels with rhos
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(rho0    ,0,0.0, iloc+1,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(rho0    ,0,0.0, iloc,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      // neutral a_1 channels with rhos
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      //  a_1 channels with sigmas	
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      // neutral a_1 channels with sigma
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
    }
  else
    {
      // channels with an a1- and a rho -
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+3,iloc);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+2,iloc);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+3,iloc);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+2,iloc);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc);
      newchannel->init();
      mode->addChannel(newchannel);
      // channels with a sigma and a10
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+2,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
      // channels with a1- and sigma
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+2,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+3);
      newchannel->init();
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      newchannel->init();
      mode->addChannel(newchannel);
    }
  // reset the parameters of the dummy resonance used for integration
  mode->resetIntermediate(rhot,_intmass,_intwidth);
  // reset the parameters of the resonances if using local values
  if(_localparameters)
    {
      mode->resetIntermediate(rhom,_rhomass,_rhowidth);
      mode->resetIntermediate(rhop,_rhomass,_rhowidth);
      mode->resetIntermediate(rho0,_rhomass,_rhowidth);
      mode->resetIntermediate(omega,_omegamass,_omegawidth);
      mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
      mode->resetIntermediate(rhot,_intmass,_intwidth);
    }
  // return if successful
  return kineallowed;
}

// the particles produced by the current
  PDVector FourPionNovosibirskCurrent::particles(int icharge, unsigned int imode,int iq,
						 int ia)
{
  PDVector output(4);
  if(icharge==3)
    {
      if(imode==1)
	{
	  output[0]=getParticleData(ParticleID::piplus);
	  output[1]=getParticleData(ParticleID::piplus);
	  output[2]=getParticleData(ParticleID::piminus);
	}
      else
	{
	  output[0]=getParticleData(ParticleID::piplus);
	  output[1]=getParticleData(ParticleID::pi0);
	  output[2]=getParticleData(ParticleID::pi0);
	}
    }
  else if(icharge==-3)
    {
      if(imode==1)
	{
	  output[0]=getParticleData(ParticleID::piminus);
	  output[1]=getParticleData(ParticleID::piminus);
	  output[2]=getParticleData(ParticleID::piplus);
	}
      else
	{
	  output[0]=getParticleData(ParticleID::piminus);
	  output[1]=getParticleData(ParticleID::pi0);
	  output[2]=getParticleData(ParticleID::pi0);
	}
    }
  output[3]=getParticleData(ParticleID::pi0);
  return output;
}

 
// the hadronic currents    
vector<LorentzPolarizationVector> 
FourPionNovosibirskCurrent::current(bool vertex, const int imode, const int ichan,
				    Energy & scale,const ParticleVector & decay) const
{
  LorentzPolarizationVector output;
  double fact(1.);
  // construct the spininfo objects if needed
  if(vertex)
    {
      for(unsigned int ix=0;ix<decay.size();++ix)
	{
	  SpinPtr temp=new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true));
	  decay[ix]->spinInfo(temp);
	}
    }
  // the momenta of the particles
  Lorentz5Momentum q1(decay[3]->momentum()),q2(decay[1]->momentum()),
    q3(decay[0]->momentum()),q4(decay[2]->momentum());
  Lorentz5Momentum Q=q1+q2+q3+q4;Q.rescaleMass();
  scale = Q.mass();
  // decide which decay mode
  // three charged pions
  if(imode==1)
    {
      // momenta of the particles
      LorentzPolarizationVector veca1rho,vecomega,veca1sig;
      if(ichan<0)
	{
	  // a_1 rho current
	  veca1rho = 
	    t1(q1,q2,q3,q4)+t1(q3,q2,q1,q4)+t1(q1,q3,q2,q4)
	    +t1(q3,q1,q2,q4)+t1(q4,q3,q1,q2)+t1(q4,q1,q3,q2);
	  // a_1 sigma current
	  veca1sig = 
	    t2(q4,q3,q1,q2)+t2(q4,q1,q3,q2)
	    -t2(q1,q4,q3,q2)-t2(q3,q4,q1,q2);
	  // omega current
	  vecomega = 
	    t3(q1,q2,q3,q4)+t3(q3,q2,q1,q4)-t3(q1,q3,q2,q4)
	    -t3(q3,q1,q2,q4)-t3(q1,q4,q3,q2)-t3(q3,q4,q1,q2);
	}
      else if(ichan== 0){vecomega = t3(q1,q4,q3,q2);}
      else if(ichan== 1){vecomega = t3(q3,q4,q1,q2);}
      else if(ichan== 2){vecomega = t3(q1,q2,q3,q4);}
      else if(ichan== 3){vecomega = t3(q3,q2,q1,q4);}
      else if(ichan== 4){vecomega = t3(q1,q3,q2,q4);}
      else if(ichan== 5){vecomega = t3(q3,q1,q2,q4);}
      else if(ichan== 6){veca1rho = t1(q4,q1,q3,q2);}
      else if(ichan== 7){veca1rho = t1(q4,q3,q1,q2);}
      else if(ichan== 8){veca1rho = t1(q1,q2,q3,q4);}
      else if(ichan== 9){veca1rho = t1(q3,q2,q1,q4);}
      else if(ichan==10){veca1rho = t1(q1,q3,q2,q4);}
      else if(ichan==11){veca1rho = t1(q3,q1,q2,q4);}
      else if(ichan==12){veca1sig = t2(q4,q1,q3,q2);}
      else if(ichan==13){veca1sig = t2(q4,q3,q1,q2);}
      else if(ichan==14){veca1sig = t2(q1,q4,q3,q2);}
      else if(ichan==15){veca1sig = t2(q3,q4,q1,q2);}
      // final manipulations
      veca1rho+=veca1sig;
      veca1rho*=gFunction(Q.mass2(),1);;
      vecomega*=gFunction(Q.mass2(),2);
      output=vecomega+veca1rho;
      // this is 1/sqrt(2) for identical particles
      fact*=0.7071067811865476;
    }
  else if(imode==0)
    {
      // momenta of the particles
      LorentzPolarizationVector veca1rho,veca1sig;
      if(ichan<0)
	{
	  // a_1 rho current
	  veca1rho=
	    t1(q2,q3,q1,q4)+t1(q2,q4,q1,q3)+t1(q3,q2,q1,q4)
	    +t1(q3,q4,q1,q2)+t1(q4,q2,q1,q3)+t1(q4,q3,q1,q2);
	  // a_1 sigma current
	  veca1sig=
	    t2(q2,q1,q3,q4)+t2(q3,q1,q2,q4)+t2(q4,q1,q3,q2)
	    -t2(q1,q2,q3,q4)-t2(q1,q3,q2,q4)-t2(q1,q4,q3,q2);
	}
      else if(ichan== 0){veca1rho = t1(q2,q3,q1,q4);}
      else if(ichan== 1){veca1rho = t1(q2,q4,q1,q3);}
      else if(ichan== 2){veca1rho = t1(q3,q2,q1,q4);}
      else if(ichan== 3){veca1rho = t1(q3,q4,q1,q2);}
      else if(ichan== 4){veca1rho = t1(q4,q2,q1,q3);}
      else if(ichan== 5){veca1rho = t1(q4,q3,q1,q2);}
      else if(ichan== 6){veca1sig = t2(q2,q1,q3,q4);}
      else if(ichan== 7){veca1sig = t2(q3,q1,q2,q4);}
      else if(ichan== 8){veca1sig = t2(q4,q1,q3,q2);}
      else if(ichan== 9){veca1sig = t2(q1,q2,q3,q4);}
      else if(ichan==10){veca1sig = t2(q1,q3,q2,q4);}
      else if(ichan==11){veca1sig = t2(q1,q4,q3,q2);}
      // add them up 
      output=veca1rho+veca1sig;
      output*=gFunction(Q.mass2(),0);
      // this is sqrt(1/3!) for identical particles
      fact*=0.40824829046386631;
    }     
  else
    {throw DecayIntegratorError() << "Unknown decay mode in the " 
				  << "FourPionNovosibirskCurrent::"
				  << "hadronCurrent()" << Exception::abortnow;}  
  output*=fact*Q.mass2();
  vector<LorentzPolarizationVector> temp; 
  temp.push_back(output);
  return temp;
}

bool FourPionNovosibirskCurrent::accept(vector<int> id)
{
  bool allowed=false;
  // check four products
  if(id.size()!=4){return false;}
  int npiminus=0,npiplus=0,npi0=0;
  for(unsigned int ix=0;ix<id.size();++ix)
    {
      if(id[ix]==ParticleID:: piplus){++npiplus;}
      else if(id[ix]==ParticleID::piminus){++npiminus;}
      else if(id[ix]==ParticleID::pi0){++npi0;}
    }
  if(npiminus==2&&npiplus==1&&npi0==1){allowed=true;}
  else if(npiminus==1&&npi0==3){allowed=true;}
  else if(npiplus==2&&npiminus==1&&npi0==1){allowed=true;}
  else if(npiplus==1&&npi0==3){allowed=true;}
  return allowed;
}

// the decay mode
unsigned int FourPionNovosibirskCurrent::decayMode(vector<int> idout)
{
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix)
    {if(abs(idout[ix])==ParticleID::piplus){++npi;}}
  if(npi==3){return 1;}
  return 0;
}


// output the information for the database
void FourPionNovosibirskCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/FourPionNovosibirskCurrent " << fullName() << " \n";
  output << "set " << fullName() << ":rhoMass "    << _rhomass/GeV << "\n";
  output << "set " << fullName() << ":a1Mass  "    << _a1mass/GeV  << "\n";
  output << "set " << fullName() << ":sigmaMass  " << _sigmamass/GeV  << "\n";
  output << "set " << fullName() << ":omegaMass  " << _omegamass/GeV  << "\n";
  output << "set " << fullName() << ":rhoWidth "    << _rhowidth/GeV << "\n";
  output << "set " << fullName() << ":a1Width  "    << _a1width/GeV  << "\n";
  output << "set " << fullName() << ":sigmaWidth  " << _sigmawidth/GeV  << "\n";
  output << "set " << fullName() << ":omegaWidth  " << _omegawidth/GeV  << "\n";
  output << "set " << fullName() << ":IntegrationMass "  << _intmass/GeV  << "\n";
  output << "set " << fullName() << ":IntegrationWidth " << _intwidth/GeV  << "\n";
  output << "set " << fullName() << ":SigmaMagnitude "  <<  _zmag << "\n";
  output << "set " << fullName() << ":SigmaPhase " << _zphase  << "\n";
  output << "set " << fullName() << ":Lambda2 "  <<  _lambda2/GeV2 << "\n";
  output << "set " << fullName() << ":LocalParameters " <<  _localparameters << "\n";
  output << "set " << fullName() << ":Initializea1 " <<  _initializea1 << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix)
    {output << "insert " << fullName() << ":a1RunningWidth " << ix 
	    << " " << _a1runwidth[ix] << "\n";}
  for(unsigned int ix=0;ix<_a1runq2.size();++ix)
    {output << "insert " << fullName() << ":a1RunningQ2 " << ix 
	    << " " << _a1runq2[ix] << "\n";}
}
 
} 

// the functions for the integrands of the a_1 width
namespace Herwig {
using namespace Genfun;

FUNCTION_OBJECT_IMP(FourPionDefaultMatrixElement)

FourPionDefaultMatrixElement::FourPionDefaultMatrixElement(Ptr<Herwig::FourPionNovosibirskCurrent>::pointer in)
  {_decayer=in;}

FourPionDefaultMatrixElement::~FourPionDefaultMatrixElement() {}
  
FourPionDefaultMatrixElement::FourPionDefaultMatrixElement(const FourPionDefaultMatrixElement & right) 
  {}
  
unsigned int FourPionDefaultMatrixElement::dimensionality() const {return 7;}

double FourPionDefaultMatrixElement::operator ()(const Argument & a) const 
{
  return _decayer->a1MatrixElement(a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
}
  
}

