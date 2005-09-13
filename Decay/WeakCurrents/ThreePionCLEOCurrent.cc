// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreePionCLEOCurrent class.
//

#include "ThreePionCLEOCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreePionCLEOCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

ThreePionCLEOCurrent::~ThreePionCLEOCurrent() {}

void ThreePionCLEOCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhomass << _rhowidth << _prhocc << _prhoc0 << _f2mass << _f2width << _pf2cc 
     << _pf200 << _f0mass << _f0width << _pf0cc << _pf000 << _sigmamass << _sigmawidth
     << _psigmacc << _psigma00 << _mpi0 << _mpic << _fpi << _fact 
     << _rhomagP << _rhophaseP 
     << _rhocoupP << _rhomagD << _rhophaseD << _rhocoupD <<_f2mag << _f2phase << _f2coup 
     << _f0mag << _f0phase << _f0coup << _sigmamag << _sigmaphase << _sigmacoup 
     << _localparameters 
     <<_a1mass<< _a1width <<  _a1runwidth << _a1runq2 <<  _initializea1
     << _mKstar << _mK << _gammk;
}

void ThreePionCLEOCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhomass >> _rhowidth >> _prhocc >> _prhoc0 >> _f2mass >> _f2width >> _pf2cc 
     >> _pf200 >> _f0mass >> _f0width >> _pf0cc >> _pf000 >> _sigmamass >> _sigmawidth
     >> _psigmacc >> _psigma00 >> _mpi0 >> _mpic >> _fpi >> _fact 
     >> _rhomagP >> _rhophaseP 
     >> _rhocoupP >> _rhomagD >> _rhophaseD >> _rhocoupD>>_f2mag >> _f2phase >> _f2coup 
     >> _f0mag >> _f0phase >> _f0coup >> _sigmamag >> _sigmaphase >> _sigmacoup 
     >> _localparameters 
     >> _a1mass >> _a1width >>  _a1runwidth >> _a1runq2 >>  _initializea1
     >> _mKstar >> _mK >> _gammk;
}

ClassDescription<ThreePionCLEOCurrent> ThreePionCLEOCurrent::initThreePionCLEOCurrent;
// Definition of the static class description member.

void ThreePionCLEOCurrent::Init() {

  static ClassDocumentation<ThreePionCLEOCurrent> documentation
    ("The \\classname{ThreePionCLEOCurrent} class performs the decay of the"
     " tau to three pions using the currents from CLEO");
  
  static ParVector<ThreePionCLEOCurrent,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &ThreePionCLEOCurrent::_rhomass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ThreePionCLEOCurrent,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &ThreePionCLEOCurrent::_rhowidth,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &ThreePionCLEOCurrent::_f2mass, GeV, 1.275*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &ThreePionCLEOCurrent::_f2width, GeV, 0.185*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_0Mass
    ("f_0Mass",
     "The mass of the f_0 meson",
     &ThreePionCLEOCurrent::_f0mass, GeV, 1.186*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_0Width
    ("f_0Width",
     "The width of the f_0 meson",
     &ThreePionCLEOCurrent::_f0width, GeV, 0.350*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma meson",
     &ThreePionCLEOCurrent::_sigmamass, GeV, 0.860*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma meson",
     &ThreePionCLEOCurrent::_sigmawidth, GeV, 0.880*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacea1Mass
    ("a1Mass",
     "The mass of the a_1 meson",
     &ThreePionCLEOCurrent::_a1mass, GeV, 1.331*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacea1Width
    ("a1Width",
     "The width of the a_1 meson",
     &ThreePionCLEOCurrent::_a1width, GeV, 0.814*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceKaonMass
    ("KaonMass",
     "The mass of the kaon",
     &ThreePionCLEOCurrent::_mK, GeV, 0.496*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceKStarMass
    ("KStarMass",
     "The mass of the k* meson",
     &ThreePionCLEOCurrent::_mKstar, GeV, 0.894*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfaceKaonCoupling
    ("KaonCoupling",
     "The relative coupling for the kaon in the a_1 running width",
     &ThreePionCLEOCurrent::_gammk, 3.32, 0.0, 10.0,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant",
     &ThreePionCLEOCurrent::_fpi, MeV, 130.7*MeV/sqrt(2.), 0.0*MeV, 500.0*MeV,
     false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhomagP
    ("RhoPWaveMagnitude",
     "The magnitude of the couplings for the p-wave rho currents",
     &ThreePionCLEOCurrent::_rhomagP,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhophaseP
    ("RhoPWavePhase",
     "The phase of the couplings for the p-wave rho currents",
     &ThreePionCLEOCurrent::_rhophaseP,
     0, 0, 0, -2.*pi, 2.*pi, false, false, true);

  static ParVector<ThreePionCLEOCurrent,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &ThreePionCLEOCurrent::_rhomagD,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &ThreePionCLEOCurrent::_rhophaseD,
     0, 0, 0, -2.*pi, 2.*pi, false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &ThreePionCLEOCurrent::_f0phase, 0.54*pi, -2.*pi, 2.*pi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &ThreePionCLEOCurrent::_f2phase, 0.56*pi,-2.*pi, 2.*pi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &ThreePionCLEOCurrent::_sigmaphase, 0.23*pi, -2.*pi, 2.*pi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef0Magnitude
    ("f0Magnitude",
     "The magnitude of the f_0 scalar current",
     &ThreePionCLEOCurrent::_f0mag, 0.77, 0.0, 10,
     false, false, true);


  static Parameter<ThreePionCLEOCurrent,InvEnergy2> interfacef2Magnitude
    ("f2Magnitude",
     "The magnitude of the f_2 tensor current",
     &ThreePionCLEOCurrent::_f2mag, 1./GeV2, 0.71/GeV2, 0./GeV2, 10./GeV2,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacesigmaMagnitude
    ("sigmaMagnitude",
     "The magnitude of the sigma scalar current",
     &ThreePionCLEOCurrent::_sigmamag, 2.1, 0.0, 10,
     false, false, true);

  static Switch<ThreePionCLEOCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &ThreePionCLEOCurrent::_localparameters, true, false, false);
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

  static ParVector<ThreePionCLEOCurrent,double> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &ThreePionCLEOCurrent::_a1runwidth,
     0, 0, 0, 0, 10000000, false, false, true);
  
  static ParVector<ThreePionCLEOCurrent,double> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &ThreePionCLEOCurrent::_a1runq2,
     0, 0, 0, 0, 10000000, false, false, true);

  static Switch<ThreePionCLEOCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &ThreePionCLEOCurrent::_initializea1, false, false, false);
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

}

// modes handled by this class
bool ThreePionCLEOCurrent::acceptMode(int imode) const{return imode>=0&&imode<=1;}

// calculate the form-factors  
void ThreePionCLEOCurrent::
 calculateFormFactors(const int ichan, const int imode,
		      Energy2 q2, Energy2 s1, Energy2 s2, Energy2 s3,
		      Complex & F1, Complex & F2, Complex & F3, Complex & F4,
		      Complex & F5) const
{
  F1=0.;F2=0.;F3=0.;F4=0.;F5=0.;
  // calculate the form factors without the a_1 piece
  CLEOFormFactor(imode,ichan,q2,s1,s2,s3,F1,F2,F3);
  // change sign of the f_2 term
  F2=-F2;
  // multiply by the a_1 factor
  Complex a1fact=a1BreitWigner(q2)*_fact;
  F1*=a1fact;F2*=a1fact;F3*=a1fact;
}

void ThreePionCLEOCurrent::CLEOFormFactor(int imode,int ichan,
					  Energy2 q2,Energy2 s1, Energy2 s2, Energy2 s3,
					  Complex & F1, Complex & F2, Complex & F3) const
{
  if(imode==0)
    {
      // compute the breit wigners we need
      Complex rhos1bw[3],rhos2bw[3],f0bws1,sigbws1,f2bws1,f0bws2,sigbws2,f2bws2;
      for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix)
	{
	  rhos1bw[ix]=rhoBreitWigner(ix,s1,0);
	  rhos2bw[ix]=rhoBreitWigner(ix,s2,0);
	}
      f0bws1  =f0BreitWigner(s1,0);
      sigbws1 =sigmaBreitWigner(s1,0);
      f2bws1  =f2BreitWigner(s1,0);
      f0bws2  =f0BreitWigner(s2,0);
      sigbws2 =sigmaBreitWigner(s2,0);
      f2bws2  =f2BreitWigner(s2,0);
      if(ichan<0)
	{
	  // the p-wave rho terms
	  for(unsigned int ix=0;ix<_rhocoupP.size();++ix)
	    {F1-=_rhocoupP[ix]*rhos1bw[ix];F2-=_rhocoupP[ix]*rhos2bw[ix];}
	  // the D-wave rho terms
	  Energy2 Dfact1=1./3.*(s1-s3);
	  Energy2 Dfact2=1./3.*(s2-s3);
	  for(unsigned int ix=0;ix<_rhocoupD.size();++ix)
	    {
	      F1-=Dfact1*_rhocoupD[ix]*rhos2bw[ix];
	      F2-=Dfact2*_rhocoupD[ix]*rhos1bw[ix];
	      F3-=_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]);
	    }
	  // the scalar terms
	  F1-=2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
	  F2-=2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
	  F3+=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	    +2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
	  // the tensor terms
	  Complex sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
	  Complex sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
	  F1+=_f2coup*(0.5*(s3-s2)*f2bws1-sfact2);
	  F2+=_f2coup*(0.5*(s3-s1)*f2bws2-sfact1);
	  F3+=_f2coup*(-sfact1+sfact2);
	}
      else if(ichan%2==0&&ichan<=4)
	{
	  unsigned int ires=ichan/2;
	  Energy2 Dfact2=1./3.*(s2-s3);
	  if(ires<_rhocoupP.size()){F1-=_rhocoupP[ires]*rhos1bw[ires];}
	  if(ires<_rhocoupD.size())
	    {F2-=Dfact2*_rhocoupD[ires]*rhos1bw[ires];
	    F3-=_rhocoupD[ires]*Dfact2*rhos1bw[ires];}
	}
      else if(ichan%2==1&&ichan<=5)
	{
	  unsigned int ires=(ichan-1)/2;
	  Energy2 Dfact1=1./3.*(s1-s3);
	  if(ires<_rhocoupP.size()){F2-=_rhocoupP[ires]*rhos2bw[ires];}
	  if(ires<_rhocoupD.size())
	    {F1-=Dfact1*_rhocoupD[ires]*rhos2bw[ires];
	    F3+=_rhocoupD[ires]*Dfact1*rhos2bw[ires];}
	}
      else if(ichan==6){F2-=2./3.*_sigmacoup*sigbws1;F3-=2./3.*_sigmacoup*sigbws1;}
      else if(ichan==7){F1-=2./3.*_sigmacoup*sigbws2;F3+=2./3.*_sigmacoup*sigbws2;}
      else if(ichan==8)
	{
	  Complex sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
	  F1+=_f2coup*0.5*(s3-s2)*f2bws1;F2-=_f2coup*sfact1;F3-=_f2coup*sfact1;
	}
      else if(ichan==9)
	{
	  Complex sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
	  F1-=_f2coup*sfact2;F2+=_f2coup*0.5*(s3-s1)*f2bws2;F3+=_f2coup*sfact2;
	}
      else if(ichan==10){F2-=2./3.*_f0coup*f0bws1;F3-=2./3.*_f0coup*f0bws1;}
      else if(ichan==11){F1-=2./3.*_f0coup*f0bws2;F3+=2./3.*_f0coup*f0bws2;}
    }
  // calculate the pi0 pi0 pi+ factor
  else if(imode==1)
    {
      // compute the breit wigners we need
      Complex rhos1bw[3],rhos2bw[3],f0bw,sigbw,f2bw;
      for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix)
	{
	  rhos1bw[ix]=rhoBreitWigner(ix,s1,1);
	  rhos2bw[ix]=rhoBreitWigner(ix,s2,1);
	}
      f0bw  =f0BreitWigner(s3,1);
      sigbw =sigmaBreitWigner(s3,1);
      f2bw  =f2BreitWigner(s3,1);
      if(ichan<0)
	{
	  // the p-wave rho terms
	  for(unsigned int ix=0;ix<_rhocoupP.size();++ix)
	    {
	      F1+=_rhocoupP[ix]*rhos1bw[ix];
	      F2+=_rhocoupP[ix]*rhos2bw[ix];
	    }
	  // the D-wave rho terms
	  Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
	  Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
	  for(unsigned int ix=0;ix<_rhocoupD.size();++ix)
	    {
	      F1+=Dfact1*_rhocoupD[ix]*rhos2bw[ix];
	      F2+=Dfact2*_rhocoupD[ix]*rhos1bw[ix];
	      F3+=_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]);
	    }
	  // the scalar terms
	  Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
	  F1+=scalar;F2+=scalar;
	  // the tensor terms
	  Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
	  F1+=Dfact3;F2+=Dfact3;
	  F3-=0.5*_f2coup*(s1-s2)*f2bw;
	}
      else if(ichan%2==0&&ichan<=4)
	{
	  unsigned int ires=ichan/2;
	  if(ires<_rhocoupP.size()){F1+=_rhocoupP[ires]*rhos1bw[ires];}
	  Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
	  if(ires<_rhocoupD.size())
	    {F2+=Dfact2*_rhocoupD[ires]*rhos1bw[ires];
	    F3+=_rhocoupD[ires]*Dfact2*rhos1bw[ires];}
	}
      else if(ichan%2==1&&ichan<=5)
	{
	  unsigned int ires=(ichan-1)/2;
	  if(ires<_rhocoupP.size()){F2+=_rhocoupP[ires]*rhos2bw[ires];}
	  Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
	  if(ires<_rhocoupD.size())
	    {F1+=Dfact1*_rhocoupD[ires]*rhos2bw[ires];
	    F3-=_rhocoupD[ires]*Dfact1*rhos2bw[ires];}
	}
      else if(ichan==6){F1+=2./3.*_sigmacoup*sigbw;F2+=2./3.*_sigmacoup*sigbw;}
      else if(ichan==7)
	{
	  Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
	  F1+=Dfact3;F2+=Dfact3;
	  F3-=0.5*_f2coup*(s1-s2)*f2bw;
	  
	}
      else if(ichan==8){F1+=2./3.*_f0coup*f0bw;F2+=2./3.*_f0coup*f0bw;}
    }
  else
    {throw DecayIntegratorError() << "ThreePionCLEOCurrent Unknown Decay" << imode
				  << Exception::abortnow;}
  // identical particle factors
  double fact=0.70710678;
  F1*=fact;F2*=fact;F3*=fact;
}

// complete the construction of the decay mode for integration
bool ThreePionCLEOCurrent::createMode(int icharge, unsigned int imode,
				      DecayPhaseSpaceModePtr mode,
				      unsigned int iloc,unsigned int ires,
				      DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  bool kineallowed=true;
  if(!acceptMode(imode)){return false;}
  int iq(0),ia(0);
  PDVector extpart=particles(1,imode,iq,ia);
  Energy min(0.);
  for(unsigned int ix=0;ix<extpart.size();++ix)
    {min+=extpart[0]->mass()-extpart[0]->widthLoCut();}
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  // pointers to the particles we need
  tPDPtr a1m = getParticleData(ParticleID::a_1minus);
  // the different rho resonances
  tPDPtr rhom[3];
  if(icharge==-3)
    {
      rhom[0] = getParticleData(-213);
      rhom[1] = getParticleData(-100213);
      rhom[2] = getParticleData(-30213);
    }
  else if(icharge==3)
    {
      rhom[0] = getParticleData(213);
      rhom[1] = getParticleData(100213);
      rhom[2] = getParticleData(30213);
    }
  else
    {return false;}
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),
		    getParticleData(30113)};
  // the sigma
  tPDPtr sigma = getParticleData(9000221);
  // the f_2
  tPDPtr f2=getParticleData(225);
  // the f_0
  tPDPtr f0=getParticleData(10221);
  // set up the integration channels
  DecayPhaseSpaceChannelPtr newchannel;
  if(imode==0)
    {
      for(unsigned int ix=0;ix<3;++ix)
	{
	  // the neutral rho channels
	  // first channel
	  newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rho0[ix],0,0.0,iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  // interchanged channel
	  newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rho0[ix],0,0.0,iloc,iloc+2);
	  mode->addChannel(newchannel);      
	}
      // the sigma channels
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(sigma,0,0.0,iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // interchanged channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(sigma,0,0.0,iloc,iloc+2);
      mode->addChannel(newchannel);
      // the f_2 channels
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(f2,0,0.0,iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // interchanged channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(f2,0,0.0,iloc,iloc+2);
      mode->addChannel(newchannel);
      // the f_0 channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(f0,0,0.0,iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // interchanged channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(f0,0,0.0,iloc,iloc+2);
      mode->addChannel(newchannel);
    }
  else
    {
      for(unsigned int ix=0;ix<3;++ix)
	{
	  // first rho+ channel
	  newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rhom[ix],0,0.0,iloc+2,iloc+1);
	  mode->addChannel(newchannel);
	  // second rho+ channel
	  newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhom[ix],0,0.0,iloc+2,iloc);
	  mode->addChannel(newchannel);
	}
      // the sigma channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(sigma,0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
      //  the f_2  channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(f2,0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
      // the f_0 channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(f0,0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
    }
  return kineallowed;
  if(_localparameters)
    {
      for(unsigned int iy=0;iy<_rhomass.size();++iy)
	{
	  mode->resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
	  mode->resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
	}
      mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
      mode->resetIntermediate(f2,_f2mass,_f2width);
      mode->resetIntermediate(f0,_f0mass,_f0width);
      mode->resetIntermediate(a1m,_a1mass,_a1width);
    }
}

PDVector ThreePionCLEOCurrent::particles(int icharge, unsigned int imode,int iq,int ia)
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
  else
    {extpart.resize(0);}
  // conjugate the particles if needed
  if(icharge==3&&extpart.size()>0)
    {for(unsigned int ix=0;ix<3;++ix)
	{if(extpart[ix]->CC()){extpart[ix]=extpart[ix]->CC();}}}
  // return the answer
  return extpart;
}

void ThreePionCLEOCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::ThreePionCLEOCurrent " << fullName() << " \n";}
  for(unsigned int ix=0;ix<_rhomass.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoMasses " << ix 
		<< " " << _rhomass[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":RhoMasses " << ix 
		<< " " << _rhomass[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoWidths " << ix 
		<< " " << _rhowidth[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":RhoWidths " << ix 
		<< " " << _rhowidth[ix] << "\n";}
    }
  output << "set " << fullName() << ":f_2Mass " << _f2mass/GeV << "\n";
  output << "set " << fullName() << ":f_2Width " << _f2width/GeV << "\n";
  output << "set " << fullName() << ":f_0Mass " << _f0mass/GeV << "\n";
  output << "set " << fullName() << ":f_0Width " << _f0width/GeV << "\n";
  output << "set " << fullName() << ":sigmaMass " << _sigmamass/GeV << "\n";
  output << "set " << fullName() << ":sigmaWidth " << _sigmawidth/GeV << "\n";
  output << "set " << fullName() << ":a1Mass " << _a1mass/GeV << "\n";
  output << "set " << fullName() << ":a1Width " <<_a1width /GeV << "\n";
  output << "set " << fullName() << ":KaonMass " << _mK/GeV << "\n";
  output << "set " << fullName() << ":KStarMass " << _mKstar/GeV << "\n";
  output << "set " << fullName() << ":KaonCoupling " << _gammk << "\n";
  output << "set " << fullName() << ":Fpi " << _fpi/MeV << "\n";
  for(unsigned int ix=0;ix<_rhomagP.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoPWaveMagnitude " << ix 
		<< " " << _rhomagP[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":RhoPWaveMagnitude " << ix 
		<< " " << _rhomagP[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_rhophaseP.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoPWavePhase " << ix 
		<< " " << _rhophaseP[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":RhoPWavePhase " << ix 
		<< " " << _rhophaseP[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_rhomagD.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoDWaveMagnitude " << ix 
		<< " " << _rhomagD[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":RhoDWaveMagnitude " << ix 
		<< " " << _rhomagD[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_rhophaseD.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoDWavePhase " << ix 
		<< " " << _rhophaseD[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":RhoDWavePhase " << ix 
		<< " " << _rhophaseD[ix] << "\n";}
    }
  output << "set " << fullName() << ":f0Phase " << _f0phase << "\n";
  output << "set " << fullName() << ":f2Phase " <<_f2phase  << "\n";
  output << "set " << fullName() << ":sigmaPhase " <<_sigmaphase  << "\n";
  output << "set " << fullName() << ":f0Magnitude " << _f0mag << "\n";
  output << "set " << fullName() << ":f2Magnitude " << _f2mag*GeV2 << "\n";
  output << "set " << fullName() << ":sigmaMagnitude " <<_sigmamag  << "\n";
  output << "set " << fullName() << ":LocalParameters " << _localparameters << "\n";
  output << "set " << fullName() << ":Initializea1 " <<_initializea1  << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix)
    {
      if(ix<200)
	{output << "set    " << fullName() << ":a1RunningWidth " << ix 
		<< " " << _a1runwidth[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":a1RunningWidth " << ix 
		<< " " << _a1runwidth[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix)
    {
      if(ix<200)
	{output << "set    " << fullName() << ":a1RunningQ2 " << ix 
		<< " " << _a1runq2[ix] << "\n";}
      else
	{output << "insert " << fullName() << ":a1RunningQ2 " << ix 
		<< " " << _a1runq2[ix] << "\n";}
    }
  ThreeMesonCurrentBase::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}

// the functions for the integrands of the a_1 width
namespace Herwig {
using namespace Genfun;
FUNCTION_OBJECT_IMP(ThreePionCLEOa1MatrixElement)

ThreePionCLEOa1MatrixElement::
ThreePionCLEOa1MatrixElement(int mode,Ptr<Herwig::ThreePionCLEOCurrent>::pointer in)
{_mode=mode,_decayer=in;}
  
  
ThreePionCLEOa1MatrixElement::~ThreePionCLEOa1MatrixElement() {
}
  
ThreePionCLEOa1MatrixElement::
ThreePionCLEOa1MatrixElement(const ThreePionCLEOa1MatrixElement & right) {  }
  
unsigned int ThreePionCLEOa1MatrixElement::dimensionality() const {return 7;}

double ThreePionCLEOa1MatrixElement::operator ()(const Argument & a) const 
 {return _decayer->a1MatrixElement(_mode,a[0],a[1],a[2],a[3],a[4],a[5],a[6]);}

}
