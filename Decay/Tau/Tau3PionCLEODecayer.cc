// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau3PionCLEODecayer class.
//

#include "Tau3PionCLEODecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Tau3PionCLEODecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig{
using namespace ThePEG;
  
Tau3PionCLEODecayer::~Tau3PionCLEODecayer() {}
  

void Tau3PionCLEODecayer::persistentOutput(PersistentOStream & os) const {
  os << _rhomass << _rhowidth << _prhocc << _prhoc0 << _f2mass << _f2width << _pf2cc 
     << _pf200 << _f0mass << _f0width << _pf0cc << _pf000 << _sigmamass << _sigmawidth
     << _psigmacc << _psigma00 << _mpi0 << _mpic << _fpi << _fact 
     << _rhomagP << _rhophaseP 
     << _rhocoupP << _rhomagD << _rhophaseD << _rhocoupD <<_f2mag << _f2phase << _f2coup 
     << _f0mag << _f0phase << _f0coup << _sigmamag << _sigmaphase << _sigmacoup 
     << _localparameters << _onechan << _threechan << _onewgts << _threewgts << _onemax  
     << _threemax <<_a1mass<< _a1width <<  _a1runwidth << _a1runq2 <<  _initializea1
     << _mKstar << _mK << _gammk;;
}
  
void Tau3PionCLEODecayer::persistentInput(PersistentIStream & is, int) {
  is >> _rhomass >> _rhowidth >> _prhocc >> _prhoc0 >> _f2mass >> _f2width >> _pf2cc 
     >> _pf200 >> _f0mass >> _f0width >> _pf0cc >> _pf000 >> _sigmamass >> _sigmawidth
     >> _psigmacc >> _psigma00 >> _mpi0 >> _mpic >> _fpi >> _fact 
     >> _rhomagP >> _rhophaseP 
     >> _rhocoupP >> _rhomagD >> _rhophaseD >> _rhocoupD>>_f2mag >> _f2phase >> _f2coup 
     >> _f0mag >> _f0phase >> _f0coup >> _sigmamag >> _sigmaphase >> _sigmacoup 
     >> _localparameters >> _onechan >> _threechan  >> _onewgts >> _threewgts >> _onemax 
     >> _threemax>> _a1mass >> _a1width >>  _a1runwidth >> _a1runq2 >>  _initializea1
     >> _mKstar >> _mK >> _gammk;
}
  
ClassDescription<Tau3PionCLEODecayer> Tau3PionCLEODecayer::initTau3PionCLEODecayer;
  // Definition of the static class description member.
  
void Tau3PionCLEODecayer::Init() {
  
  static ClassDocumentation<Tau3PionCLEODecayer> documentation
    ("The \\classname{Tau3PionCLEODecayer} class performs the decay of the"
     " tau to three pions using the currents from CLEO");
  
  static ParVector<Tau3PionCLEODecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &Tau3PionCLEODecayer::_rhomass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<Tau3PionCLEODecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &Tau3PionCLEODecayer::_rhowidth,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &Tau3PionCLEODecayer::_f2mass, GeV, 1.275*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &Tau3PionCLEODecayer::_f2width, GeV, 0.185*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacef_0Mass
    ("f_0Mass",
     "The mass of the f_0 meson",
     &Tau3PionCLEODecayer::_f0mass, GeV, 1.186*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacef_0Width
    ("f_0Width",
     "The width of the f_0 meson",
     &Tau3PionCLEODecayer::_f0width, GeV, 0.350*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma meson",
     &Tau3PionCLEODecayer::_sigmamass, GeV, 0.860*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma meson",
     &Tau3PionCLEODecayer::_sigmawidth, GeV, 0.880*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacea1Mass
    ("a1Mass",
     "The mass of the a_1 meson",
     &Tau3PionCLEODecayer::_a1mass, GeV, 1.331*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfacea1Width
    ("a1Width",
     "The width of the a_1 meson",
     &Tau3PionCLEODecayer::_a1width, GeV, 0.814*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfaceKaonMass
    ("KaonMass",
     "The mass of the kaon",
     &Tau3PionCLEODecayer::_mK, GeV, 0.496*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfaceKStatMass
    ("KStarMass",
     "The mass of the k* meson",
     &Tau3PionCLEODecayer::_mKstar, GeV, 0.894*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfaceKaonCoupling
    ("KaonCoupling",
     "The relative coupling for the kaon in the a_1 running width",
     &Tau3PionCLEODecayer::_gammk, 3.32, 0.0, 10.0,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant",
     &Tau3PionCLEODecayer::_fpi, MeV, 130.7*MeV/sqrt(2.), 0.0*MeV, 10.0*MeV,
     false, false, true);

  static ParVector<Tau3PionCLEODecayer,double> interfacerhomagP
    ("RhoPWaveMagnitude",
     "The magnitude of the couplings for the p-wave rho currents",
     &Tau3PionCLEODecayer::_rhomagP,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<Tau3PionCLEODecayer,double> interfacerhophaseP
    ("RhoPWavePhase",
     "The phase of the couplings for the p-wave rho currents",
     &Tau3PionCLEODecayer::_rhophaseP,
     0, 0, 0, 0, 2.*pi, false, false, true);

  static ParVector<Tau3PionCLEODecayer,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &Tau3PionCLEODecayer::_rhomagD,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<Tau3PionCLEODecayer,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &Tau3PionCLEODecayer::_rhophaseD,
     0, 0, 0, 0, 2.*pi, false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &Tau3PionCLEODecayer::_f0phase, 0.54*pi, 0.0, 2.*pi,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &Tau3PionCLEODecayer::_f2phase, 0.56*pi, 0.0, 2.*pi,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &Tau3PionCLEODecayer::_sigmaphase, 0.23*pi, 0.0, 2.*pi,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfacef0Magnitude
    ("f0Magnitude",
     "The magnitude of the f_0 scalar current",
     &Tau3PionCLEODecayer::_f0mag, 0.77, 0.0, 10,
     false, false, true);


  static Parameter<Tau3PionCLEODecayer,InvEnergy2> interfacef2Magnitude
    ("f2Magnitude",
     "The magnitude of the f_2 tensor current",
     &Tau3PionCLEODecayer::_f2mag, 1./GeV2, 0.71/GeV2, 0./GeV2, 10./GeV2,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfacesigmaMagnitude
    ("sigmaMagnitude",
     "The magnitude of the sigma scalar current",
     &Tau3PionCLEODecayer::_sigmamag, 2.1, 0.0, 10,
     false, false, true);

  static Switch<Tau3PionCLEODecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &Tau3PionCLEODecayer::_localparameters, true, false, false);
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

  static ParVector<Tau3PionCLEODecayer,bool> interfaceonechan
    ("OneChargedChannels",
     "The channels to use for the integration of the decay a_1^+->pi+pi0pi0",
     &Tau3PionCLEODecayer::_onechan,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<Tau3PionCLEODecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &Tau3PionCLEODecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<Tau3PionCLEODecayer,bool> interfacethreechan
    ("ThreeChargedChannels",
     "The channels to use for the integration of the decay a_1^+->pi+pi+pi-",
     &Tau3PionCLEODecayer::_threechan,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<Tau3PionCLEODecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &Tau3PionCLEODecayer::_threewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &Tau3PionCLEODecayer::_onemax,1.23756E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<Tau3PionCLEODecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &Tau3PionCLEODecayer::_threemax, 1.38754E3, 0.0, 10000.0,
     false, false, true);

  static ParVector<Tau3PionCLEODecayer,double> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &Tau3PionCLEODecayer::_a1runwidth,
     0, 0, 0, 0, 100000, false, false, true);
  
  static ParVector<Tau3PionCLEODecayer,double> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &Tau3PionCLEODecayer::_a1runq2,
     0, 0, 0, 0, 100000, false, false, true);

  static Switch<Tau3PionCLEODecayer,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &Tau3PionCLEODecayer::_initializea1, false, false, false);
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
bool Tau3PionCLEODecayer::acceptMode(int imode) const{return imode>=1&&imode<=2;}
  
// calculate the form-factors
void Tau3PionCLEODecayer::calculateFormFactors(const int imode, const int ichan,
					       Energy2 q2, Energy2 s1,
					       Energy2 s2, Energy2 s3,
					       Complex & F1,Complex & F2,
					       Complex & F3,Complex & F4,
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
  
  
void Tau3PionCLEODecayer::CLEOFormFactor(int imode,int ichan,
					 Energy2 q2,Energy2 s1, Energy2 s2, Energy2 s3,
					 Complex & F1, Complex & F2, Complex & F3) const
{
  if(imode==1)
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
  else if(imode==2)
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
      else if(ichan%2==0&&ichan>=12&&ichan<=16)
	{
	  unsigned int ires=(ichan-12)/2;
	  if(ires<_rhocoupP.size()){F1+=_rhocoupP[ires]*rhos1bw[ires];}
	  Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
	  if(ires<_rhocoupD.size())
	    {F2+=Dfact2*_rhocoupD[ires]*rhos1bw[ires];
	    F3+=_rhocoupD[ires]*Dfact2*rhos1bw[ires];}
	}
      else if(ichan%2==1&&ichan>=3&&ichan<=17)
	{
	  unsigned int ires=(ichan-13)/2;
	  if(ires<_rhocoupP.size()){F2+=_rhocoupP[ires]*rhos2bw[ires];}
	  Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
	  if(ires<_rhocoupD.size())
	    {F1+=Dfact1*_rhocoupD[ires]*rhos2bw[ires];
	    F3-=_rhocoupD[ires]*Dfact1*rhos2bw[ires];}
	}
      else if(ichan==18){F1+=2./3.*_sigmacoup*sigbw;F2+=2./3.*_sigmacoup*sigbw;}
      else if(ichan==19)
	{
	  Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
	  F1+=Dfact3;F2+=Dfact3;
	  F3-=0.5*_f2coup*(s1-s2)*f2bw;
	  
	}
      else if(ichan==20){F1+=2./3.*_f0coup*f0bw;F2+=2./3.*_f0coup*f0bw;}
    }
  else
    {throw DecayIntegratorError() << "Tau3PionCLEODecayer Unknown Decay" 
				  << Exception::abortnow;}
  // identical particle factors
  double fact=0.70710678;
  F1*=fact;F2*=fact;F3*=fact;
}

// mapping of mode to integration channel
int Tau3PionCLEODecayer::phaseSpaceMode(int in) const 
{
  if(in==2){return 0;}
  return in;
}


}

// the functions for the integrands of the a_1 width
namespace Herwig {
using namespace Genfun;

FUNCTION_OBJECT_IMP(Tau3PionCLEOa1MatrixElement)

Tau3PionCLEOa1MatrixElement::
Tau3PionCLEOa1MatrixElement(int mode,Ptr<Herwig::Tau3PionCLEODecayer>::pointer in)
{_mode=mode,_decayer=in;}
  
  
Tau3PionCLEOa1MatrixElement::~Tau3PionCLEOa1MatrixElement() {
}
  
Tau3PionCLEOa1MatrixElement::
Tau3PionCLEOa1MatrixElement(const Tau3PionCLEOa1MatrixElement & right) {  }
  
unsigned int Tau3PionCLEOa1MatrixElement::dimensionality() const {return 4;}

double Tau3PionCLEOa1MatrixElement::operator ()(const Argument & a) const 
{return _decayer->a1MatrixElement(_mode,a[0],a[1],a[2],a[3]);}

}
