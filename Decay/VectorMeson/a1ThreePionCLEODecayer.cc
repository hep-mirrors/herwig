// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1ThreePionCLEODecayer class.
//

#include "a1ThreePionCLEODecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "a1ThreePionCLEODecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"

namespace Herwig{
using namespace ThePEG;
using namespace Helicity;
using namespace ThePEG::Helicity;

a1ThreePionCLEODecayer::~a1ThreePionCLEODecayer() {}

int a1ThreePionCLEODecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  int id(dm.parent()->id());
  if(dm.products().size()!=3){return imode;}
  // check the pions
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int idtemp,npi0(0),npiplus(0),npiminus(0);
  for( ; pit!=pend;++pit)
    {
      idtemp=(**pit).id();
      if(idtemp==ParticleID::piplus){++npiplus;}
      else if(idtemp==ParticleID::piminus){++npiminus;}
      else if(idtemp==ParticleID::pi0){++npi0;}
    }
  // a_1+ decay modes
  if(id==ParticleID::a_1plus)
    {
      cc=false;
      if(npiplus==1&&npi0==2){imode=1;}
      else if(npiplus==2&&npiminus==1){imode=3;}
    }
  // a_1- modes
  else if(id==ParticleID::a_1minus)
    {
      cc=true;
      if(npiminus==1&&npi0==2){imode=1;}
      else if(npiminus==2&&npiplus==1){imode=3;}
    }
  // a_0 modes
  else if(id==ParticleID::a_10)
    {
      cc=false;
      if(npiminus==1&&npiplus==1&&npi0==1){imode=0;}
      else if(npi0==3){imode=2;}
    }
  return imode;
}
  
void a1ThreePionCLEODecayer::persistentOutput(PersistentOStream & os) const {
  os << _rhomass << _rhowidth << _prhocc << _prhoc0 << _f2mass << _f2width << _pf2cc 
     << _pf200 << _f0mass << _f0width << _pf0cc << _pf000 << _sigmamass << _sigmawidth
     << _psigmacc << _psigma00 << _mpi0 << _mpic << _coupling << _rhomagP << _rhophaseP 
     << _rhocoupP << _rhomagD << _rhophaseD << _rhocoupD <<_f2mag << _f2phase << _f2coup 
     << _f0mag << _f0phase << _f0coup << _sigmamag << _sigmaphase << _sigmacoup 
     << _localparameters << _zerowgts << _onewgts << _twowgts << _threewgts 
     << _zeromax << _onemax << _twomax << _threemax;
}
  
void a1ThreePionCLEODecayer::persistentInput(PersistentIStream & is, int) {
  is >> _rhomass >> _rhowidth >> _prhocc >> _prhoc0 >> _f2mass >> _f2width >> _pf2cc
     >> _pf200 >> _f0mass >> _f0width >> _pf0cc >> _pf000 >> _sigmamass >> _sigmawidth 
     >> _psigmacc >> _psigma00 >> _mpi0 >> _mpic >> _coupling >> _rhomagP >> _rhophaseP
     >> _rhocoupP >> _rhomagD >> _rhophaseD >> _rhocoupD>>_f2mag >> _f2phase >> _f2coup
     >> _f0mag >> _f0phase >> _f0coup >> _sigmamag >> _sigmaphase >> _sigmacoup
     >> _localparameters >> _zerowgts >> _onewgts >> _twowgts >> _threewgts
     >> _zeromax >> _onemax >> _twomax >> _threemax;
}
  
ClassDescription<a1ThreePionCLEODecayer> a1ThreePionCLEODecayer::inita1ThreePionCLEODecayer;
// Definition of the static class description member.
  
void a1ThreePionCLEODecayer::Init() {
  
  static ClassDocumentation<a1ThreePionCLEODecayer> documentation
    ("The a1ThreePionCLEODecayer class performs the decay of the "
     "a_1 to three pions using the model of CLEO");

  static ParVector<a1ThreePionCLEODecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1ThreePionCLEODecayer::_rhomass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1ThreePionCLEODecayer::_rhowidth,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &a1ThreePionCLEODecayer::_f2mass, GeV, 1.275*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &a1ThreePionCLEODecayer::_f2width, GeV, 0.185*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_0Mass
    ("f_0Mass",
     "The mass of the f_0 meson",
     &a1ThreePionCLEODecayer::_f0mass, GeV, 1.186*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_0Width
    ("f_0Width",
     "The width of the f_0 meson",
     &a1ThreePionCLEODecayer::_f0width, GeV, 0.350*GeV, 0.0*GeV, 1.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma meson",
     &a1ThreePionCLEODecayer::_sigmamass, GeV, 0.860*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma meson",
     &a1ThreePionCLEODecayer::_sigmawidth, GeV, 0.880*GeV, 0.0*GeV, 2.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1ThreePionCLEODecayer::_coupling, 1./GeV, 51.197275/GeV, -0./GeV, 1000./GeV,
     false, false, false);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhomagP
    ("RhoPWaveMagnitude",
     "The magnitude of the couplings for the p-wave rho currents",
     &a1ThreePionCLEODecayer::_rhomagP,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhophaseP
    ("RhoPWavePhase",
     "The phase of the couplings for the p-wave rho currents",
     &a1ThreePionCLEODecayer::_rhophaseP,
     0, 0, 0, -pi, pi, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhomagD,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhophaseD,
     0, 0, 0, -pi, pi, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &a1ThreePionCLEODecayer::_f0phase, 0.54*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &a1ThreePionCLEODecayer::_f2phase, 0.56*pi, -pi, pi,
     false, false, Interface::limited);


  static Parameter<a1ThreePionCLEODecayer,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &a1ThreePionCLEODecayer::_sigmaphase, 0.23*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef0Magnitude
    ("f0Magnitude",
     "The magnitude of the f_0 scalar current",
     &a1ThreePionCLEODecayer::_f0mag, 0.77, 0.0, 10,
     false, false, true);


  static Parameter<a1ThreePionCLEODecayer,InvEnergy2> interfacef2Magnitude
    ("f2Magnitude",
     "The magnitude of the f_2 tensor current",
     &a1ThreePionCLEODecayer::_f2mag, 1./GeV2, 0.71/GeV2, 0./GeV2, 10./GeV2,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacesigmaMagnitude
    ("sigmaMagnitude",
     "The magnitude of the sigma scalar current",
     &a1ThreePionCLEODecayer::_sigmamag, 2.1, 0.0, 10,
     false, false, true);

  static Switch<a1ThreePionCLEODecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1ThreePionCLEODecayer::_localparameters, true, false, false);
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

  static ParVector<a1ThreePionCLEODecayer,double> interfacezerowgts
    ("AllNeutralWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zerowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1ThreePionCLEODecayer::_threewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceZeroMax
    ("ZeroMax",
     "The maximum weight for the integration fo the channel a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zeromax, 0.0716349E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onemax,1.23756E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twomax,2.43819E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1ThreePionCLEODecayer::_threemax, 1.38754E3, 0.0, 10000.0,
     false, false, true);
}

// hadronic current
vector<LorentzPolarizationVector> 
a1ThreePionCLEODecayer::decayCurrent(const bool vertex,
				     const int ichan,const Particle & inpart,
				     const ParticleVector &outpart) const
{
  // momentum of the incoming particle
  Lorentz5Momentum Q=inpart.momentum();
  Energy2 q2=Q.mass2();
  // construct the spin info objects if needed
  if(vertex)
    {for(unsigned int ix=0;ix<outpart.size();++ix)
	{outpart[ix]->spinInfo(new_ptr(ScalarSpinInfo(outpart[ix]->momentum(),true)));}}
  // identify the mesons
  unsigned int iloc[3]={iloc[0]=0,iloc[1]=1,iloc[2]=2};
  if(imode()==1){iloc[0]=1;iloc[1]=2;iloc[2]=0;}
  // calculate the invariants and form factors
  Complex F1=0.,F2=0.,F3=0.;
  Lorentz5Momentum ps1=outpart[iloc[1]]->momentum()+outpart[iloc[2]]->momentum();
  Lorentz5Momentum ps2=outpart[iloc[0]]->momentum()+outpart[iloc[2]]->momentum();
  Lorentz5Momentum ps3=outpart[iloc[0]]->momentum()+outpart[iloc[1]]->momentum();
  ps1.rescaleMass();ps2.rescaleMass();ps3.rescaleMass();
  Energy s1=ps1.mass2(),s2=ps2.mass2(),s3=ps3.mass2();
  formFactors(imode(),ichan,q2,s1,s2,s3,F1,F2,F3);
  // use the form-factors to compute the current
  LorentzPolarizationVector output=
     F1*outpart[iloc[1]]->momentum()-F1*outpart[iloc[2]]->momentum()
    -F2*outpart[iloc[2]]->momentum()+F2*outpart[iloc[0]]->momentum()
    +F3*outpart[iloc[0]]->momentum()-F3*outpart[iloc[1]]->momentum();
  vector<LorentzPolarizationVector> temp;temp.push_back(output);
  return temp;
}

// matrix element for the running a_1 width
double a1ThreePionCLEODecayer::threeBodyMatrixElement(int iopt,Energy2 q2, Energy2 s3,
						      Energy2 s2,Energy2 s1,
						      Energy m1,Energy m2,Energy m3)
{
  Energy m12=m1*m1,m22=m2*m2,m32=m3*m3;
  // calculate the form factors
  Complex F1=0.,F2=0.,F3=0.;
  formFactors(iopt,-1,q2,s1,s2,s3,F1,F2,F3);
  // analytic calculation of the matrix element
  double dot1=( F1*conj(F1)*(2.*m22+2.*m32-s1)+F2*conj(F2)*(2.*m12+2.*m32-s2)
	       +F3*conj(F3)*(2.*m12+2.*m22-s3)-F1*conj(F2)*( s1+s2-s3-4.*m32)
	       +F1*conj(F3)*( s1-s2+s3-4.*m22)-F2*conj(F3)*(-s1+s2+s3-4.*m12)).real();
  Complex dot2 = 0.5*(F1*(s3-m32-s2+m22)-F2*(s1-m12-s3+m32)+F3*(s2-m22-s1+m12));
  return (-dot1+(dot2*conj(dot2)).real()/q2)/3.;
}

WidthCalculatorBasePtr 
a1ThreePionCLEODecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit){if(abs((**pit).id())==ParticleID::piplus){++ncharged;}}
  // integrator to perform the integral
  vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
  vector<int> intype;intype.push_back(2);intype.push_back(3);
  Energy mrho=getParticleData(ParticleID::rhoplus)->mass();
  Energy wrho=getParticleData(ParticleID::rhoplus)->width();
  vector<double> inmass;inmass.push_back(mrho);inmass.push_back(mrho);
  vector<double> inwidth;inwidth.push_back(wrho);inwidth.push_back(wrho);
  Energy m[3];
  if(ncharged==0)     {m[0]=_mpi0;m[1]=_mpi0;m[2]=_mpi0;}
  else if(ncharged==1){m[0]=_mpi0;m[1]=_mpi0;m[2]=_mpic;}
  else if(ncharged==2){m[0]=_mpic;m[1]=_mpic;m[2]=_mpi0;}
  else                {m[0]=_mpic;m[1]=_mpic;m[2]=_mpic;}
  return new_ptr(ThreeBodyAllOnCalculator(inweights,intype,inmass,inwidth,
					  const_ptr_cast<tDecayIntegratorPtr>(this),
					  ncharged,m[0],m[1],m[2]));
}

// calculate the form factos
void a1ThreePionCLEODecayer::formFactors(int iopt,int ichan,
					 Energy2 q2,Energy2 s1,Energy2 s2,
					 Energy2 s3,Complex & F1,
					 Complex& F2,Complex& F3) const
{
  double fact=_coupling;
  // a_1^0 pi0 pi0 pi0 mode
  if(iopt==0)
    {
      fact*=0.40824829;
      // compute the breit wigners we need
      Complex sigbws1 = sigmaBreitWigner(s1,1);
      Complex sigbws2 = sigmaBreitWigner(s2,1);
      Complex sigbws3 = sigmaBreitWigner(s3,1);
      Complex f0bws1  = f0BreitWigner(s1,1);
      Complex f0bws2  = f0BreitWigner(s2,1);
      Complex f0bws3  = f0BreitWigner(s3,1);
      Complex f2bws1  = f2BreitWigner(s1,1);
      Complex f2bws2  = f2BreitWigner(s2,1);
      Complex f2bws3  = f2BreitWigner(s3,1);
      if(ichan<0)
	{
	  // the scalar terms
	  F1=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	    -2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
	  F2=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	    -2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
	  F3=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	    +2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
	  // the tensor terms
	  Complex Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
	  Complex Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
	  Complex Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
	  F1+=_f2coup*( 0.5*(s3-s2)*f2bws1-Dfact2+Dfact3);
	  F2+=_f2coup*( 0.5*(s3-s1)*f2bws2-Dfact1+Dfact3);
	  F3+=_f2coup*(-0.5*(s1-s2)*f2bws3-Dfact1+Dfact2);
	}
      else if(ichan==0){F2=-2./3.*_sigmacoup*sigbws1;F3=-2./3.*_sigmacoup*sigbws1;}
      else if(ichan==1){F1=-2./3.*_sigmacoup*sigbws2;F3=+2./3.*_sigmacoup*sigbws2;}
      else if(ichan==2){F1= 2./3.*_sigmacoup*sigbws3;F2= 2./3.*_sigmacoup*sigbws3;}
      else if(ichan==3)
	{
	  Complex Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
	  F1+=_f2coup*0.5*(s3-s2)*f2bws1;F2-=_f2coup*Dfact1; F3-=_f2coup*Dfact1;
	}
      else if(ichan==4)
	{
	  Complex Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
	  F2+=_f2coup*0.5*(s3-s1)*f2bws2;F1-=_f2coup*Dfact2;F3+=_f2coup*Dfact2;
	}
      else if(ichan==5)
	{
	  Complex Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
	  F3+=-_f2coup*0.5*(s1-s2)*f2bws3;F1+=_f2coup*Dfact3;F2+=_f2coup*Dfact3;
	}
      else if(ichan==6){F2=-2./3.*_f0coup*f0bws1;F3=-2./3.*_f0coup*f0bws1;}
      else if(ichan==7){F1=-2./3.*_f0coup*f0bws2;F3=+2./3.*_f0coup*f0bws2;}
      else if(ichan==8){F1= 2./3.*_f0coup*f0bws3;F2= 2./3.*_f0coup*f0bws3;}
    }
  // a_1^+ -> pi0 pi0 pi+
  else if(iopt==1)
    {
      fact*=0.70710678;
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
  // a_1^0 ->pi+pi-pi0
  else if(iopt==2)
    {
      // compute the breit wigners we need
      Complex rhos1bw[3],rhos2bw[3],f0bw,sigbw,f2bw;
      for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix)
	{
	  rhos1bw[ix]=rhoBreitWigner(ix,s1,1);
	  rhos2bw[ix]=rhoBreitWigner(ix,s2,1);
	}
      f0bw  =f0BreitWigner(s3,0);
      sigbw =sigmaBreitWigner(s3,0);
      f2bw  =f2BreitWigner(s3,0);
      if(ichan<0)
	{
	  // the p-wave rho terms
	  for(unsigned int ix=0;ix<_rhocoupP.size();++ix)
	    {
	      F1+=_rhocoupP[ix]*rhos1bw[ix];
	      F2+=_rhocoupP[ix]*rhos2bw[ix];
	    }
	  // the D-wave rho terms
	  Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
	  Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
	  for(unsigned int ix=0;ix<_rhocoupD.size();++ix)
	    {
	      F1+=Dfact1*_rhocoupD[ix]*rhos2bw[ix];
	      F2+=Dfact2*_rhocoupD[ix]*rhos1bw[ix];
	      F3+=_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]);
	    }
	  // the scalar terms
	  Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
	  F1+=scalar;
	  F2+=scalar;
	  // the tensor terms
	  Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
	  F1+=Dfact3;F2+=Dfact3;
	  F3-=0.5*_f2coup*(s1-s2)*f2bw;
	}
      else if(ichan%2==0&&ichan<=4)
	{
	  unsigned int ires=ichan/2;
	  if(ires<_rhocoupP.size()){F1+=_rhocoupP[ires]*rhos1bw[ires];}
	  Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
	  if(ires<_rhocoupD.size())
	    {F2+=Dfact2*_rhocoupD[ires]*rhos1bw[ires];
	    F3+=_rhocoupD[ires]*Dfact2*rhos1bw[ires];}
	}
      else if(ichan%2==1&&ichan<=5)
	{
	  unsigned int ires=(ichan-1)/2;
	  if(ires<_rhocoupP.size()){F2+=_rhocoupP[ires]*rhos2bw[ires];}
	  Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
	  if(ires<_rhocoupD.size())
	    {F1+=Dfact1*_rhocoupD[ires]*rhos2bw[ires];
	    F3-=_rhocoupD[ires]*-Dfact1*rhos2bw[ires];}
	}
      else if(ichan==6){F1+=2./3.*_sigmacoup*sigbw;F2+=2./3.*_sigmacoup*sigbw;}
      else if(ichan==7)
	{
	  Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
	  F1+=Dfact3;F2+=Dfact3;
	  F3-=0.5*_f2coup*(s1-s2)*f2bw;
	}
      else if(ichan==8){F1+=2./3.*_f0coup*f0bw;F2+=2./3.*_f0coup*f0bw;}
    }
  // a_1^+ -> pi+ pi+ pi- mode
  else
    {
      fact*=0.70710678;
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
  F1*=fact;F2*=fact;F3*=fact;
} 
// output the setup information for the particle database
void a1ThreePionCLEODecayer::dataBaseOutput(ofstream & output,
					    bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  VectorMesonDecayerBase::dataBaseOutput(output,false);
  // masses and widths of the intermediate particles
  output << "set " << fullName() << ":f_2Mass "    << _f2mass/GeV     << "\n";
  output << "set " << fullName() << ":f_2Width "   << _f2width/GeV    << "\n";
  output << "set " << fullName() << ":f_0Mass "    << _f0mass/GeV     << "\n";
  output << "set " << fullName() << ":f_0Width "   << _f0width/GeV    << "\n";
  output << "set " << fullName() << ":sigmaMass "  << _sigmamass/GeV  << "\n";
  output << "set " << fullName() << ":sigmaWidth " << _sigmawidth/GeV << "\n";
  for(unsigned int ix=0;ix<_rhomass.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoMasses " << ix << " " 
		<< _rhomass[ix] << endl;}
      else
	{output << "insert " << fullName() << ":RhoMasses " << ix << " " 
		<< _rhomass[ix] << endl;}
    }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoWidths " << ix << " " 
		<< _rhowidth[ix] << endl;}
      else
	{output << "insert " << fullName() << ":RhoWidths " << ix << " " 
		<< _rhowidth[ix] << endl;}
    }
  // couplings and phases for different channels
  output << "set " << fullName() << ":f0Phase " << _f0phase << "\n";
  output << "set " << fullName() << ":f2Phase " << _f2phase<< "\n";
  output << "set " << fullName() << ":sigmaPhase " << _sigmaphase<< "\n";
  output << "set " << fullName() << ":f0Magnitude " << _f0mag<< "\n";
  output << "set " << fullName() << ":f2Magnitude " << _f2mag*GeV << "\n";
  output << "set " << fullName() << ":sigmaMagnitude " << _sigmamag/GeV << "\n";
  output << "set " << fullName() << ":Coupling " << _coupling*GeV << "\n";
  for(unsigned int ix=0;ix<_rhomagP.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoPWaveMagnitude " << ix << " " 
		<< _rhomagP[ix] << endl;}
      else
	{output << "insert " << fullName() << ":RhoPWaveMagnitude " << ix << " " 
		<< _rhomagP[ix] << endl;}
    }
  for(unsigned int ix=0;ix<_rhophaseP.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoPWavePhase " << ix << " " 
		<< _rhophaseP[ix] << endl;}
      else
	{output << "insert " << fullName() << ":RhoPWavePhase " << ix << " " 
		<< _rhophaseP[ix] << endl;}
    }  
  for(unsigned int ix=0;ix<_rhomagD.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoDWaveMagnitude " << ix << " " 
		<< _rhomagD[ix] << endl;}
      else
	{output << "insert " << fullName() << ":RhoDWaveMagnitude " << ix << " " 
		<< _rhomagD[ix] << endl;}
    }
  for(unsigned int ix=0;ix<_rhophaseD.size();++ix)
    {
      if(ix<2)
	{output << "set    " << fullName() << ":RhoDWavePhase " << ix << " " 
		<< _rhophaseD[ix] << endl;}
      else
	{output << "insert " << fullName() << ":RhoDWavePhase " << ix << " " 
		<< _rhophaseD[ix] << endl;}
    }
  // use local values of the masses etc.
  output << "set " << fullName() << ":LocalParameters " << _localparameters << "\n";
  // integration weights for the different channels
  for(unsigned int ix=0;ix<_zerowgts.size();++ix)
    {
      output << "set " << fullName() << ":AllNeutralWeights " 
	     << ix << " " << _zerowgts[ix] << "\n";
    }
  for(unsigned int ix=0;ix<_onewgts.size();++ix)
    {
      output << "set " << fullName() << ":OneChargedWeights " 
	     << ix << " " << _onewgts[ix] << "\n";
    }
  for(unsigned int ix=0;ix<_twowgts.size();++ix)
    {
      output << "set " << fullName() << ":TwoChargedWeights " 
	     << ix << " " << _twowgts[ix] << "\n";
    }
  for(unsigned int ix=0;ix<_threewgts.size();++ix)
    {
      output << "set " << fullName() << ":ThreeChargedWeights " 
	     << ix << " " << _threewgts[ix] << "\n";
    }
  // maximum weights for the different  channels
  output << "set " << fullName() << ":ZeroMax "  << _zeromax  << "\n";
  output << "set " << fullName() << ":OneMax "   << _onemax   << "\n";
  output << "set " << fullName() << ":TwoMax "   << _twomax   << "\n";
  output << "set " << fullName() << ":ThreeMax " << _threemax << "\n";
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
} 
