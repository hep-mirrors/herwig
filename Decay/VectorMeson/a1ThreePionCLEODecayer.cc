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

namespace Herwig{
using namespace ThePEG;
using namespace Helicity;
using namespace ThePEG::Helicity;

a1ThreePionCLEODecayer::~a1ThreePionCLEODecayer() {}
  
bool a1ThreePionCLEODecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  int id=dm.parent()->id();
  if((id==ParticleID::a_1plus||id==ParticleID::a_1minus||id==ParticleID::a_10)&&
     dm.products().size()==3)
    { 
      ParticleMSet::const_iterator pit  = dm.products().begin();
      ParticleMSet::const_iterator pend = dm.products().end();
      int idtemp,npi0=0,npiplus=0,npiminus=0;
      for( ; pit!=pend;++pit)
	{
	  idtemp=(**pit).id();
	  if(idtemp==ParticleID::piplus){++npiplus;}
	  else if(idtemp==ParticleID::piminus){++npiminus;}
	  else if(idtemp==ParticleID::pi0){++npi0;}
	}
      // a_1+ decay modes
      if(id==ParticleID::a_1plus&&
	 ((npiplus==2&&npiminus==1)||(npiplus==1&&npi0==2)))
	{allowed=true;}
      else if(id==ParticleID::a_1minus&&
	      ((npiminus==2&&npiplus==1)||(npiminus==1&&npi0==2)))
	{allowed=true;}
      else if(id==ParticleID::a_10&&
	      (npiminus==1&&npiplus==1&&npi0==1)||npi0==3)
	{allowed=true;}   
    }
  return allowed;
}
  
ParticleVector a1ThreePionCLEODecayer::decay(const DecayMode & dm,
					     const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  // work out which mode we are doing
  int id=parent.id(),nplus=0;
  if(id==ParticleID::a_10||id==ParticleID::a_1plus){id=ParticleID::piplus;}
  else{id=ParticleID::piminus;}
  for(unsigned int ix=0;ix<children.size();++ix)
    {if(children[ix]->id()==id){++nplus;}}
  int imode=0;
  if(nplus==0){imode=0;}
  else if(nplus==2){imode=3;}
  else if(nplus==1&&parent.id()==ParticleID::a_10){imode=2;}
  else{imode=1;}
  // perform the decay
  generate(true,imode,parent,children);
  return children;
}
  
void a1ThreePionCLEODecayer::persistentOutput(PersistentOStream & os) const {
  os << _rhomass << _rhowidth << _prhocc << _prhoc0 << _f2mass << _f2width << _pf2cc 
     << _pf200 << _f0mass << _f0width << _pf0cc << _pf000 << _sigmamass << _sigmawidth
     << _psigmacc << _psigma00 << _mpi0 << _mpic << _coupling << _rhomagP << _rhophaseP 
     << _rhocoupP << _rhomagD << _rhophaseD << _rhocoupD <<_f2mag << _f2phase << _f2coup 
     << _f0mag << _f0phase << _f0coup << _sigmamag << _sigmaphase << _sigmacoup 
     << _localparameters << _zerochan << _onechan << _twochan << _threechan << _zerowgts 
     << _onewgts << _twowgts << _threewgts << _zeromax << _onemax << _twomax 
     << _threemax;
}
  
void a1ThreePionCLEODecayer::persistentInput(PersistentIStream & is, int) {
  is >> _rhomass >> _rhowidth >> _prhocc >> _prhoc0 >> _f2mass >> _f2width >> _pf2cc 
     >> _pf200 >> _f0mass >> _f0width >> _pf0cc >> _pf000 >> _sigmamass >> _sigmawidth
     >> _psigmacc >> _psigma00 >> _mpi0 >> _mpic >> _coupling >> _rhomagP >> _rhophaseP 
     >> _rhocoupP >> _rhomagD >> _rhophaseD >> _rhocoupD>>_f2mag >> _f2phase >> _f2coup 
     >> _f0mag >> _f0phase >> _f0coup >> _sigmamag >> _sigmaphase >> _sigmacoup 
     >> _localparameters >> _zerochan >> _onechan >> _twochan >> _threechan >> _zerowgts 
     >> _onewgts >> _twowgts >> _threewgts >> _zeromax >> _onemax >> _twomax 
     >> _threemax;
}
  
ClassDescription<a1ThreePionCLEODecayer> a1ThreePionCLEODecayer::inita1ThreePionCLEODecayer;
// Definition of the static class description member.
  
void a1ThreePionCLEODecayer::Init() {
  
  static ClassDocumentation<a1ThreePionCLEODecayer> documentation
    ("The \\classname{a1ThreePionCLEODecayer} class performs the decay of the "
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
     0, 0, 0, 0, 2.*pi, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhomagD,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhophaseD,
     0, 0, 0, 0, 2.*pi, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &a1ThreePionCLEODecayer::_f0phase, 0.54*pi, 0.0, 2.*pi,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &a1ThreePionCLEODecayer::_f2phase, 0.56*pi, 0.0, 2.*pi,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &a1ThreePionCLEODecayer::_sigmaphase, 0.23*pi, 0.0, 2.*pi,
     false, false, true);

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

  static ParVector<a1ThreePionCLEODecayer,bool> interfacezerochan
    ("AllNeutralChannels",
     "The channels to use for the integration of the decay a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zerochan,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacezerowgts
    ("AllNeutralWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zerowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,bool> interfaceonechan
    ("OneChargedChannels",
     "The channels to use for the integration of the decay a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onechan,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,bool> interfacetwochan
    ("TwoChargedChannels",
     "The channels to use for the integration of the decay a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twochan,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,bool> interfacethreechan
    ("ThreeChargedChannels",
     "The channels to use for the integration of the decay a_1^+->pi+pi+pi-",
     &a1ThreePionCLEODecayer::_threechan,
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
a1ThreePionCLEODecayer::decayCurrent(const bool vertex, const int imode,
				     const int ichan,const Particle & inpart,
				     const ParticleVector &outpart) const
{
  //overall factor
  double fact=_coupling;
  // momentum of the incoming particle
  Lorentz5Momentum Q=inpart.momentum();
  Energy2 q2=Q.mass2();
  // construct the spin info objects if needed
  if(vertex)
    {
      for(unsigned int ix=0;ix<outpart.size();++ix)
	{
	  SpinPtr stemp= new_ptr(ScalarSpinInfo(outpart[ix]->momentum(),true));
	  outpart[ix]->spinInfo(stemp);
	}
    }
  // identify the mesons
  int npi0=0,npiplus=0,npiminus=0,idtemp;
  unsigned int iloc[3],ipi0[3],ipim[3],ipip[3];
  for(unsigned int ix=0; ix<outpart.size();++ix)
    {
      idtemp=outpart[ix]->id();
      if(idtemp==ParticleID::piplus){ipip[npiplus]=ix;++npiplus;}
      else if(idtemp==ParticleID::piminus){ipim[npiminus]=ix;++npiminus;}
      else if(idtemp==ParticleID::pi0){ipi0[npi0]=ix;++npi0;}
    }
  // work out which decay mode we are doing
  idtemp=inpart.id();
  Complex F1=0.,F2=0.,F3=0.;
  // a_1^0->pi0pi0pi0
  if(imode==0)
    {
      fact*=0.40824829;
      Lorentz5Momentum ps1=outpart[1]->momentum()+outpart[2]->momentum();
      Lorentz5Momentum ps2=outpart[0]->momentum()+outpart[2]->momentum();
      Lorentz5Momentum ps3=outpart[0]->momentum()+outpart[1]->momentum();
      ps1.rescaleMass();ps2.rescaleMass();ps3.rescaleMass();
      iloc[0]=ipi0[0];iloc[1]=ipi0[1];iloc[2]=ipi0[2];
      Energy2 s1=ps1.mass2(),s2=ps2.mass2(),s3=ps3.mass2();
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
      else if(ichan==21){F2=-2./3.*_sigmacoup*sigbws1;F3=-2./3.*_sigmacoup*sigbws1;}
      else if(ichan==22){F1=-2./3.*_sigmacoup*sigbws2;F3=+2./3.*_sigmacoup*sigbws2;}
      else if(ichan==23){F1= 2./3.*_sigmacoup*sigbws3;F2= 2./3.*_sigmacoup*sigbws3;}
      else if(ichan==24)
	{
	  Complex Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
	  F1+=_f2coup*0.5*(s3-s2)*f2bws1;F2-=_f2coup*Dfact1; F3-=_f2coup*Dfact1;
	}
      else if(ichan==25)
	{
	  Complex Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
	  F2+=_f2coup*0.5*(s3-s1)*f2bws2;F1-=_f2coup*Dfact2;F3+=_f2coup*Dfact2;
	}
      else if(ichan==26)
	{
	  Complex Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
	  F3+=-_f2coup*0.5*(s1-s2)*f2bws3;F1+=_f2coup*Dfact3;F2+=_f2coup*Dfact3;
	}
      else if(ichan==27){F2=-2./3.*_f0coup*f0bws1;F3=-2./3.*_f0coup*f0bws1;}
      else if(ichan==28){F1=-2./3.*_f0coup*f0bws2;F3=+2./3.*_f0coup*f0bws2;}
      else if(ichan==29){F1= 2./3.*_f0coup*f0bws3;F2= 2./3.*_f0coup*f0bws3;}
    }
  // a_1^+ ->pi0pi0pi+
  else if(imode==1)
    {
      fact*=0.70710678;
      if(idtemp==ParticleID::a_1minus){ipip[0]=ipim[0];}
      iloc[0]=ipi0[0];iloc[1]=ipi0[1];iloc[2]=ipip[0];
      // compute the breit-wigners we need
      Lorentz5Momentum ps1=outpart[ipi0[1]]->momentum()+outpart[ipip[0]]->momentum();
      Lorentz5Momentum ps2=outpart[ipi0[0]]->momentum()+outpart[ipip[0]]->momentum();
      Lorentz5Momentum ps3=outpart[ipi0[0]]->momentum()+outpart[ipi0[1]]->momentum();
      ps1.rescaleMass();ps2.rescaleMass();ps3.rescaleMass();
      Energy2 s1=ps1.mass2(),s2=ps2.mass2(),s3=ps3.mass2();
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
  // a_1^0 ->pi+pi-pi0
  else if(imode==2)
    {
      iloc[0]=ipip[0];iloc[1]=ipim[0];iloc[2]=ipi0[0];
      Lorentz5Momentum ps1=outpart[ipim[0]]->momentum()+outpart[ipi0[0]]->momentum();
      Lorentz5Momentum ps2=outpart[ipip[0]]->momentum()+outpart[ipi0[0]]->momentum();
      Lorentz5Momentum ps3=outpart[ipip[0]]->momentum()+outpart[ipim[0]]->momentum();
      ps1.rescaleMass();ps2.rescaleMass();ps3.rescaleMass();
      Energy2 s1=ps1.mass2(),s2=ps2.mass2(),s3=ps3.mass2();
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
      else if(ichan%2==0&&ichan>=30&&ichan<=34)
	{
	  unsigned int ires=(ichan-30)/2;
	  if(ires<_rhocoupP.size()){F1+=_rhocoupP[ires]*rhos1bw[ires];}
	  Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
	  if(ires<_rhocoupD.size())
	    {F2+=Dfact2*_rhocoupD[ires]*rhos1bw[ires];
	    F3+=_rhocoupD[ires]*Dfact2*rhos1bw[ires];}
	}
      else if(ichan%2==1&&ichan>=31&&ichan<=35)
	{
	  unsigned int ires=(ichan-31/2);
	  if(ires<_rhocoupP.size()){F2+=_rhocoupP[ires]*rhos2bw[ires];}
	  Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
	  if(ires<_rhocoupD.size())
	    {F1+=Dfact1*_rhocoupD[ires]*rhos2bw[ires];
	    F3-=_rhocoupD[ires]*-Dfact1*rhos2bw[ires];}
	}
      else if(ichan==36){F1+=2./3.*_sigmacoup*sigbw;F2+=2./3.*_sigmacoup*sigbw;}
      else if(ichan==37)
	{
	  Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
	  F1+=Dfact3;F2+=Dfact3;
	  F3-=0.5*_f2coup*(s1-s2)*f2bw;
	}
      else if(ichan==38){F1+=2./3.*_f0coup*f0bw;F2+=2./3.*_f0coup*f0bw;}
    }
  // a_1^+ ->pi+pi+pi-
  else if(imode==3)
    {
      fact*=0.70710678;
      // change the order  if needed
      if(idtemp==ParticleID::a_1minus)
	{npi0=ipip[0];ipip[0]=ipim[0];ipip[1]=ipim[1];ipim[0]=npi0;}
      iloc[0]=ipip[0];iloc[1]=ipip[1];iloc[2]=ipim[0];
      // compute the breit-wigners we need
      Lorentz5Momentum ps1=outpart[ipip[1]]->momentum()+outpart[ipim[0]]->momentum();
      Lorentz5Momentum ps2=outpart[ipip[0]]->momentum()+outpart[ipim[0]]->momentum();
      Lorentz5Momentum ps3=outpart[ipip[0]]->momentum()+outpart[ipip[1]]->momentum();
      ps1.rescaleMass();ps2.rescaleMass();ps3.rescaleMass();
      Energy2 s1=ps1.mass2(),s2=ps2.mass2(),s3=ps3.mass2();
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
  // use the form-factors to compute the current
  LorentzPolarizationVector output=
     F1*outpart[iloc[1]]->momentum()-F1*outpart[iloc[2]]->momentum()
    -F2*outpart[iloc[2]]->momentum()+F2*outpart[iloc[0]]->momentum()
    +F3*outpart[iloc[0]]->momentum()-F3*outpart[iloc[1]]->momentum();
  output*=fact;
  vector<LorentzPolarizationVector> temp;temp.push_back(output);
  return temp;
}
}
  
