// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1ThreePionDecayer class.
//

#include "a1ThreePionDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "a1ThreePionDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"

namespace Herwig{

using namespace ThePEG;
using namespace Helicity;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::ScalarSpinInfo;
  
a1ThreePionDecayer::~a1ThreePionDecayer() {}
  
bool a1ThreePionDecayer::accept(const DecayMode & dm) const {
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
  
ParticleVector a1ThreePionDecayer::decay(const DecayMode & dm,
					 const Particle & parent) const {
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit){if(abs((**pit).id())==ParticleID::piplus){++ncharged;}}
  bool cc = parent.id()==ParticleID::a_1minus;
  return generate(true,cc,ncharged,parent);
}
  
  
void a1ThreePionDecayer::persistentOutput(PersistentOStream & os) const {
  os << _rhomass << _rhowidth << _prho << _hm2 << _rhoD << _dhdq2m2 <<  _sigmamass
     << _sigmawidth << _psigma << _mpi << _mpi2 << _lambda2 << _a1mass2 << _zsigma 
     << _rhocoupling << _coupling << _localparameters << _zerowgts << _onewgts 
     << _twowgts << _threewgts << _zeromax
     << _onemax << _twomax << _threemax << _coupling;
}
  
void a1ThreePionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _rhomass >> _rhowidth >> _prho >> _hm2 >> _rhoD >> _dhdq2m2 >>  _sigmamass
     >> _sigmawidth >> _psigma >> _mpi >> _mpi2 >> _lambda2 >> _a1mass2 >> _zsigma 
     >> _rhocoupling >> _coupling >> _localparameters >> _zerowgts >> _onewgts 
     >> _twowgts >> _threewgts >> _zeromax
     >> _onemax >> _twomax >> _threemax >> _coupling;
}

ClassDescription<a1ThreePionDecayer> a1ThreePionDecayer::inita1ThreePionDecayer;
// Definition of the static class description member.
  
void a1ThreePionDecayer::Init() {
    
  static ClassDocumentation<a1ThreePionDecayer> documentation
    ("The a1ThreePionDecayer class is designed to decay the a_1 "
     "resonance to three pions using a model based on that used in the modelling "
     "of tau->4 pions.");


  static Switch<a1ThreePionDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1ThreePionDecayer::_localparameters, true, false, false);
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

  static Parameter<a1ThreePionDecayer,double> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1ThreePionDecayer::_coupling, 80.76, 0.0, 1000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of the mass scale squared to use in the form-factor",
     &a1ThreePionDecayer::_lambda2, GeV2, 1.2*GeV2, 0.0001*GeV2, 10.0*GeV2,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy2> interfacea1mass2
    ("a1mass2",
     "The local value of the square of the a_1 mass",
     &a1ThreePionDecayer::_a1mass2, GeV2, 1.5129*GeV2, 0.5*GeV2, 10.0*GeV2,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy> interfaceSigmaMass
    ("SigmaMass",
     "The local value of the sigma mass",
     &a1ThreePionDecayer::_sigmamass, GeV, 0.8*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy> interfaceSigmaWidth
    ("SigmaWidth",
     "The local value of the sigma width",
     &a1ThreePionDecayer::_sigmawidth, GeV, 0.8*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacezerowgts
    ("AllNeutralWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi0pi0pi0",
     &a1ThreePionDecayer::_zerowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1ThreePionDecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1ThreePionDecayer::_twowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1ThreePionDecayer::_threewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceZeroMax
    ("ZeroMax",
     "The maximum weight for the integration fo the channel a_1^0->pi0pi0pi0",
     &a1ThreePionDecayer::_zeromax, 107.793, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1ThreePionDecayer::_onemax, 1088.96, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1ThreePionDecayer::_twomax, 1750.73, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1ThreePionDecayer::_threemax, 739.334, 0.0, 10000.0,
     false, false, true);

  static ParVector<a1ThreePionDecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1ThreePionDecayer::_rhomass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1ThreePionDecayer::_rhowidth,
     0, 0, 0, -10000, 10000, false, false, true);

}

// hadronic current
vector<LorentzPolarizationVector> 
a1ThreePionDecayer::decayCurrent(const bool vertex, const int ichan, 
				 const Particle & inpart,
				 const ParticleVector &outpart) const
{
  // momentum of the incoming particle
  Lorentz5Momentum Q=inpart.momentum();
  // vector for the output
  LorentzPolarizationVector output=LorentzPolarizationVector();
  // construct the spin info objects if needed
  if(vertex)
    {for(unsigned int ix=0;ix<outpart.size();++ix)
	{outpart[ix]->spinInfo(new_ptr(ScalarSpinInfo(outpart[ix]->momentum(),true)));}}
  // identify the mesons
  int npi0=0,npiplus=0,npiminus=0,ipi0[3],ipim[3],ipip[3],idtemp;
  for(unsigned int ix=0; ix<outpart.size();++ix)
    {
      idtemp=outpart[ix]->id();
      if(idtemp==ParticleID::piplus){ipip[npiplus]=ix;++npiplus;}
      else if(idtemp==ParticleID::piminus){ipim[npiminus]=ix;++npiminus;}
      else if(idtemp==ParticleID::pi0){ipi0[npi0]=ix;++npi0;}
    }
  // work out which decay mode we are doing
  idtemp=inpart.id();
  // a_10 -> pi0pi0pi0
  if(imode()==0)
    {
      // the momenta
      Lorentz5Momentum pa=outpart[1]->momentum()+outpart[2]->momentum();
      pa.rescaleMass();
      Lorentz5Momentum pb=outpart[0]->momentum()+outpart[2]->momentum();
      pb.rescaleMass();
      Lorentz5Momentum pc=outpart[0]->momentum()+outpart[1]->momentum();
      pc.rescaleMass();
      //the breit-wigners
      Complex sig1=sigmaBreitWigner(pa.mass2());
      Complex sig2=sigmaBreitWigner(pb.mass2());
      Complex sig3=sigmaBreitWigner(pc.mass2());
      // compute the vector
      if(ichan<0)
	{
	  output= 
	    sig1*(outpart[0]->momentum())
	    +sig2*(outpart[1]->momentum())
	    +sig3*(outpart[2]->momentum());
	}
      else if(ichan==15)
	{output=sig1*(outpart[0]->momentum());}
      else if(ichan==16)
	{output=sig2*(outpart[1]->momentum());}
      else if(ichan==17)
	{output=sig3*(outpart[2]->momentum());}
      else
	{cerr << "Unknown channel for a_1^0->pi0pi0pi0 in " 
	      << "a1ThreePionDecayer::decayCurrent" << endl;
	output=LorentzPolarizationVector();}
      // the coupling z and identical particle factor
      output *=_zsigma*0.4082482904638631*Q.mass2();
    }
  // a_1+ -> pi0pi0pi+
  else if(imode()==1)
    {
      if(idtemp==ParticleID::a_1minus){ipip[0]=ipim[0];}
      Lorentz5Momentum pa=outpart[ipip[0]]->momentum()+outpart[ipi0[1]]->momentum();
      pa.rescaleMass();
      Lorentz5Momentum pb=outpart[ipip[0]]->momentum()+outpart[ipi0[0]]->momentum();
      pb.rescaleMass();
      Lorentz5Momentum pc=outpart[ipi0[0]]->momentum()+outpart[ipi0[1]]->momentum();
      pc.rescaleMass();
      // scalar propagator
      Complex sig1 = sigmaBreitWigner(pc.mass2());
      // sigma terms
      if(ichan<0||ichan==14)
	{output = _zsigma*Q.mass2()*sig1*(outpart[ipip[0]]->momentum());}
      // the rho terms
      complex<Energy2> dot01=Q*outpart[ipi0[0]]->momentum();
      complex<Energy2> dot02=Q*outpart[ipi0[1]]->momentum();
      complex<Energy2> dot03=Q*outpart[ipip[0]]->momentum();
      for(int ix=0,N=_rhocoupling.size();ix<N;++ix)
	{
	  Complex rho1=_rhocoupling[ix]*rhoBreitWigner(pa.mass2(),ix);
	  Complex rho2=_rhocoupling[ix]*rhoBreitWigner(pb.mass2(),ix);
	  if(ichan<0||ichan==8+2*ix)
	    {output +=rho1*(dot03*(outpart[ipi0[1]]->momentum())
			    -dot02*(outpart[ipip[0]]->momentum()));}
	  if(ichan<0||ichan==9+2*ix)
	    {output +=rho2*(dot03*(outpart[ipi0[0]]->momentum())
			    -dot01*(outpart[ipip[0]]->momentum()));}
	}
      // the identical particle factor
      output *=0.7071067811865476;
    }
  // a_10->pi+pi-pi0
  else if(imode()==2)
    {
      // the momenta
      Lorentz5Momentum pa=outpart[ipi0[0]]->momentum()+outpart[ipim[0]]->momentum();
      pa.rescaleMass();
      Lorentz5Momentum pb=outpart[ipi0[0]]->momentum()+outpart[ipip[0]]->momentum();
      pb.rescaleMass();
      Lorentz5Momentum pc=outpart[ipip[0]]->momentum()+outpart[ipim[0]]->momentum();
      pc.rescaleMass();
      // the sigma terms
      Complex sig1=sigmaBreitWigner(pc.mass2());
      if(ichan<0||ichan==24)
	{output = _zsigma*Q.mass2()*sig1*outpart[ipi0[0]]->momentum();}
      // rho terms
      complex<Energy2> dot01=Q*outpart[ipip[0]]->momentum();
      complex<Energy2> dot02=Q*outpart[ipim[0]]->momentum();
      complex<Energy2> dot03=Q*outpart[ipi0[0]]->momentum();
      for(int ix=0,N=_rhocoupling.size();ix<N;++ix)
	{
	  Complex rho1=_rhocoupling[ix]*rhoBreitWigner(pa.mass2(),ix);
	  Complex rho2=_rhocoupling[ix]*rhoBreitWigner(pb.mass2(),ix);
	  if(ichan<0||ichan==18+2*ix)
	    {output+=rho1*(dot03*(outpart[ipim[0]]->momentum())
			   -dot02*(outpart[ipi0[0]]->momentum()));}
	  if(ichan<0||ichan==19+2*ix)
	    {output+=rho2*(dot03*(outpart[ipip[0]]->momentum())
			   -dot01*(outpart[ipi0[0]]->momentum()));}
	}
    }
  // a1+ -> pi+pi+pi-
  else if(imode()==3)
    {
      // change the order  if needed
      if(idtemp==ParticleID::a_1minus)
	{npi0=ipip[0];ipip[0]=ipim[0];ipip[1]=ipim[1];ipim[0]=npi0;}
      // momenta of the intermediates
      Lorentz5Momentum pa=outpart[ipip[1]]->momentum()+outpart[ipim[0]]->momentum();
      pa.rescaleMass();
      Lorentz5Momentum pb=outpart[ipip[0]]->momentum()+outpart[ipim[0]]->momentum();
      pb.rescaleMass();
      Lorentz5Momentum pc=outpart[ipip[0]]->momentum()+outpart[ipip[1]]->momentum();
      pc.rescaleMass();
      // the scalar propagators 
      Complex sig1=sigmaBreitWigner(pa.mass2());
      Complex sig2=sigmaBreitWigner(pb.mass2());
      // sigma terms
      if(ichan<0||ichan==6)
	{output+=sig1*(outpart[ipip[0]]->momentum());}
      if(ichan<0||ichan==7)
	{output+=sig2*(outpart[ipip[1]]->momentum());}
      output *=_zsigma*Q.mass2();
      // rho terms
      complex<Energy2> dot01=Q*outpart[ipip[0]]->momentum();
      complex<Energy2> dot02=Q*outpart[ipip[1]]->momentum();
      complex<Energy2> dot03=Q*outpart[ipim[0]]->momentum();
      for(int ix=0,N=_rhocoupling.size();ix<N;++ix)
	{
	  Complex rho1 = _rhocoupling[ix]*rhoBreitWigner(pa.mass2(),ix);
	  Complex rho2 = _rhocoupling[ix]*rhoBreitWigner(pb.mass2(),ix);
	  if(ichan<0||ichan==2*ix)
	    {output-=rho1*( dot03*(outpart[ipip[1]]->momentum())
			    -dot02*(outpart[ipim[0]]->momentum()));}
	  if(ichan<0||ichan==2*ix+1)
	    {output-=rho2*( dot03*(outpart[ipip[0]]->momentum())
			    -dot01*(outpart[ipim[0]]->momentum()));}
	}
      // the identical particle factor
      output *=0.7071067811865476;
    }
  // form-factor
  output*=a1FormFactor(Q.mass2())*_coupling/(Q.mass()*_rhomass[0]*_rhomass[0]);
  // return the answer
  return vector<LorentzPolarizationVector>(1,output);
}

}
  
