// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau2MesonDecayer class.
//
//  Author: Peter Richardson
//

#include "Tau2MesonDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Handlers/StepHandler.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;
  
Tau2MesonDecayer::~Tau2MesonDecayer() {}
  
bool Tau2MesonDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  // can we handle this mode?
  // check the parent is tau+-/ and three decay products
  if(abs(dm.parent()->id())==15 && dm.products().size() == 3)
    {
      ParticleMSet::const_iterator pit  = dm.products().begin();
      ParticleMSet::const_iterator pend = dm.products().end();
      int id[2]; int idnu=0; int ic=0,idtemp;
      for( ; pit!=pend;++pit)
	{
	  idtemp=(**pit).id();
	  if(abs(idtemp)==16)
	    {idnu=idtemp;}
	  else
	    {id[ic]=idtemp;++ic;}
	}
      if(ic==2 && idnu*dm.parent()->id()>0)
	{
	  // pion modes
	  if((idnu ==ParticleID::nu_tau    && id[0]==ParticleID::piminus &&
	      id[1]==ParticleID::pi0) ||
	     (idnu ==ParticleID::nu_tau    && id[0]==ParticleID::pi0 &&
	      id[1]==ParticleID::piminus) ||
	     (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::piplus &&
	      id[1]==ParticleID::pi0) ||
	     (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::pi0 &&
	      id[1]==ParticleID::piplus))
	    {allowed=true;}
	  // single charged kaon
	  else if((idnu ==ParticleID::nu_tau    && id[0]==ParticleID::Kminus &&
		   id[1]==ParticleID::pi0) ||
		  (idnu ==ParticleID::nu_tau    && id[0]==ParticleID::pi0 &&
		   id[1]==ParticleID::Kminus) ||
		  (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::Kplus &&
		   id[1]==ParticleID::pi0) ||
		  (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::pi0 &&
		   id[1]==ParticleID::Kplus))
	    {allowed=true;}
	  // single neutral kaon
	  else if((idnu ==ParticleID::nu_tau    && id[0]==ParticleID::piminus &&
		   id[1]==ParticleID::Kbar0) ||
		  (idnu ==ParticleID::nu_tau    && id[0]==ParticleID::Kbar0 &&
		   id[1]==ParticleID::piminus) ||
		  (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::piplus &&
		   id[1]==ParticleID::K0) ||
		  (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::K0 &&
		   id[1]==ParticleID::piplus))
	    {allowed=true;}
	  // two kaons
	  else if((idnu ==ParticleID::nu_tau    && id[0]==ParticleID::Kminus &&
		   id[1]==ParticleID::K0) ||
		  (idnu ==ParticleID::nu_tau    && id[0]==ParticleID::K0 &&
		   id[1]==ParticleID::Kminus) ||
		  (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::Kplus &&
		   id[1]==ParticleID::Kbar0) ||
		  (idnu ==ParticleID::nu_taubar && id[0]==ParticleID::Kbar0 &&
		   id[1]==ParticleID::Kplus))
	    {allowed=true;}
	}
    }
  return allowed;
}
  
ParticleVector Tau2MesonDecayer::decay(const DecayMode & dm,
				       const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  int ikaon=1,id;
  unsigned int nkaon=0;
  for(unsigned int ix=0;ix<children.size();++ix)
    {
      id=children[ix]->id();
      if(id==ParticleID::Kbar0 || id==ParticleID::K0){ikaon=2;++nkaon;}
      else if (id==ParticleID::Kminus || id==ParticleID::Kplus){ikaon=3;++nkaon;}
    }
  if(nkaon==2){ikaon=4;}
  // perform the decay
  generate(true,ikaon-1,parent,children);
  return children;
}
  
void Tau2MesonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _pimodel << _Kmodel << _piwgt << _kwgt << _pimax << _K0max << _Kplusmax
     << _KKmax << _pichannels << _K0channels << _Kpluschannels << _KKchannels
     << _piwgts << _K0wgts << _Kpluswgts << _KKwgts << _rhoparameters
     << _Kstarparameters << _rhomasses << _rhowidths << _Kstarmasses << _Kstarwidths
     <<  _mass << _width << _mass2 << _massw << _massa <<_massb << _mom << _dhdq2
     << _hm2 << _dparam;
}
  
void Tau2MesonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _pimodel >> _Kmodel >> _piwgt >> _kwgt >> _pimax >> _K0max >> _Kplusmax
     >> _KKmax >> _pichannels >> _K0channels >> _Kpluschannels >> _KKchannels
     >> _piwgts >> _K0wgts >> _Kpluswgts >> _KKwgts >> _rhoparameters
     >> _Kstarparameters >> _rhomasses >> _rhowidths >> _Kstarmasses >> _Kstarwidths
     >> _mass >> _width >> _mass2 >> _massw >> _massa >> _massb >> _mom >> _dhdq2
     >> _hm2 >> _dparam;
}
  
ClassDescription<Tau2MesonDecayer> Tau2MesonDecayer::initTau2MesonDecayer;
  // Definition of the static class description member.
  
void Tau2MesonDecayer::Init() {
  
  static ParVector<Tau2MesonDecayer,double> interfacePiWeights
    ("PiWeights",
     "The weights of the different channels for the multichannel integration of"
     " tau -> nu pi pi",
     &Tau2MesonDecayer::_piwgts,
     0, 0, 0, 0., 1., false, false, true);
  
  static ParVector<Tau2MesonDecayer,double> interfaceK0Weights
    ("K0Weights",
     "The weights of the different channels for the multichannel integration of"
     " tau -> nu K0 pi",
     &Tau2MesonDecayer::_K0wgts,
     0, 0, 0, 0., 1., false, false, true);
  
  static ParVector<Tau2MesonDecayer,double> interfaceKplusWeights
    ("KplusWeights",
     "The weights of the different channels for the multichannel integration of"
     " tau -> nu K- pi",
     &Tau2MesonDecayer::_Kpluswgts,
     0, 0, 0, 0., 1., false, false, true);
  
  static ParVector<Tau2MesonDecayer,double> interfaceKKWeights
    ("KKWeights",
     "The weights of the different channels for the multichannel integration of"
     " tau -> nu K K",
     &Tau2MesonDecayer::_KKwgts,
     0, 0, 0, 0., 1., false, false, true);
  
  static Parameter<Tau2MesonDecayer,double> interfacePiMaximum
    ("PiMaximum",
     "The maximum weight for the integration of the decay tau -> nu pi pi",
     &Tau2MesonDecayer::_pimax, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau2MesonDecayer,double> interfaceK0Maximum
    ("K0Maximum",
     "The maximum weight for the integration of the decay tau -> nu K0 pi",
     &Tau2MesonDecayer::_K0max, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau2MesonDecayer,double> interfaceKplusMaximum
    ("KplusMaximum",
     "The maximum weight for the integration of the decay tau -> nu K- pi",
     &Tau2MesonDecayer::_Kplusmax, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau2MesonDecayer,double> interfaceKKMaximum
    ("KKMaximum",
     "The maximum weight for the integration of the decay tau -> nu K K",
     &Tau2MesonDecayer::_KKmax, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  
  static ParVector<Tau2MesonDecayer,double> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the decay tau -> nu pi pi",
     &Tau2MesonDecayer::_rhomasses,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<Tau2MesonDecayer,double> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the decay tau -> nu pi pi",
     &Tau2MesonDecayer::_rhowidths,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<Tau2MesonDecayer,double> interfaceKstarMasses
    ("KstarMasses",
     "The masses of the different rho resonances for the decay tau -> nu pi pi",
     &Tau2MesonDecayer::_Kstarmasses,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<Tau2MesonDecayer,double> interfaceKstarWidths
    ("KstarWidths",
     "The widths of the different rho resonances for the decay tau -> nu pi pi",
     &Tau2MesonDecayer::_Kstarwidths,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static Switch<Tau2MesonDecayer,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values for the rho meson masses and widths",
     &Tau2MesonDecayer::_rhoparameters, true, false, false);
  static SwitchOption interfaceRhoParameterstrue
    (interfaceRhoParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceRhoParametersParticleData
    (interfaceRhoParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static Switch<Tau2MesonDecayer,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values for the Kstar meson masses and widths",
     &Tau2MesonDecayer::_Kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static ParVector<Tau2MesonDecayer,double> interfacepiwgt
    ("PiWeight",
     "The weights of the different resonances for the decay tau -> nu pi pi",
     &Tau2MesonDecayer::_piwgt,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<Tau2MesonDecayer,double> interfacekwgt
    ("KWeight",
     "The weights of the different resonances for the decay tau -> nu K pi",
     &Tau2MesonDecayer::_kwgt,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Switch<Tau2MesonDecayer,int> interfacePiModel
    ("PiModel",
     "The model to use for the propagator for the pion modes.",
     &Tau2MesonDecayer::_pimodel, 0, false, false);
  static SwitchOption interfacePiModelKuhn
    (interfacePiModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfacePiModelGounaris
    (interfacePiModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);
  
  static Switch<Tau2MesonDecayer,int> interfaceKModel
    ("KModel",
     "The model to use for the propagator for the kaon modes.",
     &Tau2MesonDecayer::_pimodel, 0, false, false);
  static SwitchOption interfaceKModelKuhn
    (interfaceKModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfaceKModelGounaris
    (interfaceKModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);

  static ClassDocumentation<Tau2MesonDecayer> documentation
    ("The \\classname{Tau2MesonDecayer} class is designed to implement the decay "
     "of the tau to two scalar mesons using the models of either Kuhn and "
     "Santamaria (Z. Phys. C48, 445 (1990)) or Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968).  The mixing parameters are taken from "
     "Phys. Rev. D61:112002,2000 (CLEO), although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.");
  
}

// the hadronic current for this decay mode
vector<LorentzPolarizationVector> 
Tau2MesonDecayer::hadronCurrent(bool vertex, const int imode, const int ichan,
				const Particle & inpart,
				const ParticleVector & outpart) const
{
  // storage for the current
  vector<LorentzPolarizationVector> temp;
  // find the mesons
  int nmes=0,imes[3];
  for(unsigned int ix=0;ix<outpart.size();++ix)
    {if(abs(outpart[ix]->id())!=ParticleID::nu_tau){imes[nmes]=ix;++nmes;}}
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff=outpart[imes[0]]->momentum()-outpart[imes[1]]->momentum();
  Lorentz5Momentum psum=outpart[imes[0]]->momentum()+outpart[imes[1]]->momentum();
  // mass2 of vector intermediate state
  Energy2 q2=psum.m2();
  // calculate the current
  complex<double> FPI=0,denom=0.,vect[4];
  // pion mode
  if(imode==0)
    {
      if(ichan<0)
	{
	  for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix)
	    {
	      FPI+=_piwgt[ix]*BreitWigner(q2,_pimodel,0,ix);
	      denom+=_piwgt[ix];
	    }
	}
      else
	{
	  unsigned int ix=ichan/2; denom=1.;
	  if(ix<_piwgt.size()&&ix<3)
	    {FPI=_piwgt[ix]*BreitWigner(q2,_pimodel,0,ix);}
	}
      FPI = FPI/denom*sqrt(2.0*SM().CKM(0,0));
      // compute the current
      for(unsigned int ix=0;ix<4;++ix){vect[ix]=FPI*pdiff[ix];}
    }
  // single kaon modes
  else if(imode==1||imode==2)
    {
      if(ichan<0)
	{
	  for(unsigned int ix=0;ix<_kwgt.size()&&ix<3;++ix)
	    {
	      FPI+=_kwgt[ix]*BreitWigner(q2,_Kmodel,1,ix);
	      denom+=_kwgt[ix];
	    }
	}
      else
	{
	  unsigned int ix=0;
	  if(imode==1){ix=(ichan-7)/2;}
	  else if(imode==2){ix=(ichan-6)/2;}
	  denom=1.;
	  if(ix<_kwgt.size()&&ix<3)
	    {FPI=_kwgt[ix]*BreitWigner(q2,_Kmodel,1,ix);}
	}
      if(imode==1){FPI = FPI/denom*sqrt(4.*SM().CKM(0,1)/3.);}
      else{FPI = FPI/denom*sqrt(2.*SM().CKM(0,1)/3.);}
      double dot = psum.dot(pdiff)/q2;
      // compute the current
      for(unsigned int ix=0;ix<4;++ix){vect[ix]=FPI*(pdiff[ix]-psum[ix]*dot);}
    }
  // two kaon modes
  else
    {
      if(ichan<0)
	{
	  for(unsigned int ix=0;ix<_piwgt.size()&&ix<3;++ix)
	    {
	      FPI+=_piwgt[ix]*BreitWigner(q2,_pimodel,0,ix);
	      denom+=_piwgt[ix];
	    }
	}
      else
	{
	  unsigned int ix=(ichan-1)/2; denom=1.;
	  if(ix<_piwgt.size()&&ix<3)
	    {FPI=_piwgt[ix]*BreitWigner(q2,_pimodel,0,ix);}
	}
      FPI = FPI/denom*sqrt(2.0*SM().CKM(0,0));
      double dot = psum.dot(pdiff)/q2;
      // compute the current
      for(unsigned int ix=0;ix<4;++ix){vect[ix]=FPI*(pdiff[ix]-psum[ix]*dot);}
    }
  // insert the current
  temp.push_back(LorentzPolarizationVector(vect[0],vect[1],vect[2],vect[3]));
  // set up the spininfo for the decay products
  if(vertex)
    {
      for(unsigned int ix=0;ix<2;++ix)
	{
	  SpinPtr stemp = new_ptr(ScalarSpinInfo(outpart[imes[ix]]->momentum(),true));
	  outpart[imes[ix]]->spinInfo(stemp);
	}
    }
  // return the answer
  return temp;
}

}
