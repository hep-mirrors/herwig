// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau2MesonPhoton class.
//
//  Author: Peter Richardson
//

#include "Tau2MesonPhotonDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::VectorSpinPtr;
using ThePEG::Helicity::tcVectorSpinPtr;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;
  
Tau2MesonPhotonDecayer::~Tau2MesonPhotonDecayer() {}
  
bool Tau2MesonPhotonDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  // can we handle this mode?
  int idtau=dm.parent()->id();
  // check the parent is tau+-/ and four decay products
  if(abs(dm.parent()->id())==15 && dm.products().size() == 4)
    {
      int npiminus=0,npiplus=0,npi0=0,ngamma=0;
      ParticleMSet::const_iterator pit  = dm.products().begin();
      ParticleMSet::const_iterator pend = dm.products().end();
      int idnu=0,idtemp=0;
      for( ; pit!=pend;++pit)
	{
	  idtemp=(**pit).id();
	  if(abs(idtemp)==16){idnu=idtemp;}
	  else if(idtemp==ParticleID:: piplus){++npiplus;}
	  else if(idtemp==ParticleID::piminus){++npiminus;}
	  else if(idtemp==ParticleID::gamma){++ngamma;}
	  else if(idtemp==ParticleID::pi0){++npi0;}
	}
      // is it O.K
      if((idtau==ParticleID::tauminus&&idnu==ParticleID::nu_tau&&
	  npiminus==1 && npi0==1 && ngamma==1)||
	 (idtau==ParticleID::tauplus&&idnu==ParticleID::nu_taubar&&
	  npiplus==1&& npi0==1&&ngamma==1)){allowed=true;}
    }
  return allowed;
}
 
ParticleVector Tau2MesonPhotonDecayer::decay(const DecayMode & dm,
					     const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  // perform the decay
  generate(true,0,parent,children);
  return children;
}
  
  
void Tau2MesonPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _grho << _grhoomegapi << _maxwgt << _on << _weights << _resweights
     << _rhoparameters << _rhomasses << _rhowidths
     << _omegaparameters << _omegamass << _omegawidth ;}
  
void Tau2MesonPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _grho >> _grhoomegapi >> _maxwgt >> _on >> _weights >> _resweights
     >> _rhoparameters >> _rhomasses >> _rhowidths
     >> _omegaparameters >> _omegamass >> _omegawidth ;}
  
ClassDescription<Tau2MesonPhotonDecayer> Tau2MesonPhotonDecayer::initTau2MesonPhotonDecayer;
  // Definition of the static class description member.
  
void Tau2MesonPhotonDecayer::Init() {
  
  static ParVector<Tau2MesonPhotonDecayer,double> interfacereswgt
    ("Weights",
     "The weights of the different resonances for the decay tau -> nu pi pi gamma",
     &Tau2MesonPhotonDecayer::_resweights,
     0, 0, 0, -1000, 1000, false, false, true);
  
  
  static ParVector<Tau2MesonPhotonDecayer,double> interfaceintwgt
    ("IntegrationWeights",
     "The weights of the integration channels for the decay tau -> nu pi pi gamma",
     &Tau2MesonPhotonDecayer::_weights,
     0, 0, 0, 0, 1.0, false, false, true);
  
  
  static Parameter<Tau2MesonPhotonDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the integration",
     &Tau2MesonPhotonDecayer::_maxwgt, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  
  static Switch<Tau2MesonPhotonDecayer,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values for the rho meson masses and widths",
     &Tau2MesonPhotonDecayer::_rhoparameters, true, false, false);

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
  
  static Switch<Tau2MesonPhotonDecayer,bool> interfaceomegaParameters
    ("omegaParameters",
     "Use local values for the omega meson masses and widths",
     &Tau2MesonPhotonDecayer::_omegaparameters, true, false, false);
  static SwitchOption interfaceomegaParameterstrue
    (interfaceomegaParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceomegaParametersParticleData
    (interfaceomegaParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static ParVector<Tau2MesonPhotonDecayer,double> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the different rho resonances for the decay tau ->  pi pi photon",
     &Tau2MesonPhotonDecayer::_rhomasses,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static ParVector<Tau2MesonPhotonDecayer,double> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances for the decay tau -> nu pi pi photon",
     &Tau2MesonPhotonDecayer::_rhowidths,
     0, 0, 0, -10000, 10000, false, false, true);
  
  static Parameter<Tau2MesonPhotonDecayer,double> interfaceomegamass
    ("omegamass",
     "The mass of the omega for the decay tau- -> pi pi photon",
     &Tau2MesonPhotonDecayer::_omegamass, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau2MesonPhotonDecayer,double> interfaceomegawidth
    ("omegawidth",
     "The width of the omega for the decay tau- -> pi pi photon",
     &Tau2MesonPhotonDecayer::_omegawidth, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static ClassDocumentation<Tau2MesonPhotonDecayer> documentation
    ("The \\classname{Tau2MesonPhotonDecayer} class implements the decay "
     "tau+/- -> pi+/- pi0 gamma via an omega.");
  
  static Parameter<Tau2MesonPhotonDecayer,Energy2> interfacegrho
    ("grho",
     "The rho meson decay constant.",
     &Tau2MesonPhotonDecayer::_grho, GeV2, 0.11238947*GeV2, -1.0e12*GeV2, 1.0e12*GeV2,
     false, false, false);

  static Parameter<Tau2MesonPhotonDecayer,InvEnergy> interfacegrhoomegapi
    ("grhoomegapi",
     "The rho-omega-pi coupling",
     &Tau2MesonPhotonDecayer::_grhoomegapi, 1./GeV, 12.924/GeV,
     -1.0e12*1./GeV, 1.0e12*1./GeV,
     false, false, false);

}
  
// the hadronic current for this decay mode
vector<LorentzPolarizationVector>  
Tau2MesonPhotonDecayer::hadronCurrent(bool vertex,const int imode, const int ichan,
				      const Particle & inpart,
				      const ParticleVector & outpart) const
{
  vector<LorentzPolarizationVector> temp;
  complex<double> vect[4];
  // locate the particles
  unsigned int ipi0=0,ipic=0,igamma=0,id;
  Lorentz5Momentum pout;
  for(unsigned int ix=0,N=outpart.size();ix<N;++ix)
    {
      id=abs(outpart[ix]->id());
      if(id==ParticleID::pi0){ipi0=ix;}
      else if(id==ParticleID::piplus){ipic=ix;}
      else if(id==ParticleID::gamma){igamma=ix;}
      if(id!=ParticleID::nu_tau){pout+=outpart[ix]->momentum();}
    }
  // overall hadronic mass
  Energy2 q2 = pout.m2();
  // mass of the omega
  pout = outpart[igamma]->momentum()+outpart[ipi0]->momentum();
  Energy s2  = pout.m2();
  // compute the prefactor
  complex<double> prefactor=-FFunction(0.,-1)*FFunction(q2,ichan)*
    sqrt(0.5*SM().CKM(0,0))*sqrt(4.*pi*SM().alphaEM())*BreitWigner(s2,10);
  // dot products which don't depend on the polarization vector
  Energy2 dot12 = outpart[igamma]->momentum().dot(outpart[ipi0]->momentum());
  Energy2 dot13 = outpart[igamma]->momentum().dot(outpart[ipic]->momentum());
  Energy2 dot23 = outpart[ipi0]->momentum().dot(outpart[ipic]->momentum());
  Energy2 mpi2  = outpart[ipic]->mass()*outpart[ipic]->mass();
  // construct the spininfomation objects for the decay products
  VectorSpinPtr hwtemp;
  if(vertex)
    {
      SpinPtr spi0=new_ptr(ScalarSpinInfo(outpart[ipi0]->momentum(),true));
      outpart[ipi0]->spinInfo(spi0);
      SpinPtr spic=new_ptr(ScalarSpinInfo(outpart[ipic]->momentum(),true));
      outpart[ipic]->spinInfo(spic);
      SpinPtr sgamma=new_ptr(VectorSpinInfo(outpart[igamma]->momentum(),true));
      outpart[igamma]->spinInfo(sgamma);
      hwtemp= dynamic_ptr_cast<VectorSpinPtr>(sgamma);
    }
  // loop over the photon polarizations
  VectorWaveFunction photon(outpart[igamma]->momentum(),outpart[igamma]->dataPtr(),
			    outgoing);
  complex<double> dote2,dote3,coeffa,coeffb,coeffc;
  for(int ihel=-1;ihel<2;++ihel)
    {
      if(ihel!=0)
	{
	  photon.reset(ihel);
	  if(vertex){hwtemp->setBasisState(ihel,photon.Wave());}
	  // obtain the dot products we need
	  dote2 =
	    +photon.t()*outpart[ipi0]->momentum().e()
	    -photon.x()*outpart[ipi0]->momentum().px()
	    -photon.y()*outpart[ipi0]->momentum().py()
	    -photon.z()*outpart[ipi0]->momentum().pz();
	  dote3 =
	    +photon.t()*outpart[ipic]->momentum().e()
	    -photon.x()*outpart[ipic]->momentum().px()
	    -photon.y()*outpart[ipic]->momentum().py()
	    -photon.z()*outpart[ipic]->momentum().pz();
	  // now compute the coefficients
	  coeffa = mpi2*dot13-dot12*(dot23-dot13);
	  coeffb = dote2*dot13-dote3*dot12;
	  coeffc = dote2*dot23-dote3*(mpi2+dot12);
	  // finally compute the current
	  for(unsigned int ix=0;ix<4;++ix)
	    {vect[ix]=
		 photon(ix)*coeffa
	       -outpart[ipi0]->momentum()[ix]*coeffb
	       +outpart[igamma]->momentum()[ix]*coeffc;
	    vect[ix]*=prefactor;
	    }
	  temp.push_back(LorentzPolarizationVector(vect[0],vect[1],
						   vect[2],vect[3]));
	}
      else
	{
	  LorentzPolarizationVector veczero=LorentzPolarizationVector();
	  temp.push_back(veczero);
	  if(vertex){hwtemp->setBasisState(ihel,photon.Wave());}
	}
    }
  return temp;
}

}
