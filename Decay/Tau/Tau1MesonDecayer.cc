// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau1MesonDecayer class.
//
//  Author: Peter Richardson
//

#include "Tau1MesonDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity//ScalarSpinInfo.h"
#include "ThePEG/Helicity//VectorSpinInfo.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::VectorSpinPtr;
using ThePEG::Helicity::tcSpinInfoPtr;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

Tau1MesonDecayer::~Tau1MesonDecayer() {}
  
bool Tau1MesonDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  // can we handle this mode?
  // check the parent is tau+-/ and only two decay products
  if(abs(dm.parent()->id())==15 && dm.products().size() == 2)
    {
      ParticleMSet::const_iterator pit = dm.products().begin();
      int id1=(**pit).id();
      ++pit;
      int id2=(**pit).id();
      // pion modes
      if((id1==ParticleID::nu_tau    && id2==ParticleID::piminus) ||
	 (id1==ParticleID::nu_taubar && id2==ParticleID::piplus ) ||
	 (id2==ParticleID::nu_tau    && id1==ParticleID::piminus) ||
	 (id2==ParticleID::nu_taubar && id1==ParticleID::piplus ))
	{allowed=true;}
      else if((id1==ParticleID::nu_tau    && id2==ParticleID::Kminus) ||
	      (id1==ParticleID::nu_taubar && id2==ParticleID::Kplus ) ||
	      (id2==ParticleID::nu_tau    && id1==ParticleID::Kminus) ||
	      (id2==ParticleID::nu_taubar && id1==ParticleID::Kplus ))
	{allowed=true;}
      else if((id1==ParticleID::nu_tau    && id2==ParticleID::rhominus) ||
	      (id1==ParticleID::nu_taubar && id2==ParticleID::rhoplus ) ||
	      (id2==ParticleID::nu_tau    && id1==ParticleID::rhominus) ||
	      (id2==ParticleID::nu_taubar && id1==ParticleID::rhoplus ))
	{allowed=true;}
      else if((id1==ParticleID::nu_tau    && id2==ParticleID::Kstarminus) ||
	      (id1==ParticleID::nu_taubar && id2==ParticleID::Kstarplus ) ||
	      (id2==ParticleID::nu_tau    && id1==ParticleID::Kstarminus) ||
	      (id2==ParticleID::nu_taubar && id1==ParticleID::Kstarplus ))
	{allowed=true;}
    }
  return allowed;
}
  
ParticleVector Tau1MesonDecayer::decay(const DecayMode & dm,
				       const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  // perform the decay
  int id;
  unsigned int imode=0;
  for(unsigned int ix=0;ix<children.size();++ix)
    {
      id=children[ix]->id();
      if     (id==ParticleID::piminus   || id==ParticleID::piplus    ){imode=0;}
      else if(id==ParticleID::Kminus    || id==ParticleID::Kplus     ){imode=1;}
      else if(id==ParticleID::rhominus  || id==ParticleID::rhoplus   ){imode=2;}
      else if(id==ParticleID::Kstarplus || id==ParticleID::Kstarminus){imode=3;}
    }
  generate(false,imode,parent,children);
  return children;
}
  
  
void Tau1MesonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _fpi << _fk  << _grho << _gkstar <<  _piwgt <<  _kwgt;
}
  
void Tau1MesonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _fpi >> _fk  >> _grho >> _gkstar >>  _piwgt >>  _kwgt;
}
  
ClassDescription<Tau1MesonDecayer> Tau1MesonDecayer::initTau1MesonDecayer;
// Definition of the static class description member.
  
void Tau1MesonDecayer::Init() {
  
  static ClassDocumentation<Tau1MesonDecayer> documentation
    ("The \\classname{Tau1MesonDecayer} class implements the hadronic "
     "currents for the decay of the tau to a single meson, eg pi, K, rho.");
   
  static Parameter<Tau1MesonDecayer,double> interfacePiWeight
    ("PiWeight",
     "Maximum weight for the integration of the tau -> nu pi decay.",
     &Tau1MesonDecayer::_piwgt, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau1MesonDecayer,double> interfaceKaonWeight
    ("KaonWeight",
     "Maximum weight for the integration of the tau -> nu Kaon decay.",
     &Tau1MesonDecayer::_kwgt, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau1MesonDecayer,double> interfaceRhoWeight
    ("RhoWeight",
     "Maximum weight for the integration of the tau -> nu rho decay.",
     &Tau1MesonDecayer::_rhowgt, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau1MesonDecayer,double> interfaceKstarWeight
    ("KstarWeight",
     "Maximum weight for the integration of the tau -> nu Kstar decay.",
     &Tau1MesonDecayer::_kstarwgt, 1.0, -1.0e12, 1.0e12,
     false, false, false);
  
  static Parameter<Tau1MesonDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &Tau1MesonDecayer::_fpi, MeV, 130.7*MeV, -1.0e12*MeV, 1.0e12*MeV,
     false, false, false);
  
  static Parameter<Tau1MesonDecayer,Energy> interfacefk
    ("fk",
     "The Kaon decay constant",
     &Tau1MesonDecayer::_fk, MeV, 159.8*MeV, -1.0e12*MeV, 1.0e12*MeV,
     false, false, false);
  
  static Parameter<Tau1MesonDecayer,Energy2> interfacegrho
    ("grho",
     "The rho meson decay constant.",
     &Tau1MesonDecayer::_grho, GeV2, 0.115276*GeV2, -1.0e12*GeV2, 1.0e12*GeV2,
     false, false, false);
  
  static Parameter<Tau1MesonDecayer,Energy2> interfacegkstar
    ("gkstar",
     "The kstar meson decay constant.",
     &Tau1MesonDecayer::_gkstar, GeV2, 0.115276*GeV2, -1.0e12*GeV2, 1.0e12*GeV2,
     false, false, false);
  
}
  
// calculate the hadronic current
vector<LorentzPolarizationVector> 
Tau1MesonDecayer::hadronCurrent(bool vertex, const int imode, const int ichan,
				const Particle & inpart,
				const ParticleVector & decay) const 
{
  vector<LorentzPolarizationVector> temp;
  LorentzPolarizationVector vtemp;
  // find the decaying meson and perform the calculation
  int imes=0,id;
  for(unsigned int ix=0;ix<decay.size();++ix)
    {
      id = decay[ix]->id();
      if(id!=ParticleID::nu_tau&&id!=ParticleID::nu_taubar){imes=ix;}
    }
  // pseudo-scalar meson decay
  if(imode==0||imode==1)
    {
      Energy fact;
      if(imode==0){fact=_fpi*sqrt(SM().CKM(0,0));}
      else{fact=_fk*sqrt(SM().CKM(0,1));}
      temp.push_back(fact*(decay[imes]->momentum()));
      // set up the spin information for the particle
      if(vertex)
	{
	  SpinPtr stemp = new_ptr(ScalarSpinInfo(decay[imes]->momentum(),true));
	  decay[imes]->spinInfo(stemp);
	}
    }
  else
    {
      Energy2 fact;
      if(imode==2){fact=_grho*sqrt(2.*SM().CKM(0,0));}
      else{fact=_gkstar*sqrt(2.*SM().CKM(0,1));}
      // set up the spin information for the particle
      VectorSpinPtr hwtemp;
      if(vertex)
	{
	  SpinPtr vtemp = new_ptr(VectorSpinInfo(decay[imes]->momentum(),true));
	  decay[imes]->spinInfo(vtemp);
	  hwtemp= dynamic_ptr_cast<VectorSpinPtr>(vtemp);
	}
      // calculate the vectors
      VectorWaveFunction wave=VectorWaveFunction(decay[imes]->momentum(),
					      decay[imes]->dataPtr(),-1,outgoing);
      for(int iy=-1;iy<2;++iy)
	{
	  temp.push_back(fact*wave.Wave());
	  if(vertex){hwtemp->setBasisState(iy,wave.Wave());}
	  if(iy<1){wave.reset(iy+1);}
	}
    }
  return temp;
}

}
