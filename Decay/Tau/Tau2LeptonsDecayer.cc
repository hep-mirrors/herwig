// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tau2LeptonsDecayer class.
//
//  Author: Peter Richardson 
//

#include "Tau2LeptonsDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"

namespace Herwig {

using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::FermionSpinPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

Tau2LeptonsDecayer::~Tau2LeptonsDecayer() {}
  
bool Tau2LeptonsDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  // check the parent is tau+-/ and three decay products
  if(abs(dm.parent()->id())==15 && dm.products().size() == 3)
    {
      ParticleMSet::const_iterator pit  = dm.products().begin();
      ParticleMSet::const_iterator pend = dm.products().end();
      bool found[5]={false,false,false},particle[5]={false,false,false,false,false};
      int idtemp;
      for( ; pit!=pend;++pit)
	{
	  idtemp=(**pit).id();
	  if(abs(idtemp)==ParticleID::nu_tau)
	    {found[0]=true; if(idtemp>0){particle[0]=true;}}
	  else if(abs(idtemp)==ParticleID::eminus)
	    {found[1]=true; if(idtemp>0){particle[1]=true;}}
	  else if(abs(idtemp)==ParticleID::nu_e)
	    {found[2]=true; if(idtemp<0){particle[2]=true;}}
	  else if(abs(idtemp)==ParticleID::muminus)
	    {found[3]=true; if(idtemp>0){particle[3]=true;}}
	  else if(abs(idtemp)==ParticleID::nu_mu)
	    {found[4]=true; if(idtemp<0){particle[4]=true;}}
	}
      // check right particles
      // electron mode
      if(found[0]&&found[1]&&found[2]&&
	 ((dm.parent()->id()== 15&& particle[0]&& particle[1]&& particle[2])||
	  (dm.parent()->id()==-15&&!particle[0]&&!particle[1]&&!particle[2])))
	{allowed=true;}
      // muon mode
      else if(found[0]&&found[3]&&found[4]&&
	      ((dm.parent()->id()== 15&& particle[0]&& particle[3]&& particle[4])||
	       (dm.parent()->id()==-15&&!particle[0]&&!particle[3]&&!particle[4])))
	{allowed=true;}
    }
  return allowed;
}
  
ParticleVector Tau2LeptonsDecayer::decay(const DecayMode & dm,
					 const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  int imode=0;
  // decide if electron or muon mode
  for(unsigned int ix=0;ix<children.size();++ix)
    {if(abs(children[ix]->id())==ParticleID::muminus){++imode;}}
  // generate the decay
  generate(false,imode,parent,children);
  return children;
}
    
void Tau2LeptonsDecayer::persistentOutput(PersistentOStream & os) const 
{os << _electronwgt << _muonwgt;}
  
void Tau2LeptonsDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _electronwgt >> _muonwgt;}

ClassDescription<Tau2LeptonsDecayer> Tau2LeptonsDecayer::initTau2LeptonsDecayer;
// Definition of the static class description member.
  
void Tau2LeptonsDecayer::Init() {
    
  static ClassDocumentation<Tau2LeptonsDecayer> documentation
    ("The \\classname{Tau2LeptonsDecayer} class is designed to handle the "
     "leptonic decays of the tau.");

  static Parameter<Tau2LeptonsDecayer,double> interfaceElectronWeight
    ("ElectronWeight",
     "Maximum weight for the integration of the electron channel",
     &Tau2LeptonsDecayer::_electronwgt, 1.215E-9, -1.0e12, 1.0e12,
     false, false, false);

  static Parameter<Tau2LeptonsDecayer,double> interfaceMuonWeight
    ("MuonWeight",
     "Maximum weight for the integration of the muon channel",
     &Tau2LeptonsDecayer::_muonwgt, 1.065E-9, -1.0e12, 1.0e12,
     false, false, false);
  
}
  
 // current for the decay
vector<LorentzPolarizationVector> 
Tau2LeptonsDecayer::hadronCurrent(bool vertex,const int imode, const int ichan, 
				  const Particle & inpart,
				  const ParticleVector & decay) const
{
  // storage for the currents
  vector<LorentzPolarizationVector> temp;
  // locate the momenta of the lepton and neutrino
  int idtau=inpart.id();
  unsigned int inu=0,ilep=0;
  for(unsigned int ix=0;ix<decay.size();++ix)
    {
      if(abs(decay[ix]->id())==ParticleID::eminus||
	 abs(decay[ix]->id())==ParticleID::muminus){ilep=ix;}
      else if(abs(decay[ix]->id())==ParticleID::nu_e||
	      abs(decay[ix]->id())==ParticleID::nu_mu){inu=ix;}
    }
  // construct the spin information objects for the  decay products
  FermionSpinPtr lepspin,nuspin;
  if(vertex)
    {
      SpinPtr slep=new_ptr(FermionSpinInfo(decay[ilep]->momentum(),true));
      decay[ilep]->spinInfo(slep);
      lepspin=dynamic_ptr_cast<FermionSpinPtr>(slep);
      SpinPtr snu =new_ptr(FermionSpinInfo(decay[inu ]->momentum(),true));
      decay[inu ]->spinInfo(snu );
      nuspin =dynamic_ptr_cast<FermionSpinPtr>(snu);
    }
  // lepton wavefunctions for the different helicities
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  if(idtau>0)
    {
      SpinorWaveFunction nu =SpinorWaveFunction(decay[inu]->momentum(),
						decay[inu]->dataPtr(),outgoing);
      SpinorBarWaveFunction lep=SpinorBarWaveFunction(decay[ilep]->momentum(),
						      decay[ilep]->dataPtr(),outgoing);
      for(int ix=-1;ix<2;ix+=2)
	{
	  nu.reset(ix);wave.push_back(nu.Wave());
	  nuspin->setBasisState(ix,wave[(ix+1)/2]);
	  lep.reset(ix);wavebar.push_back(lep.Wave());
	  lepspin->setBasisState(ix,wavebar[(ix+1)/2].bar());
	}
    }
  else
    {
      SpinorWaveFunction lep=SpinorWaveFunction(decay[ilep]->momentum(),
						decay[ilep]->dataPtr(),outgoing);
      SpinorBarWaveFunction nu=SpinorBarWaveFunction(decay[inu]->momentum(),
						     decay[inu]->dataPtr(),outgoing);
      for(int ix=-1;ix<2;ix+=2)
	{
	  lep.reset(ix);wave.push_back(lep.Wave());
	  lepspin->setBasisState(ix,wave[(ix+1)/2]);
	  nu.reset(ix);wavebar.push_back(nu.Wave());
	  nuspin->setBasisState(ix,wavebar[(ix+1)/2].bar());
	}
    }
  // now compute the currents
  Complex vec[4], ii(0.,1.);
  temp.resize(4); int iloc=0;
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(unsigned int iy=0;iy<2;++iy)
	{
	  // calculate the current
	  if(wave[ix].Rep()==HaberDRep&&wavebar[iy].Rep()==HaberDRep)
	    {
	      Complex s2m4=wave[ix].s2()-wave[ix].s4();
	      Complex s1m3=wave[ix].s1()-wave[ix].s3();
	      vec[0] =   (-wavebar[iy].s1()*s2m4-wavebar[iy].s2()*s1m3
			  -wavebar[iy].s3()*s2m4-wavebar[iy].s4()*s1m3);
	      vec[1] =ii*(+wavebar[iy].s1()*s2m4-wavebar[iy].s2()*s1m3
			  +wavebar[iy].s3()*s2m4-wavebar[iy].s4()*s1m3);
	      vec[2] =   (-wavebar[iy].s1()*s1m3+wavebar[iy].s2()*s2m4
			  -wavebar[iy].s3()*s1m3+wavebar[iy].s4()*s2m4);
	      vec[3] =   (+wavebar[iy].s1()*s1m3+wavebar[iy].s2()*s2m4
			  +wavebar[iy].s3()*s1m3+wavebar[iy].s4()*s2m4);
	    }
	  else if(wave[ix].Rep()==HELASDRep&&wavebar[iy].Rep()==HELASDRep)
	    {
	      Complex s3s1=wavebar[iy].s3()*wave[ix].s1();
	      Complex s3s2=wavebar[iy].s3()*wave[ix].s2();
	      Complex s4s1=wavebar[iy].s4()*wave[ix].s1();
	      Complex s4s2=wavebar[iy].s4()*wave[ix].s2();
	      vec[0] =   -2.*(s3s2+s4s1);
	      vec[1] = ii*2.*(s3s2-s4s1);
	      vec[2] =   -2.*(s3s1-s4s2);
	      vec[3] =    2.*(s3s1+s4s2);
	    }
	  // location in the vector
	  if(inu<ilep)
	    {
	      if(idtau>0){iloc=2*ix+iy;}
	      else{iloc=2*iy+ix;}
	    }
	  else
	    {
	      if(idtau>0){iloc=2*iy+ix;}
	      else{iloc=2*ix+iy;}
	    }
	  // add it to the vector
	  temp[iloc]=LorentzPolarizationVector(vec[0],vec[1],vec[2],vec[3]);
	}
    }
  // return the answer
  return temp;
}

}
