// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarLeptonNeutrinoDecayer class.
//

#include "PScalarLeptonNeutrinoDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarLeptonNeutrinoDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::SpinorWaveFunction;
using Herwig::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::tcFermionSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::defaultDRep;
using ThePEG::Helicity::HELASDRep;

PScalarLeptonNeutrinoDecayer::~PScalarLeptonNeutrinoDecayer() {}

bool PScalarLeptonNeutrinoDecayer::accept(const DecayMode & dm) const {
  // must be two outgoing particles
  bool allowed=false;
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0=dm.parent()->id(),id0bar;
  if(dm.parent()->CC()){id0bar=-id0;}
  else{id0bar=id0;}
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id,ilep=4;
  for(;pit!=dm.products().end();++pit)
    {
      id=abs((**pit).id());
      if(id>=11&&id<=16){if(id%2==0){ilep=(id-10)/2;}}
    }
  // is the mode allowed
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if((id0   ==_incoming[ix]||id0bar==_incoming[ix])&&ilep<=_leptons[ix])
	{allowed=true;return allowed;}
    }
  return allowed;
}

ParticleVector PScalarLeptonNeutrinoDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // find the id of the leptons
  int ilep=4,id0=parent.id(),id0bar,id;
  ParticleMSet::const_iterator pit = dm.products().begin();
  for(;pit!=dm.products().end();++pit)
    {
      id=abs((**pit).id());
      if(id>=11&&id<=16){if(id%2==0){ilep=(id-10)/2;}}
    }
  // id's for the conjugate mode
  if(getParticleData(id0)->CC()){id0bar=-id0;}
  else{id0bar=id0;}
  // find the channel we need
  int ichan=-1;
  bool found=false; unsigned int ix=0;
  do
    {
      if(id0   ==_incoming[ix]||id0bar==_incoming[ix])
	{found=true;ichan+=ilep;}
      else{ichan+=_leptons[ix];}
      ++ix;
    }
  while (!found&&ix<_incoming.size());
  if(!found)
    {throw DecayIntegratorError() << "Unknown decay mode in PScalarLeptonNeutrino"
				  << "Decayer::decay" << Exception::abortnow;}
  // generate the decay
  bool cc = (id0==-_incoming[ix]);
  return generate(true,cc,ichan,parent);
}


void PScalarLeptonNeutrinoDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _decayconstant << _leptons << _maxweighte << _maxweightmu 
     << _maxweighttau << _GF;
}

void PScalarLeptonNeutrinoDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _decayconstant >> _leptons >> _maxweighte >> _maxweightmu 
     >> _maxweighttau >> _GF;
}

ClassDescription<PScalarLeptonNeutrinoDecayer> PScalarLeptonNeutrinoDecayer::initPScalarLeptonNeutrinoDecayer;
// Definition of the static class description member.

void PScalarLeptonNeutrinoDecayer::Init() {

  static ClassDocumentation<PScalarLeptonNeutrinoDecayer> documentation
    ("The \\classname{PScalarLeptonNeutrinoDecayer} class is the base class for"
     " the implementation of leptonic decay of a pseudoscalar meson to a lepton"
     "and a neutrino.");

  static ParVector<PScalarLeptonNeutrinoDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarLeptonNeutrinoDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,int> interfaceLepton
    ("Leptons",
     "The allowed outgoing leptons, 1 is only electrons, 2 is electrons and muons and "
     "3 is all lepton flavours",
     &PScalarLeptonNeutrinoDecayer::_leptons,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,double> interfaceMaxWeighte
    ("MaxWeightElectron",
     "The maximum weight for the electron decay mode",
     &PScalarLeptonNeutrinoDecayer::_maxweighte,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,double> interfaceMaxWeightmu
    ("MaxWeightMuon",
     "The maximum weight for the muon decay mode",
     &PScalarLeptonNeutrinoDecayer::_maxweightmu,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarLeptonNeutrinoDecayer,double> interfaceMaxWeighttau
    ("MaxWeightTau",
     "The maximum weight for the tau decay mode",
     &PScalarLeptonNeutrinoDecayer::_maxweighttau,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<PScalarLeptonNeutrinoDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &PScalarLeptonNeutrinoDecayer::_GF, 1./GeV2, 1.16639E-5/GeV2,
     -1.0e12*1./GeV2, 1.0e12*1./GeV2,
     false, false, false);

  static ParVector<PScalarLeptonNeutrinoDecayer,Energy> interfaceDecayConstant
    ("DecayConstant",
     "The decay constant for the incoming pseudoscaalr meson.",
     &PScalarLeptonNeutrinoDecayer::_decayconstant,
     0, 0, 0, 0*GeV, 10*GeV, false, false, true);

}

double PScalarLeptonNeutrinoDecayer::me2(bool vertex, const int ichan,
					 const Particle & inpart,
					 const ParticleVector & decay) const
{
  int icoup=0,id=abs(inpart.id()),idferm=0;
  // work out which decay constant to use
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {if(id==abs(_incoming[ix])){icoup=ix;}}
  // check if the decay particle has spin info 
  tcScalarSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcScalarSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  if(inspin)
    {inspin->decayed(true);}
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for th incoming particle"
				  << " in SemiLeptonicScalarScalarDecayer::me2()" 
				  << Exception::abortnow;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // find the particles
  unsigned int iferm=0,ianti=0;
  for(unsigned int ix=0;ix<decay.size();++ix)
    {
      id=decay[ix]->id();
      if(id<=-11&&id>=-16){ianti=ix;}
      else if(id>=11&&id<=16){iferm=ix;idferm=id;}
    }
  // construct the spininfo's of the outgoing particles
  FermionSpinPtr fspin,aspin;SpinPtr temp;
  if(vertex)
    {
      temp=new_ptr(FermionSpinInfo(decay[iferm]->momentum(),true));
      decay[iferm]->spinInfo(temp);
      fspin=dynamic_ptr_cast<FermionSpinPtr>(temp);
      temp=new_ptr(FermionSpinInfo(decay[ianti]->momentum(),true));
      decay[ianti]->spinInfo(temp);
      aspin=dynamic_ptr_cast<FermionSpinPtr>(temp);
    }
  // spinors for the lepton and neutrino
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wbar;
  SpinorBarWaveFunction fwave=SpinorBarWaveFunction(decay[iferm]->momentum(),
						    decay[iferm]->dataPtr(),outgoing);
  SpinorWaveFunction awave=SpinorWaveFunction(decay[ianti]->momentum(),
					      decay[ianti]->dataPtr(),outgoing);
  for(int ix=-1;ix<2;ix+=2)
    {
      awave.reset(ix);wave.push_back(awave.Wave());
      if(vertex){aspin->setBasisState(ix,wave[(ix+1)/2]);}
      fwave.reset(ix);wbar.push_back(fwave.Wave());
      if(vertex){fspin->setBasisState(ix,wbar[(ix+1)/2].bar());}

    }
  // the prefactor
  double pre;
  if(idferm%2==0){pre=decay[ianti]->mass();}
  else{pre=decay[iferm]->mass();}
  pre*=_decayconstant[icoup]*_GF;
  // compute the matrix element
  vector<int> ispin(2,2);
  DecayMatrixElement newME(1,ispin);
  Complex ii(0.,1.),term;
  ispin.resize(3);ispin[0]=0;
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(unsigned int iy=0;iy<2;++iy)
	{
	  if(idferm%2==0)
	    {
	      if(defaultDRep==HELASDRep)
		{term=2.*(wbar[iy].s3()*wave[ix].s3()+wbar[iy].s4()*wave[ix].s4());}
	      else
		{term=(wbar[iy].s1()+wbar[iy].s3())*(wave[ix].s1()+wave[ix].s3())
		     +(wbar[iy].s2()+wbar[iy].s4())*(wave[ix].s2()+wave[ix].s4());}
	    }
	  else
	    {
	      if(defaultDRep==HELASDRep)
		{term=2.*(wbar[iy].s1()*wave[ix].s1()+wbar[iy].s2()*wave[ix].s2());}
	      else
		{term=(wbar[iy].s1()-wbar[iy].s3())*(wave[ix].s1()-wave[ix].s3())
		     +(wbar[iy].s2()-wbar[iy].s4())*(wave[ix].s2()-wave[ix].s4());}
	    }
	  ispin[iferm+1]=2*iy-1;
	  ispin[ianti+1]=2*ix-1;
	  newME(ispin)=pre*term;
	}
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  double me=newME.contract(rhoin).real();
  return 0.5*me/inpart.mass()/inpart.mass();
}

}
