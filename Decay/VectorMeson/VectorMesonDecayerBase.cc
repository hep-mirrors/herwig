// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonDecayerBase class.
//
//  Author: Peter Richardson
//


#include "VectorMesonDecayerBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::RhoDMatrix;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

VectorMesonDecayerBase::~VectorMesonDecayerBase() {}
  
bool VectorMesonDecayerBase::accept(const DecayMode & dm) const {
  return false;
}
  
ParticleVector VectorMesonDecayerBase::decay(const DecayMode & dm,
					     const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  return children;
}
  
  
void VectorMesonDecayerBase::persistentOutput(PersistentOStream & os) const {
}
  
void VectorMesonDecayerBase::persistentInput(PersistentIStream & is, int) {
}
  
AbstractClassDescription<VectorMesonDecayerBase> VectorMesonDecayerBase::initVectorMesonDecayerBase;
// Definition of the static class description member.
  
void VectorMesonDecayerBase::Init() {
  
  static ClassDocumentation<VectorMesonDecayerBase> documentation
    ("The \\classname{VectorMesonDecayerBase} class is the base class for the "
     "Herwig++ decays of vector mesons.");
  
}

double VectorMesonDecayerBase::me2(bool vertex, 
				   const int imode, const int ichan,
				   const Particle & inpart,
				   const ParticleVector & decay) const
{
  // first check if the decaying particle has spin information
  tcVectorSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcVectorSpinPtr>(inpart.spinInfo());}
  vector<LorentzPolarizationVector> invec;
  RhoDMatrix rhoin(3);rhoin.average();
  // if the spin info object exists use it
  if(inspin&&inpart.spinInfo())
    {
      for(int ix=-1;ix<2;++ix){invec.push_back(inspin->getDecayBasisState(ix));}
      inspin->decay();
      rhoin = inspin->rhoMatrix();
    }
  // if has spin info but not the right type
  else if(inpart.spinInfo())
    {
      cerr << "wrong type of spin info for the incoming particle " 
	   << "in VectorMesonDecayerBase::me2()" << endl;
    }
  // if no spin info create it
  else
    {
      SpinPtr newspin=new_ptr(VectorSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcVectorSpinPtr>(newspin);
      inspin->decayed(true);
      VectorWaveFunction temp=VectorWaveFunction(inpart.momentum(),inpart.dataPtr(),
						 incoming);
      for(int ix=-1;ix<2;++ix)
	{
	  temp.reset(ix);
	  invec.push_back(temp.Wave());
	  inspin->setDecayState(ix,invec[ix+1]);
	}
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // calculate the decay current
  vector<LorentzPolarizationVector> current=decayCurrent(vertex,imode,ichan,
							 inpart,decay);
  // work out the mapping for the current vector
  vector<int> constants(decay.size()+1), ispin(decay.size()),ihel(decay.size()+1);
  int itemp=1;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      itemp*=ispin[ix];constants[ix]=itemp;
    }
  constants[decay.size()]=1;
  // compute the matrix element
  DecayMatrixElement newME(3,ispin);
  for(unsigned int hhel=0;hhel<current.size();++hhel)
    {
      // map the index for the hadrons to a helicity state
      for(unsigned int ix=decay.size();ix>0;--ix)
	{
	  ihel[ix]=(hhel%constants[ix-1])/constants[ix]-int(ispin[ix-1]/2);
	  if(ispin[ix-1]%2==0&&ihel[ix]>-0&&ispin[ix-1]!=0){++ihel[ix];}
	}
      // loop over the helicities of the incoming vector meson
      for(int thel=-1;thel<2;++thel)
	{
	  ihel[0]=thel;
	  newME(ihel)= invec[thel+1]*current[hhel];
	}
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}

// the hadronic currents (does nothing
vector<LorentzPolarizationVector> 
VectorMesonDecayerBase::decayCurrent(const bool vertex, 
				     const int imode, const int ichan,
				     const Particle & inpart,
				     const ParticleVector & decay) const
{vector<LorentzPolarizationVector> temp; return temp;}
}
