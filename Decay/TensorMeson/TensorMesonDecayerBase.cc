// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonDecayerBase class.
//

#include "TensorMesonDecayerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonDecayerBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/TensorSpinInfo.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"


namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::tcTensorSpinPtr;
using ThePEG::Helicity::TensorSpinInfo;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::incoming;


TensorMesonDecayerBase::~TensorMesonDecayerBase() {}

bool TensorMesonDecayerBase::accept(const DecayMode & dm) const {
  return false;
}

ParticleVector TensorMesonDecayerBase::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  return children;
}


  void TensorMesonDecayerBase::persistentOutput(PersistentOStream & os) const {}

  void TensorMesonDecayerBase::persistentInput(PersistentIStream & is, int) {}

AbstractClassDescription<TensorMesonDecayerBase> TensorMesonDecayerBase::initTensorMesonDecayerBase;
// Definition of the static class description member.

void TensorMesonDecayerBase::Init() {

  static ClassDocumentation<TensorMesonDecayerBase> documentation
    ("The \\classname{TensorMesonDecayerBase} class is the base class"
     " for the implementation of the decays of tensor mesons.");

}

// matrix elememt for the process
double TensorMesonDecayerBase::me2(bool vertex, const int ichan,
				   const Particle & inpart,
				   const ParticleVector & decay) const
{
  // first check if the decaying particle has spin information
  tcTensorSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin=dynamic_ptr_cast<tcTensorSpinPtr>(inpart.spinInfo());}
  vector<LorentzTensor> inten;
  RhoDMatrix rhoin(5);rhoin.average();
  // if the spin info object exists use it
  if(inspin)
    {
      for(int ix=-2;ix<3;++ix){inten.push_back(inspin->getDecayBasisState(ix));}
      inspin->decay();
      rhoin=inspin->rhoMatrix();
    }
  else if(inpart.spinInfo())
    {
      throw TensorDecayerError() << "Wrong type of spinInfo for the incoming particle "
				 << "in TensorMesonDecayerBase::me2()" 
				 << Exception::abortnow;
    }
  // if no spinInfo then create it
  else
    { 
      SpinPtr newspin=new_ptr(TensorSpinInfo(inpart.momentum(),true));
      inspin =dynamic_ptr_cast<tcTensorSpinPtr>(newspin);
      inspin->decayed(true);
      TensorWaveFunction temp=TensorWaveFunction(inpart.momentum(),
						 inpart.dataPtr(),incoming);
      for(int ix=-2;ix<3;++ix)
	{
	  temp.reset(ix);
	  inten.push_back(temp.Wave());
	  inspin->setDecayState(ix,inten[ix+2]);
	}
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // calculate the decay tensor
  vector<LorentzTensor> tensor=decayTensor(vertex,ichan,inpart,decay);
  // work out the mapping for the tensor vector
  vector<int> constants(decay.size()+1),
    ispin(decay.size()),ihel(decay.size()+1);
  int itemp=1;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      itemp*=ispin[ix];constants[ix]=itemp;
    }
  constants[decay.size()]=1;
  // compute the matrix element
  DecayMatrixElement newME(5,ispin);
  for(unsigned int hhel=0;hhel<tensor.size();++hhel)
    {
      // map the index for the hadrons to a helicity state
      for(unsigned int ix=decay.size();ix>0;--ix)
	{
	  ihel[ix]=(hhel%constants[ix-1])/constants[ix]-int(ispin[ix-1]/2);
	  if(ispin[ix-1]%2==0&&ihel[ix]>-0&&ispin[ix-1]!=0){++ihel[ix];}
	}
      // loop over the helicities of the incoming vector meson
      for(int thel=-2;thel<3;++thel)
	{ihel[0]=thel;newME(ihel)= inten[thel+2]*tensor[hhel];}
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}
}
