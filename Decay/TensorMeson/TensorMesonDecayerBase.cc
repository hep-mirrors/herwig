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
  vector<LorentzTensor> inten;
  // wave functions etc for the incoming particle
  RhoDMatrix rhoin(PDT::Spin2);rhoin.average();
  TensorWaveFunction(inten,rhoin,const_ptr_cast<tPPtr>(&inpart),incoming,
		     true,false,vertex);
  // calculate the decay tensor
  vector<LorentzTensor> tensor(decayTensor(vertex,ichan,inpart,decay));
  // work out the mapping for the tensor vector
  vector<unsigned int> constants(decay.size()+1),ihel(decay.size()+1);
  vector<PDT::Spin> ispin(decay.size());
  int itemp(1);
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      itemp*=ispin[ix];constants[ix]=itemp;
    }
  constants[decay.size()]=1;
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin2,ispin);
  for(unsigned int hhel=0;hhel<tensor.size();++hhel)
    {
      // map the index for the hadrons to a helicity state
      for(unsigned int ix=decay.size();ix>0;--ix)
	{ihel[ix]=(hhel%constants[ix-1])/constants[ix];}
      // loop over the helicities of the incoming vector meson
      for(ihel[0]=0;ihel[0]<5;++ihel[0])
	{newME(ihel)=inten[ihel[0]]*tensor[hhel];}
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}
}
