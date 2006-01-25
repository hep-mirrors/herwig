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
//#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::RhoDMatrix;
using Helicity::VectorWaveFunction;
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
    ("The VectorMesonDecayerBase class is the base class for the "
     "Herwig++ decays of vector mesons.");
  
}

double VectorMesonDecayerBase::me2(bool vertex, const int ichan,
				   const Particle & inpart,
				   const ParticleVector & decay) const
{
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // calculate the decay current
  vector<LorentzPolarizationVector> current=decayCurrent(vertex,ichan,
							 inpart,decay);
  // work out the mapping for the current vector
  vector<unsigned int> constants(decay.size()+1),ihel(decay.size()+1);
  vector<PDT::Spin> ispin(decay.size());
  int itemp=1;
  for(int ix=int(decay.size()-1);ix>=0;--ix)
    {
      ispin[ix]=decay[ix]->data().iSpin();
      itemp*=ispin[ix];constants[ix]=itemp;
    }
  constants[decay.size()]=1;
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,ispin);
  for(unsigned int hhel=0;hhel<current.size();++hhel)
    {
      // map the index for the hadrons to a helicity state
      for(unsigned int ix=decay.size();ix>0;--ix)
	{ihel[ix]=(hhel%constants[ix-1])/constants[ix];}
      // loop over the helicities of the incoming vector meson
      for(ihel[0]=0;ihel[0]<3;++ihel[0])
	{newME(ihel)=invec[ihel[0]]*current[hhel];}
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}

void VectorMesonDecayerBase::dataBaseOutput(ofstream & os,bool header) const
{DecayIntegrator::dataBaseOutput(os,header);}

}
