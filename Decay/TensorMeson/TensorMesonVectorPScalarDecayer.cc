// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorPScalarDecayer class.
//

#include "TensorMesonVectorPScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonVectorPScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::tcVectorSpinPtr;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::EpsFunction;


TensorMesonVectorPScalarDecayer::~TensorMesonVectorPScalarDecayer() {}

bool TensorMesonVectorPScalarDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(id0==_incoming[ix])
	{
	  if((id1==_outgoingV[ix]&&id2==_outgoingP[ix])||
	     (id2==_outgoingV[ix]&&id1==_outgoingP[ix])){allowed=true;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}

ParticleVector TensorMesonVectorPScalarDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  int imode=-1;
  int id=parent.id(),id1=children[0]->id(),id2=children[1]->id();
  unsigned int ix=0;
  do 
    {
      if(id==_incoming[ix])
	{
	  if((id1==_outgoingV[ix]&&id2==_outgoingP[ix])||
	     (id2==_outgoingV[ix]&&id1==_outgoingP[ix])){imode=ix;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  generate(false,imode,parent,children);
  return children;
}

void TensorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoingV << _outgoingP << _maxweight << _coupling;}

void TensorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> _coupling;}

ClassDescription<TensorMesonVectorPScalarDecayer> TensorMesonVectorPScalarDecayer::initTensorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void TensorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorPScalarDecayer> documentation
    ("The \\classname{TensorMesonVectorPScalarDecayer} class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceOutcomingV
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1particle",
     &TensorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceOutcomingP
    ("OutgoingScalar",
     "The PDG code for the outgoing pseudoscalar meson",
     &TensorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMesonVectorPScalarDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// the hadronic tensor
vector<LorentzTensor> 
TensorMesonVectorPScalarDecayer::decayTensor(const bool vertex,
					     const int imode,
					     const int, 
					     const Particle & inpart,
					     const ParticleVector & decay) const
{
  // storage for the tensor
  vector<LorentzTensor> temp;
  // work out which is the vector
  unsigned int ivec,isca;
  if(decay[0]->id()==_outgoingV[imode]){ivec=0;isca=1;}
  else{ivec=1;isca=0;}
  // set up the spin information for ther decay products
  vector<LorentzPolarizationVector> outvec;
  tcVectorSpinPtr vspin;
  if(vertex)
    {
      SpinPtr stemp=new_ptr(ScalarSpinInfo(decay[isca]->momentum(),true));
      decay[isca]->spinInfo(stemp);
      stemp=new_ptr(VectorSpinInfo(decay[ivec]->momentum(),true));
      vspin= dynamic_ptr_cast<tcVectorSpinPtr>(stemp);
      decay[ivec]->spinInfo(stemp);
    }
  // calculate the current
  VectorWaveFunction vwave(decay[ivec]->momentum(),decay[ivec]->dataPtr(),outgoing);
  LorentzPolarizationVector eps;
  for(int ix=-1;ix<2;++ix)
    {
      if(ix==0&&_outgoingV[imode]==ParticleID::gamma)
	{
	  if(vertex){vspin->setBasisState(ix,LorentzPolarizationVector());}
	  temp.push_back(LorentzTensor());
	}
      else
	{
	  vwave.reset(ix);
	  if(vertex){vspin->setBasisState(ix,vwave.Wave());}
	  eps=_coupling[imode]*EpsFunction::product(decay[ivec]->momentum(),vwave.Wave(),
						    decay[isca]->momentum());
	  temp.push_back(LorentzTensor(decay[isca]->momentum(),eps));
	}
    }
  // return the answer
  return temp;
}

}
