// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorPScalarDecayer class.
//

#include "VectorMesonVectorPScalarDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonVectorPScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::tcVectorSpinPtr;
using Helicity::VectorWaveFunction;
using Helicity::EpsFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

VectorMesonVectorPScalarDecayer::~VectorMesonVectorPScalarDecayer() {}

bool VectorMesonVectorPScalarDecayer::accept(const DecayMode & dm) const {
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

ParticleVector VectorMesonVectorPScalarDecayer::decay(const DecayMode & dm,
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


void VectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const 
{os << _incoming << _outgoingV << _outgoingP << _maxweight << _coupling;}

void VectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> _coupling;}

ClassDescription<VectorMesonVectorPScalarDecayer> VectorMesonVectorPScalarDecayer::initVectorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void VectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorPScalarDecayer> documentation
    ("The\\classname{VectorMesonVectorPScalarDecayer} class is designed for the "
     "decay of a vector meson to another vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1 particle",
     &VectorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceOutcomingPScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &VectorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonVectorPScalarDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// the hadronic currents 
vector<LorentzPolarizationVector>  
VectorMesonVectorPScalarDecayer::decayCurrent(const bool vertex, const int imode,
					      const int, 
					      const Particle & inpart,
					      const ParticleVector & decay) const
{
  // storage for the current
  vector<LorentzPolarizationVector> temp;
  // work out which of the decay products is the vector and which is the scalar
  unsigned int ivec=1,isca=0;
  bool photon=_outgoingV[imode]==ParticleID::gamma;
  if(decay[0]->id()==_outgoingV[imode]){ivec=0;isca=1;}
  tcVectorSpinPtr vecsp;
  // set up the spin information for the decay products
  if(vertex)
    {
      // scalar
      SpinPtr stemp=new_ptr(ScalarSpinInfo(decay[isca]->momentum(),true));
      decay[isca]->spinInfo(stemp);
      // vector
      stemp =new_ptr(VectorSpinInfo(decay[ivec]->momentum(),true));
      decay[ivec]->spinInfo(stemp);
    }
  // calculate the currents
  VectorWaveFunction vwave=VectorWaveFunction(decay[ivec]->momentum(),
					      decay[ivec]->dataPtr(),outgoing);
  for(int ix=-1;ix<2;++ix)
    {
      if(ix==0&&photon){temp.push_back(LorentzPolarizationVector());}
      else
	{
	  vwave.reset(ix);
	  temp.push_back(_coupling[imode]*
			 EpsFunction::product(inpart.momentum(),vwave.Wave(),
					      decay[ivec]->momentum()));
	}
    }
  return temp;
}

}
