// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorVectorDecayer class.
//

#include "TensorMesonVectorVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonVectorVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig{
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::outgoing;

TensorMesonVectorVectorDecayer::~TensorMesonVectorVectorDecayer() {}

bool TensorMesonVectorVectorDecayer::accept(const DecayMode & dm) const {
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
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){allowed=true;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}

ParticleVector TensorMesonVectorVectorDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  int imode=-1;
  int id=parent.id(),id1=children[0]->id(),id2=children[1]->id();
  unsigned int ix=0;
  do 
    {
      if(id==_incoming[ix])
	{
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  generate(false,imode,parent,children);
  return children;
}


void TensorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const 
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;}

void TensorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;}

ClassDescription<TensorMesonVectorVectorDecayer> TensorMesonVectorVectorDecayer::initTensorMesonVectorVectorDecayer;
// Definition of the static class description member.

void TensorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorVectorDecayer> documentation
    ("The \\classname{TensorMesonVectorVectorDecayer} class performs the"
     " decay of a tensor meson to two scalar mesons.");

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMesonVectorVectorDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &TensorMesonVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &TensorMesonVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMesonVectorVectorDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMesonVectorVectorDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// the hadronic tensor
vector<LorentzTensor> 
TensorMesonVectorVectorDecayer::decayTensor(const bool vertex,
					const int imode,
					const int, 
					const Particle & inpart,
					const ParticleVector & decay) const
{
  // storage for the tensor
  vector<LorentzTensor> temp;
  // set up the spin information for the decay products
  tcVectorSpinPtr vspin[2];
  if(vertex)
    {
      for(unsigned int ix=0;ix<decay.size();++ix)
	{
	  SpinPtr stemp=new_ptr(VectorSpinInfo(decay[ix]->momentum(),true));
	  decay[ix]->spinInfo(stemp);
	  vspin[ix]=dynamic_ptr_cast<tcVectorSpinPtr>(stemp);
	}
    }
  // compute the polarization vectors for the vectors
  VectorWaveFunction vwave[2]={VectorWaveFunction(decay[0]->momentum(),
						  decay[0]->dataPtr(),outgoing),
			       VectorWaveFunction(decay[1]->momentum(),
						  decay[1]->dataPtr(),outgoing)};
  LorentzPolarizationVector pol[2][3];
  int id[2]={decay[0]->id(),decay[1]->id()};
  for(unsigned int iy=0;iy<2;++iy)
    {
      for(int ix=-1;ix<2;++ix)
	{
	  if(id[iy]==ParticleID::gamma&&ix==0)
	    {pol[iy][ix+1]=LorentzPolarizationVector();}
	  else
	    {
	      vwave[iy].reset(ix);
	      pol[iy][ix+1]=vwave[iy].Wave();
	    }
	}
    }
  // compute some useful dot products etc
  Complex p1eps2[3],p2eps1[3];
  Energy2 p1p2=decay[0]->momentum()*decay[1]->momentum();
  for(unsigned int ix=0;ix<3;++ix)
    {
      p1eps2[ix]=pol[1][ix]*decay[0]->momentum();
      p2eps1[ix]=pol[0][ix]*decay[1]->momentum();
    }
  // compute some useful tensors to save CPU
  LorentzTensor tp1p2=LorentzTensor(decay[0]->momentum(),decay[1]->momentum())
                     +LorentzTensor(decay[1]->momentum(),decay[0]->momentum());
  LorentzTensor met(-1.,0. ,0. ,0.,0. ,-1.,0. ,0.,0. ,0. ,-1.,0.,0. ,0. ,0. ,1. );
  LorentzTensor tp1eps2[3],tp2eps1[3];
  for(unsigned int ix=0;ix<3;++ix)
    {
      tp1eps2[ix]=LorentzTensor(decay[0]->momentum(),pol[1][ix])
	         +LorentzTensor(pol[1][ix],decay[0]->momentum());
      tp2eps1[ix]=LorentzTensor(decay[1]->momentum(),pol[0][ix])
	         +LorentzTensor(pol[0][ix],decay[1]->momentum()); 
    }
  // main loop to compute the tensors
  Complex e1e2;
  LorentzTensor output;
  for(unsigned int ix=0;ix<3;++ix)
    {
      for(unsigned int iy=0;iy<3;++iy)
	{
	  e1e2=pol[0][ix]*pol[1][iy];
	  output=e1e2*tp1p2-p2eps1[ix]*tp1eps2[iy]-p1eps2[iy]*tp2eps1[ix]
	    +p1p2*( LorentzTensor(pol[0][ix],pol[1][iy])
		   +LorentzTensor(pol[1][iy],pol[0][ix]))
	    +(p2eps1[ix]*p1eps2[iy]-e1e2*p1p2)*met;
	  output*=_coupling[imode];
	  temp.push_back(output);
	}
    }
  // return the answer
  return temp;
}
}
