// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMeson2PScalarDecayer class.
//

#include "TensorMeson2PScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMeson2PScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;

TensorMeson2PScalarDecayer::~TensorMeson2PScalarDecayer() {}

bool TensorMeson2PScalarDecayer::accept(const DecayMode & dm) const {
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

ParticleVector TensorMeson2PScalarDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  int imode=-1;
  unsigned int ix=0;
  int id=parent.id();
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
  return generate(false,false,imode,parent);
}


void TensorMeson2PScalarDecayer::persistentOutput(PersistentOStream & os) const 
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;}

void TensorMeson2PScalarDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;}

ClassDescription<TensorMeson2PScalarDecayer> TensorMeson2PScalarDecayer::initTensorMeson2PScalarDecayer;
// Definition of the static class description member.

void TensorMeson2PScalarDecayer::Init() {

  static ClassDocumentation<TensorMeson2PScalarDecayer> documentation
    ("The \\classname{TensorMeson2PScalarDecayer} class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMeson2PScalarDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMeson2PScalarDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMeson2PScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);
}

// the hadronic tensor
vector<LorentzTensor> 
TensorMeson2PScalarDecayer::decayTensor(const bool vertex,const int, 
					const Particle & inpart,
					const ParticleVector & decay) const
{
  // storage for the tensor
  vector<LorentzTensor> temp;
  // set up the spin information for the decay products
  if(vertex)
    {
      for(unsigned int ix=0;ix<decay.size();++ix)
	{
	  SpinPtr stemp= new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true));
	  decay[ix]->spinInfo(stemp);
	}
    }
  // calcluate the current
  temp.push_back(_coupling[imode()]/inpart.mass()*
		 LorentzTensor(decay[0]->momentum(),decay[1]->momentum()));
  return temp;
}

bool TensorMeson2PScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					       double & coupling) const
{
  int imode=-1;
  int id=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
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
  coupling=_coupling[imode]*dm.parent()->mass();
  mecode=7;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}
}
