// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2MesonDecayer class.
//

#include "VectorMeson2MesonDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::tcSpinInfoPtr;

VectorMeson2MesonDecayer::~VectorMeson2MesonDecayer() {}
  
bool VectorMeson2MesonDecayer::accept(const DecayMode & dm) const {
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
  
ParticleVector VectorMeson2MesonDecayer::decay(const DecayMode & dm,
					       const Particle & parent) const 
{
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
   
void VectorMeson2MesonDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;}
  
void VectorMeson2MesonDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;}
  
ClassDescription<VectorMeson2MesonDecayer> VectorMeson2MesonDecayer::initVectorMeson2MesonDecayer;
  // Definition of the static class description member.

void VectorMeson2MesonDecayer::Init() {
  
  static ClassDocumentation<VectorMeson2MesonDecayer> documentation
    ("The \\classname{VectorMeson2MesonDecayer} class is designed to implement "
     "the decay of vector mesons to 2 scalar mesons via a current which is the "
     "difference of the momenta of the two scalars. The order of the scalar meson "
     "momenta does not matter as it only changes the sign of the matrix element.");

  static ParVector<VectorMeson2MesonDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson2MesonDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &VectorMeson2MesonDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &VectorMeson2MesonDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMeson2MesonDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMeson2MesonDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);
  
}

// the hadronic currents    
vector<LorentzPolarizationVector>  
VectorMeson2MesonDecayer::decayCurrent(const bool vertex, const int imode, const int, 
				       const Particle & inpart,
				       const ParticleVector & decay) const
{
  // storage for the current
  vector<LorentzPolarizationVector> temp;
  // setup the spininfomation for the decay products
  if(vertex)
    {
      for(unsigned int ix=0;ix<decay.size();++ix)
	{
	  SpinPtr stemp= new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true));
	  decay[ix]->spinInfo(stemp);
	}
    }
  // calculate the current
  Lorentz5Momentum ptemp=_coupling[imode]*(decay[0]->momentum()-decay[1]->momentum());
  temp.push_back(LorentzPolarizationVector(ptemp));
  return temp;
}

}
