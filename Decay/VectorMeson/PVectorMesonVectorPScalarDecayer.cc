// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PVectorMesonVectorPScalarDecayer class.
//

#include "PVectorMesonVectorPScalarDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PVectorMesonVectorPScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

PVectorMesonVectorPScalarDecayer::~PVectorMesonVectorPScalarDecayer() {}

bool PVectorMesonVectorPScalarDecayer::accept(const DecayMode & dm) const {
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

ParticleVector PVectorMesonVectorPScalarDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode=-1;
  int id=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
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
  bool cc=false;
  return generate(false,cc,imode,parent);
}


void PVectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoingV << _outgoingP << _maxweight << _coupling;}

void PVectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> _coupling;}

ClassDescription<PVectorMesonVectorPScalarDecayer> PVectorMesonVectorPScalarDecayer::initPVectorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void PVectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<PVectorMesonVectorPScalarDecayer> documentation
    ("The\\classname{PVectorMesonVectorPScalarDecayer} class is designed for the "
     "decay of a pseudovector meson to a vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PVectorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1 particle",
     &PVectorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceOutcomingPScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &PVectorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PVectorMesonVectorPScalarDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PVectorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}


// the hadronic currents 
vector<LorentzPolarizationVector>  
PVectorMesonVectorPScalarDecayer::decayCurrent(const bool vertex, const int, 
					      const Particle & inpart,
					      const ParticleVector & decay) const
{
  // storage for the current
  vector<LorentzPolarizationVector> temp;
  // work out which of the decay products is the vector and which is the scalar
  unsigned int ivec=1,isca=0;
  bool photon=_outgoingV[imode()]==ParticleID::gamma;
  if(decay[0]->id()==_outgoingV[imode()]){ivec=0;isca=1;}
  tVectorSpinPtr vecsp;
  // set up the spin information for the decay products
  if(vertex)
    {
      // scalar
      decay[isca]->spinInfo(new_ptr(ScalarSpinInfo(decay[isca]->momentum(),true)));
      // vector
      SpinPtr temp=new_ptr(VectorSpinInfo(decay[ivec]->momentum(),true));
      vecsp = dynamic_ptr_cast<tVectorSpinPtr>(temp);
      decay[ivec]->spinInfo(temp);
    }
  VectorWaveFunction vwave=VectorWaveFunction(decay[ivec]->momentum(),
					      decay[ivec]->dataPtr(),outgoing);
  // calculate the currents
  Energy2 p0dotpv=inpart.momentum()*decay[ivec]->momentum();
  Complex epsdot=0.;
  LorentzPolarizationVector output;
  for(int ix=-1;ix<2;++ix)
    {
      if(ix==0&&photon)
	{
	  temp.push_back(LorentzPolarizationVector());
	  if(vertex){vecsp->setBasisState(ix,LorentzPolarizationVector());}
	}
      else
	{
	  vwave.reset(ix);
	  //if(vertex){vecsp->setBasisState(ix,vwave.Wave());}
	  epsdot=vwave.Wave()*inpart.momentum();
	  output = p0dotpv*vwave.Wave()-epsdot*decay[ivec]->momentum();
	  output *=_coupling[imode()]/inpart.mass();
	  temp.push_back(output);
	}
    }
  return temp;
 }

bool PVectorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,
						     int & mecode,
						     double & coupling) const
{
  // work out which mode it is for the coupling
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
	  if((id1==_outgoingV[ix]&&id2==_outgoingP[ix])||
	     (id2==_outgoingV[ix]&&id1==_outgoingP[ix])){imode=ix;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling = _coupling[imode]*dm.parent()->mass();  
  mecode = 4;
  return id1==_outgoingV[imode]&&id2==_outgoingP[imode];
}
}
