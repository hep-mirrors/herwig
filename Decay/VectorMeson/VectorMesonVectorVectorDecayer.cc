// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorVectorDecayer class.
//

#include "VectorMesonVectorVectorDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonVectorVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig{
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::tVectorSpinPtr;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

VectorMesonVectorVectorDecayer::~VectorMesonVectorVectorDecayer() {}

bool VectorMesonVectorVectorDecayer::accept(const DecayMode & dm) const {
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

ParticleVector VectorMesonVectorVectorDecayer::decay(const DecayMode & dm,
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
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  bool cc=false;
  return generate(false,cc,imode,parent);
}


void VectorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;}

void VectorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;}

ClassDescription<VectorMesonVectorVectorDecayer> VectorMesonVectorVectorDecayer::initVectorMesonVectorVectorDecayer;
// Definition of the static class description member.

void VectorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorVectorDecayer> documentation
    ("The\\classname{VectorMesonVectorVectorDecayer} class is designed for the "
     "decay of a vector meson to two vector particles, either photons or other "
     "vector mesons.");

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonVectorVectorDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceOutgoing1
    ("Outgoing1",
     "The PDG code for the first outgoing particle",
     &VectorMesonVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceOutgoing2
    ("Outgoing2",
     "The PDG code for the second outgoing particle",
     &VectorMesonVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonVectorVectorDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonVectorVectorDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// the hadronic currents
vector<LorentzPolarizationVector>  
VectorMesonVectorVectorDecayer::decayCurrent(const bool vertex, const int, 
					     const Particle & inpart,
					     const ParticleVector & decay) const
{
  // storage of the current
  vector<LorentzPolarizationVector> temp;
  tVectorSpinPtr vec[2];
  // set up the spin information for the decay products
  if(vertex)
    {
      SpinPtr stemp;
      // first particle
      stemp=new_ptr(VectorSpinInfo(decay[0]->momentum(),true));
      vec[0] = dynamic_ptr_cast<tVectorSpinPtr>(stemp);
      decay[0]->spinInfo(stemp);
      // second particle
      stemp=new_ptr(VectorSpinInfo(decay[1]->momentum(),true));
      vec[1] = dynamic_ptr_cast<tVectorSpinPtr>(stemp);
      decay[1]->spinInfo(stemp);
    }
  // calculate the polarization vectors
  VectorWaveFunction vwave[2]={VectorWaveFunction(decay[0]->momentum(),
						  decay[0]->dataPtr(),outgoing),
			       VectorWaveFunction(decay[1]->momentum(),
						  decay[1]->dataPtr(),outgoing)};
  vector<LorentzPolarizationVector> eps[2];eps[0].resize(3);eps[1].resize(3);
  for(unsigned int iout=0;iout<2;++iout)
    {
      if(decay[iout]->id()==ParticleID::gamma)
	{
	  for(int ix=-1;ix<2;ix+=2)
	    {
	      vwave[iout].reset(ix);
	      eps[iout][ix+1]=vwave[iout].Wave();
	    }
	  eps[iout][1]=LorentzPolarizationVector();
	  if(vertex){vec[iout]->setBasisState(1,LorentzPolarizationVector());}
	}
      else
	{
	  for(int ix=-1;ix<2;++ix)
	    {
	      vwave[iout].reset(ix);
	      eps[iout][ix+1]=vwave[iout].Wave();
	      vec[iout]->setBasisState(ix,eps[iout][ix+1]);
	    }
	}
    }
  // work out the dot product we need for the current
  Complex p1p2=(decay[0]->momentum())*(decay[1]->momentum());
  Complex p1eps2[3],p2eps1[3];
  Complex eps1eps2;
  for(unsigned int ix=0;ix<3;++ix)
    {
      p1eps2[ix]=eps[1][ix]*(decay[0]->momentum());
      p2eps1[ix]=eps[0][ix]*(decay[1]->momentum());
    }
  // now compute the current
  Lorentz5Momentum pdiff = decay[0]->momentum()-decay[1]->momentum();
  LorentzPolarizationVector curr;
  Complex m12 = decay[0]->mass()*decay[0]->mass();
  Complex m22 = decay[1]->mass()*decay[1]->mass();
  double fact = 2.*_coupling[imode()]/(inpart.mass()*inpart.mass()*inpart.mass());
  for(unsigned int ipol1=0;ipol1<3;++ipol1)
    {
      for(unsigned int ipol2=0;ipol2<3;++ipol2)
	{
	  eps1eps2=eps[0][ipol1]*eps[1][ipol2];
	  curr = p1eps2[ipol2]*p2eps1[ipol1]*pdiff
	    +p1eps2[ipol2]*m22*eps[0][ipol1]
	    -p2eps1[ipol1]*m12*eps[1][ipol2]
	    +eps1eps2*(
		       -p1p2*pdiff
		       +m12*decay[1]->momentum()
		       -m22*decay[0]->momentum());
	  curr *=fact;
	  temp.push_back(curr);
	}
    }
  return temp;
} 

bool VectorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling = _coupling[imode]; 
  mecode = 5;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode]; 
}

}
