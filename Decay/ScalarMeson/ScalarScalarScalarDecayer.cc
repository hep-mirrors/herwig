// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarScalarScalarDecayer class.
//

#include "ScalarScalarScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarScalarScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::tcScalarSpinPtr;

ScalarScalarScalarDecayer::~ScalarScalarScalarDecayer() {}

bool ScalarScalarScalarDecayer::accept(const DecayMode & dm) const {
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
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0&&
	 ((_outgoing1[ix]==id1&&_outgoing2[ix]==id2)||
	  (_outgoing1[ix]==id2&&_outgoing2[ix]==id1))){allowed=true;}
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector ScalarScalarScalarDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // workout which mode we are doing
  int imode=-1;
  int id=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id)
	{
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  return generate(false,false,imode,parent);
}


void ScalarScalarScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void ScalarScalarScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<ScalarScalarScalarDecayer> ScalarScalarScalarDecayer::initScalarScalarScalarDecayer;
// Definition of the static class description member.

void ScalarScalarScalarDecayer::Init() {

  static ClassDocumentation<ScalarScalarScalarDecayer> documentation
    ("The \\classname{ScalarScalarScalarDecayer} class is designed for the"
     " decay of a scalar meson to two scalar mesons including off-shell effects");

  static ParVector<ScalarScalarScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarScalarScalarDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &ScalarScalarScalarDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &ScalarScalarScalarDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,Energy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarScalarScalarDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarScalarScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

double ScalarScalarScalarDecayer::me2(bool vertex, const int ichan,
				   const Particle & inpart,
				   const ParticleVector & decay) const
{
  // check if the decaying particle has spin info
  tcScalarSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcScalarSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  if(inspin)
    {inspin->decayed(true);}
  else if(inpart.spinInfo())
    {cerr << "wrong type of spin info for the incoming particle " 
	   << "in ScalarScalarScalarDecayer::me2()" << endl;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // set up the spin info for the outgoing particles
  tcScalarSpinPtr outspin[2];
  if(vertex)
    {
      for(unsigned int ix=0;ix<2;++ix)
	{
	  SpinPtr temp=new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true));
	  outspin[ix]= dynamic_ptr_cast<tcScalarSpinPtr>(temp);
	  decay[ix]->spinInfo(temp);
	}
    }
  // now compute the matrix element
  DecayMatrixElement newME(1,1,1);
  double fact = _coupling[imode()]/inpart.mass();
  newME(0,0,0)=fact;
  ME(newME);
  return fact*fact;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarScalarScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
      if(_incoming[ix]==id)
	{
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  coupling=_coupling[imode]/dm.parent()->mass();
  itype = 6;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}

}
