// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarPScalarVectorDecayer class.
//

#include "PScalarPScalarVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarPScalarVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::outgoing;

PScalarPScalarVectorDecayer::~PScalarPScalarVectorDecayer() {}

bool PScalarPScalarVectorDecayer::accept(const DecayMode & dm) const {
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
	 ((_outgoingP[ix]==id1&&_outgoingV[ix]==id2)||
	  (_outgoingP[ix]==id2&&_outgoingV[ix]==id1))){allowed=true;}
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector PScalarPScalarVectorDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
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
	  if((id1==_outgoingP[ix]&&id2==_outgoingV[ix])||
	     (id2==_outgoingP[ix]&&id1==_outgoingV[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  return generate(false,false,imode,parent);
}


void PScalarPScalarVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingP << _outgoingV << _maxweight;
}

void PScalarPScalarVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingP >> _outgoingV >> _maxweight;
}

ClassDescription<PScalarPScalarVectorDecayer> PScalarPScalarVectorDecayer::initPScalarPScalarVectorDecayer;
// Definition of the static class description member.

void PScalarPScalarVectorDecayer::Init() {

  static ClassDocumentation<PScalarPScalarVectorDecayer> documentation
    ("The \\classname{PScalarPScalarVectorDecayer} class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarPScalarVectorDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceOutgoingScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing pseudoscalar meson",
     &PScalarPScalarVectorDecayer::_outgoingP,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceOutgoingVector
    ("OutgoingVector",
     "The PDG code for the outgoing vector meson",
     &PScalarPScalarVectorDecayer::_outgoingV,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarPScalarVectorDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarPScalarVectorDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}


double PScalarPScalarVectorDecayer::me2(bool vertex, const int ichan,
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
	   << "in PScalarPScalarVectorDecayer::me2()" << endl;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // set up the spin info for the outgoing particles
  tcVectorSpinPtr outspin;
  if(vertex)
    {
      SpinPtr temp=new_ptr(ScalarSpinInfo(decay[0]->momentum(),true));
      decay[0]->spinInfo(temp);
      temp==new_ptr(VectorSpinInfo(decay[1]->momentum(),true));
      outspin= dynamic_ptr_cast<tcVectorSpinPtr>(temp);
      decay[1]->spinInfo(temp);
    }
  // setup the vector wavefunction
  VectorWaveFunction vtemp=VectorWaveFunction(decay[1]->momentum(),
					      decay[1]->dataPtr(),outgoing);
  // calculate the matrix element
  DecayMatrixElement newME(1,1,3);
  Lorentz5Momentum psum = inpart.momentum()+decay[0]->momentum();
  for(int ix=-1;ix<2;++ix)
    {
      vtemp.reset(ix);
      newME(0,0,ix)=_coupling[imode()]*(vtemp.Wave()*psum);
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  double me=newME.contract(rhoin).real();
  return me;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarPScalarVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
	  if((id1==_outgoingP[ix]&&id2==_outgoingV[ix])||
	     (id2==_outgoingP[ix]&&id1==_outgoingV[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  coupling=_coupling[imode];
  itype=10;
  return id1==_outgoingP[imode]&&id2==_outgoingV[imode];
}
}
