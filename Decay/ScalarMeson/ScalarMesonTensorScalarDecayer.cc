// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonTensorScalarDecayer class.
//

#include "ScalarMesonTensorScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarMesonTensorScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/TensorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzTensor;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::tcTensorSpinPtr;
using ThePEG::Helicity::TensorSpinInfo;
using Herwig::Helicity::TensorWaveFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::outgoing;


ScalarMesonTensorScalarDecayer::~ScalarMesonTensorScalarDecayer() {}

bool ScalarMesonTensorScalarDecayer::accept(const DecayMode & dm) const {
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
	 ((_outgoingT[ix]==id1&&_outgoingS[ix]==id2)||
	  (_outgoingT[ix]==id2&&_outgoingS[ix]==id1))){allowed=true;}
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector ScalarMesonTensorScalarDecayer::decay(const DecayMode & dm,
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
	  if((id1==_outgoingT[ix]&&id2==_outgoingS[ix])||
	     (id2==_outgoingT[ix]&&id1==_outgoingS[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  return generate(false,false,imode,parent);
}


void ScalarMesonTensorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingT << _outgoingS << _maxweight;
}

void ScalarMesonTensorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingT >> _outgoingS >> _maxweight;
}

ClassDescription<ScalarMesonTensorScalarDecayer> ScalarMesonTensorScalarDecayer::initScalarMesonTensorScalarDecayer;
// Definition of the static class description member.

void ScalarMesonTensorScalarDecayer::Init() {

  static ClassDocumentation<ScalarMesonTensorScalarDecayer> documentation
    ("The \\classname{ScalarMesonTensorScalarDecayer} class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarMesonTensorScalarDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceOutcomingT
    ("OutgoingTensor",
     "The PDG code for the outgoing tensor",
     &ScalarMesonTensorScalarDecayer::_outgoingT,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceOutcomingS
    ("OutgoingScalar",
     "The PDG code for the outgoing scalar",
     &ScalarMesonTensorScalarDecayer::_outgoingS,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarMesonTensorScalarDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarMesonTensorScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);
}


double ScalarMesonTensorScalarDecayer::me2(bool vertex, const int ichan,
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
	   << "in PScalarVectorVectorDecayer::me2()" << endl;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // set up the spin info for the outgoing particles
  tcTensorSpinPtr tenspin;
  if(vertex)
    {
      SpinPtr temp=new_ptr(TensorSpinInfo(decay[0]->momentum(),true));
      tenspin = dynamic_ptr_cast<tcTensorSpinPtr>(temp);
      decay[0]->spinInfo(temp);
      temp = new_ptr(ScalarSpinInfo(decay[1]->momentum(),true));
      decay[1]->spinInfo(temp);
    }
  // wavefunction for the outgoing tensor
  TensorWaveFunction temp=TensorWaveFunction(decay[0]->momentum(),
					     decay[0]->dataPtr(),outgoing);
  // calculate the matrix element
  DecayMatrixElement newME(1,5,1);
  Complex fact = _coupling[imode()]/inpart.mass();
  LorentzPolarizationVector vtemp;
  for(int ix=-2;ix<3;++ix)
    {
      temp.reset(ix);
      if(vertex){tenspin->setBasisState(ix,temp.Wave());}
      vtemp = temp.Wave()*(inpart.momentum()); 
      newME(0,ix,0) = fact*decay[1]->momentum()*vtemp;
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  double me=newME.contract(rhoin).real();
  return me;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarMesonTensorScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
	  if((id1==_outgoingT[ix]&&id2==_outgoingS[ix])||
	     (id2==_outgoingT[ix]&&id1==_outgoingS[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  coupling=_coupling[imode]*dm.parent()->mass();
  itype = 11;
  return id1==_outgoingT[imode]&&id2==_outgoingS[imode];
}

}
