// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorVectorDecayer class.
//

#include "PScalarVectorVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarVectorVectorDecayer.tcc"
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
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::EpsFunction;
using Herwig::Helicity::Direction;
using Herwig::Helicity::outgoing;

PScalarVectorVectorDecayer::~PScalarVectorVectorDecayer() {}

bool PScalarVectorVectorDecayer::accept(const DecayMode & dm) const {
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

ParticleVector PScalarVectorVectorDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  // workout which mode we are doing
  int imode=-1;
  int id=parent.id(),id1=children[0]->id(),id2=children[1]->id();
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
  generate(false,imode,parent,children);
  return children;
}


void PScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void PScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<PScalarVectorVectorDecayer> PScalarVectorVectorDecayer::initPScalarVectorVectorDecayer;
// Definition of the static class description member.

void PScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<PScalarVectorVectorDecayer> documentation
    ("The \\classname{PScalarVectorVectorDecayer} class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<PScalarVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarVectorVectorDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &PScalarVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &PScalarVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarVectorVectorDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarVectorVectorDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

double PScalarVectorVectorDecayer::me2(bool vertex, 
				   const int imode, const int ichan,
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
  tcVectorSpinPtr outspin[2];
  if(vertex)
    {
      for(unsigned int ix=0;ix<2;++ix)
	{
	  SpinPtr temp=new_ptr(VectorSpinInfo(decay[ix]->momentum(),true));
	  outspin[ix]= dynamic_ptr_cast<tcVectorSpinPtr>(temp);
	  decay[ix]->spinInfo(temp);
	}
    }
  // calculate the wavefunctions for the outgoing vectors
  LorentzPolarizationVector wave[3][2];
  for(unsigned int ix=0;ix<2;++ix)
    {
      VectorWaveFunction temp=VectorWaveFunction(decay[ix]->momentum(),
						 decay[ix]->dataPtr(),outgoing);
      for(int iy=-1;iy<2;++iy)
	{
	  if(iy==0&&decay[ix]->id()==ParticleID::gamma)
	    {wave[iy+1][ix]=LorentzPolarizationVector();}
	  else
	    {temp.reset(iy);wave[iy+1][ix]=temp.Wave();}
	  if(vertex){outspin[ix]->setBasisState(iy,wave[iy+1][ix]);}
	}
    }
  // now compute the matrix element
  DecayMatrixElement newME(1,3,3);
  for(unsigned int ix=0;ix<3;++ix)
    {
      for(unsigned int iy=0;iy<3;++iy)
	{
	  newME(0,ix-1,iy-1)=_coupling[imode]*
	    EpsFunction::product(wave[ix][0],decay[1]->momentum(),wave[iy][1])
	    *decay[0]->momentum();
	}
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  double me=newME.contract(rhoin).real();
  return me;
}

}
