// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWZDecayer class.
//

#include "SMWZDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SMWZDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::tcFermionSpinPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::RhoDMatrix;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

SMWZDecayer::~SMWZDecayer() {}

bool SMWZDecayer::accept(const DecayMode & dm) const {
  bool allowed(false);
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  // Z to quarks or leptons
  if(id0==ParticleID::Z0)
    {if(id1==-id2&&(abs(id1)<=5||(abs(id1)>=11&&abs(id1)<=16))){allowed=true;}}
  // W to quarks
  else if(abs(id0)==ParticleID::Wplus)
    {
      int idd,idu;
      if(abs(id1)%2==1&&abs(id2)%2==0)
	{
	  idd=abs(id1);idu=abs(id2);
	  if((id1<0&&id2>0&&id0==ParticleID::Wminus)||
	     (id1>0&&id2<0&&id0==ParticleID::Wplus)){return allowed;}
	}
      else if(abs(id1)%2==0&&abs(id2)%2==1)
	{
	  idd=abs(id2);idu=abs(id1);
	  if((id2<0&&id1>0&&id0==ParticleID::Wminus)||
	     (id2>0&&id1<0&&id0==ParticleID::Wplus)){return allowed;}
	}
      else{return allowed;}
      if(idd<6&&idu<6){allowed=true;}
      else if(idd==idu-1&&idd>=11&&idd<=15){allowed=true;}
    }
  return allowed;
}

ParticleVector SMWZDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // id's of the decaying particles
  int id0=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  int imode=-1;
  if(id0==ParticleID::Z0)
    {
      if(abs(id1)<6){imode=abs(id1)-1;}
      else if(abs(id1)>=11&&abs(id1)<=16){imode=abs(id1)-6;}
    }
  else if(abs(id0)==ParticleID::Wplus)
    {
      int idd,idu;
      if(abs(id1)%2==1){idd=abs(id1);idu=abs(id2);}
      else{idd=abs(id2);idu=abs(id1);}
      if(idd<=5){imode=idd+idu/2+9;}
      else{imode=(idd-1)/2+12;}
    }
  bool cc = parent.id()==ParticleID::Wminus;
  return generate(false,cc,imode,parent);
}


void SMWZDecayer::persistentOutput(PersistentOStream & os) const {
  os << _Wvertex << _Zvertex << _Zquarkwgt << _Wquarkwgt << _Zleptonwgt << _Wleptonwgt;
}

void SMWZDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _Wvertex >> _Zvertex >> _Zquarkwgt >> _Wquarkwgt >> _Zleptonwgt >> _Wleptonwgt;
}

ClassDescription<SMWZDecayer> SMWZDecayer::initSMWZDecayer;
// Definition of the static class description member.

void SMWZDecayer::Init() {

  static ClassDocumentation<SMWZDecayer> documentation
    ("The \\classname{SMWZDecayer} class is the implementation of the decay"
     " of the W and Z bosons to the Standard Model fermions.");

  static ParVector<SMWZDecayer,double> interfaceZquarkMax
    ("ZquarkMax",
     "The maximum weight for the decay of the Z to quarks",
     &SMWZDecayer::_Zquarkwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZDecayer,double> interfaceWquarkMax
    ("WquarkMax",
     "The maximum weight for the decay of the W to quarks",
     &SMWZDecayer::_Wquarkwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZDecayer,double> interfaceZleptonMax
    ("ZleptonMax",
     "The maximum weight for the decay of the Z to leptons",
     &SMWZDecayer::_Zleptonwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZDecayer,double> interfaceWleptonMax
    ("WleptonMax",
     "The maximum weight for the decay of the W to leptons",
     &SMWZDecayer::_Wleptonwgt,
     0, 0, 0, -10000, 10000, false, false, true);

}


// return the matrix element squared
double SMWZDecayer::me2(bool vertex, const int ichan, const Particle & inpart,
			    const ParticleVector & decay) const 
{
  // check if the incoming particle has a spin info 
  tcVectorSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcVectorSpinPtr>(inpart.spinInfo());}
  RhoDMatrix rhoin(3);rhoin.average();
  VectorWaveFunction inwave[3];
  // if the spin info object exists use it
  if(inspin&&inpart.spinInfo())
    {
      for(int ix=-1;ix<2;++ix)
	{inwave[ix+1]=VectorWaveFunction(inpart.momentum(),inpart.dataPtr(),
					 inspin->getDecayBasisState(ix),incoming);}
      inspin->decay();
      rhoin = inspin->rhoMatrix();
    }
  else
    {
      // if has spin info but not the right type issue warning and throw away
      if(inpart.spinInfo())
	{throw DecayIntegratorError() << "Wrong type of spin info for the "
				      << "incoming particle in SMWZDecayer::me2()" 
				      << Exception::warning;}
      SpinPtr newspin=new_ptr(VectorSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcVectorSpinPtr>(newspin);
      inspin->decayed(true);
      VectorWaveFunction temp=VectorWaveFunction(inpart.momentum(),inpart.dataPtr(),
						 incoming);
      for(int ix=-1;ix<2;++ix)
	{
	  temp.reset(ix);
	  inwave[ix+1]=temp;
	  inspin->setDecayState(ix,temp.Wave());
	}
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // construct the spinors for the outgoing particles
  int iferm,ianti;
  if(decay[0]->id()<0){iferm=1;ianti=0;}
  else{iferm=0;ianti=1;}
  SpinorWaveFunction awave = SpinorWaveFunction(decay[ianti]->momentum(),
						decay[ianti]->dataPtr(),outgoing);
  SpinorBarWaveFunction fwave = SpinorBarWaveFunction(decay[iferm]->momentum(),
						      decay[iferm]->dataPtr(),outgoing);
  // spin info for the outgoing particles
  FermionSpinPtr ferm,anti;
  if(vertex)
    {
      ferm = new_ptr(FermionSpinInfo(decay[iferm]->momentum(),true));
      anti = new_ptr(FermionSpinInfo(decay[ianti]->momentum(),true));
      decay[iferm]->spinInfo(ferm);
      decay[ianti]->spinInfo(anti);
    }
  // compute the matrix element
  DecayMatrixElement newme(3,2,2);
  Energy2 scale=inpart.mass()*inpart.mass();
  for(int ifm=-1;ifm<2;ifm+=2)
    {
      fwave.reset(ifm);
      if(vertex){ferm->setBasisState(ifm,fwave.Wave().bar());}
      for(int ia=-1;ia<2;ia+=2)
	{
	  awave.reset(ia);
	  if(vertex&&ifm==-1){anti->setBasisState(ia,awave.Wave());}
	  for(int vhel=-1;vhel<2;++vhel)
	    {
	      if(inpart.id()==ParticleID::Z0)
		{
		  if(iferm>ianti){newme(vhel,ia,ifm)=
		      _Zvertex->evaluate(scale,awave,fwave,inwave[vhel+1]);}
		  else{newme(vhel,ifm,ia)=
		      _Zvertex->evaluate(scale,awave,fwave,inwave[vhel+1]);}
		}
	      else
		{
		  if(iferm>ianti){newme(vhel,ia,ifm)=
		      _Wvertex->evaluate(scale,awave,fwave,inwave[vhel+1]);}
		  else{newme(vhel,ifm,ia)=
		      _Wvertex->evaluate(scale,awave,fwave,inwave[vhel+1]);}
		}
	    }
	}
    }
  ME(newme);
  double output=(newme.contract(rhoin)).real()/scale;
  if(abs(decay[0]->id())<=6){output*=3.;}
  return output;
}

}
