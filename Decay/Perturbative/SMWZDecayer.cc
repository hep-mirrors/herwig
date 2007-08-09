// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWZDecayer class.
//

#include "SMWZDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMWZDecayer::SMWZDecayer() 
  : _zquarkwgt(5,0.),
    _zleptonwgt(6,0.),
    _wquarkwgt(6,0.),
    _wleptonwgt(3,0.)
{
   _zquarkwgt[0]  = 0.488029;
   _zquarkwgt[1]  = 0.378461;
   _zquarkwgt[2]  = 0.488019;
   _zquarkwgt[3]  = 0.378027;
   _zquarkwgt[4]  = 0.483207;
   _zleptonwgt[0] = 0.110709;
   _zleptonwgt[1] = 0.220276;
   _zleptonwgt[2] = 0.110708;
   _zleptonwgt[3] = 0.220276;
   _zleptonwgt[4] = 0.110458;
   _zleptonwgt[5] = 0.220276;
   _wquarkwgt[0]  = 1.01596;
   _wquarkwgt[1]  = 0.0537308;
   _wquarkwgt[2]  = 0.0538085;
   _wquarkwgt[3]  = 1.01377;
   _wquarkwgt[4]  = 1.45763e-05;
   _wquarkwgt[5]  = 0.0018143;
   _wleptonwgt[0] = 0.356594;
   _wleptonwgt[1] = 0.356593;
   _wleptonwgt[2] = 0.356333;
   // intermediates
   generateIntermediates(false);
}

void SMWZDecayer::doinit() throw(InitException) 
{
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(hwsm)
    {
      _wvertex = hwsm->vertexFFW();
      _zvertex = hwsm->vertexFFZ();
      // make sure they are initialized
      _wvertex->init();
      _zvertex->init();
    }
  else
    {throw InitException();}
  // now set up the decay modes
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3);
  vector<double> wgt(0);
  // the Z decay modes
  extpart[0]=getParticleData(ParticleID::Z0);
  // loop over the  quarks and the leptons
  int ix,istep=0,iy;
  for( ;istep<11;istep+=10)
    {
      for(ix=1;ix<7;++ix)
	{
	  iy=istep+ix;
	  if(iy!=6)
	    {
	      // check that the combination of particles is allowed
	      if(_zvertex->allowed(-iy,iy,ParticleID::Z0))
		{
		  extpart[1] = getParticleData(-iy);
		  extpart[2] = getParticleData( iy);
		  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
		  if(iy<=6){addMode(mode,_zquarkwgt.at(ix-1),wgt);}
		  else{addMode(mode,_zleptonwgt.at(iy-11),wgt);}
		}
	      else
		{throw InitException() << "SMWZDecayer::doinit() the Z vertex" 
				       << "cannot handle all the modes" 
				       << Exception::abortnow;}
	    }
	}
    }
  // and the W modes
  extpart[0]=getParticleData(ParticleID::Wplus);
  // loop for the quarks
  unsigned int iz=0;
  for(ix=1;ix<6;ix+=2)
    {
      for(iy=2;iy<6;iy+=2)
	{
	  // check that the combination of particles is allowed
	  if(_wvertex->allowed(-ix,iy,ParticleID::Wminus))
	    {
	      extpart[1] = getParticleData(-ix);
	      extpart[2] = getParticleData( iy);
	      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	      addMode(mode,_wquarkwgt[iz],wgt);
	      ++iz;
	    }
	  else
	    {throw InitException() << "SMWZDecayer::doinit() the W vertex" 
				   << "cannot handle all the quark modes" 
				   << Exception::abortnow;}
	}
    }
  for(ix=11;ix<17;ix+=2)
    {
      // check that the combination of particles is allowed
      if(_wvertex->allowed(-ix,ix+1,ParticleID::Wminus))
	{
	  extpart[1] = getParticleData(-ix);
	  extpart[2] = getParticleData(ix+1);
	  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	  addMode(mode,_wleptonwgt[(ix-11)/2],wgt);
	}
	  else
	    {throw InitException() << "SMWZDecayer::doinit() the W vertex" 
				   << "cannot handle all the lepton modes" 
				   << Exception::abortnow;}
    }
}

int SMWZDecayer::modeNumber(bool & cc,tcPDPtr parent, 
			    const PDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  int id0=parent->id();
  PDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  // Z to quarks or leptons
  cc =false;
  if(id0==ParticleID::Z0)
    {
      if(abs(id1)<6&&id1==-id2){imode=abs(id1)-1;}
      else if(abs(id1)>=11&&abs(id1)<=16&&id1==-id2){imode=abs(id1)-6;}
      cc =false;
    }
  // W to quarks and leptons
  else if(abs(id0)==ParticleID::Wplus)
    {
      int idd(0),idu(0);
      if(abs(id1)%2==1&&abs(id2)%2==0)
	{idd=abs(id1);idu=abs(id2);}
      else if(abs(id1)%2==0&&abs(id2)%2==1)
	{idd=abs(id2);idu=abs(id1);}
      if(idd==0&&idu==0)
	{return imode;}
      else if(idd<=5)
	{imode=idd+idu/2+9;}
      else
	{imode=(idd-1)/2+12;}
      cc=(id0==ParticleID::Wminus);
    }
  return imode;
}

void SMWZDecayer::persistentOutput(PersistentOStream & os) const {
  os << _wvertex << _zvertex << _zquarkwgt << _wquarkwgt << _zleptonwgt << _wleptonwgt;
}

void SMWZDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _wvertex >> _zvertex >> _zquarkwgt >> _wquarkwgt >> _zleptonwgt >> _wleptonwgt;
}

ClassDescription<SMWZDecayer> SMWZDecayer::initSMWZDecayer;
// Definition of the static class description member.

void SMWZDecayer::Init() {

  static ClassDocumentation<SMWZDecayer> documentation
    ("The SMWZDecayer class is the implementation of the decay"
     " of the W and Z bosons to the Standard Model fermions.");

  static ParVector<SMWZDecayer,double> interfaceZquarkMax
    ("ZquarkMax",
     "The maximum weight for the decay of the Z to quarks",
     &SMWZDecayer::_zquarkwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZDecayer,double> interfaceWquarkMax
    ("WquarkMax",
     "The maximum weight for the decay of the W to quarks",
     &SMWZDecayer::_wquarkwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZDecayer,double> interfaceZleptonMax
    ("ZleptonMax",
     "The maximum weight for the decay of the Z to leptons",
     &SMWZDecayer::_zleptonwgt,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWZDecayer,double> interfaceWleptonMax
    ("WleptonMax",
     "The maximum weight for the decay of the W to leptons",
     &SMWZDecayer::_wleptonwgt,
     0, 0, 0, -10000, 10000, false, false, true);

}


// return the matrix element squared
double SMWZDecayer::me2(bool vertex, const int, const Particle & inpart,
			const ParticleVector & decay) const 
{
  // get/calculate the spin info for the decaying particle
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<VectorWaveFunction> inwave;
  VectorWaveFunction(inwave,rhoin,const_ptr_cast<tPPtr>(&inpart),incoming,
		     true,false,vertex);
  // construct the spinors for the outgoing particles
  int iferm,ianti;
  if(decay[0]->id()<0){iferm=1;ianti=0;}
  else{iferm=0;ianti=1;}
  vector<SpinorWaveFunction   > awave;
  vector<SpinorBarWaveFunction> fwave;
  SpinorWaveFunction   (awave,decay[ianti],outgoing,true,vertex);
  SpinorBarWaveFunction(fwave,decay[iferm],outgoing,true,vertex);
  // compute the matrix element
  DecayMatrixElement newme(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  Energy2 scale(inpart.mass()*inpart.mass());
  unsigned int ifm,ia,vhel;
  for(ifm=0;ifm<2;++ifm)
    {
      for(ia=0;ia<2;++ia)
	{
	  for(vhel=0;vhel<3;++vhel)
	    {
	      if(inpart.id()==ParticleID::Z0)
		{
		  if(iferm>ianti){newme(vhel,ia,ifm)=
		      _zvertex->evaluate(scale,awave[ia],fwave[ifm],inwave[vhel]);}
		  else{newme(vhel,ifm,ia)=
		      _zvertex->evaluate(scale,awave[ia],fwave[ifm],inwave[vhel]);}
		}
	      else
		{
		  if(iferm>ianti){newme(vhel,ia,ifm)=
		      _wvertex->evaluate(scale,awave[ia],fwave[ifm],inwave[vhel]);}
		  else{newme(vhel,ifm,ia)=
		      _wvertex->evaluate(scale,awave[ia],fwave[ifm],inwave[vhel]);}
		}
	    }
	}
    }
  ME(newme);
  double output=(newme.contract(rhoin)).real()*UnitRemoval::E2/scale;
  if(abs(decay[0]->id())<=6){output*=3.;}
  if(decay[0]->hasColour())
    {decay[0]->antiColourNeighbour(decay[1]);}
  else if(decay[1]->hasColour())
    {decay[1]->antiColourNeighbour(decay[0]);}
  return output;
}
