// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiPiDecayer class.
//

#include "EtaPiPiPiDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EtaPiPiPiDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "Herwig++/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig++/PDT/OneOffShellCalculator.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;

EtaPiPiPiDecayer::~EtaPiPiPiDecayer() {}

bool EtaPiPiPiDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  // check three outgoing particles
  if(dm.products().size()!=3){return false;}
  unsigned int npi0(0),npip(0),npim(0); int id,iother(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {
      id=(**pit).id();
      if(id==ParticleID::piplus){++npip;}
      else if(id==ParticleID::piminus){++npim;}
      else if(id==ParticleID::pi0&&npi0<2){++npi0;}
      else{iother=id;}
    }
  if(!(npi0==2||(npip==1&&npim==1))){return allowed;}
  if(npi0==1){iother=ParticleID::pi0;}
  id=dm.parent()->id();
  bool charged=npi0<2;
  unsigned int ix=0;
  do 
    {
      allowed=(id==_incoming[ix]&&iother==_outgoing[ix]&&charged==_charged[ix]);
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector EtaPiPiPiDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int idout(0),id;
  unsigned int npi0(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {
      id=(**pit).id();
      if(id==ParticleID::pi0&&npi0<2){++npi0;}
      else if(id!=ParticleID::piplus&&id!=ParticleID::piminus){idout=id;}
    }
  if(npi0==1){idout=ParticleID::pi0;}
  bool charged=npi0<2;
  int imode=-1;
  id=parent.id();
  unsigned int ix=0;
  do 
    {
      if(id==_incoming[ix]&&idout==_outgoing[ix]&&_charged[ix]==charged){imode=ix;}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  if(imode<0){throw DecayIntegratorError() << "Unknown mode in EtaPiPiPiDecayer::decay()"
					   << Exception::runerror;}
  bool cc=false;
  return generate(false,cc,imode,parent);

}

void EtaPiPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing << _charged << _prefactor << _a << _b << _c  
     << _maxweight;
}

void EtaPiPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing >> _charged >> _prefactor >> _a >> _b >> _c 
     >> _maxweight;
}

ClassDescription<EtaPiPiPiDecayer> EtaPiPiPiDecayer::initEtaPiPiPiDecayer;
// Definition of the static class description member.

void EtaPiPiPiDecayer::Init() {

  static ClassDocumentation<EtaPiPiPiDecayer> documentation
    ("The \\classname{EtaPiPiPiDecayer} class");

}

double EtaPiPiPiDecayer::me2(bool vertex,const int,const Particle & inpart,
			     const ParticleVector & decay) const
{
  // construct spin info objects (this is pretty much a waste of time)
  // check if the decay particle has spin info 
  tcScalarSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcScalarSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  if(inspin)
    {inspin->decayed(true);}
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for th incoming particle"
	 			  << " in EtaPiGammaGammaDecayer::me2()" 
	 			  << Exception::abortnow;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  if(vertex)
    {for(unsigned int ix=0;ix<decay.size();++ix)
	{decay[ix]->spinInfo(new_ptr(ScalarSpinInfo(decay[ix]->momentum(),true)));}}
  // calculate the matrix element
  // compute the variables we need
  Lorentz5Momentum ps=inpart.momentum()-decay[2]->momentum();ps.rescaleMass();
  Lorentz5Momentum pu=inpart.momentum()-decay[0]->momentum();pu.rescaleMass();
  Lorentz5Momentum pt=inpart.momentum()-decay[1]->momentum();pt.rescaleMass();
  Energy2 s =ps.mass2();
  Energy2 u =pu.mass2();
  Energy2 t =pt.mass2();
  Energy m34=decay[0]->mass()+decay[1]->mass();
  Energy msum=decay[2]->mass()+m34;
  Energy Q=inpart.mass()-msum;
  Energy2 Mmm2=(inpart.mass()-decay[2]->mass())*(inpart.mass()-decay[2]->mass());
  // compute the variables
  double x = 0.5*sqrt(3.)*(u-t)/inpart.mass()/Q;
  double y = 0.5*msum/inpart.mass()*(Mmm2-s)/m34/Q-1;
  double x2=x*x,y2=y*y;
  double me;
  me=_prefactor[imode()]*(1+_a[imode()]*y+_b[imode()]*y2+_c[imode()]*x2);
  DecayMatrixElement newME(1,1,1,1);
  newME(0,0,0,0)=sqrt(me);
  ME(newME);
  return me;
}

double EtaPiPiPiDecayer::threeBodydGammads(int imodeb,Energy q2, Energy2 s,
					   Energy m1,Energy m2,Energy m3)
{
  Energy q=sqrt(q2);
  Energy m34=m1+m2;
  Energy msum=m34+m3;
  Energy Q=q-msum;
  Energy2 Mmm2=(q-m3)*(q-m3);
  Energy2 m12=m1*m1,m22=m2*m2,m32=m3*m3;
  double y = 0.5*msum/q*(Mmm2-s)/m34/Q-1;
  double y2=y*y;
  InvEnergy2 xfact=0.5*sqrt(3.)/q/Q;
  Energy2 xc=q2+m12+m22+m32-s;
  Energy rs=sqrt(s);
  Energy e2star = 0.5*(s-m12+m22)/rs;
  Energy e3star = 0.5*(q2-s-m32)/rs;
  Energy e2sm=sqrt(e2star*e2star-m22);
  Energy e3sm=sqrt(e3star*e3star-m32);
  Energy2 a = 2*e2star*e3star+m22+m32;
  Energy2 b = 2*e2sm*e3sm;
  double output=2*b*(1+_a[imodeb]*y+_b[imodeb]*y2+_c[imodeb]*xfact*xfact*(xc*xc))
    +_c[imodeb]*(-8.*xfact*xfact*xc*a*b
		 +4.*2*b*(3.*a*a+b*b)/3.*xfact*xfact);
  return output*_prefactor[imodeb]/256./pi/pi/pi/q2/q;
}


WidthCalculatorBasePtr 
EtaPiPiPiDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  int idout(0),id;
  unsigned int npi0(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {
      id=(**pit).id();
      if(id==ParticleID::pi0&&npi0<2){++npi0;}
      else if(id!=ParticleID::piplus&&id!=ParticleID::piminus){idout=id;}
    }
  if(npi0==1){idout=ParticleID::pi0;}
  bool charged=npi0<2;
  int imode=-1;
  id=dm.parent()->id();
  unsigned int ix=0;
  do 
    {
      if(id==_incoming[ix]&&idout==_outgoing[ix]&&_charged[ix]==charged){imode=ix;}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  Energy mpi;
  if(charged){mpi=getParticleData(ParticleID::piplus)->mass();}
  else{mpi=getParticleData(ParticleID::pi0)->mass();}
  Energy m[3]={mpi,mpi,getParticleData(_outgoing[imode])->mass()};
  tDecayIntegratorPtr decayer=const_ptr_cast<tDecayIntegratorPtr>(this);
  WidthCalculatorBasePtr 
    temp=new_ptr(ThreeBodyAllOn1IntegralCalculator(1,-1000.,0.0,decayer,imode,
						   m[0],m[1],m[2]));
  if(_outgoing[imode]==ParticleID::eta)
    {
      tcGenericMassGeneratorPtr test;
      tGenericMassGeneratorPtr massptr;
      if(getParticleData(_outgoing[imode])->massGenerator())
	{
	  test=dynamic_ptr_cast<tcGenericMassGeneratorPtr>
	    (getParticleData(_outgoing[imode])->massGenerator());
	  massptr=const_ptr_cast<tGenericMassGeneratorPtr>(test);
	}
      if(massptr)
	{
	  massptr->init();
	  return new_ptr(OneOffShellCalculator(3,temp,massptr,0.));
	}
    }
  return temp;
} 
}
