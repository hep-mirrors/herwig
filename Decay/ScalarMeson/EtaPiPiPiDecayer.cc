// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiPiDecayer class.
//

#include "EtaPiPiPiDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EtaPiPiPiDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig++/PDT/OneOffShellCalculator.h"

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::ScalarWaveFunction;

EtaPiPiPiDecayer::EtaPiPiPiDecayer() 
{
  // eta to pi+pi-pi0
  _incoming.push_back(221);_outgoing.push_back(111);_charged.push_back(true);
  _prefactor.push_back(0.0404509);
  _a.push_back(-1.17);_b.push_back(0.21);_c.push_back(0.06);
  _maxweight.push_back(1.32);
  // eta to pi0pi0pi0
  _incoming.push_back(221);_outgoing.push_back(111);_charged.push_back(false);
  _prefactor.push_back(0.0883547);
  _a.push_back(0.);_b.push_back(-0.062);_c.push_back(-0.062);
  _maxweight.push_back(1.33);
  // eta' to pi+pi-pi0
  _incoming.push_back(331);_outgoing.push_back(111);_charged.push_back(true);
  _prefactor.push_back(0.037165);
  _a.push_back(-3.08);_b.push_back(0.13);_c.push_back(0.62);
  _maxweight.push_back(0.0227363);
  // eta' to pi0pi0pi0
  _incoming.push_back(331);_outgoing.push_back(111);_charged.push_back(false);
  _prefactor.push_back(0.016203);
  _a.push_back(0.0);_b.push_back(-0.86);_c.push_back(-0.86);
  _maxweight.push_back(2.26);
  // eta' to pi+pi-eta
  _incoming.push_back(331);_outgoing.push_back(221);_charged.push_back(true);
  _prefactor.push_back(46.47);
  _a.push_back(-0.093);_b.push_back(-0.059);_c.push_back(-0.003);
  _maxweight.push_back(1.30);
  // eta' to pi0pi0eta
  _incoming.push_back(331);_outgoing.push_back(221);_charged.push_back(false);
  _prefactor.push_back(19.408225);
  _a.push_back(-0.105);_b.push_back(-0.065);_c.push_back(-0.004);
  _maxweight.push_back(1.30);
  // initial size of the arrays
  _initsize=_maxweight.size();
}

void EtaPiPiPiDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize(_incoming.size());
  if(isize!=_outgoing.size()||isize!=_prefactor.size()||
     isize!=_charged.size()||isize!=_a.size()||
     isize!=_b.size()||isize!=_c.size()||isize!=_maxweight.size())
    {throw InitException() << "Inconsistent parameters in EtaPiPiPiDecayer::doinit()"
			   << Exception::runerror;}
  // external particles for the modes
  PDVector extneut(4),extcharged(4);
  extneut[1]    = getParticleData(ParticleID::pi0);
  extneut[2]    = getParticleData(ParticleID::pi0);
  extcharged[1] = getParticleData(ParticleID::piplus);
  extcharged[2] = getParticleData(ParticleID::piminus);
  tPDPtr sigma(getParticleData(9000221));
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr newchannel;
  vector<double> dummyweights(1,1.);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extneut[0]    = getParticleData(_incoming[ix]);
      extcharged[0] = getParticleData(_incoming[ix]);
      extneut[3]    = getParticleData(_outgoing[ix]);
      extcharged[3] = getParticleData(_outgoing[ix]);
      if(_charged[ix])
	{
	  // the pi+pi- mode
	  mode = new_ptr(DecayPhaseSpaceMode(extcharged,this));
	  newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	  newchannel->addIntermediate(extcharged[0],0, 0.0,-1,3);
	  newchannel->addIntermediate(sigma,1,0.0, 1,2);
	  mode->addChannel(newchannel);
	}
      else
	{
	  // the pi0pi0 mode
	  mode = new_ptr(DecayPhaseSpaceMode(extneut,this));
	  newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
	  newchannel->addIntermediate(extneut[0],0, 0.0,-1,3);
	  newchannel->addIntermediate(sigma,1,0.0, 1,2);
	  mode->addChannel(newchannel);
	}
      addMode(mode,_maxweight[ix],dummyweights);
    }
}

EtaPiPiPiDecayer::~EtaPiPiPiDecayer() {}

bool EtaPiPiPiDecayer::accept(const DecayMode & dm) const {
  bool allowed(false);
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
  bool charged(npi0<2);
  unsigned int ix(0);
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
  int idout(0),id,imode(-1);
  unsigned int npi0(0),ix(0);
  ParticleMSet::const_iterator pit(dm.products().begin());
  for( ;pit!=dm.products().end();++pit)
    {
      id=(**pit).id();
      if(id==ParticleID::pi0&&npi0<2){++npi0;}
      else if(id!=ParticleID::piplus&&id!=ParticleID::piminus){idout=id;}
    }
  if(npi0==1){idout=ParticleID::pi0;}
  bool charged(npi0<2);
  id=parent.id();
  do 
    {
      if(id==_incoming[ix]&&idout==_outgoing[ix]&&_charged[ix]==charged){imode=ix;}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  if(imode<0){throw DecayIntegratorError() << "Unknown mode in EtaPiPiPiDecayer::decay()"
					   << Exception::runerror;}
  bool cc(false);
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
    ("The \\classname{EtaPiPiPiDecayer} class performs the decay of a scalar meson to"
     " two pions and another meson using a simple paramterisation of the dalitz plot.");

  static ParVector<EtaPiPiPiDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming particle",
     &EtaPiPiPiDecayer::_incoming, -1, 0,  0, 1000000,
     false, false, true);

  static ParVector<EtaPiPiPiDecayer,int> interfaceOutgoing
    ("Outgoing",
     "The PDG code of the outgoing particle",
     &EtaPiPiPiDecayer::_outgoing, -1, 0,  0, 1000000,
     false, false, true);

  static ParVector<EtaPiPiPiDecayer,bool> interfaceCharged
    ("Charged",
     "Whether the pions or charged or neutral",
     &EtaPiPiPiDecayer::_charged,  -1,false, 0, 0,
     false, false, false);

  static ParVector<EtaPiPiPiDecayer,double> interfacePrefactor
    ("Prefactor",
     "The prefactor for the decay to get the correct partial width",
     &EtaPiPiPiDecayer::_prefactor, -1,1.0,  0, 0,
     false, false, false);

  static ParVector<EtaPiPiPiDecayer,double> interfacea
    ("a",
     "The a parameter for the dalitz plot",
     &EtaPiPiPiDecayer::_a, -1, 0.0,  -10.0, 10.0,
     false, false, true);

  static ParVector<EtaPiPiPiDecayer,double> interfaceb
    ("b",
     "The b parameter for the dalitz plot",
     &EtaPiPiPiDecayer::_b, -1, 0.0,  -10.0, 10.0,
     false, false, true);

  static ParVector<EtaPiPiPiDecayer,double> interfacec
    ("c",
     "The c parameter for the dalitz plot",
     &EtaPiPiPiDecayer::_c, -1, 0.0,  -10.0, 10.0,
     false, false, true);

  static ParVector<EtaPiPiPiDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &EtaPiPiPiDecayer::_maxweight,
     0, 0, 0, 0., 200., false, false, true);

}

double EtaPiPiPiDecayer::me2(bool vertex,const int,const Particle & inpart,
			     const ParticleVector & decay) const
{
  // workaround for gcc 3.2.3 bug
  // construct spin info objects (this is pretty much a waste of time)
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  for(unsigned int ix=0;ix<decay.size();++ix)
    //ALB {ScalarWaveFunction(decay[ix],outgoing,true,vertex);}
    {PPtr mytemp = decay[ix]; ScalarWaveFunction(mytemp,outgoing,true,vertex);}

  // calculate the matrix element
  // compute the variables we need
  Lorentz5Momentum ps(inpart.momentum()-decay[2]->momentum());ps.rescaleMass();
  Lorentz5Momentum pu(inpart.momentum()-decay[0]->momentum());pu.rescaleMass();
  Lorentz5Momentum pt(inpart.momentum()-decay[1]->momentum());pt.rescaleMass();
  Energy2 s(ps.mass2()),u(pu.mass2()),t(pt.mass2());
  Energy m34(decay[0]->mass()+decay[1]->mass());
  Energy msum(decay[2]->mass()+m34);
  Energy Q(inpart.mass()-msum);
  Energy2 Mmm2((inpart.mass()-decay[2]->mass())*(inpart.mass()-decay[2]->mass()));
  // compute the variables
  double x(0.5*sqrt(3.)*(u-t)/inpart.mass()/Q),x2(x*x);
  double y(0.5*msum/inpart.mass()*(Mmm2-s)/m34/Q-1),y2(y*y);
  double me(_prefactor[imode()]*(1+_a[imode()]*y+_b[imode()]*y2+_c[imode()]*x2));
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newME(0,0,0,0)=sqrt(me);
  ME(newME);
  return me;
}

double EtaPiPiPiDecayer::threeBodydGammads(int imodeb,Energy q2, Energy2 s,
					   Energy m1,Energy m2,Energy m3)
{
  Energy q(sqrt(q2)),m34(m1+m2),msum(m34+m3),Q(q-msum);
  Energy2 Mmm2((q-m3)*(q-m3)),m12(m1*m1),m22(m2*m2),m32(m3*m3);
  double y(0.5*msum/q*(Mmm2-s)/m34/Q-1),y2(y*y);
  InvEnergy2 xfact=0.5*sqrt(3.)/q/Q;
  Energy2 xc(q2+m12+m22+m32-s);
  Energy rs(sqrt(s)),e2star(0.5*(s-m12+m22)/rs),e3star(0.5*(q2-s-m32)/rs);
  Energy e2sm(sqrt(e2star*e2star-m22)),e3sm(sqrt(e3star*e3star-m32));
  Energy2 a(2*e2star*e3star+m22+m32),b(2*e2sm*e3sm);
  double output=2*b*(1+_a[imodeb]*y+_b[imodeb]*y2+_c[imodeb]*xfact*xfact*(xc*xc))
    +_c[imodeb]*(-8.*xfact*xfact*xc*a*b
		 +4.*2*b*(3.*a*a+b*b)/3.*xfact*xfact);
  return output*_prefactor[imodeb]/256./pi/pi/pi/q2/q;
}


WidthCalculatorBasePtr 
EtaPiPiPiDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  int idout(0),id,imode(-1);
  unsigned int npi0(0),ix(0);
  ParticleMSet::const_iterator pit(dm.products().begin());
  for( ;pit!=dm.products().end();++pit)
    {
      id=(**pit).id();
      if(id==ParticleID::pi0&&npi0<2){++npi0;}
      else if(id!=ParticleID::piplus&&id!=ParticleID::piminus){idout=id;}
    }
  if(npi0==1){idout=ParticleID::pi0;}
  bool charged(npi0<2);
  id=dm.parent()->id();
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
    temp(new_ptr(ThreeBodyAllOn1IntegralCalculator(1,-1000.,0.0,decayer,imode,
						   m[0],m[1],m[2])));
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
  
void EtaPiPiPiDecayer::dataBaseOutput(ofstream & output,
				      bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming   " << ix << " "
		 << _incoming[ix]   << "\n";
	  output << "set " << fullName() << ":Outgoing  " << ix << " "
		 << _outgoing[ix]  << "\n";
	  output << "set " << fullName() << ":Charged " << ix << " "
		 << _charged[ix]  << "\n";
	  output << "set " << fullName() << ":Prefactor " << ix << " "
		 << _prefactor[ix]  << "\n";
	  output << "set " << fullName() << ":a " << ix << " "
		 << _a[ix]  << "\n";
	  output << "set " << fullName() << ":b " << ix << " "
		 << _b[ix]  << "\n";
	  output << "set " << fullName() << ":c " << ix << " "
		 << _c[ix]  << "\n";
	  output << "set " << fullName() << ":MaxWeight  " << ix << " "
		 << _maxweight[ix]  << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming   " << ix << " "
		 << _incoming[ix]   << "\n";
	  output << "insert " << fullName() << ":Outgoing  " << ix << " "
		 << _outgoing[ix]  << "\n";
	  output << "insert " << fullName() << ":Charged " << ix << " "
		 << _charged[ix]  << "\n";
	  output << "insert " << fullName() << ":Prefactor " << ix << " "
		 << _prefactor[ix]  << "\n";
	  output << "insert " << fullName() << ":a " << ix << " "
		 << _a[ix]  << "\n";
	  output << "insert " << fullName() << ":b " << ix << " "
		 << _b[ix]  << "\n";
	  output << "insert " << fullName() << ":c " << ix << " "
		 << _c[ix]  << "\n";
	  output << "insert " << fullName() << ":MaxWeight  " << ix << " "
		 << _maxweight[ix]  << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
