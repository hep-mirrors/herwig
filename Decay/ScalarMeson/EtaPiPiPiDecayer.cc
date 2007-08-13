// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiPiDecayer class.
//
#include "EtaPiPiPiDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig++/PDT/OneOffShellCalculator.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

EtaPiPiPiDecayer::EtaPiPiPiDecayer() 
  : _incoming(6), _outgoing(6), _charged(6), _prefactor(6),
    _a(6), _b(6), _c(6), _maxweight(6) {
  // eta to pi+pi-pi0
  _incoming[0] = 221; _outgoing[0] = 111; _charged[0] = true; 
  _prefactor[0] = 0.0404509; _maxweight[0] = 1.32;
  _a[0] = -1.17; _b[0] = 0.21; _c[0] = 0.06; 
  // eta to pi0pi0pi0
  _incoming[1] = 221; _outgoing[1] = 111; _charged[1] = false; 
  _prefactor[1] = 0.0883547; _maxweight[1] = 1.33; 
  _a[1] = 0.; _b[1] = -0.062; _c[1] = -0.062; 
  // eta' to pi+pi-pi0
  _incoming[2] = 331; _outgoing[2] = 111; _charged[2] = true; 
  _prefactor[2] = 0.037165; _maxweight[2] = 0.0227363;
  _a[2] = -3.08; _b[2] = 0.13; _c[2] = 0.62; 
  // eta' to pi0pi0pi0
  _incoming[3] = 331; _outgoing[3] = 111; _charged[3] = false; 
  _prefactor[3] = 0.016203; _maxweight[3] = 2.26; 
  _a[3] = 0.0; _b[3] = -0.86; _c[3] = -0.86; 
  // eta' to pi+pi-eta
  _incoming[4] = 331; _outgoing[4] = 221; _charged[4] = true; 
  _prefactor[4] = 46.47; _maxweight[4] = 1.30;
  _a[4] = -0.093; _b[4] = -0.059; _c[4] = -0.003; 
  // eta' to pi0pi0eta
  _incoming[5] = 331; _outgoing[5] = 221; _charged[5] = false; 
  _prefactor[5] = 19.408225; _maxweight[5] = 1.30;
  _a[5] = -0.105; _b[5] = -0.065; _c[5] = -0.004; 
  // initial size of the arrays
  _initsize=_maxweight.size();
  // intermediates
  generateIntermediates(false);
}
 
void EtaPiPiPiDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize(_incoming.size());
  if(isize!=_outgoing.size()||isize!=_prefactor.size()||
     isize!=_charged.size()||isize!=_a.size()||
     isize!=_b.size()||isize!=_c.size()||isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in EtaPiPiPiDecayer::doinit()"
			  << Exception::runerror;
  // external particles for the modes
  PDVector extneut(4),extcharged(4);
  extneut[1]    = getParticleData(ParticleID::pi0);
  extneut[2]    = getParticleData(ParticleID::pi0);
  extcharged[1] = getParticleData(ParticleID::piplus);
  extcharged[2] = getParticleData(ParticleID::piminus);
  tPDPtr rho(getParticleData(113));
  DecayPhaseSpaceModePtr mode;
  DecayPhaseSpaceChannelPtr newchannel;
  vector<double> dummyweights(1,1.);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extneut[0]    = getParticleData(_incoming[ix]);
    extcharged[0] = getParticleData(_incoming[ix]);
    extneut[3]    = getParticleData(_outgoing[ix]);
    extcharged[3] = getParticleData(_outgoing[ix]);
    if(_charged[ix]) {
      // the pi+pi- mode
      mode = new_ptr(DecayPhaseSpaceMode(extcharged,this));
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extcharged[0],0, 0.0,-1,3);
      newchannel->addIntermediate(rho,1,0.0, 1,2);
      mode->addChannel(newchannel);
    }
    else {
      // the pi0pi0 mode
      mode = new_ptr(DecayPhaseSpaceMode(extneut,this));
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extneut[0],0, 0.0,-1,3);
      newchannel->addIntermediate(rho,1,0.0, 1,2);
      mode->addChannel(newchannel);
    }
    addMode(mode,_maxweight[ix],dummyweights);
  }
  resetIntermediate(rho,600.*MeV,600.*MeV);
}

int EtaPiPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
				 const PDVector & children) const {
  if(children.size()!=3) return -1;
  unsigned int npi0(0),npip(0),npim(0); int id,iother(0);
  PDVector::const_iterator pit = children.begin();
  for( ;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)           ++npip;
    else if(id==ParticleID::piminus)     ++npim;
    else if(id==ParticleID::pi0&&npi0<2) ++npi0;
    else iother=id;
  }
  bool charged;
  if(npim==1&&npip==1) {
    charged=true;
    if(npi0==1) iother=ParticleID::pi0;
  }
  else if(npi0==2) charged=false;
  else return -1;
  // find the mode
  id=parent->id();
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id==_incoming[ix]&&iother==_outgoing[ix]&&_charged[ix]==charged) 
      imode=ix;
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
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
    ("The EtaPiPiPiDecayer class performs the decay of a scalar meson to"
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
			     const ParticleVector & decay) const {
  // workaround for gcc 3.2.3 bug
  // construct spin info objects (this is pretty much a waste of time)
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    //ALB ScalarWaveFunction(decay[ix],outgoing,true,vertex);
    PPtr mytemp = decay[ix]; ScalarWaveFunction(mytemp,outgoing,true,vertex);
  }
  // calculate the matrix element
  // compute the variables we need
  Lorentz5Momentum ps(inpart.momentum()-decay[2]->momentum());ps.rescaleMass();
  Lorentz5Momentum pu(inpart.momentum()-decay[0]->momentum());pu.rescaleMass();
  Lorentz5Momentum pt(inpart.momentum()-decay[1]->momentum());pt.rescaleMass();
  Energy2 s(ps.mass2()),u(pu.mass2()),t(pt.mass2());
  Energy m34(0.5*(decay[0]->mass()+decay[1]->mass()));
  Energy msum(decay[2]->mass()+2.*m34);
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

InvEnergy EtaPiPiPiDecayer::threeBodydGammads(const int imodeb, const Energy2 q2,
					   const  Energy2 s, const Energy m1,
					   const Energy m2, const Energy m3) const {
  Energy q(sqrt(q2)),m34(m1+m2),msum(m34+m3),Q(q-msum);
  Energy2 Mmm2((q-m3)*(q-m3)),m12(m1*m1),m22(m2*m2),m32(m3*m3);
  double y(0.5*msum/q*(Mmm2-s)/m34/Q-1),y2(y*y);
  InvEnergy2 xfact=0.5*sqrt(3.)/q/Q;
  Energy2 xc(q2+m12+m22+m32-s);
  Energy rs(sqrt(s)),e2star(0.5*(s-m12+m22)/rs),e3star(0.5*(q2-s-m32)/rs);
  Energy e2sm(sqrt(e2star*e2star-m22)),e3sm(sqrt(e3star*e3star-m32));
  Energy2 a(2*e2star*e3star+m22+m32),b(2*e2sm*e3sm);
  Energy2 output=2*b*(1+_a[imodeb]*y+_b[imodeb]*y2+_c[imodeb]*xfact*xfact*(xc*xc))
    +_c[imodeb]*(-8.*xfact*xfact*xc*a*b
		 +4.*2*b*(3.*a*a+b*b)/3.*xfact*xfact);
  using Constants::pi;
  return output*_prefactor[imodeb]/256./pi/pi/pi/q2/q;
}


WidthCalculatorBasePtr 
EtaPiPiPiDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  int idout(0),id,imode(-1);
  unsigned int npi0(0),ix(0);
  ParticleMSet::const_iterator pit(dm.products().begin());
  for( ;pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::pi0&&npi0<2)                            ++npi0;
    else if(id!=ParticleID::piplus&&id!=ParticleID::piminus) idout=id;
  }
  if(npi0==1) idout=ParticleID::pi0;
  bool charged(npi0<2);
  id=dm.parent()->id();
  do {
    if(id==_incoming[ix]&&idout==_outgoing[ix]&&_charged[ix]==charged) 
      imode=ix;
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  Energy mpi;
  if(charged){mpi=getParticleData(ParticleID::piplus)->mass();}
  else{mpi=getParticleData(ParticleID::pi0)->mass();}
  Energy m[3]={mpi,mpi,getParticleData(_outgoing[imode])->mass()};
  WidthCalculatorBasePtr 
    temp(new_ptr(ThreeBodyAllOn1IntegralCalculator<EtaPiPiPiDecayer>
		 (1,-1000.*MeV,0.0*MeV,0.0,*this,imode,m[0],m[1],m[2])));
  if(_outgoing[imode]==ParticleID::eta) {
    tcGenericMassGeneratorPtr test;
    tGenericMassGeneratorPtr massptr;
    if(getParticleData(_outgoing[imode])->massGenerator()) {
      test=dynamic_ptr_cast<tcGenericMassGeneratorPtr>
	(getParticleData(_outgoing[imode])->massGenerator());
      massptr=const_ptr_cast<tGenericMassGeneratorPtr>(test);
    }
    if(massptr) {
      massptr->init();
      return new_ptr(OneOffShellCalculator(3,temp,massptr,0.*MeV));
    }
  }
  return temp;
} 
  
void EtaPiPiPiDecayer::dataBaseOutput(ofstream & output,
				      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
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
    else {
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
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
