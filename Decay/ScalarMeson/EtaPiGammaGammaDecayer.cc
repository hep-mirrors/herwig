// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiGammaGammaDecayer class.
//

#include "EtaPiGammaGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EtaPiGammaGammaDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::outgoing;
using Helicity::incoming;
using Helicity::VectorWaveFunction;
using Helicity::ScalarWaveFunction;

void EtaPiGammaGammaDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // set rho parameters if needed
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!_localparameters)
    {
      _rhomass  = rho->mass();
      _rhowidth = rho->width();
    }
  // constant for the running rho width
  _mpi=getParticleData(ParticleID::pi0)->mass();
  Energy pcm =Kinematics::pstarTwoBodyDecay(_rhomass,_mpi,_mpi);
  _rhoconst=_rhomass*_rhomass*_rhowidth/(pcm*pcm*pcm);
  // set the prefactors
  double conv(sqrt(4.*pi*SM().alphaEM()));
  conv *=_fpi*_fpi*_grho/_rhomass/_rhomass;
  InvEnergy pre(2.*sqrt(3.)/9.*_grhoomega*_grhoomega*conv*conv);
  double fact[2];
  // constants for eta
  fact[0] = _ratiofpif8*cos(_theta)-sqrt(2.)*_ratiofpif0*sin(_theta);
  // constants for eta'
  fact[1] = _ratiofpif8*sin(_theta)+sqrt(2.)*_ratiofpif0*cos(_theta);
  for(unsigned int ix=0;ix<2;++ix)
    {
      _Dconst[ix]=fact[ix]*pre;
      _Econst[ix]=fact[ix]*pre;
    }
  // set up the phsae space for the decays
  tPDPtr eta[2]={getParticleData(ParticleID::eta),getParticleData(ParticleID::etaprime)};
  PDVector extpart;extpart.resize(4);
  extpart[1] = getParticleData(ParticleID::pi0);
  extpart[2] = getParticleData(ParticleID::gamma);
  extpart[3] = getParticleData(ParticleID::gamma);
  vector<double> dummyweights(2,0.5);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<2;++ix)
    {
      extpart[0] = eta[ix];
      mode = new DecayPhaseSpaceMode(extpart,this);
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
      newchannel->addIntermediate(rho,0,0.0, 1,3);
      mode->addChannel(newchannel);
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
      newchannel->addIntermediate(rho,0,0.0, 1,2);
      mode->addChannel(newchannel);
      if(ix==0){addMode(mode,_etamax,dummyweights);}
      else if(ix==1){addMode(mode,_etapmax,dummyweights);}
    }
}

EtaPiGammaGammaDecayer::~EtaPiGammaGammaDecayer() {}

bool EtaPiGammaGammaDecayer::accept(const DecayMode & dm) const {
  bool allowed=false;
  // check three outgoing particles
  if(dm.products().size()!=3){return false;}
  // check the incoming particle is an eta or eta prime
  int id=dm.parent()->id();
  if(id==ParticleID::eta||id==ParticleID::etaprime)
    {
      ParticleMSet::const_iterator pit = dm.products().begin();
      unsigned int npi0=0,ngamma=0;
      for( ;pit!=dm.products().end();++pit)
	{
	  id=(**pit).id();
	  if(id==ParticleID::pi0){++npi0;}
	  else if(id==ParticleID::gamma){++ngamma;}
	}
      if(npi0==1&&ngamma==2){allowed=true;}
    }
  return allowed;
}

ParticleVector EtaPiGammaGammaDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int id(parent.id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  bool cc(false);
  return generate(false,cc,imode,parent);
}


void EtaPiGammaGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _grhoomega << _fpi << _grho << _rhomass << _rhowidth << _localparameters 
     << _ratiofpif8 << _ratiofpif0 << _theta << _etamax << _etapmax << _rhoconst << _mpi;
  for(unsigned int ix=0;ix<2;++ix){os << _Dconst[ix] << _Econst[ix];}
}

void EtaPiGammaGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _grhoomega >> _fpi >> _grho >> _rhomass >> _rhowidth >> _localparameters 
     >> _ratiofpif8 >> _ratiofpif0 >> _theta >> _etamax >> _etapmax >> _rhoconst >> _mpi;
  for(unsigned int ix=0;ix<2;++ix){is >> _Dconst[ix] >> _Econst[ix];}
}

ClassDescription<EtaPiGammaGammaDecayer> EtaPiGammaGammaDecayer::initEtaPiGammaGammaDecayer;
// Definition of the static class description member.

void EtaPiGammaGammaDecayer::Init() {

  static ClassDocumentation<EtaPiGammaGammaDecayer> documentation
    ("The \\classname{EtaPiGammaGammaDecayer} class implements a VMD model for the"
     " decay of the eta or etaprime to a pion and two photons.");

  static Parameter<EtaPiGammaGammaDecayer,InvEnergy> interfacegrhoomega
    ("grhoomega",
     "The couping of the rho, omega and a pion",
     &EtaPiGammaGammaDecayer::_grhoomega, 1./GeV, 12.924/GeV, 0.0/GeV, 100./GeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant",
     &EtaPiGammaGammaDecayer::_fpi, MeV, 130.7*MeV, 0.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfacegrho
    ("grho",
     "Rho decay constant",
     &EtaPiGammaGammaDecayer::_grho, 5.9, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho meson",
     &EtaPiGammaGammaDecayer::_rhomass, MeV, 771.1*MeV, 500.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho meson",
     &EtaPiGammaGammaDecayer::_rhowidth, MeV, 149.2*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceRatioFpiF8
    ("RatioFpiF8",
     "The ratio of the decay constant Fpi to F8",
     &EtaPiGammaGammaDecayer::_ratiofpif8, 1./1.3, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceRatioFpiF0
    ("RatioFpiF0",
     "The ratio of the decay constant Fpi to F0",
     &EtaPiGammaGammaDecayer::_ratiofpif0, 1./1.04, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceTheta
    ("Theta",
     "The eta etaprime mixing angle",
     &EtaPiGammaGammaDecayer::_theta, -pi/9., -pi, pi,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaMax
    ("EtaMax",
     "THe maximum weight for the eta decay",
     &EtaPiGammaGammaDecayer::_etamax, 1.35, -1.0e12, 1.0e12,
     false, false, false);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaPrimeMax
    ("EtaPrimeMax",
     "THe maximum weight for the eta prime decay",
     &EtaPiGammaGammaDecayer::_etapmax, 0.006, -1.0e12, 1.0e12,
     false, false, false);

  static Switch<EtaPiGammaGammaDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the parameters",
     &EtaPiGammaGammaDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);
}

double EtaPiGammaGammaDecayer::me2(bool vertex, const int,const Particle & inpart,
				   const ParticleVector & decay) const
{
  unsigned int ix,iy;
  // workaround for gcc 3.2.3 bug
  // spin info of the decaying particle
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);

  // spin info and wavefunctions for outgoing particles
  vector<LorentzPolarizationVector> vwave[2];

  // workaround for gcc 3.2.3 bug
  for(ix=0;ix<2;++ix)
    //ALB  {VectorWaveFunction(vwave[ix],decay[ix+1],outgoing,true,true,vertex);}
  //ALB ScalarWaveFunction(decay[0],outgoing,true,vertex);
    {
      vector<LorentzPolarizationVector> mytempLPV;
      VectorWaveFunction(mytempLPV,decay[ix+1],outgoing,true,true,vertex);
      vwave[ix]=mytempLPV; 
    }
  PPtr mytemp = decay[0];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);

  // dot products we need
  Energy2 q1dotq2(decay[1]->momentum()*decay[2]->momentum()),
    pdotq1(inpart.momentum()*decay[1]->momentum()),
    pdotq2(inpart.momentum()*decay[2]->momentum());
  complex<Energy> e1dotq2[3],e1dotp[3],e2dotq1[3],e2dotp[3];
  for(ix=0;ix<3;++ix)
    {
      if(ix==1)
	{e1dotq2[ix]=0.;e1dotp[ix] =0.;e2dotq1[ix]=0.;e2dotp[ix] =0.;}
      else
	{
	  e1dotq2[ix] =vwave[0][ix]*decay[2]->momentum();
	  e1dotp[ix]  =vwave[0][ix]*inpart.momentum();
	  e2dotq1[ix] =vwave[1][ix]*decay[1]->momentum();
	  e2dotp[ix]  =vwave[1][ix]*inpart.momentum();
	}
    }
  // the momentum dependent pieces of the matrix element
  Complex ii(0.,1.);
  Energy2 mpi2(decay[0]->mass()*decay[0]->mass()),meta2(inpart.mass()*inpart.mass()),
    mrho2(_rhomass*_rhomass),
    t(mpi2+2.*((decay[0]->momentum())*(decay[1]->momentum()))),
    u(mpi2+2.*((decay[0]->momentum())*(decay[2]->momentum())));
  Energy q(sqrt(t)),pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
  complex<Energy> tgamma(ii*pcm*pcm*pcm*_rhoconst/q);
  q=sqrt(u);pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy> ugamma(ii*pcm*pcm*pcm*_rhoconst/q);
  complex<InvEnergy2> prop1(1./(mrho2-t-tgamma)),prop2(1./(mrho2-u-ugamma));
  Complex Dfact(_Dconst[imode()]*(prop1*(pdotq2-meta2)+prop2*(pdotq1-meta2)));
  complex<InvEnergy2> Efact(_Econst[imode()]*(prop1+prop2));
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin1,PDT::Spin1);
  Complex e1dote2;
  for(ix=0;ix<3;++ix)
    {
      for(iy=0;iy<3;++iy)
	{
	  if(ix==1||iy==1){newME(0,0,ix,iy)=0.;}
	  else
	    {
	      e1dote2=vwave[0][ix]*vwave[1][iy];
	      newME(0,0,ix,iy) = 
		Dfact*(e1dote2*q1dotq2-e1dotq2[ix]*e2dotq1[iy])
	       -Efact*(-e1dote2*pdotq1*pdotq2-e1dotp[ix]*e2dotp[iy]*q1dotq2
			+e1dotq2[ix]*e2dotp[iy]*pdotq1+e1dotp[ix]*e2dotq1[iy]*pdotq2);
	    }
	}
    }
  // contract the whole thing
  ME(newME);
  RhoDMatrix rhoin(PDT::Spin0);rhoin.average();
  /*
  double me(newME.contract(rhoin).real());
  Energy M(inpart.mass()),M2(M*M);
  Energy2 s1(2.*(decay[1]->momentum()*decay[2]->momentum()));
  Energy2 s2(M2-2.*(inpart.momentum()*decay[1]->momentum()));
  Energy2 s3(M2-2.*(inpart.momentum()*decay[2]->momentum()));
  cout << "testing the matrix element " << (
   2*(2*(Dfact*conj(Dfact)).real() + 2*(Dfact*conj(Efact)).real()*M2 
      + (Efact*conj(Efact)).real()*M2*M2)*
      s1*s1 - 2*(Efact*conj(Efact)).real()*M2*s1*(M2 - s2)*
   (M2 - s3) +(Efact*conj(Efact)).real()*(M2 - s2)*(M2 - s2)*
   (M2-s3)*(M2-s3))/8. - me << endl;
  return me;
  */
  return newME.contract(rhoin).real();
}
 
double EtaPiGammaGammaDecayer::threeBodyMatrixElement(int imodeb,Energy2 q2, Energy2 s3,
						      Energy2 s2,Energy2 s1,
						      Energy m1,Energy m2,Energy m3)
{
  // compute the prefactors
  Energy2 mrho2(_rhomass*_rhomass);
  Energy q(sqrt(s3));
  Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
  Complex ii(0.,1.);
  complex<Energy> tgamma(ii*pcm*pcm*pcm*_rhoconst/q);
  q=sqrt(s2);pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy> ugamma(ii*pcm*pcm*pcm*_rhoconst/q);
  complex<InvEnergy2> prop1(1./(mrho2-s3-tgamma)),prop2(1./(mrho2-s2-ugamma));
  Energy2 pdotq2(0.5*(q2-s3)),pdotq1(0.5*(q2-s2));
  Complex Dfact(_Dconst[imodeb]*(prop1*(pdotq2-q2)+prop2*(pdotq1-q2)));
  complex<InvEnergy2> Efact(_Econst[imodeb]*(prop1+prop2));
  double D2((Dfact*conj(Dfact)).real()),E2((Efact*conj(Efact)).real()),
    ED((Efact*conj(Dfact)).real());
  return (2*(2*D2+2*ED*q2+E2*q2*q2)*s1*s1-2*E2*q2*s1*(q2-s2)*(q2-s3)
		  +E2*(q2-s2)*(q2-s2)*(q2-s3)*(q2-s3))/8.;
}

WidthCalculatorBasePtr 
EtaPiGammaGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int id(dm.parent()->id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
  Energy mrho(getParticleData(ParticleID::rhoplus)->mass());
  Energy wrho(getParticleData(ParticleID::rhoplus)->width());
  vector<double> inmass;inmass.push_back(mrho);inmass.push_back(mrho);
  vector<double> inwidth;inwidth.push_back(wrho);inwidth.push_back(wrho);
  vector<int> intype;intype.push_back(1);intype.push_back(2);
  tcDecayIntegratorPtr decayer=this;
  return new_ptr(ThreeBodyAllOnCalculator(inweights,intype,inmass,inwidth,
					  const_ptr_cast<tDecayIntegratorPtr>(decayer),
					  imode,_mpi,0.,0.));
}

void EtaPiGammaGammaDecayer::dataBaseOutput(ofstream & output) const
{
  output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":Points " << _npoint << "\n";
  output << "set " << fullName() << ":grhoomega " << _grhoomega*GeV << "\n";
  output << "set " << fullName() << ":Fpi " << _fpi/MeV  << "\n";
  output << "set " << fullName() << ":grho " << _grho << "\n";
  output << "set " << fullName() << ":RhoMass " << _rhomass/MeV << "\n";
  output << "set " << fullName() << ":RhoWidth " << _rhowidth/MeV << "\n";
  output << "set " << fullName() << ":RatioFpiF8 " << _ratiofpif8 << "\n";
  output << "set " << fullName() << ":RatioFpiF0 " << _ratiofpif0 << "\n";
  output << "set " << fullName() << ":Theta " << _theta  << "\n";
  output << "set " << fullName() << ":EtaMax " << _etamax << "\n";
  output << "set " << fullName() << ":EtaPrimeMax " << _etapmax << "\n";
  output << "set " << fullName() << ":LocalParameters " << _localparameters << "\n";

  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
