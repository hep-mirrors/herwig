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
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::VectorSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;
using Helicity::VectorWaveFunction;

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
  int id=parent.id(),imode=1;
  if(id==ParticleID::eta){imode=0;}
  bool cc=false;
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
     &EtaPiGammaGammaDecayer::_theta, -pi/9., 0.0, 2.*pi,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaMax
    ("EtaMax",
     "THe maximum weight for the eta decay",
     &EtaPiGammaGammaDecayer::_etamax, 1.2232e-06, -1.0e12, 1.0e12,
     false, false, false);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaPrimeMax
    ("EtaPrimeMax",
     "THe maximum weight for the eta prime decay",
     &EtaPiGammaGammaDecayer::_etapmax, 0.00105024, -1.0e12, 1.0e12,
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
  VectorSpinPtr vspin[2];
  if(vertex)
    {
      decay[0]->spinInfo(new_ptr(ScalarSpinInfo(decay[0]->momentum(),true)));
      vspin[0]=new_ptr(VectorSpinInfo(decay[1]->momentum(),true));
      decay[1]->spinInfo(vspin[0]);
      vspin[1]=new_ptr(VectorSpinInfo(decay[2]->momentum(),true));
      decay[2]->spinInfo(vspin[1]);
    }
  // compute the polarization vectors for the outgoing particles
  VectorWaveFunction vwave1=VectorWaveFunction(decay[1]->momentum(),
					      decay[1]->dataPtr(),outgoing);
  VectorWaveFunction vwave2=VectorWaveFunction(decay[2]->momentum(),
					      decay[2]->dataPtr(),outgoing);
  LorentzPolarizationVector wave1[3],wave2[3];
  // dot products we need
  Energy2 q1dotq2 = decay[1]->momentum()*decay[2]->momentum();
  Energy2 pdotq1  = inpart.momentum()*decay[1]->momentum();
  Energy2 pdotq2  = inpart.momentum()*decay[2]->momentum();
  complex<Energy> e1dotq2[3],e1dotp[3],e2dotq1[3],e2dotp[3];
  for(int ix=-1;ix<2;++ix)
    {
      if(ix==0)
	{
	  wave1[ix+1]=LorentzPolarizationVector();
	  wave2[ix+1]=LorentzPolarizationVector();
	  e1dotq2[ix+1]=0.;
	  e1dotp[ix+1] =0.;
	  e2dotq1[ix+1]=0.;
	  e2dotp[ix+1] =0.;
	}
      else
	{
	  vwave1.reset(ix);wave1[ix+1]=vwave1.Wave();
	  vwave2.reset(ix);wave2[ix+1]=vwave2.Wave();
	  e1dotq2[ix+1] =wave1[ix+1]*decay[2]->momentum();
	  e1dotp[ix+1]  =wave1[ix+1]*inpart.momentum();
	  e2dotq1[ix+1] =wave2[ix+1]*decay[1]->momentum();
	  e2dotp[ix+1]  =wave2[ix+1]*inpart.momentum();
	}
      if(vertex)
	{
	  vspin[0]->setBasisState(ix,wave1[ix+1]);
	  vspin[1]->setBasisState(ix,wave2[ix+1]);
	}
    }
  // the momentum dependent pieces of the matrix element
  Energy2 mpi2=decay[0]->mass()*decay[0]->mass();
  Energy2 meta2=inpart.mass()*inpart.mass();
  Energy2 mrho2=_rhomass*_rhomass;
  Energy2 t = mpi2+2.*((decay[0]->momentum())*(decay[1]->momentum()));
  Energy2 u = mpi2+2.*((decay[0]->momentum())*(decay[2]->momentum()));
  Energy q=sqrt(t);
  Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy> tgamma=Complex(0.,1.)*pcm*pcm*pcm*_rhoconst/q;
  q=sqrt(u);
  pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy> ugamma=Complex(0.,1.)*pcm*pcm*pcm*_rhoconst/q;
  complex<InvEnergy2> prop1=1./(mrho2-t-tgamma);
  complex<InvEnergy2> prop2=1./(mrho2-u-ugamma);
  Complex Dfact = _Dconst[imode()]*(prop1*(pdotq2-meta2)+prop2*(pdotq1-meta2));
  complex<InvEnergy2> Efact= _Econst[imode()]*(prop1+prop2);
  // compute the matrix element
  DecayMatrixElement newME(1,1,3,3);
  Complex e1dote2;
  for(int ihel1=-1;ihel1<2;++ihel1)
    {
      for(int ihel2=-1;ihel2<2;++ihel2)
	{
	  if(ihel1==0||ihel2==0){newME(0,0,ihel1,ihel2)=0.;}
	  else
	    {
	      e1dote2=wave1[ihel1+1]*wave2[ihel2+1];
	      newME(0,0,ihel1,ihel2) = 
		Dfact*(e1dote2*q1dotq2-e1dotq2[ihel1+1]*e2dotq1[ihel2+1])
		-Efact*(-e1dote2*pdotq1*pdotq2-e1dotp[ihel1+1]*e2dotp[ihel2+1]*q1dotq2
			+e1dotq2[ihel1+1]*e2dotp[ihel2+1]*pdotq1
			+e1dotp[ihel1+1]*e2dotq1[ihel2+1]*pdotq2);
	    }
	}
    }
  // contract the whole thing
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  return newME.contract(rhoin).real();
}
 
double EtaPiGammaGammaDecayer::threeBodyMatrixElement(int imodeb,Energy2 q2, Energy2 s3,
						      Energy2 s2,Energy2 s1,
						      Energy m1,Energy m2,Energy m3)
{
  // compute the prefactors
  Energy2 mrho2=_rhomass*_rhomass;
  Energy q=sqrt(s3);
  Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy> tgamma=Complex(0.,1.)*pcm*pcm*pcm*_rhoconst/q;
  q=sqrt(s2);
  pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy> ugamma=Complex(0.,1.)*pcm*pcm*pcm*_rhoconst/q;
  complex<InvEnergy2> prop1=1./(mrho2-s3-tgamma);
  complex<InvEnergy2> prop2=1./(mrho2-s2-ugamma);
  Energy2 pdotq2=0.5*(q2-s3);
  Energy2 pdotq1=0.5*(q2-s2);
  Complex Dfact = _Dconst[imodeb]*(prop1*(pdotq2-q2)+prop2*(pdotq1-q2));
  complex<InvEnergy2> Efact= _Econst[imodeb]*(prop1+prop2);
  double D2 = (Dfact*conj(Dfact)).real();
  double E2 = (Efact*conj(Efact)).real();
  double ED = (Efact*conj(Dfact)).real();
  double output =(2*(2*D2+2*ED*q2+E2*q2*q2)*s1*s1-2*E2*q2*s1*(q2-s2)*(q2-s3)
		  +E2*(q2-s2)*(q2-s2)*(q2-s3)*(q2-s3))/8.;
  return output;
}


WidthCalculatorBasePtr 
EtaPiGammaGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int id=dm.parent()->id(),imode=1;
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
  Energy mrho=getParticleData(ParticleID::rhoplus)->mass();
  Energy wrho=getParticleData(ParticleID::rhoplus)->width();
  vector<double> inmass;inmass.push_back(mrho);inmass.push_back(mrho);
  vector<double> inwidth;inwidth.push_back(wrho);inwidth.push_back(wrho);
  vector<int> intype;intype.push_back(1);intype.push_back(2);
  tcDecayIntegratorPtr decayer=this;
  WidthCalculatorBasePtr 
    output(new_ptr(ThreeBodyAllOnCalculator(inweights,intype,inmass,inwidth,
					    const_ptr_cast<tDecayIntegratorPtr>(decayer),
					    imode,_mpi,0.,0.)));
  return output;
}

}
