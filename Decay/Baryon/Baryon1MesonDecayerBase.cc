// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Baryon1MesonDecayerBase class.
//

#include "Baryon1MesonDecayerBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Baryon1MesonDecayerBase.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcFermionSpinPtr;
using ThePEG::Helicity::tcRSFermionSpinPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::RSSpinorWaveFunction;
using Helicity::RSSpinorBarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

Baryon1MesonDecayerBase::~Baryon1MesonDecayerBase() {}

void Baryon1MesonDecayerBase::persistentOutput(PersistentOStream & os) const {
}

void Baryon1MesonDecayerBase::persistentInput(PersistentIStream & is, int) {
}

AbstractClassDescription<Baryon1MesonDecayerBase> Baryon1MesonDecayerBase::initBaryon1MesonDecayerBase;
// Definition of the static class description member.

void Baryon1MesonDecayerBase::Init() {

  static ClassDocumentation<Baryon1MesonDecayerBase> documentation
    ("The \\classname{Baryon1MesonDecayerBase} class is the base class for"
     " the decays of the baryons to a baryon and a pseudoscalar or vector meson.");

}

// return the matrix element squared for a given mode and phase-space channel
// (inherited from DecayIntegrator and implemented here)
double Baryon1MesonDecayerBase::me2(bool vertex, const int ichan,
				    const Particle & inpart,
				    const ParticleVector & decay) const
{
  double me(0.);
  // decide which matrix element we are doing
  // incoming spin-1/2 particle
  if(inpart.dataPtr()->iSpin()==2)
    {
      // decay to spin-1/2 particle
      if(decay[0]->dataPtr()->iSpin()==2)
	{
	  // scalar meson
	  if(decay[1]->dataPtr()->iSpin()==1)
	    {me=halfHalfScalar(vertex,ichan,inpart,decay);}
	  // vector meson
	  else if(decay[1]->dataPtr()->iSpin()==3)
	    {me=halfHalfVector(vertex,ichan,inpart,decay);}
	  else
	    {throw DecayIntegratorError() << "Unknown outgoing meson spin in "
					  << "Baryon1MesonDecayerBase::me2()" 
					  << Exception::abortnow; }
	}
      // decay to spin-3/2 particle
      else if(decay[0]->dataPtr()->iSpin()==4)
	{
	  // scalar meson
	  if(decay[1]->dataPtr()->iSpin()==1)
	    {me=halfThreeHalfScalar(vertex,ichan,inpart,decay);}
	  // vector meson
	  else if(decay[1]->dataPtr()->iSpin()==3)
	    {me=halfThreeHalfVector(vertex,ichan,inpart,decay);}
	  else
	    {throw DecayIntegratorError() << "Unknown outgoing meson spin in "
					  << "Baryon1MesonDecayerBase::me2()" 
					  << Exception::abortnow; }
	}
      // unknown
      else
	{throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
				       << "Baryon1MesonDecayerBase::me2()" 
				       << Exception::abortnow;}
    }
  // incoming spin-3/2 particle
  else if(inpart.dataPtr()->iSpin()==4)
    {
      // decay to spin-1/2 particle
      if(decay[0]->dataPtr()->iSpin()==2)
	{
	  // scalar meson
	  if(decay[1]->dataPtr()->iSpin()==1)
	    {me=threeHalfHalfScalar(vertex,ichan,inpart,decay);}
	  // vector meson
	  else if(decay[1]->dataPtr()->iSpin()==3)
	    {me=threeHalfHalfVector(vertex,ichan,inpart,decay);}
	  else
	    {throw DecayIntegratorError() << "Unknown outgoing meson spin in "
					  << "Baryon1MesonDecayerBase::me2()" 
					  << Exception::abortnow; }
	}
      // decay to spin-3/2 particle
      else if(decay[0]->dataPtr()->iSpin()==4)
	{
	  // scalar meson
	  if(decay[1]->dataPtr()->iSpin()==1)
	    {me=threeHalfThreeHalfScalar(vertex,ichan,inpart,decay);}
	  else
	    {throw DecayIntegratorError() << "Unknown outgoing meson spin in "
					  << "Baryon1MesonDecayerBase::me2()" 
					  << Exception::abortnow; }
	}
      // unknown
      else
	{throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
				       << "Baryon1MesonDecayerBase::me2()" 
				       << Exception::abortnow;}
    }
  // unknown
  else
    {
      throw DecayIntegratorError() << "Unknown incoming spin in "
				   << "Baryon1MesonDecayerBase::me2()" 
				   << Exception::abortnow;
    }
  return me;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a pseudoscalar meson
double Baryon1MesonDecayerBase::
halfHalfScalar(bool vertex, const int ichan,const Particle & inpart,
	       const ParticleVector & decay) const
{
  // check if the decay particle has spin info 
  tcFermionSpinPtr inspin;
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  sp.push_back(inspin->getDecayBasisState(-1));
	  sp.push_back(inspin->getDecayBasisState( 1));
	}
      else
	{
	  sbar.push_back(inspin->getDecayBasisState(-1).bar());
	  sbar.push_back(inspin->getDecayBasisState( 1).bar());
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in Baryon1MesonDecayerBase::halfHalfScalar()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(FermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  SpinorWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sp.push_back(temp.Wave());
	  temp.reset(1);sp.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
		{inspin->setDecayState(ix,sp[(ix+1)/2]);}}
	  
	}
      else
	{
	  SpinorBarWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sbar.push_back(temp.Wave());
	  temp.reset(1);sbar.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
	      {inspin->setDecayState(ix,sbar[(ix+1)/2].bar());}}
	}
    }
  tcFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr temp(new_ptr(FermionSpinInfo(decay[0]->momentum(),true)));
      outspin=dynamic_ptr_cast<tcFermionSpinPtr>(temp);
      decay[0]->spinInfo(temp);
    }
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      SpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
				 -1,outgoing,dirac);
      sbar.push_back(temp.Wave());
      temp.reset(1);sbar.push_back(temp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2)
	    {outspin->setBasisState(ix,sbar[(ix+1)/2].bar());}}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      SpinorWaveFunction temp = SpinorWaveFunction(decay[0]->momentum(),
						   decay[0]->dataPtr(),-1,
						   outgoing,dirac);
      sp.push_back(temp.Wave());
      temp.reset(1);sp.push_back(temp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2)
	    {outspin->setBasisState(ix,sp[(ix+1)/2]);}}
    }
  // construct the spinInfo for the scalar
  decay[1]->spinInfo(new_ptr(ScalarSpinInfo(decay[1]->momentum(),true)));
  // get the couplings
  Complex A,B;
  halfHalfScalarCoupling(imode(),A,B);
  Complex left,right,meout;
  // coupling for an incoming particle
  if(inpart.id()>0){left=(A-B);right=(A+B);}
  // coupling for an incoming antiparticle
  else{left=conj(A+B);right=conj(A-B);}
  // calculate the matrix element
  vector<int> ispin;
  ispin.push_back(decay[0]->dataPtr()->iSpin());
  ispin.push_back(decay[1]->dataPtr()->iSpin());
  DecayMatrixElement newME(2,ispin);
  ispin.push_back(0);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix)
    {
      for(iy=0;iy<2;++iy)
	{
	  // low energy conventions
	  if(sp[iy].Rep()==HaberDRep&&sbar[ix].Rep()==HaberDRep)
	    {
	      meout = left*((sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
			    +(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4()))
		+right*( (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
			 +(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4()));
	      meout*=0.5;
	    }
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    {
	      meout=  left*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  // mixed conventions
	  else
	    {
	      sp[iy].changeRep(HELASDRep);
	      sbar[ix].changeRep(HELASDRep);
	      meout=  left*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  if(decay[0]->id()>0){ispin[0]=2*iy-1;ispin[1]=2*ix-1;}
	  else{ispin[0]=2*ix-1;ispin[1]=2*iy-1;}
	  newME(ispin)=meout;
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  double me=(newME.contract(temp)).real()/inpart.mass()/inpart.mass();
  /*
  // calculate the matrix element using the KK results
  Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
  Energy2 qplus  = (m1+m2)*(m1+m2)-m3*m3;
  Energy2 qminus = (m1-m2)*(m1-m2)-m3*m3;
  Energy rqplus=sqrt(qplus),rqminus=sqrt(qminus); 
  Complex h1,h2;
  // the amplitudes
  h1 =  2.*rqplus *A;
  h2 = -2.*rqminus*B;
  cout << "testing A " << 0.25*(h1*conj(h1)+h2*conj(h2))/(me*m1*m1) << endl;
  */
  return me;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfHalfVector(bool vertex, const int ichan,const Particle & inpart,
	       const ParticleVector & decay) const
{
  // check if the outgoing meson is really a photon
  bool photon=decay[1]->id()==ParticleID::gamma;
  // check if the decay particle has spin info 
  tcFermionSpinPtr inspin;
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  sp.push_back(inspin->getDecayBasisState(-1));
	  sp.push_back(inspin->getDecayBasisState( 1));
	}
      else
	{
	  sbar.push_back(inspin->getDecayBasisState(-1).bar());
	  sbar.push_back(inspin->getDecayBasisState( 1).bar());
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in Baryon1MesonDecayerBase::halfHalfVector()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(FermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  SpinorWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sp.push_back(temp.Wave());
	  temp.reset(1);sp.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
		{inspin->setDecayState(ix,sp[(ix+1)/2]);}}
	}
      else
	{
	  SpinorBarWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sbar.push_back(temp.Wave());
	  temp.reset(1);sbar.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
	      {inspin->setDecayState(ix,sbar[(ix+1)/2].bar());}}
	}
    }
  tcFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr temp(new_ptr(FermionSpinInfo(decay[0]->momentum(),true)));
      outspin=dynamic_ptr_cast<tcFermionSpinPtr>(temp);
      decay[0]->spinInfo(temp);
    }
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      SpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
				 -1,outgoing,dirac);
      sbar.push_back(temp.Wave());
      temp.reset(1);sbar.push_back(temp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2)
	    {outspin->setBasisState(ix,sbar[(ix+1)/2].bar());}}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      SpinorWaveFunction temp = SpinorWaveFunction(decay[0]->momentum(),
						   decay[0]->dataPtr(),-1,
						   outgoing,dirac);
      sp.push_back(temp.Wave());
      temp.reset(1);sp.push_back(temp.Wave());
      if(vertex)
	{for(int ix=-1;ix<2;ix+=2)
		{outspin->setBasisState(ix,sp[(ix+1)/2]);}}
    }
  // construct the wavefunction and spin info for the vector
  VectorWaveFunction vtemp(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  SpinPtr stemp=new_ptr(VectorSpinInfo(decay[1]->momentum(),true));
  tcVectorSpinPtr vspin=dynamic_ptr_cast<tcVectorSpinPtr>(stemp);
  decay[1]->spinInfo(stemp);
  // get the couplings
  Complex A1,A2,B1,B2;
  halfHalfVectorCoupling(imode(),A1,A2,B1,B2);
  Complex lS,rS,lV,rV,scalar,ii(0.,1.);
  // couplings for an incoming particle
  if(inpart.id()>0){lS=(A2-B2);rS=(A2+B2);lV=(A1-B1);rV=(A1+B1);}
  else{lS=conj(A2+B2);rS=conj(A2-B2);lV=conj(A1-B1);rV=conj(A1+B1);}
  // wavefunctions for the outgoing vector
  vector<LorentzPolarizationVector> eps;
  for(int ix=-1;ix<2;++ix)
    {
      if(!(photon&&ix==0))
	{
	  vtemp.reset(ix);
	  if(vertex){vspin->setBasisState(ix,vtemp.Wave());}
	  eps.push_back(vtemp.Wave());
	}
      else
	{
	  if(vertex){vspin->setBasisState(ix,LorentzPolarizationVector());}
	  eps.push_back(LorentzPolarizationVector());
	}
    }
  // calculate the matrix element
  // decide which type of mode to do
  vector<int> ispin;
  ispin.push_back(decay[0]->dataPtr()->iSpin());
  ispin.push_back(decay[1]->dataPtr()->iSpin());
  Energy msum=inpart.mass()+decay[0]->mass();
  DecayMatrixElement newME(2,ispin);
  ispin.resize(3);
  LorentzPolarizationVector svec;
  Complex s2m4,s1m3,s1p3,s2p4,s3s2,s4s1,s1s4,s2s3,s3s1,s4s2,s1s3,s2s4,prod;
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(unsigned int iy=0;iy<2;++iy)
	{
	  // spinor part of the matrix element
	  // low energy conventions
	  if(sp[iy].Rep()==HaberDRep&&sbar[ix].Rep()==HaberDRep)
	    {
	      // scalar like term
	      scalar = lS*( (sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
			    +(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4()))
		+rS*( (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
		      +(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4()));
	      scalar*=0.5;
	      // vector like term
	      s2m4=sp[iy].s2()-sp[iy].s4();
	      s1m3=sp[iy].s1()-sp[iy].s3();
	      s1p3=sp[iy].s1()+sp[iy].s3();
	      s2p4=sp[iy].s2()+sp[iy].s4();
	      svec[0] =   0.5*lV*(-sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  -sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+0.5*rV*(+sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			 -sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[1] =ii*0.5*lV*(+sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  +sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+ii*0.5*rV*(-sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			    +sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[2] =   0.5*lV*(-sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  -sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3-sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3+sbar[ix].s4()*s2p4);
	      svec[3] =   0.5*lV*(+sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  +sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3+sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3-sbar[ix].s4()*s2p4);
	    }
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    {
	      // scalar like term
	      scalar=  lS*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+rS*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	      // vector like term
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      svec[0] =    -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] = ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =    -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =     lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	    }
	  // mixed conventions
	  else
	    {
	      // scalar like term
	      sp[iy].changeRep(HELASDRep);
	      sbar[ix].changeRep(HELASDRep);
	      scalar=  lS*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+rS*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      // vector term
	      svec[0] =   -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] =ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =   -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =    lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	    }
	  if(decay[0]->id()>0){ispin[0]=2*iy-1;ispin[1]=2*ix-1;}
	  else{ispin[0]=2*ix-1;ispin[1]=2*iy-1;}
	  for(unsigned int iz=0;iz<3;++iz)
	    {
	      ispin[2]=iz-1;
	      prod=eps[iz]*inpart.momentum();
	      prod/=msum;
	      newME(ispin)=svec*eps[iz]+prod*scalar;
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  double me=(newME.contract(temp)).real()/inpart.mass()/inpart.mass();
  /*
  // calculate the matrix element using the KK results
  Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
  Energy2 qplus  = (m1+m2)*(m1+m2)-m3*m3;
  Energy2 qminus = (m1-m2)*(m1-m2)-m3*m3;
  Energy rqplus=sqrt(qplus),rqminus=sqrt(qminus); 
  Complex h1,h2,h3,h4;
  Energy pcm = sqrt((m1*m1-(m2+m3)*(m2+m3))*(m1*m1-(m2-m3)*(m2-m3)))/2./m1;
  A2 /=msum;
  B2 /=msum;
  // the amplitudes
  h1 =  2.*sqrt(2.)*rqplus *B1;
  h2 = -2.*sqrt(2.)*rqminus*A1;
  h3 = 2./m3*(rqplus *(m1-m2)*B1-rqminus*m1*pcm*B2);
  h4 = 2./m3*(rqminus*(m1+m2)*A1+rqplus *m1*pcm*A2);
  cout << "testing B " << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)
		    +h4*conj(h4))/(me*m1*m1) << endl;
  */
  return me;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
halfThreeHalfScalar(bool vertex, const int ichan,const Particle & inpart,
	       const ParticleVector & decay) const
{
  // check if the decay particle has spin info 
  tcFermionSpinPtr inspin;
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  sp.push_back(inspin->getDecayBasisState(-1));
	  sp.push_back(inspin->getDecayBasisState( 1));
	}
      else
	{
	  sbar.push_back(inspin->getDecayBasisState(-1).bar());
	  sbar.push_back(inspin->getDecayBasisState( 1).bar());
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in Baryon1MesonDecayerBase::"
				  << "halfThreeHalfScalar()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(FermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  SpinorWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sp.push_back(temp.Wave());
	  temp.reset(1);sp.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
		{inspin->setDecayState(ix,sp[(ix+1)/2]);}}
	  
	}
      else
	{
	  SpinorBarWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sbar.push_back(temp.Wave());
	  temp.reset(1);sbar.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
	      {inspin->setDecayState(ix,sbar[(ix+1)/2].bar());}}
	}
    }
  tcRSFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr temp(new_ptr(RSFermionSpinInfo(decay[0]->momentum(),true)));
      outspin=dynamic_ptr_cast<tcRSFermionSpinPtr>(temp);
      decay[0]->spinInfo(temp);
    }
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      RSSpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),-2,
				   outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sbar.push_back(temp.Wave().dot(inpart.momentum()));
	      if(vertex){outspin->setBasisState(ix,temp.Wave().bar());}
	    }
	}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      RSSpinorWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),-2,
				outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sp.push_back(temp.Wave().dot(inpart.momentum()));
	      if(vertex){outspin->setBasisState(ix,temp.Wave());}
	    }
	}
    }
  // construct the spinInfo for the scalar
  decay[1]->spinInfo(new_ptr(ScalarSpinInfo(decay[1]->momentum(),true)));
  // get the couplings
  Complex A,B;
  Energy msum=inpart.mass()+decay[0]->mass();
  halfThreeHalfScalarCoupling(imode(),A,B);
  Complex left,right,meout;
  // incoming particle
  if(inpart.id()>0){left=(A-B)/msum;right=(A+B)/msum;}
  // incoming anti-particle
  else{left=conj(A+B)/msum;right=conj(A-B)/msum;}
  // compute the matrix element
  vector<int> ispin;
  ispin.push_back(decay[0]->dataPtr()->iSpin());
  ispin.push_back(decay[1]->dataPtr()->iSpin());
  DecayMatrixElement newME(2,ispin);
  ispin.push_back(0);
  unsigned int ix,iy,ixa,iya;
  for(ixa=0;ixa<2;++ixa)
    {
      for(iya=0;iya<4;++iya)
	{
	  if(decay[0]->id()>0){ix=iya;iy=ixa;}
	  else{ix=ixa;iy=iya;}
	  // low energy conventions
	  if(sp[iy].Rep()==HaberDRep&&sbar[ix].Rep()==HaberDRep)
	    {
	      meout = 
		left*( (sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
		       +(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4()))
		+right*( (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
			 +(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4()));
	      meout*=0.5;
	    } 
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    {
	      meout=  
		left*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  // mixed conventions
	  else
	    {
	      sp[iy].changeRep(HELASDRep);
	      sbar[ix].changeRep(HELASDRep);
	      meout=  
		left*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  ispin[0]=2*ixa-1;ispin[1]=iya-2;if(ispin[1]>=0){++ispin[1];}
	  newME(ispin)=meout;
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  double me=(newME.contract(temp)).real()/inpart.mass()/inpart.mass();
  /*
  // calculate the matrix element using the KK results
  Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
  Energy2 qplus  = (m1+m2)*(m1+m2)-m3*m3;
  Energy2 qminus = (m1-m2)*(m1-m2)-m3*m3;
  Energy rqplus=sqrt(qplus),rqminus=sqrt(qminus); 
  Complex h1,h2;
  Energy pcm = sqrt((m1*m1-(m2+m3)*(m2+m3))*(m1*m1-(m2-m3)*(m2-m3)))/2./m1;
  A /=msum;
  B /=msum;
  // the amplitudes
  h1 =-2.*sqrt(2./3.)*pcm*m1/m2*rqminus*B;
  h2 = 2.*sqrt(2./3.)*pcm*m1/m2*rqplus *A;
  cout << "testing C " << 0.25*(h1*conj(h1)+h2*conj(h2))/(me*m1*m1) << endl;
  */
  return me;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfThreeHalfVector(bool vertex, const int ichan,const Particle & inpart,
	       const ParticleVector & decay) const
{
  // check if the decay particle has spin info 
  tcFermionSpinPtr inspin;
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector <LorentzSpinorBar> sbar;
  vector<LorentzRSSpinor> RSsp;
  vector<LorentzRSSpinorBar> RSsbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  sp.push_back(inspin->getDecayBasisState(-1));
	  sp.push_back(inspin->getDecayBasisState( 1));
	}
      else
	{
	  sbar.push_back(inspin->getDecayBasisState(-1).bar());
	  sbar.push_back(inspin->getDecayBasisState( 1).bar());
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in Baryon1MesonDecayerBase::" 
				  << "halfThreeHalfVector()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(FermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  SpinorWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sp.push_back(temp.Wave());
	  temp.reset(1);sp.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
		{inspin->setDecayState(ix,sp[(ix+1)/2]);}}
	  
	}
      else
	{
	  SpinorBarWaveFunction temp(inpart.momentum(),inpart.dataPtr(),-1,incoming);
	  sbar.push_back(temp.Wave());
	  temp.reset(1);sbar.push_back(temp.Wave());
	  if(vertex)
	    {for(int ix=-1;ix<2;ix+=2)
	      {inspin->setDecayState(ix,sbar[(ix+1)/2].bar());}}
	}
    }
  tcRSFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr temp(new_ptr(RSFermionSpinInfo(decay[0]->momentum(),true)));
      outspin=dynamic_ptr_cast<tcRSFermionSpinPtr>(temp);
      decay[0]->spinInfo(temp);
    }
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      RSSpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),-2,
				   outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sbar.push_back(temp.Wave().dot(inpart.momentum()));
	      RSsbar.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave().bar());}
	    }
	}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      RSSpinorWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),-2,
				outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sp.push_back(temp.Wave().dot(inpart.momentum()));
	      RSsp.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave());}
	    }
	}
    }
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3;
  halfThreeHalfVectorCoupling(imode(),A1,A2,A3,B1,B2,B3);
  Energy msum=inpart.mass()+decay[0]->mass();
  Complex lS,rS,lV,rV,scalar,ii(0.,1.),left,right;
  // incoming particle
  if(inpart.id()>0)
    {
      lS=(A3-B3);rS=(A3+B3);
      lV=(A2-B2);rV=(A2+B2);
      left=(A1-B1);right=(A1+B1);
    }
  // incoming anti-particle
  else
    {
      lS=conj(A3+B3);rS=conj(A3-B3);
      lV=conj(A2-B2);rV=conj(A2+B2);
      left=conj(A1+B1);right=conj(A1-B1);
    }
  // construct the wavefunction and spin info for the vector
  VectorWaveFunction vtemp(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  SpinPtr itemp=new_ptr(VectorSpinInfo(decay[1]->momentum(),true));
  tcVectorSpinPtr vspin=dynamic_ptr_cast<tcVectorSpinPtr>(itemp);
  decay[1]->spinInfo(itemp);
  // wavefunctions for the outgoing vector
  vector<LorentzPolarizationVector> eps;
  for(int ix=-1;ix<2;++ix)
    {
      vtemp.reset(ix);
      if(vertex){vspin->setBasisState(ix,vtemp.Wave());}
      eps.push_back(vtemp.Wave());
    }
  // compute the matrix element
  vector<int> ispin;
  ispin.push_back(decay[0]->dataPtr()->iSpin());
  ispin.push_back(decay[1]->dataPtr()->iSpin());
  DecayMatrixElement newME(2,ispin);
  ispin.resize(3);
  LorentzPolarizationVector svec;
  Complex s2m4,s1m3,s1p3,s2p4,s3s2,s4s1,s1s4,s2s3,s3s1,s4s2,s1s3,s2s4,prod,meout;
  unsigned int ix,iy,ixa,iya,iz;
  LorentzSpinor stemp;
  LorentzSpinorBar sbtemp;
  for(iya=0;iya<4;++iya)
    {
      ispin[1]=iya-2;if(ispin[1]>=0){++ispin[1];}
      // piece where the vector-spinor is dotted with the momentum of the
      // incoming fermion
      for(ixa=0;ixa<2;++ixa)
	{
	  if(decay[0]->id()>0){ix=iya;iy=ixa;}
	  else{ix=ixa;iy=iya;}
	  if(sp[iy].Rep()==HaberDRep&&sbar[ix].Rep()==HaberDRep)
	    {
	      // scalar like term
	      scalar = lS*( (sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
			    +(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4()))
		      +rS*( (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
			    +(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4()));
	      scalar*=0.5;
	      // vector like term
	      s2m4=sp[iy].s2()-sp[iy].s4();
	      s1m3=sp[iy].s1()-sp[iy].s3();
	      s1p3=sp[iy].s1()+sp[iy].s3();
	      s2p4=sp[iy].s2()+sp[iy].s4();
	      svec[0] =   0.5*lV*(-sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  -sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+0.5*rV*(+sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			 -sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[1] =ii*0.5*lV*(+sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  +sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+ii*0.5*rV*(-sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			    +sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[2] =   0.5*lV*(-sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  -sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3-sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3+sbar[ix].s4()*s2p4);
	      svec[3] =   0.5*lV*(+sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  +sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3+sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3-sbar[ix].s4()*s2p4);
	    }
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    {
	      // scalar like term
	      scalar=  lS*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		      +rS*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	      // vector like term
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      svec[0] =    -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] = ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =    -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =     lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	    }
	  // mixed conventions
	  else
	    {
	      // scalar like term
	      sp[iy].changeRep(HELASDRep);
	      sbar[ix].changeRep(HELASDRep);
	      scalar=  lS*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+rS*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      // vector term
	      svec[0] =   -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] =ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =   -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =    lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	    }
	  ispin[0]=2*ixa-1;
	  for(iz=0;iz<3;++iz)
	    {
	      ispin[2]=iz-1;
	      prod=eps[iz]*inpart.momentum();
	      prod/=msum;
	      newME(ispin)+=(svec*eps[iz]+prod*scalar)/msum;
	    }
	}
      // the piece where the vector spinor is dotted with the polarization vector
      for(unsigned int iz=0;iz<3;++iz)
	{
	  ispin[2]=iz-1;
	  if(decay[0]->id()>0){sbtemp=RSsbar[iya].dot(eps[iz]);}
	  else{stemp=RSsp[iya].dot(eps[iz]);}
	  for(ixa=0;ixa<2;++ixa)
	    {
	      ispin[0]=2*ixa-1;
	      if(decay[0]->id()>0){stemp=sp[ixa];}
	      else{sbtemp=sbar[ixa];}
	      // low energy conventions
	      if(stemp.Rep()==HaberDRep&&sbtemp.Rep()==HaberDRep)
		{
		  meout = left*((sbtemp.s1()-sbtemp.s3())*(stemp.s1()-stemp.s3())
				+(sbtemp.s2()-sbtemp.s4())*(stemp.s2()-stemp.s4()))
		    +right*( (sbtemp.s1()+sbtemp.s3())*(stemp.s1()+stemp.s3())
			     +(sbtemp.s2()+sbtemp.s4())*(stemp.s2()+stemp.s4()));
		  meout*=0.5;
		}
	      // high energy conventions
	      else if(stemp.Rep()==HELASDRep&&sbtemp.Rep()==HELASDRep)
		{
		  meout=  left*(sbtemp.s1()*stemp.s1()+sbtemp.s2()*stemp.s2())
		    +right*(sbtemp.s3()*stemp.s3()+sbtemp.s4()*stemp.s4());
		}
	      // mixed conventions
	      else
		{
		  stemp.changeRep(HELASDRep);
		  sbtemp.changeRep(HELASDRep);
		  meout=  left*(sbtemp.s1()*stemp.s1()+sbtemp.s2()*stemp.s2())
		    +right*(sbtemp.s3()*stemp.s3()+sbtemp.s4()*stemp.s4());
		}
	      newME(ispin)+=meout;
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  double me=(newME.contract(temp)).real()/inpart.mass()/inpart.mass();
  /*
  // calculate the matrix element using the KK results
  Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
  Energy2 qplus  = (m1+m2)*(m1+m2)-m3*m3;
  Energy2 qminus = (m1-m2)*(m1-m2)-m3*m3;
  Energy rqplus=sqrt(qplus),rqminus=sqrt(qminus); 
  Complex h1,h2,h3,h4,h5,h6;
  Energy pcm = sqrt((m1*m1-(m2+m3)*(m2+m3))*(m1*m1-(m2-m3)*(m2-m3)))/2./m1;
  A2 /=msum;
  B2 /=msum;
  A3 /=(msum*msum);
  B3 /=(msum*msum);
  // the amplitudes
  h1 = -2.*rqplus *A1;
  h2 =  2.*rqminus*B1;
  h3 = -2./sqrt(3.)*rqplus *(A1-qminus/m2*A2);
  h4 =  2./sqrt(3.)*rqminus*(B1-qplus /m2*B2);
  h5 =  2.*sqrt(2./3.)/m2/m3*(-rqplus *(0.5*(m1*m1-m2*m2-m3*m3)*A1
					+0.5*qminus*(m1+m2)*A2+m1*m1*pcm*pcm*A3));
  h6 =  2.*sqrt(2./3.)/m2/m3*( rqminus*(0.5*(m1*m1-m2*m2-m3*m3)*B1
					-0.5*qplus *(m1-m2)*B2+m1*m1*pcm*pcm*A3));
  cout << "testing D " << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)
		    +h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/(me*m1*m1) << endl;
  */
  return me;
}


// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
threeHalfHalfScalar(bool vertex, const int ichan,const Particle & inpart,
		    const ParticleVector & decay) const
{
  // check if the decay particle has spin info 
  tcRSFermionSpinPtr inspin;
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcRSFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{sp.push_back(inspin->getDecayBasisState(ix).dot(decay[0]->momentum()));}
	    }
	}
      else
	{
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{sbar.push_back(inspin->getDecayBasisState(ix).bar().dot(decay[0]->momentum()));}
	    }
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in Baryon1MesonDecayerBase::"
				  << "halfThreeHalfScalar()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(RSFermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcRSFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  RSSpinorWaveFunction temp(inpart.momentum(),inpart.dataPtr(),incoming);
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{
		  temp.reset(ix);
		  sp.push_back(temp.Wave().dot(decay[0]->momentum()));
		  if(vertex){inspin->setDecayState(ix,temp.Wave());}
		}
	    }
	}
      else
	{
	  RSSpinorBarWaveFunction temp(inpart.momentum(),inpart.dataPtr(),incoming);
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{
		  temp.reset(ix);
		  sbar.push_back(temp.Wave().dot(decay[0]->momentum()));
		  if(vertex){inspin->setDecayState(ix,temp.Wave().bar());}
		}
	    }
	}
    }
  tcFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr temp(new_ptr(FermionSpinInfo(decay[0]->momentum(),true)));
      outspin=dynamic_ptr_cast<tcFermionSpinPtr>(temp);
      decay[0]->spinInfo(temp);
    }
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      SpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
				 outgoing,dirac);
      for(int ix=-1;ix<2;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sbar.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave().bar());}
	    }
	}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      SpinorWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
			      outgoing,dirac);
      for(int ix=-1;ix<2;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sp.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave());}
	    }
	}
    }
  // construct the spinInfo for the scalar
  decay[1]->spinInfo(new_ptr(ScalarSpinInfo(decay[1]->momentum(),true)));
  // get the couplings
  Complex A,B;
  Energy msum=inpart.mass()+decay[0]->mass();
  halfThreeHalfScalarCoupling(imode(),A,B);
  Complex left,right,meout;
  // incoming particle
  if(inpart.id()>0){left=(A-B)/msum;right=(A+B)/msum;}
  // incoming anti-particle
  else{left=conj(A+B)/msum;right=conj(A-B)/msum;}
  // compute the matrix element
  vector<int> ispin;
  ispin.push_back(decay[0]->dataPtr()->iSpin());
  ispin.push_back(decay[1]->dataPtr()->iSpin());
  DecayMatrixElement newME(4,ispin);
  ispin.push_back(0);
  unsigned int ix,iy,ixa,iya;
  Complex mesum(0.);
  for(ixa=0;ixa<2;++ixa)
    {
      for(iya=0;iya<4;++iya)
	{
	  if(decay[0]->id()>0){iy=iya;ix=ixa;}
	  else{iy=ixa;ix=iya;}
	  // low energy conventions
	  if(sp[iy].Rep()==HaberDRep&&sbar[ix].Rep()==HaberDRep)
	    {
	      meout = 
		left*( (sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
		       +(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4()))
		+right*( (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
			 +(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4()));
	      meout*=0.5;
	    } 
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    {
	      //	      cout << "testing doing the helas bit " << ix << "  " << iy << endl;
	      //cout << "testing bar  " 
	      //   << sbar[ix][0] << " " 
	      //   << sbar[ix][1] << " " 
	      //	   << sbar[ix][2] << " " 
	      //   << sbar[ix][3] << " " << endl;
	      //cout << "testing spin " 
	      //   << sp[iy][0] << " " 
	      //	   << sp[iy][1] << " " 
	      //   << sp[iy][2] << " " 
	      //	   << sp[iy][3] << " " << endl;
	      meout=  
		left*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  // mixed conventions
	  else
	    {
	      sp[iy].changeRep(HELASDRep);
	      sbar[ix].changeRep(HELASDRep);
	      meout=  
		left*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  ispin[0]=iya-2;ispin[1]=2*ixa-1;if(ispin[0]>=0){++ispin[0];}
	  //cout << "testing " 
	  //    << ispin[0] << "  " 
	  //    << ispin[1] << "  " 
	  //    << ispin[2] << "  " << meout << endl;
	  newME(ispin)=meout;
	  mesum +=meout*conj(meout);
	}
    }
  //  cout << "testing the matrix element " << mesum/4./inpart.mass()/inpart.mass() << endl;
  // store the matrix element
  ME(newME);
  // return the answer
  double me=(newME.contract(temp)).real()/inpart.mass()/inpart.mass();
  // cout << "testing the matrix element " << me << " " << inpart.mass() << " " 
  //    << decay[0]->mass() << " " << decay[1]->mass() << endl;
  return me;
}

// matrix element for the decay of a spin-3/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::threeHalfThreeHalfScalar(bool vertex, const int ichan,
							 const Particle & inpart,
							 const ParticleVector & decay) const
{
  // check if the decay particle has spin info 
  tcRSFermionSpinPtr inspin;
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcRSFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  vector<LorentzRSSpinor> Rsp;
  vector<LorentzRSSpinorBar> Rsbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{
		  Rsp.push_back(inspin->getDecayBasisState(ix));
		  sp.push_back(Rsp.back().dot(decay[0]->momentum()));
		}
	    }
	}
      else
	{
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{
		  Rsbar.push_back(inspin->getDecayBasisState(ix).bar());
		  sbar.push_back(Rsbar.back().dot(decay[0]->momentum()));
		}
	    }
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in Baryon1MesonDecayerBase::"
				  << "halfThreeHalfScalar()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(RSFermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcRSFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  RSSpinorWaveFunction temp(inpart.momentum(),inpart.dataPtr(),incoming);
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{
		  temp.reset(ix);
		  sp.push_back(temp.Wave().dot(decay[0]->momentum()));
		  Rsp.push_back(temp.Wave());
		  if(vertex){inspin->setDecayState(ix,temp.Wave());}
		}
	    }
	}
      else
	{
	  RSSpinorBarWaveFunction temp(inpart.momentum(),inpart.dataPtr(),incoming);
	  for(int ix=-2;ix<3;++ix)
	    {
	      if(ix!=0)
		{
		  temp.reset(ix);
		  sbar.push_back(temp.Wave().dot(decay[0]->momentum()));
		  Rsbar.push_back(temp.Wave());
		  if(vertex){inspin->setDecayState(ix,temp.Wave().bar());}
		}
	    }
	}
    }
  tcRSFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr temp(new_ptr(RSFermionSpinInfo(decay[0]->momentum(),true)));
      outspin=dynamic_ptr_cast<tcRSFermionSpinPtr>(temp);
      decay[0]->spinInfo(temp);
    }
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      RSSpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
				   outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      Rsbar.push_back(temp.Wave());
	      sbar.push_back(temp.Wave().dot(inpart.momentum()));
	      if(vertex){outspin->setBasisState(ix,temp.Wave().bar());}
	    }
	}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      RSSpinorWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
			      outgoing,dirac);
      for(int ix=-2;ix<3;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      Rsp.push_back(temp.Wave());
	      sp.push_back(temp.Wave().dot(inpart.momentum()));
	      if(vertex){outspin->setBasisState(ix,temp.Wave());}
	    }
	}
    }
  // construct the spinInfo for the scalar
  decay[1]->spinInfo(new_ptr(ScalarSpinInfo(decay[1]->momentum(),true)));
  // get the couplings
  Complex A1,B1,A2,B2;
  Energy msum=inpart.mass()+decay[0]->mass();
  threeHalfThreeHalfScalarCoupling(imode(),A1,A2,B1,B2);
  Complex left1,right1,left2,right2,meout;
  // incoming particle
  if(inpart.id()>0)
    {
      left1=(A1-B1);           right1=(A1+B1);
      left2=(A2-B2)/msum/msum; right2=(A2+B2)/msum/msum;
    }
  // incoming anti-particle
  else
    {
      left1=(A1+B1);           right1=(A1-B1);
      left2=(A2+B2)/msum/msum; right2=(A2-B2)/msum/msum;
    }
  // compute the matrix element
  vector<int> ispin;
  ispin.push_back(decay[0]->dataPtr()->iSpin());
  ispin.push_back(decay[1]->dataPtr()->iSpin());
  DecayMatrixElement newME(4,ispin);
  ispin.push_back(0);
  unsigned int ix,iy,ixa,iya,iz;
  Complex mesum(0.);
  for(ixa=0;ixa<4;++ixa)
    {
      for(iya=0;iya<4;++iya)
	{
	  if(decay[0]->id()>0){iy=iya;ix=ixa;}
	  else{iy=ixa;ix=iya;}
	  // low energy conventions
	  if(sp[iy].Rep()==HaberDRep&&sbar[ix].Rep()==HaberDRep)
	    {
	      meout=
		left1*( (Rsbar[ix](3,0)-Rsbar[ix](3,2))*(Rsp[iy](3,0)-Rsp[iy](3,2))
			+(Rsbar[ix](3,1)-Rsbar[ix](3,3))*(Rsp[iy](3,1)-Rsp[iy](3,3)))
		+right1*( (Rsbar[ix](3,0)+Rsbar[ix](3,2))*(Rsp[iy](3,0)+Rsp[iy](3,2))
			  +(Rsbar[ix](3,1)+Rsbar[ix](3,3))*(Rsp[iy](3,1)+Rsp[iy](3,3)));
	      // the first piece
	      for(iz=0;iz<3;++iz)
		{
		  meout-=
		    left1*( (Rsbar[ix](iz,0)-Rsbar[ix](iz,2))*(Rsp[iy](iz,0)-Rsp[iy](iz,2))
			    +(Rsbar[ix](iz,1)-Rsbar[ix](iz,3))*(Rsp[iy](iz,1)-Rsp[iy](iz,3)))
		    +right1*( (Rsbar[ix](iz,0)+Rsbar[ix](iz,2))*(Rsp[iy](iz,0)+Rsp[iy](iz,2))
			      +(Rsbar[ix](iz,1)+Rsbar[ix](iz,3))*(Rsp[iy](iz,1)+Rsp[iy](iz,3)));
		}
	      // the second piece
	      meout+=
		left2*( (sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
			+(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4()))
		+right2*( (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
			  +(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4()));
	      meout*=0.5;
	    }
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    {
	      //  first piece 
	      meout = 
		left1*(Rsbar[ix](3,0)*Rsp[iy](3,0)+Rsbar[ix](3,1)*Rsp[iy](3,1))
		+right1*(Rsbar[ix](3,2)*Rsp[iy](3,2)+Rsbar[ix](3,3)*Rsp[iy](3,3));
	      for(iz=0;iz<3;++iz)
		{
		  meout-=
		    left1*(Rsbar[ix](iz,0)*Rsp[iy](iz,0)+Rsbar[ix](iz,1)*Rsp[iy](iz,1))
		    +right1*(Rsbar[ix](iz,2)*Rsp[iy](iz,2)+Rsbar[ix](iz,3)*Rsp[iy](iz,3));
		}
	      // second piece
	      meout+=
		left2*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right2*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  // mixed conventions
	  else
	    {
	      sp[iy].changeRep(HELASDRep);
	      sbar[ix].changeRep(HELASDRep);
	      Rsp[iy].changeRep(HELASDRep);
	      Rsbar[ix].changeRep(HELASDRep);
	      //  first piece 
	      meout = 
		left1*(Rsbar[ix](3,0)*Rsp[iy](3,0)+Rsbar[ix](3,1)*Rsp[iy](3,1))
		+right1*(Rsbar[ix](3,2)*Rsp[iy](3,2)+Rsbar[ix](3,3)*Rsp[iy](3,3));
	      for(iz=0;iz<3;++iz)
		{
		  meout-=
		    left1*(Rsbar[ix](iz,0)*Rsp[iy](iz,0)+Rsbar[ix](iz,1)*Rsp[iy](iz,1))
		    +right1*(Rsbar[ix](iz,2)*Rsp[iy](iz,2)+Rsbar[ix](iz,3)*Rsp[iy](iz,3));
		}
	      // second piece
	      meout+=
		left2*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+right2*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	    }
	  ispin[0]=iya-2;if(ispin[0]>=0){++ispin[0];}
	  ispin[1]=ixa-2;if(ispin[1]>=0){++ispin[1];}
	  newME(ispin)=meout;
	  mesum +=meout*conj(meout);
	}
    }  
  // store the matrix element
  ME(newME);
  // return the answer
  double me=(newME.contract(temp)).real()/inpart.mass()/inpart.mass();
  return me;
}

// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a vector meson

double Baryon1MesonDecayerBase::
threeHalfHalfVector(bool vertex, const int ichan,const Particle & inpart,
 const ParticleVector & decay) const
{
  bool photon=decay[1]->id()==ParticleID::gamma;
  // check if the decay particle has spin info 
  tcRSFermionSpinPtr inspin;
  RhoDMatrix temp(inpart.dataPtr()->iSpin()); temp.average();
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcRSFermionSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  vector<LorentzRSSpinor> RSsp;
  vector<LorentzRSSpinorBar> RSsbar;
  if(inspin)
    {
      inspin->decay();
      temp=inspin->rhoMatrix();
      if(inpart.id()>0)
	{
	  for(int iz=-2;iz<3;++iz)
	    {
	      if(iz!=0)
		{
		  RSsp.push_back(inspin->getDecayBasisState(iz));
		  sp.push_back(RSsp.back().dot(decay[0]->momentum()));
		}
	    }
	}
      else
	{
	  for(int iz=-2;iz<3;++iz)
	    {
	      if(iz!=0)
		{
		  RSsbar.push_back(inspin->getDecayBasisState(iz).bar());
		  sbar.push_back(RSsbar.back().dot(decay[0]->momentum()));
		}
	    }
	}
    }
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in Baryon1MesonDecayerBase::" 
				  << "halfThreeHalfVector()" 
				  << Exception::abortnow;}
  // spininfo does not exist create it
  else
    {
      SpinPtr newspin=new_ptr(RSFermionSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcRSFermionSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
      if(inpart.id()>0)
	{
	  RSSpinorWaveFunction temp(inpart.momentum(),inpart.dataPtr(),incoming);
	  for(int iz=-2;iz<3;++iz)
	    {
	      if(iz!=0)
		{
		  temp.reset(iz);
		  RSsp.push_back(temp.Wave());
		  sp.push_back(temp.Wave().dot(decay[0]->momentum()));
		  if(vertex){inspin->setDecayState(iz,temp.Wave());}
		}
	    }
	}
      else
	{
	  RSSpinorBarWaveFunction temp(inpart.momentum(),inpart.dataPtr(),incoming);
	  for(int iz=-2;iz<3;++iz)
	    {
	      if(iz!=0)
		{
		  temp.reset(iz);
		  RSsbar.push_back(temp.Wave());
		  sbar.push_back(temp.Wave().dot(decay[0]->momentum()));
		  if(vertex){inspin->setDecayState(iz,temp.Wave().bar());}
		}
	    }
	}
    }
  tcFermionSpinPtr outspin;
  // construct the spin info for the outgoing baryon
  if(vertex)
    {
      SpinPtr temp(new_ptr(FermionSpinInfo(decay[0]->momentum(),true)));
      outspin=dynamic_ptr_cast<tcFermionSpinPtr>(temp);
      decay[0]->spinInfo(temp);
    }
  // calculate the wavefunction for the outgoing particles
  if(decay[0]->id()>0)
    {
      DiracRep dirac=sp[0].Rep();
      SpinorBarWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
				 outgoing,dirac);
      for(int ix=-1;ix<2;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sbar.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave().bar());}
	    }
	}
    }
  else
    {
      DiracRep dirac=sbar[0].Rep();
      SpinorWaveFunction temp(decay[0]->momentum(),decay[0]->dataPtr(),
			      outgoing,dirac);
      for(int ix=-1;ix<2;++ix)
	{
	  if(ix!=0)
	    {
	      temp.reset(ix);
	      sp.push_back(temp.Wave());
	      if(vertex){outspin->setBasisState(ix,temp.Wave());}
	    }
	}
    }
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3;
  halfThreeHalfVectorCoupling(imode(),A1,A2,A3,B1,B2,B3);
  Energy msum=inpart.mass()+decay[0]->mass();
  Complex lS,rS,lV,rV,scalar,ii(0.,1.),left,right;
  // incoming particle
  if(inpart.id()>0)
    {
      lS=(A3-B3);rS=(A3+B3);
      lV=(A2-B2);rV=(A2+B2);
      left=(A1-B1);right=(A1+B1);
    }
  // incoming anti-particle
  else
    {
      lS=conj(A3+B3);rS=conj(A3-B3);
      lV=conj(A2-B2);rV=conj(A2+B2);
      left=conj(A1+B1);right=conj(A1-B1);
    }
  // construct the wavefunction and spin info for the vector
  VectorWaveFunction vtemp(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  SpinPtr itemp=new_ptr(VectorSpinInfo(decay[1]->momentum(),true));
  tcVectorSpinPtr vspin=dynamic_ptr_cast<tcVectorSpinPtr>(itemp);
  decay[1]->spinInfo(itemp);
  // wavefunctions for the outgoing vector
  vector<LorentzPolarizationVector> eps;
  for(int ix=-1;ix<2;++ix)
    {
      if(!(photon&&ix==0))
	{
	  vtemp.reset(ix);
	  if(vertex){vspin->setBasisState(ix,vtemp.Wave());}
	  eps.push_back(vtemp.Wave());
	}
      else
	{
	  if(vertex){vspin->setBasisState(ix,LorentzPolarizationVector());}
	  eps.push_back(LorentzPolarizationVector());
	}
    }
  // compute the matrix element
  vector<int> ispin;
  ispin.push_back(decay[0]->dataPtr()->iSpin());
  ispin.push_back(decay[1]->dataPtr()->iSpin());
  DecayMatrixElement newME(4,ispin);
  ispin.resize(3);
  LorentzPolarizationVector svec;
  Complex s2m4,s1m3,s1p3,s2p4,s3s2,s4s1,s1s4,s2s3,s3s1,s4s2,s1s3,s2s4,prod,meout;
  unsigned int ix,iy,ixa,iya,iz;
  LorentzSpinor stemp;
  LorentzSpinorBar sbtemp;
  for(iya=0;iya<4;++iya)
    {
      ispin[0]=iya-2;if(ispin[0]>=0){++ispin[0];}
      for(ixa=0;ixa<2;++ixa)
	{
	  if(decay[0]->id()>0){iy=iya;ix=ixa;}
	  else{iy=ixa;ix=iya;}
	  if(sp[iy].Rep()==HaberDRep&&sbar[ix].Rep()==HaberDRep)
	    {
	      // scalar like term
	      scalar = lS*( (sbar[ix].s1()-sbar[ix].s3())*(sp[iy].s1()-sp[iy].s3())
			    +(sbar[ix].s2()-sbar[ix].s4())*(sp[iy].s2()-sp[iy].s4()))
		      +rS*( (sbar[ix].s1()+sbar[ix].s3())*(sp[iy].s1()+sp[iy].s3())
			    +(sbar[ix].s2()+sbar[ix].s4())*(sp[iy].s2()+sp[iy].s4()));
	      scalar*=0.5;
	      // vector like term
	      s2m4=sp[iy].s2()-sp[iy].s4();
	      s1m3=sp[iy].s1()-sp[iy].s3();
	      s1p3=sp[iy].s1()+sp[iy].s3();
	      s2p4=sp[iy].s2()+sp[iy].s4();
	      svec[0] =   0.5*lV*(-sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  -sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+0.5*rV*(+sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			 -sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[1] =ii*0.5*lV*(+sbar[ix].s1()*s2m4-sbar[ix].s2()*s1m3
				  +sbar[ix].s3()*s2m4-sbar[ix].s4()*s1m3)
		+ii*0.5*rV*(-sbar[ix].s1()*s2p4+sbar[ix].s2()*s1p3
			    +sbar[ix].s3()*s2p4-sbar[ix].s4()*s1p3);
	      svec[2] =   0.5*lV*(-sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  -sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3-sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3+sbar[ix].s4()*s2p4);
	      svec[3] =   0.5*lV*(+sbar[ix].s1()*s1m3+sbar[ix].s2()*s2m4
				  +sbar[ix].s3()*s1m3+sbar[ix].s4()*s2m4)
		+0.5*rV*(+sbar[ix].s1()*s1p3+sbar[ix].s2()*s2p4
			 -sbar[ix].s3()*s1p3-sbar[ix].s4()*s2p4);
	    }
	  // high energy conventions
	  else if(sp[iy].Rep()==HELASDRep&&sbar[ix].Rep()==HELASDRep)
	    {
	      // scalar like term
	      scalar=  lS*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		      +rS*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	      // vector like term
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      svec[0] =    -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] = ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =    -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =     lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	    }
	  // mixed conventions
	  else
	    {
	      // scalar like term
	      sp[iy].changeRep(HELASDRep);
	      sbar[ix].changeRep(HELASDRep);
	      scalar=  lS*(sbar[ix].s1()*sp[iy].s1()+sbar[ix].s2()*sp[iy].s2())
		+rS*(sbar[ix].s3()*sp[iy].s3()+sbar[ix].s4()*sp[iy].s4());
	      s3s2=sbar[ix].s3()*sp[iy].s2();s4s1=sbar[ix].s4()*sp[iy].s1();
	      s1s4=sbar[ix].s1()*sp[iy].s4();s2s3=sbar[ix].s2()*sp[iy].s3();
	      s3s1=sbar[ix].s3()*sp[iy].s1();s4s2=sbar[ix].s4()*sp[iy].s2();
	      s1s3=sbar[ix].s1()*sp[iy].s3();s2s4=sbar[ix].s2()*sp[iy].s4();
	      // vector term
	      svec[0] =   -lV*(s3s2+s4s1)+rV*(s1s4+s2s3);
	      svec[1] =ii*(lV*(s3s2-s4s1)-rV*(s1s4-s2s3));
	      svec[2] =   -lV*(s3s1-s4s2)+rV*(s1s3-s2s4);
	      svec[3] =    lV*(s3s1+s4s2)+rV*(s1s3+s2s4);
	    }
	  ispin[1]=2*ixa-1;
	  for(iz=0;iz<3;++iz)
	    {
	      ispin[2]=iz-1;
	      prod=eps[iz]*decay[0]->momentum();
	      prod/=msum;
	      newME(ispin)+=(svec*eps[iz]+prod*scalar)/msum;
	    }
	}
      // the piece where the vector spinor is dotted with the polarization vector
      for(unsigned int iz=0;iz<3;++iz)
	{
	  ispin[2]=iz-1;
	  if(decay[0]->id()>0){stemp=RSsp[iya].dot(eps[iz]);}
	  else{sbtemp=RSsbar[iya].dot(eps[iz]);}
	  for(ixa=0;ixa<2;++ixa)
	    {
	      ispin[1]=2*ixa-1;
	      if(decay[0]->id()>0){sbtemp=sbar[ixa];}
	      else{stemp=sp[ixa];}
	      // low energy conventions
	      if(stemp.Rep()==HaberDRep&&sbtemp.Rep()==HaberDRep)
		{
		  meout = left*((sbtemp.s1()-sbtemp.s3())*(stemp.s1()-stemp.s3())
				+(sbtemp.s2()-sbtemp.s4())*(stemp.s2()-stemp.s4()))
		    +right*( (sbtemp.s1()+sbtemp.s3())*(stemp.s1()+stemp.s3())
			     +(sbtemp.s2()+sbtemp.s4())*(stemp.s2()+stemp.s4()));
		  meout*=0.5;
		}
	      // high energy conventions
	      else if(stemp.Rep()==HELASDRep&&sbtemp.Rep()==HELASDRep)
		{
		  meout=  left*(sbtemp.s1()*stemp.s1()+sbtemp.s2()*stemp.s2())
		    +right*(sbtemp.s3()*stemp.s3()+sbtemp.s4()*stemp.s4());
		}
	      // mixed conventions
	      else
		{
		  stemp.changeRep(HELASDRep);
		  sbtemp.changeRep(HELASDRep);
		  meout=  left*(sbtemp.s1()*stemp.s1()+sbtemp.s2()*stemp.s2())
		    +right*(sbtemp.s3()*stemp.s3()+sbtemp.s4()*stemp.s4());
		}
	      newME(ispin)+=meout;
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  double me=(newME.contract(temp)).real()/inpart.mass()/inpart.mass();
  /*
  // calculate the matrix element using the KK results
  Energy m1(inpart.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
  Energy2 qplus  = (m1+m2)*(m1+m2)-m3*m3;
  Energy2 qminus = (m1-m2)*(m1-m2)-m3*m3;
  Energy rqplus=sqrt(qplus),rqminus=sqrt(qminus); 
  Complex h1,h2,h3,h4;
  Energy pcm = sqrt((m1*m1-(m2+m3)*(m2+m3))*(m1*m1-(m2-m3)*(m2-m3)))/2./m1;
  A2 /=msum;
  B2 /=msum;
  A3 /=(msum*msum);
  B3 /=(msum*msum);
  // the amplitudes
  h1 = -2.*rqplus *A1;
  h2 =  2.*rqminus*B1;
  h3 = -2./sqrt(3.)*rqplus *(A1+qminus/m1*A2);
  h4 =  2./sqrt(3.)*rqminus*(B1+qplus /m1*B2);
  cout << "testing the matrix element " << me << " " << m1 << " " << m2 << endl;
  cout << "testing D " << 0.125*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)
		    +h4*conj(h4))/(me*m1*m1) << endl;
*/
  return me;
}


void Baryon1MesonDecayerBase::halfHalfScalarCoupling(int,Complex&,Complex&) const
 {throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfScalarCoupling() called from base class this must be implemented in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::halfHalfVectorCoupling(int,Complex&,Complex&,
						     Complex&,Complex&) const
 {throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfVectorCoupling() called from base class this must be implemented in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling(int,Complex&,
							  Complex&) const
 {throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling() called from base class this must be implemented in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling(int,Complex&,Complex&,
							  Complex&,Complex&,
							  Complex&,Complex&) const
 {throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling() called from base class this must be implemented in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::threeHalfThreeHalfScalarCoupling(int,Complex&,
							       Complex&,Complex&,
							       Complex&) const
 {throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfThreeHalfScalarCoupling() called from base class this must be implemented in the inheriting class" << Exception::abortnow;}


}



