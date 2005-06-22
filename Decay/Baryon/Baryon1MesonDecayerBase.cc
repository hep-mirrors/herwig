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
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::SpinorWaveFunction;
  using Herwig::Helicity::ScalarWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::RSSpinorWaveFunction;
using Helicity::RSSpinorBarWaveFunction;
using Helicity::VectorWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

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
      else{throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
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
	  else{throw DecayIntegratorError() << "Unknown outgoing meson spin in "
					    << "Baryon1MesonDecayerBase::me2()" 
					    << Exception::abortnow; }
	}
      // decay to spin-3/2 particle
      else if(decay[0]->dataPtr()->iSpin()==4)
	{
	  // scalar meson
	  if(decay[1]->dataPtr()->iSpin()==1)
	    {me=threeHalfThreeHalfScalar(vertex,ichan,inpart,decay);}
	  else{throw DecayIntegratorError() << "Unknown outgoing meson spin in "
					    << "Baryon1MesonDecayerBase::me2()" 
					    << Exception::abortnow; }
	}
      // unknown
      else{throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
					<< "Baryon1MesonDecayerBase::me2()" 
					<< Exception::abortnow;}
    }
  // unknown
  else{throw DecayIntegratorError() << "Unknown incoming spin in "
				    << "Baryon1MesonDecayerBase::me2()" 
				    << Exception::abortnow;}
  return me;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a pseudoscalar meson
double Baryon1MesonDecayerBase::
halfHalfScalar(bool vertex, const int ichan,const Particle & inpart,
	       const ParticleVector & decay) const
{
  // set up the spins and calculate the spinors
  RhoDMatrix temp;
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inpart.id()>0)
    {
      SpinorWaveFunction(sp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      SpinorBarWaveFunction(sbar,decay[0],outgoing,true,vertex);
    }
  else
    {
      SpinorBarWaveFunction(sbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      SpinorWaveFunction(sp,decay[0],outgoing,true,vertex);
    }
  Herwig::Helicity::ScalarWaveFunction(decay[1],Herwig::Helicity::outgoing,true,vertex);
  // get the couplings
  Complex A,B;
  halfHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),A,B);
  Complex left,right,meout;
  // coupling for an incoming particle
  if(inpart.id()>0){left=(A-B);right=(A+B);}
  // coupling for an incoming antiparticle
  else{left=conj(A+B);right=conj(A-B);}
  // calculate the matrix element
  DecayMatrixElement newME(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  vector<unsigned int> ispin(3,0);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix)
    {
      for(iy=0;iy<2;++iy)
	{
	  if(decay[0]->id()>0){ispin[0]=iy;ispin[1]=ix;}
	  else{ispin[0]=ix;ispin[1]=iy;}
	  newME(ispin)=sp[iy].generalScalar(sbar[ix],left,right);
	}
    }
  // store the matrix element
  ME(newME);
  return (newME.contract(temp)).real()/inpart.mass()/inpart.mass();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfHalfVector(bool vertex, const int ichan,const Particle & inpart,
	       const ParticleVector & decay) const
{
  // check if the outgoing meson is really a photon
  bool photon=decay[1]->id()==ParticleID::gamma;
  // set up the spins and calculate the spinors
  RhoDMatrix temp;
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  vector<LorentzPolarizationVector> eps;
  if(inpart.id()>0)
    {
      SpinorWaveFunction(sp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      SpinorBarWaveFunction(sbar,decay[0],outgoing,true,vertex);
    }
  else
    {
      SpinorBarWaveFunction(sbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      SpinorWaveFunction(sp,decay[0],outgoing,true,vertex);
    }
  // construct the wavefunction and spin info for the vector
  VectorWaveFunction(eps,decay[1],outgoing,true,photon,vertex);
  //for(unsigned int ix=0;ix<eps.size();++ix)
  //  {eps[ix]=LorentzPolarizationVector(decay[1]->momentum());}
  // get the couplings
  Complex A1,A2,B1,B2;
  halfHalfVectorCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			 A1,A2,B1,B2);
  Complex lS,rS,lV,rV,scalar,ii(0.,1.);
  // couplings for an incoming particle
  if(inpart.id()>0){lS=(A2-B2);rS=(A2+B2);lV=(A1-B1);rV=(A1+B1);}
  else{lS=conj(A2+B2);rS=conj(A2-B2);lV=conj(A1-B1);rV=conj(A1+B1);}
  // calculate the matrix element
  // decide which type of mode to do
  Energy msum(inpart.mass()+decay[0]->mass());
  DecayMatrixElement newME(PDT::Spin1Half,decay[0]->dataPtr()->iSpin(),
			   decay[1]->dataPtr()->iSpin());
  vector<unsigned int> ispin(3);
  LorentzPolarizationVector svec;
  Complex s2m4,s1m3,s1p3,s2p4,s3s2,s4s1,s1s4,s2s3,s3s1,s4s2,s1s3,s2s4,prod;
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix)
    {
      for(iy=0;iy<2;++iy)
	{
	  // scalar like piece
	  scalar = sp[iy].generalScalar(sbar[ix],lS,rS);
	  // vector like piece
	  svec   = sp[iy].generalCurrent(sbar[ix],lV,rV);
	  if(decay[0]->id()>0){ispin[0]=iy;ispin[1]=ix;}
	  else{ispin[0]=ix;ispin[1]=iy;}
	  for(ispin[2]=0;ispin[2]<3;++ispin[2])
	    {
	      ispin[2]=ispin[2];
	      prod=eps[ispin[2]]*inpart.momentum();
	      prod/=msum;
	      newME(ispin)=svec*eps[ispin[2]]+prod*scalar;
	      //	      cout << "testing ME " 
	      //	   << ispin[0] << " " << ispin[1] << " " << ispin[2] << " " 
	      //	   << newME(ispin) << endl;
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  return (newME.contract(temp)).real()/inpart.mass()/inpart.mass();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::halfThreeHalfScalar(bool vertex, const int ichan,
						    const Particle & inpart,
						    const ParticleVector & decay) const
{
  unsigned int ix,iy,ixa,iya;
  // set up the spins and calculate the spinors
  RhoDMatrix temp;
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inpart.id()>0)
    {
      SpinorWaveFunction(sp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      vector<LorentzRSSpinorBar> Rsbar;
      RSSpinorBarWaveFunction(Rsbar,decay[0],outgoing,true,vertex);
      sbar.resize(Rsbar.size());
      for(ix=0;ix<Rsbar.size();++ix){sbar[ix]=Rsbar[ix].dot(inpart.momentum());}
    }
  else
    {
      SpinorBarWaveFunction(sbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      vector<LorentzRSSpinor> Rsp;
      RSSpinorWaveFunction(Rsp,decay[0],outgoing,true,vertex);
      sp.resize(Rsp.size());
      for(ix=0;ix<Rsp.size();++ix){sp[ix]=Rsp[ix].dot(inpart.momentum());}
    }
  // construct the spinInfo for the scalar
  ScalarWaveFunction(decay[1],outgoing,true,vertex);
  // get the couplings
  Complex A,B,left,right;
  Energy msum(inpart.mass()+decay[0]->mass());
  halfThreeHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A,B);
  // incoming particle
  if(inpart.id()>0){left=(A-B)/msum;right=(A+B)/msum;}
  // incoming anti-particle
  else{left=conj(A+B)/msum;right=conj(A-B)/msum;}
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin0);
  vector<unsigned int> ispin(3,0);
  for(ixa=0;ixa<2;++ixa)
    {
      for(iya=0;iya<4;++iya)
	{
	  if(decay[0]->id()>0){ix=iya;iy=ixa;}
	  else{ix=ixa;iy=iya;}
	  // low energy conventions
	  ispin[0]=ixa;ispin[1]=iya;
	  newME(ispin)=sp[iy].generalScalar(sbar[ix],left,right);
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  return (newME.contract(temp)).real()/inpart.mass()/inpart.mass();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfThreeHalfVector(bool vertex, const int ichan,const Particle & inpart,
	       const ParticleVector & decay) const
{
  // check if the outgoing meson is really a photon
  bool photon=decay[1]->id()==ParticleID::gamma;
  // set up the spins and calculate the spinors
  RhoDMatrix temp;
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  vector<LorentzRSSpinor> RSsp;
  vector<LorentzRSSpinorBar> RSsbar;
  vector<LorentzPolarizationVector> eps;
  if(inpart.id()>0)
    {
      SpinorWaveFunction(sp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      RSSpinorBarWaveFunction(RSsbar,decay[0],outgoing,true,vertex);
      sbar.resize(RSsbar.size());
      for(unsigned int ix=0;ix<RSsbar.size();++ix)
	{sbar[ix]=RSsbar[ix].dot(inpart.momentum());}
    }
  else
    {
      SpinorBarWaveFunction(sbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      RSSpinorWaveFunction(RSsp,decay[0],outgoing,true,vertex);
      sp.resize(RSsp.size());
      for(unsigned int ix=0;ix<RSsp.size();++ix)
	{sp[ix]=RSsp[ix].dot(inpart.momentum());}
    }
  // construct the wavefunction and spin info for the vector
  VectorWaveFunction(eps,decay[1],outgoing,true,photon,vertex);
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3;
  halfThreeHalfVectorCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(inpart.mass()+decay[0]->mass());
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
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin1);
  vector<unsigned int> ispin(3);
  LorentzPolarizationVector svec;
  Complex prod,meout;
  unsigned int ix,iy,ixa,iya,iz;
  LorentzSpinor stemp;
  LorentzSpinorBar sbtemp;
  for(iya=0;iya<4;++iya)
    {
      ispin[1]=iya;
      // piece where the vector-spinor is dotted with the momentum of the
      // incoming fermion
      for(ixa=0;ixa<2;++ixa)
	{
	  if(decay[0]->id()>0){ix=iya;iy=ixa;}
	  else{ix=ixa;iy=iya;}
	  scalar = sp[iy].generalScalar(sbar[ix],lS,rS);
	  svec   = sp[iy].generalCurrent(sbar[ix],lV,rV);
	  ispin[0]=ixa;
	  for(iz=0;iz<3;++iz)
	    {
	      ispin[2]=iz;
	      prod=eps[iz]*inpart.momentum();
	      prod/=msum;
	      newME(ispin)+=(svec*eps[iz]+prod*scalar)/msum;
	    }
	}
      // the piece where the vector spinor is dotted with the polarization vector
      for(unsigned int iz=0;iz<3;++iz)
	{
	  ispin[2]=iz;
	  if(decay[0]->id()>0){sbtemp=RSsbar[iya].dot(eps[iz]);}
	  else{stemp=RSsp[iya].dot(eps[iz]);}
	  for(ixa=0;ixa<2;++ixa)
	    {
	      ispin[0]=ixa;
	      if(decay[0]->id()>0){stemp=sp[ixa];}
	      else{sbtemp=sbar[ixa];}
	      meout = stemp.generalScalar(sbtemp,left,right);
	      newME(ispin)+=meout;
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  return (newME.contract(temp)).real()/inpart.mass()/inpart.mass();
}


// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
threeHalfHalfScalar(bool vertex, const int ichan,const Particle & inpart,
		    const ParticleVector & decay) const
{
  // set up the spins and calculate the spinors
  RhoDMatrix temp;
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  if(inpart.id()>0)
    {
      vector<LorentzRSSpinor> Rsp;
      RSSpinorWaveFunction(Rsp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      SpinorBarWaveFunction(sbar,decay[0],outgoing,true,vertex);
      sp.resize(Rsp.size());
      for(unsigned int ix=0;ix<Rsp.size();++ix)
	{sp[ix]=Rsp[ix].dot(decay[0]->momentum());}
    }
  else
    {
      vector<LorentzRSSpinorBar> Rsbar;
      RSSpinorBarWaveFunction(Rsbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      SpinorWaveFunction(sp,decay[0],outgoing,true,vertex);
      sbar.resize(Rsbar.size());
      for(unsigned int ix=0;ix<Rsbar.size();++ix)
	{sbar[ix]=Rsbar[ix].dot(decay[0]->momentum());}
    }
  // construct the spinInfo for the scalar
  ScalarWaveFunction(decay[1],outgoing,true,vertex);
  // get the couplings
  Complex A,B;
  Energy msum=inpart.mass()+decay[0]->mass();
  threeHalfHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A,B);
  Complex left,right;
  // incoming particle
  if(inpart.id()>0){left=(A-B)/msum;right=(A+B)/msum;}
  // incoming anti-particle
  else{left=conj(A+B)/msum;right=conj(A-B)/msum;}
  // compute the matrix element
  vector<unsigned int> ispin(3,0);
  DecayMatrixElement newME(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin0);
  unsigned int ix,iy,ixa,iya;
  for(ixa=0;ixa<2;++ixa)
    {
      for(iya=0;iya<4;++iya)
	{
	  if(decay[0]->id()>0){iy=iya;ix=ixa;}
	  else{iy=ixa;ix=iya;}
	  ispin[0]=iya;ispin[1]=ixa;
	  newME(ispin) = sp[iy].generalScalar(sbar[ix],left,right);
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  return (newME.contract(temp)).real()/inpart.mass()/inpart.mass();
}

// matrix element for the decay of a spin-3/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::threeHalfThreeHalfScalar(bool vertex, const int ichan,
							 const Particle & inpart,
							 const ParticleVector & decay) const
{
  RhoDMatrix temp;
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  vector<LorentzRSSpinor> Rsp;
  vector<LorentzRSSpinorBar> Rsbar;
  if(inpart.id()>0)
    {
      RSSpinorWaveFunction(Rsp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
      RSSpinorBarWaveFunction(Rsbar,decay[0],outgoing,true,vertex);
      for(unsigned int ix=0;ix<4;++ix)
	{
	  sp.push_back(Rsp[ix].dot(decay[0]->momentum()));
	  sbar.push_back(Rsbar[ix].dot(inpart.momentum()));
	}
    }
  else
    {
      RSSpinorBarWaveFunction(Rsbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			      vertex);
      RSSpinorWaveFunction(Rsp,decay[0],outgoing,true,vertex);
      for(unsigned int ix=0;ix<4;++ix)
	{
	  sbar.push_back(Rsbar[ix].dot(decay[0]->momentum()));
	  sp.push_back(Rsp[ix].dot(inpart.momentum()));
	}
    }
  // construct the spinInfo for the scalar
  ScalarWaveFunction(decay[1],outgoing,true,vertex);
  // get the couplings
  Complex A1,B1,A2,B2;
  Energy msum(inpart.mass()+decay[0]->mass());
  threeHalfThreeHalfScalarCoupling(imode(),inpart.mass(),decay[0]->mass(),
				   decay[1]->mass(),A1,A2,B1,B2);
  Complex left1,right1,left2,right2;
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
  DecayMatrixElement newME(PDT::Spin3Half,PDT::Spin3Half,PDT::Spin0);
  vector<unsigned int> ispin(3,0);
  unsigned int ix,iy,ixa,iya;
  for(ixa=0;ixa<4;++ixa)
    {
      for(iya=0;iya<4;++iya)
	{
	  if(decay[0]->id()>0){iy=iya;ix=ixa;}
	  else{iy=ixa;ix=iya;}
	  ispin[0]=iya;
	  ispin[1]=ixa;
	  newME(ispin)=Rsp[iy].generalScalar(Rsbar[ix],left1,right1)
	               +sp[iy].generalScalar( sbar[ix],left2,right2);;
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  return (newME.contract(temp)).real()/inpart.mass()/inpart.mass();
}

// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a vector meson

double Baryon1MesonDecayerBase::
threeHalfHalfVector(bool vertex, const int ichan,const Particle & inpart,
 const ParticleVector & decay) const
{
  // check if the outgoing meson is really a photon
  bool photon=decay[1]->id()==ParticleID::gamma;
  // set up the spins and calculate the spinors
  RhoDMatrix temp;
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  vector<LorentzRSSpinor> RSsp;
  vector<LorentzRSSpinorBar> RSsbar;
  vector<LorentzPolarizationVector> eps;
  if(inpart.id()>0)
    {
      RSSpinorWaveFunction(RSsp,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			   vertex);
      SpinorBarWaveFunction(sbar,decay[0],outgoing,true,vertex);
      sp.resize(RSsp.size());
      for(unsigned int ix=0;ix<RSsp.size();++ix)
	{sp[ix]=RSsp[ix].dot(decay[0]->momentum());}
    }
  else
    {
      RSSpinorBarWaveFunction(RSsbar,temp,const_ptr_cast<tPPtr>(&inpart),incoming,true,
			    vertex);
      SpinorWaveFunction(sp,decay[0],outgoing,true,vertex);
      sbar.resize(RSsbar.size());
      for(unsigned int ix=0;ix<RSsbar.size();++ix)
	{sbar[ix]=RSsbar[ix].dot(decay[0]->momentum());}
    }
  // construct the wavefunction and spin info for the vector
  VectorWaveFunction(eps,decay[1],outgoing,true,photon,vertex);
  //for(unsigned int ix=0;ix<eps.size();++ix)
  //  {eps[ix]=LorentzPolarizationVector(decay[1]->momentum());}
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3,prod,meout;
  threeHalfHalfVectorCoupling(imode(),inpart.mass(),decay[0]->mass(),decay[1]->mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(inpart.mass()+decay[0]->mass());
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
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin1);
  vector<unsigned int> ispin(3);
  LorentzPolarizationVector svec;
  unsigned int ix,iy,ixa,iya,iz;
  LorentzSpinor stemp;
  LorentzSpinorBar sbtemp;
  for(iya=0;iya<4;++iya)
    {
      ispin[0]=iya;
      for(ixa=0;ixa<2;++ixa)
	{
	  if(decay[0]->id()>0){iy=iya;ix=ixa;}
	  else{iy=ixa;ix=iya;}
	  scalar = sp[iy].generalScalar( sbar[ix],lS,rS);
	  svec   = sp[iy].generalCurrent(sbar[ix],lV,rV);
	  ispin[1]=ixa;
	  for(iz=0;iz<3;++iz)
	    {
	      ispin[2]=iz;
	      prod=eps[iz]*decay[0]->momentum();
	      prod/=msum;
	      newME(ispin)+=(svec*eps[iz]+prod*scalar)/msum;
	    }
	}
      // the piece where the vector spinor is dotted with the polarization vector
      for(iz=0;iz<3;++iz)
	{
	  ispin[2]=iz;
	  if(decay[0]->id()>0){stemp=RSsp[iya].dot(eps[iz]);}
	  else{sbtemp=RSsbar[iya].dot(eps[iz]);}
	  for(ixa=0;ixa<2;++ixa)
	    {
	      ispin[1]=ixa;
	      if(decay[0]->id()>0){sbtemp=sbar[ixa];}
	      else{stemp=sp[ixa];}
	      meout = stemp.generalScalar(sbtemp,left,right);
	      newME(ispin)+=meout;
	    }
	}
    }
  // store the matrix element
  ME(newME);
  // return the answer
  return (newME.contract(temp)).real()/inpart.mass()/inpart.mass();
}


void Baryon1MesonDecayerBase::halfHalfScalarCoupling(int,Energy,Energy,Energy,
						     Complex&,Complex&) const
{throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfScalarCoupling()"
			      << " called from base class this must be implemented"
			      << " in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::halfHalfVectorCoupling(int,Energy,Energy,Energy,
						     Complex&,Complex&,
						     Complex&,Complex&) const
{throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfVectorCoupling()" 
			      << " called from base class this must be implemented " 
			      << "in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&) const
{throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling"
			      << "() called from base class this must be implemented"
			      << " in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&,
							  Complex&,Complex&,
							  Complex&,Complex&) const
{throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling"
			      << "() called from base class this must be implemented "
			      << "in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::threeHalfHalfScalarCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&) const
{throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfHalfScalarCoupling"
			      << "() called from base class this must be implemented"
			      << " in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::threeHalfHalfVectorCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&,
							  Complex&,Complex&,
							  Complex&,Complex&) const
{throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfHalfVectorCoupling"
			      << "() called from base class this must be implemented "
			      << "in the inheriting class" << Exception::abortnow;}

void Baryon1MesonDecayerBase::threeHalfThreeHalfScalarCoupling(int,Energy,Energy,
							       Energy,Complex&,
							       Complex&,Complex&,
							       Complex&) const
{throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfThreeHalfScalar"
			      << "Coupling() called from base class this must be "
			      << "implemented in the inheriting class" 
			      << Exception::abortnow;}

bool Baryon1MesonDecayerBase::twoBodyMEcode(const DecayMode & dm,int & mecode,
					    double & coupling) const
{
  coupling=1.;
  unsigned int inspin(dm.parent()->iSpin()),outspin,outmes;
  ParticleMSet::const_iterator pit(dm.products().begin());
  bool order; 
  if((**pit).iSpin()%2==0)
    {
      order=true;
      outspin=(**pit).iSpin();
      ++pit;outmes=(**pit).iSpin();
    }
  else
    {
      order=false;
      outmes=(**pit).iSpin();++pit;
      outspin=(**pit).iSpin();
    }
  mecode=-1;
  if(inspin==2)
    {
      if(outspin==2){if(outmes==1){mecode=101;}else{mecode=102;}}
      else if(outspin==4){if(outmes==1){mecode=103;}else{mecode=104;}}
    }
  else if(inspin==4)
    {
      if(outspin==2){if(outmes==1){mecode=105;}else{mecode=106;}}
      else if(outspin==4){if(outmes==1){mecode=107;}else{mecode=108;}}
    }
  return order;
}

}



