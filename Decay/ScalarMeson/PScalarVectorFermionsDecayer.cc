// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorFermionsDecayer class.
//

#include "PScalarVectorFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarVectorFermionsDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "Herwig++/PDT/ThreeBodyAllOn1IntegralCalculator.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VectorWaveFunction;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::EpsFunction;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::tcVectorSpinPtr;
using ThePEG::Helicity::tcFermionSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::defaultDRep;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

PScalarVectorFermionsDecayer::~PScalarVectorFermionsDecayer() {}

bool PScalarVectorFermionsDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  // must be three outgoing particles
  if(dm.products().size()!=3){return allowed;}
  // ids of the particles
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit; int id2=(**pit).id();
  ++pit; int id3=(**pit).id();
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0)
      {
	if((id1==_outgoingV[ix]&&id2==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id1==_outgoingV[ix]&&id3==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id2==_outgoingV[ix]&&id1==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id2==_outgoingV[ix]&&id3==_outgoingf[ix]&&id1==_outgoinga[ix])||
	   (id3==_outgoingV[ix]&&id1==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id3==_outgoingV[ix]&&id2==_outgoingf[ix]&&id1==_outgoinga[ix]))
	  {allowed=true;}
      }
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector PScalarVectorFermionsDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // workout which mode we are doing
  int imode=-1;
  int id0=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit; int id2=(**pit).id();
  ++pit; int id3=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0)
	{
	  if((id1==_outgoingV[ix]&&id2==_outgoingf[ix]&&id3==_outgoinga[ix])||
	     (id1==_outgoingV[ix]&&id3==_outgoingf[ix]&&id2==_outgoinga[ix])||
	     (id2==_outgoingV[ix]&&id1==_outgoingf[ix]&&id3==_outgoinga[ix])||
	     (id2==_outgoingV[ix]&&id3==_outgoingf[ix]&&id1==_outgoinga[ix])||
	     (id3==_outgoingV[ix]&&id1==_outgoingf[ix]&&id2==_outgoinga[ix])||
	     (id3==_outgoingV[ix]&&id2==_outgoingf[ix]&&id1==_outgoinga[ix]))
	    {imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  return generate(false,false,imode,parent);
}


void PScalarVectorFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingV << _outgoingf << _outgoinga << _maxweight
     << _includeVMD << _VMDid << _VMDmass << _VMDwidth;
}

void PScalarVectorFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingV >> _outgoingf >> _outgoinga >> _maxweight
     >> _includeVMD >> _VMDid >> _VMDmass >> _VMDwidth;
}

ClassDescription<PScalarVectorFermionsDecayer> PScalarVectorFermionsDecayer::initPScalarVectorFermionsDecayer;
// Definition of the static class description member.

void PScalarVectorFermionsDecayer::Init() {

  static ClassDocumentation<PScalarVectorFermionsDecayer> documentation
    ("The \\classname{PScalarVectorFermionsDecayer} class is designed"
     " for the decay of a pseudoscalar meson to a photon and a"
     "fermion-antifermion pair");

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarVectorFermionsDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingV
    ("OutgoingVector",
     "The PDG code for the outgoing pseudoscalar",
     &PScalarVectorFermionsDecayer::_outgoingV,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingA
    ("OutgoingAntiFermion",
     "The PDG code for the outgoing antifermion",
     &PScalarVectorFermionsDecayer::_outgoinga,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarVectorFermionsDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarVectorFermionsDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceIncludeVMD
    ("IncludeVMD",
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalarVectorFermionsDecayer::_includeVMD,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &PScalarVectorFermionsDecayer::_VMDid,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,double> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &PScalarVectorFermionsDecayer::_VMDmass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,double> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &PScalarVectorFermionsDecayer::_VMDwidth,
     0, 0, 0, -10000, 10000, false, false, true);

}

double PScalarVectorFermionsDecayer::me2(bool vertex, const int ichan,
					 const Particle & inpart,
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
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in PScalarVectorFermionsDecayer::me2()" 
				  << Exception::abortnow;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // identify the outgoing particles
  unsigned int ivec=0,iferm=0,ianti=0;int id;
  for(unsigned int ix=0;ix<decay.size();++ix)
    {
      id=decay[ix]->id();
      if(id==_outgoingV[imode()]){ivec=ix;}
      else if(id==_outgoingf[imode()]){iferm=ix;}
      else if(id==_outgoinga[imode()]){ianti=ix;}
    }
  // set up the spin info for the outgoing particles
  tcVectorSpinPtr vspin;
  tcFermionSpinPtr fspin,aspin;
  if(vertex)
    {
      // for the vector
      SpinPtr temp=new_ptr(VectorSpinInfo(inpart.momentum(),true));
      vspin= dynamic_ptr_cast<tcVectorSpinPtr>(temp);
      decay[ivec]->spinInfo(temp);
      // for the fermion
      temp=new_ptr(FermionSpinInfo(decay[iferm]->momentum(),true));
      decay[iferm]->spinInfo(temp);
      fspin = dynamic_ptr_cast<tcFermionSpinPtr>(temp);
      // for the antifermion
      temp=new_ptr(FermionSpinInfo(decay[ianti]->momentum(),true));
      decay[ianti]->spinInfo(temp);
      aspin =dynamic_ptr_cast<tcFermionSpinPtr>(temp);
    }
  // vectors containing the spinors
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  // calculate the spinor and antispinor
  SpinorBarWaveFunction fwave = SpinorBarWaveFunction(decay[iferm]->momentum(),
						      decay[iferm]->dataPtr(),outgoing);
  SpinorWaveFunction awave=SpinorWaveFunction(decay[ianti]->momentum(),
					      decay[ianti]->dataPtr(),outgoing);
  for(int ix=-1;ix<2;ix+=2)
    {
      // spinor for the fermion
      fwave.reset(ix);
      wavebar.push_back(fwave.Wave());
      if(vertex){fspin->setBasisState(ix,wavebar[(ix+1)/2].bar());}
      // spinorbar for the antifermion
      awave.reset(ix);
      wave.push_back(awave.Wave());
      if(vertex){aspin->setBasisState(ix,wave[(ix+1)/2]);}
    }    
  // calculate the wavefunctions for the outgoing vectors
  LorentzPolarizationVector vwave[3];
  VectorWaveFunction temp=VectorWaveFunction(decay[ivec]->momentum(),
					     decay[ivec]->dataPtr(),outgoing);
  for(int iy=-1;iy<2;++iy)
    {
      if(iy==0&&decay[ivec]->id()==ParticleID::gamma)
	{vwave[iy+1]=LorentzPolarizationVector();}
      else
	{temp.reset(iy);vwave[iy+1]=temp.Wave();}
      if(vertex){vspin->setBasisState(iy,vwave[iy+1]);}
    }
  // now compute the matrix element
  Complex ii(0.,1.);
  Lorentz5Momentum pff=decay[iferm]->momentum()+decay[ianti]->momentum();
  pff.rescaleMass();
  Energy2 mff2=pff.mass()*pff.mass();
  // compute the prefactor
  Complex pre=_coupling[imode()];
  pre /= mff2;
  // the VMD factor
  if(_includeVMD[imode()]>0)
    {
      Energy2 mrho2=_VMDmass[imode()]*_VMDmass[imode()];
      Energy2 mwrho=_VMDmass[imode()]*_VMDwidth[imode()];
      pre*= (-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
    }
  LorentzPolarizationVector eps,fcurrent;
  vector<int> ispin(3);ispin[ivec]=3;ispin[iferm]=2;ispin[ianti]=2;
  // compute the matrix element
  DecayMatrixElement newME(1,ispin);
  ispin.resize(4);ispin[0]=0;
  Complex s1s4,s2s3,s3s2,s4s1,s1s3,s2s4,s3s1,s4s2;
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(unsigned int iy=0;iy<2;++iy)
	{
	  s1s4 = wavebar[iy].s1()*wave[ix].s4();
	  s2s3 = wavebar[iy].s2()*wave[ix].s3();
	  s3s2 = wavebar[iy].s3()*wave[ix].s2();
	  s4s1 = wavebar[iy].s4()*wave[ix].s1();
	  s1s3 = wavebar[iy].s1()*wave[ix].s3();
	  s2s4 = wavebar[iy].s2()*wave[ix].s4();
	  s3s1 = wavebar[iy].s3()*wave[ix].s1();
	  s4s2 = wavebar[iy].s4()*wave[ix].s2();
	  // calculate the current
	  if(defaultDRep==HaberDRep)
	    {
	      fcurrent[0] =       s1s4+s2s3-s3s2-s4s1;
	      fcurrent[1] =  -ii*(s1s4-s2s3-s3s2+s4s1);
	      fcurrent[2] =       s1s3-s2s4-s3s1+s4s2;
	      fcurrent[3] = 
		+wavebar[iy].s1()*wave[ix].s1()+wavebar[iy].s2()*wave[ix].s2()
		-wavebar[iy].s3()*wave[ix].s3()-wavebar[iy].s4()*wave[ix].s4();
	    }
	  else
	    {
	      fcurrent[0] =      s1s4+s2s3-s3s2-s4s1;
	      fcurrent[1] = -ii*(s1s4-s2s3-s3s2+s4s1);
	      fcurrent[2] =      s1s3-s2s4-s3s1+s4s2;
	      fcurrent[3] =      s1s3+s2s4+s3s1+s4s2;
	    }
	  // compute the current for this part
	  eps = EpsFunction::product(decay[ivec]->momentum(),pff,fcurrent);
	  ispin[iferm+1]=2*iy-1;ispin[ianti+1]=2*ix-1;
	  for(unsigned int iz=0;iz<3;++iz)
	    {ispin[ivec+1]=iz-1;newME(ispin)=pre*(vwave[iz]*eps);}
	}	  
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  double me=newME.contract(rhoin).real();
  return me;
}

// method to return an object to calculate the 3 or higher body partial width
WidthCalculatorBasePtr 
PScalarVectorFermionsDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int imode=-1;
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit; int id2=(**pit).id();
  ++pit; int id3=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0)
      {
	if((id1==_outgoingV[ix]&&id2==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id1==_outgoingV[ix]&&id3==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id2==_outgoingV[ix]&&id1==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id2==_outgoingV[ix]&&id3==_outgoingf[ix]&&id1==_outgoinga[ix])||
	   (id3==_outgoingV[ix]&&id1==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id3==_outgoingV[ix]&&id2==_outgoingf[ix]&&id1==_outgoinga[ix]))
	  {imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // get the masses we need
  Energy m[3]={getParticleData(_outgoingV[imode])->mass(),
	       getParticleData(_outgoingf[imode])->mass(),
	       getParticleData(_outgoinga[imode])->mass()};
  return new_ptr(ThreeBodyAllOn1IntegralCalculator(3,-1000.,-0.9,
						   const_ptr_cast<tDecayIntegratorPtr>(this),
						   imode,m[0],m[1],m[2]));
}


double PScalarVectorFermionsDecayer::threeBodydGammads(int imodeb,Energy2 q2,
						       Energy2 mff2, Energy m1, 
						       Energy m2, Energy m3)
{
  // the masses of the external particles
  Energy q=sqrt(q2);
  Energy2 m12=m1*m1;
  Energy2 m22=m2*m2;
  Energy2 m32=m3*m3;
  // calculate the prefactor
  Complex pre=_coupling[imodeb],ii(0.,1.);
  pre /= mff2;
  // the VMD factor
  if(_includeVMD[imodeb]>0)
    {
      Energy2 mrho2=_VMDmass[imodeb]*_VMDmass[imodeb];
      Energy2 mwrho=_VMDmass[imodeb]*_VMDwidth[imodeb];
      pre*= (-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
    }
  double factor=real(pre*conj(pre));
  // compute the pieces from the integration limits
  Energy mff=sqrt(mff2);
  Energy e2star = 0.5*(mff2-m32+m22)/mff;
  Energy e1star = 0.5*(q2-mff2-m12)/mff;
  Energy e1sm = sqrt(e1star*e1star-m12);
  Energy e2sm = sqrt(e2star*e2star-m22);
  Energy2 a = 2*e1star*e2star+m12+m22;
  Energy2 b = 2*e1sm*e2sm;
  // term independent of s3
  double me = 2*b*(2*(m12*(mff2*mff2 + 4*mff2*m2*m3 -(m22 - m32)*(m22 - m32)) + 
		      2*m2*(m12 +m22)*m3*(-mff2 +m22 + q2))
		   +(m12 +m22)*(m12 +m22)*(-mff2 +m22 - 2*m2*m3 - m32)
		   -(mff2 +m22 + 2*m2*m3 - m32)*(-mff2 +m22 + q2)*(-mff2 +m22 + q2));
  // linear term
  me+= 2.*a*b*(2*(-mff2*mff2 - (m22 - m32)*(m12 - q2) + 
		  mff2*(m12 + m22 + m32 + q2)));
  // quadratic term
  me+=-4.*mff2*b*(3.*a*a+b*b)/3.;
  me*=-factor;
  // phase space factors
  return me/256./pi/pi/pi/q2/q;
}
}


