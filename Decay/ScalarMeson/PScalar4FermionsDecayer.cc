// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalar4FermionsDecayer class.
//

#include "PScalar4FermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalar4FermionsDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::EpsFunction;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::tcScalarSpinPtr;
using ThePEG::Helicity::tcFermionSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::defaultDRep;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

PScalar4FermionsDecayer::~PScalar4FermionsDecayer() {}

bool PScalar4FermionsDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  // must be four outgoing particles
  if(dm.products().size()!=4){return allowed;}
  // get the id's of the outgoing particles
  int id[4]; bool done[4]; unsigned int ix(0),iy(0);
  // ids of the particles
  int id0(dm.parent()->id()),idtemp(-1),idl1(-1),idl2(-1),idt[2];
  ParticleMSet::const_iterator pit = dm.products().begin();
  for ( ;pit!=dm.products().end();++pit)
    {id[ix]=(**pit).id();done[ix]=false;++ix;}
  // find the two lepton pairs
  // find the first fermion
  ix=0;
  do 
    {if(id[ix]>0&!done[ix]){done[ix]=true;idtemp=id[ix];}++ix;}
  while(ix<4&&idtemp<0);
  if(idtemp<0){return false;}
  // find its antiparticle
  ix=0;
  do 
    {if(id[ix]==-idtemp&!done[ix]){done[ix]=true;idl1=idtemp;}++ix;}
  while(ix<4&&idl1<0);
  if(idl1<0){return false;}
  idtemp=-1;
  for(ix=0;ix<4;++ix){if(!done[ix]){idt[iy]=id[ix];++iy;}}
  if(idt[0]==-idt[1]){idl2=abs(idt[0]);}
  if(idl2<0){return false;}
  // loop over the modes and see if this is one of them
  ix=0;
  do
    {
      if(_incoming[ix]==id0)
	{
	  if((idl1==_outgoing1[ix]&&idl2==_outgoing2[ix])||
	     (idl2==_outgoing1[ix]&&idl1==_outgoing2[ix])){allowed=true;}
	}
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector PScalar4FermionsDecayer::decay(const DecayMode & dm,
					      const Particle & parent) const {
  // workout which mode we are doing
  int imode=-1;
  int id0=parent.id(),idl[2],id; unsigned int iy(0);
  // find the ids of the outgoing fermions
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit)
    {id=(**pit).id();if(id>0){idl[iy]=id;++iy;}}
  iy=0;
  do
    {
      if(_incoming[iy]==id0)
	{
	  if((idl[0]==_outgoing1[iy]&&idl[1]==_outgoing2[iy])||
	     (idl[1]==_outgoing1[iy]&&idl[0]==_outgoing2[iy])){imode=iy;}
	}
      ++iy;
    }
  while(imode<0&&iy<_incoming.size());
  // perform the decay
  return generate(false,false,imode,parent);
}

void PScalar4FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoing1 << _outgoing2 << _maxweight 
     << _includeVMD << _VMDid << _VMDmass << _VMDwidth;
}

void PScalar4FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight 
     >> _includeVMD >> _VMDid >> _VMDmass >> _VMDwidth;
}

ClassDescription<PScalar4FermionsDecayer> PScalar4FermionsDecayer::initPScalar4FermionsDecayer;
// Definition of the static class description member.

void PScalar4FermionsDecayer::Init() {

  static ClassDocumentation<PScalar4FermionsDecayer> documentation
    ("The \\classname{PScalar4FermionsDecayer} class is designed for the decay"
     " of a pseudosclar meson to four fermions. It is intended for the decay of"
     "the pion to two electron-positron pairs.");

  static ParVector<PScalar4FermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalar4FermionsDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceOutcoming1
    ("Outgoing1",
     "The PDG code for the first outgoing fermion",
     &PScalar4FermionsDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceOutcoming2
    ("Outgoing2",
     "The PDG code for the first outgoing fermion",
     &PScalar4FermionsDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalar4FermionsDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalar4FermionsDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceIncludeVMD
    ("IncludeVMD",
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalar4FermionsDecayer::_includeVMD,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &PScalar4FermionsDecayer::_VMDid,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,double> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &PScalar4FermionsDecayer::_VMDmass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,double> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &PScalar4FermionsDecayer::_VMDwidth,
     0, 0, 0, -10000, 10000, false, false, true);

}

double PScalar4FermionsDecayer::me2(bool vertex, const int ichan,
				    const Particle & inpart,
				    const ParticleVector & decay) const
{
  bool identical=(_outgoing1[imode()]==_outgoing2[imode()]);
  // check if the decay particle has spin info 
  tcScalarSpinPtr inspin;
  if(inpart.spinInfo())
    {inspin = dynamic_ptr_cast<tcScalarSpinPtr>(inpart.spinInfo());}
  // if the spin info object exists use it
  if(inspin)
    {inspin->decayed(true);}
  else if(inpart.spinInfo())
    {throw DecayIntegratorError() << "Wrong type of spin info for the incoming particle"
				  << " in PScalar4FermionsDecayer::me2()" 
				  << Exception::abortnow;}
  else
    {
      SpinPtr newspin=new_ptr(ScalarSpinInfo(inpart.momentum(),true));
      inspin = dynamic_ptr_cast<tcScalarSpinPtr>(newspin);
      inspin->decayed(true);
      const_ptr_cast<tPPtr>(&inpart)->spinInfo(newspin);
    }
  // set up the spin info for the outgoing particles
  tcFermionSpinPtr fspin[2],aspin[2];
  if(vertex)
    {
      for(unsigned int ix=0;ix<2;++ix)
	{
	  // for the fermions
	  SpinPtr temp;
	  temp=new_ptr(FermionSpinInfo(decay[2*ix]->momentum(),true));
	  decay[2*ix]->spinInfo(temp);
	  fspin[ix] = dynamic_ptr_cast<tcFermionSpinPtr>(temp);
	  // for the antifermions
	  temp=new_ptr(FermionSpinInfo(decay[2*ix+1]->momentum(),true));
	  decay[2*ix+1]->spinInfo(temp);
	  aspin[ix] =dynamic_ptr_cast<tcFermionSpinPtr>(temp);
	}
    }
  // vectors containing the spinors
  vector<LorentzSpinor> wave[2];
  vector<LorentzSpinorBar> wavebar[2];
  // calculate the spinor and antispinor
  SpinorBarWaveFunction fwave[2] =
    {SpinorBarWaveFunction(decay[0]->momentum(),decay[0]->dataPtr(),outgoing),
     SpinorBarWaveFunction(decay[2]->momentum(),decay[2]->dataPtr(),outgoing)};
  SpinorWaveFunction awave[2] = 
    {SpinorWaveFunction(decay[1]->momentum(),decay[1]->dataPtr(),outgoing),
     SpinorWaveFunction(decay[3]->momentum(),decay[3]->dataPtr(),outgoing)};
  for(unsigned int ix=0;ix<2;++ix)
    {
      for(int iy=-1;iy<2;iy+=2)
	{
	  // spinor for the fermions
	  fwave[ix].reset(iy);
	  wavebar[ix].push_back(fwave[ix].Wave());
	  if(vertex){fspin[ix]->setBasisState(iy,wavebar[ix][(iy+1)/2].bar());}
	  // spinorbar for the antifermion
	  awave[ix].reset(iy);
	  wave[ix].push_back(awave[ix].Wave());
	  if(vertex){aspin[ix]->setBasisState(iy,wave[ix][(iy+1)/2]);}
	}
    }
  // momenta of the outgoing photons
  Lorentz5Momentum momentum[4];
  momentum[0]=decay[0]->momentum()+decay[1]->momentum();
  momentum[1]=decay[2]->momentum()+decay[3]->momentum();
  momentum[0].rescaleMass();
  momentum[1].rescaleMass();
  if(identical)
    {
      momentum[2]=decay[2]->momentum()+decay[1]->momentum();
      momentum[3]=decay[0]->momentum()+decay[3]->momentum();
      momentum[2].rescaleMass();
      momentum[3].rescaleMass();
    }
  // compute the currents for the two leptonic decays
  LorentzPolarizationVector current[4][2][2];
  Complex ii(0.,1.);
  Complex s1s4,s2s3,s3s2,s4s1,s1s3,s2s4,s3s1,s4s2;
  unsigned int it;
  for(unsigned int iz=0;iz<2;++iz)
    {
      if(iz==0){it=1;}
      else{it=0;}
      for(unsigned int ix=0;ix<2;++ix)
	{
	  for(unsigned int iy=0;iy<2;++iy)
	    {
	      // first two currents
	      s1s4 = wavebar[iz][iy].s1()*wave[iz][ix].s4();
	      s2s3 = wavebar[iz][iy].s2()*wave[iz][ix].s3();
	      s3s2 = wavebar[iz][iy].s3()*wave[iz][ix].s2();
	      s4s1 = wavebar[iz][iy].s4()*wave[iz][ix].s1();
	      s1s3 = wavebar[iz][iy].s1()*wave[iz][ix].s3();
	      s2s4 = wavebar[iz][iy].s2()*wave[iz][ix].s4();
	      s3s1 = wavebar[iz][iy].s3()*wave[iz][ix].s1();
	      s4s2 = wavebar[iz][iy].s4()*wave[iz][ix].s2();
	      // calculate the current
	      if(defaultDRep==HaberDRep)
		{
		  current[iz][iy][ix][0] =       s1s4+s2s3-s3s2-s4s1;
		  current[iz][iy][ix][1] =  -ii*(s1s4-s2s3-s3s2+s4s1);
		  current[iz][iy][ix][2] =       s1s3-s2s4-s3s1+s4s2;
		  current[iz][iy][ix][3] = 
		    +wavebar[iz][iy].s1()*wave[iz][ix].s1()
		    +wavebar[iz][iy].s2()*wave[iz][ix].s2()
		    -wavebar[iz][iy].s3()*wave[iz][ix].s3()
		    -wavebar[iz][iy].s4()*wave[iz][ix].s4();
		}
	      else
		{
		  current[iz][iy][ix][0] =      s1s4+s2s3-s3s2-s4s1;
		  current[iz][iy][ix][1] = -ii*(s1s4-s2s3-s3s2+s4s1);
		  current[iz][iy][ix][2] =      s1s3-s2s4-s3s1+s4s2;
		  current[iz][iy][ix][3] =      s1s3+s2s4+s3s1+s4s2;
		}
	      // the second two currents      
	      if(identical)
		{
		  s1s4 = wavebar[iz][iy].s1()*wave[it][ix].s4();
		  s2s3 = wavebar[iz][iy].s2()*wave[it][ix].s3();
		  s3s2 = wavebar[iz][iy].s3()*wave[it][ix].s2();
		  s4s1 = wavebar[iz][iy].s4()*wave[it][ix].s1();
		  s1s3 = wavebar[iz][iy].s1()*wave[it][ix].s3();
		  s2s4 = wavebar[iz][iy].s2()*wave[it][ix].s4();
		  s3s1 = wavebar[iz][iy].s3()*wave[it][ix].s1();
		  s4s2 = wavebar[iz][iy].s4()*wave[it][ix].s2();
		  // calculate the current
		  if(defaultDRep==HaberDRep)
		    {
		      current[iz+2][iy][ix][0] =       s1s4+s2s3-s3s2-s4s1;
		      current[iz+2][iy][ix][1] =  -ii*(s1s4-s2s3-s3s2+s4s1);
		      current[iz+2][iy][ix][2] =       s1s3-s2s4-s3s1+s4s2;
		      current[iz+2][iy][ix][3] = 
			+wavebar[it][iy].s1()*wave[iz][ix].s1()
			+wavebar[it][iy].s2()*wave[iz][ix].s2()
			-wavebar[it][iy].s3()*wave[iz][ix].s3()
			-wavebar[it][iy].s4()*wave[iz][ix].s4();
		    }
		  else
		    {
		      current[iz+2][iy][ix][0] =      s1s4+s2s3-s3s2-s4s1;
		      current[iz+2][iy][ix][1] = -ii*(s1s4-s2s3-s3s2+s4s1);
		      current[iz+2][iy][ix][2] =      s1s3-s2s4-s3s1+s4s2;
		      current[iz+2][iy][ix][3] =      s1s3+s2s4+s3s1+s4s2;
		    }
		}
	    }
	}
    }
  // invariants
  Energy m12(momentum[0].mass()*momentum[0].mass());
  Energy m34(momentum[1].mass()*momentum[1].mass()),m14(0.),m23(0.);
  Complex prop1(1./m12/m34),prop2(0.);
  if(identical)
    {
      m14=momentum[2].mass()*momentum[2].mass();
      m23=momentum[3].mass()*momentum[3].mass();
      prop2=1./m14/m23;
    }
  // the VMD factor if needed
  if(_includeVMD[imode()]>0)
    {
      Energy2 mrho2=_VMDmass[imode()]*_VMDmass[imode()];
      Energy2 mwrho=_VMDmass[imode()]*_VMDwidth[imode()];
      prop1*= 
	(-mrho2+ii*mwrho)/(m12-mrho2+ii*mwrho)*
	(-mrho2+ii*mwrho)/(m34-mrho2+ii*mwrho);
      if(identical)
	{prop2*= 
	    (-mrho2+ii*mwrho)/(m14-mrho2+ii*mwrho)*
	    (-mrho2+ii*mwrho)/(m23-mrho2+ii*mwrho);}
    }
  // prefactor
  Complex pre=_coupling[imode()]*4.*pi*SM().alphaEM()*inpart.mass(),diag;
  // now compute the matrix element
  LorentzPolarizationVector eps;
  vector<int> ispin(4,2);
  DecayMatrixElement newME(1,ispin);
  ispin.resize(5);ispin[0]=0;
  for(unsigned int if1=0;if1<2;++if1)
    {
      ispin[1]=2*if1-1;
      for(unsigned int if2=0;if2<2;++if2)
	{
	  ispin[2]=2*if2-1;
	  for(unsigned int ia1=0;ia1<2;++ia1)
	    {
	      ispin[3]=2*ia1-1;
	      for(unsigned int ia2=0;ia2<2;++ia2)
		{
		  ispin[4]=2*ia2-1;
		  // the first diagram
		  eps = EpsFunction::product(current[0][if1][ia1],momentum[1],
					     current[1][if2][ia2]);
		  diag = prop1*(eps*momentum[0]);
		  // exchanged diagram if identical particles
		  //  (sign due normal ordering) 
	          if(identical)
		    {
		      eps = EpsFunction::product(current[2][if1][ia2],momentum[3],
						 current[3][if2][ia1]);
		      diag-= prop2*(eps*momentum[2]);
		    }
		  newME(ispin)=pre*diag;
		}
	    }
	}
    }
  ME(newME);
  RhoDMatrix rhoin=RhoDMatrix(1);rhoin.average();
  double me=newME.contract(rhoin).real();
  if(identical){me*=0.25;}
  return me;
}
}
