// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2FermionDecayer class.
//

#include "VectorMeson2FermionDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson2FermionDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::FermionSpinPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::defaultDRep;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

VectorMeson2FermionDecayer::~VectorMeson2FermionDecayer() {}

bool VectorMeson2FermionDecayer::accept(const DecayMode & dm) const 
{
  // is this mode allowed
  bool allowed=false;
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0&&
	 ((_outgoing1[ix]==id1&&_outgoing2[ix]==id2)||
	  (_outgoing1[ix]==id2&&_outgoing2[ix]==id1))){allowed=true;}
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector VectorMeson2FermionDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  // workout which mode we are doing
  int imode=-1;
  int id=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id)
	{
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  bool cc=false;
  return generate(false,cc,imode,parent);
}


void VectorMeson2FermionDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void VectorMeson2FermionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<VectorMeson2FermionDecayer> VectorMeson2FermionDecayer::initVectorMeson2FermionDecayer;
// Definition of the static class description member.

void VectorMeson2FermionDecayer::Init() {

  static ClassDocumentation<VectorMeson2FermionDecayer> documentation
    ("The \\classname{VectorMeson2FermionDecayer} class is designed for the decay "
     "of vectro mesons to fermions. It is mainly used for the decay of vector mesons "
     "to electrons and muons.");

  static ParVector<VectorMeson2FermionDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson2FermionDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the outgoing fermion",
     &VectorMeson2FermionDecayer::_outgoing1,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing anti-fermion",
     &VectorMeson2FermionDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMeson2FermionDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMeson2FermionDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}


// the hadronic currents    
vector<LorentzPolarizationVector>  
VectorMeson2FermionDecayer::decayCurrent(const bool vertex, const int, 
					 const Particle & inpart,
					 const ParticleVector & decay) const
  {
    // prefactor
    double pre=_coupling[imode()]/inpart.mass();
    // storage for the current
    vector<LorentzPolarizationVector> temp;
    // find the fermion and the antifermion
    unsigned int iferm=0,ianti=0;
    for(unsigned int ix=0;ix<decay.size();++ix)
      {
	if(decay[ix]->id()>0){iferm=ix;}
	else if(decay[ix]->id()<0){ianti=ix;}
      }
    // construct the spin information objects for the  decay products
    FermionSpinPtr fspin,aspin;
    if(vertex)
      {
	SpinPtr sferm=new_ptr(FermionSpinInfo(decay[iferm]->momentum(),true));
	decay[iferm]->spinInfo(sferm);
	fspin = dynamic_ptr_cast<FermionSpinPtr>(sferm);
	SpinPtr santi=new_ptr(FermionSpinInfo(decay[ianti]->momentum(),true));
	decay[ianti]->spinInfo(santi);
	aspin =dynamic_ptr_cast<FermionSpinPtr>(santi);
      }
    // vectors containing the spinors
    vector<LorentzSpinor> wave;
    vector<LorentzSpinorBar> wavebar;
    // calculate the spinor and antispinor
    SpinorWaveFunction fwave = SpinorWaveFunction(decay[iferm]->momentum(),
						  decay[iferm]->dataPtr(),outgoing);
    SpinorBarWaveFunction awave=SpinorBarWaveFunction(decay[ianti]->momentum(),
						      decay[ianti]->dataPtr(),outgoing);
    for(int ix=-1;ix<2;ix+=2)
      {
	// spinor for the fermion
	fwave.reset(ix);
	wave.push_back(fwave.Wave());
	fspin->setBasisState(ix,wave[(ix+1)/2]);
	// spinorbar for the antifermion
	awave.reset(ix);
	wavebar.push_back(awave.Wave());
	aspin->setBasisState(ix,wavebar[(ix+1)/2].bar());
      }
    // now compute the currents
    Complex ii(0.,1.);
    temp.resize(4);
    unsigned int iloc;
    LorentzPolarizationVector vec;
    for(unsigned int ix=0;ix<2;++ix)
      {
	for(unsigned int iy=0;iy<2;++iy)
	  {
	    Complex s1s4 = wavebar[iy].s1()*wave[ix].s4();
	    Complex s2s3 = wavebar[iy].s2()*wave[ix].s3();
	    Complex s3s2 = wavebar[iy].s3()*wave[ix].s2();
	    Complex s4s1 = wavebar[iy].s4()*wave[ix].s1();
	    Complex s1s3 = wavebar[iy].s1()*wave[ix].s3();
	    Complex s2s4 = wavebar[iy].s2()*wave[ix].s4();
	    Complex s3s1 = wavebar[iy].s3()*wave[ix].s1();
	    Complex s4s2 = wavebar[iy].s4()*wave[ix].s2();
	    // calculate the current
	    if(defaultDRep==HaberDRep)
	      {
		vec[0] =       s1s4+s2s3-s3s2-s4s1;
		vec[1] =  -ii*(s1s4-s2s3-s3s2+s4s1);
		vec[2] =       s1s3-s2s4-s3s1+s4s2;
		vec[3] = 
		  +wavebar[iy].s1()*wave[ix].s1()+wavebar[iy].s2()*wave[ix].s2()
		  -wavebar[iy].s3()*wave[ix].s3()-wavebar[iy].s4()*wave[ix].s4();
	      }
	    else
	      {
		vec[0] =      s1s4+s2s3-s3s2-s4s1;
		vec[1] = -ii*(s1s4-s2s3-s3s2+s4s1);
		vec[2] =      s1s3-s2s4-s3s1+s4s2;
		vec[3] =      s1s3+s2s4+s3s1+s4s2;
	      }
	    // location in the vector
	    if(iferm<ianti){iloc=2*ix+iy;}
	    else{iloc=2*iy+ix;}
	    // add it to the vector
	    temp[iloc]=pre*vec;
	  }
      }
    // return the answer
    return temp;
  }

bool VectorMeson2FermionDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						double & coupling) const
{
  int imode=-1;
  int id=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id)
	{
	  if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  coupling=_coupling[imode];
  mecode=2;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
 }
}

