// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonPScalarFermionsDecayer class.
//
//  Author: Peter Richardson
//

#include "VectorMesonPScalarFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonPScalarFermionsDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::EpsFunction;
using ThePEG::Helicity::FermionSpinPtr;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::defaultDRep;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;


VectorMesonPScalarFermionsDecayer::~VectorMesonPScalarFermionsDecayer() {}

bool VectorMesonPScalarFermionsDecayer::accept(const DecayMode & dm) const {
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
      cout << "testing " << ix << "  " << _incoming.size() << "   " << id0 
	   << "    " << _incoming[ix] << endl;
      if(_incoming[ix]==id0)
      {
	cout << "matched incoming id" << id1 << " " << id2 << " " << id3 << endl;
	cout << _incoming[ix] << " -> " << _outgoingP[ix] << "   " 
	     << _outgoingf[ix] << "    " << _outgoinga[ix] << "   " << endl;
	if((id1==_outgoingP[ix]&&id2==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id1==_outgoingP[ix]&&id3==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id2==_outgoingP[ix]&&id1==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id2==_outgoingP[ix]&&id3==_outgoingf[ix]&&id1==_outgoinga[ix])||
	   (id3==_outgoingP[ix]&&id1==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id3==_outgoingP[ix]&&id2==_outgoingf[ix]&&id1==_outgoinga[ix]))
	  {allowed=true;}
      }
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector VectorMesonPScalarFermionsDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  // workout which mode we are doing
  int imode=-1;
  int id0=parent.id(),id1=children[0]->id(),id2=children[1]->id(),id3=children[2]->id();
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0)
      {
	if((id1==_outgoingP[ix]&&id2==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id1==_outgoingP[ix]&&id3==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id2==_outgoingP[ix]&&id1==_outgoingf[ix]&&id3==_outgoinga[ix])||
	   (id2==_outgoingP[ix]&&id3==_outgoingf[ix]&&id1==_outgoinga[ix])||
	   (id3==_outgoingP[ix]&&id1==_outgoingf[ix]&&id2==_outgoinga[ix])||
	   (id3==_outgoingP[ix]&&id2==_outgoingf[ix]&&id1==_outgoinga[ix]))
	  {imode=ix;}
	}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  generate(false,imode,parent,children);
  return children;
}


void VectorMesonPScalarFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingP << _outgoingf << _outgoinga << _maxweight
     << _includeVMD << _VMDid << _VMDmass << _VMDwidth;
}

void VectorMesonPScalarFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingP >> _outgoingf >> _outgoinga >> _maxweight
     >> _includeVMD >> _VMDid >> _VMDmass >> _VMDwidth;
}

ClassDescription<VectorMesonPScalarFermionsDecayer> VectorMesonPScalarFermionsDecayer::initVectorMesonPScalarFermionsDecayer;
// Definition of the static class description member.

void VectorMesonPScalarFermionsDecayer::Init() {

  static ClassDocumentation<VectorMesonPScalarFermionsDecayer> documentation
    ("The \\classname{VectorMesonPScalarFermionsDecayer} class is designed to "
     "perform the decay of a vector meson to a pseudoscalar meson and a "
     "fermion-antifermion pair.");





  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonPScalarFermionsDecayer::_incoming,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingP
    ("FirstPseudoScalar",
     "The PDG code for the outgoing pseudoscalar",
     &VectorMesonPScalarFermionsDecayer::_outgoingP,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingF
    ("SecondOutgoing",
     "The PDG code for the outgoing fermion",
     &VectorMesonPScalarFermionsDecayer::_outgoingf,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingA
    ("SecondOutgoing",
     "The PDG code for the outgoing antifermion",
     &VectorMesonPScalarFermionsDecayer::_outgoinga,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonPScalarFermionsDecayer::_coupling,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonPScalarFermionsDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceIncludeVMD
    ("IncludeVMD",
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &VectorMesonPScalarFermionsDecayer::_includeVMD,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &VectorMesonPScalarFermionsDecayer::_VMDid,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &VectorMesonPScalarFermionsDecayer::_VMDmass,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &VectorMesonPScalarFermionsDecayer::_VMDwidth,
     0, 0, 0, -10000, 10000, false, false, true);

}


// the hadronic currents    
vector<LorentzPolarizationVector>  
VectorMesonPScalarFermionsDecayer::decayCurrent(const bool vertex, const int imode,
						const int, const Particle & inpart,
						const ParticleVector & decay) const
  {
    // prefactor
    Complex pre=_coupling[imode];
    // storage for the current
    vector<LorentzPolarizationVector> temp;
    // find the outgoing particles
    unsigned int iferm=0,ianti=0,isca=0;int id;
    for(unsigned int ix=0;ix<decay.size();++ix)
      {
	id = decay[ix]->id();
	if(id==_outgoingP[ix]){isca=ix;}
	else if(id>0){iferm=ix;}
	else if(id<0){ianti=ix;}
      }    
    // construct the spin information objects for the  decay products
    FermionSpinPtr fspin,aspin;
    if(vertex)
      {
	decay[isca]->spinInfo(new_ptr(ScalarSpinInfo(decay[isca]->momentum(),true)));
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
    LorentzPolarizationVector fcurrent;
    Lorentz5Momentum pff=decay[iferm]->momentum()+decay[ianti]->momentum();
    pff.rescaleMass();
    // the factor for the off-shell photon
    Energy2 mff2=pff.mass()*pff.mass();
    pre /= mff2;
    // the VMD factor
    if(_includeVMD[imode]>0)
      {
	Energy2 mrho2=_VMDmass[imode]*_VMDmass[imode];
	Energy2 mwrho=_VMDmass[imode]*_VMDwidth[imode];
	pre*= (-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
      }
    unsigned int iloc;
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
	    // location in the vector
	    if(iferm<ianti){iloc=2*ix+iy;}
	    else{iloc=2*iy+ix;}
	    // add it to the vector
	    temp[iloc]=pre*EpsFunction::product(inpart.momentum(),pff,fcurrent);
	  }
      }
    // return the answer
    return temp;
  }
}
