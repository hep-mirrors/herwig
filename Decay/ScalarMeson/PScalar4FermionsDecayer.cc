// -*- C++ -*-
//
// PScalar4FermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalar4FermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      _maxweight[ix] = mode(ix)->maxWeight();
  }
}

PScalar4FermionsDecayer::PScalar4FermionsDecayer() 
  : _coupling(1,0.025159062/GeV), _incoming(1,111), _outgoing1(1,11), 
    _outgoing2(1,11), _maxweight(1,0.000234211), 
    _includeVMD(1,2),_VMDid(1,113), _VMDmass(1,0.7758*GeV), 
    _VMDwidth(1,0.1503*GeV), _initsize(1) {
  // intermediates
  generateIntermediates(false);
}

void PScalar4FermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoing1.size() || isize!=_outgoing2.size()||
     isize!=_maxweight.size() || isize!=_includeVMD.size()|| isize!=_VMDid.size()    ||
     isize!=_VMDmass.size()  || isize!=_VMDwidth.size())
    throw InitException() << "Inconsistent parameters in PScalar4FermionsDecayer"
			  << Exception::abortnow;
  // create the integration channels for each mode 
  tPDVector extpart(5);
  tPDPtr gamma=getParticleData(ParticleID::gamma);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  vector<double> wgt;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    wgt.resize(1);wgt[0]=1.;
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData( _outgoing1[ix]);
    extpart[2] = getParticleData(-_outgoing1[ix]);
    extpart[3] = getParticleData( _outgoing2[ix]);
    extpart[4] = getParticleData(-_outgoing2[ix]);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    // first channel always need this
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,-2);
    newchannel->addIntermediate(gamma     ,1,-1.1, 1,2);
    newchannel->addIntermediate(gamma     ,1,-1.1, 3,4);
    mode->addChannel(newchannel);
    if(_outgoing1[ix]==_outgoing2[ix]) {
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,-2);
      newchannel->addIntermediate(gamma     ,1,-1.1, 3,2);
      newchannel->addIntermediate(gamma     ,1,-1.1, 1,4);
      mode->addChannel(newchannel);
      wgt.resize(2);wgt[0]=0.5;wgt[1]=0.5;
    }
    else{wgt.resize(1);wgt[0]=1.;}
    addMode(mode,_maxweight[ix],wgt);
  }
  // set up the values for the VMD factor if needed (copy the default mass and width 
  //                                                 into the array)
  for(unsigned ix=0;ix<isize;++ix) {
    if(_includeVMD[ix]==1) {
      _VMDmass[ix]=getParticleData(_VMDid[ix])->mass();
      _VMDwidth[ix]=getParticleData(_VMDid[ix])->width();
    }
  }
}

int PScalar4FermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  // must be four outgoing particles
  if(children.size()!=4) return -1;
  // get the id's of the outgoing particles
  int id[4]={0,0,0,0}; 
  bool done[4]={false,false,false,false}; 
  unsigned int ix(0),iy(0);
  // ids of the particles
  int id0(parent->id()),idtemp(-1),idl1(-1),idl2(-1),idt[2];
  tPDVector::const_iterator pit = children.begin();
  for ( ;pit!=children.end();++pit) {
    id[ix]=(**pit).id();
    done[ix]=false;
    ++ix;
  }
  // find the two lepton pairs
  // find the first fermion
  ix=0;
  do {
    if( id[ix]>0 && !done[ix] ) {
      done[ix]=true;
      idtemp=id[ix];
    }
    ++ix;
  }
  while(ix<4&&idtemp<0);
  if(idtemp<0) return -1;
  // find its antiparticle
  ix=0;
  do {
    if( id[ix]==-idtemp && !done[ix] ) {
      done[ix]=true;
      idl1=idtemp;
    }
    ++ix;
  } while( ix<4 && idl1<0 );
  if(idl1<0) return -1;
  // find the second particle antiparticle pair
  for(ix=0;ix<4;++ix) {
    if(!done[ix]) {
      idt[iy]=id[ix];
      ++iy;
    }
  }
  if(idt[0]==-idt[1]) idl2=abs(idt[0]);
  if(idl2<0) return -1;
  // loop over the modes and see if this is one of them
  ix=0;
  int imode(-1);
  do {
    if(_incoming[ix]==id0) {
      if((idl1==_outgoing1[ix]&&idl2==_outgoing2[ix])||
	 (idl2==_outgoing1[ix]&&idl1==_outgoing2[ix])) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
}

void PScalar4FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/MeV) 
     << _incoming << _outgoing1 << _outgoing2 << _maxweight 
     << _includeVMD << _VMDid 
     << ounit(_VMDmass,MeV) << ounit(_VMDwidth,MeV);
}

void PScalar4FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/MeV) 
     >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight 
     >> _includeVMD >> _VMDid 
     >> iunit(_VMDmass,MeV) >> iunit(_VMDwidth,MeV);
}

ClassDescription<PScalar4FermionsDecayer> 
PScalar4FermionsDecayer::initPScalar4FermionsDecayer;
// Definition of the static class description member.

void PScalar4FermionsDecayer::Init() {

  static ClassDocumentation<PScalar4FermionsDecayer> documentation
    ("The PScalar4FermionsDecayer class is designed for the decay"
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
     "The PDG code for the second outgoing fermion",
     &PScalar4FermionsDecayer::_outgoing2,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalar4FermionsDecayer::_coupling,
     1/MeV, 0, ZERO, -10000/MeV, 10000/MeV, false, false, true);

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
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<PScalar4FermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &PScalar4FermionsDecayer::_VMDid,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<PScalar4FermionsDecayer,Energy> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &PScalar4FermionsDecayer::_VMDmass,
     1.*MeV, -1, ZERO, -10000*MeV, 10000*MeV, false, false, true);

  static ParVector<PScalar4FermionsDecayer,Energy> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &PScalar4FermionsDecayer::_VMDwidth,
     1.*MeV, -1, ZERO, -10000*MeV, 10000*MeV, false, false, true);

}

double PScalar4FermionsDecayer::me2(const int,
				    const Particle & inpart,
				    const ParticleVector & decay,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
					 PDT::Spin1Half,PDT::Spin1Half)));
  bool identical((_outgoing1[imode()]==_outgoing2[imode()]));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    // set up the spin information for the decay products
    for(unsigned int ix=0;ix<2;++ix) {
      SpinorBarWaveFunction::
	constructSpinInfo(_wavebar[ix],decay[2*ix  ],outgoing,true);
      SpinorWaveFunction::
	constructSpinInfo(_wave[ix]   ,decay[2*ix+1],outgoing,true);
    }
    return 0.;
  }
  // calculate the spinors
  for(unsigned int ix=0;ix<2;++ix) {
    SpinorBarWaveFunction::
      calculateWaveFunctions(_wavebar[ix],decay[2*ix  ],outgoing);
    SpinorWaveFunction::
      calculateWaveFunctions(_wave[ix]   ,decay[2*ix+1],outgoing);
  }
  // momenta of the outgoing photons
  Lorentz5Momentum momentum[4];
  momentum[0]=decay[0]->momentum()+decay[1]->momentum();momentum[0].rescaleMass();
  momentum[1]=decay[2]->momentum()+decay[3]->momentum();momentum[1].rescaleMass();
  if(identical) {
    momentum[2]=decay[2]->momentum()+decay[1]->momentum();momentum[2].rescaleMass();
    momentum[3]=decay[0]->momentum()+decay[3]->momentum();momentum[3].rescaleMass();
  }
  // compute the currents for the two leptonic decays
  LorentzPolarizationVectorE current[4][2][2];
  unsigned int it,ix,iy,iz;
  for(iz=0;iz<2;++iz) {
    it = iz==0 ? 1 : 0;
    for(ix=0;ix<2;++ix) {
      for(iy=0;iy<2;++iy) {
	current[iz][iy][ix] = _wave[iz][ix].vectorCurrent(_wavebar[iz][iy]);
	// the second two currents      
	if(identical)
	  current[iz+2][iy][ix] = _wave[it][ix].vectorCurrent(_wavebar[iz][iy]);
      }
    }
  }
  // invariants
  Energy2 m12(sqr(momentum[0].mass()));
  Energy2 m34(sqr(momentum[1].mass()));
  Energy2 m14(ZERO), m23(ZERO);
  complex<InvEnergy4> prop1(1./m12/m34),prop2(0./sqr(MeV2));
  Complex ii(0.,1.);
  if(identical) {
    m14=momentum[2].mass()*momentum[2].mass();
    m23=momentum[3].mass()*momentum[3].mass();
    prop2=1./m14/m23;
  }
  // the VMD factor if needed
  if(_includeVMD[imode()]>0) {
    Energy2 mrho2(_VMDmass[imode()]*_VMDmass[imode()]);
    Energy2 mwrho(_VMDmass[imode()]*_VMDwidth[imode()]);
    prop1 = prop1*(-mrho2+ii*mwrho)/(m12-mrho2+ii*mwrho)*
                  (-mrho2+ii*mwrho)/(m34-mrho2+ii*mwrho);
    if(identical) {
      prop2 = prop2*(-mrho2+ii*mwrho)/(m14-mrho2+ii*mwrho)*
	            (-mrho2+ii*mwrho)/(m23-mrho2+ii*mwrho);
    }
  }
  // prefactor
  Complex pre(_coupling[imode()]*4.*Constants::pi
	      *SM().alphaEM()*inpart.mass());
  Complex diag;
  // now compute the matrix element
  LorentzVector<complex<Energy3> > eps;
  vector<unsigned int> ispin(5,0);
  for(ispin[1]=0; ispin[1]<2;++ispin[1]) {
    for(ispin[2]=0;ispin[2]<2;++ispin[2]) {
      for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
	for(ispin[4]=0;ispin[4]<2;++ispin[4]) {
	  // the first diagram
	  eps = epsilon(current[0][ispin[1]][ispin[2]],momentum[1],
			current[1][ispin[3]][ispin[4]]);
	  diag = Complex(prop1*(eps*momentum[0]));
	  // exchanged diagram if identical particles
	  //  (sign due normal ordering) 
	  if(identical) {
	    eps = epsilon(current[2][ispin[1]][ispin[4]],momentum[3],
			  current[3][ispin[3]][ispin[2]]);
	    diag-= Complex(prop2*(eps*momentum[2]));
	  }
	  (*ME())(ispin)=pre*diag;
	}
      }
    }
  }
  double me=ME()->contract(_rho).real();
  if(identical) me *= 0.25;
  return me;
}

// output the setup info for the particle database
void PScalar4FermionsDecayer::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming   " << ix << " " 
	     << _incoming[ix]   << "\n";
      output << "newdef " << name() << ":Outgoing1  " << ix << " " 
	     << _outgoing1[ix]  << "\n";
      output << "newdef " << name() << ":Outgoing2  " << ix << " " 
	     << _outgoing2[ix]  << "\n";
      output << "newdef " << name() << ":Coupling   " << ix << " " 
	     << _coupling[ix]*MeV   << "\n";
      output << "newdef " << name() << ":MaxWeight  " << ix << " " 
	     << _maxweight[ix]  << "\n";
      output << "newdef " << name() << ":IncludeVMD " << ix << " " 
	     << _includeVMD[ix] << "\n";
      output << "newdef " << name() << ":VMDID      " << ix << " " 
	     << _VMDid[ix]      << "\n";
      output << "newdef " << name() << ":VMDmass    " << ix << " " 
	     << _VMDmass[ix]/MeV    << "\n";
      output << "newdef " << name() << ":VMDwidth   " << ix << " " 
	     << _VMDwidth[ix]/MeV   << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming   " << ix << " " 
	     << _incoming[ix]   << "\n";
      output << "insert " << name() << ":Outgoing1  " << ix << " " 
	     << _outgoing1[ix]  << "\n";
      output << "insert " << name() << ":Outgoing2  " << ix << " " 
	     << _outgoing2[ix]  << "\n";
      output << "insert " << name() << ":Coupling   " << ix << " " 
		 << _coupling[ix]*MeV   << "\n";
      output << "insert " << name() << ":MaxWeight  " << ix << " " 
	     << _maxweight[ix]  << "\n";
      output << "insert " << name() << ":IncludeVMD " << ix << " " 
	     << _includeVMD[ix] << "\n";
      output << "insert " << name() << ":VMDID      " << ix << " " 
	     << _VMDid[ix]      << "\n";
      output << "insert " << name() << ":VMDmass    " << ix << " " 
	     << _VMDmass[ix]/MeV    << "\n";
      output << "insert " << name() << ":VMDwidth   " << ix << " " 
	     << _VMDwidth[ix]/MeV   << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
