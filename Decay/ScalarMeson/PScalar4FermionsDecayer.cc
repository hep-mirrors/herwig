// -*- C++ -*-
//
// PScalar4FermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalar4FermionsDecayer class.
//

#include "PScalar4FermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
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
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      maxweight_[ix] = mode(ix)->maxWeight();
  }
}

PScalar4FermionsDecayer::PScalar4FermionsDecayer() {
  // intermediates
  generateIntermediates(false);
}

void PScalar4FermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!=maxweight_.size() || isize!=includeVMD_.size()|| isize!=VMDid_.size()    ||
     isize!=VMDmass_.size()  || isize!=VMDwidth_.size())
    throw InitException() << "Inconsistent parameters in PScalar4FermionsDecayer"
  			  << Exception::abortnow;
  // create the integration channels for each mode
  tPDPtr gamma=getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in = getParticleData(incoming_[ix]);
    tPDVector out={getParticleData( outgoing_[ix].first),
		   getParticleData(-outgoing_[ix].first),
		   getParticleData( outgoing_[ix].second),
		   getParticleData(-outgoing_[ix].second)};
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,gamma,0,gamma,1,1,1,2,2,3,2,4));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    newChannel.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel);
    PhaseSpaceChannel newChannel2((PhaseSpaceChannel(mode),0,gamma,0,gamma,1,3,1,2,2,1,2,4));
    newChannel2.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    newChannel2.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel2);
    addMode(mode);
  }
  // set up the values for the VMD factor if needed (copy the default mass and width 
  //                                                 into the array)
  for(unsigned ix=0;ix<isize;++ix) {
    if(includeVMD_[ix]==1) {
      VMDmass_[ix]=getParticleData(VMDid_[ix])->mass();
      VMDwidth_[ix]=getParticleData(VMDid_[ix])->width();
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
    if(incoming_[ix]==id0) {
      if((idl1==outgoing_[ix].first&&idl2==outgoing_[ix].second)||
	 (idl2==outgoing_[ix].first&&idl1==outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void PScalar4FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/MeV) 
     << incoming_ << outgoing_ << maxweight_ 
     << includeVMD_ << VMDid_ 
     << ounit(VMDmass_,MeV) << ounit(VMDwidth_,MeV);
}

void PScalar4FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/MeV) 
     >> incoming_ >> outgoing_ >> maxweight_ 
     >> includeVMD_ >> VMDid_ 
     >> iunit(VMDmass_,MeV) >> iunit(VMDwidth_,MeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalar4FermionsDecayer,DecayIntegrator>
describeHerwigPScalar4FermionsDecayer("Herwig::PScalar4FermionsDecayer", "HwSMDecay.so");

void PScalar4FermionsDecayer::Init() {

  static ClassDocumentation<PScalar4FermionsDecayer> documentation
    ("The PScalar4FermionsDecayer class is designed for the decay"
     " of a pseudosclar meson to four fermions. It is intended for the decay of"
     "the pion to two electron-positron pairs.");

  static Command<PScalar4FermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, 1st fermion, 2nd fermion, coupling(1/GeV), includeVMD,"
     " VMD id, VMD mass(GeV), VMD width(GeV) and max weight for a decay."
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalar4FermionsDecayer::setUpDecayMode, false);

  static Deleted<PScalar4FermionsDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceOutcoming1
    ("Outgoing1","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceOutcoming2
    ("Outgoing2","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceIncludeVMD
    ("IncludeVMD","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceVMDID
    ("VMDID","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceVMDmass
    ("VMDmass","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceVMDwidth
    ("VMDwidth","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

}
void PScalar4FermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<2;++ix) {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_[ix],decay[2*ix  ],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave_[ix]   ,decay[2*ix+1],outgoing,true);
  }
}

double PScalar4FermionsDecayer::me2(const int,const Particle & part,
				    const tPDVector &,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
  					 PDT::Spin1Half,PDT::Spin1Half)));
  bool identical((outgoing_[imode()].first==outgoing_[imode()].second));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the spinors
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix].resize(2);
    wave_   [ix].resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel) {
      wavebar_[ix][ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[2*ix  ],ihel,Helicity::outgoing);
      wave_   [ix][ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[2*ix+1],ihel,Helicity::outgoing);
    }
  }
  // // momenta of the outgoing photons
  Lorentz5Momentum poff[4];
  poff[0]=momenta[0]+momenta[1];
  poff[0].rescaleMass();
  poff[1]=momenta[2]+momenta[3];
  poff[1].rescaleMass();
  if(identical) {
    poff[2]=momenta[2]+momenta[1];
    poff[2].rescaleMass();
    poff[3]=momenta[0]+momenta[3];
    poff[3].rescaleMass();
  }
  // compute the currents for the two leptonic decays
  LorentzPolarizationVectorE current[4][2][2];
  unsigned int it,ix,iy,iz;
  for(iz=0;iz<2;++iz) {
    it = iz==0 ? 1 : 0;
    for(ix=0;ix<2;++ix) {
      for(iy=0;iy<2;++iy) {
  	current[iz][iy][ix] = wave_[iz][ix].vectorCurrent(wavebar_[iz][iy]);
  	// the second two currents      
  	if(identical)
  	  current[iz+2][iy][ix] = wave_[it][ix].vectorCurrent(wavebar_[iz][iy]);
      }
    }
  }
  // invariants
  Energy2 m12(sqr(poff[0].mass()));
  Energy2 m34(sqr(poff[1].mass()));
  Energy2 m14(ZERO), m23(ZERO);
  complex<InvEnergy4> prop1(1./m12/m34),prop2(0./sqr(MeV2));
  Complex ii(0.,1.);
  if(identical) {
    m14=poff[2].mass()*poff[2].mass();
    m23=poff[3].mass()*poff[3].mass();
    prop2=1./m14/m23;
  }
  // the VMD factor if needed
  if(includeVMD_[imode()]>0) {
    Energy2 mrho2(VMDmass_[imode()]*VMDmass_[imode()]);
    Energy2 mwrho(VMDmass_[imode()]*VMDwidth_[imode()]);
    prop1 = prop1*(-mrho2+ii*mwrho)/(m12-mrho2+ii*mwrho)*
                  (-mrho2+ii*mwrho)/(m34-mrho2+ii*mwrho);
    if(identical) {
      prop2 = prop2*(-mrho2+ii*mwrho)/(m14-mrho2+ii*mwrho)*
  	            (-mrho2+ii*mwrho)/(m23-mrho2+ii*mwrho);
    }
  }
  // prefactor
  Complex pre(coupling_[imode()]*4.*Constants::pi
  	      *SM().alphaEM()*part.mass());
  Complex diag;
  // now compute the matrix element
  LorentzVector<complex<Energy3> > eps;
  vector<unsigned int> ispin(5,0);
  for(ispin[1]=0; ispin[1]<2;++ispin[1]) {
    for(ispin[2]=0;ispin[2]<2;++ispin[2]) {
      for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
  	for(ispin[4]=0;ispin[4]<2;++ispin[4]) {
  	  // the first diagram
  	  eps = epsilon(current[0][ispin[1]][ispin[2]],poff[1],
  			current[1][ispin[3]][ispin[4]]);
  	  diag = Complex(prop1*(eps*poff[0]));
  	  // exchanged diagram if identical particles
  	  //  (sign due normal ordering) 
  	  if(identical) {
  	    eps = epsilon(current[2][ispin[1]][ispin[4]],poff[3],
  			  current[3][ispin[3]][ispin[2]]);
  	    diag-= Complex(prop2*(eps*poff[2]));
  	  }
  	  (*ME())(ispin)=pre*diag;
  	}
      }
    }
  }
  double me=ME()->contract(rho_).real();
  if(identical) me *= 0.25;
  return me;
}

// output the setup info for the particle database
void PScalar4FermionsDecayer::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << includeVMD_[ix] << " "
	   << VMDid_[ix] << " " << VMDmass_[ix]/GeV << " "
	   << VMDwidth_[ix]/GeV << " " << maxweight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string PScalar4FermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // vmd option  
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDinclude = stoi(stype);
  // vmd id
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDid = stoi(stype);
  // vmd mass
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double VMDmass = stof(stype);
  // vmd width
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double VMDwidth = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  includeVMD_.push_back(VMDinclude);
  VMDid_.push_back(VMDid);
  VMDmass_.push_back(VMDmass*GeV);
  VMDwidth_.push_back(VMDwidth*GeV);
  maxweight_.push_back(wgt);
  // success
  return "";
}
