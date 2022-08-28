// -*- C++ -*-
//
// VectorMesonPScalarFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonPScalarFermionsDecayer class.
//
//  Author: Peter Richardson
//

#include "VectorMesonPScalarFermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonPScalarFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      maxweight_[ix] = mode(ix)->maxWeight();
      weight_[ix]    = mode(ix)->channels()[1].weight();
    }
  }
}

VectorMesonPScalarFermionsDecayer::VectorMesonPScalarFermionsDecayer() {
  // intermediates
  generateIntermediates(false);
}

void VectorMesonPScalarFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxweight_.size()|| isize!=includeVMD_.size()||
     isize!=VMDid_.size()     || isize!=VMDmass_.size()  || isize!=VMDwidth_.size()||
     isize!=weight_.size())
    throw InitException() << "Inconsistent parameters in VectorMesonPScalar"
			  << "FermionsDecayer" << Exception::abortnow;
  // create the integration channel for each mode
  tPDPtr gamma(getParticleData(ParticleID::gamma)),rho;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    rho=getParticleData(VMDid_[ix]);
    tPDPtr in     =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second),
		     getParticleData(-outgoing_[ix].second)};
    PhaseSpaceModePtr newmode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    // photon channel
    PhaseSpaceChannel newChannel ((PhaseSpaceChannel(newmode),0,gamma,0,1,1,2,1,3));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    newChannel.weight(1.-weight_[ix]);
    newmode->addChannel(newChannel);
    // vmd channel
    PhaseSpaceChannel newChannel2((PhaseSpaceChannel(newmode),0,rho,0,1,1,2,1,3));
    newChannel2.weight(weight_[ix]);
    newmode->addChannel(newChannel2);
    addMode(newmode);
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

int VectorMesonPScalarFermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
						  const tPDVector & children) const {
  int imode(-1);
  // must be three outgoing particles
  if(children.size()!=3){return imode;}
  // ids of the particles
  int id0(parent->id()),idf[2]={0,0},ids(0);
  unsigned int nf(0);
  tPDVector::const_iterator pit = children.begin();
  for( ;pit!=children.end();++pit) {
    if((**pit).iSpin()==PDT::Spin0) ids=(**pit).id();
    else {
      idf[nf]=(**pit).id();
      ++nf;
    }
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==ids) {
      if((idf[0]==outgoing_[ix].second&&idf[1]==-outgoing_[ix].second)||
	 (idf[1]==outgoing_[ix].second&&idf[0]==-outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  // perform the decay
  cc=false;
  return imode;
}

void VectorMesonPScalarFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/GeV) << incoming_ << outgoing_
     << maxweight_ << weight_ << includeVMD_ << VMDid_ << ounit(VMDmass_,GeV) 
     << ounit(VMDwidth_,GeV);
}

void VectorMesonPScalarFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/GeV) >> incoming_ >> outgoing_
     >> maxweight_ >> weight_ >> includeVMD_ >> VMDid_ >> iunit(VMDmass_,GeV) 
     >> iunit(VMDwidth_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonPScalarFermionsDecayer,DecayIntegrator>
describeHerwigVectorMesonPScalarFermionsDecayer("Herwig::VectorMesonPScalarFermionsDecayer", "HwVMDecay.so");

void VectorMesonPScalarFermionsDecayer::Init() {

  static ClassDocumentation<VectorMesonPScalarFermionsDecayer> documentation
    ("The VectorMesonPScalarFermionsDecayer class is designed to "
     "perform the decay of a vector meson to a pseudoscalar meson and a "
     "fermion-antifermion pair.");

  static Command<VectorMesonPScalarFermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, fermion, coupling(1/GeV), includeVMD,"
     " VMD id, VMD mass(GeV), VMD width(GeV) and max weight and weight of second channel for a decay."
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &VectorMesonPScalarFermionsDecayer::setUpDecayMode, false);

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceOutcomingP
    ("OutgoingPseudoScalar","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceOutcomingF
    ("OutgoingFermion","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceOutcomingA
    ("OutgoingAntiFermion","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceWeight
    ("Weight","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceIncludeVMD
    ("IncludeVMD","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceVMDID
    ("VMDID","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceVMDmass
    ("VMDmass","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceVMDwidth
    ("VMDwidth","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

}

void VectorMesonPScalarFermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[1],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[2],outgoing,true);
}

double VectorMesonPScalarFermionsDecayer::me2(const int, const Particle & part,
					      const tPDVector & ,
					      const vector<Lorentz5Momentum> & momenta,
					      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,
					 PDT::Spin1Half,PDT::Spin1Half)));
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[1],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[2],ix,Helicity::outgoing);
  }
  // the factor for the off-shell photon
  Lorentz5Momentum pff(momenta[1]+momenta[2]);
  pff.rescaleMass();
  Energy2 mff2(pff.mass2());
  // prefactor
  complex<InvEnergy3> pre(coupling_[imode()]/mff2);
  Complex ii(0.,1.);
  // the VMD factor
  if(includeVMD_[imode()]>0) {
    Energy2 mrho2(VMDmass_[imode()]*VMDmass_[imode()]);
    Energy2 mwrho(VMDmass_[imode()]*VMDwidth_[imode()]);
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  // calculate the matrix element
  LorentzPolarizationVector temp;
  unsigned int ix,iy,iz;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      temp=pre*epsilon(part.momentum(),pff,
				    wave_[ix].vectorCurrent(wavebar_[iy]));
      for(iz=0;iz<3;++iz) 
	(*ME())(iz,0,iy,ix)=temp.dot(vectors_[iz]); 
    }
  }
  // code for the spin averaged me for testing only
  // Energy  m[4]={part.mass(),momenta[0].mass(),
  // 		momenta[1].mass(),momenta[2].mass()};
  // Energy2 m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
  // Lorentz5Momentum p12=momenta[0]+momenta[1];p12.rescaleMass();
  // Energy2 m122(p12.mass2());
  // cout << "testing the matrix element " 
  //      << -1./3.*(pre*conj(pre)).real()*(-2*m122*m122*mff2 - mff2*mff2*mff2 + 
  // 		   m2[1]*(2*m2[2]*m2[3] - 2*m2[3]*m2[3] + 
  // 			  m2[1]*(m2[2] - 2*m[2]*m[3] - m2[3])) - 
  // 		   2*m[2]*(m2[2]*m[2] - 2*m2[1]*m[3] - m[2]*m2[3])*
  // 		   m2[0] - (m2[2] + 2*m[2]*m[3] - m2[3])*
  // 		   m2[0]*m2[0] + mff2*mff2*
  // 		   (2*m2[1] + (m[2] - m[3])*(m[2] - m[3]) + 2*m2[0]) - 
  // 		   mff2*(m2[1]*m2[1] + 2*m2[1]*m[2]*(m[2] - 2*m[3]) + 
  // 			 2*m2[2]*m2[3] - 2*(2*m[2] - m[3])*m[3]*m2[0] + 
  // 			 m2[0]*m2[0]) + 2*m122*
  // 		   (-mff2*mff2 - (m[2] - m[3])*(m[2] + m[3])*(m[1] - m[0])*
  // 		    (m[1] + m[0]) + mff2*
  // 		    (m2[1] + m2[2] + m2[3] + m2[0])))
  //      << endl;
  // return the answer
  return ME()->contract(rho_).real();
}

WidthCalculatorBasePtr 
VectorMesonPScalarFermionsDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int imode(-1);
  // ids of the particles
  int id0(dm.parent()->id());
  int idf[2] = {0,0};
  int ids(0);
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit) {
    if((**pit).iSpin()==PDT::Spin0){ids=(**pit).id();}
    else{idf[nf]=(**pit).id();++nf;}
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==ids) {
      if((idf[0]==outgoing_[ix].second&&idf[1]==-outgoing_[ix].second)||
	 (idf[1]==outgoing_[ix].second&&idf[0]==-outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  // get the masses we need
  Energy m[3]={getParticleData( outgoing_[imode].first )->mass(),
	       getParticleData( outgoing_[imode].second)->mass(),
	       getParticleData(-outgoing_[imode].second)->mass()};
  return 
    new_ptr(ThreeBodyAllOn1IntegralCalculator<VectorMesonPScalarFermionsDecayer>
	    (3,-1000.*MeV,-0.8*MeV,-0.8,*this,imode,m[0],m[1],m[2]));
} 

InvEnergy VectorMesonPScalarFermionsDecayer::
threeBodydGammads(int imodeb, const Energy2 q2, const Energy2 mff2, 
		  const  Energy m1, const Energy m2, const  Energy m3) const {
  // the masses of the external particles
  Energy q(sqrt(q2));
  Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
  // prefactor
  complex<InvEnergy3> pre(coupling_[imodeb]/mff2);
  Complex ii(0.,1.);
  // the VMD factor
  if(includeVMD_[imodeb]>0) {
    Energy2 mrho2(VMDmass_[imodeb]*VMDmass_[imodeb]),
      mwrho(VMDmass_[imodeb]*VMDwidth_[imodeb]);
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  InvEnergy6 factor(real(pre*conj(pre)));
  // compute the pieces from the integration limits
  Energy mff(sqrt(mff2)),e2star(0.5*(mff2-m32+m22)/mff),e1star(0.5*(q2-mff2-m12)/mff),
    e1sm(sqrt(e1star*e1star-m12)),e2sm(sqrt(e2star*e2star-m22));
  Energy2 a(2*e1star*e2star+m12+m22),b(2*e1sm*e2sm);
  // term independent of s3
  Energy8 me = 2*b*(-mff2*mff2*mff2 +m12*m12*(m22 - 2*m2*m3 - m32) - 
		   2*m22*(m22 - m32)*q2 -(m22 + 2*m2*m3 - m32)*q2*q2 + 
		   mff2*mff2*(2*m12 +(m2-m3)*(m2-m3)+2*q2) + 2*m12*m3*
		   ((m22-m32)*m3 + 2*m2*q2) - 
		   mff2*(m12*m12 + 2*m12*m2*(m2 - 2*m3) + 2*m22*m32 - 
			 2*(2*m2 - m3)*m3*q2 + q2*q2));
  // linear term
  me+= 2.*a*b*(2*(-mff2*mff2-(m22-m32)*(m12-q2)+mff2*(m12+m22+m32+q2)));
  // quadratic term
  me+=2*b*(3.*a*a+b*b)/3.*(-2*mff2);
  // phase space factors
  using Constants::pi;
  return -factor * me/768./pi/pi/pi/q2/q;
}

// output the setup information for the particle database
void VectorMesonPScalarFermionsDecayer::dataBaseOutput(ofstream & output,
						       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second << "  " 
	   << coupling_[ix]*GeV << " " << includeVMD_[ix] << " "
	   << VMDid_[ix] << " " << VMDmass_[ix]/GeV    << " "
	   << VMDwidth_[ix]/GeV << " " << maxweight_[ix]  << " "
	   << weight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}


string VectorMesonPScalarFermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 0";
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
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt2 = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  includeVMD_.push_back(VMDinclude);
  VMDid_.push_back(VMDid);
  VMDmass_.push_back(VMDmass*GeV);
  VMDwidth_.push_back(VMDwidth*GeV);
  maxweight_.push_back(wgt);
  weight_.push_back(wgt2);
  // success
  return "";
}
