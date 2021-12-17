// -*- C++ -*-
//
// PScalarVectorFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorFermionsDecayer class.
//

#include "PScalarVectorFermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarVectorFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      maxweight_[ix]=mode(ix)->maxWeight();
  }
}

PScalarVectorFermionsDecayer::PScalarVectorFermionsDecayer() {
  // intermediates
  generateIntermediates(false);
}

void PScalarVectorFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!=maxweight_.size() || isize!=includeVMD_.size()||
     isize!=VMDid_.size()     || isize!=VMDmass_.size()   || isize!=VMDwidth_.size())
    throw InitException() << "Inconsistent parameters in PScalarVectorFermionsDecayer"
  			  << Exception::abortnow;
  // create the integration channel for each mode 
  tPDPtr gamma(getParticleData(ParticleID::gamma));
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second),
		     getParticleData(-outgoing_[ix].second)};
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,gamma,0,1,1,2,1,3));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel);
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

int PScalarVectorFermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  int imode(-1);
  // must be three outgoing particles
  if(children.size()!=3) return imode;
  // ids of the particles
  int id0(parent->id()),idf[2]={0,0},idv(0);
  unsigned int nf(0);
  tPDVector::const_iterator pit = children.begin();
  for( ;pit!=children.end();++pit) {
    if((**pit).iSpin()==PDT::Spin1) {
      idv=(**pit).id();
    }
    else {
      idf[nf]=(**pit).id();
      ++nf;
    }
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==idv)
      {if((idf[0]==outgoing_[ix].second&&idf[1]==-outgoing_[ix].second)||
	  (idf[1]==outgoing_[ix].second&&idf[0]==-outgoing_[ix].second)){imode=ix;}}
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void PScalarVectorFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/MeV) 
     << incoming_ << outgoing_ << maxweight_
     << includeVMD_ << VMDid_ 
     << ounit(VMDmass_,MeV) << ounit(VMDwidth_,MeV);
}

void PScalarVectorFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/MeV) 
     >> incoming_ >> outgoing_ >> maxweight_
     >> includeVMD_ >> VMDid_ 
     >> iunit(VMDmass_,MeV) >> iunit(VMDwidth_,MeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarVectorFermionsDecayer,DecayIntegrator>
describeHerwigPScalarVectorFermionsDecayer("Herwig::PScalarVectorFermionsDecayer", "HwSMDecay.so");

void PScalarVectorFermionsDecayer::Init() {

  static ClassDocumentation<PScalarVectorFermionsDecayer> documentation
    ("The PScalarVectorFermionsDecayer class is designed"
     " for the decay of a pseudoscalar meson to a photon and a"
     "fermion-antifermion pair");

  static Command<PScalarVectorFermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, fermion, coupling(1/GeV), includeVMD,"
     " VMD id, VMD mass(GeV), VMD width(GeV) and max weight for a decay."
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalarVectorFermionsDecayer::setUpDecayMode, false);
  
  static Deleted<PScalarVectorFermionsDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceOutcomingV
    ("OutgoingVector","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceOutcomingF
    ("OutgoingFermion","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceOutcomingA
    ("OutgoingAntiFermion","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceIncludeVMD
    ("IncludeVMD","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceVMDID
    ("VMDID","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceVMDmass
    ("VMDmass","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceVMDwidth
    ("VMDwidth","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

}

void PScalarVectorFermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  // set up the spin information for the decay products
  VectorWaveFunction::
    constructSpinInfo(vectors_,decay[0],outgoing,true,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[1],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[2],outgoing,true);
}

double PScalarVectorFermionsDecayer::me2(const int,const Particle & part,
					 const tPDVector &,
					 const vector<Lorentz5Momentum> & momenta,
					 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1Half,
					 PDT::Spin1Half)));
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate polarization vector
  vectors_.resize(3);
  for(unsigned int ihel=0;ihel<3;ihel+=2) {
    vectors_[ihel] = HelicityFunctions::polarizationVector(-momenta[0],ihel,Helicity::outgoing);
  }
  // calculate the spinors
  wavebar_.resize(2);
  wave_   .resize(2);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wavebar_[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[1],ihel,Helicity::outgoing);
    wave_   [ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[2],ihel,Helicity::outgoing);
  }
  // now compute the matrix element
  Complex ii(0.,1.);
  Lorentz5Momentum pff(momenta[1]+momenta[2]);
  pff.rescaleMass();
  Energy2 mff2(pff.mass()*pff.mass());
  // compute the prefactor
  complex<InvEnergy3> pre(coupling_[imode()]/mff2);
  // the VMD factor
  if(includeVMD_[imode()]>0) {
    Energy2 mrho2=VMDmass_[imode()]*VMDmass_[imode()];
    Energy2 mwrho=VMDmass_[imode()]*VMDwidth_[imode()];
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  LorentzVector<complex<Energy3> > eps;
  LorentzVector<complex<Energy> > fcurrent;
  // compute the matrix element
  vector<unsigned int> ispin(4);ispin[0]=0;
  for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
    for(ispin[2]=0;ispin[2]<2;++ispin[2]) {
      fcurrent = wave_[ispin[3]].vectorCurrent(wavebar_[ispin[2]]);
      // compute the current for this part
      eps = epsilon(momenta[0],pff,fcurrent);
      for(ispin[1]=0;ispin[1]<3;++ispin[1]) {
	(*ME())(ispin) = Complex(pre *vectors_[ispin[1]].dot(eps));
      }
    }	  
  }
  double me = ME()->contract(rho_).real();
  //   //code to test the matrix element against the analytic result
  //   Energy   m[4]={part.mass(),momenta[0].mass(),momenta[1].mass(),momenta[2].mass()};
  //   Energy2 m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
  //   Lorentz5Momentum p12=momenta[0]+momenta[1];p12.rescaleMass();
  //   Energy2 m122(p12.mass2());
  //   Complex output( ((pre*conj(pre)).real()*(
  // 				 -2*m122*m122*mff2 - mff2*mff2*mff2 + 
  // 				 m2[1]*(2*m2[2]*m2[3] - 2*m2[3]*m2[3] + 
  // 					m2[1]*(m2[2] - 2*m[2]*m[3] - m2[3])) - 
  // 				 2*m[2]*(m[2]*m2[2] - 2*m2[1]*m[3] - m[2]*m2[3])*
  // 				 m2[0] - (m2[2] + 2*m[2]*m[3] - m2[3])*
  // 				 m2[0]*m2[0] +mff2*mff2*
  // 				 (2*m2[1] + (m[2] - m[3])*(m[2] - m[3]) + 2*m2[0]) - 
  // 				 mff2*(m2[1]*m2[1] + 2*m2[1]*m[2]*(m[2] - 2*m[3]) + 
  // 				       2*m2[2]*m2[3] - 2*(2*m[2] - m[3])*m[3]*m2[0] + 
  // 				       m2[0]*m2[0]) + 2*m122*
  // 				 (-mff2*mff2 - (m2[2] - m2[3])*
  // 				  (m2[1] - m2[0]) + 
  // 				  mff2*(m2[1] + m2[2] + m2[3] + 
  // 					m2[0])))));
  //   cout << "testing the matrix element " 
  //        << real(output) << " " << me << " " << test2 << endl;
  return me;
}

// method to return an object to calculate the 3 or higher body partial width
WidthCalculatorBasePtr 
PScalarVectorFermionsDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int imode(-1);
  // ids of the particles
  int id0(dm.parent()->id()),idf[2]={0,0},idv(0);
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit) {
    if((**pit).iSpin()==PDT::Spin1){idv=(**pit).id();}
    else{idf[nf]=(**pit).id();++nf;}
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==idv) {
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
    new_ptr(ThreeBodyAllOn1IntegralCalculator<PScalarVectorFermionsDecayer>
	    (3,-1000.*MeV,-0.9*MeV,-0.9,*this,imode,m[0],m[1],m[2]));
}

InvEnergy PScalarVectorFermionsDecayer::threeBodydGammads(const int imodeb,
							  const Energy2 q2, 
							  const Energy2 mff2, 
							  const Energy m1,
							  const Energy m2,
							  const Energy m3) const {
  // the masses of the external particles
  Energy q=sqrt(q2);
  Energy2 m12=m1*m1;
  Energy2 m22=m2*m2;
  Energy2 m32=m3*m3;
  // calculate the prefactor
  Complex ii(0.,1.);
  complex<InvEnergy3> pre = coupling_[imodeb] / mff2;
  // the VMD factor
  if(includeVMD_[imodeb]>0) {
    Energy2 mrho2=VMDmass_[imodeb]*VMDmass_[imodeb];
    Energy2 mwrho=VMDmass_[imodeb]*VMDwidth_[imodeb];
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  InvEnergy6 factor=real(pre*conj(pre));
  // compute the pieces from the integration limits
  Energy mff=sqrt(mff2);
  Energy e2star = 0.5*(mff2-m32+m22)/mff;
  Energy e1star = 0.5*(q2-mff2-m12)/mff;
  Energy e1sm = sqrt(e1star*e1star-m12);
  Energy e2sm = sqrt(e2star*e2star-m22);
  Energy2 a = 2*e1star*e2star+m12+m22;
  Energy2 b = 2*e1sm*e2sm;
  // term independent of s3
  Energy8 me = 2*b*(2*(m12*(mff2*mff2 + 4*mff2*m2*m3 -(m22 - m32)*(m22 - m32)) + 
		       2*m2*(m12 +m22)*m3*(-mff2 +m22 + q2))
		    +(m12 +m22)*(m12 +m22)*(-mff2 +m22 - 2*m2*m3 - m32)
		    -(mff2 +m22 + 2*m2*m3 - m32)*(-mff2 +m22 + q2)*(-mff2 +m22 + q2));
  // linear term
  me+= 2.*a*b*(2*(-mff2*mff2 - (m22 - m32)*(m12 - q2) + 
		  mff2*(m12 + m22 + m32 + q2)));
  // quadratic term
  me+=-4.*mff2*b*(3.*a*a+b*b)/3.;

  // phase space factors
  using Constants::pi;
  return -factor * me/256./pi/pi/pi/q2/q;
}

// output the setup information for the particle database
void PScalarVectorFermionsDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*GeV << includeVMD_[ix] << " "
	   << VMDid_[ix] << " " << VMDmass_[ix]/GeV << " " << VMDwidth_[ix]/GeV << " "
	   << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string PScalarVectorFermionsDecayer::setUpDecayMode(string arg) {
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
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
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
