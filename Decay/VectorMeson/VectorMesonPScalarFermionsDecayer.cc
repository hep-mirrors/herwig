// -*- C++ -*-
//
// VectorMesonPScalarFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonPScalarFermionsDecayer class.
//
//  Author: Peter Richardson
//

#include "VectorMesonPScalarFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
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

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonPScalarFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix) {
      _maxweight[ix] = mode(ix)->maxWeight();
      _weight[ix]    = mode(ix)->channelWeight(1);
    }
  }
}

VectorMesonPScalarFermionsDecayer::VectorMesonPScalarFermionsDecayer() 
  : _coupling(6), _incoming(6), _outgoingP(6), _outgoingf(6), _outgoinga(6), 
    _maxweight(6), _weight(6), _includeVMD(6), _VMDid(6), _VMDmass(6), 
    _VMDwidth(6) {
  // omega -> pi e+e- /mu+mu-
  _incoming[0] =  223; _outgoingP[0] =  111; 
  _outgoingf[0] = 11; _outgoinga[0] = -11; 
  _coupling[0] = 0.2179/GeV; _maxweight[0] = 4.0; _weight[0] = 0.; 
  _includeVMD[0] = 2; _VMDid[0] = 113; 
  _VMDmass[0] = 0.7758*GeV; _VMDwidth[0] = 0.1503*GeV; 
  _incoming[1] =  223; _outgoingP[1] =  111; 
  _outgoingf[1] = 13; _outgoinga[1] = -13; 
  _coupling[1] = 0.2179/GeV; _maxweight[1] = 2.8; _weight[1] = 0.; 
  _includeVMD[1] = 2; _VMDid[1] = 113; 
  _VMDmass[1] = 0.7758*GeV; _VMDwidth[1] = 0.1503*GeV; 
  // phi -> eta e+e-/mu+mu-
  _incoming[2] =  333; _outgoingP[2] =  221; 
  _outgoingf[2] = 11; _outgoinga[2] = -11; 
  _coupling[2] = 0.0643/GeV; _maxweight[2] = 3.7; _weight[2] = 0.; 
  _includeVMD[2] = 2; _VMDid[2] = 113; 
  _VMDmass[2] = 0.7758*GeV; _VMDwidth[2] = 0.1503*GeV; 
  _incoming[3] =  333; ;_outgoingP[3] =  221; 
  _outgoingf[3] = 13; _outgoinga[3] = -13; 
  _coupling[3] = 0.0643/GeV; _maxweight[3] = 2.8; _weight[3] = 0.; 
  _includeVMD[3] = 2; _VMDid[3] = 113; 
  _VMDmass[3] = 0.7758*GeV; _VMDwidth[3] = 0.1503*GeV; 
  // phi -> pi e+e-/mu+mu-
  _incoming[4] =  333; ;_outgoingP[4] =  111; 
  _outgoingf[4] = 11; _outgoinga[4] = -11; 
  _coupling[4] = 0.0120094/GeV; _maxweight[4] = 4.9; _weight[4] = 0.21; 
  _includeVMD[4] = 2; _VMDid[4] = 113; 
  _VMDmass[4] = 0.7758*GeV; _VMDwidth[4] = 0.1503*GeV; 
  _incoming[5] =  333; ;_outgoingP[5] =  111; 
  _outgoingf[5] = 13; _outgoinga[5] = -13; 
  _coupling[5] = 0.0120094/GeV; _maxweight[5] = 3.2; _weight[5] = 0.33; 
  _includeVMD[5] = 2; _VMDid[5] = 113; 
  _VMDmass[5] = 0.7758*GeV; _VMDwidth[5] = 0.1503*GeV; 
  // the initial size of the arrays
  _initsize=_coupling.size();
  // intermediates
  generateIntermediates(false);
}

void VectorMesonPScalarFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize(_coupling.size());
  if(isize!=_incoming.size()  || isize!=_outgoingP.size()|| isize!=_outgoingf.size()||
     isize!=_outgoinga.size() || isize!=_maxweight.size()|| isize!=_includeVMD.size()||
     isize!=_VMDid.size()     || isize!=_VMDmass.size()  || isize!=_VMDwidth.size()||
     isize!=_weight.size())
    throw InitException() << "Inconsistent parameters in VectorMesonPScalar"
			  << "FermionsDecayer" << Exception::abortnow;
  // create the integration channel for each mode 
  tPDVector extpart(4);
  tPDPtr gamma(getParticleData(ParticleID::gamma)),rho;
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr newmode;
  vector<double> wgt(2,1.);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    rho=getParticleData(_VMDid[ix]);
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoingP[ix]);
    extpart[2] = getParticleData(_outgoingf[ix]);
    extpart[3] = getParticleData(_outgoinga[ix]);
    newmode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    // photon channel
    newchannel=new_ptr(DecayPhaseSpaceChannel(newmode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(gamma     ,1,-1.1, 2,3);
    newmode->addChannel(newchannel);
    wgt[0]=1.-_weight[ix];
    // vmd channel
    newchannel=new_ptr(DecayPhaseSpaceChannel(newmode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho       ,0, 0.0, 2,3);
    newmode->addChannel(newchannel);
    wgt[1]=_weight[ix];
    addMode(newmode,_maxweight[ix],wgt);
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
    if(_incoming[ix]==id0&&_outgoingP[ix]==ids) {
      if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	 (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  cc=false;
  return imode;
}

void VectorMesonPScalarFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/GeV) << _incoming << _outgoingP << _outgoingf << _outgoinga 
     << _maxweight << _weight << _includeVMD << _VMDid << ounit(_VMDmass,GeV) 
     << ounit(_VMDwidth,GeV);
}

void VectorMesonPScalarFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/GeV) >> _incoming >> _outgoingP >> _outgoingf >> _outgoinga 
     >> _maxweight >> _weight >> _includeVMD >> _VMDid >> iunit(_VMDmass,GeV) 
     >> iunit(_VMDwidth,GeV);
}

ClassDescription<VectorMesonPScalarFermionsDecayer> 
VectorMesonPScalarFermionsDecayer::initVectorMesonPScalarFermionsDecayer;
// Definition of the static class description member.

void VectorMesonPScalarFermionsDecayer::Init() {

  static ClassDocumentation<VectorMesonPScalarFermionsDecayer> documentation
    ("The VectorMesonPScalarFermionsDecayer class is designed to "
     "perform the decay of a vector meson to a pseudoscalar meson and a "
     "fermion-antifermion pair.");

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonPScalarFermionsDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingP
    ("OutgoingPseudoScalar",
     "The PDG code for the outgoing pseudoscalar",
     &VectorMesonPScalarFermionsDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingF
    ("OutgoingFermion",
     "The PDG code for the outgoing fermion",
     &VectorMesonPScalarFermionsDecayer::_outgoingf,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceOutcomingA
    ("OutgoingAntiFermion",
     "The PDG code for the outgoing antifermion",
     &VectorMesonPScalarFermionsDecayer::_outgoinga,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonPScalarFermionsDecayer::_coupling,
     1/MeV, 0, ZERO, -10000000/MeV, 10000000/MeV, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonPScalarFermionsDecayer::_maxweight,
     0, 0, 0, 0., 200.0, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,double> interfaceWeight
    ("Weight",
     "The weight for vector meson phasse space channel",
     &VectorMesonPScalarFermionsDecayer::_weight,
     0, 0, 0, 0., 200.0, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceIncludeVMD
    ("IncludeVMD",
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &VectorMesonPScalarFermionsDecayer::_includeVMD,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,int> interfaceVMDID
    ("VMDID",
     "The PDG code for the particle to be used for the VMD factor.",
     &VectorMesonPScalarFermionsDecayer::_VMDid,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,Energy> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &VectorMesonPScalarFermionsDecayer::_VMDmass,
     MeV, 0, ZERO, -10000000*MeV, 10000000*MeV, false, false, true);

  static ParVector<VectorMesonPScalarFermionsDecayer,Energy> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &VectorMesonPScalarFermionsDecayer::_VMDwidth,
     MeV, 0, ZERO, -10000000*MeV, 10000000*MeV, false, false, true);

}

double VectorMesonPScalarFermionsDecayer::me2(const int,
					      const Particle & inpart,
					      const ParticleVector & decay,
					      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,
					 PDT::Spin1Half,PDT::Spin1Half)));
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[1],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[2],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[1],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[2],outgoing);
  // the factor for the off-shell photon
  Lorentz5Momentum pff(decay[1]->momentum()+decay[2]->momentum());
  pff.rescaleMass();
  Energy2 mff2(pff.mass2());
  // prefactor
  complex<InvEnergy3> pre(_coupling[imode()]/mff2);
  Complex ii(0.,1.);
  // the VMD factor
  if(_includeVMD[imode()]>0) {
    Energy2 mrho2(_VMDmass[imode()]*_VMDmass[imode()]);
    Energy2 mwrho(_VMDmass[imode()]*_VMDwidth[imode()]);
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  // calculate the matrix element
  LorentzPolarizationVector temp;
  unsigned int ix,iy,iz;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      temp=pre*epsilon(inpart.momentum(),pff,
				    _wave[ix].vectorCurrent(_wavebar[iy]));
      for(iz=0;iz<3;++iz) 
	(*ME())(iz,0,iy,ix)=temp.dot(_vectors[iz]); 
    }
  }
  /* code for the spin averaged me for testing only
  Energy  m[4]={inpart.mass(),decay[0]->mass(),decay[1]->mass(),decay[2]->mass()};
  Energy m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
  Lorentz5Momentum p12=decay[0]->momentum()+decay[1]->momentum();p12.rescaleMass();
  Energy m122(p12.mass2());
  cout << "testing the matrix element " 
       << -1./3.*(pre*conj(pre)).real()*(-2*m122*m122*mff2 - mff2*mff2*mff2 + 
		   m2[1]*(2*m2[2]*m2[3] - 2*m2[3]*m2[3] + 
			  m2[1]*(m2[2] - 2*m[2]*m[3] - m2[3])) - 
		   2*m[2]*(m2[2]*m[2] - 2*m2[1]*m[3] - m[2]*m2[3])*
		   m2[0] - (m2[2] + 2*m[2]*m[3] - m2[3])*
		   m2[0]*m2[0] + mff2*mff2*
		   (2*m2[1] + (m[2] - m[3])*(m[2] - m[3]) + 2*m2[0]) - 
		   mff2*(m2[1]*m2[1] + 2*m2[1]*m[2]*(m[2] - 2*m[3]) + 
			 2*m2[2]*m2[3] - 2*(2*m[2] - m[3])*m[3]*m2[0] + 
			 m2[0]*m2[0]) + 2*m122*
		   (-mff2*mff2 - (m[2] - m[3])*(m[2] + m[3])*(m[1] - m[0])*
		    (m[1] + m[0]) + mff2*
		    (m2[1] + m2[2] + m2[3] + m2[0])))
       << endl;
   */
  // return the answer
  return ME()->contract(_rho).real();
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
    if(_incoming[ix]==id0&&_outgoingP[ix]==ids) {
      if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	 (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  // get the masses we need
  Energy m[3]={getParticleData(_outgoingP[imode])->mass(),
	       getParticleData(_outgoingf[imode])->mass(),
	       getParticleData(_outgoinga[imode])->mass()};
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
  complex<InvEnergy3> pre(_coupling[imodeb]/mff2);
  Complex ii(0.,1.);
  // the VMD factor
  if(_includeVMD[imodeb]>0) {
    Energy2 mrho2(_VMDmass[imodeb]*_VMDmass[imodeb]),
      mwrho(_VMDmass[imodeb]*_VMDwidth[imodeb]);
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
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming   " << ix << "  " 
	     << _incoming[ix]   << "\n";
      output << "newdef " << name() << ":OutgoingPseudoScalar  " 
	     << ix << "  " << _outgoingP[ix]  << "\n";
      output << "newdef " << name() << ":OutgoingFermion  " 
	     << ix << "  " << _outgoingf[ix]  << "\n";
      output << "newdef " << name() << ":OutgoingAntiFermion " 
	     << ix << "  " << _outgoinga[ix]  << "\n";
      output << "newdef " << name() << ":Coupling   " << ix << "  " 
	     << _coupling[ix]*MeV   << "\n";
      output << "newdef " << name() << ":MaxWeight  " << ix << "  " 
	     << _maxweight[ix]  << "\n";
      output << "newdef " << name() << ":Weight  "    << ix << "  " 
	     << _weight[ix]  << "\n";
      output << "newdef " << name() << ":IncludeVMD " << ix << "  " 
	     << _includeVMD[ix] << "\n";
      output << "newdef " << name() << ":VMDID      " << ix << "  " 
	     << _VMDid[ix]      << "\n";
      output << "newdef " << name() << ":VMDmass    " << ix << "  " 
	     << _VMDmass[ix]/MeV    << "\n";
      output << "newdef " << name() << ":VMDwidth   " << ix << "  " 
	     << _VMDwidth[ix]/MeV   << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming   " << ix << "  " 
	     << _incoming[ix]   << "\n";
      output << "insert " << name() << ":OutgoingPseudoScalar  " 
	     << ix << "  " << _outgoingP[ix]  << "\n";
      output << "insert " << name() << ":OutgoingFermion  " 
	     << ix << "  " << _outgoingf[ix]  << "\n";
      output << "insert " << name() << ":OutgoingAntiFermion " 
	     << ix << "  " << _outgoinga[ix]  << "\n";
      output << "insert " << name() << ":Coupling   " << ix << "  " 
	     << _coupling[ix]*MeV   << "\n";
      output << "insert " << name() << ":Weight  "    << ix << "  " 
	     << _weight[ix]  << "\n";
      output << "insert " << name() << ":IncludeVMD " << ix << "  " 
	     << _includeVMD[ix] << "\n";
      output << "insert " << name() << ":VMDID      " << ix << "  " 
	     << _VMDid[ix]      << "\n";
      output << "insert " << name() << ":VMDmass    " << ix << "  " 
	     << _VMDmass[ix]/MeV    << "\n";
      output << "insert " << name() << ":VMDwidth   " << ix << "  " 
	     << _VMDwidth[ix]/MeV   << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
