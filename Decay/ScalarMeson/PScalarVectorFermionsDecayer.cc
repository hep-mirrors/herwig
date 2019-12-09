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
#include "ThePEG/Interface/ParVector.h"
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
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      _maxweight[ix]=mode(ix)->maxWeight();
  }
}

PScalarVectorFermionsDecayer::PScalarVectorFermionsDecayer() 
  : _coupling(5), _incoming(5), _outgoingV(5), _outgoingf(5), 
    _outgoinga(5),_maxweight(5), _includeVMD(5), _VMDid(5),
    _VMDmass(5), _VMDwidth(5) {
  // pi0 -> gamma e+e-
  _incoming[0] = 111;_outgoingV[0] =  22;
  _outgoingf[0] = 11;_outgoinga[0] = -11;
  _coupling[0] = 0.00761872/GeV;_maxweight[0] = 0.027;
  _includeVMD[0] = 2;_VMDid[0] = 113;
  _VMDmass[0] = 0.7758*GeV;_VMDwidth[0] = 0.1503*GeV;
  // eta -> gamma e+e-/mu+/mu-
  _incoming[1] = 221;_outgoingV[1] =  22;
  _outgoingf[1] = 11;_outgoinga[1] = -11;
  _coupling[1] = 0.007554164/GeV;_maxweight[1] = 2.8;
  _includeVMD[1] = 2;_VMDid[1] = 113;
  _VMDmass[1] = 0.7758*GeV;_VMDwidth[1] = 0.1503*GeV;
  _incoming[2] = 221;_outgoingV[2] =  22;
  _outgoingf[2] = 13;_outgoinga[2] = -13;
  _coupling[2] = 0.007554164/GeV;_maxweight[2] = 2.1;
  _includeVMD[2] = 2;_VMDid[2] = 113;
  _VMDmass[2] = 0.7758*GeV;_VMDwidth[2] = 0.1503*GeV;
  // eta' -> gamma e+e-/mu+mu-
  _incoming[3] = 331;_outgoingV[3] =  22;
  _outgoingf[3] = 11;_outgoinga[3] = -11;
  _coupling[3] = 0.0104/GeV;_maxweight[3] = 5.2;
  _includeVMD[3] = 2;_VMDid[3] = 113;
  _VMDmass[3] = 0.7758*GeV;_VMDwidth[3] = 0.1503*GeV;
  _incoming[4] = 331;_outgoingV[4] =  22;
  _outgoingf[4] = 13;_outgoinga[4] = -13;
  _coupling[4] = 0.0104/GeV;_maxweight[4] = 3.0;
  _includeVMD[4] = 2;_VMDid[4] = 113;
  _VMDmass[4] = 0.7758*GeV;_VMDwidth[4] = 0.1503*GeV;
  // initial size of the arrays
  _initsize = _incoming.size();
  // intermediates
  generateIntermediates(false);
}

void PScalarVectorFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoingV.size()|| isize!=_outgoingf.size()||
     isize!=_outgoinga.size() || isize!=_maxweight.size()|| isize!=_includeVMD.size()||
     isize!=_VMDid.size()     || isize!=_VMDmass.size()  || isize!=_VMDwidth.size())
    throw InitException() << "Inconsistent parameters in PScalarVectorFermionsDecayer"
			  << Exception::abortnow;
  // create the integration channel for each mode 
  tPDVector extpart(4);
  tPDPtr gamma(getParticleData(ParticleID::gamma));
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  vector<double> wgt(1,1.);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoingV[ix]);
    extpart[2] = getParticleData(_outgoingf[ix]);
    extpart[3] = getParticleData(_outgoinga[ix]);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(gamma     ,1,-1.1, 2,3);
    mode->addChannel(newchannel);
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
    if(_incoming[ix]==id0&&_outgoingV[ix]==idv)
      {if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	  (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])){imode=ix;}}
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
}

void PScalarVectorFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/MeV) 
     << _incoming << _outgoingV << _outgoingf << _outgoinga << _maxweight
     << _includeVMD << _VMDid 
     << ounit(_VMDmass,MeV) << ounit(_VMDwidth,MeV);
}

void PScalarVectorFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/MeV) 
     >> _incoming >> _outgoingV >> _outgoingf >> _outgoinga >> _maxweight
     >> _includeVMD >> _VMDid 
     >> iunit(_VMDmass,MeV) >> iunit(_VMDwidth,MeV);
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

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarVectorFermionsDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingV
    ("OutgoingVector",
     "The PDG code for the outgoing pseudoscalar",
     &PScalarVectorFermionsDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingF
    ("OutgoingFermion",
     "The PDG code for the outgoing fermion",
     &PScalarVectorFermionsDecayer::_outgoingf,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,int> interfaceOutcomingA
    ("OutgoingAntiFermion",
     "The PDG code for the outgoing antifermion",
     &PScalarVectorFermionsDecayer::_outgoinga,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarVectorFermionsDecayer::_coupling,
     1/MeV, 0, ZERO, -10000000/MeV, 10000000/MeV, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarVectorFermionsDecayer::_maxweight,
     0, 0, 0, 0.0, 100., false, false, true);

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
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,Energy> interfaceVMDmass
    ("VMDmass",
     "The mass to use for the particle in the VMD factor",
     &PScalarVectorFermionsDecayer::_VMDmass,
     MeV, 0, ZERO, ZERO, 10000.*MeV, false, false, true);

  static ParVector<PScalarVectorFermionsDecayer,Energy> interfaceVMDwidth
    ("VMDwidth",
     "The width to use for the particle in the VMD factor",
     &PScalarVectorFermionsDecayer::_VMDwidth,
     MeV, 0, ZERO, ZERO, 10000.*MeV, false, false, true);

}

double PScalarVectorFermionsDecayer::me2(const int,
					 const Particle & inpart,
					 const ParticleVector & decay,
					 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1Half,
					 PDT::Spin1Half)));
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    // set up the spin information for the decay products
    VectorWaveFunction::
      constructSpinInfo(_vectors,decay[0],outgoing,true,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[1],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[2],outgoing,true);
    return 0.;
  }
  // calculate the spinors and polarization vectors
  VectorWaveFunction::
    calculateWaveFunctions(_vectors,decay[0],outgoing,true);
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[1],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[2],outgoing);
  // now compute the matrix element
  Complex ii(0.,1.);
  Lorentz5Momentum pff(decay[1]->momentum()+decay[2]->momentum());
  pff.rescaleMass();
  Energy2 mff2(pff.mass()*pff.mass());
  // compute the prefactor
  complex<InvEnergy3> pre(_coupling[imode()]/mff2);
  // the VMD factor
  if(_includeVMD[imode()]>0) {
    Energy2 mrho2=_VMDmass[imode()]*_VMDmass[imode()];
    Energy2 mwrho=_VMDmass[imode()]*_VMDwidth[imode()];
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  LorentzVector<complex<Energy3> > eps;
  LorentzVector<complex<Energy> > fcurrent;
  // compute the matrix element
  vector<unsigned int> ispin(4);ispin[0]=0;
  for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
    for(ispin[2]=0;ispin[2]<2;++ispin[2]) {
      fcurrent = _wave[ispin[3]].vectorCurrent(_wavebar[ispin[2]]);
      // compute the current for this part
      eps = epsilon(decay[0]->momentum(),pff,fcurrent);
      for(ispin[1]=0;ispin[1]<3;++ispin[1]) {
	(*ME())(ispin) = Complex(pre *_vectors[ispin[1]].dot(eps));
      }
    }	  
  }
  double me = ME()->contract(_rho).real();
//   //code to test the matrix element against the analytic result
//   Energy   m[4]={inpart.mass(),decay[0]->mass(),decay[1]->mass(),decay[2]->mass()};
//   Energy2 m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
//   Lorentz5Momentum p12=decay[0]->momentum()+decay[1]->momentum();p12.rescaleMass();
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
    if(_incoming[ix]==id0&&_outgoingV[ix]==idv) {
      if((idf[0]==_outgoingf[ix]&&idf[1]==_outgoinga[ix])||
	 (idf[1]==_outgoingf[ix]&&idf[0]==_outgoinga[ix])) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  // get the masses we need
  Energy m[3]={getParticleData(_outgoingV[imode])->mass(),
	       getParticleData(_outgoingf[imode])->mass(),
	       getParticleData(_outgoinga[imode])->mass()};
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
  complex<InvEnergy3> pre = _coupling[imodeb] / mff2;
  // the VMD factor
  if(_includeVMD[imodeb]>0) {
    Energy2 mrho2=_VMDmass[imodeb]*_VMDmass[imodeb];
    Energy2 mwrho=_VMDmass[imodeb]*_VMDwidth[imodeb];
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
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming   " << ix << "  " 
	     << _incoming[ix]   << "\n";
      output << "newdef " << name() << ":OutgoingVector  " 
	     << ix << "  " << _outgoingV[ix]  << "\n";
      output << "newdef " << name() << ":OutgoingFermion  " 
	     << ix << "  " << _outgoingf[ix]  << "\n";
      output << "newdef " << name() << ":OutgoingAntiFermion " 
	     << ix << "  " << _outgoinga[ix]  << "\n";
      output << "newdef " << name() << ":Coupling   " << ix << "  " 
	     << _coupling[ix]*MeV   << "\n";
      output << "newdef " << name() << ":MaxWeight  " << ix << "  " 
	     << _maxweight[ix]  << "\n";
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
      output << "insert " << name() << ":OutgoingVector  " 
	     << ix << "  " << _outgoingV[ix]  << "\n";
      output << "insert " << name() << ":OutgoingFermion  " 
	     << ix << "  " << _outgoingf[ix]  << "\n";
      output << "insert " << name() << ":OutgoingAntiFermion " 
	     << ix << "  " << _outgoinga[ix]  << "\n";
      output << "insert " << name() << ":Coupling   " << ix << "  " 
	     << _coupling[ix]*MeV   << "\n";
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
