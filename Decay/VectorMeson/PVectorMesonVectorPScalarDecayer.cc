// -*- C++ -*-
//
// PVectorMesonVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PVectorMesonVectorPScalarDecayer class.
//

#include "PVectorMesonVectorPScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
 
void PVectorMesonVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

void PVectorMesonVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingV.size()||isize!=_outgoingP.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters in "
			  << "PVectorMesonVectorPScalarDecayer::doinit()" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoingV[ix]);
    extpart[2]=getParticleData(_outgoingP[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

PVectorMesonVectorPScalarDecayer::PVectorMesonVectorPScalarDecayer() 
  : _coupling(67), _incoming(67), _outgoingV(67), _outgoingP(67), 
    _maxweight(67) {
  // decay mode h'_1 to K K*
  _incoming[0] =  10333; _outgoingV[0] =  313; _outgoingP[0] = -311; 
  _coupling[0] = 4.889/GeV; _maxweight[0] = 11.; 
  _incoming[1] =  10333; _outgoingV[1] =  323; _outgoingP[1] = -321; 
  _coupling[1] = 4.156/GeV; _maxweight[1] = 11.; 
  // decay mode h_1 to rho
  _incoming[2] = 10223; _outgoingV[2] =  113; _outgoingP[2] =  111; 
  _coupling[2] = 4.411/GeV; _maxweight[2] = 6.; 
  _incoming[3] = 10223; _outgoingV[3] = -213; _outgoingP[3] =  211; 
  _coupling[3] = 4.411/GeV; _maxweight[3] = 6.; 
  // decay mode b_1 to omega pi
  _incoming[4] =  10113; _outgoingV[4] =  223; _outgoingP[4] =  111; 
  _coupling[4] = 3.862/GeV; _maxweight[4] = 7.5; 
  _incoming[5] =  10213; _outgoingV[5] =  223; _outgoingP[5] =  211; 
  _coupling[5] = 3.862/GeV; _maxweight[5] = 7.5; 
  // decay mode b_1^+ to gamma pi
  _incoming[6] =  10213; _outgoingV[6] =   22; _outgoingP[6] =  211; 
  _coupling[6] = 0.195/GeV; _maxweight[6] = 2.; 
  // decay mode D'_s1 to D*K
  _incoming[7] =  20433; _outgoingV[7] =  413; _outgoingP[7] =  311; 
  _coupling[7] = 0.161/GeV; _maxweight[7] = 2.5; 
  _incoming[8] =  20433; _outgoingV[8] =  423; _outgoingP[8] =  321; 
  _coupling[8] = 0.161/GeV; _maxweight[8] = 2.5; 
  // decay mode B'_s1 to B*K
  _incoming[9] =  20533; _outgoingV[9] =  513; _outgoingP[9] = -311; 
  _coupling[9] = 0.389/GeV; _maxweight[9] = 2.1; 
  _incoming[10] =  20533; _outgoingV[10] =  523; _outgoingP[10] = -321; 
  _coupling[10] = 0.389/GeV; _maxweight[10] = 2.1; 
  // decay mode K_1 to rho K
  _incoming[11] =  10323; _outgoingV[11] =  213; _outgoingP[11] =  311; 
  _coupling[11] = 4.98/GeV; _maxweight[11] = 3.5; 
  _incoming[12] =  10323; _outgoingV[12] =  113; _outgoingP[12] =  321; 
  _coupling[12] = 3.40/GeV; _maxweight[12] = 3.5; 
  _incoming[13] =  10313; _outgoingV[13] = -213; _outgoingP[13] =  321; 
  _coupling[13] = 4.87/GeV; _maxweight[13] = 3.5; 
  _incoming[14] =  10313; _outgoingV[14] =  113; _outgoingP[14] =  311; 
  _coupling[14] = 3.55/GeV; _maxweight[14] = 3.5; 
  // decay mode K'_1 to rho K
  _incoming[15] =  20323; _outgoingV[15] =  213; _outgoingP[15] =  311; 
  _coupling[15] = 0.97/GeV; _maxweight[15] = 6.5; 
  _incoming[16] =  20323; _outgoingV[16] =  113; _outgoingP[16] =  321; 
  _coupling[16] = 0.69/GeV; _maxweight[16] = 6.5; 
  _incoming[17] =  20313; _outgoingV[17] = -213; _outgoingP[17] =  321; 
  _coupling[17] = 0.97/GeV; _maxweight[17] = 6.5; 
  _incoming[18] =  20313; _outgoingV[18] =  113; _outgoingP[18] =  311; 
  _coupling[18] = 0.707/GeV; _maxweight[18] = 6.5; 
  // decay mode K_1 to omega K
  _incoming[19] =  10323; _outgoingV[19] =  223; _outgoingP[19] =  321; 
  _coupling[19] = 4.76/GeV; _maxweight[19] = 7.5; 
  _incoming[20] =  10313; _outgoingV[20] =  223; _outgoingP[20] =  311; 
  _coupling[20] = 6.0/GeV; _maxweight[20] = 7.5; 
  // decay mode K'_1 to omega K
  _incoming[21] =  20323; _outgoingV[21] =  223; _outgoingP[21] =  321; 
  _coupling[21] = 0.600/GeV; _maxweight[21] = 8.; 
  _incoming[22] =  20313; _outgoingV[22] =  223; _outgoingP[22] =  311; 
  _coupling[22] = 0.600/GeV; _maxweight[22] = 8.; 
  // decay mode K_1 to K* pi
  _incoming[23] =  10323; _outgoingP[23] =  211; _outgoingV[23] =  313; 
  _coupling[23] = 0.941/GeV; _maxweight[23] = 8.5; 
  _incoming[24] =  10323; _outgoingP[24] =  111; _outgoingV[24] =  323; 
  _coupling[24] = 0.656/GeV; _maxweight[24] = 8.5; 
  _incoming[25] =  10313; _outgoingP[25] = -211; _outgoingV[25] =  323; 
  _coupling[25] = 0.932/GeV; _maxweight[25] = 8.5; 
  _incoming[26] =  10313; _outgoingP[26] =  111; _outgoingV[26] =  313; 
  _coupling[26] = 0.658/GeV; _maxweight[26] = 8.5; 
  // decay mode K'_1 to K* pi
  _incoming[27] =  20323; _outgoingP[27] =  211; _outgoingV[27] =  313; 
  _coupling[27] = 2.845/GeV; _maxweight[27] = 12.; 
  _incoming[28] =  20323; _outgoingP[28] =  111; _outgoingV[28] =  323; 
  _coupling[28] = 1.99/GeV; _maxweight[28] = 12.; 
  _incoming[29] =  20313; _outgoingP[29] = -211; _outgoingV[29] =  323; 
  _coupling[29] = 2.84/GeV; _maxweight[29] = 12.; 
  _incoming[30] =  20313; _outgoingP[30] =  111; _outgoingV[30] =  313; 
  _coupling[30] = 2.00/GeV; _maxweight[30] = 12.; 
  // decaymode D_1 to D* pi
  _incoming[31] =  10423; _outgoingP[31] = -211; _outgoingV[31] =  413; 
  _coupling[31] = 0.489/GeV; _maxweight[31] = 3.; 
  _incoming[32] =  10423; _outgoingP[32] =  111; _outgoingV[32] =  423; 
  _coupling[32] = 0.347/GeV; _maxweight[32] = 3.; 
  _incoming[33] =  10413; _outgoingP[33] =  211; _outgoingV[33] =  423; 
  _coupling[33] = 0.542/GeV; _maxweight[33] = 3.; 
  _incoming[34] =  10413; _outgoingP[34] =  111; _outgoingV[34] =  413; 
  _coupling[34] = 0.383/GeV; _maxweight[34] = 3.; 
  // decaymode D'_1 to D* pi
  _incoming[35] =  20423; _outgoingP[35] = -211; _outgoingV[35] =  413; 
  _coupling[35] = 1.933/GeV; _maxweight[35] = 3.; 
  _incoming[36] =  20423; _outgoingP[36] =  111; _outgoingV[36] =  423; 
  _coupling[36] = 1.367/GeV; _maxweight[36] = 3.; 
  _incoming[37] =  20413; _outgoingP[37] =  211; _outgoingV[37] =  423; 
  _coupling[37] = 1.926/GeV; _maxweight[37] = 3.; 
  _incoming[38] =  20413; _outgoingP[38] =  111; _outgoingV[38] =  413; 
  _coupling[38] = 1.367/GeV; _maxweight[38] = 3.; 
  // decaymode B_1 to B* pi
  _incoming[39] =  10523; _outgoingP[39] =  211; _outgoingV[39] =  513; 
  _coupling[39] = 0.130/GeV; _maxweight[39] = 2.2; 
  _incoming[40] =  10523; _outgoingP[40] =  111; _outgoingV[40] =  523; 
  _coupling[40] = 0.0924/GeV; _maxweight[40] = 2.2; 
  _incoming[41] =  10513; _outgoingP[41] = -211; _outgoingV[41] =  523; 
  _coupling[41] = 0.130/GeV; _maxweight[41] = 2.2; 
  _incoming[42] =  10513; _outgoingP[42] =  111; _outgoingV[42] =  513; 
  _coupling[42] = 0.0924/GeV; _maxweight[42] = 2.2; 
  // decaymode B'_1 to B* pi
  _incoming[43] =  20523; _outgoingP[43] =  211; _outgoingV[43] =  513; 
  _coupling[43] = 0.445/GeV; _maxweight[43] = 2.2; 
  _incoming[44] =  20523; _outgoingP[44] =  111; _outgoingV[44] =  523; 
  _coupling[44] = 0.314/GeV; _maxweight[44] = 2.2; 
  _incoming[45] =  20513; _outgoingP[45] = -211; _outgoingV[45] =  523; 
  _coupling[45] = 0.445/GeV; _maxweight[45] = 2.2; 
  _incoming[46] =  20513; _outgoingP[46] =  111; _outgoingV[46] =  513; 
  _coupling[46] = 0.314/GeV; _maxweight[46] = 2.2; 
  // decaymode D_s1 to D* pi
  _incoming[47] =  10433; _outgoingP[47] =  111; _outgoingV[47] =  433; 
  _coupling[47] = 0.022/GeV; _maxweight[47] = 2.5; 
  // decaymode D_s1 to D gamma
  _incoming[48] =  10433; _outgoingP[48] =  431; _outgoingV[48] = 22; 
  _coupling[48] = 0.0587/GeV; _maxweight[48] = 2.1; 
  // decaymode B_s1 to B gamma
  _incoming[49] =  10533; _outgoingP[49] =  531; _outgoingV[49] = 22; 
  _coupling[49] = 0.142/GeV; _maxweight[49] = 2.;
  // decaymode B_s1 to B* pi
  _incoming[50] =  10533; _outgoingP[50] =  111; _outgoingV[50] =  533; 
  _coupling[50] = 0.0074/GeV; _maxweight[50] = 2.1; 
  // decaymode B_c1 to B_c gamma
  _incoming[51] =  10543; _outgoingP[51] =  541; _outgoingV[51] = 22; 
  _coupling[51] = 0.0759/GeV; _maxweight[51] = 2.2; 
  // decaymode B'_c1 to B_c gamma
  _incoming[52] =  20543; _outgoingP[52] =  541; _outgoingV[52] = 22; 
  _coupling[52] = 0.175/GeV; _maxweight[52] = 2.2; 
  // decaymode h_c to eta_c gamma
  _incoming[53] = 10443; _outgoingP[53] = 441; _outgoingV[53] = 22; 
  _coupling[53] = 0.329/GeV; _maxweight[53] = 4.; 
  // decaymode h_b to eta_b gamma
  _incoming[54] = 10553; _outgoingP[54] = 551; _outgoingV[54] = 22; 
  _coupling[54] = 0.0356/GeV; _maxweight[54] = 3.5; 
  // a_1 to K* K
  _incoming[55] = 20213; _outgoingP[55] = -311; _outgoingV[55] = 323; 
  _coupling[55] = 3.42/GeV; _maxweight[55] = 2.5; 
  _incoming[56] = 20213; _outgoingP[56] = 321; _outgoingV[56] = -313; 
  _coupling[56] = 3.42/GeV; _maxweight[56] = 2.5; 
  _incoming[57] = 20113; _outgoingP[57] = 321; _outgoingV[57] = -323; 
  _coupling[57] = 3.42/GeV; _maxweight[57] = 4.0; 
  _incoming[58] = 20113; _outgoingP[58] = 311; _outgoingV[58] = -313; 
  _coupling[58] = 3.42/GeV; _maxweight[58] = 4.0; 
  // a_1 to gamma pi
  _incoming[59] =  20113; _outgoingP[59] =  111; _outgoingV[59] = 22; 
  _coupling[59] = 0.01/GeV; _maxweight[59] = 2.; 
  _incoming[60] =  20213; _outgoingP[60] =  211; _outgoingV[60] = 22; 
  _coupling[60] = 0.01/GeV; _maxweight[60] = 2.; 
  // f'_1 to K* K
  _incoming[61] = 20333; _outgoingP[61] = 321; _outgoingV[61] = -323; 
  _coupling[61] = 1.637/GeV; _maxweight[61] = 7.; 
  _incoming[62] = 20333; _outgoingP[62] = 311; _outgoingV[62] = -313; 
  _coupling[62] = 1.737/GeV; _maxweight[62] = 7.; 
  // decay mode K_1 to gamma K
  _incoming[63] =  10313; _outgoingV[63] =  22; _outgoingP[63] =  311; 
  _coupling[63] = 0.119/GeV; _maxweight[63] = 7.5; 
  // decay mode K'_1 to gamma K
  _incoming[64] =  20313; _outgoingV[64] =  22; _outgoingP[64] =  311; 
  _coupling[64] = 0.220/GeV; _maxweight[64] = 8.;  
  // decaymode B_s1 to B K
  _incoming[65] =  10533; _outgoingP[65] =  -311; _outgoingV[65] = 513; 
  _coupling[65] = 0.0418/GeV; _maxweight[65] = 2.;
  _incoming[66] =  10533; _outgoingP[66] =  -321; _outgoingV[66] = 523; 
  _coupling[66] = 0.0373/GeV; _maxweight[66] = 2.;  
  // initial size of the arrays
  _initsize = _coupling.size();
  // intermediates
  generateIntermediates(false);
}

int PVectorMesonVectorPScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
						 const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==_incoming[ix]) {
      if((id1   ==_outgoingV[ix]&&id2   ==_outgoingP[ix])||
	 (id2   ==_outgoingV[ix]&&id1   ==_outgoingP[ix])) imode=ix;
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if((id1bar==_outgoingV[ix]&&id2bar==_outgoingP[ix])||
	 (id2bar==_outgoingV[ix]&&id1bar==_outgoingP[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}


void PVectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoingV << _outgoingP << _maxweight << ounit(_coupling,1/GeV);
}

void PVectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> iunit(_coupling,1/GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PVectorMesonVectorPScalarDecayer,DecayIntegrator>
describeHerwigPVectorMesonVectorPScalarDecayer("Herwig::PVectorMesonVectorPScalarDecayer", "HwVMDecay.so");

void PVectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<PVectorMesonVectorPScalarDecayer> documentation
    ("The PVectorMesonVectorPScalarDecayer class is designed for the "
     "decay of a pseudovector meson to a vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PVectorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1 particle",
     &PVectorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceOutcomingPScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &PVectorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PVectorMesonVectorPScalarDecayer::_coupling,
     1/GeV, 0, ZERO, ZERO, 100./GeV, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PVectorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, 0., 10000., false, false, true);

}


double PVectorMesonVectorPScalarDecayer::me2(const int,
					     const Particle & inpart,
					     const ParticleVector& decay,
					     MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  // is the vector massless
  bool photon(_outgoingV[imode()]==ParticleID::gamma);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors[0],_rho,
						const_ptr_cast<tPPtr>(&inpart),
						incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    VectorWaveFunction::constructSpinInfo(_vectors[1],decay[0],
					  outgoing,true,photon);
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  VectorWaveFunction::calculateWaveFunctions(_vectors[1],decay[0],outgoing,photon);
  // compute the matrix element
  Energy2 p0dotpv(inpart.momentum()*decay[0]->momentum());
  complex<Energy> epsdot(ZERO);
  InvEnergy2 pre(_coupling[imode()]/inpart.mass());
  for(unsigned ix=0;ix<3;++ix) {
    if(ix==1&&photon) {
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=0.;
    }
    else {
      epsdot=_vectors[1][ix]*inpart.momentum();
      for(unsigned int iy=0;iy<3;++iy)
	(*ME())(iy,ix,0)=Complex(pre*_vectors[0][iy].dot(p0dotpv*_vectors[1][ix]
							 -epsdot*decay[0]->momentum()));
    }
  }
  // test of the matrix element
//   double me = newME.contract(rhoin).real();
//   Energy pcm=Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					   decay[1]->mass());
//   double test = sqr(_coupling[imode()])/3.*(2.*sqr(pcm)+3.*sqr(decay[0]->mass()));
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return ME()->contract(_rho).real();
}

bool PVectorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,
						     int & mecode,
						     double & coupling) const {
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0); bool order(false);
  int imode(-1);
  do {
    if(id==_incoming[ix]) {
      if(id1   ==_outgoingV[ix]&&id2   ==_outgoingP[ix]) {
	imode=ix;
	order=true;
      }
      if(id2   ==_outgoingV[ix]&&id1   ==_outgoingP[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoingV[ix]&&id2bar==_outgoingP[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoingV[ix]&&id1bar==_outgoingP[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling = _coupling[imode]*dm.parent()->mass();  
  mecode = 4;
  return order;
}

void PVectorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " "
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":OutgoingVector " << ix << " "
	     << _outgoingV[ix] << "\n";
      output << "newdef " << name() << ":OutgoingPScalar " << ix << " "
	     << _outgoingP[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " "
	     << _coupling[ix]*GeV << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming "  << ix << " "
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":OutgoingVector " << ix << " "
	     << _outgoingV[ix] << "\n";
      output << "insert " << name() << ":OutgoingPScalar " << ix << " "
	     << _outgoingP[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " "
	     << _coupling[ix]*GeV << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
