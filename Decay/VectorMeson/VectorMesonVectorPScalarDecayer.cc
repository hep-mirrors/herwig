// -*- C++ -*-
//
// VectorMesonVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorPScalarDecayer class.
//

#include "VectorMesonVectorPScalarDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix) {
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMesonVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingV.size()||isize!=_outgoingP.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "VectorMesonVectorPScalarDecayer" << Exception::runerror;
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

VectorMesonVectorPScalarDecayer::VectorMesonVectorPScalarDecayer() 
  : _coupling(73), _incoming(73), _outgoingV(73), _outgoingP(73), _maxweight(73) {
  // intermediates
  generateIntermediates(false);
  // rho -> gamma pi modes
  _incoming[0] =  113; _outgoingV[0] =  22; _outgoingP[0] =  111; 
  _coupling[0] = 0.2527/GeV; _maxweight[0] = 1.7; 
  _incoming[1] =  213; _outgoingV[1] =  22; _outgoingP[1] =  211; 
  _coupling[1] = 0.2210/GeV; _maxweight[1] = 1.7; 
  // rho  -> gamma eta mode
  _incoming[2] =  113; _outgoingV[2] =  22; _outgoingP[2] =  221; 
  _coupling[2] = 0.492/GeV; _maxweight[2] = 1.7; 
  // omega -> gamma pi 
  _incoming[3] =  223; _outgoingV[3] =  22; _outgoingP[3] =  111; 
  _coupling[3] = 0.7279465/GeV; _maxweight[3] = 1.7; 
  // omega -> gamma eta
  _incoming[4] =  223; _outgoingV[4] =  22; _outgoingP[4] =  221; 
  _coupling[4] = 0.143/GeV; _maxweight[4] = 1.7; 
  // phi -> gamma pi 
  _incoming[5] =  333; _outgoingV[5] =  22; _outgoingP[5] =  111; 
  _coupling[5] = 0.0397/GeV; _maxweight[5] = 1.7; 
  // phi -> gamma eta
  _incoming[6] =  333; _outgoingV[6] =  22; _outgoingP[6] =  221; 
  _coupling[6] = 0.212/GeV; _maxweight[6] = 1.7; 
  // phi -> gamma eta'
  _incoming[7] =  333; _outgoingV[7] =  22; _outgoingP[7] =  331; 
  _coupling[7] = 0.219/GeV; _maxweight[7] = 1.8; 
  // phi -> omega pi
  _incoming[8] =  333; _outgoingV[8] = 223; _outgoingP[8] =  111; 
  _coupling[8] = 0.0417/GeV; _maxweight[8] = 1.9; 
  // phi' -> K* K
  _incoming[9] =  100333; _outgoingV[9] =  323; _outgoingP[9] = -321; 
  _coupling[9] = 3.934/GeV; _maxweight[9] = 5.; 
  _incoming[10] =  100333; _outgoingV[10] =  313; _outgoingP[10] = -311; 
  _coupling[10] = 4.011/GeV; _maxweight[10] = 5.; 
  // K* -> gamma K 
  _incoming[11] =  313; _outgoingV[11] =  22; _outgoingP[11] =  311; 
  _coupling[11] = 0.384/GeV; _maxweight[11] = 1.7; 
  _incoming[12] =  323; _outgoingV[12] =  22; _outgoingP[12] =  321; 
  _coupling[12] = 0.253/GeV; _maxweight[12] = 1.7; 
  // d* decay
  _incoming[13] =  423; _outgoingV[13] = 22; _outgoingP[13] =  421; 
  _coupling[13] = 0.616/GeV; _maxweight[13] = 1.7; 
  _incoming[14] =  413; _outgoingV[14] = 22; _outgoingP[14] =  411; 
  _coupling[14] = 0.152/GeV; _maxweight[14] = 1.7; 
  // D_s* decays
  _incoming[15] =  433; _outgoingV[15] = 22; _outgoingP[15] =  431; 
  _coupling[15] = 0.764/GeV; _maxweight[15] = 1.7; 
  // B_s* decays
  _incoming[16] =  533; _outgoingV[16] = 22; _outgoingP[16] =  531; 
  _coupling[16] = 0.248/GeV; _maxweight[16] = 1.7; 
  // B_c* decays
  _incoming[17] =  543; _outgoingV[17] = 22; _outgoingP[17] =  541; 
  _coupling[17] = 0.266/GeV; _maxweight[17] = 1.7; 
  // B* decay
  _incoming[18] =  523; _outgoingV[18] = 22; _outgoingP[18] =  521; 
  _coupling[18] = 0.553/GeV; _maxweight[18] = 1.7; 
  _incoming[19] =  513; _outgoingV[19] = 22; _outgoingP[19] =  511; 
  _coupling[19] = 0.310/GeV; _maxweight[19] = 1.7; 
  // rho'' eta rho
  _incoming[20] =  30113; _outgoingV[20] =  113; _outgoingP[20] = 221; 
  _coupling[20] = 2.663/GeV; _maxweight[20] = 4.; 
  _incoming[21] =  30213; _outgoingV[21] =  213; _outgoingP[21] = 221; 
  _coupling[21] = 2.663/GeV; _maxweight[21] = 4.; 
  // rho '' K* K
  _incoming[22] =  30113; _outgoingV[22] =  323; _outgoingP[22] = -321; 
  _coupling[22] = 0.894/GeV; _maxweight[22] = 5.; 
  _incoming[23] =  30113; _outgoingV[23] =  313; _outgoingP[23] = -311; 
  _coupling[23] = 0.908/GeV; _maxweight[23] = 5.; 
  _incoming[24] =  30213; _outgoingV[24] =  323; _outgoingP[24] = -311; 
  _coupling[24] = 1.265/GeV; _maxweight[24] = 5.; 
  _incoming[25] =  30213; _outgoingV[25] = -313; _outgoingP[25] =  321; 
  _coupling[25] = 1.273/GeV; _maxweight[25] = 5.; 
  // omega'' rho pi
  _incoming[26] =  30223; _outgoingV[26] =  213; _outgoingP[26] = -211; 
  _coupling[26] = 2.996/GeV; _maxweight[26] = 3.5; 
  _incoming[27] =  30223; _outgoingV[27] =  113; _outgoingP[27] =  111; 
  _coupling[27] = 2.996/GeV; _maxweight[27] = 3.5; 
  // omega' rho pi
  _incoming[28] =  100223; _outgoingV[28] =  213; _outgoingP[28] = -211; 
  _coupling[28] = 4.507/GeV; _maxweight[28] = 4.; 
  _incoming[29] =  100223; _outgoingV[29] =  113; _outgoingP[29] =  111; 
  _coupling[29] = 4.507/GeV; _maxweight[29] = 4.; 
  // K*''->K* pi decays
  _incoming[30] =  30313; _outgoingV[30] =  323; _outgoingP[30] = -211; 
  _coupling[30] = 3.36/GeV;  _maxweight[30] = 4.; 
  _incoming[31] =  30313; _outgoingV[31] =  313; _outgoingP[31] =  111; 
  _coupling[31] = 2.38/GeV;  _maxweight[31] = 4.; 
  _incoming[32] =  30323; _outgoingV[32] =  313; _outgoingP[32] =  211; 
  _coupling[32] = 3.36/GeV;  _maxweight[32] = 4.; 
  _incoming[33] =  30323; _outgoingV[33] =  323; _outgoingP[33] =  111; 
  _coupling[33] = 2.38/GeV;  _maxweight[33] = 4.; 
  // K*''->K rho decays
  _incoming[34] =  30313; _outgoingP[34] =  321; _outgoingV[34] = -213; 
  _coupling[34] = 4.159/GeV;  _maxweight[34] = 3.; 
  _incoming[35] =  30313; _outgoingP[35] =  311; _outgoingV[35] =  113; 
  _coupling[35] = 2.939/GeV;  _maxweight[35] = 3.; 
  _incoming[36] =  30323; _outgoingP[36] =  311; _outgoingV[36] =  213; 
  _coupling[36] = 4.159/GeV;  _maxweight[36] = 3.; 
  _incoming[37] =  30323; _outgoingP[37] =  321; _outgoingV[37] =  113; 
  _coupling[37] = 2.939/GeV;  _maxweight[37] = 3.; 
  // K*' decays
  _incoming[38] =  100313; _outgoingV[38] =  323; _outgoingP[38] = -211; 
  _coupling[38] = 9.469/GeV;  _maxweight[38] = 6.; 
  _incoming[39] =  100313; _outgoingV[39] =  313; _outgoingP[39] =  111; 
  _coupling[39] = 6.781/GeV;  _maxweight[39] = 6.; 
  _incoming[40] =  100323; _outgoingV[40] =  313; _outgoingP[40] =  211; 
  _coupling[40] = 9.469/GeV;  _maxweight[40] = 6.; 
  _incoming[41] =  100323; _outgoingV[41] =  323; _outgoingP[41] =  111; 
  _coupling[41] = 6.781/GeV;  _maxweight[41] = 6.; 
  // J/psi -> gamma eta_c decay
  _incoming[42] = 443; _outgoingV[42] = 22; _outgoingP[42] = 441; 
  _coupling[42] = 0.149/GeV; _maxweight[42] = 30.; 
  // J/psi -> gamma eta' decay
  _incoming[43] = 443; _outgoingV[43] = 22; _outgoingP[43] = 331; 
  _coupling[43] = 0.00250/GeV; _maxweight[43] = 1.7; 
  // J/psi -> rho pi decay
  _incoming[44] = 443; _outgoingV[44] = 213; _outgoingP[44] = -211; 
  _coupling[44] = 0.00274/GeV; _maxweight[44] = 3.3; 
  _incoming[45] = 443; _outgoingV[45] = 113; _outgoingP[45] = 111; 
  _coupling[45] = 0.00274/GeV; _maxweight[45] = 3.3; 
  // J/psi -> K* K decay
  _incoming[46] = 443; _outgoingV[46] = 323; _outgoingP[46] = -321; 
  _coupling[46] = 0.00180/GeV; _maxweight[46] = 6.; 
  _incoming[47] = 443; _outgoingV[47] = 313; _outgoingP[47] = -311; 
  _coupling[47] = 0.00180/GeV; _maxweight[47] = 6.; 
  // J/psi -> omega eta decay
  _incoming[48] = 443; _outgoingV[48] = 223; _outgoingP[48] = 221; 
  _coupling[48] = 0.00154/GeV; _maxweight[48] = 6.; 
  // J/psi -> gamma eta decay
  _incoming[49] = 443; _outgoingV[49] = 22; _outgoingP[49] = 221; 
  _coupling[49] = 0.00103/GeV; _maxweight[49] = 1.7; 
  // J/psi -> phi eta decay
  _incoming[50] = 443; _outgoingV[50] = 333; _outgoingP[50] = 221; 
  _coupling[50] = 0.00110/GeV; _maxweight[50] = 5.5; 
  // J/psi -> phi eta' decay
  _incoming[51] = 443; _outgoingV[51] = 333; _outgoingP[51] = 331; 
  _coupling[51] = 0.00085/GeV; _maxweight[51] = 5.5; 
  // J/psi -> omega pi0 decay
  _incoming[52] = 443; _outgoingV[52] = 223; _outgoingP[52] = 111; 
  _coupling[52] = 0.00073/GeV; _maxweight[52] = 6.; 
  // J/psi -> rho0 eta decay
  _incoming[53] = 443; _outgoingV[53] = 113; _outgoingP[53] = 221; 
  _coupling[53] = 0.00054/GeV; _maxweight[53] = 3.; 
  // J/psi -> rho0 eta' decay
  _incoming[54] = 443; _outgoingV[54] = 113; _outgoingP[54] = 331; 
  _coupling[54] = 0.00045/GeV; _maxweight[54] = 3.; 
  // J/psi -> omega eta' decay
  _incoming[55] = 443; _outgoingV[55] = 223; _outgoingP[55] = 331; 
  _coupling[55] = 0.00058/GeV; _maxweight[55] = 6.; 
  // J/psi -> gamma pi0 decay
  _incoming[56] = 443; _outgoingV[56] = 22; _outgoingP[56] = 111; 
  _coupling[56] = 0.000177/GeV; _maxweight[56] = 1.7; 
  // psi(2s)-> j/psi eta decay
  _incoming[57] = 100443; _outgoingV[57] = 443; _outgoingP[57] = 221; 
  _coupling[57] = 0.230/GeV; _maxweight[57] = 1.7; 
  // psi(2s)-> gamma eta_c decay
  _incoming[58] = 100443; _outgoingV[58] = 22; _outgoingP[58] = 441; 
  _coupling[58] = 0.0114/GeV; _maxweight[58] = 3.4; 
  // psi(2s)-> gamma eta' decay
  _incoming[59] = 100443; _outgoingV[59] = 22; _outgoingP[59] = 331; 
  _coupling[59] = 0.00062/GeV; _maxweight[59] = 1.7; 
  // psi(2s)-> J/psi pi0 decay
  _incoming[60] = 100443; _outgoingV[60] = 443; _outgoingP[60] = 111; 
  _coupling[60] = 0.0106/GeV; _maxweight[60] = 1.7; 
  // psi(1d)-> gamma eta_c decay
  _incoming[61] = 30443; _outgoingV[61] = 443; _outgoingP[61] = 221; 
  _coupling[61] = 0.135/GeV; _maxweight[61] = 1.7;
  _incoming[62] = 30443; _outgoingV[62] = 333; _outgoingP[62] = 221; 
  _coupling[62] = 0.0076/GeV; _maxweight[62] = 5.5; 
  // psi(2s) K* K
  _incoming[63] = 100443; _outgoingV[63] = 323; _outgoingP[63] = -321; 
  _coupling[63] = 0.00021/GeV; _maxweight[63] = 6.5; 
  _incoming[64] = 100443; _outgoingV[64] = 313; _outgoingP[64] = -311; 
  _coupling[64] = 0.00054/GeV; _maxweight[64] = 6.5; 
  // psi(2s) -> phi eta decay
  _incoming[65] = 100443; _outgoingV[65] = 333; _outgoingP[65] = 221; 
  _coupling[65] = 0.00029/GeV; _maxweight[65] = 5.5; 
  // psi(2s) -> phi eta' decay
  _incoming[66] = 100443; _outgoingV[66] = 333; _outgoingP[66] = 331; 
  _coupling[66] = 0.00033/GeV; _maxweight[66] = 5.5; 
  // psi(2s) -> rho pi decay
  _incoming[67] = 100443; _outgoingV[67] = 213; _outgoingP[67] = -211; 
  _coupling[67] = 0.00017/GeV; _maxweight[67] = 3.5; 
  _incoming[68] = 100443; _outgoingV[68] = 113; _outgoingP[68] = 111; 
  _coupling[68] = 0.00017/GeV; _maxweight[68] = 3.5; 
  // psi(2s) -> omega eta' decay
  _incoming[69] = 100443; _outgoingV[69] = 223; _outgoingP[69] = 331; 
  _coupling[69] = 0.00032/GeV; _maxweight[69] = 6.;
  // psi(2s) -> rho0 eta decay
  _incoming[70] = 100443; _outgoingV[70] = 113; _outgoingP[70] = 221; 
  _coupling[70] = 0.00025/GeV; _maxweight[70] = 3.5; 
  // psi(2s) -> omega pi0 decay
  _incoming[71] = 100443; _outgoingV[71] = 223; _outgoingP[71] = 111; 
  _coupling[71] = 0.00022/GeV; _maxweight[71] = 6.; 
  // psi(2s) -> rho0 eta' decay
  _incoming[72] = 100443; _outgoingV[72] = 113; _outgoingP[72] = 331; 
  _coupling[72] = 0.00026/GeV; _maxweight[72] = 3.5; 
  // initial size of the vectors for the database output
  _initsize=_incoming.size();
}

int VectorMesonVectorPScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(idbar==_incoming[ix]) {
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


void VectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << _incoming << _outgoingV << _outgoingP << _maxweight << ounit(_coupling,1/MeV);
}

void VectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> iunit(_coupling,1/MeV);
}

ClassDescription<VectorMesonVectorPScalarDecayer> 
VectorMesonVectorPScalarDecayer::initVectorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void VectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorPScalarDecayer> documentation
    ("The VectorMesonVectorPScalarDecayer class is designed for the "
     "decay of a vector meson to another vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1 particle",
     &VectorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceOutcomingPScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &VectorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonVectorPScalarDecayer::_coupling,
     1/MeV, 0, ZERO, -10000000/MeV, 10000000/MeV, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, -10000000, 10000000, false, false, true);

}

double VectorMesonVectorPScalarDecayer::me2(const int,
					    const Particle & inpart,
					    const ParticleVector & decay,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
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
  LorentzPolarizationVector vtemp;
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1&&photon) {
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=0.;
    }
    else {
      vtemp=_coupling[imode()]/inpart.mass()*
	epsilon(inpart.momentum(),_vectors[1][ix],decay[0]->momentum());
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=_vectors[0][iy].dot(vtemp);
    }
  }
  // test of the matrix element
//   double me = newME.contract(rhoin).real();
//   Energy pcm=Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					   decay[1]->mass());
//   double test = sqr(_coupling[imode()]*pcm)*2./3.;
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return ME()->contract(_rho).real();
}

bool VectorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						    double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
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
  mecode = 1;
  return order;
}

// output the setup info for the particle database
void VectorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
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
	     << _coupling[ix]*MeV << "\n";
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
	     << _coupling[ix]*MeV << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
