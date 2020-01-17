// -*- C++ -*-
//
// ScalarScalarScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarScalarScalarDecayer class.
//

#include "ScalarScalarScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarScalarScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix) )_maxweight[ix] = mode(ix)->maxWeight();
  }
}

ScalarScalarScalarDecayer::ScalarScalarScalarDecayer() 
  : _incoming(78), _outgoing1(78), _outgoing2(78), 
    _coupling(78), _maxweight(78) {
  // f_0(980) to pi pi
  _incoming[0] = 9010221; _outgoing1[0] = 111; _outgoing2[0] =  111; 
  _coupling[0] = 1.66*GeV; _maxweight[0] = 1.05; 
  _incoming[1] = 9010221; _outgoing1[1] = 211; _outgoing2[1] = -211; 
  _coupling[1] = 2.35*GeV; _maxweight[1] = 1.05; 
  // f_0(980) to K K 
  _incoming[2] = 9010221; _outgoing1[2] = 321; _outgoing2[2] = -321; 
  _coupling[2] = 1.02*GeV; _maxweight[2] = 1.05; 
  _incoming[3] = 9010221; _outgoing1[3] = 311; _outgoing2[3] = -311; 
  _coupling[3] = 1.02*GeV; _maxweight[3] = 1.05; 
  // f_0(1370) to pi pi
  _incoming[4] = 10221; _outgoing1[4] = 111; _outgoing2[4] = 111; 
  _coupling[4] = 0.745*GeV; _maxweight[4] = 1.05; 
  _incoming[5] = 10221; _outgoing1[5] = 211; _outgoing2[5] = -211; 
  _coupling[5] = 1.054*GeV; _maxweight[5] = 1.05; 
  // f_0(1370) to pi' pi
  _incoming[6] = 10221; _outgoing1[6] = 100111; _outgoing2[6] = 111; 
  _coupling[6] = 5.027*GeV; _maxweight[6] = 2.1; 
  _incoming[7] = 10221; _outgoing1[7] = 100211; _outgoing2[7] = -211; 
  _coupling[7] = 5.027*GeV; _maxweight[7] = 2.1; 
  // f_0(1370) to K K 
  _incoming[8] = 10221; _outgoing1[8] = 321; _outgoing2[8] = -321; 
  _coupling[8] = 0.886*GeV; _maxweight[8] = 1.1; 
  _incoming[9] = 10221; _outgoing1[9] = 311; _outgoing2[9] = -311; 
  _coupling[9] = 0.886*GeV; _maxweight[9] = 1.1; 
  // f_0(1710) to pi pi
  _incoming[10] = 10331; _outgoing1[10] = 111; _outgoing2[10] = 111; 
  _coupling[10] = 0.503*GeV; _maxweight[10] = 1.05; 
  _incoming[11] = 10331; _outgoing1[11] = 211; _outgoing2[11] = -211; 
  _coupling[11] = 0.711*GeV; _maxweight[11] = 1.05; 
  // f_0(1710) to K K 
  _incoming[12] = 10331; _outgoing1[12] = 321; _outgoing2[12] = -321; 
  _coupling[12] = 2.096*GeV; _maxweight[12] = 1.05; 
  _incoming[13] = 10331; _outgoing1[13] = 311; _outgoing2[13] = -311; 
  _coupling[13] = 2.096*GeV; _maxweight[13] = 1.05; 
  // sigma to pi pi
  _incoming[14] = 9000221; _outgoing1[14] = 111; _outgoing2[14] = 111; 
  _coupling[14] = 3.654*GeV; _maxweight[14] = 1.05; 
  _incoming[15] = 9000221; _outgoing1[15] = 211; _outgoing2[15] = -211; 
  _coupling[15] = 5.178*GeV; _maxweight[15] = 1.05; 
  // a_0 to eta pi
  _incoming[16] =  9000111; _outgoing1[16] = 221; _outgoing2[16] =  111; 
  _coupling[16] = 3.33*GeV; _maxweight[16] = 1.1; 
  _incoming[17] =  9000211; _outgoing1[17] = 221; _outgoing2[17] =  211; 
  _coupling[17] = 3.33*GeV; _maxweight[17] = 1.1; 
  // a_0 to K K
  _incoming[18] =  9000111; _outgoing1[18] = 321; _outgoing2[18] = -321; 
  _coupling[18] = 2.54*GeV; _maxweight[18] = 1.05; 
  _incoming[19] =  9000111; _outgoing1[19] = 311; _outgoing2[19] = -311; 
  _coupling[19] = 2.54*GeV; _maxweight[19] = 1.05; 
  _incoming[20] =  9000211; _outgoing1[20] = 321; _outgoing2[20] = -311; 
  _coupling[20] = 3.59*GeV; _maxweight[20] = 1.05; 
  // a'_0 to eta pi
  _incoming[21] =  10111; _outgoing1[21] = 221; _outgoing2[21] =  111; 
  _coupling[21] = 1.357*GeV; _maxweight[21] = 1.05; 
  _incoming[22] =  10211; _outgoing1[22] = 221; _outgoing2[22] =  211; 
  _coupling[22] = 1.357*GeV; _maxweight[22] = 1.05; 
  // a'_0 to eta' pi
  _incoming[23] =  10111; _outgoing1[23] = 331; _outgoing2[23] =  111; 
  _coupling[23] = 0.995*GeV; _maxweight[23] = 1.6; 
  _incoming[24] =  10211; _outgoing1[24] = 331; _outgoing2[24] =  211; 
  _coupling[24] = 0.995*GeV; _maxweight[24] = 1.6; 
  // a'_0 to K K
  _incoming[25] =  10111; _outgoing1[25] = 321; _outgoing2[25] = -321; 
  _coupling[25] = 0.950*GeV; _maxweight[25] = 1.05; 
  _incoming[26] =  10111; _outgoing1[26] = 311; _outgoing2[26] = -311; 
  _coupling[26] = 0.950*GeV; _maxweight[26] = 1.05; 
  _incoming[27] =  10211; _outgoing1[27] = 321; _outgoing2[27] = -311; 
  _coupling[27] = 1.344*GeV; _maxweight[27] = 1.05; 
  // f_0(1370) to eta eta  
  _incoming[28] =  10221; _outgoing1[28] = 221; _outgoing2[28] = 221; 
  _coupling[28] = 0.235*GeV; _maxweight[28] = 1.1; 
  // f_0(1710) to eta eta  
  _incoming[29] =  10331; _outgoing1[29] = 221; _outgoing2[29] = 221; 
  _coupling[29] = 2.189*GeV; _maxweight[29] = 1.2; 
  // f_0(1370) to sigma sigma  
  _incoming[30] =  10221; _outgoing1[30] = 9000221; _outgoing2[30] = 9000221; 
  _coupling[30] = 21.46*GeV; _maxweight[30] = 7.; 
  // K_0* to K pi
  _incoming[31] =  10311; _outgoing1[31] =  311; _outgoing2[31] =  111; 
  _coupling[31] = 2.837*GeV; _maxweight[31] = 1.05; 
  _incoming[32] =  10311; _outgoing1[32] =  321; _outgoing2[32] = -211; 
  _coupling[32] = 4.000*GeV; _maxweight[32] = 1.05; 
  _incoming[33] =  10321; _outgoing1[33] =  321; _outgoing2[33] =  111; 
  _coupling[33] = 2.837*GeV; _maxweight[33] = 1.05; 
  _incoming[34] =  10321; _outgoing1[34] =  311; _outgoing2[34] =  211; 
  _coupling[34] = 4.000*GeV; _maxweight[34] = 1.05; 
  // D_0* to D pi
  _incoming[35] =  10411; _outgoing1[35] =  411; _outgoing2[35] =  111; 
  _coupling[35] = 5.472*GeV; _maxweight[35] = 1.05; 
  _incoming[36] =  10411; _outgoing1[36] =  421; _outgoing2[36] =  211; 
  _coupling[36] = 7.714*GeV; _maxweight[36] = 1.05; 
  _incoming[37] =  10421; _outgoing1[37] =  421; _outgoing2[37] =  111; 
  _coupling[37] = 5.447*GeV; _maxweight[37] = 1.05; 
  _incoming[38] =  10421; _outgoing1[38] =  411; _outgoing2[38] = -211; 
  _coupling[38] = 7.818*GeV; _maxweight[38] = 1.05; 
  // B_0* to B pi
  _incoming[39] =  10511; _outgoing1[39] =  511; _outgoing2[39] =  111; 
  _coupling[39] = 9.698*GeV; _maxweight[39] = 1.05; 
  _incoming[40] =  10511; _outgoing1[40] =  521; _outgoing2[40] = -211; 
  _coupling[40] = 13.71*GeV; _maxweight[40] = 1.05; 
  _incoming[41] =  10521; _outgoing1[41] =  521; _outgoing2[41] =  111; 
  _coupling[41] = 9.698*GeV; _maxweight[41] = 1.05; 
  _incoming[42] =  10521; _outgoing1[42] =  511; _outgoing2[42] =  211; 
  _coupling[42] = 13.71*GeV; _maxweight[42] = 1.05; 
  // K' to K_0* pi
  _incoming[43] =  100311; _outgoing1[43] =  10311; _outgoing2[43] =  111; 
  _coupling[43] = 6.595*GeV; _maxweight[43] = 2.; 
  _incoming[44] =  100311; _outgoing1[44] =  10321; _outgoing2[44] = -211; 
  _coupling[44] = 9.445*GeV; _maxweight[44] = 2.; 
  _incoming[45] =  100321; _outgoing1[45] =  10321; _outgoing2[45] =  111; 
  _coupling[45] = 6.595*GeV; _maxweight[45] = 2.; 
  _incoming[46] =  100321; _outgoing1[46] =  10311; _outgoing2[46] =  211; 
  _coupling[46] = 9.445*GeV; _maxweight[46] = 2.; 
  // D_s0* to D_s pi
  _incoming[47] =  10431; _outgoing1[47] =  431; _outgoing2[47] =  111; 
  _coupling[47] = 0.103*GeV; _maxweight[47] = 1.05; 
  // B_s0* to B_s pi
  _incoming[48] =  10531; _outgoing1[48] =  531; _outgoing2[48] =  111; 
  _coupling[48] = 8.314*GeV; _maxweight[48] = 1.05; 
  // eta'' to a_0 pi
  _incoming[49] = 100221; _outgoing1[49] = 9000111; _outgoing2[49] = 111; 
  _coupling[49] = 2.057*GeV; _maxweight[49] = 2.; 
  _incoming[50] = 100221; _outgoing1[50] = 9000211; _outgoing2[50] = -211; 
  _coupling[50] = 2.057*GeV; _maxweight[50] = 2.; 
  // eta''' to a_0 pi
  _incoming[51] = 9020221; _outgoing1[51] = 9000111; _outgoing2[51] = 111; 
  _coupling[51] = 1.470*GeV; _maxweight[51] = 2.; 
  _incoming[52] = 9020221; _outgoing1[52] = 9000211; _outgoing2[52] = -211; 
  _coupling[52] = 1.470*GeV; _maxweight[52] = 2.; 
  // eta''' to sigma eta
  _incoming[53] = 9020221; _outgoing1[53] = 221; _outgoing2[53] = 9000221; 
  _coupling[53] = 4.051*GeV; _maxweight[53] = 2.; 
  // eta'' to sigma eta
  _incoming[54] = 100221; _outgoing1[54] = 9000221; _outgoing2[54] = 221; 
  _coupling[54] = 4.316*GeV; _maxweight[54] = 2.; 
  // chi_0c decays to K K 
  _incoming[55] = 10441; _outgoing1[55] = 321; _outgoing2[55] = -321; 
  _coupling[55] = 0.104*GeV; _maxweight[55] = 1.05; 
  _incoming[56] = 10441; _outgoing1[56] = 311; _outgoing2[56] = -311; 
  _coupling[56] = 0.104*GeV; _maxweight[56] = 1.05; 
  // chi_0c decays to pi pi
  _incoming[57] = 10441; _outgoing1[57] = 211; _outgoing2[57] = -211; 
  _coupling[57] = 0.093*GeV; _maxweight[57] = 1.05; 
  _incoming[58] = 10441; _outgoing1[58] = 111; _outgoing2[58] = 111; 
  _coupling[58] = 0.066*GeV; _maxweight[58] = 1.05; 
  // chi_0c decays to eta eta
  _incoming[59] = 10441; _outgoing1[59] = 221; _outgoing2[59] = 221; 
  _coupling[59] = 0.064*GeV; _maxweight[59] = 1.05; 
  // f_0(1500) to pipi
  _incoming[60] = 9030221; _outgoing1[60] = 211; _outgoing2[60] = -211; 
  _coupling[60] = 1.398*GeV; _maxweight[60] = 1.05; 
  _incoming[61] = 9030221; _outgoing1[61] = 111; _outgoing2[61] = 111; 
  _coupling[61] = 0.989*GeV; _maxweight[61] = 1.05; 
  // f_0(1500) to sigma sigma
  _incoming[62] = 9030221; _outgoing1[62] = 9000221; _outgoing2[62] = 9000221; 
  _coupling[62] = 6.079*GeV; _maxweight[62] = 7.; 
  // f_0(1500) to eta eta
  _incoming[63] = 9030221; _outgoing1[63] = 221; _outgoing2[63] = 221; 
  _coupling[63] = 0.809*GeV; _maxweight[63] = 1.05; 
  // f_0(1500) to eta eta'
  _incoming[64] = 9030221; _outgoing1[64] = 221; _outgoing2[64] = 331; 
  _coupling[64] = 2.844*GeV; _maxweight[64] = 4.5; 
  // f_0(1500) to K K
  _incoming[65] = 9030221; _outgoing1[65] = 321; _outgoing2[65] = -321; 
  _coupling[65] = 0.686*GeV; _maxweight[65] = 1.05; 
  _incoming[66] = 9030221; _outgoing1[66] = 311; _outgoing2[66] = -311; 
  _coupling[66] = 0.686*GeV; _maxweight[66] = 1.05; 
  // f_0(1500) to pi' pi
  _incoming[67] = 9030221; _outgoing1[67] = 100111; _outgoing2[67] = 111; 
  _coupling[67] = 2.615*GeV; _maxweight[67] = 2.1; 
  _incoming[68] = 9030221; _outgoing1[68] = 100211;_outgoing2[68] = -211; 
  _coupling[68] = 2.615*GeV; _maxweight[68] = 2.1;
  // kappa to K pi
  _incoming[69] =  9000311; _outgoing1[69] =  311; _outgoing2[69] =  111; 
  _coupling[69] = 3.834*GeV; _maxweight[69] = 1.05; 
  _incoming[70] =  9000311; _outgoing1[70] =  321; _outgoing2[70] = -211; 
  _coupling[70] = 5.406*GeV; _maxweight[70] = 1.05; 
  _incoming[71] =  9000321; _outgoing1[71] =  321; _outgoing2[71] =  111; 
  _coupling[71] = 3.834*GeV; _maxweight[71] = 1.05; 
  _incoming[72] =  9000321; _outgoing1[72] =  311; _outgoing2[72] =  211; 
  _coupling[72] = 5.406*GeV; _maxweight[72] = 1.05; 
  // chi_0c decays to K*_0 K*_0 
  _incoming[73] = 10441; _outgoing1[73] = 10321; _outgoing2[73] = -10321; 
  _coupling[73] = 0.104*GeV; _maxweight[73] = 1.05; 
  _incoming[74] = 10441; _outgoing1[74] = 10311; _outgoing2[74] = -10311; 
  _coupling[74] = 0.104*GeV; _maxweight[74] = 1.05; 
  // B*_s0 decays
  _incoming[75] = 10531; _outgoing1[75] = 511; _outgoing2[75] = -311; 
  _coupling[75] = 12.17*GeV; _maxweight[75] = 1.05; 
  _incoming[76] = 10531; _outgoing1[76] = 521; _outgoing2[76] = -321; 
  _coupling[76] = 12.17*GeV; _maxweight[76] = 1.05;  
  // chi_0c decays to f_0 f_0
  _incoming[77] = 10441; _outgoing1[77] = 9010221; _outgoing2[77] = 9010221; 
  _coupling[77] = 0.084*GeV; _maxweight[77] = 1.05; 

  // initial size
  _initsize = _coupling.size();
  // intermediates
  generateIntermediates(false);
}

void ScalarScalarScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(_coupling.size());
  if(isize!=_incoming.size()  || isize!=_outgoing1.size()||
     isize!=_outgoing2.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in ScalarScalarScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt;
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoing1[ix]);
    extpart[2] = getParticleData(_outgoing2[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int ScalarScalarScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id0(parent->id());
  int id0bar = parent->CC() ? parent->CC()->id() : id0;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  // loop over the modes and see if this is one of them
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id0   ==_incoming[ix]) {
      if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	 (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==_incoming[ix]&&imode<0) {
      if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	 (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void ScalarScalarScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,MeV) << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void ScalarScalarScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,MeV) >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarScalarScalarDecayer,DecayIntegrator>
describeHerwigScalarScalarScalarDecayer("Herwig::ScalarScalarScalarDecayer", "HwSMDecay.so");

void ScalarScalarScalarDecayer::Init() {

  static ClassDocumentation<ScalarScalarScalarDecayer> documentation
    ("The ScalarScalarScalarDecayer class is designed for the"
     " decay of a scalar meson to two scalar mesons including off-shell effects");

  static ParVector<ScalarScalarScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarScalarScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &ScalarScalarScalarDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &ScalarScalarScalarDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,Energy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarScalarScalarDecayer::_coupling,
     MeV, 0, ZERO, ZERO, 1000000.*MeV, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarScalarScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

double ScalarScalarScalarDecayer::me2(const int,
				      const Particle & inpart,
				      const ParticleVector & decay,
				      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return 0.;
  }
  double fact(_coupling[imode()]/inpart.mass());
  (*ME())(0,0,0) = fact;
  return sqr(fact);
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarScalarScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
  unsigned int ix(0); bool order(true);
  do {
    if(id   ==_incoming[ix]) {
      if(id1==_outgoing1[ix]&&id2==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoing1[ix]&&id1==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode]/dm.parent()->mass();
  itype = 6;
  return order;
}

// output the setup information for the particle database
void ScalarScalarScalarDecayer::dataBaseOutput(ofstream & output,
					       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "newdef " << name() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]/MeV << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "insert " << name() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]/MeV << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
