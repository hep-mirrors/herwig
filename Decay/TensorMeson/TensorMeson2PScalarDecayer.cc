// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMeson2PScalarDecayer class.
//

#include "TensorMeson2PScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

TensorMeson2PScalarDecayer::TensorMeson2PScalarDecayer() 
  : _incoming(48), _outgoing1(48), _outgoing2(48), 
    _coupling(48), _maxweight(48) {
  // a_2 -> eta pi
  _incoming[0] = 115; _outgoing1[0] =  221; _outgoing2[0] = 111; 
  _coupling[0] = 10.90/GeV; _maxweight[0] = 2.02; 
  _incoming[1] = 215; _outgoing1[1] =  221; _outgoing2[1] = 211; 
  _coupling[1] = 10.90/GeV; _maxweight[1] = 2.00; 
  // a_2 -> eta' pi
  _incoming[2] = 115; _outgoing1[2] =  331; _outgoing2[2] = 111; 
  _coupling[2] = 9.92/GeV; _maxweight[2] = 2.2; 
  _incoming[3] = 215; _outgoing1[3] =  331; _outgoing2[3] = 211; 
  _coupling[3] = 9.92/GeV; _maxweight[3] = 2.1; 
  // a_2 -> K K
  _incoming[4] = 115; _outgoing1[4] =  311; _outgoing2[4] = -311; 
  _coupling[4] = 7.36/GeV; _maxweight[4] = 2.; 
  _incoming[5] = 115; _outgoing1[5] =  321; _outgoing2[5] = -321; 
  _coupling[5] = 7.36/GeV; _maxweight[5] = 2.02; 
  _incoming[6] = 215; _outgoing1[6] =  321; _outgoing2[6] = -311; 
  _coupling[6] = 10.41/GeV; _maxweight[6] = 2.; 
  // f_2 -> pi pi
  _incoming[7] = 225; _outgoing1[7] =  211; _outgoing2[7] = -211; 
  _coupling[7] = 18.73/GeV; _maxweight[7] = 2.; 
  _incoming[8] = 225; _outgoing1[8] =  111; _outgoing2[8] = 111; 
  _coupling[8] = 13.24/GeV; _maxweight[8] = 2.02; 
  // f_2 -> eta eta
  _incoming[9] = 225; _outgoing1[9] =  221; _outgoing2[9] = 221; 
  _coupling[9] = 8.854/GeV; _maxweight[9] = 2.15; 
  // f_2 -> KK
  _incoming[10] = 225; _outgoing1[10] =  321; _outgoing2[10] = -321; 
  _coupling[10] = 11.03/GeV; _maxweight[10] = 2.; 
  _incoming[11] = 225; _outgoing1[11] =  311; _outgoing2[11] = -311; 
  _coupling[11] = 11.38/GeV; _maxweight[11] = 2.; 
  // f'_2 -> KK
  _incoming[12] = 335; _outgoing1[12] =  321; _outgoing2[12] = -321; 
  _coupling[12] = 14.65/GeV; _maxweight[12] = 2.02; 
  _incoming[13] = 335; _outgoing1[13] =  311; _outgoing2[13] = -311; 
  _coupling[13] = 14.65/GeV; _maxweight[13] = 2.; 
  // f'_2 -> eta eta
  _incoming[14] = 335; _outgoing1[14] =  221; _outgoing2[14] = 221; 
  _coupling[14] = 9.15/GeV; _maxweight[14] = 2.02; 
  // f'_2 -> pi pi 
  _incoming[15] = 335; _outgoing1[15] =  211; _outgoing2[15] = -211; 
  _coupling[15] = 0.860/GeV; _maxweight[15] = 2.02; 
  _incoming[16] = 335; _outgoing1[16] =  111; _outgoing2[16] = 111; 
  _coupling[16] = 0.608/GeV; _maxweight[16] = 2.02; 
  // K_2 -> K eta
  _incoming[17] =  325; _outgoing1[17] =  321; _outgoing2[17] = 221; 
  _coupling[17] = 1.52/GeV; _maxweight[17] = 2.2; 
  _incoming[18] =  315; _outgoing1[18] =  311; _outgoing2[18] = 221; 
  _coupling[18] = 1.52/GeV; _maxweight[18] = 2.; 
  // K_2 -> K pi
  _incoming[19] =  325; _outgoing1[19] =  321; _outgoing2[19] =  111; 
  _coupling[19] = 8.30/GeV;_maxweight[19] = 2.; 
  _incoming[20] =  325; _outgoing1[20] =  311; _outgoing2[20] =  211; 
  _coupling[20] = 11.74/GeV; _maxweight[20] = 2.; 
  _incoming[21] =  315; _outgoing1[21] =  311; _outgoing2[21] =  111; 
  _coupling[21] = 8.68/GeV; _maxweight[21] = 2.; 
  _incoming[22] =  315; _outgoing1[22] =  321; _outgoing2[22] = -211; 
  _coupling[22] = 12.28/GeV; _maxweight[22] = 2.02; 
  // B_2 -> B pi
  _incoming[23] =  525; _outgoing1[23] =  521; _outgoing2[23] = 111; 
  _coupling[23] = 25.62/GeV; _maxweight[23] = 2.05; 
  _incoming[24] =  525; _outgoing1[24] =  511; _outgoing2[24] = 211; 
  _coupling[24] = 36.24/GeV; _maxweight[24] = 2.; 
  _incoming[25] =  515; _outgoing1[25] =  511; _outgoing2[25] = 111; 
  _coupling[25] = 25.62/GeV; _maxweight[25] = 2.05; 
  _incoming[26] =  515; _outgoing1[26] =  521; _outgoing2[26] = -211; 
  _coupling[26] = 36.24/GeV; _maxweight[26] = 2.02; 
  // D_2 -> D pi
  _incoming[27] =  425; _outgoing1[27] =  421; _outgoing2[27] = 111; 
  _coupling[27] = 13.17/GeV; _maxweight[27] = 2.02; 
  _incoming[28] =  425; _outgoing1[28] =  411; _outgoing2[28] = -211; 
  _coupling[28] = 18.62/GeV; _maxweight[28] = 2.; 
  _incoming[29] =  415; _outgoing1[29] =  411; _outgoing2[29] = 111; 
  _coupling[29] = 14.00/GeV; _maxweight[29] = 2.02; 
  _incoming[30] =  415; _outgoing1[30] =  421; _outgoing2[30] = 211; 
  _coupling[30] = 19.80/GeV; _maxweight[30] = 2.10; 
  // D_s2
  _incoming[31] =  435; _outgoing1[31] =  421; _outgoing2[31] =  321; 
  _coupling[31] = 23.39/GeV; _maxweight[31] = 2.02; 
  _incoming[32] =  435; _outgoing1[32] =  411; _outgoing2[32] =  311; 
  _coupling[32] = 23.39/GeV; _maxweight[32] = 2.02; 
  // B_s2
  _incoming[33] =  535; _outgoing1[33] =  521; _outgoing2[33] = -321; 
  _coupling[33] = 45.88/GeV; _maxweight[33] = 2.08; 
  _incoming[34] =  535; _outgoing1[34] =  511; _outgoing2[34] = -311; 
  _coupling[34] = 45.88/GeV; _maxweight[34] = 2.; 
  // chi_c2 to pi pi
  _incoming[35] =  445; _outgoing1[35] =  211; _outgoing2[35] = -211; 
  _coupling[35] = 0.0226/GeV; _maxweight[35] = 2.02; 
  _incoming[36] =  445; _outgoing1[36] =  111; _outgoing2[36] =  111; 
  _coupling[36] = 0.0159/GeV; _maxweight[36] = 1.62; 
  // chi_c2 to K K
  _incoming[37] =  445; _outgoing1[37] =  321; _outgoing2[37] = -321; 
  _coupling[37] = 0.0187/GeV; _maxweight[37] = 2.02; 
  _incoming[38] =  445; _outgoing1[38] =  311; _outgoing2[38] = -311; 
  _coupling[38] = 0.0187/GeV; _maxweight[38] = 2.02; 
  // f_2 to sigma sigma
  _incoming[39] = 225; _outgoing1[39] = 9000221; _outgoing2[39] = 9000221; 
  _coupling[39] = 99.27/GeV; _maxweight[39] = 30.; 
  // pi_2 to sigma pi
  _incoming[40] = 10115; _outgoing1[40] = 9000221; _outgoing2[40] = 111; 
  _coupling[40] = 99.27/GeV; _maxweight[40] = 30.; 
  _incoming[41] = 10215; _outgoing1[41] = 9000221; _outgoing2[41] = 211; 
  _coupling[41] = 99.27/GeV; _maxweight[41] = 30.; 
  // eta_2 to a_0 pi
  _incoming[42] = 10225; _outgoing1[42] = 9000111; _outgoing2[42] = 111; 
  _coupling[42] = 99.27/GeV; _maxweight[42] = 30.; 
  _incoming[43] = 10225; _outgoing1[43] = 9000211; _outgoing2[43] = -211; 
  _coupling[43] = 99.27/GeV; _maxweight[43] = 30.; 
  // eta'_2 to a_0 pi
  _incoming[44] = 10335; _outgoing1[44] = 9000111; _outgoing2[44] = 111; 
  _coupling[44] = 99.27/GeV; _maxweight[44] = 30.; 
  _incoming[45] = 10335; _outgoing1[45] = 9000211; _outgoing2[45] = -211; 
  _coupling[45] = 99.27/GeV; _maxweight[45] = 30.; 
  // chi_c2(2P) to D D
  _incoming[46] =  100445; _outgoing1[46] =  411; _outgoing2[46] = -411; 
  _coupling[46] = 0.0226/GeV; _maxweight[46] = 2.02; 
  _incoming[47] =  100445; _outgoing1[47] =  421; _outgoing2[47] = -421; 
  _coupling[47] = 0.0159/GeV; _maxweight[47] = 1.62; 
  // initial size of the vectors
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void TensorMeson2PScalarDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoing1.size()||isize!=_outgoing2.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters TensorMeson2PScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoing1[ix]);
    extpart[2]=getParticleData(_outgoing2[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int TensorMeson2PScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==_incoming[ix]) {
      if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	 (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])) imode=ix;
    }
    if(idbar==_incoming[ix]) {
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

void TensorMeson2PScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << _incoming << _outgoing1 << _outgoing2 << _maxweight 
     << ounit(_coupling,1/GeV);
}

void TensorMeson2PScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight 
     >> iunit(_coupling,1/GeV);
}

ClassDescription<TensorMeson2PScalarDecayer> 
TensorMeson2PScalarDecayer::initTensorMeson2PScalarDecayer;
// Definition of the static class description member.

void TensorMeson2PScalarDecayer::Init() {

  static ClassDocumentation<TensorMeson2PScalarDecayer> documentation
    ("The TensorMeson2PScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMeson2PScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMeson2PScalarDecayer::_coupling,
     1/GeV, 0, 0/GeV, 0/GeV, 100./GeV, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMeson2PScalarDecayer::_maxweight,
     0, 0, 0, 0., 100000., false, false, true);
}

// matrix elememt for the process
double TensorMeson2PScalarDecayer::me2(bool vertex, const int,
				       const Particle & inpart,
				       const ParticleVector & decay) const {
  vector<LorentzTensor<double> > inten;
  // wave functions etc for the incoming particle
  RhoDMatrix rhoin(PDT::Spin2);rhoin.average();
  TensorWaveFunction(inten,rhoin,const_ptr_cast<tPPtr>(&inpart),incoming,
		     true,false,vertex);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix) {
    // workaround for gcc 3.2.3 bug
    //ALB {ScalarWaveFunction(decay[ix],outgoing,true,vertex);}
    PPtr mytemp=decay[ix]; 
    ScalarWaveFunction(mytemp,outgoing,true,vertex);
  }
  // calculate the matrix element
  DecayMatrixElement newME(PDT::Spin2,PDT::Spin0,PDT::Spin0);
  for(unsigned int ix=0;ix<5;++ix) {
    newME(ix,0,0) = _coupling[imode()]/inpart.mass()*
      ((inten[ix]*decay[1]->momentum())*decay[0]->momentum());
  }
  ME(newME);
//   // test of the answer
//   double me = newME.contract(rhoin).real();
//   Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					     decay[1]->mass());
//   double test = Energy4(pow<4,1>(2*pcm))*sqr( _coupling[imode()]/inpart.mass())/120.;
//   cout << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " " 
//        << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return newME.contract(rhoin).real();
}

bool TensorMeson2PScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling=_coupling[imode]*dm.parent()->mass();
  mecode=7;
  return order;
}

void TensorMeson2PScalarDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "set " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "set " << fullName() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "set " << fullName() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "set " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV << "\n";
      output << "set " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << fullName() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "insert " << fullName() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "insert " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV << "\n";
      output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
