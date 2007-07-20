// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2FermionDecayer class.
//

#include "VectorMeson2FermionDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson2FermionDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoingf.size()||
     isize!=_outgoinga.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "FermionDecayer::doiin() " << Exception::runerror;
  // set up the integration channels
  vector<double> wgt(0);
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoingf[ix]);
    extpart[2]=getParticleData(_outgoinga[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

VectorMeson2FermionDecayer::VectorMeson2FermionDecayer() 
  : _coupling(42), _incoming(42), _outgoingf(42), _outgoinga(42), _maxweight(42) {
  // don't include intermediates
  generateIntermediates(false);
  // rho -> e+e-, mu+mu
  _incoming[0] = 113; _outgoingf[0] = 11; _outgoinga[0] = -11; 
  _coupling[0] = 18.524E-3;   _maxweight[0] = 2.; 
  _incoming[1] = 113; _outgoingf[1] = 13; _outgoinga[1] = -13; 
  _coupling[1] = 18.35E-3;    _maxweight[1] = 2.; 
  // omega -> e+e-, mu+mu-
  _incoming[2] = 223; _outgoingf[2] = 11; _outgoinga[2] = -11; 
  _coupling[2] = 5.387E-3;  _maxweight[2] = 2.; 
  _incoming[3] = 223; _outgoingf[3] = 13; _outgoinga[3] = -13; 
  _coupling[3] = 5.387E-3;  _maxweight[3] = 2.; 
  // phi -> e+e-, mu+mu-
  _incoming[4] = 333; _outgoingf[4] = 11; _outgoinga[4] = -11; 
  _coupling[4] = 6.852E-3;  _maxweight[4] = 2.; 
  _incoming[5] = 333; _outgoingf[5] = 13; _outgoinga[5] = -13; 
  _coupling[5] = 6.852E-3;  _maxweight[5] = 2.; 
  // psi(1d) to leptons
  _incoming[6] = 30443; _outgoingf[6] = 11; _outgoinga[6] = -11; 
  _coupling[6] = 1.611E-3;  _maxweight[6] = 2.; 
  _incoming[7] = 30443; _outgoingf[7] = 13; _outgoinga[7] = -13; 
  _coupling[7] = 1.611E-3;  _maxweight[7] = 2.; 
  _incoming[8] = 30443; _outgoingf[8] = 15; _outgoinga[8] = -15;  
  _coupling[8] = 1.611E-3;  _maxweight[8] = 2.; 
  // J/psi decay
  _incoming[9] = 443; _outgoingf[9] = 11; _outgoinga[9] = -11; 
  _coupling[9] = 8.088E-3;  _maxweight[9] = 2.; 
  _incoming[10] = 443; _outgoingf[10] = 13; _outgoinga[10] = -13; 
  _coupling[10] = 8.088E-3;  _maxweight[10] = 2.; 
  // psi2s to leptons
  _incoming[11] = 100443; _outgoingf[11] = 11; _outgoinga[11] = -11; 
  _coupling[11] = 4.645E-3;  _maxweight[11] = 2.; 
  _incoming[12] = 100443; _outgoingf[12] = 13; _outgoinga[12] = -13; 
  _coupling[12] = 4.645E-3;  _maxweight[12] = 2.; 
  _incoming[13] = 100443; _outgoingf[13] = 15; _outgoinga[13] = -15; 
  _coupling[13] = 4.645E-3;  _maxweight[13] = 2.; 
  // upsilon to leptons
  _incoming[14] = 553; _outgoingf[14] = 11; _outgoinga[14] = -11; 
  _coupling[14] = 2.290E-3;  _maxweight[14] = 2.; 
  _incoming[15] = 553; _outgoingf[15] = 13; _outgoinga[15] = -13; 
  _coupling[15] = 2.290E-3;  _maxweight[15] = 2.; 
  _incoming[16] = 553; _outgoingf[16] = 15; _outgoinga[16] = -15; 
  _coupling[16] = 2.290E-3;  _maxweight[16] = 2.; 
  // upsilon 2s to leptons
  _incoming[17] = 100553; _outgoingf[17] = 11; _outgoinga[17] = -11; 
  _coupling[17] = 1.466E-3;  _maxweight[17] = 2.; 
  _incoming[18] = 100553; _outgoingf[18] = 13; _outgoinga[18] = -13; 
  _coupling[18] = 1.466E-3;  _maxweight[18] = 2.; 
  _incoming[19] = 100553; _outgoingf[19] = 15; _outgoinga[19] = -15; 
  _coupling[19] = 1.466E-3;  _maxweight[19] = 2.; 
  // upsilon 3s to leptons
  _incoming[20] = 200553; _outgoingf[20] = 11; _outgoinga[20] = -11; 
  _coupling[20] = 1.316E-3;  _maxweight[20] = 2.; 
  _incoming[21] = 200553; _outgoingf[21] = 13; _outgoinga[21] = -13; 
  _coupling[21] = 1.316E-3;  _maxweight[21] = 2.; 
  _incoming[22] = 200553; _outgoingf[22] = 15; _outgoinga[22] = -15; 
  _coupling[22] = 1.316E-3;  _maxweight[22] = 2.; 
  // upsilon 4s to leptons
  _incoming[23] = 300553; _outgoingf[23] = 11; _outgoinga[23] = -11; 
  _coupling[23] = 1.411E-3;  _maxweight[23] = 2.; 
  _incoming[24] = 300553; _outgoingf[24] = 13; _outgoinga[24] = -13; 
  _coupling[24] = 1.411E-3;  _maxweight[24] = 2.; 
  _incoming[25] = 300553; _outgoingf[25] = 15; _outgoinga[25] = -15; 
  _coupling[25] = 1.411E-3;  _maxweight[25] = 2.; 
  // baryonic jpsi decays
  // to neutrons and proton
  _incoming[26] = 443; _outgoingf[26] = 2212; _outgoinga[26] = -2212; 
  _maxweight[26] = 2.; _coupling[26] = 1.581E-3; 
  _incoming[27] = 443; _outgoingf[27] = 2112; _outgoinga[27] = -2112; 
  _maxweight[27] = 2.; _coupling[27] = 1.581E-3; 
  // to sigma's
  _incoming[28] = 443; _outgoingf[28] = 3112; _outgoinga[28] = -3112; 
  _maxweight[28] = 2.; _coupling[28] = 1.307E-3; 
  _incoming[29] = 443; _outgoingf[29] = 3212; _outgoinga[29] = -3212; 
  _maxweight[29] = 2.; _coupling[29] = 1.307E-3; 
  _incoming[30] = 443; _outgoingf[30] = 3222; _outgoinga[30] = -3222; 
  _maxweight[30] = 2.; _coupling[30] = 1.307E-3; 
  // to Xi's 
  _incoming[31] = 443; _outgoingf[31] = 3322; _outgoinga[31] = -3322; 
  _maxweight[31] = 2.; _coupling[31] = 1.183E-3; 
  _incoming[32] = 443; _outgoingf[32] = 3312; _outgoinga[32] = -3312; 
  _maxweight[32] = 2.; _coupling[32] = 1.183E-3; 
  // to lambda
  _incoming[33] = 443; _outgoingf[33] = 3122; _outgoinga[33] = -3122; 
  _maxweight[33] = 2.; _coupling[33] = 1.284E-3; 
  // baryonic psi(2s) decays
  // to neutrons and protons
  _incoming[34] = 100443; _outgoingf[34] = 2212; _outgoinga[34] = -2212; 
  _maxweight[34] = 2.; _coupling[34] = 7.822E-4; 
  _incoming[35] = 100443; _outgoingf[35] = 2112; _outgoinga[35] = -2112; 
  _maxweight[35] = 2.; _coupling[35] = 7.822E-4; 
  // to sigma's
  _incoming[36] = 100443; _outgoingf[36] = 3112; _outgoinga[36] = -3112; 
  _maxweight[36] = 2.; _coupling[36] = 6.120E-4; 
  _incoming[37] = 100443; _outgoingf[37] = 3212; _outgoinga[37] = -3212; 
  _maxweight[37] = 2.; _coupling[37] = 6.120E-4; 
  _incoming[38] = 100443; _outgoingf[38] = 3222; _outgoinga[38] = -3222; 
  _maxweight[38] = 2.; _coupling[38] = 6.120E-4; 
  // to Xi's
  _incoming[39] = 100443; _outgoingf[39] = 3322; _outgoinga[39] = -3322; 
  _maxweight[39] = 2.; _coupling[39] = 5.544E-4; 
  _incoming[40] = 100443; _outgoingf[40] = 3312; _outgoinga[40] = -3312; 
  _maxweight[40] = 2.; _coupling[40] = 5.544E-4; 
  // to lambda
  _incoming[41] = 100443; _outgoingf[41] = 3122; _outgoinga[41] = -3122; 
  _maxweight[41] = 2.; _coupling[41] = 7.432E-4;   
  // set the initial size
  _initsize=_incoming.size();
}

int VectorMeson2FermionDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
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
    if(_incoming[ix]==id   ) {
      if((id1   ==_outgoingf[ix]&&id2   ==_outgoinga[ix])||
	 (id2   ==_outgoingf[ix]&&id1   ==_outgoinga[ix])) imode=ix;
    }
    if(_incoming[ix]==idbar) {
      if((id1bar==_outgoingf[ix]&&id2bar==_outgoinga[ix])||
	 (id2bar==_outgoingf[ix]&&id1bar==_outgoinga[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  return imode;
}

void VectorMeson2FermionDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingf << _outgoinga << _maxweight;
}

void VectorMeson2FermionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingf >> _outgoinga >> _maxweight;
}

ClassDescription<VectorMeson2FermionDecayer> 
VectorMeson2FermionDecayer::initVectorMeson2FermionDecayer;
// Definition of the static class description member.

void VectorMeson2FermionDecayer::Init() {

  static ClassDocumentation<VectorMeson2FermionDecayer> documentation
    ("The VectorMeson2FermionDecayer class is designed for the decay "
     "of vectro mesons to fermions. It is mainly used for the decay of vector mesons "
     "to electrons and muons.");

  static ParVector<VectorMeson2FermionDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson2FermionDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,int> interfaceOutcoming1
    ("OutgoingFermion",
     "The PDG code for the outgoing fermion",
     &VectorMeson2FermionDecayer::_outgoingf,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,int> interfaceOutcoming2
    ("OutgoingAntiFermion",
     "The PDG code for the second outgoing anti-fermion",
     &VectorMeson2FermionDecayer::_outgoinga,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMeson2FermionDecayer::_coupling,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMeson2FermionDecayer::_maxweight,
     0, 0, 0, -10000000, 10000000, false, false, true);

}

double VectorMeson2FermionDecayer::me2(bool vertex, const int,
				   const Particle & inpart,
				   const ParticleVector & decay) const {
  // wavefunctions of the decaying particle
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(_outgoingf[imode()]!=decay[iferm]->id()){iferm=1;ianti=0;}
  // construct the spin information objects for the  decay products
  vector<LorentzSpinor<SqrtEnergy> > wave;
  vector<LorentzSpinorBar<SqrtEnergy> > wavebar;
  SpinorBarWaveFunction(wavebar,decay[iferm],outgoing,true,vertex);
  SpinorWaveFunction(   wave   ,decay[ianti],outgoing,true,vertex);
  // prefactor
  InvEnergy pre(_coupling[imode()]/inpart.mass());
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  // now compute the currents
  LorentzPolarizationVector temp;
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      temp = pre*wave[ix].vectorCurrent(wavebar[iy]);
      for(unsigned int iz=0;iz<3;++iz) {
	if(iferm>ianti) newME(iz,ix,iy)=invec[iz].dot(temp);
	else            newME(iz,iy,ix)=invec[iz].dot(temp);
      }
    }
  }
  ME(newME);
  // test of the matrix element
//   double me = newME.contract(rhoin).real();
//   double test = 4.*sqr(_coupling[imode()])/3.*
//     (1.+2.*sqr(decay[0]->mass()/inpart.mass()));
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return newME.contract(rhoin).real();
}

bool VectorMeson2FermionDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==_incoming[ix]) {
      if(id1==_outgoingf[ix]&&id2==_outgoinga[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoingf[ix]&&id1==_outgoinga[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoingf[ix]&&id2bar==_outgoinga[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoingf[ix]&&id1bar==_outgoinga[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode];
  mecode=2;
  return order;
}

// output the setup information for the particle database
void VectorMeson2FermionDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "set " << fullName() << ":Incoming " << ix << " "
	     << _incoming[ix] << "\n";
      output << "set " << fullName() << ":OutgoingFermion " << ix << " "
	     << _outgoingf[ix] << "\n";
      output << "set " << fullName() << ":OutgoingAntiFermion "  << ix << " "
	     << _outgoinga[ix] << "\n";
      output << "set " << fullName() << ":Coupling " << ix << " "
	     << _coupling[ix] << "\n";
      output << "set " << fullName() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << fullName() << ":Incoming " << ix << " "
	     << _incoming[ix] << "\n";
      output << "insert " << fullName() << ":OutgoingFermion "  << ix << " "
	     << _outgoingf[ix] << "\n";
      output << "insert " << fullName() << ":OutgoingAntiFermion "  << ix << " "
	     << _outgoinga[ix] << "\n";
      output << "insert " << fullName() << ":Coupling " << ix << " "
	     << _coupling[ix] << "\n";
      output << "insert " << fullName() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
