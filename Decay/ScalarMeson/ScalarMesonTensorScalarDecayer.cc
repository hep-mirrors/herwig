// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonTensorScalarDecayer class.
//

#include "ScalarMesonTensorScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

ScalarMesonTensorScalarDecayer::ScalarMesonTensorScalarDecayer() {
  _incoming.push_back( 411);_outgoingT.push_back(225);_outgoingS.push_back( 211);
  _coupling.push_back(7.24E-7/GeV);_maxweight.push_back(0.006);
  // initial size of the arrays
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

inline void ScalarMesonTensorScalarDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoingT.size()||
     isize!=_outgoingS.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in "
			  << "ScalarMesonTensorScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt;
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoingT[ix]);
    extpart[2] = getParticleData(_outgoingS[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int ScalarMesonTensorScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  // must be two outgoing particles
  if(dm.products().size()!=2) return -1;
  // ids of the particles
  int id0(dm.parent()->id());
  int id0bar = dm.parent()->CC() ? dm.parent()->CC()->id() : id0;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id0   ==_incoming[ix]) {
      if((id1   ==_outgoingT[ix]&&id2   ==_outgoingS[ix])||
	 (id2   ==_outgoingT[ix]&&id1   ==_outgoingS[ix])) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==_incoming[ix]&&imode<0) {
      if((id1bar==_outgoingT[ix]&&id2bar==_outgoingS[ix])||
	 (id2bar==_outgoingT[ix]&&id1bar==_outgoingS[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void ScalarMesonTensorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/GeV) << _incoming << _outgoingT << _outgoingS << _maxweight;
}

void ScalarMesonTensorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/GeV) >> _incoming >> _outgoingT >> _outgoingS >> _maxweight;
}

ClassDescription<ScalarMesonTensorScalarDecayer> 
ScalarMesonTensorScalarDecayer::initScalarMesonTensorScalarDecayer;
// Definition of the static class description member.

void ScalarMesonTensorScalarDecayer::Init() {

  static ClassDocumentation<ScalarMesonTensorScalarDecayer> documentation
    ("The ScalarMesonTensorScalarDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarMesonTensorScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceOutcomingT
    ("OutgoingTensor",
     "The PDG code for the outgoing tensor",
     &ScalarMesonTensorScalarDecayer::_outgoingT,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceOutcomingS
    ("OutgoingScalar",
     "The PDG code for the outgoing scalar",
     &ScalarMesonTensorScalarDecayer::_outgoingS,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarMesonTensorScalarDecayer::_coupling,
     1/GeV, 0, 0/GeV, 0./GeV, 100./GeV, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarMesonTensorScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);
}

double ScalarMesonTensorScalarDecayer::me2(bool vertex, const int,
					   const Particle & inpart,
					   const ParticleVector & decay) const {
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // set up the spin info for the outgoing particles
  vector<LorentzTensor<double> > twave;
  TensorWaveFunction(twave,decay[0],outgoing,true,false,vertex);
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr mytemp = decay[1];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);
  // calculate the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin2,PDT::Spin0);
  InvEnergy2 fact(_coupling[imode()]/inpart.mass());
  LorentzPolarizationVectorE vtemp;
  for(unsigned int ix=0;ix<5;++ix) {
    vtemp = twave[ix]*inpart.momentum(); 
    newME(0,ix,0) = fact * decay[1]->momentum().dot(vtemp);
  }
  ME(newME);
  RhoDMatrix rhoin(PDT::Spin0);
  rhoin.average();
  // test of the matrix element
//   double me=newME.contract(rhoin).real();
//   Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					     decay[1]->mass());
//   double test = 2.*pow<4,1>(pcm)*sqr(_coupling[imode()]*inpart.mass())/
//     3./pow<4,1>(decay[0]->mass());
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << (me-test)/(me+test) << "\n";
  // output the answer
  return newME.contract(rhoin).real();
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarMesonTensorScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
						   double & coupling) const {
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); 
  bool order(false);
  int imode(-1);
  do {
    if(id   ==_incoming[ix]) {
      if(id1==_outgoingT[ix]&&id2==_outgoingS[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoingT[ix]&&id1==_outgoingS[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoingT[ix]&&id2bar==_outgoingS[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoingT[ix]&&id1bar==_outgoingS[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode]*dm.parent()->mass();
  itype = 11;
  return order;
}

// output the setup information for the particle database
void ScalarMesonTensorScalarDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "set " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "set " << fullName() << ":OutgoingTensor " << ix << " " 
	     << _outgoingT[ix] << "\n";
      output << "set " << fullName() << ":OutgoingScalar " << ix << " " 
	     << _outgoingS[ix] << "\n";
      output << "set " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix]*MeV << "\n";
      output << "set " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << fullName() << ":OutgoingTensor " << ix << " " 
	     << _outgoingT[ix] << "\n";
      output << "insert " << fullName() << ":OutgoingScalar " << ix << " " 
	     << _outgoingS[ix] << "\n";
      output << "insert " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix]*MeV << "\n";
      output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
