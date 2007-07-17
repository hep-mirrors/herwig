// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorScalarDecayer class.
//
#include "VectorMesonVectorScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

VectorMesonVectorScalarDecayer::VectorMesonVectorScalarDecayer() 
  : _coupling(15), _incoming(15), _outgoingV(15), 
    _outgoingS(15), _maxweight(15) {
  // decay of the phi to the a_0 and f_0 and a photon
  _incoming[0] = 333; _outgoingV[0] = 22; _outgoingS[0] = 9000111; 
  _coupling[0] = 0.154/GeV; _maxweight[0] = 17.; 
  _incoming[1] = 333; _outgoingV[1] = 22; _outgoingS[1] = 9010221; 
  _coupling[1] = 0.267/GeV; _maxweight[1] = 14.; 
  // Jpsi decayers
  _incoming[2] = 443; _outgoingV[2] = 22; _outgoingS[2] = 10331; 
  _coupling[2] = 0.00207/GeV; _maxweight[2] = 3.; 
  _incoming[3] = 443; _outgoingV[3] = 223; _outgoingS[3] = 10331; 
  _coupling[3] = 0.00144/GeV; _maxweight[3] = 9.; 
  _incoming[4] = 443; _outgoingV[4] = 333; _outgoingS[4] = 10331; 
  _coupling[4] = 0.00127/GeV; _maxweight[4] = 9.; 
  _incoming[5] = 443; _outgoingV[5] = 333; _outgoingS[5] = 9010221; 
  _coupling[5] = 0.00070/GeV; _maxweight[5] = 12.; 
  _incoming[6] = 443; _outgoingV[6] = 223; _outgoingS[6] = 9010221; 
  _coupling[6] = 0.00048/GeV; _maxweight[6] = 13.; 
  // upsilon(2s)
  _incoming[7] = 100553; _outgoingV[7] = 22; _outgoingS[7] = 10551; 
  _coupling[7] = 0.122/GeV; _maxweight[7] = 2.5; 
  // upsilon(3s)
  _incoming[8] = 200553; _outgoingV[8] = 22; _outgoingS[8] = 110551; 
  _coupling[8] = 0.174/GeV; _maxweight[8] = 2.5; 
  // psi2s decays
  _incoming[9] = 100443; _outgoingV[9] = 22; _outgoingS[9] = 10441; 
  _coupling[9] = 0.229/GeV; _maxweight[9] = 5.; 
  _incoming[10] = 100443; _outgoingV[10] = 22; _outgoingS[10] = 331; 
  _coupling[10] = 0.0464/GeV; _maxweight[10] = 2.1; 
  _incoming[11] = 100443; _outgoingV[11] = 22; _outgoingS[11] = 10331; 
  _coupling[11] = 0.000695/GeV; _maxweight[11] = 2.5; 
  _incoming[12] = 100443; _outgoingV[12] = 333; _outgoingS[12] = 9010221; 
  _coupling[12] = 0.000530/GeV; _maxweight[12] = 10.; 
  // rho' to sigma rho
  _incoming[13] = 100213; _outgoingV[13] = 213; _outgoingS[13] = 9000221; 
  _coupling[13] = 0.174/GeV; _maxweight[13] = 2.5; 
  _incoming[14] = 100113; _outgoingV[14] = 113; _outgoingS[14] = 9000221; 
  _coupling[14] = 0.174/GeV; _maxweight[14] = 2.5; 
  // initial size of the arrays
  _initsize = _coupling.size();
  // intermediates
  generateIntermediates(false);
}

void VectorMesonVectorScalarDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingV.size()||isize!=_outgoingS.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters in "
			  << "VectorMesonVectorScalarDecayer::doinit()" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoingV[ix]);
    extpart[2]=getParticleData(_outgoingS[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int VectorMesonVectorScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  // must be two outgoing particles
  if(dm.products().size()!=2) return -1;
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
  int imode(-1);
  do {
    if(id   ==_incoming[ix]) {
      if((id1   ==_outgoingV[ix]&&id2   ==_outgoingS[ix])||
	 (id2   ==_outgoingV[ix]&&id1   ==_outgoingS[ix])) imode=ix;
    }
    if(idbar==_incoming[ix]) {
      if((id1bar==_outgoingV[ix]&&id2bar==_outgoingS[ix])||
	 (id2bar==_outgoingV[ix]&&id1bar==_outgoingS[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void VectorMesonVectorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoingV << _outgoingS << _maxweight << ounit(_coupling,1/GeV);
}

void VectorMesonVectorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoingV >> _outgoingS >> _maxweight >> iunit(_coupling,1/GeV);
}

ClassDescription<VectorMesonVectorScalarDecayer> 
VectorMesonVectorScalarDecayer::initVectorMesonVectorScalarDecayer;
// Definition of the static class description member.

void VectorMesonVectorScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorScalarDecayer> documentation
    ("The VectorMesonVectorScalarDecayer class is designed for the "
     "decay of a vector meson to a vector meson, or the photon, and a "
     "scalar meson.");

  static ParVector<VectorMesonVectorScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonVectorScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1 particle",
     &VectorMesonVectorScalarDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorScalarDecayer,int> interfaceOutcomingScalar
    ("OutgoingScalar",
     "The PDG code for the outgoing spin-0 particle",
     &VectorMesonVectorScalarDecayer::_outgoingS,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonVectorScalarDecayer::_coupling,
     1/GeV, 0, 0/GeV, 0./GeV, 100./GeV, false, false, true);

  static ParVector<VectorMesonVectorScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonVectorScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);
}

double VectorMesonVectorScalarDecayer::me2(bool vertex, const int,
					   const Particle & inpart,
					   const ParticleVector & decay) const {
  // wavefunction for the decaying particle
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // check for photons
  bool photon=_outgoingV[imode()]==ParticleID::gamma;
  // set up the spin information for the decay products
  vector<LorentzPolarizationVector> vout;
  VectorWaveFunction(vout,decay[0],outgoing,true,photon,vertex);
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr myvout = decay[1];
  ScalarWaveFunction(myvout,outgoing,true,vertex);
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  Energy2 p0dotpv(inpart.momentum()*decay[0]->momentum());
  complex<Energy> epsdot(0.*MeV);
  InvEnergy2 pre(_coupling[imode()]/inpart.mass());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1&&photon) {
      for(unsigned int iy=0;iy<3;++iy) newME(iy,ix,0)=0.;
    }
    else {
      epsdot=vout[ix]*inpart.momentum();
      for(unsigned int iy=0;iy<3;++iy) {
	newME(iy,ix,0)=pre*invec[iy].dot(p0dotpv*vout[ix]-
					 epsdot*decay[0]->momentum());
      }
    }
  }
  ME(newME);
  // test of the matrix element
//   double me = newME.contract(rhoin).real();
//   Energy pcm=Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					   decay[1]->mass());
//   double test = sqr(_coupling[imode()])/3.*(2.*sqr(pcm)+3.*sqr(decay[0]->mass()));
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return newME.contract(rhoin).real();
}

bool VectorMesonVectorScalarDecayer::twoBodyMEcode(const DecayMode & dm,
						     int & mecode,
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
    if(id==_incoming[ix]) {
      if(id1   ==_outgoingV[ix]&&id2   ==_outgoingS[ix]) {
	imode=ix;
	order=true;
      }
      if(id2   ==_outgoingV[ix]&&id1   ==_outgoingS[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoingV[ix]&&id2bar==_outgoingS[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoingV[ix]&&id1bar==_outgoingS[ix]) {
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

void VectorMesonVectorScalarDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "set " << fullName() << ":Incoming " << ix << " "
	     << _incoming[ix] << "\n";
      output << "set " << fullName() << ":OutgoingVector " << ix << " "
	     << _outgoingV[ix] << "\n";
      output << "set " << fullName() << ":OutgoingScalar " << ix << " "
	     << _outgoingS[ix] << "\n";
      output << "set " << fullName() << ":Coupling " << ix << " "
	     << _coupling[ix]*MeV << "\n";
      output << "set " << fullName() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << fullName() << ":Incoming "  << ix << " "
	     << _incoming[ix] << "\n";
      output << "insert " << fullName() << ":OutgoingVector " << ix << " "
	     << _outgoingV[ix] << "\n";
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
