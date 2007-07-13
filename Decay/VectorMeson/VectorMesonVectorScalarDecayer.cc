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

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::ScalarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::incoming;
using Helicity::outgoing;

void VectorMesonVectorScalarDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingV.size()||isize!=_outgoingS.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    {throw InitException() << "Inconsistent parameters in "
			   << "VectorMesonVectorScalarDecayer::doinit()" 
			   << Exception::abortnow;}
  // set up the integration channels
  vector<double> wgt(0);
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0]=getParticleData(_incoming[ix]);
      extpart[1]=getParticleData(_outgoingV[ix]);
      extpart[2]=getParticleData(_outgoingS[ix]);
      mode=new DecayPhaseSpaceMode(extpart,this);
      addMode(mode,_maxweight[ix],wgt);
    }
}

VectorMesonVectorScalarDecayer::VectorMesonVectorScalarDecayer() 
{
  // intermediates
  generateIntermediates(false);
  // reserve sizes of the vectors
  _incoming.reserve(20);_outgoingV.reserve(20);_outgoingS.reserve(20);
  _coupling.reserve(20);_maxweight.reserve(20);
  // decay of the phi to the a_0 and f_0 and a photon
  _incoming.push_back(333);_outgoingV.push_back(22);_outgoingS.push_back(9000111);
  _coupling.push_back(0.154/GeV);_maxweight.push_back(17.);
  _incoming.push_back(333);_outgoingV.push_back(22);_outgoingS.push_back(9010221);
  _coupling.push_back(0.267/GeV);_maxweight.push_back(14.);
  // Jpsi decayers
  _incoming.push_back(443);_outgoingV.push_back(22);_outgoingS.push_back(10331);
  _coupling.push_back(0.00207/GeV);_maxweight.push_back(3.);
  _incoming.push_back(443);_outgoingV.push_back(223);_outgoingS.push_back(10331);
  _coupling.push_back(0.00144/GeV);_maxweight.push_back(9.);
  _incoming.push_back(443);_outgoingV.push_back(333);_outgoingS.push_back(10331);
  _coupling.push_back(0.00127/GeV);_maxweight.push_back(9.);
  _incoming.push_back(443);_outgoingV.push_back(333);_outgoingS.push_back(9010221);
  _coupling.push_back(0.00070/GeV);_maxweight.push_back(12.);
  _incoming.push_back(443);_outgoingV.push_back(223);_outgoingS.push_back(9010221);
  _coupling.push_back(0.00048/GeV);_maxweight.push_back(13.);
  // upsilon(2s)
  _incoming.push_back(100553);_outgoingV.push_back(22);_outgoingS.push_back(10551);
  _coupling.push_back(0.122/GeV);_maxweight.push_back(2.5);
  // upsilon(3s)
  _incoming.push_back(200553);_outgoingV.push_back(22);_outgoingS.push_back(110551);
  _coupling.push_back(0.174/GeV);_maxweight.push_back(2.5);
  // psi2s decays
  _incoming.push_back(100443);_outgoingV.push_back(22);_outgoingS.push_back(10441);
  _coupling.push_back(0.229/GeV);_maxweight.push_back(5.);
  _incoming.push_back(100443);_outgoingV.push_back(22);_outgoingS.push_back(331);
  _coupling.push_back(0.0464/GeV);_maxweight.push_back(2.1);
  _incoming.push_back(100443);_outgoingV.push_back(22);_outgoingS.push_back(10331);
  _coupling.push_back(0.000695/GeV);_maxweight.push_back(2.5);
  _incoming.push_back(100443);_outgoingV.push_back(333);_outgoingS.push_back(9010221);
  _coupling.push_back(0.000530/GeV);_maxweight.push_back(10.);
  // rho' to sigma rho
  _incoming.push_back(100213);_outgoingV.push_back(213);_outgoingS.push_back(9000221);
  _incoming.push_back(100113);_outgoingV.push_back(113);_outgoingS.push_back(9000221);
  _coupling.push_back(0.174/GeV);_maxweight.push_back(2.5);
  _coupling.push_back(0.174/GeV);_maxweight.push_back(2.5);
  // initial size of the arrays
  _initsize = _coupling.size();
}

VectorMesonVectorScalarDecayer::~VectorMesonVectorScalarDecayer() {}

int VectorMesonVectorScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0);
  cc=false;
  do 
    {
      if(id   ==_incoming[ix])
	{if((id1   ==_outgoingV[ix]&&id2   ==_outgoingS[ix])||
	    (id2   ==_outgoingV[ix]&&id1   ==_outgoingS[ix])){imode=ix;}}
      if(idbar==_incoming[ix])
	{if((id1bar==_outgoingV[ix]&&id2bar==_outgoingS[ix])||
	    (id2bar==_outgoingV[ix]&&id1bar==_outgoingS[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void VectorMesonVectorScalarDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoingV << _outgoingS << _maxweight << ounit(_coupling,1/GeV);}

void VectorMesonVectorScalarDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoingV >> _outgoingS >> _maxweight >> iunit(_coupling,1/GeV);}

ClassDescription<VectorMesonVectorScalarDecayer> VectorMesonVectorScalarDecayer::initVectorMesonVectorScalarDecayer;
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
  // return the answer
  return newME.contract(rhoin).real();
}

bool VectorMesonVectorScalarDecayer::twoBodyMEcode(const DecayMode & dm,
						     int & mecode,
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
  do 
    {
      if(id==_incoming[ix])
	{
	  if(id1   ==_outgoingV[ix]&&id2   ==_outgoingS[ix]){imode=ix;order=true;}
	  if(id2   ==_outgoingV[ix]&&id1   ==_outgoingS[ix]){imode=ix;order=false;}
	}
      if(idbar==_incoming[ix]&&imode<0)
	{
	  if(id1bar==_outgoingV[ix]&&id2bar==_outgoingS[ix]){imode=ix;order=true;}
	  if(id2bar==_outgoingV[ix]&&id1bar==_outgoingS[ix]){imode=ix;order=false;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling = _coupling[imode]*dm.parent()->mass();  
  mecode = 4;
  return order;
}

void VectorMesonVectorScalarDecayer::dataBaseOutput(ofstream & output,
						    bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
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
      else
	{
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
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
