// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarPScalarVectorDecayer class.
//

#include "PScalarPScalarVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarPScalarVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

PScalarPScalarVectorDecayer::PScalarPScalarVectorDecayer() 
{
  // decay pi' to rho pi
  _incoming.push_back( 100111);_outgoingV.push_back( 213);_outgoingP.push_back(-211);
  _incoming.push_back( 100211);_outgoingV.push_back( 213);_outgoingP.push_back( 111);
  _incoming.push_back( 100211);_outgoingV.push_back( 113);_outgoingP.push_back( 211);
  _coupling.push_back(3.57);_maxweight.push_back(4.5);
  _coupling.push_back(3.57);_maxweight.push_back(4.5);
  _coupling.push_back(3.57);_maxweight.push_back(4.5);
  // K' to K rho
  _incoming.push_back( 100311);_outgoingP.push_back( 311);_outgoingV.push_back( 113);
  _incoming.push_back( 100321);_outgoingP.push_back( 321);_outgoingV.push_back( 113);
  _incoming.push_back( 100311);_outgoingP.push_back( 321);_outgoingV.push_back(-213);
  _incoming.push_back( 100321);_outgoingP.push_back( 311);_outgoingV.push_back( 213);
  _coupling.push_back(1.);_maxweight.push_back(4.);
  _coupling.push_back(1.);_maxweight.push_back(4.);
  _coupling.push_back(1.41);_maxweight.push_back(4.);
  _coupling.push_back(1.41);_maxweight.push_back(4.);
  // K' to K* pi
  _incoming.push_back( 100311);_outgoingV.push_back( 313);_outgoingP.push_back( 111);
  _incoming.push_back( 100321);_outgoingV.push_back( 323);_outgoingP.push_back( 111);
  _incoming.push_back( 100311);_outgoingV.push_back( 323);_outgoingP.push_back(-211);
  _incoming.push_back( 100321);_outgoingV.push_back( 313);_outgoingP.push_back( 211);
  _coupling.push_back(1.55);_maxweight.push_back(2.);
  _coupling.push_back(1.55);_maxweight.push_back(2.);
  _coupling.push_back(2.19);_maxweight.push_back(2.);
  _coupling.push_back(2.19);_maxweight.push_back(2.);
  // eta (1475) to K* K
  _incoming.push_back( 100331);_outgoingV.push_back( 323);_outgoingP.push_back(-321);
  _incoming.push_back( 100331);_outgoingV.push_back( 313);_outgoingP.push_back(-311);
  _coupling.push_back(2.80);_maxweight.push_back(3.5);
  _coupling.push_back(2.80);_maxweight.push_back(3.5);
  // eta (1475) to K* K
  _incoming.push_back( 9020221);_outgoingV.push_back( 323);_outgoingP.push_back(-321);
  _incoming.push_back( 9020221);_outgoingV.push_back( 313);_outgoingP.push_back(-311);
  _coupling.push_back(0.98);_maxweight.push_back(4.);
  _coupling.push_back(0.98);_maxweight.push_back(4.);
  // initial size of the arrays
  _initsize=_incoming.size();
}

PScalarPScalarVectorDecayer::~PScalarPScalarVectorDecayer() {}

bool PScalarPScalarVectorDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id()),id0bar(id0);
  if(dm.parent()->CC()){id0bar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0);
  do
    {
      if(id0   ==_incoming[ix])
	{if((id1   ==_outgoingP[ix]&&id2   ==_outgoingV[ix])||
	    (id2   ==_outgoingP[ix]&&id1   ==_outgoingV[ix])){allowed=true;}}
      if(id0bar==_incoming[ix]&&!allowed)
	{if((id1bar==_outgoingP[ix]&&id2bar==_outgoingV[ix])||
	    (id2bar==_outgoingP[ix]&&id1bar==_outgoingV[ix])){allowed=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}

ParticleVector PScalarPScalarVectorDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode(-1);
  int id(parent.id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0);
  bool cc(false);
  do 
    {
      if(id   ==_incoming[ix])
	{if((id1   ==_outgoingP[ix]&&id2   ==_outgoingV[ix])||
	    (id2   ==_outgoingP[ix]&&id1   ==_outgoingV[ix])){imode=ix;}}
      if(idbar==_incoming[ix])
	{if((id1bar==_outgoingP[ix]&&id2bar==_outgoingV[ix])||
	    (id2bar==_outgoingP[ix]&&id1bar==_outgoingV[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void PScalarPScalarVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingP << _outgoingV << _maxweight;
}

void PScalarPScalarVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingP >> _outgoingV >> _maxweight;
}

ClassDescription<PScalarPScalarVectorDecayer> PScalarPScalarVectorDecayer::initPScalarPScalarVectorDecayer;
// Definition of the static class description member.

void PScalarPScalarVectorDecayer::Init() {

  static ClassDocumentation<PScalarPScalarVectorDecayer> documentation
    ("The \\classname{PScalarPScalarVectorDecayer} class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarPScalarVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceOutgoingScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing pseudoscalar meson",
     &PScalarPScalarVectorDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceOutgoingVector
    ("OutgoingVector",
     "The PDG code for the outgoing vector meson",
     &PScalarPScalarVectorDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarPScalarVectorDecayer::_coupling,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarPScalarVectorDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

double PScalarPScalarVectorDecayer::me2(bool vertex, const int ichan,
					const Particle & inpart,
					const ParticleVector & decay) const
{
  // set up spins for the incoming particles
  ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  // set up the spin info for the outgoing particles
  vector<LorentzPolarizationVector> eps;
  ScalarWaveFunction(decay[0],outgoing,true,vertex);
  VectorWaveFunction(eps,decay[1],outgoing,true,false,vertex);
  // calculate the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin1);
  Lorentz5Momentum psum(inpart.momentum()+decay[0]->momentum());
  for(unsigned int ix=0;ix<3;++ix)
    {newME(0,0,ix)=_coupling[imode()]/inpart.mass()*(eps[ix]*psum);}
  ME(newME);
  RhoDMatrix rhoin(PDT::Spin0);rhoin.average();
  return newME.contract(rhoin).real();
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarPScalarVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & mecode,
					       double & coupling) const
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
  unsigned int ix(0); bool order;
  do 
    {
      if(id   ==_incoming[ix])
	{
	  if(id1==_outgoingP[ix]&&id2==_outgoingV[ix]){imode=ix;order=true;}
	  if(id2==_outgoingP[ix]&&id1==_outgoingV[ix]){imode=ix;order=false;}
	}
      if(idbar==_incoming[ix]&&imode<0)
	{
	  if(id1bar==_outgoingP[ix]&&id2bar==_outgoingV[ix]){imode=ix;order=true;}
	  if(id2bar==_outgoingP[ix]&&id1bar==_outgoingV[ix]){imode=ix;order=false;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode];
  mecode=10;
  return order;
}

// output the setup information for the particle database
void PScalarPScalarVectorDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":Points " << _npoint << "\n";
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming " << ix << " " 
		 << _incoming[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingPScalar " << ix << " " 
		 << _outgoingP[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingVector " << ix << " " 
		 << _outgoingV[ix] << "\n";
	  output << "set " << fullName() << ":Coupling " << ix << " " 
		 << _coupling[ix] << "\n";
	  output << "set " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming " << ix << " " 
		 << _incoming[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingPScalar " << ix << " " 
		 << _outgoingP[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingVector " << ix << " " 
		 << _outgoingV[ix] << "\n";
	  output << "insert " << fullName() << ":Coupling " << ix << " " 
		 << _coupling[ix] << "\n";
	  output << "insert " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	}
    }
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
