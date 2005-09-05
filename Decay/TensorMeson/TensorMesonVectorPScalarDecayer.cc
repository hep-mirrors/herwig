// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorPScalarDecayer class.
//

#include "TensorMesonVectorPScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonVectorPScalarDecayer.tcc"
#endif
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::outgoing;
using Herwig::Helicity::EpsFunction;

TensorMesonVectorPScalarDecayer::TensorMesonVectorPScalarDecayer() 
{
  // a_2 -> rho pi
  _incoming.push_back( 115);_outgoingV.push_back( 213);_outgoingP.push_back(-211);
  _incoming.push_back( 215);_outgoingV.push_back( 113);_outgoingP.push_back( 211);
  _incoming.push_back( 215);_outgoingV.push_back( 213);_outgoingP.push_back( 111);
  _coupling.push_back(21.1/GeV2);_maxweight.push_back(10.);
  _coupling.push_back(21.1/GeV2);_maxweight.push_back(9.);
  _coupling.push_back(21.1/GeV2);_maxweight.push_back(9.);
  // a_2+/- -> gamma pi+/-
  _incoming.push_back( 215);_outgoingV.push_back( 22);_outgoingP.push_back( 211);
  _coupling.push_back(0.551/GeV2);_maxweight.push_back(2.);
  // k_2 -> K_2 omega
  _incoming.push_back( 315);_outgoingV.push_back(223);_outgoingP.push_back( 311);
  _incoming.push_back( 325);_outgoingV.push_back(223);_outgoingP.push_back( 321);
  _coupling.push_back(11.66/GeV2);_maxweight.push_back(17.);
  _coupling.push_back(11.66/GeV2);_maxweight.push_back(20.5);
  // k_2+/- -> K+/- gamma
  _incoming.push_back( 325);_outgoingV.push_back(22);_outgoingP.push_back( 321);
  _coupling.push_back(0.553/GeV2);_maxweight.push_back(2.2);
  // B_c2 -> B_c gamma
  _incoming.push_back( 545);_outgoingV.push_back(22);_outgoingP.push_back( 541);
  _coupling.push_back(0.651/GeV2);_maxweight.push_back(2.);
  // K_2 -> K rho
  _incoming.push_back( 325);_outgoingV.push_back( 113);_outgoingP.push_back( 321);
  _incoming.push_back( 325);_outgoingV.push_back( 213);_outgoingP.push_back( 311);
  _incoming.push_back( 315);_outgoingV.push_back( 113);_outgoingP.push_back( 311);
  _incoming.push_back( 315);_outgoingV.push_back(-213);_outgoingP.push_back( 321);
  _coupling.push_back(10.14/GeV2);_maxweight.push_back(9.);
  _coupling.push_back(14.33/GeV2);_maxweight.push_back(9.);
  _coupling.push_back(10.14/GeV2);_maxweight.push_back(9.);
  _coupling.push_back(14.33/GeV2);_maxweight.push_back(9.);
  // K_2 -> K* pi 
  _incoming.push_back( 325);_outgoingV.push_back( 323);_outgoingP.push_back( 111);
  _incoming.push_back( 325);_outgoingV.push_back( 313);_outgoingP.push_back( 211);
  _incoming.push_back( 315);_outgoingV.push_back( 313);_outgoingP.push_back( 111);
  _incoming.push_back( 315);_outgoingV.push_back( 323);_outgoingP.push_back(-211);
  _coupling.push_back(9.627/GeV2);_maxweight.push_back(13);
  _coupling.push_back(13.62/GeV2);_maxweight.push_back(11);
  _coupling.push_back(9.627/GeV2);_maxweight.push_back(8.);
  _coupling.push_back(13.62/GeV2);_maxweight.push_back(8.);
  // D_2 -> D* pi 
  _incoming.push_back( 425);_outgoingV.push_back( 423);_outgoingP.push_back( 111);
  _incoming.push_back( 425);_outgoingV.push_back( 413);_outgoingP.push_back(-211);
  _incoming.push_back( 415);_outgoingV.push_back( 413);_outgoingP.push_back( 111);
  _incoming.push_back( 415);_outgoingV.push_back( 423);_outgoingP.push_back( 211);
  _coupling.push_back(5.931/GeV2);_maxweight.push_back(2.2);
  _coupling.push_back(8.387/GeV2);_maxweight.push_back(2.4);
  _coupling.push_back(5.931/GeV2);_maxweight.push_back(2.4);
  _coupling.push_back(8.387/GeV2);_maxweight.push_back(2.);
  // D_s2 -> D* K
  _incoming.push_back( 435);_outgoingV.push_back( 423);_outgoingP.push_back( 321);
  _incoming.push_back( 435);_outgoingV.push_back( 413);_outgoingP.push_back( 311);
  _coupling.push_back(13.10/GeV2);_maxweight.push_back(2.2);
  _coupling.push_back(13.10/GeV2);_maxweight.push_back(2.5);
  // B_2 -> B* pi 
  _incoming.push_back( 525);_outgoingV.push_back( 523);_outgoingP.push_back( 111);
  _incoming.push_back( 525);_outgoingV.push_back( 513);_outgoingP.push_back( 211);
  _incoming.push_back( 515);_outgoingV.push_back( 513);_outgoingP.push_back( 111);
  _incoming.push_back( 515);_outgoingV.push_back( 523);_outgoingP.push_back(-211);
  _coupling.push_back(4.384/GeV2);_maxweight.push_back(2.1);
  _coupling.push_back(6.199/GeV2);_maxweight.push_back(2.1);
  _coupling.push_back(4.384/GeV2);_maxweight.push_back(2.1);
  _coupling.push_back(6.199/GeV2);_maxweight.push_back(2.1);
  // D_s2
  _incoming.push_back( 435);_outgoingV.push_back( 423);_outgoingP.push_back( 321);
  _incoming.push_back( 435);_outgoingV.push_back( 413);_outgoingP.push_back( 311);
  _coupling.push_back(13.09/GeV2);_maxweight.push_back(2.2);
  _coupling.push_back(13.09/GeV2);_maxweight.push_back(2.5);
  // B_s2
  _incoming.push_back( 535);_outgoingV.push_back( 523);_outgoingP.push_back(-321);
  _incoming.push_back( 535);_outgoingV.push_back( 513);_outgoingP.push_back(-311);
  _coupling.push_back(8.12/GeV2);_maxweight.push_back(2.4);
  _coupling.push_back(8.12/GeV2);_maxweight.push_back(2.1);
  // upsilon_2(1d) to chi_b gamma
  _incoming.push_back(20555);_outgoingV.push_back(  22);_outgoingP.push_back(10551);
  _coupling.push_back(8.12/GeV2);_maxweight.push_back(2.4);
  // initial size of the arrays
  _initsize=_incoming.size();
}

TensorMesonVectorPScalarDecayer::~TensorMesonVectorPScalarDecayer() {}

bool TensorMesonVectorPScalarDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
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

ParticleVector TensorMesonVectorPScalarDecayer::decay(const DecayMode & dm,
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

void TensorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoingV << _outgoingP << _maxweight << _coupling;}

void TensorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> _coupling;}

ClassDescription<TensorMesonVectorPScalarDecayer> TensorMesonVectorPScalarDecayer::initTensorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void TensorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorPScalarDecayer> documentation
    ("The \\classname{TensorMesonVectorPScalarDecayer} class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceOutcomingV
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1particle",
     &TensorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceOutcomingP
    ("OutgoingScalar",
     "The PDG code for the outgoing pseudoscalar meson",
     &TensorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,InvEnergy2> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMesonVectorPScalarDecayer::_coupling,
     0, 0, 0, 0., 100./GeV2, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// the hadronic tensor
vector<LorentzTensor> 
TensorMesonVectorPScalarDecayer::decayTensor(const bool vertex,const int, 
					     const Particle & inpart,
					     const ParticleVector & decay) const
{
  // storage for the tensor
  vector<LorentzTensor> temp;
  vector<LorentzPolarizationVector> vwave;
  // check for photons
  bool photon(_outgoingV[imode()]==ParticleID::gamma);
  // work out which is the vector

  // workaround for gcc 3.2.3 bug
  // set up the spin information for ther decay products
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr mytemp = decay[1];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);

  VectorWaveFunction(vwave,decay[0],outgoing,true,photon,vertex);
  double fact(_coupling[imode()]/inpart.mass());
  for(unsigned int ix=0;ix<3;++ix)
    {
      if(ix==1&&photon){temp.push_back(LorentzTensor());}
      else
	{
	  temp.push_back(LorentzTensor(decay[1]->momentum(),
				       fact*EpsFunction::product(decay[0]->momentum(),
								 vwave[ix],
								 decay[1]->momentum())));
	}
    }
  // return the answer
  return temp;
}

bool TensorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  unsigned int ix(0); bool order(false);
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
  coupling=_coupling[imode]*dm.parent()->mass()*dm.parent()->mass();
  mecode=8;
  return order;
}

void TensorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output) const
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
	  output << "set " << fullName() << ":OutgoingVector " << ix << " " 
		 << _outgoingV[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingScalar " << ix << " " 
		 << _outgoingP[ix] << "\n";
	  output << "set " << fullName() << ":Coupling " << ix << " " 
		 << _coupling[ix] << "\n";
	  output << "set " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming " << ix << " " 
		 << _incoming[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingVector " << ix << " " 
		 << _outgoingV[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingScalar " << ix << " " 
		 << _outgoingP[ix] << "\n";
	  output << "insert " << fullName() << ":Coupling " << ix << " " 
		 << _coupling[ix] << "\n";
	  output << "insert " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	}
    }
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
