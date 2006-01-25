// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonPVectorPScalarDecayer class.
//

#include "VectorMesonPVectorPScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonPVectorPScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::ScalarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::outgoing;

VectorMesonPVectorPScalarDecayer::VectorMesonPVectorPScalarDecayer() 
{
  // Jpsi to K_1 K
  _incoming.push_back(443);_outgoingA.push_back( 20313);_outgoingP.push_back(-311);
  _incoming.push_back(443);_outgoingA.push_back( 20323);_outgoingP.push_back(-321);
  _coupling.push_back(0.00114/GeV);_maxweight.push_back(12.);
  _coupling.push_back(0.00114/GeV);_maxweight.push_back(12.);
  // Jpsi to b_1 pi
  _incoming.push_back(443);_outgoingA.push_back( 10213);_outgoingP.push_back(-211);
  _incoming.push_back(443);_outgoingA.push_back( 10113);_outgoingP.push_back( 111);
  _coupling.push_back(0.00106/GeV);_maxweight.push_back(10.5);
  _coupling.push_back(0.00106/GeV);_maxweight.push_back(10.5);
  // psi(2s) to K_1 K
  _incoming.push_back(100443);_outgoingA.push_back( 10313);_outgoingP.push_back(-311);
  _incoming.push_back(100443);_outgoingA.push_back( 10323);_outgoingP.push_back(-321);
  _coupling.push_back(0.000898/GeV);_maxweight.push_back(12.);
  _coupling.push_back(0.000898/GeV);_maxweight.push_back(12.);
  // psi(2s) to b_1 pi
  _incoming.push_back(100443);_outgoingA.push_back( 10213);_outgoingP.push_back(-211);
  _incoming.push_back(100443);_outgoingA.push_back( 10113);_outgoingP.push_back( 111);
  _coupling.push_back(0.000464/GeV);_maxweight.push_back(10.5);
  _coupling.push_back(0.000464/GeV);_maxweight.push_back(10.5);
  // rho'' decays
  // to h_1
  _incoming.push_back( 30213);_outgoingA.push_back( 10223);_outgoingP.push_back(211);
  _incoming.push_back( 30113);_outgoingA.push_back( 10223);_outgoingP.push_back(111);
  _coupling.push_back(1.41/GeV);_maxweight.push_back(5.5);
  _coupling.push_back(1.41/GeV);_maxweight.push_back(5.5);
  // to a_1
  _incoming.push_back( 30213);_outgoingA.push_back( 20213);_outgoingP.push_back( 111);
  _incoming.push_back( 30213);_outgoingA.push_back( 20113);_outgoingP.push_back( 211);
  _incoming.push_back( 30113);_outgoingA.push_back( 20213);_outgoingP.push_back(-211);
  _coupling.push_back(1.29/GeV);_maxweight.push_back(4.);
  _coupling.push_back(1.29/GeV);_maxweight.push_back(4.);
  _coupling.push_back(1.29/GeV);_maxweight.push_back(4.);
  //  rho' decays
  // to h_1
  _incoming.push_back( 100213);_outgoingA.push_back( 10223);_outgoingP.push_back(211);
  _incoming.push_back( 100113);_outgoingA.push_back( 10223);_outgoingP.push_back(111);
  _coupling.push_back(1.95/GeV);_maxweight.push_back(5.);
  _coupling.push_back(1.95/GeV);_maxweight.push_back(5.);
  // to a_1
  _incoming.push_back( 100213);_outgoingA.push_back( 20213);_outgoingP.push_back(111);
  _incoming.push_back( 100213);_outgoingA.push_back( 20113);_outgoingP.push_back(211);
  _incoming.push_back( 100113);_outgoingA.push_back( 20213);_outgoingP.push_back(-211);
  _coupling.push_back(3.73/GeV);_maxweight.push_back(4.);
  _coupling.push_back(3.73/GeV);_maxweight.push_back(4.);
  _coupling.push_back(3.73/GeV);_maxweight.push_back(4.);
  // omega' to b pi
  _incoming.push_back(100223);_outgoingA.push_back( 10113);_outgoingP.push_back( 111);
  _incoming.push_back(100223);_outgoingA.push_back( 10213);_outgoingP.push_back(-211);
  _coupling.push_back(1.64/GeV);_maxweight.push_back(7.);
  _coupling.push_back(1.64/GeV);_maxweight.push_back(6.);
  // initial size of the arrays
  _initsize = _coupling.size();
}

VectorMesonPVectorPScalarDecayer::~VectorMesonPVectorPScalarDecayer() {}

bool VectorMesonPVectorPScalarDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id()),id0bar(id0);
  if(dm.parent()->CC()){id0bar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0);
  do
    {
      if(id0   ==_incoming[ix])
	{if((id1   ==_outgoingA[ix]&&id2   ==_outgoingP[ix])||
	    (id2   ==_outgoingA[ix]&&id1   ==_outgoingP[ix])){allowed=true;}}
      if(id0bar==_incoming[ix]&&!allowed)
	{if((id1bar==_outgoingA[ix]&&id2bar==_outgoingP[ix])||
	    (id2bar==_outgoingA[ix]&&id1bar==_outgoingP[ix])){allowed=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}

ParticleVector VectorMesonPVectorPScalarDecayer::decay(const DecayMode & dm,
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
	{if((id1   ==_outgoingA[ix]&&id2   ==_outgoingP[ix])||
	    (id2   ==_outgoingA[ix]&&id1   ==_outgoingP[ix])){imode=ix;}}
      if(idbar==_incoming[ix])
	{if((id1bar==_outgoingA[ix]&&id2bar==_outgoingP[ix])||
	    (id2bar==_outgoingA[ix]&&id1bar==_outgoingP[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void VectorMesonPVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const {
os << _incoming << _outgoingA << _outgoingP << _maxweight << _coupling;}

void VectorMesonPVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) {
is >> _incoming >> _outgoingA >> _outgoingP >> _maxweight >> _coupling;}

ClassDescription<VectorMesonPVectorPScalarDecayer> VectorMesonPVectorPScalarDecayer::initVectorMesonPVectorPScalarDecayer;
// Definition of the static class description member.

void VectorMesonPVectorPScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonPVectorPScalarDecayer> documentation
    ("The VectorMesonPVectorPScalarDecayer class is designed for the "
     "decay of a vector meson to a pseudovector meson and a "
     "pseudoscalar meson.");

  static ParVector<VectorMesonPVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonPVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingPVector",
     "The PDG code for the outgoing spin-1 particle",
     &VectorMesonPVectorPScalarDecayer::_outgoingA,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,int> interfaceOutcomingScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &VectorMesonPVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonPVectorPScalarDecayer::_coupling,
     0, 0, 0, 0./GeV, 100./GeV, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonPVectorPScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);
}

// the hadronic currents 
vector<LorentzPolarizationVector>  
VectorMesonPVectorPScalarDecayer::decayCurrent(const bool vertex, const int, 
					      const Particle & inpart,
					      const ParticleVector & decay) const
{
  // storage for the current
  vector<LorentzPolarizationVector> temp;
  // set up the spin information for the decay products
  VectorWaveFunction(temp,decay[0],outgoing,true,false,vertex);

  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr mytemp = decay[1];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);

  // calculate the currents
  Energy2 p0dotpv(inpart.momentum()*decay[0]->momentum());
  Complex epsdot(0.),pre(_coupling[imode()]/inpart.mass());
  for(unsigned int ix=0;ix<3;++ix)
    {
      epsdot=temp[ix]*inpart.momentum();
      temp[ix]=pre*(p0dotpv*temp[ix]-epsdot*decay[0]->momentum());
    }
  return temp;
 }

bool VectorMesonPVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,
						     int & mecode,
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
      if(id==_incoming[ix])
	{
	  if(id1   ==_outgoingA[ix]&&id2   ==_outgoingP[ix]){imode=ix;order=true;}
	  if(id2   ==_outgoingA[ix]&&id1   ==_outgoingP[ix]){imode=ix;order=false;}
	}
      if(idbar==_incoming[ix]&&imode<0)
	{
	  if(id1bar==_outgoingA[ix]&&id2bar==_outgoingP[ix]){imode=ix;order=true;}
	  if(id2bar==_outgoingA[ix]&&id1bar==_outgoingP[ix]){imode=ix;order=false;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling = _coupling[imode]*dm.parent()->mass();  
  mecode = 4;
  return order;
}
void VectorMesonPVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						      bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  VectorMesonDecayerBase::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming " << ix << " "
		 << _incoming[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingPVector " << ix << " "
		 << _outgoingA[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingPScalar " << ix << " "
		 << _outgoingP[ix] << "\n";
	  output << "set " << fullName() << ":Coupling " << ix << " "
		 << _coupling[ix] << "\n";
	  output << "set " << fullName() << ":MaxWeight " << ix << " "
		 << _maxweight[ix] << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming "  << ix << " "
		 << _incoming[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingPVector " << ix << " "
		 << _outgoingA[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingPScalar " << ix << " "
		 << _outgoingP[ix] << "\n";
	  output << "insert " << fullName() << ":Coupling " << ix << " "
		 << _coupling[ix] << "\n";
	  output << "insert " << fullName() << ":MaxWeight " << ix << " "
		 << _maxweight[ix] << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
