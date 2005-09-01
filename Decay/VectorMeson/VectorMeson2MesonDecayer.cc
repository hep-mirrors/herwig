// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2MesonDecayer class.
//

#include "VectorMeson2MesonDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::outgoing;

VectorMeson2MesonDecayer::VectorMeson2MesonDecayer() 
{
  // particles and couplings for the different modes
  // rho -> pi pi
  _incoming.push_back( 113);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _incoming.push_back( 213);_outgoing1.push_back( 111);_outgoing2.push_back( 211);
  _coupling.push_back(6.);_coupling.push_back(6.);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // rho' -> pi pi
  _incoming.push_back( 100113);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _incoming.push_back( 100213);_outgoing1.push_back( 111);_outgoing2.push_back( 211);
  _coupling.push_back(3.428);_coupling.push_back(3.428);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // rho'' -> pi pi
  _incoming.push_back( 30113);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _incoming.push_back( 30213);_outgoing1.push_back( 111);_outgoing2.push_back( 211);
  _coupling.push_back(1.611);_coupling.push_back(1.611);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // rho'' -> K K
  _incoming.push_back( 30113);_outgoing1.push_back( 321);_outgoing2.push_back(-321);
  _incoming.push_back( 30113);_outgoing1.push_back( 311);_outgoing2.push_back(-311);
  _incoming.push_back( 30213);_outgoing1.push_back( 321);_outgoing2.push_back(-311);
  _coupling.push_back(0.294);_coupling.push_back(0.294);_coupling.push_back(0.416);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // rho'' -> pi' pi
  _incoming.push_back( 30113);_outgoing1.push_back(100211);_outgoing2.push_back(-211);
  _incoming.push_back( 30213);_outgoing1.push_back(100111);_outgoing2.push_back( 211);
  _incoming.push_back( 30213);_outgoing1.push_back(111);_outgoing2.push_back( 100211);
  _coupling.push_back(7.630);_coupling.push_back(7.630);_coupling.push_back(7.630);
  _maxweight.push_back(5.);_maxweight.push_back(5.);_maxweight.push_back(5.);
  // rho' -> pi' pi
  _incoming.push_back( 100113);_outgoing1.push_back(100211);_outgoing2.push_back(-211);
  _incoming.push_back( 100213);_outgoing1.push_back(100111);_outgoing2.push_back( 211);
  _incoming.push_back( 100213);_outgoing1.push_back(111);_outgoing2.push_back( 100211);
  _coupling.push_back(28.6);_coupling.push_back(28.6);_coupling.push_back(28.6); 
  _maxweight.push_back(5.);_maxweight.push_back(5.);_maxweight.push_back(5.);
  // omega -> pi pi
  _incoming.push_back(223);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _coupling.push_back(0.1847);_maxweight.push_back(2.);
  // K* decays
  _incoming.push_back( 313);_outgoing1.push_back( 321);_outgoing2.push_back(-211);
  _incoming.push_back( 313);_outgoing1.push_back( 311);_outgoing2.push_back( 111);
  _incoming.push_back( 323);_outgoing1.push_back( 311);_outgoing2.push_back( 211);
  _incoming.push_back( 323);_outgoing1.push_back( 321);_outgoing2.push_back( 111);
  _coupling.push_back(4.57);_coupling.push_back(3.23);
  _coupling.push_back(4.57);_coupling.push_back(3.23);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // K*' decays
  _incoming.push_back( 100313);_outgoing1.push_back( 321);_outgoing2.push_back(-211);
  _incoming.push_back( 100313);_outgoing1.push_back( 311);_outgoing2.push_back( 111);
  _incoming.push_back( 100323);_outgoing1.push_back( 311);_outgoing2.push_back( 211);
  _incoming.push_back( 100323);_outgoing1.push_back( 321);_outgoing2.push_back( 111);
  _coupling.push_back(1.296);_coupling.push_back(0.916);
  _coupling.push_back(1.296);_coupling.push_back(0.916);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // K*'' decays
  _incoming.push_back( 30313);_outgoing1.push_back( 321);_outgoing2.push_back(-211);
  _incoming.push_back( 30313);_outgoing1.push_back( 311);_outgoing2.push_back( 111);
  _incoming.push_back( 30323);_outgoing1.push_back( 311);_outgoing2.push_back( 211);
  _incoming.push_back( 30323);_outgoing1.push_back( 321);_outgoing2.push_back( 111);
  _coupling.push_back(3.114);_coupling.push_back(2.201);
  _coupling.push_back(3.114);_coupling.push_back(2.201);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // phi decays
  _incoming.push_back( 333);_outgoing1.push_back( 321);_outgoing2.push_back(-321);
  _incoming.push_back( 333);_outgoing1.push_back( 311);_outgoing2.push_back(-311);
  _incoming.push_back( 333);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _coupling.push_back(4.48);_coupling.push_back(4.59);_coupling.push_back(8.986E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // phi' decays
  _incoming.push_back( 100333);_outgoing1.push_back( 321);_outgoing2.push_back(-321);
  _incoming.push_back( 100333);_outgoing1.push_back( 311);_outgoing2.push_back(-311);
  _coupling.push_back(0.912);_coupling.push_back(0.918);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // excited psi decays
  _incoming.push_back(30443);_outgoing1.push_back( 411);_outgoing2.push_back(-411);
  _incoming.push_back(30443);_outgoing1.push_back( 421);_outgoing2.push_back(-421);
  _coupling.push_back(13.375);_maxweight.push_back(2.);
  _coupling.push_back(13.375);_maxweight.push_back(2.);
  // D* decays
  _incoming.push_back( 423);_outgoing1.push_back( 421);_outgoing2.push_back(111);
  _incoming.push_back( 413);_outgoing1.push_back( 411);_outgoing2.push_back(111);
  _incoming.push_back( 413);_outgoing1.push_back( 421);_outgoing2.push_back(211);
  _coupling.push_back(6.366);_maxweight.push_back(2.);
  _coupling.push_back(6.370);_maxweight.push_back(2.);
  _coupling.push_back(9.019);_maxweight.push_back(2.);
  // D_s* decays
  _incoming.push_back( 433);_outgoing1.push_back( 431);_outgoing2.push_back(111);
  _coupling.push_back(5.635);_maxweight.push_back(2.);
  // K_1 decays to K*_0 pion
  _incoming.push_back( 10323);_outgoing1.push_back( 10321);_outgoing2.push_back( 111);
  _incoming.push_back( 10323);_outgoing1.push_back( 10311);_outgoing2.push_back( 211);
  _incoming.push_back( 10313);_outgoing1.push_back( 10311);_outgoing2.push_back( 111);
  _incoming.push_back( 10313);_outgoing1.push_back( 10321);_outgoing2.push_back(-211);
  _coupling.push_back(20.366);_maxweight.push_back(12.5);
  _coupling.push_back(28.802);_maxweight.push_back(12.5);
  _coupling.push_back(20.366);_maxweight.push_back(12.5);
  _coupling.push_back(28.802);_maxweight.push_back(12.5);
  // K_1 decays to f(1370) kaon
  _incoming.push_back( 10323);_outgoing1.push_back( 321);_outgoing2.push_back( 10221);
  _incoming.push_back( 10313);_outgoing1.push_back( 311);_outgoing2.push_back( 10221);
  _coupling.push_back(38.8);_maxweight.push_back(6.);
  _coupling.push_back(38.8);_maxweight.push_back(6.);
  // K'_1 decays to f(1370) kaon
  _incoming.push_back( 20323);_outgoing1.push_back( 321);_outgoing2.push_back( 10221);
  _incoming.push_back( 20313);_outgoing1.push_back( 311);_outgoing2.push_back( 10221);
  _coupling.push_back(23.34);_maxweight.push_back(6.);
  _coupling.push_back(23.34);_maxweight.push_back(6.);
  // upsilon(4s)
  _incoming.push_back(300553);_outgoing1.push_back(521);_outgoing2.push_back(-521);
  _incoming.push_back(300553);_outgoing1.push_back(511);_outgoing2.push_back(-511);
  _coupling.push_back(23.653);_maxweight.push_back(2.);
  _coupling.push_back(23.653);_maxweight.push_back(2.);
  // jpsi to pions
  _incoming.push_back(443);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _coupling.push_back(2.568E-3);_maxweight.push_back(2.);
  // jpsi to kaons
  _incoming.push_back(443);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _coupling.push_back(1.111E-3);_maxweight.push_back(2.);
  _incoming.push_back(443);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _coupling.push_back(0.873E-3);_maxweight.push_back(2.);
  // psi(2s) to pions
  _incoming.push_back(100443);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _coupling.push_back(0.963E-3);_maxweight.push_back(2.);
  // psi(2s) to kaons
  _incoming.push_back(100443);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _coupling.push_back(0.817E-3);_maxweight.push_back(2.);
  _incoming.push_back(100443);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _coupling.push_back(0.817E-3);_maxweight.push_back(2.);
  // f_1 to a_0 pi
  _incoming.push_back(20223);_outgoing1.push_back( 9000111);_outgoing2.push_back( 111);
  _coupling.push_back(3.035);_maxweight.push_back(10.);
  _incoming.push_back(20223);_outgoing1.push_back( 9000211);_outgoing2.push_back(-211);
  _coupling.push_back(3.035);_maxweight.push_back(10.);
  _incoming.push_back(20223);_outgoing1.push_back(-9000211);_outgoing2.push_back( 211);
  _coupling.push_back(3.035);_maxweight.push_back(10.);
  // f'_1 to a_0 pi
  _incoming.push_back(20333);_outgoing1.push_back( 9000111);_outgoing2.push_back( 111);
  _coupling.push_back(0.954);_maxweight.push_back(10.);
  _incoming.push_back(20333);_outgoing1.push_back( 9000211);_outgoing2.push_back(-211);
  _coupling.push_back(0.954);_maxweight.push_back(10.);
  _incoming.push_back(20333);_outgoing1.push_back(-9000211);_outgoing2.push_back( 211);
  _coupling.push_back(0.954);_maxweight.push_back(10.);
  // initial size of the vectors for the database output
  _initsize=_incoming.size();
}

VectorMeson2MesonDecayer::~VectorMeson2MesonDecayer() {}
  
bool VectorMeson2MesonDecayer::accept(const DecayMode & dm) const {
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
	{if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	    (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])){allowed=true;}}
      if(id0bar==_incoming[ix]&&!allowed)
	{if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	    (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])){allowed=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}
  
ParticleVector VectorMeson2MesonDecayer::decay(const DecayMode & dm,
					       const Particle & parent) const 
{
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
	{if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	    (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])){imode=ix;}}
      if(idbar==_incoming[ix])
	{if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	    (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}
   
void VectorMeson2MesonDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;}
  
void VectorMeson2MesonDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;}
  
ClassDescription<VectorMeson2MesonDecayer> VectorMeson2MesonDecayer::initVectorMeson2MesonDecayer;
  // Definition of the static class description member.

void VectorMeson2MesonDecayer::Init() {
  
  static ClassDocumentation<VectorMeson2MesonDecayer> documentation
    ("The \\classname{VectorMeson2MesonDecayer} class is designed to implement "
     "the decay of vector mesons to 2 scalar mesons via a current which is the "
     "difference of the momenta of the two scalars. The order of the scalar meson "
     "momenta does not matter as it only changes the sign of the matrix element.");

  static ParVector<VectorMeson2MesonDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson2MesonDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &VectorMeson2MesonDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &VectorMeson2MesonDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMeson2MesonDecayer::_coupling,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMeson2MesonDecayer::_maxweight,
     0, 0, 0, -10000000, 10000000, false, false, true);
  
}

// the hadronic currents    
vector<LorentzPolarizationVector>  
VectorMeson2MesonDecayer::decayCurrent(const bool vertex, const int, 
				       const Particle & inpart,
				       const ParticleVector & decay) const
{
  // setup the spininfomation for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix)

    // workaround for gcc 3.2.3 bug
    //ALB {ScalarWaveFunction(decay[ix],outgoing,true,vertex);}
    {PPtr mytemp = decay[ix]; ScalarWaveFunction(mytemp,outgoing,true,vertex);}

  // calculate the current
  return vector<LorentzPolarizationVector>(1,_coupling[imode()]/inpart.mass()*
					   (decay[0]->momentum()-decay[1]->momentum()));
}
 
bool VectorMeson2MesonDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
	  if(id1==_outgoing1[ix]&&id2==_outgoing2[ix]){imode=ix;order=true;}
	  if(id2==_outgoing1[ix]&&id1==_outgoing2[ix]){imode=ix;order=false;}
	}
      if(idbar==_incoming[ix]&&imode<0)
	{
	  if(id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix]){imode=ix;order=true;}
	  if(id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix]){imode=ix;order=false;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode];
  mecode=0;
  return order;
}

// output the setup information for the particle database
void VectorMeson2MesonDecayer::dataBaseOutput(ofstream & output)
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
	  output << "set " << fullName() << ":FirstOutgoing " << ix << " " 
		 << _outgoing1[ix] << "\n";
	  output << "set " << fullName() << ":SecondOutgoing " << ix << " " 
		 << _outgoing2[ix] << "\n";
	  output << "set " << fullName() << ":Coupling " << ix << " " 
		 << _coupling[ix] << "\n";
	  output << "set " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming " << ix << " " 
		 << _incoming[ix] << "\n";
	  output << "insert " << fullName() << ":FirstOutgoing " << ix << " " 
		 << _outgoing1[ix] << "\n";
	  output << "insert " << fullName() << ":SecondOutgoing " << ix << " " 
		 << _outgoing2[ix] << "\n";
	  output << "insert " << fullName() << ":Coupling " << ix << " " 
		 << _coupling[ix] << "\n";
	  output << "insert " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	}
    }
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
