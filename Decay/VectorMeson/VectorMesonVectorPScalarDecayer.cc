// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorPScalarDecayer class.
//

#include "VectorMesonVectorPScalarDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonVectorPScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VectorWaveFunction;
using Helicity::ScalarWaveFunction;
using Helicity::EpsFunction;
using Helicity::outgoing;

VectorMesonVectorPScalarDecayer::VectorMesonVectorPScalarDecayer() 
{
  // rho -> gamma pi modes
  _incoming.push_back( 113);_outgoingV.push_back( 22);_outgoingP.push_back( 111);
  _incoming.push_back( 213);_outgoingV.push_back( 22);_outgoingP.push_back( 211);
  _coupling.push_back(0.2527/GeV);_maxweight.push_back(2.);
  _coupling.push_back(0.2210/GeV);_maxweight.push_back(2.);
  // rho  -> gamma eta mode
  _incoming.push_back( 113);_outgoingV.push_back( 22);_outgoingP.push_back( 221);
  _coupling.push_back(0.492/GeV);_maxweight.push_back(2.);
  // omega -> gamma pi 
  _incoming.push_back( 223);_outgoingV.push_back( 22);_outgoingP.push_back( 111);
  _coupling.push_back(0.7279465/GeV);_maxweight.push_back(2.);
  // omega -> gamma eta
  _incoming.push_back( 223);_outgoingV.push_back( 22);_outgoingP.push_back( 221);
  _coupling.push_back(0.143/GeV);_maxweight.push_back(2.);
  // phi -> gamma pi 
  _incoming.push_back( 333);_outgoingV.push_back( 22);_outgoingP.push_back( 111);
  _coupling.push_back(0.0397/GeV);_maxweight.push_back(2.);
  // phi -> gamma eta
  _incoming.push_back( 333);_outgoingV.push_back( 22);_outgoingP.push_back( 221);
  _coupling.push_back(0.212/GeV);_maxweight.push_back(2.);
  // phi -> gamma eta'
  _incoming.push_back( 333);_outgoingV.push_back( 22);_outgoingP.push_back( 331);
  _coupling.push_back(0.219/GeV);_maxweight.push_back(2.5);
  // phi -> omega pi
  _incoming.push_back( 333);_outgoingV.push_back(223);_outgoingP.push_back( 111);
  _coupling.push_back(0.0417/GeV);_maxweight.push_back(2.5);
  // phi' -> K* K
  _incoming.push_back( 100333);_outgoingV.push_back( 323);_outgoingP.push_back(-321);
  _coupling.push_back(3.934/GeV);_maxweight.push_back(5.);
  _incoming.push_back( 100333);_outgoingV.push_back( 313);_outgoingP.push_back(-311);
  _coupling.push_back(4.011/GeV);_maxweight.push_back(5.);
  // K* -> gamma K 
  _incoming.push_back( 313);_outgoingV.push_back( 22);_outgoingP.push_back( 311);
  _incoming.push_back( 323);_outgoingV.push_back( 22);_outgoingP.push_back( 321);
  _coupling.push_back(0.384/GeV);_maxweight.push_back(2.);
  _coupling.push_back(0.253/GeV);_maxweight.push_back(2.);
  // d* decay
  _incoming.push_back( 423);_outgoingV.push_back(22);_outgoingP.push_back( 421);
  _coupling.push_back(0.616/GeV);_maxweight.push_back(2.);
  _incoming.push_back( 413);_outgoingV.push_back(22);_outgoingP.push_back( 411);
  _coupling.push_back(0.152/GeV);_maxweight.push_back(2.);
  // D_s* decays
  _incoming.push_back( 433);_outgoingV.push_back(22);_outgoingP.push_back( 431);
  _coupling.push_back(0.189/GeV);_maxweight.push_back(2.);
  // B_s* decays
  _incoming.push_back( 533);_outgoingV.push_back(22);_outgoingP.push_back( 531);
  _coupling.push_back(0.235/GeV);_maxweight.push_back(2.);
  // B_c* decays
  _incoming.push_back( 543);_outgoingV.push_back(22);_outgoingP.push_back( 541);
  _coupling.push_back(0.1025/GeV);_maxweight.push_back(2.);
  // B* decay
  _incoming.push_back( 523);_outgoingV.push_back(22);_outgoingP.push_back( 521);
  _coupling.push_back(0.553/GeV);_maxweight.push_back(2.);
  _incoming.push_back( 513);_outgoingV.push_back(22);_outgoingP.push_back( 511);
  _coupling.push_back(0.310/GeV);_maxweight.push_back(2.);
  // rho'' eta rho
  _incoming.push_back( 30113);_outgoingV.push_back( 113);_outgoingP.push_back(221);
  _incoming.push_back( 30213);_outgoingV.push_back( 213);_outgoingP.push_back(221);
  _coupling.push_back(2.600/GeV);_maxweight.push_back(5.);
  _coupling.push_back(2.600/GeV);_maxweight.push_back(5.);
  // rho '' K* K
  _incoming.push_back( 30113);_outgoingV.push_back( 323);_outgoingP.push_back(-321);
  _incoming.push_back( 30113);_outgoingV.push_back( 313);_outgoingP.push_back(-311);
  _incoming.push_back( 30213);_outgoingV.push_back( 323);_outgoingP.push_back(-311);
  _incoming.push_back( 30213);_outgoingV.push_back(-313);_outgoingP.push_back( 321);
  _coupling.push_back(1.357/GeV);_maxweight.push_back(3.);
  _coupling.push_back(1.378/GeV);_maxweight.push_back(3.);
  _coupling.push_back(1.932/GeV);_maxweight.push_back(3.);
  _coupling.push_back(1.932/GeV);_maxweight.push_back(3.);
  // omega'' rho pi
  _incoming.push_back( 30223);_outgoingV.push_back( 213);_outgoingP.push_back(-211);
  _incoming.push_back( 30223);_outgoingV.push_back( 113);_outgoingP.push_back( 111);
  _coupling.push_back(2.996/GeV);_maxweight.push_back(3.);
  _coupling.push_back(2.996/GeV);_maxweight.push_back(3.);
  // omega' rho pi
  _incoming.push_back( 100223);_outgoingV.push_back( 213);_outgoingP.push_back(-211);
  _incoming.push_back( 100223);_outgoingV.push_back( 113);_outgoingP.push_back( 111);
  _coupling.push_back(4.507/GeV);_maxweight.push_back(5.);
  _coupling.push_back(4.507/GeV);_maxweight.push_back(5.);
  // K*''->K* pi decays
  _incoming.push_back( 30313);_outgoingV.push_back( 323);_outgoingP.push_back(-211);
  _incoming.push_back( 30313);_outgoingV.push_back( 313);_outgoingP.push_back( 111);
  _incoming.push_back( 30323);_outgoingV.push_back( 313);_outgoingP.push_back( 211);
  _incoming.push_back( 30323);_outgoingV.push_back( 323);_outgoingP.push_back( 111);
  _coupling.push_back(3.36/GeV);_coupling.push_back(2.38/GeV);
  _coupling.push_back(3.36/GeV);_coupling.push_back(2.38/GeV);
  _maxweight.push_back(5.);_maxweight.push_back(5.);
  _maxweight.push_back(5.);_maxweight.push_back(5.);
  // K*''->K rho decays
  _incoming.push_back( 30313);_outgoingP.push_back( 321);_outgoingV.push_back(-213);
  _incoming.push_back( 30313);_outgoingP.push_back( 311);_outgoingV.push_back( 113);
  _incoming.push_back( 30323);_outgoingP.push_back( 311);_outgoingV.push_back( 213);
  _incoming.push_back( 30323);_outgoingP.push_back( 321);_outgoingV.push_back( 113);
  _coupling.push_back(4.159/GeV);_coupling.push_back(2.939/GeV);
  _coupling.push_back(4.159/GeV);_coupling.push_back(2.939/GeV);
  _maxweight.push_back(5.);_maxweight.push_back(5.);
  _maxweight.push_back(5.);_maxweight.push_back(5.);
  // K*' decays
  _incoming.push_back( 100313);_outgoingV.push_back( 323);_outgoingP.push_back(-211);
  _incoming.push_back( 100313);_outgoingV.push_back( 313);_outgoingP.push_back( 111);
  _incoming.push_back( 100323);_outgoingV.push_back( 313);_outgoingP.push_back( 211);
  _incoming.push_back( 100323);_outgoingV.push_back( 323);_outgoingP.push_back( 111);
  _coupling.push_back(9.469/GeV);_coupling.push_back(6.781/GeV);
  _coupling.push_back(9.469/GeV);_coupling.push_back(6.781/GeV);
  _maxweight.push_back(5.);_maxweight.push_back(5.);
  _maxweight.push_back(5.);_maxweight.push_back(5.);
  // J/psi -> gamma eta_c decay
  _incoming.push_back(443);_outgoingV.push_back(22);_outgoingP.push_back(441);
  _coupling.push_back(0.160/GeV);_maxweight.push_back(20.);
  // J/psi -> gamma eta' decay
  _incoming.push_back(443);_outgoingV.push_back(22);_outgoingP.push_back(331);
  _coupling.push_back(0.00236/GeV);_maxweight.push_back(2.);
  // J/psi -> rho pi decay
  _incoming.push_back(443);_outgoingV.push_back(213);_outgoingP.push_back(-211);
  _incoming.push_back(443);_outgoingV.push_back(113);_outgoingP.push_back(111);
  _coupling.push_back(0.00234/GeV);_maxweight.push_back(4.);
  _coupling.push_back(0.00234/GeV);_maxweight.push_back(4.);
  // J/psi -> K* K decay
  _incoming.push_back(443);_outgoingV.push_back(323);_outgoingP.push_back(-321);
  _incoming.push_back(443);_outgoingV.push_back(313);_outgoingP.push_back(-311);
  _coupling.push_back(0.00177/GeV);_maxweight.push_back(7.);
  _coupling.push_back(0.00177/GeV);_maxweight.push_back(7.);
  // J/psi -> omega eta decay
  _incoming.push_back(443);_outgoingV.push_back(223);_outgoingP.push_back(221);
  _coupling.push_back(0.00145/GeV);_maxweight.push_back(7.);
  // J/psi -> gamma eta decay
  _incoming.push_back(443);_outgoingV.push_back(22);_outgoingP.push_back(221);
  _coupling.push_back(0.00095/GeV);_maxweight.push_back(2.);
  // J/psi -> phi eta decay
  _incoming.push_back(443);_outgoingV.push_back(333);_outgoingP.push_back(221);
  _coupling.push_back(0.00102/GeV);_maxweight.push_back(7.);
  // J/psi -> phi eta' decay
  _incoming.push_back(443);_outgoingV.push_back(333);_outgoingP.push_back(331);
  _coupling.push_back(0.00084/GeV);_maxweight.push_back(6.5);
  // J/psi -> omega pi0 decay
  _incoming.push_back(443);_outgoingV.push_back(223);_outgoingP.push_back(111);
  _coupling.push_back(0.00069/GeV);_maxweight.push_back(7.);
  // J/psi -> rho0 eta decay
  _incoming.push_back(443);_outgoingV.push_back(113);_outgoingP.push_back(221);
  _coupling.push_back(0.00054/GeV);_maxweight.push_back(4.);
  // J/psi -> rho0 eta' decay
  _incoming.push_back(443);_outgoingV.push_back(113);_outgoingP.push_back(331);
  _coupling.push_back(0.00045/GeV);_maxweight.push_back(3.5);
  // J/psi -> omega eta' decay
  _incoming.push_back(443);_outgoingV.push_back(223);_outgoingP.push_back(331);
  _coupling.push_back(0.00054/GeV);_maxweight.push_back(7.);
  // J/psi -> gamma pi0 decay
  _incoming.push_back(443);_outgoingV.push_back(22);_outgoingP.push_back(111);
  _coupling.push_back(0.00019/GeV);_maxweight.push_back(2.);
  // psi(2s)-> j/psi eta decay
  _incoming.push_back(100443);_outgoingV.push_back(443);_outgoingP.push_back(221);
  _coupling.push_back(0.213/GeV);_maxweight.push_back(2.);
  // psi(2s)-> gamma eta_c decay
  _incoming.push_back(100443);_outgoingV.push_back(22);_outgoingP.push_back(441);
  _coupling.push_back(0.0108/GeV);_maxweight.push_back(3.5);
  // psi(2s)-> gamma eta' decay
  _incoming.push_back(100443);_outgoingV.push_back(22);_outgoingP.push_back(331);
  _coupling.push_back(0.00057/GeV);_maxweight.push_back(2.);
  // psi(2s)-> J/psi pi0 decay
  _incoming.push_back(100443);_outgoingV.push_back(443);_outgoingP.push_back(111);
  _coupling.push_back(0.00846/GeV);_maxweight.push_back(2.);
  // initial size of the vectors for the database output
  _initsize=_incoming.size();
}

VectorMesonVectorPScalarDecayer::~VectorMesonVectorPScalarDecayer() {}

bool VectorMesonVectorPScalarDecayer::accept(const DecayMode & dm) const {
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
	{if((id1   ==_outgoingV[ix]&&id2   ==_outgoingP[ix])||
	    (id2   ==_outgoingV[ix]&&id1   ==_outgoingP[ix])){allowed=true;}}
      if(id0bar==_incoming[ix]&&!allowed)
	{if((id1bar==_outgoingV[ix]&&id2bar==_outgoingP[ix])||
	    (id2bar==_outgoingV[ix]&&id1bar==_outgoingP[ix])){allowed=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}

ParticleVector VectorMesonVectorPScalarDecayer::decay(const DecayMode & dm,
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
	{if((id1   ==_outgoingV[ix]&&id2   ==_outgoingP[ix])||
	    (id2   ==_outgoingV[ix]&&id1   ==_outgoingP[ix])){imode=ix;}}
      if(idbar==_incoming[ix])
	{if((id1bar==_outgoingV[ix]&&id2bar==_outgoingP[ix])||
	    (id2bar==_outgoingV[ix]&&id1bar==_outgoingP[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void VectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const 
{os << _incoming << _outgoingV << _outgoingP << _maxweight << _coupling;}

void VectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> _coupling;}

ClassDescription<VectorMesonVectorPScalarDecayer> VectorMesonVectorPScalarDecayer::initVectorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void VectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorPScalarDecayer> documentation
    ("The\\classname{VectorMesonVectorPScalarDecayer} class is designed for the "
     "decay of a vector meson to another vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1 particle",
     &VectorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,int> interfaceOutcomingPScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &VectorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonVectorPScalarDecayer::_coupling,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, -10000000, 10000000, false, false, true);

}

// the hadronic currents 
vector<LorentzPolarizationVector>  
VectorMesonVectorPScalarDecayer::decayCurrent(const bool vertex, const int, 
					      const Particle & inpart,
					      const ParticleVector & decay) const
{
  // storage for the current
  vector<LorentzPolarizationVector> temp;
  // check if the outgoing vector is a photon
  bool photon(_outgoingV[imode()]==ParticleID::gamma);

  // workaround for gcc 3.2.3 bug
  // set up the spin information for the decay products
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr mytemp = decay[1];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);

  VectorWaveFunction(temp,decay[0],outgoing,true,photon,vertex);
  for(unsigned int ix=0;ix<3;++ix)
    {
      if(ix==1&&photon){temp[ix]=LorentzPolarizationVector();}
      else
	{temp[ix]=_coupling[imode()]/inpart.mass()*
	    EpsFunction::product(inpart.momentum(),temp[ix],decay[0]->momentum());}
    }
  return temp;
}
 
bool VectorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
	  if(id1   ==_outgoingV[ix]&&id2   ==_outgoingP[ix]){imode=ix;order=true;}
	  if(id2   ==_outgoingV[ix]&&id1   ==_outgoingP[ix]){imode=ix;order=false;}
	}
      if(idbar==_incoming[ix]&&imode<0)
	{
	  if(id1bar==_outgoingV[ix]&&id2bar==_outgoingP[ix]){imode=ix;order=true;}
	  if(id2bar==_outgoingV[ix]&&id1bar==_outgoingP[ix]){imode=ix;order=false;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling = _coupling[imode]*dm.parent()->mass();  
  mecode = 1;
  return order;
}

// output the setup info for the particle database
void VectorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming " << ix << " "
		 << _incoming[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingVector " << ix << " "
		 << _outgoingV[ix] << "\n";
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
	  output << "insert " << fullName() << ":OutgoingVector " << ix << " "
		 << _outgoingV[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingPScalar " << ix << " "
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
