// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorVectorDecayer class.
//

#include "VectorMesonVectorVectorDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonVectorVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig{
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::VectorWaveFunction;
using Helicity::outgoing;

VectorMesonVectorVectorDecayer::VectorMesonVectorVectorDecayer() 
{
  // decay of rho'' to rho rho
  _incoming.push_back(30213);_outgoing1.push_back( 213);_outgoing2.push_back(113);
  _incoming.push_back(30113);_outgoing1.push_back(-213);_outgoing2.push_back(213);
  _coupling.push_back(3.26);_maxweight.push_back(25.);
  _coupling.push_back(3.26);_maxweight.push_back(25.);
  // decay of rho' to rho rho
  _incoming.push_back( 100213);_outgoing1.push_back( 213);_outgoing2.push_back(113);
  _incoming.push_back( 100113);_outgoing1.push_back(-213);_outgoing2.push_back(213);
  _coupling.push_back(14.63);_maxweight.push_back(50.);
  _coupling.push_back(14.63);_maxweight.push_back(50.);
  // initial size of the arrays
  _initsize=_incoming.size();
}

VectorMesonVectorVectorDecayer::~VectorMesonVectorVectorDecayer() {}

bool VectorMesonVectorVectorDecayer::accept(const DecayMode & dm) const {
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

ParticleVector VectorMesonVectorVectorDecayer::decay(const DecayMode & dm,
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


void VectorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;}

void VectorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;}

ClassDescription<VectorMesonVectorVectorDecayer> VectorMesonVectorVectorDecayer::initVectorMesonVectorVectorDecayer;
// Definition of the static class description member.

void VectorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorVectorDecayer> documentation
    ("The\\classname{VectorMesonVectorVectorDecayer} class is designed for the "
     "decay of a vector meson to two vector particles, either photons or other "
     "vector mesons.");

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceOutgoing1
    ("Outgoing1",
     "The PDG code for the first outgoing particle",
     &VectorMesonVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceOutgoing2
    ("Outgoing2",
     "The PDG code for the second outgoing particle",
     &VectorMesonVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonVectorVectorDecayer::_coupling,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// the hadronic currents
vector<LorentzPolarizationVector>  
VectorMesonVectorVectorDecayer::decayCurrent(const bool vertex, const int, 
					     const Particle & inpart,
					     const ParticleVector & decay) const
{
  // storage of the current
  vector<LorentzPolarizationVector> temp,eps[2];
  unsigned int ix,ipol1,ipol2;
  // set up the spin information for the decay products
  for(ix=0;ix<2;++ix)

    // workaround for gcc 3.2.3 bug
    //ALB {VectorWaveFunction(eps[ix],decay[ix],outgoing,true,
    //ALB			decay[ix]->id()==ParticleID::gamma,vertex);}
    {
      vector<LorentzPolarizationVector> mytemp; 
      VectorWaveFunction(mytemp,decay[ix],outgoing,true,
			decay[ix]->id()==ParticleID::gamma,vertex);
      eps[ix]=mytemp;
    }

  // work out the dot product we need for the current
  Complex p1p2((decay[0]->momentum())*(decay[1]->momentum()));
  Complex p1eps2[3],p2eps1[3];
  Complex eps1eps2;
  for(ix=0;ix<3;++ix)
    {
      p1eps2[ix]=eps[1][ix]*(decay[0]->momentum());
      p2eps1[ix]=eps[0][ix]*(decay[1]->momentum());
    }
  // now compute the current
  Lorentz5Momentum pdiff(decay[0]->momentum()-decay[1]->momentum());
  Complex m12(decay[0]->mass()*decay[0]->mass()),m22(decay[1]->mass()*decay[1]->mass());
  double fact(2.*_coupling[imode()]/(inpart.mass()*inpart.mass()*inpart.mass()));
  for(ipol1=0;ipol1<3;++ipol1)
    {
      for(ipol2=0;ipol2<3;++ipol2)
	{
	  eps1eps2=eps[0][ipol1]*eps[1][ipol2];
	  temp.push_back(fact*(p1eps2[ipol2]*p2eps1[ipol1]*pdiff
			       +p1eps2[ipol2]*m22*eps[0][ipol1]
			       -p2eps1[ipol1]*m12*eps[1][ipol2]
			       +eps1eps2*(-p1p2*pdiff
					  +m12*decay[1]->momentum()
					  -m22*decay[0]->momentum())));
	}
    }
  return temp;
} 

bool VectorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling = _coupling[imode]; 
  mecode = 5;
  return order; 
}

// output the setup information for the particle database
void VectorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output) const
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
	  output << "set " << fullName() << ":Outgoing1 " << ix << " " 
		 << _outgoing1[ix] << "\n";
	  output << "set " << fullName() << ":Outgoing2 " << ix << " " 
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
	  output << "insert " << fullName() << ":Outgoing1 " << ix << " " 
		 << _outgoing1[ix] << "\n";
	  output << "insert " << fullName() << ":Outgoing2 " << ix << " " 
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
