// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorVectorDecayer class.
//

#include "PScalarVectorVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PScalarVectorVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::EpsFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

PScalarVectorVectorDecayer::PScalarVectorVectorDecayer() 
{
  // decay eta -> omega gamma
  _incoming.push_back(331);_outgoing1.push_back(223);_outgoing2.push_back(22);
  _coupling.push_back(0.1412/GeV);_maxweight.push_back(1.1);
  // decay pi -> gamma gamma
  _incoming.push_back(111);_outgoing1.push_back(22);_outgoing2.push_back(22);
  _coupling.push_back(0.0178/GeV);_maxweight.push_back(1.);
  // decay eta -> gamma gamma
  _incoming.push_back(221);_outgoing1.push_back(22);_outgoing2.push_back(22);
  _coupling.push_back(0.0176/GeV);_maxweight.push_back(1.);
  // decay eta' -> gamma gamma
  _incoming.push_back(331);_outgoing1.push_back(22);_outgoing2.push_back(22);
  _coupling.push_back(0.0221/GeV);_maxweight.push_back(1.1);
  // decay eta_c -> rho rho
  _incoming.push_back(441);_outgoing1.push_back(213);_outgoing2.push_back(-213);
  _coupling.push_back(0.0494/GeV);_maxweight.push_back(2.5);
  _incoming.push_back(441);_outgoing1.push_back(113);_outgoing2.push_back( 113);
  _coupling.push_back(0.0349/GeV);_maxweight.push_back(2.5);
  // decay eta-c -> phi phi
  _incoming.push_back(441);_outgoing1.push_back(333);_outgoing2.push_back(333);
  _coupling.push_back(0.0215/GeV);_maxweight.push_back(6.);
  // decay eta-c -> gamma gamma
  _incoming.push_back(441);_outgoing1.push_back(22);_outgoing2.push_back(22);
  _coupling.push_back(0.00531/GeV);_maxweight.push_back(1.);
  // decay eta_c -> K* K*
  _incoming.push_back(441);_outgoing1.push_back(323);_outgoing2.push_back(-323);
  _coupling.push_back(0.0242/GeV);_maxweight.push_back(5.);
  _incoming.push_back(441);_outgoing1.push_back(313);_outgoing2.push_back(-313);
  _coupling.push_back(0.0242/GeV);_maxweight.push_back(5.);
  // initial size of the vectors
  _initsize = _incoming.size();
}

PScalarVectorVectorDecayer::~PScalarVectorVectorDecayer() {}

bool PScalarVectorVectorDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do
    {
      if(_incoming[ix]==id0)
	{if((_outgoing1[ix]==id1&&_outgoing2[ix]==id2)||
	    (_outgoing1[ix]==id2&&_outgoing2[ix]==id1)){allowed=true;}}
      ++ix;
    }
  while(!allowed&&ix<_incoming.size());
  return allowed;
}

ParticleVector PScalarVectorVectorDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  // workout which mode we are doing
  int imode(-1);
  int id(parent.id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  do
    {
      if(_incoming[ix]==id)
	{if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  // perform the decay
  return generate(false,false,imode,parent);
}


void PScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void PScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<PScalarVectorVectorDecayer> PScalarVectorVectorDecayer::initPScalarVectorVectorDecayer;
// Definition of the static class description member.

void PScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<PScalarVectorVectorDecayer> documentation
    ("The \\classname{PScalarVectorVectorDecayer} class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<PScalarVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &PScalarVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &PScalarVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarVectorVectorDecayer::_coupling,
     0, 0, 0, 0., 10000/GeV, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 200., false, false, true);
}

double PScalarVectorVectorDecayer::me2(bool vertex, const int ichan,
				   const Particle & inpart,
				   const ParticleVector & decay) const
{
  ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  // set up the spin info for the outgoing particles
  bool photon[2]={false,false};
  vector<LorentzPolarizationVector> wave[2];
  for(unsigned int ix=0;ix<2;++ix)
    {
      if(decay[ix]->id()==ParticleID::gamma){photon[ix]=true;}
      VectorWaveFunction(wave[ix],decay[ix],outgoing,true,photon[ix],vertex);
    }
  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1,PDT::Spin1);
  double fact(_coupling[imode()]/inpart.mass());
  unsigned int ix,iy;
  for(ix=0;ix<3;++ix)
    {
      for(iy=0;iy<3;++iy)
	{
	  newME(0,ix,iy)=fact*
	    EpsFunction::product(wave[0][ix],decay[1]->momentum(),wave[1][iy])
	    *decay[0]->momentum();
	}
    }
  ME(newME);
  RhoDMatrix rhoin(RhoDMatrix(PDT::Spin0));rhoin.average();
  return newME.contract(rhoin).real();
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
					       double & coupling) const
{
  int imode(-1);
  int id(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  do
    {
      if(_incoming[ix]==id)
	{if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	     (id2==_outgoing1[ix]&&id1==_outgoing2[ix])){imode=ix;}}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  coupling=_coupling[imode]*dm.parent()->mass();
  itype = 3;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}


// output the setup info for the particle database
void PScalarVectorVectorDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  output << "set " << fullName() << ":Iteration " << _niter << "\n";
  output << "set " << fullName() << ":Ntry " << _ntry << "\n";
  output << "set " << fullName() << ":Points " << _npoint << "\n";
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":Incoming   " << ix << " "
		 << _incoming[ix]   << "\n";
	  output << "set " << fullName() << ":FirstOutgoing  " << ix << " "
		 << _outgoing1[ix]  << "\n";
	  output << "set " << fullName() << ":SecondOutgoing " << ix << " "
		 << _outgoing2[ix]  << "\n";
	  output << "set " << fullName() << ":Coupling   " << ix << " "
		 << _coupling[ix]   << "\n";
	  output << "set " << fullName() << ":MaxWeight  " << ix << " "
		 << _maxweight[ix]  << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming   " << ix << " "
		 << _incoming[ix]   << "\n";
	  output << "insert " << fullName() << ":FirstOutgoing  " << ix << " "
		 << _outgoing1[ix]  << "\n";
	  output << "insert " << fullName() << ":SecondOutgoing " << ix << " "
		 << _outgoing2[ix]  << "\n";
	  output << "insert " << fullName() << ":Coupling   " << ix << " "
		 << _coupling[ix]   << "\n";
	  output << "insert " << fullName() << ":MaxWeight  " << ix << " "
		 << _maxweight[ix]  << "\n";
	}
    }
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}
