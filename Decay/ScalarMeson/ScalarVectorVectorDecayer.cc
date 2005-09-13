// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarVectorVectorDecayer class.
//

#include "ScalarVectorVectorDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarVectorVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::VectorWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;


inline ScalarVectorVectorDecayer::ScalarVectorVectorDecayer() 
{
  // f_0(1370) to rho rho
  _incoming.push_back(10221);_outgoing1.push_back(113);_outgoing2.push_back(113);
  _incoming.push_back(10221);_outgoing1.push_back(213);_outgoing2.push_back(-213);
  _coupling.push_back(11.26/GeV);_maxweight.push_back(20.);
  _coupling.push_back(15.92/GeV);_maxweight.push_back(20.);
  // f_0(1500) to rho rho
  _incoming.push_back(9030221);_outgoing1.push_back(113);_outgoing2.push_back(113);
  _incoming.push_back(9030221);_outgoing1.push_back(213);_outgoing2.push_back(-213);
  _coupling.push_back(1.691/GeV);_maxweight.push_back(20.);
  _coupling.push_back(2.391/GeV);_maxweight.push_back(20.);
  // chi_c0 decays
  _incoming.push_back(10441);_outgoing1.push_back(443);_outgoing2.push_back(22);
  _coupling.push_back(1./GeV);_maxweight.push_back(1.);
  _incoming.push_back(10441);_outgoing1.push_back(323);_outgoing2.push_back(-323);
  _incoming.push_back(10441);_outgoing1.push_back(313);_outgoing2.push_back(-313);
  _coupling.push_back(1./GeV);_maxweight.push_back(1.);
  _coupling.push_back(1./GeV);_maxweight.push_back(1.);
  _incoming.push_back(10441);_outgoing1.push_back(333);_outgoing2.push_back(333);
  _coupling.push_back(1./GeV);_maxweight.push_back(1.);
  _incoming.push_back(10441);_outgoing1.push_back(22);_outgoing2.push_back(22);
  _coupling.push_back(1./GeV);_maxweight.push_back(1.);
}

ScalarVectorVectorDecayer::~ScalarVectorVectorDecayer() {}

bool ScalarVectorVectorDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // check that at least some modes exist
  if(_incoming.size()==0){return false;}
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

ParticleVector ScalarVectorVectorDecayer::decay(const DecayMode & dm,
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


void ScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void ScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<ScalarVectorVectorDecayer> ScalarVectorVectorDecayer::initScalarVectorVectorDecayer;
// Definition of the static class description member.

void ScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<ScalarVectorVectorDecayer> documentation
    ("The \\classname{ScalarVectorVectorDecayer} class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<ScalarVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &ScalarVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &ScalarVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarVectorVectorDecayer::_coupling,
     0, 0, 0, 0., 10000/GeV, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 500000., false, false, true);

}

double ScalarVectorVectorDecayer::me2(bool vertex, const int ichan,
				      const Particle & inpart,
				      const ParticleVector & decay) const
{
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // set up the spin info for the outgoing particles
  bool photon[2]={false,false};
  vector<LorentzPolarizationVector> wave[2];
  for(unsigned int ix=0;ix<2;++ix)
    {
      if(decay[ix]->id()==ParticleID::gamma){photon[ix]=true;}

      // workaround for gcc 3.2.3 bug
      //ALB VectorWaveFunction(wave[ix],decay[ix],outgoing,true,photon[ix],vertex);
      vector<LorentzPolarizationVector> mytempLPV; 
      VectorWaveFunction(mytempLPV,decay[ix],outgoing,true,photon[ix],vertex);
      wave[ix]=mytempLPV;
    }
  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin1,PDT::Spin1);
  double fact(_coupling[imode()]/inpart.mass());
  Energy p1p2(decay[0]->momentum()*decay[1]->momentum());
  unsigned int ix,iy;
  for(ix=0;ix<3;++ix)
    {for(iy=0;iy<3;++iy)
	{newME(0,ix,iy)=fact*(p1p2*wave[0][ix]*wave[1][iy]-
			      (wave[1][iy]*decay[0]->momentum())*
			      (wave[0][ix]*decay[1]->momentum()));}}
  ME(newME);
  RhoDMatrix rhoin(PDT::Spin0); rhoin.average();
  return newME.contract(rhoin).real();
}

// output the setup info for the particle database
void ScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
					       bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
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
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
  itype = 12;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}

}
