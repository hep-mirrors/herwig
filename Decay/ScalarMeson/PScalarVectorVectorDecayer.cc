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
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using ThePEG::Helicity::ScalarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

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
  // intermediates
  generateIntermediates(false);
}

PScalarVectorVectorDecayer::~PScalarVectorVectorDecayer() {}

int PScalarVectorVectorDecayer::modeNumber(bool &,const DecayMode & dm) const
{
  int imode(-1);
  if(dm.products().size()!=2){return imode;}
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
  return imode;
}

void PScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/GeV) 
     << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void PScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/GeV)
     >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<PScalarVectorVectorDecayer> PScalarVectorVectorDecayer::initPScalarVectorVectorDecayer;
// Definition of the static class description member.

void PScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<PScalarVectorVectorDecayer> documentation
    ("The PScalarVectorVectorDecayer class is designed for"
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
     1/GeV, 0, 0/GeV, 0./GeV, 10000/GeV, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 200., false, false, true);
}

double PScalarVectorVectorDecayer::me2(bool vertex, const int,
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
  InvEnergy2 fact(_coupling[imode()]/inpart.mass());
  unsigned int ix,iy;
  for(ix=0;ix<3;++ix)
    {
      for(iy=0;iy<3;++iy)
	{
	  newME(0,ix,iy)=fact*
	    Helicity::
	    epsilon(wave[0][ix],decay[1]->momentum(),wave[1][iy])
	    *decay[0]->momentum();
	}
    }
  ME(newME);
  RhoDMatrix rhoin(PDT::Spin0); rhoin.average();
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
void PScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
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
		 << _coupling[ix]*MeV   << "\n";
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
		 << _coupling[ix]*MeV   << "\n";
	  output << "insert " << fullName() << ":MaxWeight  " << ix << " "
		 << _maxweight[ix]  << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
