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
using Helicity::incoming;

void VectorMesonPVectorPScalarDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingA.size()||isize!=_outgoingP.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    {throw InitException() << "Inconsistent parameters in "
			   << "VectorMesonPVectorPScalarDecayer::doinit()" 
			   << Exception::abortnow;}
  // set up the integration channels
  vector<double> wgt(0);
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0]=getParticleData(_incoming[ix]);
      extpart[1]=getParticleData(_outgoingA[ix]);
      extpart[2]=getParticleData(_outgoingP[ix]);
      mode=new DecayPhaseSpaceMode(extpart,this);
      addMode(mode,_maxweight[ix],wgt);
    }
}

VectorMesonPVectorPScalarDecayer::VectorMesonPVectorPScalarDecayer() 
{
  // intermediates
  generateIntermediates(false);
  // reserve sizes of the vectors
  _incoming.reserve(25);_outgoingA.reserve(25);_outgoingP.reserve(25);
  _coupling.reserve(25);_maxweight.reserve(25);
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

int VectorMesonPVectorPScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const
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
	{if((id1   ==_outgoingA[ix]&&id2   ==_outgoingP[ix])||
	    (id2   ==_outgoingA[ix]&&id1   ==_outgoingP[ix])){imode=ix;}}
      if(idbar==_incoming[ix])
	{if((id1bar==_outgoingA[ix]&&id2bar==_outgoingP[ix])||
	    (id2bar==_outgoingA[ix]&&id1bar==_outgoingP[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}


void VectorMesonPVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const {
os << _incoming << _outgoingA << _outgoingP << _maxweight << ounit(_coupling,1/GeV);}

void VectorMesonPVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) {
is >> _incoming >> _outgoingA >> _outgoingP >> _maxweight >> iunit(_coupling,1/GeV);}

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
     1/GeV, 0, 0/GeV, 0./GeV, 100./GeV, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonPVectorPScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);
}

double VectorMesonPVectorPScalarDecayer::me2(bool vertex, const int,
					     const Particle & inpart,
					     const ParticleVector & decay) const
{
  // wavefunctions of the incoming particle
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // set up the spin information for the decay products
  vector<LorentzPolarizationVector> vout;
  VectorWaveFunction(vout,decay[0],outgoing,true,false,vertex);
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr myvout = decay[1];
  ScalarWaveFunction(myvout,outgoing,true,vertex);
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  Energy2 p0dotpv(inpart.momentum()*decay[0]->momentum());
  complex<Energy> epsdot(0.*MeV);
  InvEnergy2 pre(_coupling[imode()]/inpart.mass());  
  for(unsigned int ix=0;ix<3;++ix)
    {
      epsdot=vout[ix]*inpart.momentum();
      for(unsigned int iy=0;iy<3;++iy)
	{
	  newME(iy,ix,0)=pre*(p0dotpv*(vout[ix].dot(invec[iy]))-
			      epsdot*(invec[iy]*decay[0]->momentum()));
	}
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
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
  DecayIntegrator::dataBaseOutput(output,false);
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
		 << _coupling[ix]*MeV << "\n";
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
		 << _coupling[ix]*MeV << "\n";
	  output << "insert " << fullName() << ":MaxWeight " << ix << " "
		 << _maxweight[ix] << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
