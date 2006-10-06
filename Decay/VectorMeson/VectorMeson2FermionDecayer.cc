// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2FermionDecayer class.
//

#include "VectorMeson2FermionDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMeson2FermionDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::VectorWaveFunction;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::incoming;
using Helicity::outgoing;

void VectorMeson2FermionDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoingf.size()||
     isize!=_outgoinga.size() || isize!=_maxweight.size())
    {throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "FermionDecayer::doiin() " << Exception::runerror;}
  // set up the integration channels
  vector<double> wgt(0);
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0]=getParticleData(_incoming[ix]);
      extpart[1]=getParticleData(_outgoingf[ix]);
      extpart[2]=getParticleData(_outgoinga[ix]);
      mode = new DecayPhaseSpaceMode(extpart,this);
      addMode(mode,_maxweight[ix],wgt);
    }
}

VectorMeson2FermionDecayer::VectorMeson2FermionDecayer() 
{
  // don't include intermediates
  generateIntermediates(false);
  // reserve size of vectors for speed
  _incoming.reserve(45);_outgoingf.reserve(45);_outgoinga.reserve(45);
  _coupling.reserve(45);_maxweight.reserve(45);
  // rho -> e+e-, mu+mu
  _incoming.push_back(113);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(113);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(18.524E-3);_coupling.push_back(18.35E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // omega -> e+e-, mu+mu-
  _incoming.push_back(223);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(223);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(5.387E-3);_coupling.push_back(5.387E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // phi -> e+e-, mu+mu-
  _incoming.push_back(333);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(333);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(6.852E-3);_coupling.push_back(6.852E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // psi(1d) to leptons
  _incoming.push_back(30443);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(30443);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _incoming.push_back(30443);_outgoingf.push_back(15);_outgoinga.push_back(-15); 
  _coupling.push_back(1.611E-3);_coupling.push_back(1.611E-3);
  _coupling.push_back(1.611E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // J/psi decay
  _incoming.push_back(443);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(443);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _coupling.push_back(8.088E-3);_coupling.push_back(8.088E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);
  // psi2s to leptons
  _incoming.push_back(100443);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(100443);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _incoming.push_back(100443);_outgoingf.push_back(15);_outgoinga.push_back(-15);
  _coupling.push_back(4.645E-3);_coupling.push_back(4.645E-3);
  _coupling.push_back(4.645E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // upsilon to leptons
  _incoming.push_back(553);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(553);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _incoming.push_back(553);_outgoingf.push_back(15);_outgoinga.push_back(-15);
  _coupling.push_back(2.290E-3);_coupling.push_back(2.290E-3);
  _coupling.push_back(2.290E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // upsilon 2s to leptons
  _incoming.push_back(100553);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(100553);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _incoming.push_back(100553);_outgoingf.push_back(15);_outgoinga.push_back(-15);
  _coupling.push_back(1.466E-3);_coupling.push_back(1.466E-3);
  _coupling.push_back(1.466E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // upsilon 3s to leptons
  _incoming.push_back(200553);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(200553);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _incoming.push_back(200553);_outgoingf.push_back(15);_outgoinga.push_back(-15);
  _coupling.push_back(1.316E-3);_coupling.push_back(1.316E-3);
  _coupling.push_back(1.316E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // upsilon 4s to leptons
  _incoming.push_back(300553);_outgoingf.push_back(11);_outgoinga.push_back(-11);
  _incoming.push_back(300553);_outgoingf.push_back(13);_outgoinga.push_back(-13);
  _incoming.push_back(300553);_outgoingf.push_back(15);_outgoinga.push_back(-15);
  _coupling.push_back(1.411E-3);_coupling.push_back(1.411E-3);
  _coupling.push_back(1.411E-3);
  _maxweight.push_back(2.);_maxweight.push_back(2.);_maxweight.push_back(2.);
  // baryonic jpsi decays
  // to neutrons and proton
  _incoming.push_back(443);_outgoingf.push_back(2212);_outgoinga.push_back(-2212);
  _incoming.push_back(443);_outgoingf.push_back(2112);_outgoinga.push_back(-2112);
  _maxweight.push_back(2.);_coupling.push_back(1.581E-3);
  _maxweight.push_back(2.);_coupling.push_back(1.581E-3);
  // to sigma's
  _incoming.push_back(443);_outgoingf.push_back(3112);_outgoinga.push_back(-3112);
  _incoming.push_back(443);_outgoingf.push_back(3212);_outgoinga.push_back(-3212);
  _incoming.push_back(443);_outgoingf.push_back(3222);_outgoinga.push_back(-3222);
  _maxweight.push_back(2.);_coupling.push_back(1.307E-3);
  _maxweight.push_back(2.);_coupling.push_back(1.307E-3);
  _maxweight.push_back(2.);_coupling.push_back(1.307E-3);
  // to Xi's 
  _incoming.push_back(443);_outgoingf.push_back(3322);_outgoinga.push_back(-3322);
  _incoming.push_back(443);_outgoingf.push_back(3312);_outgoinga.push_back(-3312);
  _maxweight.push_back(2.);_coupling.push_back(1.183E-3);
  _maxweight.push_back(2.);_coupling.push_back(1.183E-3);
  // to lambda
  _incoming.push_back(443);_outgoingf.push_back(3122);_outgoinga.push_back(-3122);
  _maxweight.push_back(2.);_coupling.push_back(1.284E-3);
  // baryonic psi(2s) decays
  // to neutrons and protons
  _incoming.push_back(100443);_outgoingf.push_back(2212);_outgoinga.push_back(-2212);
  _incoming.push_back(100443);_outgoingf.push_back(2112);_outgoinga.push_back(-2112);
  _maxweight.push_back(2.);_coupling.push_back(7.822E-4);
  _maxweight.push_back(2.);_coupling.push_back(7.822E-4);
  // to sigma's
  _incoming.push_back(100443);_outgoingf.push_back(3112);_outgoinga.push_back(-3112);
  _incoming.push_back(100443);_outgoingf.push_back(3212);_outgoinga.push_back(-3212);
  _incoming.push_back(100443);_outgoingf.push_back(3222);_outgoinga.push_back(-3222);
  _maxweight.push_back(2.);_coupling.push_back(6.120E-4);
  _maxweight.push_back(2.);_coupling.push_back(6.120E-4);
  _maxweight.push_back(2.);_coupling.push_back(6.120E-4);
  // to Xi's
  _incoming.push_back(100443);_outgoingf.push_back(3322);_outgoinga.push_back(-3322);
  _incoming.push_back(100443);_outgoingf.push_back(3312);_outgoinga.push_back(-3312);
  _maxweight.push_back(2.);_coupling.push_back(5.544E-4);
  _maxweight.push_back(2.);_coupling.push_back(5.544E-4);
  // to lambda
  _incoming.push_back(100443);_outgoingf.push_back(3122);_outgoinga.push_back(-3122);
  _maxweight.push_back(2.);_coupling.push_back(7.432E-4);  
  // set the initial size
  _initsize=_incoming.size();
}

VectorMeson2FermionDecayer::~VectorMeson2FermionDecayer() {}

int VectorMeson2FermionDecayer::modeNumber(bool & cc,const DecayMode & dm) const
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
      if(_incoming[ix]==id   )
	{if((id1   ==_outgoingf[ix]&&id2   ==_outgoinga[ix])||
	    (id2   ==_outgoingf[ix]&&id1   ==_outgoinga[ix])){imode=ix;}}
      if(_incoming[ix]==idbar)
	{if((id1bar==_outgoingf[ix]&&id2bar==_outgoinga[ix])||
	    (id2bar==_outgoingf[ix]&&id1bar==_outgoinga[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(imode<0&&ix<_incoming.size());
  return imode;
}

void VectorMeson2FermionDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingf << _outgoinga << _maxweight;
}

void VectorMeson2FermionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingf >> _outgoinga >> _maxweight;
}

ClassDescription<VectorMeson2FermionDecayer> VectorMeson2FermionDecayer::initVectorMeson2FermionDecayer;
// Definition of the static class description member.

void VectorMeson2FermionDecayer::Init() {

  static ClassDocumentation<VectorMeson2FermionDecayer> documentation
    ("The VectorMeson2FermionDecayer class is designed for the decay "
     "of vectro mesons to fermions. It is mainly used for the decay of vector mesons "
     "to electrons and muons.");

  static ParVector<VectorMeson2FermionDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson2FermionDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,int> interfaceOutcoming1
    ("OutgoingFermion",
     "The PDG code for the outgoing fermion",
     &VectorMeson2FermionDecayer::_outgoingf,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,int> interfaceOutcoming2
    ("OutgoingAntiFermion",
     "The PDG code for the second outgoing anti-fermion",
     &VectorMeson2FermionDecayer::_outgoinga,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMeson2FermionDecayer::_coupling,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2FermionDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMeson2FermionDecayer::_maxweight,
     0, 0, 0, -10000000, 10000000, false, false, true);

}

double VectorMeson2FermionDecayer::me2(bool vertex, const int,
				   const Particle & inpart,
				   const ParticleVector & decay) const
{
  // wavefunctions of the decaying particle
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(_outgoingf[imode()]!=decay[iferm]->id()){iferm=1;ianti=0;}
  // construct the spin information objects for the  decay products
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  SpinorBarWaveFunction(wavebar,decay[iferm],outgoing,true,vertex);
  SpinorWaveFunction(   wave   ,decay[ianti],outgoing,true,vertex);
  // prefactor
  double pre(_coupling[imode()]/inpart.mass());
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half);
  // now compute the currents
  LorentzPolarizationVector temp;
  for(unsigned ix=0;ix<2;++ix)
    {
      for(unsigned iy=0;iy<2;++iy)
	{
	  temp = pre*wave[ix].vectorCurrent(wavebar[iy]);
	  for(unsigned int iz=0;iz<3;++iz)
	    {
	      if(iferm>ianti){newME(iz,ix,iy)=invec[iz]*temp;}
	      else           {newME(iz,iy,ix)=invec[iz]*temp;}
	    }
	}
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}

bool VectorMeson2FermionDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
	  if(id1==_outgoingf[ix]&&id2==_outgoinga[ix]){imode=ix;order=true;}
	  if(id2==_outgoingf[ix]&&id1==_outgoinga[ix]){imode=ix;order=false;}
	}
      if(idbar==_incoming[ix]&&imode<0)
	{
	  if(id1bar==_outgoingf[ix]&&id2bar==_outgoinga[ix]){imode=ix;order=true;}
	  if(id2bar==_outgoingf[ix]&&id1bar==_outgoinga[ix]){imode=ix;order=false;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode];
  mecode=2;
  return order;
}


// output the setup information for the particle database
void VectorMeson2FermionDecayer::dataBaseOutput(ofstream & output,
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
	  output << "set " << fullName() << ":OutgoingFermion " << ix << " "
		 << _outgoingf[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingAntiFermion "  << ix << " "
		 << _outgoinga[ix] << "\n";
	  output << "set " << fullName() << ":Coupling " << ix << " "
		 << _coupling[ix] << "\n";
	  output << "set " << fullName() << ":MaxWeight " << ix << " "
		 << _maxweight[ix] << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming " << ix << " "
		 << _incoming[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingFermion "  << ix << " "
		 << _outgoingf[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingAntiFermion "  << ix << " "
		 << _outgoinga[ix] << "\n";
	  output << "insert " << fullName() << ":Coupling " << ix << " "
		 << _coupling[ix] << "\n";
	  output << "insert " << fullName() << ":MaxWeight " << ix << " "
		 << _maxweight[ix] << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}

