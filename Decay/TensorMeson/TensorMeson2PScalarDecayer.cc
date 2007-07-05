// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMeson2PScalarDecayer class.
//

#include "TensorMeson2PScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMeson2PScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using namespace ThePEG::Helicity;

void TensorMeson2PScalarDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoing1.size()||isize!=_outgoing2.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    {throw InitException() << "Inconsistent parameters TensorMeson2PScalarDecayer" 
			   << Exception::abortnow;}
  // set up the integration channels
  vector<double> wgt(0);
  DecayPhaseSpaceModePtr mode;
  PDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0]=getParticleData(_incoming[ix]);
      extpart[1]=getParticleData(_outgoing1[ix]);
      extpart[2]=getParticleData(_outgoing2[ix]);
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      addMode(mode,_maxweight[ix],wgt);
    }
}

TensorMeson2PScalarDecayer::TensorMeson2PScalarDecayer() 
{
  // reserve size of vectors for speed
  _incoming.reserve(50);
  _outgoing1.reserve(50);
  _outgoing2.reserve(50);
  _coupling.reserve(50);
  _maxweight.reserve(50);
  // intermediates
  generateIntermediates(false);
  // a_2 -> eta pi
  _incoming.push_back(115);_outgoing1.push_back( 221);_outgoing2.push_back(111);
  _incoming.push_back(215);_outgoing1.push_back( 221);_outgoing2.push_back(211);
  _coupling.push_back(10.90/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(10.90/GeV);_maxweight.push_back(2.00);
  // a_2 -> eta' pi
  _incoming.push_back(115);_outgoing1.push_back( 331);_outgoing2.push_back(111);
  _incoming.push_back(215);_outgoing1.push_back( 331);_outgoing2.push_back(211);
  _coupling.push_back(9.92/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(9.92/GeV);_maxweight.push_back(2.1);
  // a_2 -> K K
  _incoming.push_back(115);_outgoing1.push_back( 311);_outgoing2.push_back(-311);
  _incoming.push_back(115);_outgoing1.push_back( 321);_outgoing2.push_back(-321);
  _incoming.push_back(215);_outgoing1.push_back( 321);_outgoing2.push_back(-311);
  _coupling.push_back(7.36/GeV);_maxweight.push_back(2.);
  _coupling.push_back(7.36/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(10.41/GeV);_maxweight.push_back(2.);
  // f_2 -> pi pi
  _incoming.push_back(225);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _incoming.push_back(225);_outgoing1.push_back( 111);_outgoing2.push_back(111);
  _coupling.push_back(18.73/GeV);_maxweight.push_back(2.);
  _coupling.push_back(13.24/GeV);_maxweight.push_back(2.02);
  // f_2 -> eta eta
  _incoming.push_back(225);_outgoing1.push_back( 221);_outgoing2.push_back(221);
  _coupling.push_back(8.854/GeV);_maxweight.push_back(2.15);
  // f_2 -> KK
  _incoming.push_back(225);_outgoing1.push_back( 321);_outgoing2.push_back(-321);
  _incoming.push_back(225);_outgoing1.push_back( 311);_outgoing2.push_back(-311);
  _coupling.push_back(11.03/GeV);_maxweight.push_back(2.);
  _coupling.push_back(11.38/GeV);_maxweight.push_back(2.);
  // f'_2 -> KK
  _incoming.push_back(335);_outgoing1.push_back( 321);_outgoing2.push_back(-321);
  _incoming.push_back(335);_outgoing1.push_back( 311);_outgoing2.push_back(-311);
  _coupling.push_back(14.65/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(14.65/GeV);_maxweight.push_back(2.);
  // f'_2 -> eta eta
  _incoming.push_back(335);_outgoing1.push_back( 221);_outgoing2.push_back(221);
  _coupling.push_back(9.15/GeV);_maxweight.push_back(2.02);
  // f'_2 -> pi pi 
  _incoming.push_back(335);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _incoming.push_back(335);_outgoing1.push_back( 111);_outgoing2.push_back(111);
  _coupling.push_back(0.860/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(0.608/GeV);_maxweight.push_back(2.02);
  // K_2 -> K eta
  _incoming.push_back( 325);_outgoing1.push_back( 321);_outgoing2.push_back(221);
  _incoming.push_back( 315);_outgoing1.push_back( 311);_outgoing2.push_back(221);
  _coupling.push_back(1.52/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(1.52/GeV);_maxweight.push_back(2.);
  // K_2 -> K pi
  _incoming.push_back( 325);_outgoing1.push_back( 321);_outgoing2.push_back( 111);
  _incoming.push_back( 325);_outgoing1.push_back( 311);_outgoing2.push_back( 211);
  _incoming.push_back( 315);_outgoing1.push_back( 311);_outgoing2.push_back( 111);
  _incoming.push_back( 315);_outgoing1.push_back( 321);_outgoing2.push_back(-211);
  _coupling.push_back(8.30/GeV) ;_maxweight.push_back(2.);
  _coupling.push_back(11.74/GeV);_maxweight.push_back(2.);
  _coupling.push_back(8.68/GeV) ;_maxweight.push_back(2.);
  _coupling.push_back(12.28/GeV);_maxweight.push_back(2.02);
  // B_2 -> B pi
  _incoming.push_back( 525);_outgoing1.push_back( 521);_outgoing2.push_back(111);
  _incoming.push_back( 525);_outgoing1.push_back( 511);_outgoing2.push_back(211);
  _incoming.push_back( 515);_outgoing1.push_back( 511);_outgoing2.push_back(111);
  _incoming.push_back( 515);_outgoing1.push_back( 521);_outgoing2.push_back(-211);
  _coupling.push_back(25.62/GeV);_maxweight.push_back(2.05);
  _coupling.push_back(36.24/GeV);_maxweight.push_back(2.);
  _coupling.push_back(25.62/GeV);_maxweight.push_back(2.05);
  _coupling.push_back(36.24/GeV);_maxweight.push_back(2.02);
  // D_2 -> D pi
  _incoming.push_back( 425);_outgoing1.push_back( 421);_outgoing2.push_back(111);
  _incoming.push_back( 425);_outgoing1.push_back( 411);_outgoing2.push_back(-211);
  _incoming.push_back( 415);_outgoing1.push_back( 411);_outgoing2.push_back(111);
  _incoming.push_back( 415);_outgoing1.push_back( 421);_outgoing2.push_back(211);
  _coupling.push_back(13.17/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(18.62/GeV);_maxweight.push_back(2.);
  _coupling.push_back(14.00/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(19.80/GeV);_maxweight.push_back(2.10);
  // D_s2
  _incoming.push_back( 435);_outgoing1.push_back( 421);_outgoing2.push_back( 321);
  _incoming.push_back( 435);_outgoing1.push_back( 411);_outgoing2.push_back( 311);
  _coupling.push_back(23.39/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(23.39/GeV);_maxweight.push_back(2.02);
  // B_s2
  _incoming.push_back( 535);_outgoing1.push_back( 521);_outgoing2.push_back(-321);
  _incoming.push_back( 535);_outgoing1.push_back( 511);_outgoing2.push_back(-311);
  _coupling.push_back(45.88/GeV);_maxweight.push_back(2.08);
  _coupling.push_back(45.88/GeV);_maxweight.push_back(2.);
  // chi_c2 to pi pi
  _incoming.push_back( 445);_outgoing1.push_back( 211);_outgoing2.push_back(-211);
  _incoming.push_back( 445);_outgoing1.push_back( 111);_outgoing2.push_back( 111);
  _coupling.push_back(0.0226/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(0.0159/GeV);_maxweight.push_back(1.62);
  // chi_c2 to K K
  _incoming.push_back( 445);_outgoing1.push_back( 321);_outgoing2.push_back(-321);
  _incoming.push_back( 445);_outgoing1.push_back( 311);_outgoing2.push_back(-311);
  _coupling.push_back(0.0187/GeV);_maxweight.push_back(2.02);
  _coupling.push_back(0.0187/GeV);_maxweight.push_back(2.02);
  // f_2 to sigma sigma
  _incoming.push_back(225);_outgoing1.push_back(9000221);_outgoing2.push_back(9000221);
  _coupling.push_back(99.27/GeV);_maxweight.push_back(30.);
  // pi_2 to sigma pi
  _incoming.push_back(10115);_outgoing1.push_back(9000221);_outgoing2.push_back(111);
  _incoming.push_back(10215);_outgoing1.push_back(9000221);_outgoing2.push_back(211);
  _coupling.push_back(99.27/GeV);_maxweight.push_back(30.);
  _coupling.push_back(99.27/GeV);_maxweight.push_back(30.);
  // eta_2 to a_0 pi
  _incoming.push_back(10225);_outgoing1.push_back(9000111);_outgoing2.push_back(111);
  _incoming.push_back(10225);_outgoing1.push_back(9000211);_outgoing2.push_back(-211);
  _coupling.push_back(99.27/GeV);_maxweight.push_back(30.);
  _coupling.push_back(99.27/GeV);_maxweight.push_back(30.);
  // eta'_2 to a_0 pi
  _incoming.push_back(10335);_outgoing1.push_back(9000111);_outgoing2.push_back(111);
  _incoming.push_back(10335);_outgoing1.push_back(9000211);_outgoing2.push_back(-211);
  _coupling.push_back(99.27/GeV);_maxweight.push_back(30.);
  _coupling.push_back(99.27/GeV);_maxweight.push_back(30.);
  // initial size of the vectors
  _initsize=_incoming.size();
}

TensorMeson2PScalarDecayer::~TensorMeson2PScalarDecayer() {}

int TensorMeson2PScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const
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
	{if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	    (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])){imode=ix;}}
      if(idbar==_incoming[ix])
	{if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	    (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}


void TensorMeson2PScalarDecayer::persistentOutput(PersistentOStream & os) const 
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << ounit(_coupling,1/GeV);}

void TensorMeson2PScalarDecayer::persistentInput(PersistentIStream & is, int) 
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> iunit(_coupling,1/GeV);}

ClassDescription<TensorMeson2PScalarDecayer> TensorMeson2PScalarDecayer::initTensorMeson2PScalarDecayer;
// Definition of the static class description member.

void TensorMeson2PScalarDecayer::Init() {

  static ClassDocumentation<TensorMeson2PScalarDecayer> documentation
    ("The TensorMeson2PScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMeson2PScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMeson2PScalarDecayer::_coupling,
     1/GeV, 0, 0/GeV, 0/GeV, 100./GeV, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMeson2PScalarDecayer::_maxweight,
     0, 0, 0, 0., 100000., false, false, true);
}


// matrix elememt for the process
double TensorMeson2PScalarDecayer::me2(bool vertex, const int,
				       const Particle & inpart,
				       const ParticleVector & decay) const
{
  vector<LorentzTensor<double> > inten;
  // wave functions etc for the incoming particle
  RhoDMatrix rhoin(PDT::Spin2);rhoin.average();
  TensorWaveFunction(inten,rhoin,const_ptr_cast<tPPtr>(&inpart),incoming,
		     true,false,vertex);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix)
    // workaround for gcc 3.2.3 bug
    //ALB {ScalarWaveFunction(decay[ix],outgoing,true,vertex);}
    {PPtr mytemp=decay[ix]; ScalarWaveFunction(mytemp,outgoing,true,vertex);}
  // calculate the matrix element
  DecayMatrixElement newME(PDT::Spin2,PDT::Spin0,PDT::Spin0);
  for(unsigned int ix=0;ix<5;++ix)
    {
      newME(ix,0,0) = _coupling[imode()]/inpart.mass()*
	((inten[ix]*decay[1]->momentum())*decay[0]->momentum());
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}

bool TensorMeson2PScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling=_coupling[imode]*dm.parent()->mass();
  mecode=7;
  return order;
}

void TensorMeson2PScalarDecayer::dataBaseOutput(ofstream & output,
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
	  output << "set " << fullName() << ":FirstOutgoing " << ix << " " 
		 << _outgoing1[ix] << "\n";
	  output << "set " << fullName() << ":SecondOutgoing " << ix << " " 
		 << _outgoing2[ix] << "\n";
	  output << "set " << fullName() << ":Coupling " << ix << " " 
		 << _coupling[ix]*GeV << "\n";
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
		 << _coupling[ix]*GeV << "\n";
	  output << "insert " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
