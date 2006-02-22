// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorVectorDecayer class.
//

#include "TensorMesonVectorVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TensorMesonVectorVectorDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig{
using namespace ThePEG;
using namespace ThePEG::Helicity;
using namespace Herwig::Helicity;

void TensorMesonVectorVectorDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoing1.size()||isize!=_outgoing2.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    {throw InitException() << "Inconsistent parameters TensorMesonVectorVectorDecayer" 
			   << Exception::abortnow;}
  // set up the integration channels
  vector<double> wgt(0);
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0]=getParticleData(_incoming[ix]);
      extpart[1]=getParticleData(_outgoing1[ix]);
      extpart[2]=getParticleData(_outgoing2[ix]);
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      addMode(mode,_maxweight[ix],wgt);
    }
}

TensorMesonVectorVectorDecayer::TensorMesonVectorVectorDecayer() 
{
  // intermediates
  generateIntermediates(false);
  // size of the vectors to improve speed
  _incoming.reserve(25);_outgoing1.reserve(25);_outgoing2.reserve(25);
  _coupling.reserve(25);_maxweight.reserve(25);
  // a_2 -> gamma gamma
  _incoming.push_back(115);_outgoing1.push_back( 22);_outgoing2.push_back(22);
  _coupling.push_back(0.02336/GeV);_maxweight.push_back(2.);
  // f_2 -> gamma gamma
  _incoming.push_back(225);_outgoing1.push_back( 22);_outgoing2.push_back(22);
  _coupling.push_back(0.01253/GeV);_maxweight.push_back(2.01);
  // f'_2 -> gamma gamma
  _incoming.push_back(335);_outgoing1.push_back( 22);_outgoing2.push_back(22);
  _coupling.push_back(0.00161/GeV);_maxweight.push_back(2.01);
  // chi_b(2P) decays
  _incoming.push_back(100555);_outgoing1.push_back(553);_outgoing2.push_back(223);
  _coupling.push_back(0.0118/GeV);_maxweight.push_back(2.25);
  _incoming.push_back(100555);_outgoing1.push_back(553);_outgoing2.push_back(22);
  _coupling.push_back(0.0172/GeV);_maxweight.push_back(2.02);
  _incoming.push_back(100555);_outgoing1.push_back(100553);_outgoing2.push_back(22);
  _coupling.push_back(0.145/GeV);_maxweight.push_back(2.02);
  _incoming.push_back(100555);_outgoing1.push_back(333);_outgoing2.push_back(333);
  _coupling.push_back(0.00483/GeV);_maxweight.push_back(13.5);
  // chi_c decays
  _incoming.push_back(445);_outgoing1.push_back(443);_outgoing2.push_back(22);
  _coupling.push_back(0.243/GeV);_maxweight.push_back(2.02);
  _incoming.push_back(445);_outgoing1.push_back(323);_outgoing2.push_back(-323);
  _coupling.push_back(0.00627/GeV);_maxweight.push_back(21.);
  _incoming.push_back(445);_outgoing1.push_back(313);_outgoing2.push_back(-313);
  _coupling.push_back(0.00627/GeV);_maxweight.push_back(19.);
  _incoming.push_back(445);_outgoing1.push_back(333);_outgoing2.push_back(333);
  _coupling.push_back(0.00475/GeV);_maxweight.push_back(19.);
  _incoming.push_back(445);_outgoing1.push_back(22);_outgoing2.push_back(22);
  _coupling.push_back(0.00120/GeV);_maxweight.push_back(2.);
  // chi_b(1P) decays
  _incoming.push_back(555);_outgoing1.push_back(553);_outgoing2.push_back(22);
  _coupling.push_back(0.0683/GeV);_maxweight.push_back(2.02);
  // a_2 omega rho 
  _incoming.push_back( 115);_outgoing1.push_back( 223);_outgoing2.push_back( 113);
  _incoming.push_back( 215);_outgoing1.push_back( 223);_outgoing2.push_back( 213);
  _coupling.push_back(22.2/GeV);_maxweight.push_back(51.);
  _coupling.push_back(22.2/GeV);_maxweight.push_back(51.);
  // f_2 rho rho
  _incoming.push_back( 225);_outgoing1.push_back( 113);_outgoing2.push_back( 113);
  _incoming.push_back( 225);_outgoing1.push_back( 213);_outgoing2.push_back(-213);
  _coupling.push_back(10.5/GeV);_maxweight.push_back(40.);
  _coupling.push_back(14.8/GeV);_maxweight.push_back(40.);
  // K_2-> K* rho
  _incoming.push_back( 315);_outgoing1.push_back( 113);_outgoing2.push_back( 313);
  _incoming.push_back( 315);_outgoing1.push_back(-213);_outgoing2.push_back( 323);
  _incoming.push_back( 325);_outgoing1.push_back( 113);_outgoing2.push_back( 323);
  _incoming.push_back( 325);_outgoing1.push_back( 213);_outgoing2.push_back( 313);
  _coupling.push_back(13.14/GeV);_maxweight.push_back(34.);
  _coupling.push_back(18.58/GeV);_maxweight.push_back(34.);
  _coupling.push_back(13.14/GeV);_maxweight.push_back(34.);
  _coupling.push_back(18.58/GeV);_maxweight.push_back(34.);
  // initial size of the vectors
  _initsize = _incoming.size();
}

TensorMesonVectorVectorDecayer::~TensorMesonVectorVectorDecayer() {}

int TensorMesonVectorVectorDecayer::modeNumber(bool & cc,const DecayMode & dm) const
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


void TensorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const 
{os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;}

void TensorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;}

ClassDescription<TensorMesonVectorVectorDecayer> TensorMesonVectorVectorDecayer::initTensorMesonVectorVectorDecayer;
// Definition of the static class description member.

void TensorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorVectorDecayer> documentation
    ("The TensorMesonVectorVectorDecayer class performs the"
     " decay of a tensor meson to two scalar mesons.");

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMesonVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &TensorMesonVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &TensorMesonVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMesonVectorVectorDecayer::_coupling,
     0, 0, 0, 0./GeV, 1000./GeV, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMesonVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// matrix elememt for the process
double TensorMesonVectorVectorDecayer::me2(bool vertex, const int ichan,
					   const Particle & inpart,
					   const ParticleVector & decay) const
{
  // wave functions etc for the incoming particle
  vector<LorentzTensor> inten;
  RhoDMatrix rhoin(PDT::Spin2);rhoin.average();
  TensorWaveFunction(inten,rhoin,const_ptr_cast<tPPtr>(&inpart),incoming,
		     true,false,vertex);
  // wave functions for the decay products
  vector<LorentzPolarizationVector> pol[2];
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    {
      photon[ix]=decay[ix]->id()==ParticleID::gamma;

      // workaround for gcc 3.2.3 bug
      //ALB VectorWaveFunction(pol[ix],decay[ix],outgoing,true,photon[ix],vertex);
      vector<LorentzPolarizationVector> mytempLPV;
      VectorWaveFunction(mytempLPV,decay[ix],outgoing,true,photon[ix],vertex);
      pol[ix]=mytempLPV;
    }
  // compute some useful dot products etc
  Complex p1eps2[3],p2eps1[3];
  Energy2 p1p2(decay[0]->momentum()*decay[1]->momentum());
  for(unsigned int ix=0;ix<3;++ix)
    {
      p1eps2[ix]=pol[1][ix]*decay[0]->momentum();
      p2eps1[ix]=pol[0][ix]*decay[1]->momentum();
    }
  // compute the traces and useful dot products
  Complex trace[5],pboth[5];
  LorentzPolarizationVector pleft[2][5],pright[2][5];
  for(unsigned int ix=0;ix<5;++ix)
    {
      trace[ix]=inten[ix].trace();
      pleft[0][ix] =(-decay[0]->momentum())*inten[ix];
      pleft[1][ix] =(-decay[1]->momentum())*inten[ix];
      pright[0][ix]=inten[ix]*(-decay[0]->momentum());
      pright[1][ix]=inten[ix]*(-decay[1]->momentum());
      pboth[ix]=((pleft[0][ix]+pright[0][ix])*decay[1]->momentum());
    }
  // loop to compute the matrix element
  DecayMatrixElement newME(PDT::Spin2,PDT::Spin1,PDT::Spin1);
  Complex e1e2;
  LorentzTensor te1e2;
  double fact(_coupling[imode()]/inpart.mass());
  Complex me;
  for(unsigned int ix=0;ix<3;++ix)
    {
      for(unsigned iy=0;iy<3;++iy)
	{
	  e1e2=pol[0][ix]*pol[1][iy];
	  te1e2=p1p2*(LorentzTensor(pol[0][ix],pol[1][iy])
		       +LorentzTensor(pol[1][iy],pol[0][ix]));
	  for(unsigned int inhel=0;inhel<5;++inhel)
	    {
	      me = inten[inhel]*te1e2
		-p2eps1[ix]*(pol[1][iy]*(pleft[0][inhel]+pright[0][inhel])) 
		-p1eps2[iy]*(pol[0][ix]*(pleft[1][inhel]+pright[1][inhel]))
		+pboth[inhel]*e1e2
		+(p2eps1[ix]*p1eps2[iy]-e1e2*p1p2)*trace[inhel];
	      newME(inhel,ix,iy)=fact*me;
	    }


	}
    }
  ME(newME);
  /*
  // testing the matrix element
  Energy2 m02=inpart.mass()*inpart.mass();
  Energy2 m12=decay[0]->mass()*decay[0]->mass();
  Energy2 m22=decay[1]->mass()*decay[1]->mass();
  double root = inpart.mass()*inpart.mass()*(inpart.mass()*inpart.mass()
  -2.*decay[0]->mass()*decay[0]->mass()
  -2.*decay[1]->mass()*decay[1]->mass())
  +decay[0]->mass()*decay[0]->mass()*(decay[0]->mass()*decay[0]->mass()
  -2.*decay[1]->mass()*decay[1]->mass())
  +decay[1]->mass()*decay[1]->mass()*decay[1]->mass()*decay[1]->mass();
  Energy pcm =0.5*sqrt(root)/inpart.mass();
  Energy2 pcm2=pcm*pcm;
  cout << "testing the matrix element VV " 
  << inpart.PDGName() << " -> " 
  << decay[0]->PDGName() << "  " 
  << decay[1]->PDGName() << "  " 
  << pcm/30./pi/m02/m02*_coupling[imode()]*_coupling[imode()]*
  (3.*m02*(8.*pcm2*pcm2+5.*(m12*m22+pcm2*(m12+m22)))
  -5.*(m12-m22)*(m12-m22)*pcm2) << endl;
  */
  // return the answer
  return newME.contract(rhoin).real();
}

bool TensorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling=_coupling[imode]*dm.parent()->mass();
  mecode=9;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}
void TensorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
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
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
