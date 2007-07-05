// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PVectorMesonVectorPScalarDecayer class.
//

#include "PVectorMesonVectorPScalarDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PVectorMesonVectorPScalarDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::VectorWaveFunction;
using Helicity::ScalarWaveFunction;
  using Helicity::incoming;
using Helicity::outgoing;

void PVectorMesonVectorPScalarDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingV.size()||isize!=_outgoingP.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    {throw InitException() << "Inconsistent parameters in "
			   << "PVectorMesonVectorPScalarDecayer::doinit()" 
			   << Exception::abortnow;}
  // set up the integration channels
  vector<double> wgt(0);
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {
      extpart[0]=getParticleData(_incoming[ix]);
      extpart[1]=getParticleData(_outgoingV[ix]);
      extpart[2]=getParticleData(_outgoingP[ix]);
      mode=new DecayPhaseSpaceMode(extpart,this);
      addMode(mode,_maxweight[ix],wgt);
    }
}

PVectorMesonVectorPScalarDecayer::PVectorMesonVectorPScalarDecayer() 
{
  // intermediates
  generateIntermediates(false);
  // reserve size of the vectors for speed
  _incoming.reserve(65);_outgoingV.reserve(65);_outgoingP.reserve(65);
  _coupling.reserve(65);_maxweight.reserve(65);
  // decay mode h'_1 to K K*
  _incoming.push_back( 10333);_outgoingV.push_back( 313);_outgoingP.push_back(-311);
  _incoming.push_back( 10333);_outgoingV.push_back( 323);_outgoingP.push_back(-321);
  _coupling.push_back(4.766/GeV);_maxweight.push_back(11.);
  _coupling.push_back(4.051/GeV);_maxweight.push_back(11.);
  // decay mode h_1 to rho
  _incoming.push_back(10223);_outgoingV.push_back( 113);_outgoingP.push_back( 111);
  _incoming.push_back(10223);_outgoingV.push_back(-213);_outgoingP.push_back( 211);
  _coupling.push_back(4.411/GeV);_maxweight.push_back(6.);
  _coupling.push_back(4.411/GeV);_maxweight.push_back(6.);
  // decay mode b_1 to omega pi
  _incoming.push_back( 10113);_outgoingV.push_back( 223);_outgoingP.push_back( 111);
  _incoming.push_back( 10213);_outgoingV.push_back( 223);_outgoingP.push_back( 211);
  _coupling.push_back(3.862/GeV);_maxweight.push_back(7.5);
  _coupling.push_back(3.862/GeV);_maxweight.push_back(7.5);
  // decay mode b_1^+ to gamma pi
  _incoming.push_back( 10213);_outgoingV.push_back(  22);_outgoingP.push_back( 211);
  _coupling.push_back(0.195/GeV);_maxweight.push_back(2.);
  // decay mode D'_s1 to D*K
  _incoming.push_back( 20433);_outgoingV.push_back( 413);_outgoingP.push_back( 311);
  _incoming.push_back( 20433);_outgoingV.push_back( 423);_outgoingP.push_back( 321);
  _coupling.push_back(0.161/GeV);_maxweight.push_back(2.5);
  _coupling.push_back(0.161/GeV);_maxweight.push_back(2.5);
  // decay mode B'_s1 to B*K
  _incoming.push_back( 20533);_outgoingV.push_back( 513);_outgoingP.push_back(-311);
  _incoming.push_back( 20533);_outgoingV.push_back( 523);_outgoingP.push_back(-321);
  _coupling.push_back(0.070/GeV);_maxweight.push_back(2.1);
  _coupling.push_back(0.070/GeV);_maxweight.push_back(2.1);
  // decay mode K_1 to rho K
  _incoming.push_back( 10323);_outgoingV.push_back( 213);_outgoingP.push_back( 311);
  _incoming.push_back( 10323);_outgoingV.push_back( 113);_outgoingP.push_back( 321);
  _incoming.push_back( 10313);_outgoingV.push_back(-213);_outgoingP.push_back( 321);
  _incoming.push_back( 10313);_outgoingV.push_back( 113);_outgoingP.push_back( 311);
  _coupling.push_back(4.36/GeV);_maxweight.push_back(3.5);
  _coupling.push_back(3.08/GeV);_maxweight.push_back(3.5);
  _coupling.push_back(4.36/GeV);_maxweight.push_back(3.5);
  _coupling.push_back(3.08/GeV);_maxweight.push_back(3.5);
  // decay mode K'_1 to rho K
  _incoming.push_back( 20323);_outgoingV.push_back( 213);_outgoingP.push_back( 311);
  _incoming.push_back( 20323);_outgoingV.push_back( 113);_outgoingP.push_back( 321);
  _incoming.push_back( 20313);_outgoingV.push_back(-213);_outgoingP.push_back( 321);
  _incoming.push_back( 20313);_outgoingV.push_back( 113);_outgoingP.push_back( 311);
  _coupling.push_back(0.834/GeV);_maxweight.push_back(6.5);
  _coupling.push_back(0.590/GeV);_maxweight.push_back(6.5);
  _coupling.push_back(0.834/GeV);_maxweight.push_back(6.5);
  _coupling.push_back(0.590/GeV);_maxweight.push_back(6.5);
  // decay mode K_1 to omega K
  _incoming.push_back( 10323);_outgoingV.push_back( 223);_outgoingP.push_back( 321);
  _incoming.push_back( 10313);_outgoingV.push_back( 223);_outgoingP.push_back( 311);
  _coupling.push_back(5.872/GeV);_maxweight.push_back(7.5);
  _coupling.push_back(7.906/GeV);_maxweight.push_back(7.5);
  // decay mode K'_1 to omega K
  _incoming.push_back( 20323);_outgoingV.push_back( 223);_outgoingP.push_back( 321);
  _incoming.push_back( 20313);_outgoingV.push_back( 223);_outgoingP.push_back( 311);
  _coupling.push_back(0.494/GeV);_maxweight.push_back(8.);
  _coupling.push_back(0.494/GeV);_maxweight.push_back(8.);
  // decay mode K_1 to K* pi
  _incoming.push_back( 10323);_outgoingP.push_back( 211);_outgoingV.push_back( 313);
  _incoming.push_back( 10323);_outgoingP.push_back( 111);_outgoingV.push_back( 323);
  _incoming.push_back( 10313);_outgoingP.push_back(-211);_outgoingV.push_back( 323);
  _incoming.push_back( 10313);_outgoingP.push_back( 111);_outgoingV.push_back( 313);
  _coupling.push_back(1.02/GeV);_maxweight.push_back(8.5);
  _coupling.push_back(0.72/GeV);_maxweight.push_back(8.5);
  _coupling.push_back(1.02/GeV);_maxweight.push_back(8.5);
  _coupling.push_back(0.72/GeV);_maxweight.push_back(8.5);
  // decay mode K'_1 to K* pi
  _incoming.push_back( 20323);_outgoingP.push_back( 211);_outgoingV.push_back( 313);
  _incoming.push_back( 20323);_outgoingP.push_back( 111);_outgoingV.push_back( 323);
  _incoming.push_back( 20313);_outgoingP.push_back(-211);_outgoingV.push_back( 323);
  _incoming.push_back( 20313);_outgoingP.push_back( 111);_outgoingV.push_back( 313);
  _coupling.push_back(2.83/GeV);_maxweight.push_back(12.);
  _coupling.push_back(2.01/GeV);_maxweight.push_back(12.);
  _coupling.push_back(2.83/GeV);_maxweight.push_back(12.);
  _coupling.push_back(2.01/GeV);_maxweight.push_back(12.);
  // decaymode D_1 to D* pi
  _incoming.push_back( 10423);_outgoingP.push_back(-211);_outgoingV.push_back( 413);
  _incoming.push_back( 10423);_outgoingP.push_back( 111);_outgoingV.push_back( 423);
  _incoming.push_back( 10413);_outgoingP.push_back( 211);_outgoingV.push_back( 423);
  _incoming.push_back( 10413);_outgoingP.push_back( 111);_outgoingV.push_back( 413);
  _coupling.push_back(0.473/GeV);_maxweight.push_back(3.);
  _coupling.push_back(0.333/GeV);_maxweight.push_back(3.);
  _coupling.push_back(0.570/GeV);_maxweight.push_back(3.);
  _coupling.push_back(0.403/GeV);_maxweight.push_back(3.);
  // decaymode D'_1 to D* pi
  _incoming.push_back( 20423);_outgoingP.push_back(-211);_outgoingV.push_back( 413);
  _incoming.push_back( 20423);_outgoingP.push_back( 111);_outgoingV.push_back( 423);
  _incoming.push_back( 20413);_outgoingP.push_back( 211);_outgoingV.push_back( 423);
  _incoming.push_back( 20413);_outgoingP.push_back( 111);_outgoingV.push_back( 413);
  _coupling.push_back(1.933/GeV);_maxweight.push_back(3.);
  _coupling.push_back(1.367/GeV);_maxweight.push_back(3.);
  _coupling.push_back(1.993/GeV);_maxweight.push_back(3.);
  _coupling.push_back(1.367/GeV);_maxweight.push_back(3.);
  // decaymode B_1 to B* pi
  _incoming.push_back( 10523);_outgoingP.push_back( 211);_outgoingV.push_back( 513);
  _incoming.push_back( 10523);_outgoingP.push_back( 111);_outgoingV.push_back( 523);
  _incoming.push_back( 10513);_outgoingP.push_back(-211);_outgoingV.push_back( 523);
  _incoming.push_back( 10513);_outgoingP.push_back( 111);_outgoingV.push_back( 513);
  _coupling.push_back(0.128/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(0.091/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(0.128/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(0.091/GeV);_maxweight.push_back(2.2);
  // decaymode B'_1 to B* pi
  _incoming.push_back( 20523);_outgoingP.push_back( 211);_outgoingV.push_back( 513);
  _incoming.push_back( 20523);_outgoingP.push_back( 111);_outgoingV.push_back( 523);
  _incoming.push_back( 20513);_outgoingP.push_back(-211);_outgoingV.push_back( 523);
  _incoming.push_back( 20513);_outgoingP.push_back( 111);_outgoingV.push_back( 513);
  _coupling.push_back(0.481/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(0.340/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(0.481/GeV);_maxweight.push_back(2.2);
  _coupling.push_back(0.340/GeV);_maxweight.push_back(2.2);
  // decaymode D_s1 to D* pi
  _incoming.push_back( 10433);_outgoingP.push_back( 111);_outgoingV.push_back( 433);
  _coupling.push_back(0.021/GeV);_maxweight.push_back(2.5);
  // decaymode D_s1 to D gamma
  _incoming.push_back( 10433);_outgoingP.push_back( 431);_outgoingV.push_back(22);
  _coupling.push_back(0.072/GeV);_maxweight.push_back(2.1);
  // decaymode B_s1 to B gamma
  _incoming.push_back( 10533);_outgoingP.push_back( 531);_outgoingV.push_back(22);
  _coupling.push_back(0.142/GeV);_maxweight.push_back(2.);
  // decaymode B_s1 to B* pi
  _incoming.push_back( 10533);_outgoingP.push_back( 111);_outgoingV.push_back( 533);
  _coupling.push_back(0.0074/GeV);_maxweight.push_back(2.1);
  // decaymode B_c1 to B_c gamma
  _incoming.push_back( 10543);_outgoingP.push_back( 541);_outgoingV.push_back(22);
  _coupling.push_back(0.072/GeV);_maxweight.push_back(2.2);
  // decaymode B'_c1 to B_c gamma
  _incoming.push_back( 20543);_outgoingP.push_back( 541);_outgoingV.push_back(22);
  _coupling.push_back(0.197/GeV);_maxweight.push_back(2.2);
  // decaymode h_c to eta_c gamma
  _incoming.push_back(10443);_outgoingP.push_back(441);_outgoingV.push_back(22);
  _coupling.push_back(0.364/GeV);_maxweight.push_back(4.);
  // decaymode h_b to eta_b gamma
  _incoming.push_back(10553);_outgoingP.push_back(551);_outgoingV.push_back(22);
  _coupling.push_back(0.0356/GeV);_maxweight.push_back(3.5);
  // a_1 to K* K
  _incoming.push_back(20213);_outgoingP.push_back(-311);_outgoingV.push_back(323);
  _incoming.push_back(20213);_outgoingP.push_back(321);_outgoingV.push_back(-313);
  _incoming.push_back(20113);_outgoingP.push_back(321);_outgoingV.push_back(-323);
  _incoming.push_back(20113);_outgoingP.push_back(311);_outgoingV.push_back(-313);
  _coupling.push_back(3.03/GeV);_maxweight.push_back(2.5);
  _coupling.push_back(3.03/GeV);_maxweight.push_back(2.5);
  _coupling.push_back(3.03/GeV);_maxweight.push_back(4.0);
  _coupling.push_back(3.03/GeV);_maxweight.push_back(4.0);
  // a_1 to gamma pi
  _incoming.push_back( 20113);_outgoingP.push_back( 111);_outgoingV.push_back(22);
  _incoming.push_back( 20213);_outgoingP.push_back( 211);_outgoingV.push_back(22);
  _coupling.push_back(0.01/GeV);_maxweight.push_back(2.);
  _coupling.push_back(0.01/GeV);_maxweight.push_back(2.);
  // f'_1 to K* K
  _incoming.push_back(20333);_outgoingP.push_back(321);_outgoingV.push_back(-323);
  _incoming.push_back(20333);_outgoingP.push_back(311);_outgoingV.push_back(-313);
  _coupling.push_back(1.637/GeV);_maxweight.push_back(7.);
  _coupling.push_back(1.737/GeV);_maxweight.push_back(7.);
  // decay mode K_1 to gamma K
  _incoming.push_back( 10313);_outgoingV.push_back( 22);_outgoingP.push_back( 311);
  _coupling.push_back(5.872/GeV);_maxweight.push_back(7.5);
  // decay mode K'_1 to gamma K
  _incoming.push_back( 20313);_outgoingV.push_back( 22);_outgoingP.push_back( 311);
  _coupling.push_back(0.494/GeV);_maxweight.push_back(8.);
  // initial size of the arrays
  _initsize = _coupling.size();
}

PVectorMesonVectorPScalarDecayer::~PVectorMesonVectorPScalarDecayer() {}

int PVectorMesonVectorPScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  ParticleMSet::const_iterator it;
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=(dm.parent()->CC())->id();}
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
	{if((id1   ==_outgoingV[ix]&&id2   ==_outgoingP[ix])||
	    (id2   ==_outgoingV[ix]&&id1   ==_outgoingP[ix])){imode=ix;}}
      if(idbar==_incoming[ix]&&imode<0)
	{if((id1bar==_outgoingV[ix]&&id2bar==_outgoingP[ix])||
	    (id2bar==_outgoingV[ix]&&id1bar==_outgoingP[ix])){imode=ix;cc=true;}}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}


void PVectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const
{os << _incoming << _outgoingV << _outgoingP << _maxweight << ounit(_coupling,1/GeV);}

void PVectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)
{is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> iunit(_coupling,1/GeV);}

ClassDescription<PVectorMesonVectorPScalarDecayer> PVectorMesonVectorPScalarDecayer::initPVectorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void PVectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<PVectorMesonVectorPScalarDecayer> documentation
    ("The PVectorMesonVectorPScalarDecayer class is designed for the "
     "decay of a pseudovector meson to a vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PVectorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1 particle",
     &PVectorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,int> interfaceOutcomingPScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &PVectorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PVectorMesonVectorPScalarDecayer::_coupling,
     1/GeV, 0, 0/GeV, 0./GeV, 100./GeV, false, false, true);

  static ParVector<PVectorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PVectorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, 0., 10000., false, false, true);

}


double PVectorMesonVectorPScalarDecayer::me2(bool vertex, const int,
					     const Particle & inpart,
					     const ParticleVector & decay) const
{
  // wavefunction for the decaying vector
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // is the vector massless
  bool photon(_outgoingV[imode()]==ParticleID::gamma);
  // set up the spin information for the decay products
  vector<LorentzPolarizationVector> vout;
  VectorWaveFunction(vout,decay[0],outgoing,true,photon,vertex);
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(decay[1],outgoing,true,vertex);
  PPtr mytemp = decay[1];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  Energy2 p0dotpv(inpart.momentum()*decay[0]->momentum());
  complex<Energy> epsdot(0.*MeV);
  InvEnergy2 pre(_coupling[imode()]/inpart.mass());
  for(unsigned ix=0;ix<3;++ix)
    {
      if(ix==1&&photon)
	{for(unsigned int iy=0;iy<3;++iy){newME(iy,ix,0)=0.;}}
      else
	{
	  epsdot=vout[ix]*inpart.momentum();
	  for(unsigned int iy=0;iy<3;++iy)
	    {newME(iy,ix,0)=pre*invec[iy].dot(p0dotpv*vout[ix]
					      -epsdot*decay[0]->momentum());}
	}
    }
  ME(newME);
  // return the answer
  return newME.contract(rhoin).real();
}

bool PVectorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,
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
  mecode = 4;
  return order;
}

void PVectorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
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
	  output << "set " << fullName() << ":OutgoingVector " << ix << " "
		 << _outgoingV[ix] << "\n";
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
	  output << "insert " << fullName() << ":OutgoingVector " << ix << " "
		 << _outgoingV[ix] << "\n";
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
