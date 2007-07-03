// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarScalarScalarDecayer class.
//

#include "ScalarScalarScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

ScalarScalarScalarDecayer::ScalarScalarScalarDecayer() {
  // f_0(980) to pi pi
  _incoming.push_back(9010221);_outgoing1.push_back(111);_outgoing2.push_back( 111);
  _incoming.push_back(9010221);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _coupling.push_back(2.093*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(2.961*GeV);_maxweight.push_back(1.05);
  // f_0(980) to K K 
  _incoming.push_back(9010221);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _incoming.push_back(9010221);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _coupling.push_back(5.921*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(5.921*GeV);_maxweight.push_back(1.05);
  // f_0(1370) to pi pi
  _incoming.push_back(10221);_outgoing1.push_back(111);_outgoing2.push_back(111);
  _incoming.push_back(10221);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _coupling.push_back(0.745*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(1.054*GeV);_maxweight.push_back(1.05);
  // f_0(1370) to pi' pi
  _incoming.push_back(10221);_outgoing1.push_back(100111);_outgoing2.push_back(111);
  _incoming.push_back(10221);_outgoing1.push_back(100211);_outgoing2.push_back(-211);
  _coupling.push_back(5.027*GeV);_maxweight.push_back(2.1);
  _coupling.push_back(5.027*GeV);_maxweight.push_back(2.1);
  // f_0(1370) to K K 
  _incoming.push_back(10221);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _incoming.push_back(10221);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _coupling.push_back(0.886*GeV);_maxweight.push_back(1.1);
  _coupling.push_back(0.886*GeV);_maxweight.push_back(1.1);
  // f_0(1710) to pi pi
  _incoming.push_back(10331);_outgoing1.push_back(111);_outgoing2.push_back(111);
  _incoming.push_back(10331);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _coupling.push_back(0.696*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(0.985*GeV);_maxweight.push_back(1.05);
  // f_0(1710) to K K 
  _incoming.push_back(10331);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _incoming.push_back(10331);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _coupling.push_back(2.096*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(2.096*GeV);_maxweight.push_back(1.05);
  // sigma to pi pi
  _incoming.push_back(9000221);_outgoing1.push_back(111);_outgoing2.push_back(111);
  _incoming.push_back(9000221);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _coupling.push_back(3.654*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(5.178*GeV);_maxweight.push_back(1.05);
  // a_0 to eta pi
  _incoming.push_back( 9000111);_outgoing1.push_back(221);_outgoing2.push_back( 111);
  _incoming.push_back( 9000211);_outgoing1.push_back(221);_outgoing2.push_back( 211);
  _coupling.push_back(5.307*GeV);_maxweight.push_back(1.1);
  _coupling.push_back(5.307*GeV);_maxweight.push_back(1.1);
  // a_0 to K K
  _incoming.push_back( 9000111);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _incoming.push_back( 9000111);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _incoming.push_back( 9000211);_outgoing1.push_back(321);_outgoing2.push_back(-311);
  _coupling.push_back(2.242*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(2.242*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(3.171*GeV);_maxweight.push_back(1.05);
  // a'_0 to eta pi
  _incoming.push_back( 10111);_outgoing1.push_back(221);_outgoing2.push_back( 111);
  _incoming.push_back( 10211);_outgoing1.push_back(221);_outgoing2.push_back( 211);
  _coupling.push_back(1.357*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(1.357*GeV);_maxweight.push_back(1.05);
  // a'_0 to eta' pi
  _incoming.push_back( 10111);_outgoing1.push_back(331);_outgoing2.push_back( 111);
  _incoming.push_back( 10211);_outgoing1.push_back(331);_outgoing2.push_back( 211);
  _coupling.push_back(0.995*GeV);_maxweight.push_back(1.6);
  _coupling.push_back(0.995*GeV);_maxweight.push_back(1.6);
  // a'_0 to K K
  _incoming.push_back( 10111);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _incoming.push_back( 10111);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _incoming.push_back( 10211);_outgoing1.push_back(321);_outgoing2.push_back(-311);
  _coupling.push_back(0.950*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(0.950*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(1.344*GeV);_maxweight.push_back(1.05);
  // f_0(1370) to eta eta  
  _incoming.push_back( 10221);_outgoing1.push_back(221);_outgoing2.push_back(221);
  _coupling.push_back(0.235*GeV);_maxweight.push_back(1.1);
  // f_0(1710) to eta eta  
  _incoming.push_back( 10331);_outgoing1.push_back(221);_outgoing2.push_back(221);
  _coupling.push_back(2.189*GeV);_maxweight.push_back(1.2);
  // f_0(1370) to sigma sigma  
  _incoming.push_back( 10221);_outgoing1.push_back(9000221);
  _outgoing2.push_back(9000221);
  _coupling.push_back(21.46*GeV);_maxweight.push_back(7.);
  // K_0* to K pi
  _incoming.push_back( 10311);_outgoing1.push_back( 311);_outgoing2.push_back( 111);
  _incoming.push_back( 10311);_outgoing1.push_back( 321);_outgoing2.push_back(-211);
  _incoming.push_back( 10321);_outgoing1.push_back( 321);_outgoing2.push_back( 111);
  _incoming.push_back( 10321);_outgoing1.push_back( 311);_outgoing2.push_back( 211);
  _coupling.push_back(2.837*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(4.000*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(2.837*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(4.000*GeV);_maxweight.push_back(1.05);
  // D_0* to D pi
  _incoming.push_back( 10411);_outgoing1.push_back( 411);_outgoing2.push_back( 111);
  _incoming.push_back( 10411);_outgoing1.push_back( 421);_outgoing2.push_back( 211);
  _incoming.push_back( 10421);_outgoing1.push_back( 421);_outgoing2.push_back( 111);
  _incoming.push_back( 10421);_outgoing1.push_back( 411);_outgoing2.push_back(-211);
  _coupling.push_back(5.408*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(7.623*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(5.408*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(7.623*GeV);_maxweight.push_back(1.05);
  // B_0* to B pi
  _incoming.push_back( 10511);_outgoing1.push_back( 511);_outgoing2.push_back( 111);
  _incoming.push_back( 10511);_outgoing1.push_back( 521);_outgoing2.push_back(-211);
  _incoming.push_back( 10521);_outgoing1.push_back( 521);_outgoing2.push_back( 111);
  _incoming.push_back( 10521);_outgoing1.push_back( 511);_outgoing2.push_back( 211);
  _coupling.push_back(9.327*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(13.19*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(9.327*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(13.19*GeV);_maxweight.push_back(1.05);
  // K' to K_0* pi
  _incoming.push_back( 100311);_outgoing1.push_back( 10311);_outgoing2.push_back( 111);
  _incoming.push_back( 100311);_outgoing1.push_back( 10321);_outgoing2.push_back(-211);
  _incoming.push_back( 100321);_outgoing1.push_back( 10321);_outgoing2.push_back( 111);
  _incoming.push_back( 100321);_outgoing1.push_back( 10311);_outgoing2.push_back( 211);
  _coupling.push_back(6.822*GeV);_maxweight.push_back(2.);
  _coupling.push_back(9.770*GeV);_maxweight.push_back(2.);
  _coupling.push_back(6.822*GeV);_maxweight.push_back(2.);
  _coupling.push_back(9.770*GeV);_maxweight.push_back(2.);
  // D_s0* to D_s pi
  _incoming.push_back( 10431);_outgoing1.push_back( 431);_outgoing2.push_back( 111);
  _coupling.push_back(0.103*GeV);_maxweight.push_back(1.05);
  // B_s0* to B_s pi
  _incoming.push_back( 10531);_outgoing1.push_back( 531);_outgoing2.push_back( 111);
  _coupling.push_back(8.314*GeV);_maxweight.push_back(1.05);
  // eta'' to a_0 pi
  _incoming.push_back(100221);_outgoing1.push_back(9000111);_outgoing2.push_back(111);
  _incoming.push_back(100221);_outgoing1.push_back(9000211);_outgoing2.push_back(-211);
  _coupling.push_back(0.936*GeV);_maxweight.push_back(2.);
  _coupling.push_back(0.936*GeV);_maxweight.push_back(2.);
  // eta''' to a_0 pi
  _incoming.push_back(9020221);_outgoing1.push_back(9000111);_outgoing2.push_back(111);
  _incoming.push_back(9020221);_outgoing1.push_back(9000211);_outgoing2.push_back(-211);
  _coupling.push_back(1.257*GeV);_maxweight.push_back(2.);
  _coupling.push_back(1.257*GeV);_maxweight.push_back(2.);
  // eta''' to sigma eta
  _incoming.push_back(9020221);_outgoing1.push_back(221);_outgoing2.push_back(9000221);
  _coupling.push_back(3.808*GeV);_maxweight.push_back(2.);
  // eta'' to sigma eta
  _incoming.push_back(100221);_outgoing1.push_back(9000221);_outgoing2.push_back(221);
  _coupling.push_back(4.828*GeV);_maxweight.push_back(2.);
  // chi_0c decays to K K 
  _incoming.push_back(10441);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _incoming.push_back(10441);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _coupling.push_back(0.104*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(0.104*GeV);_maxweight.push_back(1.05);
  // chi_0c decays to pi pi
  _incoming.push_back(10441);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _incoming.push_back(10441);_outgoing1.push_back(111);_outgoing2.push_back(111);
  _coupling.push_back(0.093*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(0.066*GeV);_maxweight.push_back(1.05);
  // chi_0c decays to eta eta
  _incoming.push_back(10441);_outgoing1.push_back(221);_outgoing2.push_back(221);
  _coupling.push_back(0.064*GeV);_maxweight.push_back(1.05);
  // f_0(1500) to pipi
  _incoming.push_back(9030221);_outgoing1.push_back(211);_outgoing2.push_back(-211);
  _incoming.push_back(9030221);_outgoing1.push_back(111);_outgoing2.push_back(111);
  _coupling.push_back(1.398*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(0.989*GeV);_maxweight.push_back(1.05);
  // f_0(1500) to sigma sigma
  _incoming.push_back(9030221);_outgoing1.push_back(9000221);
  _outgoing2.push_back(9000221);
  _coupling.push_back(5.902*GeV);_maxweight.push_back(7.);
  // f_0(1500) to eta eta
  _incoming.push_back(9030221);_outgoing1.push_back(221);_outgoing2.push_back(221);
  _coupling.push_back(0.809*GeV);_maxweight.push_back(1.05);
  // f_0(1500) to eta eta'
  _incoming.push_back(9030221);_outgoing1.push_back(221);_outgoing2.push_back(331);
  _coupling.push_back(2.712*GeV);_maxweight.push_back(4.5);
  // f_0(1500) to K K
  _incoming.push_back(9030221);_outgoing1.push_back(321);_outgoing2.push_back(-321);
  _incoming.push_back(9030221);_outgoing1.push_back(311);_outgoing2.push_back(-311);
  _coupling.push_back(0.686*GeV);_maxweight.push_back(1.05);
  _coupling.push_back(0.686*GeV);_maxweight.push_back(1.05);
  // f_0(1500) to pi' pi
  _incoming.push_back(9030221);_outgoing1.push_back(100111);_outgoing2.push_back(111);
  _incoming.push_back(9030221);_outgoing1.push_back(100211);_outgoing2.push_back(-211);
  _coupling.push_back(5.027*GeV);_maxweight.push_back(2.1);
  _coupling.push_back(5.027*GeV);_maxweight.push_back(2.1);
  // initial size
  _initsize = _coupling.size();
  // intermediates
  generateIntermediates(false);
}

int ScalarScalarScalarDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  int imode(-1);
  // must be two outgoing particles
  if(dm.products().size()!=2){return imode;}
  // ids of the particles
  int id0(dm.parent()->id()),id0bar(id0);
  if(dm.parent()->CC()){id0bar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  // loop over the modes and see if this is one of them
  unsigned int ix(0);
  do {
    if(id0   ==_incoming[ix])
      {if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	  (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])){imode=ix;cc=false;}}
    if(id0bar==_incoming[ix]&&imode<0)
      {if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	  (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])){imode=ix;cc=true;}}
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void ScalarScalarScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,MeV)
     << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void ScalarScalarScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,MeV)
     >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<ScalarScalarScalarDecayer> 
ScalarScalarScalarDecayer::initScalarScalarScalarDecayer;
// Definition of the static class description member.

void ScalarScalarScalarDecayer::Init() {

  static ClassDocumentation<ScalarScalarScalarDecayer> documentation
    ("The ScalarScalarScalarDecayer class is designed for the"
     " decay of a scalar meson to two scalar mesons including off-shell effects");

  static ParVector<ScalarScalarScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarScalarScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &ScalarScalarScalarDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &ScalarScalarScalarDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,Energy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarScalarScalarDecayer::_coupling,
     MeV, 0, 0*MeV, 0.0*MeV, 1000000.*MeV, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarScalarScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

double ScalarScalarScalarDecayer::me2(bool vertex, const int,
				   const Particle & inpart,
				   const ParticleVector & decay) const {
  // workaround for gcc 3.2.3 bug
  //ALB ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // set up the spin info for the outgoing particles
  //ALB for(unsigned int ix=0;ix<2;++ix){ScalarWaveFunction(decay[ix],outgoing,true,vertex);}
  for (unsigned int ix=0;ix<2;++ix) {
    PPtr mytemp = decay[ix]; 
    ScalarWaveFunction(mytemp,outgoing,true,vertex);
  }

  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0);
  double fact(_coupling[imode()]/inpart.mass());
  newME(0,0,0) = fact;
  ME(newME);
  return sqr(fact);
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarScalarScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
  unsigned int ix(0); bool order(true);
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
  coupling=_coupling[imode]/dm.parent()->mass();
  itype = 6;
  return order;
}

// output the setup information for the particle database
void ScalarScalarScalarDecayer::dataBaseOutput(ofstream & output,
					       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "set " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "set " << fullName() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "set " << fullName() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "set " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix]/MeV << "\n";
      output << "set " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << fullName() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "insert " << fullName() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "insert " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix]/MeV << "\n";
      output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
}
