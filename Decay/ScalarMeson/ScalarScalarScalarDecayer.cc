// -*- C++ -*-
//
// ScalarScalarScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarScalarScalarDecayer class.
//

#include "ScalarScalarScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarScalarScalarDecayer::doinitrun() {
  DecayIntegrator2::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix) ) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

ScalarScalarScalarDecayer::ScalarScalarScalarDecayer() {
  _incoming  = {9010221, 9010221, 9010221, 9010221,   10221,   10221,
	          10221,   10221,   10221,   10221,   10331,   10331,
	          10331,   10331, 9000221, 9000221, 9000111, 9000211,
	        9000111, 9000111, 9000211,   10111,   10211,   10111,
	          10211,   10111,   10111,   10211,   10221,   10331,
	          10221,   10311,   10311,   10321,   10321,   10411,
	          10411,   10421,   10421,   10511,   10511,   10521,
	          10521,  100311,  100311,  100321,  100321,   10431,
	          10531,  100221,  100221, 9020221, 9020221, 9020221,
	         100221,   10441,   10441,   10441,   10441,   10441,
	        9030221, 9030221, 9030221, 9030221, 9030221, 9030221,
	        9030221, 9030221, 9030221, 9000311, 9000311, 9000321,
	        9000321,   10441,   10441,   10531,   10531,   10441};
  _outgoing1 = {    111,     211,     321,     311,     111,     211,
		 100111,  100211,     321,     311,     111,     211,
		    321,     311,     111,     211,     221,     221,
		    321,     311,     321,     221,     221,     331,
		    331,     321,     311,     321,     221,     221,
		9000221,     311,     321,     321,     311,     411,
		    421,     421,     411,     511,     521,     521,
		    511,   10311,   10321,   10321,   10311,     431,
		    531, 9000111, 9000211, 9000111, 9000211,     221,
		9000221,     321,     311,     211,     111,     221,
		    211,     111, 9000221,     221,     221,     321,
		    311,  100111,  100211,     311,     321,     321,
		    311,   10321,   10311,     511,     521, 9010221};
  _outgoing2 = {    111,    -211,    -321,    -311,     111,    -211,
		    111,    -211,    -321,    -311,     111,    -211,
		   -321,    -311,     111,    -211,     111,     211,
		   -321,    -311,    -311,     111,     211,     111,
		    211,    -321,    -311,    -311,     221,     221,
		9000221,     111,    -211,     111,     211,     111,
		    211,     111,    -211,     111,    -211,     111,
		    211,     111,    -211,     111,     211,     111,
		    111,     111,    -211,     111,    -211, 9000221,
		    221,    -321,    -311,    -211,     111,     221,
		   -211,     111, 9000221,     221,     331,    -321,
		   -311,     111,    -211,     111,    -211,     111,
		    211,  -10321,  -10311,    -311,    -321, 9010221};
  _coupling  ={ 1.66*GeV,  2.35*GeV,  1.02*GeV,  1.02*GeV, 0.745*GeV, 1.054*GeV,
	       5.027*GeV, 5.027*GeV, 0.886*GeV, 0.886*GeV, 0.503*GeV, 0.711*GeV,
	       2.096*GeV, 2.096*GeV, 3.654*GeV, 5.178*GeV,  3.33*GeV,  3.33*GeV,
		2.54*GeV,  2.54*GeV,  3.59*GeV, 1.357*GeV, 1.357*GeV, 0.995*GeV,
	       0.995*GeV, 0.950*GeV, 0.950*GeV, 1.344*GeV, 0.235*GeV, 2.189*GeV,
	       21.46*GeV, 2.837*GeV, 4.000*GeV, 2.837*GeV, 4.000*GeV, 5.472*GeV,
	       7.714*GeV, 5.447*GeV, 7.818*GeV, 9.698*GeV, 13.71*GeV, 9.698*GeV,
	       13.71*GeV, 6.595*GeV, 9.445*GeV, 6.595*GeV, 9.445*GeV, 0.103*GeV,
	       8.314*GeV, 2.057*GeV, 2.057*GeV, 1.470*GeV, 1.470*GeV, 4.051*GeV,
	       4.316*GeV, 0.104*GeV, 0.104*GeV, 0.093*GeV, 0.066*GeV, 0.064*GeV,
	       1.398*GeV, 0.989*GeV, 6.079*GeV, 0.809*GeV, 2.844*GeV, 0.686*GeV,
	       0.686*GeV, 2.615*GeV, 2.615*GeV, 3.834*GeV, 5.406*GeV, 3.834*GeV,
	       5.406*GeV, 0.104*GeV, 0.104*GeV, 12.17*GeV, 12.17*GeV, 0.084*GeV};
  _maxweight={1.05, 1.05, 1.05, 1.05, 1.05, 1.05,
	       2.1,  2.1,  1.1,  1.1, 1.05, 1.05,
	      1.05, 1.05, 1.05, 1.05,  1.1,  1.1,
	      1.05, 1.05, 1.05, 1.05, 1.05,  1.6,
	       1.6, 1.05, 1.05, 1.05,  1.1,  1.2,
	       7.0, 1.05, 1.05, 1.05, 1.05, 1.05,
	      1.05, 1.05, 1.05, 1.05, 1.05, 1.05,
	      1.05,   2.,   2.,   2.,   2., 1.05,
	      1.05,   2.,   2.,   2.,   2.,   2.,
	        2., 1.05, 1.05, 1.05, 1.05, 1.05,
	      1.05, 1.05,  7.0, 1.05, 4.5,  1.05,
	      1.05,  2.1,  2.1, 1.05, 1.05, 1.05,
	      1.05, 1.05, 1.05, 1.05, 1.05, 1.05};
  // initial size
  _initsize = _coupling.size();
  // intermediates
  generateIntermediates(false);
}

void ScalarScalarScalarDecayer::doinit() {
  DecayIntegrator2::doinit();
  // check the parameters arew consistent
  unsigned int isize(_coupling.size());
  if(isize!=_incoming.size()  || isize!=_outgoing1.size()||
     isize!=_outgoing2.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in ScalarScalarScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    tPDPtr     in = getParticleData(_incoming[ix]);
    tPDVector out = {getParticleData(_outgoing1[ix]),
		     getParticleData(_outgoing2[ix])};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,_maxweight[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int ScalarScalarScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id0(parent->id());
  int id0bar = parent->CC() ? parent->CC()->id() : id0;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  // loop over the modes and see if this is one of them
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id0   ==_incoming[ix]) {
      if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	 (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==_incoming[ix]&&imode<0) {
      if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	 (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void ScalarScalarScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,MeV) << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void ScalarScalarScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,MeV) >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarScalarScalarDecayer,DecayIntegrator2>
describeHerwigScalarScalarScalarDecayer("Herwig::ScalarScalarScalarDecayer", "HwSMDecay.so");

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
     MeV, 0, ZERO, ZERO, 1000000.*MeV, false, false, true);

  static ParVector<ScalarScalarScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarScalarScalarDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

void ScalarScalarScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double ScalarScalarScalarDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> &,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
  }
  double fact(_coupling[imode()]/part.mass());
  (*ME())(0,0,0) = fact;
  return sqr(fact);
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarScalarScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
					       double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0);
  bool order(true);
  do {
    if(id   ==_incoming[ix]) {
      if(id1==_outgoing1[ix]&&id2==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoing1[ix]&&id1==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
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
  // parameters for the DecayIntegrator2 base class
  DecayIntegrator2::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "newdef " << name() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]/MeV << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "insert " << name() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]/MeV << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
