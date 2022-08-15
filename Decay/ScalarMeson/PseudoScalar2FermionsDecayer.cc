// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PseudoScalar2FermionsDecayer class.
//

#include "PseudoScalar2FermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PseudoScalar2FermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}


void PseudoScalar2FermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "FermionDecayer::doiin() " << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

PseudoScalar2FermionsDecayer::PseudoScalar2FermionsDecayer() {
  // don't include intermediates
  generateIntermediates(false);
}

int PseudoScalar2FermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
					     const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(incoming_[ix]==id   ) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(incoming_[ix]==idbar) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  return imode;
}

IBPtr PseudoScalar2FermionsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr PseudoScalar2FermionsDecayer::fullclone() const {
  return new_ptr(*this);
}

void PseudoScalar2FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << coupling_ << incoming_ << outgoing_ << maxweight_;
}

void PseudoScalar2FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> coupling_ >> incoming_ >> outgoing_ >> maxweight_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PseudoScalar2FermionsDecayer,DecayIntegrator>
describeHerwigPseudoScalar2FermionsDecayer("Herwig::PseudoScalar2FermionsDecayer", "HwSMDecay.so");

void PseudoScalar2FermionsDecayer::Init() {

  static ClassDocumentation<PseudoScalar2FermionsDecayer> documentation
    ("The PseudoScalar2FermionsDecayer class implements the decay of a pseudoscalar meson "
     "to a fermion and antifermion.");
  
  static Command<PseudoScalar2FermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, fermion, antifermion), coupling and max weight for a decay",
     &PseudoScalar2FermionsDecayer::setUpDecayMode, false);

}

void PseudoScalar2FermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  // set up the spin information for the decay products
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}


double PseudoScalar2FermionsDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // prefactor
  InvEnergy pre(coupling_[imode()]/part.mass());
  // now compute the ME
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      Complex temp = pre*wave_[ix].pseudoScalar(wavebar_[iy]);
      if(iferm>ianti) (*ME())(0,ix,iy)=temp;
      else            (*ME())(0,iy,ix)=temp;
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // double test = 2*sqr(coupling_[imode()])*(1.-sqr((momenta[0].mass()-momenta[1].mass())/part.mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool PseudoScalar2FermionsDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						double & coupling) const {
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
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode];
  mecode=2;
  return order;
}

// output the setup information for the particle database
void PseudoScalar2FermionsDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string PseudoScalar2FermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxweight_.push_back(wgt);
  // success
  return "";
}
