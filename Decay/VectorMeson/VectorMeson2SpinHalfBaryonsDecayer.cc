// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2SpinHalfBaryonsDecayer class.
//

#include "VectorMeson2SpinHalfBaryonsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson2SpinHalfBaryonsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMeson2SpinHalfBaryonsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=gm_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!= maxweight_.size()||
     isize!= phi_.size() || isize!= ge_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "BaryonsDecayer::doiin() " << Exception::runerror;
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

int VectorMeson2SpinHalfBaryonsDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

IBPtr VectorMeson2SpinHalfBaryonsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VectorMeson2SpinHalfBaryonsDecayer::fullclone() const {
  return new_ptr(*this);
}

void VectorMeson2SpinHalfBaryonsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ge_ << gm_ << phi_ << incoming_ << outgoing_ << maxweight_;
}

void VectorMeson2SpinHalfBaryonsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> ge_ >> gm_ >> phi_ >> incoming_ >> outgoing_ >> maxweight_;
}


DescribeClass<VectorMeson2SpinHalfBaryonsDecayer,DecayIntegrator>
describeHerwigVectorMeson2SpinHalfBaryonsDecayer("Herwig::VectorMeson2SpinHalfBaryonsDecayer", "HwVMDecay.so");

void VectorMeson2SpinHalfBaryonsDecayer::Init() {

  static ClassDocumentation<VectorMeson2SpinHalfBaryonsDecayer> documentation
    ("The VectorMeson2SpinHalfBaryonsDecayer class is designed for "
     "the decay of vector mesons to baryons.");
  
  static Command<VectorMeson2SpinHalfBaryonsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, baryon, antibaryon), GE, GM, phase and max weight for a decay",
     &VectorMeson2SpinHalfBaryonsDecayer::setUpDecayMode, false);
}

void VectorMeson2SpinHalfBaryonsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // outgoing fermion
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  // outgoing antifermion
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double VectorMeson2SpinHalfBaryonsDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  // initialze me
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // coefficients
  Complex GM = gm_[imode()];
  Complex GE = ge_[imode()]*exp(Complex(0.,phi_[imode()]));
  LorentzPolarizationVector c2 = -2.*outgoing[0]->mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()))*
    (GM-GE)*(momenta[iferm]-momenta[ianti]);
  // now compute the currents
  LorentzPolarizationVector temp;
  //double mesum(0.);
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      LorentzPolarizationVector temp = (GM*wave_[ix].vectorCurrent(wavebar_[iy])+c2*wave_[ix].scalar(wavebar_[iy]))/part.mass();
      for(unsigned int iz=0;iz<3;++iz) {
	if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
	//mesum += norm(vectors_[iz].dot(temp));
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // cerr << "testing decay " << part.PDGName() << " -> " << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << "\n";
  // cerr << "testing ME " << mesum/3. << " " << me << " " << 4./3.*(norm(GM)+2.*sqr(outgoing[0]->mass()/part.mass())*norm(GE)) << "\n";
  // cerr << "testing gamma " << mesum/3./8./Constants::pi*sqrt(sqr(part.mass())-4.*sqr(outgoing[0]->mass()))/MeV << "\n";
  // return the answer
  return me;
}

// output the setup information for the particle database
void VectorMeson2SpinHalfBaryonsDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << gm_[ix] << " " << ge_[ix] << " " << phi_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string VectorMeson2SpinHalfBaryonsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
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
  // get the couplings
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  ge_.push_back(stof(stype));
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  gm_.push_back(stof(stype));
  // and phase
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  phi_.push_back(stof(stype));
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
