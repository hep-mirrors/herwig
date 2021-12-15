// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2SpinThreeHalfBaryonsDecayer class.
//

#include "VectorMeson2SpinThreeHalfBaryonsDecayer.h"
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

void VectorMeson2SpinThreeHalfBaryonsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::doinit() {
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

int VectorMeson2SpinThreeHalfBaryonsDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

IBPtr VectorMeson2SpinThreeHalfBaryonsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VectorMeson2SpinThreeHalfBaryonsDecayer::fullclone() const {
  return new_ptr(*this);
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ge_ << gm_ << phi_ << incoming_ << outgoing_ << maxweight_;
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> ge_ >> gm_ >> phi_ >> incoming_ >> outgoing_ >> maxweight_;
}


DescribeClass<VectorMeson2SpinThreeHalfBaryonsDecayer,DecayIntegrator>
describeHerwigVectorMeson2SpinThreeHalfBaryonsDecayer("Herwig::VectorMeson2SpinThreeHalfBaryonsDecayer", "HwVMDecay.so");

void VectorMeson2SpinThreeHalfBaryonsDecayer::Init() {

  static ClassDocumentation<VectorMeson2SpinThreeHalfBaryonsDecayer> documentation
    ("The VectorMeson2SpinThreeHalfBaryonsDecayer class is designed for "
     "the decay of vector mesons to baryons.");
  
  static Command<VectorMeson2SpinThreeHalfBaryonsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, baryon, antibaryon), GE, GM, phase and max weight for a decay",
     &VectorMeson2SpinThreeHalfBaryonsDecayer::setUpDecayMode, false);
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // outgoing fermion
  RSSpinorBarWaveFunction::
    constructSpinInfo(wave2bar_,decay[iferm],outgoing,true);
  // outgoing antifermion
  RSSpinorWaveFunction::
    constructSpinInfo(wave2_,decay[ianti],outgoing,true);
}

double VectorMeson2SpinThreeHalfBaryonsDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  // initialze me
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin3Half,PDT::Spin3Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave2_.resize(4);
  wave2bar_.resize(4);
  wave_.resize(4);
  wavebar_.resize(4);
  RSSpinorBarWaveFunction swave(momenta[iferm],outgoing[iferm],Helicity::outgoing);
  RSSpinorWaveFunction    awave(momenta[ianti],outgoing[ianti],Helicity::outgoing);
  LorentzPolarizationVector vtemp = part.momentum()/part.mass();
  for(unsigned int ix=0;ix<4;++ix) {
    swave.reset(ix);
    awave.reset(ix);
    wave2bar_[ix] = swave.dimensionedWf();
    wavebar_ [ix] = wave2bar_[ix].dot(vtemp);
    wave2_   [ix] = awave.dimensionedWf();
    wave_    [ix] = wave2_[ix].dot(vtemp);
  }
  // coefficients
  Complex GM = gm_[imode()];
  Complex GE = ge_[imode()]*exp(Complex(0.,phi_[imode()]));
  LorentzPolarizationVector c2 = -2.*outgoing[0]->mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()))*
    (GM-GE)*(momenta[iferm]-momenta[ianti]);
  // now compute the currents
  for(unsigned ix=0;ix<4;++ix) {
    for(unsigned iy=0;iy<4;++iy) {
      // q(al)q(be) piece
      LorentzPolarizationVector temp2 = (GM*wave_[ix].vectorCurrent(wavebar_[iy])+c2*wave_[ix].scalar(wavebar_[iy]))*
	2.*part.mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()));
      // g(al)g(be) * GM-GE piece
      LorentzPolarizationVector temp3 = wave2_[ix].generalScalar(wave2bar_[iy],1.,1.)*c2/part.mass();
      // g(al)g(be) * gamma^mu
      LorentzPolarizationVector temp1(GM/part.mass()*(wave2bar_[iy](0,3)*wave2_[ix](0,0) + wave2bar_[iy](0,2)*wave2_[ix](0,1) - 
						      wave2bar_[iy](0,1)*wave2_[ix](0,2) - wave2bar_[iy](0,0)*wave2_[ix](0,3) + 
						      wave2bar_[iy](1,3)*wave2_[ix](1,0) + wave2bar_[iy](1,2)*wave2_[ix](1,1) - 
						      wave2bar_[iy](1,1)*wave2_[ix](1,2) - wave2bar_[iy](1,0)*wave2_[ix](1,3) + 
						      wave2bar_[iy](2,3)*wave2_[ix](2,0) + wave2bar_[iy](2,2)*wave2_[ix](2,1) - 
						      wave2bar_[iy](2,1)*wave2_[ix](2,2) - wave2bar_[iy](2,0)*wave2_[ix](2,3) - 
						      wave2bar_[iy](3,3)*wave2_[ix](3,0) - wave2bar_[iy](3,2)*wave2_[ix](3,1) + 
						      wave2bar_[iy](3,1)*wave2_[ix](3,2) + wave2bar_[iy](3,0)*wave2_[ix](3,3)),
				      Complex(0,1)*GM/part.mass()*(wave2bar_[iy](0,3)*wave2_[ix](0,0) - wave2bar_[iy](0,2)*wave2_[ix](0,1) - 
								   wave2bar_[iy](0,1)*wave2_[ix](0,2) + wave2bar_[iy](0,0)*wave2_[ix](0,3) + 
								   wave2bar_[iy](1,3)*wave2_[ix](1,0) - wave2bar_[iy](1,2)*wave2_[ix](1,1) - 
								   wave2bar_[iy](1,1)*wave2_[ix](1,2) + wave2bar_[iy](1,0)*wave2_[ix](1,3) + 
								   wave2bar_[iy](2,3)*wave2_[ix](2,0) - wave2bar_[iy](2,2)*wave2_[ix](2,1) - 
								   wave2bar_[iy](2,1)*wave2_[ix](2,2) + wave2bar_[iy](2,0)*wave2_[ix](2,3) - 
								   wave2bar_[iy](3,3)*wave2_[ix](3,0) + wave2bar_[iy](3,2)*wave2_[ix](3,1) + 
								   wave2bar_[iy](3,1)*wave2_[ix](3,2) - wave2bar_[iy](3,0)*wave2_[ix](3,3)),
				      GM/part.mass()*(wave2bar_[iy](0,2)*wave2_[ix](0,0) - wave2bar_[iy](0,3)*wave2_[ix](0,1) - 
						      wave2bar_[iy](0,0)*wave2_[ix](0,2) + wave2bar_[iy](0,1)*wave2_[ix](0,3) + 
						      wave2bar_[iy](1,2)*wave2_[ix](1,0) - wave2bar_[iy](1,3)*wave2_[ix](1,1) - 
						      wave2bar_[iy](1,0)*wave2_[ix](1,2) + wave2bar_[iy](1,1)*wave2_[ix](1,3) + 
						      wave2bar_[iy](2,2)*wave2_[ix](2,0) - wave2bar_[iy](2,3)*wave2_[ix](2,1) - 
						      wave2bar_[iy](2,0)*wave2_[ix](2,2) + wave2bar_[iy](2,1)*wave2_[ix](2,3) - 
						      wave2bar_[iy](3,2)*wave2_[ix](3,0) + wave2bar_[iy](3,3)*wave2_[ix](3,1) + 
						      wave2bar_[iy](3,0)*wave2_[ix](3,2) - wave2bar_[iy](3,1)*wave2_[ix](3,3)),
				      GM/part.mass()*(-wave2bar_[iy](0,2)*wave2_[ix](0,0) - wave2bar_[iy](0,3)*wave2_[ix](0,1) - 
						      wave2bar_[iy](0,0)*wave2_[ix](0,2) - wave2bar_[iy](0,1)*wave2_[ix](0,3) - 
						      wave2bar_[iy](1,2)*wave2_[ix](1,0) - wave2bar_[iy](1,3)*wave2_[ix](1,1) - 
						      wave2bar_[iy](1,0)*wave2_[ix](1,2) - wave2bar_[iy](1,1)*wave2_[ix](1,3) - 
						      wave2bar_[iy](2,2)*wave2_[ix](2,0) - wave2bar_[iy](2,3)*wave2_[ix](2,1) - 
						      wave2bar_[iy](2,0)*wave2_[ix](2,2) - wave2bar_[iy](2,1)*wave2_[ix](2,3) + 
						      wave2bar_[iy](3,2)*wave2_[ix](3,0) + wave2bar_[iy](3,3)*wave2_[ix](3,1) + 
						      wave2bar_[iy](3,0)*wave2_[ix](3,2) + wave2bar_[iy](3,1)*wave2_[ix](3,3)));
      LorentzPolarizationVector temp = temp1+temp2+temp3;
      for(unsigned int iz=0;iz<3;++iz) {
	if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
      }
    }
  }
  // double mesum = ME()->contract(RhoDMatrix(PDT::Spin1)).real();
  // generator()->log() << "testing decay " << part.PDGName() << " -> " << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << "\n";
  // generator()->log() << "testing ME " << mesum  << " " << me << " " << 1./3.*(40.*norm(GM)/9.+16.*sqr(outgoing[0]->mass()/part.mass())*norm(GE)) << "\n";
  // return the answer
  return ME()->contract(rho_).real();
}

// output the setup information for the particle database
void VectorMeson2SpinThreeHalfBaryonsDecayer::dataBaseOutput(ofstream & output,
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

string VectorMeson2SpinThreeHalfBaryonsDecayer::setUpDecayMode(string arg) {
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
  if(pData->iSpin()!=PDT::Spin3Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3Half)
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
