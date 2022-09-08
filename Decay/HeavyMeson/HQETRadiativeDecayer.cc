// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HQETRadiativeDecayer class.
//

#include "HQETRadiativeDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

using namespace Herwig;

void HQETRadiativeDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=type_     .size())
    throw InitException() << "Inconsistent parameters in HQETRadiativeDecayer"
    			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  tPDPtr gamma = getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix]),gamma};
    if(in&&out[0]&&out[1]) {
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
    // calculate coupling A
    Energy mQ = abs(incoming_[ix])/100 == 4 ?
      getParticleData(ParticleID::c)->constituentMass() :
      getParticleData(ParticleID::b)->constituentMass() ;
    double eq = ((abs(incoming_[ix])/10)%10)==2 ?  2./3. : -1./3.;
    double eQ =  abs(incoming_[ix])/100==4      ?  2./3. : -1./3.;
    coupling_.push_back(8.*sqrt(Constants::pi)*
			(eQ/4./mQ*sqrt(SM().alphaEM(sqr(mQ)))+Ch_/Lambda_*eq*sqrt(SM().alphaEM(sqr(Lambda_)))));
  }
}

void HQETRadiativeDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

IBPtr HQETRadiativeDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr HQETRadiativeDecayer::fullclone() const {
  return new_ptr(*this);
}

void HQETRadiativeDecayer::persistentOutput(PersistentOStream & os) const {
  os << Ch_ << ounit(Lambda_,GeV) << incoming_ << outgoing_ << type_ << maxWeight_ << ounit(coupling_,1./GeV);
}

void HQETRadiativeDecayer::persistentInput(PersistentIStream & is, int) {
  is >> Ch_ >> iunit(Lambda_,GeV) >> incoming_ >> outgoing_ >> type_ >> maxWeight_ >> iunit(coupling_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HQETRadiativeDecayer,DecayIntegrator>
describeHerwigHQETRadiativeDecayer("Herwig::HQETRadiativeDecayer", "HwHMDecay.so");

void HQETRadiativeDecayer::Init() {

  static ClassDocumentation<HQETRadiativeDecayer> documentation
    ("The HQETRadiativeDecayer class performs the EM decays of excited heavy mesons using HQET results.");

  static Parameter<HQETRadiativeDecayer,double> interfaceCh
     ("Ch",
      "EM coefficient for heavy meson decays",
      &HQETRadiativeDecayer::Ch_, -1.058455, -2.0, 2.0,
      false, false, Interface::limited);

  static ParVector<HQETRadiativeDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &HQETRadiativeDecayer::maxWeight_,
     0, 0, 0, 0., 100000., false, false, true);

  static Parameter<HQETRadiativeDecayer,Energy> interfacefLambda
    ("Lambda",
     "Strong decays momentum scale",
     &HQETRadiativeDecayer::Lambda_, GeV, 1.*GeV, .1*GeV, 2.*GeV, false, false, Interface::limited);

  static Command<HQETRadiativeDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(1/GeV^2) and max weight for a decay",
     &HQETRadiativeDecayer::setUpDecayMode, false);
}

int HQETRadiativeDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix]&&id2 == ParticleID::gamma)||
    	 (id2   ==outgoing_[ix]&&id1 == ParticleID::gamma)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix]&&id2bar==ParticleID::gamma)||
    	 (id2bar==outgoing_[ix]&&id1bar==ParticleID::gamma)) {
    	imode=ix;
    	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void HQETRadiativeDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  switch (part.dataPtr()->iSpin()) {
  case PDT::Spin0:
    Helicity::ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
						    Helicity::incoming,true);
    break;
  case PDT::Spin1:
    Helicity::VectorWaveFunction::constructSpinInfo(vecIn_,const_ptr_cast<tPPtr>(&part),
  						    Helicity::incoming,true,false);
    break;
  case PDT::Spin2:
    Helicity::TensorWaveFunction::constructSpinInfo(tensorIn_,const_ptr_cast<tPPtr>(&part),
						    Helicity::incoming,true,false);
    break;
  default:
    assert(false);
  }
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix) {
    switch (decay[ix]->dataPtr()->iSpin()) {
    case PDT::Spin1:
      Helicity::VectorWaveFunction::constructSpinInfo(vecOut_,decay[ix],
						      Helicity::outgoing,true,false);
      break;
    case PDT::Spin0:
      Helicity::ScalarWaveFunction::constructSpinInfo(decay[ix],Helicity::outgoing,true);
      break;
    default:
      assert(false);
    }
  }
}

// matrix elememt for the process
double HQETRadiativeDecayer::me2(const int, const Particle & part,
				 const tPDVector & ,
				 const vector<Lorentz5Momentum> & momenta,
				 MEOption meopt) const {
  if(!ME()) {
    if(abs(type_[imode()])==1 || abs(type_[imode()])==2) {
      ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin1)));
    }
  }
  // stuff for incoming particle
  if(meopt==Initialize) {
    unsigned int Type = abs(type_[imode()]);
    if(Type==1 || Type==2) {
      rho_ = RhoDMatrix(PDT::Spin1);
      Helicity::VectorWaveFunction::calculateWaveFunctions(vecIn_,rho_,const_ptr_cast<tPPtr>(&part),
  							   Helicity::incoming,false);
    }
    else {
      cerr << "Unknown decay mode: type " << type_[imode()] << " for " << part << "\n";
      assert(false);
    }
  }
  // calculate the matrix element
  // double test(0.);
  if(abs(type_[imode()])==1 || abs(type_[imode()])==2) {
    // get the polarization vectors
    vecOut_.resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix==1) continue;
      vecOut_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
    }
    // calculate ME
    InvEnergy2 fact = coupling_[imode()]*sqrt(momenta[0].mass()/part.mass())/part.mass();
    for(unsigned int ix=0;ix<3;++ix) {
      LorentzVector<complex<Energy2> > v0 = epsilon(part.momentum(),momenta[1],vecIn_[ix]);
      for(unsigned int iy=0;iy<3;++iy) {
        if(iy==1) {
          (*ME())(ix,0,iy) = 0.;
        }
        else {
	  (*ME())(ix,0,iy) = Complex(fact*vecOut_[iy].dot(v0));
        }
      }
    }
    // analytic test of the answer
    // test = 0.5*sqr(coupling_[imode()])*momenta[0].mass()/part.mass()/3.
    //   *sqr(sqr(part.mass())-sqr(momenta[0].mass()))/sqr(part.mass());
  }
  else {
    assert(false);
  }
  double output = ME()->contract(rho_).real();
  // testing
  // double ratio = (output-test)/(output+test);
  // cout << "testing matrix element for " << part.PDGName() << " -> "
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << output << " " << test << " " << ratio << " " << output/test << endl;
  // return the answer
  return output;
}

bool HQETRadiativeDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix]&&id2==ParticleID::gamma) {
  	imode=ix;
  	order=true;
      }
      if(id2==outgoing_[ix]&&id1==ParticleID::gamma) {
  	imode=ix;
  	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix]&&id2bar==ParticleID::gamma) {
  	imode=ix;
  	order=true;
      }
      if(id2bar==outgoing_[ix]&&id1bar==ParticleID::gamma) {
  	imode=ix;
  	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=abs(coupling_[imode]*dm.parent()->mass());
  mecode=18;
  return order;
}

void HQETRadiativeDecayer::dataBaseOutput(ofstream & output,
  				       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  // couplings
  output << "newdef " << name() << ":Ch     " << Ch_         << "\n";
  output << "newdef " << name() << ":Lambda " << Lambda_/GeV << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix] << " "
	   << type_[ix] << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\""
  		    << fullName() << "\";" << endl;
}

string HQETRadiativeDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(int(pData->iSpin())%2!=1)
    return "Incoming particle with id " + std::to_string(in) + "does not integer spin";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0 && pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out) + "does not have spin 0/1";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int itype = stoi(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  type_.push_back(itype);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
