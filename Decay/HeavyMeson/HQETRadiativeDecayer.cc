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

HQETRadiativeDecayer::HQETRadiativeDecayer(): Ch_(-1.058455), deltaEta_(1./43.7), Lambda_(1.*GeV),
    incoming_ ({413,423,433,                         //D* EM decay modes: VtoSV
                513,523,533}),                       //B* EM decay modes: VtoSV
    outgoingH_({411,421,431,
                511,521,531}),
    outgoingL_({22,22,22,
                22,22,22}),
    type_     ({1,  1,  1,
                2,  2,  2}),
    maxWeight_({1., 1., 1.,
                1., 1., 1.})
{}

void HQETRadiativeDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoingH_.size()||isize!=outgoingL_.size()||
     isize!=maxWeight_.size()||isize!=type_     .size())
    throw InitException() << "Inconsistent parameters in HQETRadiativeDecayer"
    			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoingH_[ix]),
		     getParticleData(outgoingL_[ix])};
    if(in&&out[0]&&out[1]) {
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
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
  os << Ch_ << deltaEta_<< ounit(Lambda_,GeV) << maxWeight_;
}

void HQETRadiativeDecayer::persistentInput(PersistentIStream & is, int) {
  is >> Ch_ >> deltaEta_>> iunit(Lambda_,GeV) >> maxWeight_;
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

  static Parameter<HQETRadiativeDecayer,double> interfaceDeltaEta
    ("DeltaEta",
     "The mixing parameter for eta-pi0 of isospin violating decays",
     &HQETRadiativeDecayer::deltaEta_, 1./43.7, 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<HQETRadiativeDecayer,Energy> interfacefLambda
    ("Lambda",
     "Strong decays momentum scale",
     &HQETRadiativeDecayer::Lambda_, GeV, 1.*GeV, .1*GeV, 2.*GeV,
     false, false, Interface::limited);
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
      if((id1   ==outgoingH_[ix]&&id2   ==outgoingL_[ix])||
    	 (id2   ==outgoingH_[ix]&&id1   ==outgoingL_[ix])) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoingH_[ix]&&id2bar==outgoingL_[ix])||
    	 (id2bar==outgoingH_[ix]&&id1bar==outgoingL_[ix])) {
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
  			      const tPDVector & outgoing,
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
  Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
					     momenta[1].mass());
  if(abs(type_[imode()])==1 || abs(type_[imode()])==2) {
    // get the polarization vectors
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix==1) continue;
      vecOut_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
    }
    // calculate coupling A
    InvEnergy A;
    Energy mQ;
    double eq(1./3.), eQ(2./3.);
    if(abs(part.id())==423 || abs(part.id())==523) {
      eq = -2./3.;
    }
    if(abs(type_[imode()])==7) { //i.e. for the c quarks
      mQ = getParticleData(ParticleID::c)->mass();
    }
    else { //i.e. for the b quarks
      eQ = -1./3.;
      mQ = getParticleData(ParticleID::b)->mass();
    }
    if(part.id()<0) {
      eq *= -1.;
      eQ *= -1.;
    }
    A = eQ/4./mQ*sqrt(SM().alphaEM(sqr(mQ)))+Ch_/Lambda_*eq*sqrt(SM().alphaEM(sqr(Lambda_)));
    //calculate ME
    InvEnergy2 fact = 4.*sqrt(8.*M_PI)*A*sqrt(momenta[0].mass()/part.mass())/part.mass();
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int iy=0;iy<3;++iy) {
        if(iy==1) {
          (*ME())(ix,iy,0) = 0.;
        }
        else {
	         (*ME())(ix,iy,0) = Complex(fact*epsilon(vecOut_[iy],momenta[1],vecIn_[ix])
                      *part.momentum());
        }
      }
    }
    // analytic test of the answer
    //test = 32.*M_PI*sqr(A)*momenta[0].mass()/part.mass()/3.
    //     *sqr(sqr(part.mass())-sqr(momenta[0].mass()))/sqr(part.mass());
  }
  else {
    assert(false);
  }
  double output = ME()->contract(rho_).real();
  // testing
  // double ratio = (output-test)/(output+test);
  // generator()->log() << "testing matrix element for " << part.PDGName() << " -> "
  // 		     << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  // 		     << output << " " << test << " " << ratio << endl;
  // isospin factors
  if(abs(outgoing[1]->id())==ParticleID::pi0) {
    output *= type_[imode()]>0 ? 0.5 : 0.125*sqr(deltaEta_);
  }
  // return the answer
  return output;
}

bool HQETRadiativeDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
				      double & coupling) const {
  // int imode(-1);
  // int id(dm.parent()->id());
  // int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  // ParticleMSet::const_iterator pit(dm.products().begin());
  // int id1((**pit).id());
  // int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  // ++pit;
  // int id2((**pit).id());
  // int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  // unsigned int ix(0); bool order(false);
  // do {
  //   if(id   ==incoming_[ix]) {
  //     if(id1==outgoingH_[ix]&&id2==outgoingL_[ix]) {
  // 	imode=ix;
  // 	order=true;
  //     }
  //     if(id2==outgoingH_[ix]&&id1==outgoingL_[ix]) {
  // 	imode=ix;
  // 	order=false;
  //     }
  //   }
  //   if(idbar==incoming_[ix]&&imode<0) {
  //     if(id1bar==outgoingH_[ix]&&id2bar==outgoingL_[ix]) {
  // 	imode=ix;
  // 	order=true;
  //     }
  //     if(id2bar==outgoingH_[ix]&&id1bar==outgoingL_[ix]) {
  // 	imode=ix;
  // 	order=false;
  //     }
  //   }
  //   ++ix;
  // }
  // while(ix<incoming_.size()&&imode<0);
  // coupling=_coupling[imode]*dm.parent()->mass();
  // mecode=7;
  // return order;
}

void HQETRadiativeDecayer::dataBaseOutput(ofstream & output,
  				       bool header) const {
                   cerr<<"gets is dataBaseOutput\n";
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  // couplings
  output << "newdef " << name() << ":Ch     " << Ch_         << "\n";
  output << "newdef " << name() << ":Lambda " << Lambda_/GeV << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\""
  		    << fullName() << "\";" << endl;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "newdef " << name() << ":MaxWeight " << ix << " " << maxWeight_[ix] << "\n";
  }
}
