// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the f1RhoPiPiDecayer class.
//

#include "f1RhoPiPiDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

f1RhoPiPiDecayer::f1RhoPiPiDecayer() :
  ga1RhoPi_(4.8*GeV), gf1a1Pi_(9.77/GeV), maxWeight_({1.,1.}) {
  generateIntermediates(true);
}

// normally not implemented but do it here to get rid of the a_1 which
// can't be on-shell
ParticleVector f1RhoPiPiDecayer::decay(const Particle & parent,
				       const tPDVector & children) const {
  ParticleVector output = DecayIntegrator::decay(parent,children);
  ParticleVector::iterator it =output.begin();
  for(;it!=output.end();++it) {
    long id = (**it).id();
    if(id==20113 || id==20213 || id==-20213)
      break;
  }
  if(it!=output.end()) {
    PPtr a1 = *it;
    output.erase(it);
    for(PPtr child : a1->children()) {
      output.push_back(child);
      a1->abandonChild(child);
    }
  }
  return output;
}


IBPtr f1RhoPiPiDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr f1RhoPiPiDecayer::fullclone() const {
  return new_ptr(*this);
}

void f1RhoPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(ga1RhoPi_,GeV) << ounit(gf1a1Pi_,1./GeV)
     << ounit(ma1_,GeV) << ounit(ga1_,GeV) << maxWeight_;
}

void f1RhoPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(ga1RhoPi_,GeV) >> iunit(gf1a1Pi_,1./GeV)
     >> iunit(ma1_,GeV) >> iunit(ga1_,GeV) >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<f1RhoPiPiDecayer,DecayIntegrator>
describeHerwigf1RhoPiPiDecayer("Herwig::f1RhoPiPiDecayer", "HwVMDecay.so");

void f1RhoPiPiDecayer::Init() {

  static ClassDocumentation<f1RhoPiPiDecayer> documentation
    ("The f1RhoPiPiDecayer class implements a simple model for "
     "f1 -> pipirho via an intermediate a_1");

  static ParVector<f1RhoPiPiDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "Maximum weights for the decays",
     &f1RhoPiPiDecayer::maxWeight_, 2, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<f1RhoPiPiDecayer,Energy> interfacega1RhoPi
    ("ga1RhoPi",
     "Coupling of the a_1 to rho pi",
     &f1RhoPiPiDecayer::ga1RhoPi_, GeV, 4.8*GeV, 0.*GeV, 20.*GeV,
     false, false, Interface::limited);
  
  static Parameter<f1RhoPiPiDecayer,InvEnergy> interfacegf1a1Pi
    ("gf1a1Pi",
     "Coupling of f_1 to a_1 pi",
     &f1RhoPiPiDecayer::gf1a1Pi_, 1./GeV, 1.0/GeV, 0.0/GeV, 10.0/GeV,
     false, false, Interface::limited);

}

void f1RhoPiPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // pointers to the particles we need as external particles
  tPDPtr f1 = getParticleData(ParticleID::f_1);
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a1m = getParticleData(ParticleID::a_1minus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  tPDPtr rhop = getParticleData(ParticleID::rhoplus);
  tPDPtr rhom = getParticleData(ParticleID::rhominus);
  tPDPtr rho0 = getParticleData(ParticleID::rho0);
  // decay mode f_1 -> pi+ pi- rho0
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(f1,{pip,pim,rho0},maxWeight_[0]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,3));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,1,1,3));
  addMode(mode);
  // decay mode f_1 -> pi- pi0 rho+
  mode = new_ptr(PhaseSpaceMode(f1,{pim,pi0,rhop},maxWeight_[1]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,1,1,2,1,3));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,1,1,3));
  addMode(mode);
  ma1_ = a1p->mass();
  ga1_ = a1p->width(); 
}

void f1RhoPiPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<maxWeight_.size();++ix)
      maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

int f1RhoPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  if(children.size()!=3 || parent->id()!=ParticleID::f_1) return -1;
  // check the pions
  int npi0(0),npiplus(0),npiminus(0),idrho(0);
  for(auto child : children) {
    int idtemp= child->id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
    else                                 idrho=idtemp;
  }
  cc = false;
  // f_1 -> pi+pi-rho0 mode
  if(idrho==113 && npiplus==1 && npiminus==1)         return 0;
  // f_1 -> pi-pi0rho+ mode
  else if  (idrho== 213 && npiminus==1 && npi0==1)    return 1;
  else if  (idrho==-213 && npiplus ==1 && npi0==1) {
    cc=true;
    return 1;
  }
  // not found
  return -1;
}

void f1RhoPiPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
  					incoming,true,false);
  // set up the spin information for the decay products
  // pions
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  // rho
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[2],outgoing,true,false);
}

double f1RhoPiPiDecayer::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  useMe();
  // polarization vectors for incoming
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  // polarization vectors for rho
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix)
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  // breit wigners
  Lorentz5Momentum pa1[2] = {part.momentum()-momenta[0],
			     part.momentum()-momenta[1]};
  for(unsigned int ix=0;ix<2;++ix) pa1[ix].rescaleMass();
  complex<InvEnergy2> bwa1[2] = {Resonance::BreitWignera1(pa1[0].mass2(),ma1_,ga1_)/sqr(ma1_),
				 Resonance::BreitWignera1(pa1[1].mass2(),ma1_,ga1_)/sqr(ma1_)};
  if(ichan>0) bwa1[ ichan == 0 ? 1 : 0 ] = ZERO;
  // compute the matrix element
  for(unsigned int ihel=0;ihel<3;++ihel) {
    LorentzVector<complex<Energy2> > pol[2] = {epsilon(vectors_[0][ihel],part.momentum(),pa1[0]),
					       epsilon(vectors_[0][ihel],part.momentum(),pa1[1])};
    for(unsigned int ohel=0;ohel<3;++ohel) {
      (*ME())(ihel,0,0,ohel) = Complex(gf1a1Pi_*ga1RhoPi_*(LorentzPolarizationVector(pol[0]*bwa1[0])-
							   LorentzPolarizationVector(pol[1]*bwa1[1])).dot(vectors_[1][ohel]));
    }
  } 
  // matrix element
  return ME()->contract(rho_).real();
}

void f1RhoPiPiDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":ga1RhoPi " << ga1RhoPi_/GeV << "\n";
  output << "newdef " << name() << ":gf1a1Pi "  << gf1a1Pi_*GeV  << "\n";
  for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
    output << "newdef    " << name() << ":maxWeight " << ix << " "
  	   << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
