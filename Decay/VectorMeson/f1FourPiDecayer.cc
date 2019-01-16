// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the f1FourPiDecayer class.
//

#include "f1FourPiDecayer.h"
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

f1FourPiDecayer::f1FourPiDecayer() :
  gRhoPiPi_(6.), ga1RhoPi_(4.8*GeV), gf1a1Pi_(9.77/GeV), maxWeight_({1.,1.}) {
  generateIntermediates(true);
}

// normally not implemented but do it here to get rid of the a_1 which
// can't be on-shell
ParticleVector f1FourPiDecayer::decay(const Particle & parent,
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

IBPtr f1FourPiDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr f1FourPiDecayer::fullclone() const {
  return new_ptr(*this);
}

void f1FourPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << gRhoPiPi_ << ounit(ga1RhoPi_,GeV) << ounit(gf1a1Pi_,1./GeV)
     << ounit(ma1_,GeV) << ounit(ga1_,GeV)
     << ounit(mrho_,GeV) << ounit(grho_,GeV) << maxWeight_;
}

void f1FourPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> gRhoPiPi_ >> iunit(ga1RhoPi_,GeV) >> iunit(gf1a1Pi_,1./GeV)
     >> iunit(ma1_,GeV) >> iunit(ga1_,GeV)
     >> iunit(mrho_,GeV) >> iunit(grho_,GeV) >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<f1FourPiDecayer,DecayIntegrator>
describeHerwigf1FourPiDecayer("Herwig::f1FourPiDecayer", "HwVMDecay.so");

void f1FourPiDecayer::Init() {

  static ClassDocumentation<f1FourPiDecayer> documentation
    ("The f1FourPiDecayer class implements a simple model for "
     "f1 -> pipirho via an intermediate a_1");

  static ParVector<f1FourPiDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "Maximum weights for the decays",
     &f1FourPiDecayer::maxWeight_, 2, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<f1FourPiDecayer,double> interfacegRhoPiPi
    ("gRhoPiPi",
     "The coupling of the rho to two pions",
     &f1FourPiDecayer::gRhoPiPi_, 6., 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<f1FourPiDecayer,Energy> interfacega1RhoPi
    ("ga1RhoPi",
     "Coupling of the a_1 to rho pi",
     &f1FourPiDecayer::ga1RhoPi_, GeV, 4.8*GeV, 0.*GeV, 20.*GeV,
     false, false, Interface::limited);
  
  static Parameter<f1FourPiDecayer,InvEnergy> interfacegf1a1Pi
    ("gf1a1Pi",
     "Coupling of f_1 to a_1 pi",
     &f1FourPiDecayer::gf1a1Pi_, 1./GeV, 1.0/GeV, 0.0/GeV, 10.0/GeV,
     false, false, Interface::limited);

}

void f1FourPiDecayer::doinit() {
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
  // decay mode f_1 -> pi+ pi- pi+ pi-
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(f1,{pip,pim,pip,pim},maxWeight_[0]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rho0,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,1,1,rho0,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,2,1,rho0,2,1,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,3,1,rho0,2,1,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rho0,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,1,1,rho0,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,4,1,rho0,2,1,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,3,1,rho0,2,1,2,2));
  addMode(mode);
  // decay mode f_1 -> pi+ pi0 pi- pi0
  mode = new_ptr(PhaseSpaceMode(f1,{pip,pi0,pim,pi0},maxWeight_[0]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rhop,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,1,1,rhop,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rhop,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,1,1,rhop,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,4,1,rhom,2,1,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,3,1,rhom,2,1,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,2,1,rhom,2,1,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,3,1,rhom,2,1,2,4));
  addMode(mode);
  // masses of intermediates
  ma1_ = a1p->mass();
  ga1_ = a1p->width(); 
  mrho_ = rhop->mass();
  grho_ = rhop->width(); 
}

void f1FourPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<maxWeight_.size();++ix)
      maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

int f1FourPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  if(children.size()!=4 || parent->id()!=ParticleID::f_1) return -1;
  // check the pions
  int npi0(0),npiplus(0),npiminus(0);
  for(auto child : children) {
    int idtemp= child->id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
  }
  cc = false;
  // f_1 -> 2pi+ 2pi-mode
  if(npiplus==2 && npiminus==2)                    return 0;
  // f_1 -> pi+ pi- pi0 pi0 mode
  else if  (npiplus ==1 && npiminus==1 && npi0==2) return 1;
  // not found
  return -1;
}

void f1FourPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vector_,const_ptr_cast<tPPtr>(&part),
  					incoming,true,false);
  // set up the spin information for the pions
  for(unsigned int ix=0;ix<4;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double f1FourPiDecayer::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  // polarization vectors for incoming
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vector_,rho_,
  					       const_ptr_cast<tPPtr>(&part),
  					       incoming,false);
  }
  // breit wigners
  Lorentz5Momentum pa1[4] = {part.momentum()-momenta[0],part.momentum()-momenta[1],
			     part.momentum()-momenta[2],part.momentum()-momenta[3]};
  complex<InvEnergy2> bwa1[4];
  for(unsigned int ix=0;ix<4;++ix) {
    pa1[ix].rescaleMass();
    bwa1[ix] = Resonance::BreitWignera1(pa1[ix].mass2(),ma1_,ga1_)/sqr(ma1_);
  }
  // compute the matrix element (as a current and then contract)
  LorentzVector<complex<InvEnergy> > current;
  // decay mode f_1 -> pi+ pi- pi+pi-
  double sym(0.5);
  if(imode()==0) {
    sym=0.25;
    complex<InvEnergy2> bwrho[4] =
      {Resonance::BreitWignerPWave((momenta[0]+momenta[1]).m2(),mrho_,grho_,momenta[0].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[0]+momenta[3]).m2(),mrho_,grho_,momenta[0].mass(),momenta[3].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[1]).m2(),mrho_,grho_,momenta[2].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[3]).m2(),mrho_,grho_,momenta[2].mass(),momenta[3].mass())/sqr(mrho_)};
    if(ichan<=0) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rho0,2,3,2,4));
      current += bwa1[0]*bwrho[3]*epsilon(part.momentum(),pa1[0],momenta[2]-momenta[3]);
    }
    if(ichan<0||ichan==1) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,1,1,rho0,2,3,2,4));
      current -= bwa1[1]*bwrho[3]*epsilon(part.momentum(),pa1[1],momenta[2]-momenta[3]);
    }
    if(ichan<0||ichan==2) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,2,1,rho0,2,1,2,4));
      current += bwa1[2]*bwrho[1]*epsilon(part.momentum(),pa1[2],momenta[0]-momenta[3]);
    }
    if(ichan<0||ichan==3) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,3,1,rho0,2,1,2,4));
      current -= bwa1[1]*bwrho[1]*epsilon(part.momentum(),pa1[1],momenta[0]-momenta[3]);
    }
    if(ichan<0||ichan==4) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rho0,2,3,2,2));
      current += bwa1[0]*bwrho[2]*epsilon(part.momentum(),pa1[0],momenta[2]-momenta[1]);
    }
    if(ichan<0||ichan==5) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,1,1,rho0,2,3,2,2));
      current -= bwa1[3]*bwrho[2]*epsilon(part.momentum(),pa1[3],momenta[2]-momenta[1]);
    }
    if(ichan<0||ichan==6) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,4,1,rho0,2,1,2,2));
      current += bwa1[2]*bwrho[0]*epsilon(part.momentum(),pa1[2],momenta[0]-momenta[1]);
    }
    if(ichan<0||ichan==7) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,3,1,rho0,2,1,2,2));
      current -= bwa1[3]*bwrho[0]*epsilon(part.momentum(),pa1[3],momenta[0]-momenta[1]);
    }
  }
  // decay mode f_1 -> pi+ pi0 pi- pi0
  else {
    complex<InvEnergy2> bwrho[4] =
      {Resonance::BreitWignerPWave((momenta[0]+momenta[1]).m2(),mrho_,grho_,momenta[0].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[0]+momenta[3]).m2(),mrho_,grho_,momenta[0].mass(),momenta[3].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[1]).m2(),mrho_,grho_,momenta[2].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[3]).m2(),mrho_,grho_,momenta[2].mass(),momenta[3].mass())/sqr(mrho_)};
    double f1 = (momenta[2].mass2()-momenta[3].mass2())/sqr(mrho_);
    double f2 = 1+f1;
    f1 = 1.-f1;
    if(ichan<=0) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rhop,2,3,2,4));
      current += bwa1[0]*bwrho[3]*epsilon(part.momentum(),pa1[0],f1*momenta[2]-f2*momenta[3]);
    }
    if(ichan<0||ichan==1) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,1,1,rhop,2,3,2,4));
      current -= bwa1[1]*bwrho[3]*epsilon(part.momentum(),pa1[1],f1*momenta[2]-f2*momenta[3]);
    }
    if(ichan<0||ichan==2) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rhop,2,3,2,2));
      current += bwa1[0]*bwrho[2]*epsilon(part.momentum(),pa1[0],f1*momenta[2]-f2*momenta[1]);
    }
    if(ichan<0||ichan==3) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,1,1,rhop,2,3,2,2));
      current -= bwa1[3]*bwrho[2]*epsilon(part.momentum(),pa1[3],f1*momenta[2]-f2*momenta[1]);
    }
    if(ichan<0||ichan==4) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,4,1,rhom,2,1,2,2));
      current += bwa1[2]*bwrho[0]*epsilon(part.momentum(),pa1[2],f1*momenta[0]-f2*momenta[1]);
    }
    if(ichan<0||ichan==5) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,3,1,rhom,2,1,2,2));
      current -= bwa1[3]*bwrho[0]*epsilon(part.momentum(),pa1[3],f1*momenta[0]-f2*momenta[1]);
    }
    if(ichan<0||ichan==6) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,2,1,rhom,2,1,2,4));
      current += bwa1[2]*bwrho[1]*epsilon(part.momentum(),pa1[2],f1*momenta[0]-f2*momenta[3]);
    }
    if(ichan<0||ichan==7) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,3,1,rhom,2,1,2,4));
      current -= bwa1[1]*bwrho[1]*epsilon(part.momentum(),pa1[1],f1*momenta[0]-f2*momenta[3]);
    }
  }
  // contract the current
  double pre = gRhoPiPi_*gf1a1Pi_*ga1RhoPi_;
  for(unsigned int ihel=0;ihel<3;++ihel)
    (*ME())(ihel,0,0,0,0) = pre*current.dot(vector_[ihel])*part.mass();
  // matrix element
  return sym*ME()->contract(rho_).real();
}

void f1FourPiDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":gRhoPiPi " << gRhoPiPi_     << "\n";
  output << "newdef " << name() << ":ga1RhoPi " << ga1RhoPi_/GeV << "\n";
  output << "newdef " << name() << ":gf1a1Pi "  << gf1a1Pi_*GeV  << "\n";
  for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
    output << "newdef    " << name() << ":maxWeight " << ix << " "
  	   << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
