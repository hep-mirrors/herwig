// -*- C++ -*-
//
// EtaPiPiPiDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiPiDecayer class.
//
#include "EtaPiPiPiDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/PDT/OneOffShellCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void EtaPiPiPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void EtaPiPiPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()||isize!=prefactor_.size()||
     isize!=charged_.size()||isize!=a_.size()||
     isize!=b_.size()||isize!=c_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in EtaPiPiPiDecayer::doinit()"
  			  << Exception::runerror;
  // external particles for the modes
  tPDVector outneut(3),outcharged(3);
  outneut[0]    = getParticleData(ParticleID::pi0);
  outneut[1]    = getParticleData(ParticleID::pi0);
  outcharged[0] = getParticleData(ParticleID::piplus);
  outcharged[1] = getParticleData(ParticleID::piminus);
  tPDPtr rho(getParticleData(113));
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr incoming = getParticleData(incoming_[ix]);
    outneut[2]    = getParticleData(outgoing_[ix]);
    outcharged[2] = getParticleData(outgoing_[ix]);
    // the pi+pi- mode
    if(charged_[ix]) {
      mode = new_ptr(PhaseSpaceMode(incoming,outcharged,maxWeight_[ix]));
    }
    // the pi0pi0 mode
    else {
      mode = new_ptr(PhaseSpaceMode(incoming,outneut,maxWeight_[ix]));
    }
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,rho,0,3,1,1,1,2));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.0);
    mode->addChannel(newChannel);
    addMode(mode);
  }
  resetIntermediate(rho,600.*MeV,600.*MeV);
}

int EtaPiPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
				 const tPDVector & children) const {
  if(children.size()!=3) return -1;
  unsigned int npi0(0),npip(0),npim(0); int id,iother(0);
  tPDVector::const_iterator pit = children.begin();
  for( ;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)           ++npip;
    else if(id==ParticleID::piminus)     ++npim;
    else if(id==ParticleID::pi0&&npi0<2) ++npi0;
    else iother=id;
  }
  bool charged;
  if(npim==1&&npip==1) {
    charged=true;
    if(npi0==1) iother=ParticleID::pi0;
  }
  else if(npi0==2) charged=false;
  else return -1;
  // find the mode
  id=parent->id();
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id==incoming_[ix]&&iother==outgoing_[ix]&&charged_[ix]==charged) 
      imode=ix;
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void EtaPiPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << charged_ << prefactor_ << a_ << b_ << c_  
     << maxWeight_;
}

void EtaPiPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> charged_ >> prefactor_ >> a_ >> b_ >> c_ 
     >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPiPiPiDecayer,DecayIntegrator>
describeHerwigEtaPiPiPiDecayer("Herwig::EtaPiPiPiDecayer", "HwSMDecay.so");

void EtaPiPiPiDecayer::Init() {

  static ClassDocumentation<EtaPiPiPiDecayer> documentation
    ("The EtaPiPiPiDecayer class performs the decay of a scalar meson to"
     " two pions and another meson using a simple paramterisation of the dalitz plot.",
     "The decay of eta to two pions follows \\cite{Beisert:2003zs,Gormley:1970qz,Tippens:2001fm}.",
     "%\\cite{Beisert:2003zs}\n"
     "\\bibitem{Beisert:2003zs}\n"
     "  N.~Beisert and B.~Borasoy,\n"
     "  %``Hadronic decays of eta and eta' with coupled channels,''\n"
     "  Nucl.\\ Phys.\\  A {\\bf 716}, 186 (2003)\n"
     "  [arXiv:hep-ph/0301058].\n"
     "  %%CITATION = NUPHA,A716,186;%%\n"
     "%\\cite{Gormley:1970qz}\n"
     "\\bibitem{Gormley:1970qz}\n"
     "  M.~Gormley, E.~Hyman, W.~Y.~Lee, T.~Nash, J.~Peoples, C.~Schultz and S.~Stein,\n"
     "   ``Experimental determination of the dalitz-plot distribution of the decays\n"
     "   eta $\\to$ pi+ pi- pi0 and eta $\\to$ pi+ pi- gamma, and the branching ratio\n"
     "  %eta $\\to$ pi+ pi- gamma/eta $\\to$ pi+,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 2}, 501 (1970).\n"
     "  %%CITATION = PHRVA,D2,501;%%\n"
     "%\\cite{Tippens:2001fm}\n"
     "\\bibitem{Tippens:2001fm}\n"
     "  W.~B.~Tippens {\\it et al.}  [Crystal Ball Collaboration],\n"
     "  %``Determination of the quadratic slope parameter in eta $\\to$ 3pi0 decay,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 87}, 192001 (2001).\n"
     "  %%CITATION = PRLTA,87,192001;%%\n"
     );

  static Command<EtaPiPiPiDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the decay mode (incoming, outgoing, charged/neutral pions, prefactor, a, b, c parameters and maximum weight",
     &EtaPiPiPiDecayer::setUpDecayMode, false);

  static Deleted<EtaPiPiPiDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<EtaPiPiPiDecayer> interfaceOutgoing
    ("Outgoing","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfaceCharged
    ("Charged","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfacePrefactor
    ("Prefactor","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfacea
    ("a","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfaceb
    ("b","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfacec
    ("c","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

}

void EtaPiPiPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}



double EtaPiPiPiDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the matrix element
  // compute the variables we need
  Lorentz5Momentum ps(part.momentum()-momenta[2]);
  ps.rescaleMass();
  Lorentz5Momentum pu(part.momentum()-momenta[0]);
  pu.rescaleMass();
  Lorentz5Momentum pt(part.momentum()-momenta[1]);
  pt.rescaleMass();
  Energy2 s(ps.mass2()),u(pu.mass2()),t(pt.mass2());
  Energy m34(0.5*(momenta[0].mass()+momenta[1].mass()));
  Energy msum(momenta[2].mass()+2.*m34);
  Energy Q(part.mass()-msum);
  Energy2 Mmm2((part.mass()-momenta[2].mass())*(part.mass()-momenta[2].mass()));
  // compute the variables
  double x(0.5*sqrt(3.)*(u-t)/part.mass()/Q),x2(x*x);
  double y(0.5*msum/part.mass()*(Mmm2-s)/m34/Q-1),y2(y*y);
  double me(prefactor_[imode()]*(1+a_[imode()]*y+b_[imode()]*y2+c_[imode()]*x2));
  if(me<0.) me=0.;
  (*ME())(0,0,0,0)=sqrt(me);
  return me;
}

InvEnergy EtaPiPiPiDecayer::threeBodydGammads(const int imodeb, const Energy2 q2,
					   const  Energy2 s, const Energy m1,
					   const Energy m2, const Energy m3) const {
  Energy q(sqrt(q2)),m34(m1+m2),msum(m34+m3),Q(q-msum);
  Energy2 Mmm2((q-m3)*(q-m3)),m12(m1*m1),m22(m2*m2),m32(m3*m3);
  double y(0.5*msum/q*(Mmm2-s)/m34/Q-1),y2(y*y);
  InvEnergy2 xfact=0.5*sqrt(3.)/q/Q;
  Energy2 xc(q2+m12+m22+m32-s);
  Energy rs(sqrt(s)),e2star(0.5*(s-m12+m22)/rs),e3star(0.5*(q2-s-m32)/rs);
  Energy e2sm(sqrt(e2star*e2star-m22)),e3sm(sqrt(e3star*e3star-m32));
  Energy2 a(2*e2star*e3star+m22+m32),b(2*e2sm*e3sm);
  Energy2 output=2*b*(1+a_[imodeb]*y+b_[imodeb]*y2+c_[imodeb]*xfact*xfact*(xc*xc))
    +c_[imodeb]*(-8.*xfact*xfact*xc*a*b
		 +4.*2*b*(3.*a*a+b*b)/3.*xfact*xfact);
  using Constants::pi;
  return output*prefactor_[imodeb]/256./pi/pi/pi/q2/q;
}


WidthCalculatorBasePtr 
EtaPiPiPiDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  int idout(0),id,imode(-1);
  unsigned int npi0(0),ix(0);
  ParticleMSet::const_iterator pit(dm.products().begin());
  for( ;pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::pi0&&npi0<2)                            ++npi0;
    else if(id!=ParticleID::piplus&&id!=ParticleID::piminus) idout=id;
  }
  if(npi0==1) idout=ParticleID::pi0;
  bool charged(npi0<2);
  id=dm.parent()->id();
  do {
    if(id==incoming_[ix]&&idout==outgoing_[ix]&&charged_[ix]==charged) 
      imode=ix;
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  Energy mpi;
  if(charged){mpi=getParticleData(ParticleID::piplus)->mass();}
  else{mpi=getParticleData(ParticleID::pi0)->mass();}
  Energy m[3]={mpi,mpi,getParticleData(outgoing_[imode])->mass()};
  WidthCalculatorBasePtr 
    temp(new_ptr(ThreeBodyAllOn1IntegralCalculator<EtaPiPiPiDecayer>
  		 (1,-1000.*MeV,ZERO,0.0,*this,imode,m[0],m[1],m[2])));
  if(outgoing_[imode]==ParticleID::eta) {
    tcGenericMassGeneratorPtr test;
    tGenericMassGeneratorPtr massptr;
    if(getParticleData(outgoing_[imode])->massGenerator()) {
      test=dynamic_ptr_cast<tcGenericMassGeneratorPtr>
  	(getParticleData(outgoing_[imode])->massGenerator());
      massptr=const_ptr_cast<tGenericMassGeneratorPtr>(test);
    }
    if(massptr) {
      massptr->init();
      return new_ptr(OneOffShellCalculator(3,temp,massptr,ZERO));
    }
  }
  return temp;
} 
  
void EtaPiPiPiDecayer::dataBaseOutput(ofstream & output,
				      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix]   << " "
	   << outgoing_[ix]  << " " << charged_[ix]  << " " << prefactor_[ix]  << " "
	   << a_[ix]  << " " << b_[ix]  << " " << c_[ix]  << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string EtaPiPiPiDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "First outgoing particle with id " + std::to_string(out) + "does not have spin 0";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  bool charge = stoi(stype);
  // get the couplings
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double a = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double b = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double c = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_ .push_back(in);
  outgoing_ .push_back(out);
  charged_  .push_back(charge);
  prefactor_.push_back(g);
  a_        .push_back(a);
  b_        .push_back(b);
  c_        .push_back(c);
  maxWeight_.push_back(wgt);
  return "";
}
