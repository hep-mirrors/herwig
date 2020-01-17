// -*- C++ -*-
//
// SMHiggsFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsFermionsDecayer class.
//

#include "SMHiggsFermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "Herwig/Utilities/Maths.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/ShowerAlpha.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMHiggsFermionsDecayer::SMHiggsFermionsDecayer() :
  CF_(4./3.), NLO_(false) {
  _maxwgt.resize(9);
  _maxwgt[0]=0.;
  _maxwgt[1]=0;		
  _maxwgt[2]=0;		
  _maxwgt[3]=0.0194397;	
  _maxwgt[4]=0.463542;	
  _maxwgt[5]=0.;		
  _maxwgt[6]=6.7048e-09; 
  _maxwgt[7]=0.00028665; 
  _maxwgt[8]=0.0809643;  
}

void SMHiggsFermionsDecayer::doinit() {
  PerturbativeDecayer::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "SMHiggsFermionsDecayer needs the StandardModel class"
			  << " to be either the Herwig one or a class inheriting"
			  << " from it";
  _hvertex = hwsm->vertexFFH();
  // make sure they are initialized
  _hvertex->init();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
  // set up the decay modes
  vector<double> wgt(0);
  unsigned int imode=0;
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  int iy;
  extpart[0]=higgs;
  for(unsigned int istep=0;istep<11;istep+=10) {
    for(unsigned ix=1;ix<7;++ix) {
      if(istep<10||ix%2!=0) {
	iy = ix+istep;
	extpart[1]=getParticleData( iy);
	extpart[2]=getParticleData(-iy);
	mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
	addMode(mode,_maxwgt[imode],wgt);
	++imode;
      }
    }
  }
//   Energy quarkMass = getParticleData(ParticleID::b )->mass();
//   Energy higgsMass = getParticleData(ParticleID::h0)->mass();
//   double mu = quarkMass/higgsMass;
//   double beta = sqrt(1.-4.*sqr(mu));
//   double beta2 = sqr(beta);
//   double aS = SM().alphaS(sqr(higgsMass));
//   double L = log((1.+beta)/(1.-beta));
//   cerr << "testing " << beta << " " << mu << "\n";
//   cerr << "testing " << aS << " " << L << "\n";
//   double fact = 
//     6.-0.75*(1.+beta2)/beta2+12.*log(mu)-8.*log(beta)
//     +(5./beta-2.*beta+0.375*sqr(1.-beta2)/beta2/beta)*L
//     +(1.+beta2)/beta*(4.*L*log(0.5*(1.+beta)/beta)
// 		      -2.*log(0.5*(1.+beta))*log(0.5*(1.-beta))
// 		      +8.*Herwig::Math::ReLi2((1.-beta)/(1.+beta))
// 		      -4.*Herwig::Math::ReLi2(0.5*(1.-beta)));
//   cerr << "testing correction " 
//        << 1.+4./3.*aS/Constants::twopi*fact
//        << "\n"; 
//   double real = 4./3.*aS/Constants::twopi*
//     (8.-0.75*(1.+beta2)/beta2+8.*log(mu)-8.*log(beta)
//      +(3./beta+0.375*sqr(1.-beta2)/pow(beta,3))*L
//      +(1.+beta2)/beta*(-0.5*sqr(L)+4.*L*log(0.5*(1.+beta))
// 		       -2.*L*log(beta)-2.*log(0.5*(1.+beta))*log(0.5*(1.-beta))
// 		       +6.*Herwig::Math::ReLi2((1.-beta)/(1.+beta))
// 		       -4.*Herwig::Math::ReLi2(0.5*(1.-beta))
// 		       -2./3.*sqr(Constants::pi)));
//   double virt = 4./3.*aS/Constants::twopi*
//     (-2.+4.*log(mu)+(2./beta-2.*beta)*L
//      +(1.+beta2)/beta*(0.5*sqr(L)-2.*L*log(beta)+2.*sqr(Constants::pi)/3.
// 		       +2.*Herwig::Math::ReLi2((1.-beta)/(1.+beta))));
//   cerr << "testing real " << real << "\n";
//   cerr << "testing virtual " << virt << "\n";
//   cerr << "testing total no mb corr " << 1.+real+virt << "\n";
//   cerr << "testing total    mb corr " << 1.+real+virt +(8./3. - 2.*log(sqr(mu)))*aS/Constants::pi << "\n";
//   InvEnergy2 Gf = 1.166371e-5/GeV2;
//   Gf = sqrt(2.)*4*Constants::pi*SM().alphaEM(sqr(higgsMass))/8./SM().sin2ThetaW()/
//     sqr(getParticleData(ParticleID::Wplus)->mass());
//   cerr << "testing GF " << Gf*GeV2 << "\n";
//   Energy LO = (3./8./Constants::pi)*sqrt(2)*sqr(quarkMass)*Gf*higgsMass*beta*beta*beta;
//   cerr << "testing LO " << LO/GeV << "\n";
//   cerr << "testing quark mass " << quarkMass/GeV << "\n";
//   cerr << "testing gamma " << (1.+real+virt)*LO/MeV << "\n";
}
  
bool SMHiggsFermionsDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  if(parent->id()!=ParticleID::h0||children.size()!=2) return false;
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if(id1==-id2&&(abs(id1)<=6||(abs(id1)>=11&&abs(id1)<=16)))
    return true;
  else
    return false;
}

ParticleVector SMHiggsFermionsDecayer::decay(const Particle & parent,
					     const tPDVector & children) const {
  // id's of the decaying particles
  tPDVector::const_iterator pit(children.begin());
  int id1((**pit).id());
  int imode=-1;
  if(abs(id1)<=6)                     imode = abs(id1)-1;
  else if(abs(id1)>=11&&abs(id1)<=16) imode = (abs(id1)-11)/2+6;
  ParticleVector output(generate(false,false,imode,parent));
  // set up the colour flow
  if(output[0]->hasColour())      output[0]->antiColourNeighbour(output[1]);
  else if(output[1]->hasColour()) output[1]->antiColourNeighbour(output[0]);
  return output;
}


void SMHiggsFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _maxwgt << _hvertex << NLO_;
}

void SMHiggsFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _maxwgt >> _hvertex >> NLO_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHiggsFermionsDecayer,PerturbativeDecayer>
describeHerwigSMHiggsFermionsDecayer("Herwig::SMHiggsFermionsDecayer", "HwPerturbativeHiggsDecay.so");

void SMHiggsFermionsDecayer::Init() {

  static ClassDocumentation<SMHiggsFermionsDecayer> documentation
    ("The SMHiggsFermionsDecayer class implements the decat of the Standard Model"
     " Higgs boson to the Standard Model fermions");

  static ParVector<SMHiggsFermionsDecayer,double> interfaceMaxWeights
    ("MaxWeights",
     "Maximum weights for the various decays",
     &SMHiggsFermionsDecayer::_maxwgt, 9, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<SMHiggsFermionsDecayer,bool> interfaceNLO
    ("NLO",
     "Whether to return the LO or NLO result",
     &SMHiggsFermionsDecayer::NLO_, false, false, false);
  static SwitchOption interfaceNLOLO
    (interfaceNLO,
     "No",
     "Leading-order result",
     false);
  static SwitchOption interfaceNLONLO
    (interfaceNLO,
     "Yes",
     "NLO result",
     true);

}

// return the matrix element squared
double SMHiggsFermionsDecayer::me2(const int, const Particle & part,
				   const ParticleVector & decay,
				   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _swave = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(_rho);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					  incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[ianti],outgoing);
  Energy2 scale(sqr(part.mass()));
  unsigned int ifm,ia;
  for(ifm=0;ifm<2;++ifm) {
    for(ia=0;ia<2;++ia) {
      if(iferm>ianti)
	(*ME())(0,ia,ifm)=_hvertex->evaluate(scale,_wave[ia],
					  _wavebar[ifm],_swave);
      else
	(*ME())(0,ifm,ia)=_hvertex->evaluate(scale,_wave[ia],
					  _wavebar[ifm],_swave);
    }
  }
  int id = abs(decay[0]->id());
  double output=(ME()->contract(_rho)).real()*UnitRemoval::E2/scale;
  if(id <=6) output*=3.;
  // test of the partial width
//   Ptr<Herwig::StandardModel>::transient_const_pointer 
//     hwsm=dynamic_ptr_cast<Ptr<Herwig::StandardModel>::transient_const_pointer>(standardModel());
//   double g2(hwsm->alphaEM(scale)*4.*Constants::pi/hwsm->sin2ThetaW());
//   Energy mass(hwsm->mass(scale,decay[0]->dataPtr())),
//     mw(getParticleData(ParticleID::Wplus)->mass());
//   double beta(sqrt(1.-4.*decay[0]->mass()*decay[0]->mass()/scale));
//   cerr << "testing alpha " << hwsm->alphaEM(scale) << "\n";
//   Energy test(g2*mass*mass*beta*beta*beta*part.mass()/32./Constants::pi/mw/mw);
//   if(abs(decay[0]->id())<=6){test *=3.;}
//   cout << "testing the answer " << output << "     " 
//        << test/GeV
//        << endl;
  // leading-order result
  if(!NLO_) return output;
  // fermion mass
  Energy particleMass = decay[0]->dataPtr()->mass();
  // check decay products coloured, otherwise return
  if(!decay[0]->dataPtr()->coloured()||
     particleMass==ZERO) return output;
  // inital masses, couplings  etc
  // higgs mass
  mHiggs_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mHiggs_));
  // reduced mass
  mu_  = particleMass/mHiggs_;
  mu2_ = sqr(mu_);
  // generate y
  double yminus = 0.; 
  double yplus  = 1.-2.*mu_*(1.-mu_)/(1.-2*mu2_);
  double y = yminus + UseRandom::rnd()*(yplus-yminus);
  //generate z for D31,2
  double v  = sqrt(sqr(2.*mu2_+(1.-2.*mu2_)*(1.-y))-4.*mu2_)/(1.-2.*mu2_)/(1.-y);
  double zplus  = (1.+v)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
  double zminus = (1.-v)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
  double z = zminus + UseRandom::rnd()*(zplus-zminus);
  // map y,z to x1,x2 for both possible emissions
  double x2 = 1. - y*(1.-2.*mu2_);
  double x1 = 1. - z*(x2-2.*mu2_);
  //get the dipoles
  InvEnergy2 D1 = dipoleSubtractionTerm( x1, x2); 
  InvEnergy2 D2 = dipoleSubtractionTerm( x2, x1); 
  InvEnergy2 dipoleSum = abs(D1) + abs(D2);
  //jacobian
  double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
  //calculate real
  Energy2 realPrefactor = 0.25*sqr(mHiggs_)*sqr(1.-2.*mu2_)
    /sqrt(calculateLambda(1,mu2_,mu2_))/sqr(Constants::twopi);
  InvEnergy2 realEmission = 4.*Constants::pi*aS_*CF_*calculateRealEmission( x1, x2);
  // calculate the virtual
  double virtualTerm = calculateVirtualTerm();
  // running mass correction
  virtualTerm += (8./3. - 2.*log(mu2_))*aS_/Constants::pi;
  //answer = (born + virtual + real)/born * LO
  output *= 1. + virtualTerm + 2.*jac*realPrefactor*(realEmission*abs(D1)/dipoleSum  - D1);
  // return the answer
  return output;
}

void SMHiggsFermionsDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the PerturbativeDecayer base class
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    os << "newdef " << name() << ":MaxWeights " << ix << " "
	   << _maxwgt[ix] << "\n";
  }
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}

void SMHiggsFermionsDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _maxwgt[ix] = mode(ix)->maxWeight();
    }
  }
}


//calculate lambda
double SMHiggsFermionsDecayer::calculateLambda(double x, double y, double z) const{
  return sqr(x)+sqr(y)+sqr(z)-2.*x*y-2.*x*z-2.*y*z;
}

//calculates the dipole subtraction term for x1, D31,2 (Dij,k),
// 2 is the spectator anti-fermion and 3 is the gluon
InvEnergy2 SMHiggsFermionsDecayer::
dipoleSubtractionTerm(double x1, double x2) const{
  InvEnergy2 commonPrefactor = CF_*8.*Constants::pi*aS_/sqr(mHiggs_);
  return commonPrefactor/(1.-x2)*
    (2.*(1.-2.*mu2_)/(2.-x1-x2)- 
     sqrt((1.-4.*mu2_)/(sqr(x2)-4.*mu2_))*
     (x2-2.*mu2_)*(2.+(x1-1.)/(x2-2.*mu2_)+2.*mu2_/(1.-x2))/(1.-2.*mu2_));
}

//return ME for real emission
InvEnergy2 SMHiggsFermionsDecayer::
calculateRealEmission(double x1, double x2) const {
  InvEnergy2 prefactor = 2./sqr(mHiggs_)/(1.-4.*mu2_);
  return prefactor*(2. + (1.-x1)/(1.-x2) + (1.-x2)/(1.-x1) 
                    + 2.*(1.-2.*mu2_)*(1.-4.*mu2_)/(1.-x1)/(1.-x2)
                    - 2.*(1.-4.*mu2_)*(1./(1.-x2)+1./(1.-x1)) 
                    - 2.*mu2_*(1.-4.*mu2_)*(1./sqr(1.-x2)+1./sqr(1.-x1)));
}

double SMHiggsFermionsDecayer::
calculateVirtualTerm() const {
  // logs and prefactors
  double beta = sqrt(1.-4.*mu2_);
  double L = log((1.+beta)/(1.-beta));
  double prefactor = CF_*aS_/Constants::twopi;
  // non-singlet piece
  double nonSingletTerm = calculateNonSingletTerm(beta, L);
  double virtualTerm = 
    -2.+4.*log(mu_)+(2./beta - 2.*beta)*L 
    + (2.-4.*mu2_)/beta*(0.5*sqr(L) - 2.*L*log(beta)
			 + 2.*Herwig::Math::ReLi2((1.-beta)/(1.+beta)) 
			 + 2.*sqr(Constants::pi)/3.);
  double iEpsilonTerm = 
    2.*(3.-sqr(Constants::pi)/2. + 0.5*log(mu2_) - 1.5*log(1.-2.*mu2_)
	-(1.-2.*mu2_)/beta*(0.5*sqr(L)+sqr(Constants::pi)/6.
			    -2.*L*log(1.-2.*mu2_))
	+ nonSingletTerm);
  return prefactor*(virtualTerm+iEpsilonTerm);
}

//non-singlet piece of I(epsilon) insertion operator
double SMHiggsFermionsDecayer::
calculateNonSingletTerm(double beta, double L) const {
  return  1.5*log(1.-2.*mu2_)  
    + (1.-2.*mu2_)/beta*(- 2.*L*log(4.*(1.-2.*mu2_)/sqr(1.+beta))+
			 + 2.*Herwig::Math::ReLi2(sqr((1.-beta)/(1.+beta)))
			 - 2.*Herwig::Math::ReLi2(2.*beta/(1.+beta)) 
			 - sqr(Constants::pi)/6.) 
    + log(1.-mu_) 
    - 2.*log(1.-2.*mu_) 
    - 2.*mu2_/(1.-2.*mu2_)*log(mu_/(1.-mu_))
    - mu_/(1.-mu_)
    + 2.*mu_*(2*mu_-1.)/(1.-2.*mu2_)
    + 0.5*sqr(Constants::pi);
}

double SMHiggsFermionsDecayer::matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
						  const ParticleVector & decay3, MEOption,
						  ShowerInteraction inter) {
  mHiggs_ = inpart.mass();
  mu_ = decay2[0]->mass()/mHiggs_;
  mu2_ = sqr(mu_);
  double x1 = 2.*decay3[0]->momentum().t()/mHiggs_;
  double x2 = 2.*decay3[1]->momentum().t()/mHiggs_;
  double pre = inter==ShowerInteraction::QCD ? CF_ : sqr(double(decay2[0]->dataPtr()->iCharge())/3.);
  return pre*calculateRealEmission(x1,x2)*4.*Constants::pi*sqr(mHiggs_);
}
