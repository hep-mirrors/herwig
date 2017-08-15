// -*- C++ -*-
//
// SMHiggsFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
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
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMHiggsFermionsDecayer::SMHiggsFermionsDecayer() :
  CF_(4./3.), pTmin_(1.*GeV), NLO_(false) {
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
  gluon_ = getParticleData(ParticleID::g);
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
  os << _maxwgt << _hvertex << NLO_
     << alphaS_ << gluon_ << ounit( pTmin_, GeV );
}

void SMHiggsFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _maxwgt >> _hvertex >> NLO_
     >> alphaS_ >> gluon_ >> iunit( pTmin_, GeV );
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
  
  static Reference<SMHiggsFermionsDecayer,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "The object calculating the strong coupling constant",
     &SMHiggsFermionsDecayer::alphaS_, false, false, true, false, false);

  static Parameter<SMHiggsFermionsDecayer, Energy> interfacePtMin
    ("minpT",
     "The pt cut on hardest emision generation",
     &SMHiggsFermionsDecayer::pTmin_, GeV, 1.*GeV, 0*GeV, 100000.0*GeV,
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
  InvEnergy2 realEmission = calculateRealEmission( x1, x2);
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

RealEmissionProcessPtr SMHiggsFermionsDecayer::
generateHardest(RealEmissionProcessPtr born) {
  // check coloured
  if(!born->bornOutgoing()[0]->dataPtr()->coloured()) return RealEmissionProcessPtr();
  assert(born->bornOutgoing().size()==2);
  // extract required info
  higgs_ = born->bornIncoming()[0];
  partons_.resize(2);
  quark_.resize(2);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    partons_[ix] = born->bornOutgoing()[ix]->dataPtr();
    quark_[ix]   = born->bornOutgoing()[ix]->momentum();
    quark_[ix].setMass(partons_[ix]->mass());
  }
  bool order = partons_[0]->id()<0;
  if(order) {
    swap(partons_[0]   ,partons_[1]   );
    swap(quark_[0]     ,quark_[1]     );
  }
  gauge_.setMass(0.*MeV);
  // Get the Higgs boson mass.
  mh2_ = (quark_[0] + quark_[1]).m2();
  mHiggs_ = sqrt(mh2_);
  aS_ = SM().alphaS(sqr(mHiggs_));
  Energy particleMass = partons_[0]->mass();
  mu_  = particleMass/mHiggs_;
  mu2_ = sqr(mu_);
  // Generate emission and set _quark[0,1] and _gauge to be the 
  // momenta of q, qbar and g after the hardest emission:
  if(!getEvent()) {
    born->pT()[ShowerInteraction::QCD] = pTmin_;
    return born;
  }
  // Ensure the energies are greater than the constituent masses:
  for (int i=0; i<2; i++) {
    if (quark_[i].e() < partons_[i]->constituentMass()) return RealEmissionProcessPtr();
  }
  if (gauge_.e()    < gluon_     ->constituentMass()) return RealEmissionProcessPtr();
  // set masses
  quark_[0].setMass( partons_[0]->mass() );
  quark_[1].setMass( partons_[1]->mass() );
  gauge_   .setMass( ZERO );
  // assign the emitter based on evolution scales
  unsigned int iemitter   = quark_[0]*gauge_ > quark_[1]*gauge_ ? 2 : 1;
  unsigned int ispectator = iemitter==1                         ? 1 : 2;
  // create new partices and insert
  PPtr hboson = higgs_->dataPtr()->produceParticle(higgs_->momentum());
  born->incoming().push_back(hboson);
  PPtr newq = partons_[0]->produceParticle(quark_[0]);
  PPtr newa = partons_[1]->produceParticle(quark_[1]);
  PPtr newg = gluon_->produceParticle(gauge_);
  // make colour connections
  newg->colourNeighbour(newq);
  newa->colourNeighbour(newg);
  // insert in output structure
  if(!order) {
    born->outgoing().push_back(newq);
    born->outgoing().push_back(newa);
  }
  else {
    born->outgoing().push_back(newa);
    born->outgoing().push_back(newq);
    swap(iemitter,ispectator);
  }
  born->outgoing().push_back(newg);
  born->emitter  (iemitter  );
  born->spectator(ispectator);
  born->emitted  (3);
  born->pT()[ShowerInteraction::QCD] = pT_;
  // return process
  born->interaction(ShowerInteraction::QCD);
  return born;
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
  InvEnergy2 prefactor = CF_*8.*Constants::pi*aS_/sqr(mHiggs_)/(1.-4.*mu2_);
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

bool SMHiggsFermionsDecayer::getEvent() {
  // max pT
  Energy pTmax = 0.5*sqrt(mh2_);
  // Define over valued y_max & y_min according to the associated pt_min cut.
  double ymax  =  acosh(pTmax/pTmin_);
  double ymin  = -ymax;
  // pt of the emmission
  pT_ = pTmax;
  // prefactor
  double overEst = 4.;
  double prefactor = overEst*alphaS_->overestimateValue()*4./3./Constants::twopi;
  // loop to generate the pt and rapidity
  bool reject;
  
  //arrays to hold the temporary  probabilities whilst the for loop progresses
  double probTemp[2][2]={{0.,0.},{0.,0.}};
  probTemp[0][0]=probTemp[0][1]=probTemp[1][0]=probTemp[1][1]=0.;
  double x1Solution[2][2] = {{0.,0.},{0.,0.}};
  double x2Solution[2][2] = {{0.,0.},{0.,0.}};
  double x3Solution[2]    =  {0.,0.};
  Energy pT[2]            =  {pTmax,pTmax};
  double yTemp[2]         =  {0.,0.};
  for(int i=0; i<2; i++) {
    do {
      // reject the emission
      reject = true;
      // generate pt
      pT[i] *= pow(UseRandom::rnd(),1./(prefactor*(ymax-ymin)));
      Energy2 pT2 = sqr(pT[i]);
      if(pT[i]<pTmin_) {
        pT[i] = -GeV;
        break;
      }
      // generate y
      yTemp[i] = ymin + UseRandom::rnd()*(ymax-ymin);
      //generate x3 & x1 from pT & y
      double x1Plus  = 1.;
      double x1Minus = 2.*mu_;
      x3Solution[i] = 2.*pT[i]*cosh(yTemp[i])/mHiggs_;
      // prefactor
      Energy2 weightPrefactor = mh2_/16./sqr(Constants::pi)/sqrt(1.-4.*mu2_);
      weightPrefactor /= prefactor;
      // calculate x1 & x2 solutions
      Energy4 discrim2 = (sqr(x3Solution[i]*mHiggs_) - 4.*pT2)*
        (mh2_*(x3Solution[i]-1.)*(4.*mu2_+x3Solution[i]-1.)-4.*mu2_*pT2);
      //check discriminant2 is > 0
      if( discrim2 < ZERO) continue;
      Energy2 discriminant = sqrt(discrim2);
      Energy2 fact1 = 3.*mh2_*x3Solution[i]-2.*mh2_+2.*pT2*x3Solution[i]-4.*pT2-mh2_*sqr(x3Solution[i]);
      Energy2 fact2 = 2.*mh2_*(x3Solution[i]-1.)-2.*pT2;
      // two solns for x1
      x1Solution[i][0] = (fact1 + discriminant)/fact2;
      x1Solution[i][1] = (fact1  - discriminant)/fact2;
      x2Solution[i][0] = 2.-x3Solution[i]-x1Solution[i][0];
      x2Solution[i][1] = 2.-x3Solution[i]-x1Solution[i][1];
      bool found = false;
      for(unsigned int j=0;j<2;++j) {
        if(x1Solution[i][j]>=x1Minus && x1Solution[i][j]<=x1Plus &&
           checkZMomenta(x1Solution[i][j], x2Solution[i][j], x3Solution[i], yTemp[i], pT[i])) {
          InvEnergy2 D1 = dipoleSubtractionTerm( x1Solution[i][j], x2Solution[i][j]); 
          InvEnergy2 D2 = dipoleSubtractionTerm( x2Solution[i][j], x1Solution[i][j]);
	  double dipoleFactor = abs(D1)/(abs(D1) + abs(D2));
	  probTemp[i][j] = weightPrefactor*pT[i]*dipoleFactor*
            calculateJacobian(x1Solution[i][j], x2Solution[i][j], pT[i])*
            calculateRealEmission(x1Solution[i][j], x2Solution[i][j]);
          
          found = true;
        }
        else {
          probTemp[i][j] = 0.;
        }
      }
      if(!found) continue;
      // alpha S piece
      double wgt = (probTemp[i][0]+probTemp[i][1])*alphaS_->value(sqr(pT[i]))/aS_;
      // matrix element weight
      reject = UseRandom::rnd()>wgt;
  }
    while(reject);
  } //end of emitter for loop
  // no emission
  if(pT[0]<ZERO&&pT[1]<ZERO) return false;
  //pick the spectator and x1 x2 values
  double x1,x2,y;
  //particle 1 emits, particle 2 spectates
  unsigned int iemit=0;
  if(pT[0]>pT[1]){ 
    pT_ = pT[0];
    y=yTemp[0];
    if(probTemp[0][0]>UseRandom::rnd()*(probTemp[0][0]+probTemp[0][1])) {
      x1 = x1Solution[0][0];
      x2 = x2Solution[0][0];
    }
    else {
      x1 = x1Solution[0][1];
      x2 = x2Solution[0][1];
    }
  }
  //particle 2 emits, particle 1 spectates
  else {
    iemit=1;
    pT_ = pT[1];
    y=yTemp[1];
    if(probTemp[1][0]>UseRandom::rnd()*(probTemp[1][0]+probTemp[1][1])) {
      x1 = x1Solution[1][0];
      x2 = x2Solution[1][0];
    }
    else {
      x1 = x1Solution[1][1];
      x2 = x2Solution[1][1];
    }
  }
  // find spectator
  unsigned int ispect = iemit == 0 ? 1 : 0;
  // Find the boost from the lab to the c.o.m with the spectator 
  // along the -z axis, and then invert it.
  LorentzRotation eventFrame( ( quark_[0] + quark_[1] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*quark_[ispect];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta() - Constants::pi );
  eventFrame.invert();
  //generation of phi
  double phi = UseRandom::rnd() * Constants::twopi;
  // spectator
  quark_[ispect].setT( 0.5*x2*mHiggs_ );
  quark_[ispect].setX( ZERO );
  quark_[ispect].setY( ZERO );
  quark_[ispect].setZ( -sqrt(0.25*mh2_*x2*x2-mh2_*mu2_) );
  // gluon
  gauge_.setT( pT_*cosh(y)  );
  gauge_.setX( pT_*cos(phi) );
  gauge_.setY( pT_*sin(phi)  );
  gauge_.setZ( pT_*sinh(y)  );
  gauge_.setMass(ZERO);
  // emitter reconstructed from gluon & spectator
  quark_[iemit] = - gauge_ - quark_[ispect];
  quark_[iemit].setT( 0.5*mHiggs_*x1 );
  // boost constructed vectors into the event frame
  quark_[0] = eventFrame * quark_[0];
  quark_[1] = eventFrame * quark_[1];
  gauge_     = eventFrame * gauge_;
  // need to reset masses because for whatever reason the boost  
  // touches the mass component of the five-vector and can make  
  // zero mass objects acquire a floating point negative mass(!).
  gauge_.setMass( ZERO );
  quark_[iemit] .setMass(partons_[iemit ]->mass());
  quark_[ispect].setMass(partons_[ispect]->mass());

  return true;
}

InvEnergy SMHiggsFermionsDecayer::calculateJacobian(double x1, double x2, Energy pT) const{
  double xPerp = abs(2.*pT/mHiggs_);
  Energy jac = mHiggs_*fabs((x1*x2-2.*mu2_*(x1+x2)+sqr(x2)-x2)/xPerp/pow(sqr(x2)-4.*mu2_,1.5));   
  return 1./jac; //jacobian as defined is dptdy=jac*dx1dx2, therefore we have to divide by it
}

bool SMHiggsFermionsDecayer::checkZMomenta(double x1, double x2, double x3, double y, Energy pT) const {
  double xPerp2 = 4.*pT*pT/mHiggs_/mHiggs_;
  static double tolerance = 1e-6; 
  bool isMomentaReconstructed = false;  

  if(pT*sinh(y)>ZERO) {
    if(abs(-sqrt(sqr(x2)-4.*mu2_)+sqrt(sqr(x3)-xPerp2) + sqrt(sqr(x1)-xPerp2 - 4.*mu2_)) <= tolerance ||
       abs(-sqrt(sqr(x2)-4.*mu2_)+sqrt(sqr(x3)-xPerp2)  - sqrt(sqr(x1)-xPerp2 - 4.*mu2_))  <= tolerance) isMomentaReconstructed=true;
  }
  else if(pT*sinh(y) < ZERO){
      if(abs(-sqrt(sqr(x2)-4.*mu2_)-sqrt(sqr(x3)-xPerp2) + sqrt(sqr(x1)-xPerp2 - 4.*mu2_)) <= tolerance ||
         abs(-sqrt(sqr(x2)-4.*mu2_)-sqrt(sqr(x3)-xPerp2)  - sqrt(sqr(x1)-xPerp2 - 4.*mu2_))  <= tolerance) isMomentaReconstructed=true;
  }
  else 
    if(abs(-sqrt(sqr(x2)-4.*mu2_)+ sqrt(sqr(x1)-xPerp2 - 4.*mu2_)) <= tolerance) isMomentaReconstructed=true;
      
  return isMomentaReconstructed;
}
  
