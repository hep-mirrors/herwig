// -*- C++ -*-
//
// SMWDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWDecayer class.
//

#include "SMWDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMWDecayer::EPS_=0.00000001;

SMWDecayer::SMWDecayer()
  : quarkWeight_(6,0.), leptonWeight_(3,0.), CF_(4./3.),
    NLO_(false) {
  quarkWeight_[0]  = 1.01596;
  quarkWeight_[1]  = 0.0537308;
  quarkWeight_[2]  = 0.0538085;
  quarkWeight_[3]  = 1.01377;
  quarkWeight_[4]  = 1.45763e-05;
  quarkWeight_[5]  = 0.0018143;
  leptonWeight_[0] = 0.356594;
  leptonWeight_[1] = 0.356593;
  leptonWeight_[2] = 0.356333;
  // intermediates
  generateIntermediates(false);
}

void SMWDecayer::doinit() {
  PerturbativeDecayer::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig StandardModel object in"
				  << "SMWDecayer::doinit()"
				  << Exception::runerror;
  FFWVertex_ = hwsm->vertexFFW();
  FFGVertex_ = hwsm->vertexFFG();
  WWWVertex_ = hwsm->vertexWWW();
  FFPVertex_ = hwsm->vertexFFP();
  // make sure they are initialized
  FFGVertex_->init();
  FFWVertex_->init();
  WWWVertex_->init();
  FFPVertex_->init();
  // now set up the decay modes
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  vector<double> wgt(0);
  // W modes
  extpart[0]=getParticleData(ParticleID::Wplus);
  // loop for the quarks
  unsigned int iz=0;
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<6;iy+=2) {
      // check that the combination of particles is allowed
      if(!FFWVertex_->allowed(-ix,iy,ParticleID::Wminus))
	throw InitException() << "SMWDecayer::doinit() the W vertex" 
			      << "cannot handle all the quark modes" 
			      << Exception::abortnow;
      extpart[1] = getParticleData(-ix);
      extpart[2] = getParticleData( iy);
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      addMode(mode,quarkWeight_[iz],wgt);
      ++iz;
    }
  }
  // loop for the leptons
  for(int ix=11;ix<17;ix+=2) {
    // check that the combination of particles is allowed
    // if(!FFWVertex_->allowed(-ix,ix+1,ParticleID::Wminus))
    //   throw InitException() << "SMWDecayer::doinit() the W vertex" 
    // 			    << "cannot handle all the lepton modes" 
    // 			    << Exception::abortnow;
    extpart[1] = getParticleData(-ix);
    extpart[2] = getParticleData(ix+1);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    addMode(mode,leptonWeight_[(ix-11)/2],wgt);
  }
  gluon_ = getParticleData(ParticleID::g);
}

int SMWDecayer::modeNumber(bool & cc,tcPDPtr parent, 
			    const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  int id0=parent->id();
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if(abs(id0)!=ParticleID::Wplus) return imode;
  int idd(0),idu(0);
  if(abs(id1)%2==1&&abs(id2)%2==0) {
    idd=abs(id1);
    idu=abs(id2);
  }
  else if(abs(id1)%2==0&&abs(id2)%2==1) {
    idd=abs(id2);
    idu=abs(id1);
  }
  if(idd==0&&idu==0) {
    return imode;
  }
  else if(idd<=5) {
    imode=idd+idu/2-2;
  }
  else {
    imode=(idd-1)/2+1;
  }
  cc= (id0==ParticleID::Wminus);
  return imode;
}

void SMWDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFWVertex_ << quarkWeight_ << leptonWeight_
     << FFGVertex_ << gluon_ << NLO_
     << WWWVertex_ << FFPVertex_;  
}

void SMWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFWVertex_ >> quarkWeight_ >> leptonWeight_
     >> FFGVertex_ >> gluon_ >> NLO_
     >> WWWVertex_ >> FFPVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMWDecayer,PerturbativeDecayer>
describeHerwigSMWDecayer("Herwig::SMWDecayer", "HwPerturbativeDecay.so");

void SMWDecayer::Init() {

  static ClassDocumentation<SMWDecayer> documentation
    ("The SMWDecayer class is the implementation of the decay"
     " of the W boson to the Standard Model fermions.");

  static ParVector<SMWDecayer,double> interfaceWquarkMax
    ("QuarkMax",
     "The maximum weight for the decay of the W to quarks",
     &SMWDecayer::quarkWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWDecayer,double> interfaceWleptonMax
    ("LeptonMax",
     "The maximum weight for the decay of the W to leptons",
     &SMWDecayer::leptonWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static Switch<SMWDecayer,bool> interfaceNLO
    ("NLO",
     "Whether to return the LO or NLO result",
     &SMWDecayer::NLO_, false, false, false);
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
double SMWDecayer::me2(const int, const Particle & part,
			const ParticleVector & decay,
			MEOption meopt) const {
  if(!ME()) 
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					  incoming,true,false);
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar_,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave_   ,decay[ianti],outgoing);
  // compute the matrix element
  Energy2 scale(sqr(part.mass()));
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      for(unsigned int vhel=0;vhel<3;++vhel) {
	if(iferm>ianti) (*ME())(vhel,ia,ifm)=
	  FFWVertex_->evaluate(scale,wave_[ia],wavebar_[ifm],vectors_[vhel]);
	else            (*ME())(vhel,ifm,ia)=
	  FFWVertex_->evaluate(scale,wave_[ia],wavebar_[ifm],vectors_[vhel]);
      }
    }
  }
  double output=(ME()->contract(rho_)).real()*UnitRemoval::E2/scale;
  if(abs(decay[0]->id())<=6) output*=3.;
  if(decay[0]->hasColour())      decay[0]->antiColourNeighbour(decay[1]);
  else if(decay[1]->hasColour()) decay[1]->antiColourNeighbour(decay[0]);
  // leading-order result
  if(!NLO_) return output;
  // check decay products coloured, otherwise return
  if(!decay[0]->dataPtr()->coloured()) return output;
  // inital masses, couplings  etc
  // W mass
  mW_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mW_));
  // reduced mass
  double mu1  = (decay[0]->dataPtr()->mass())/mW_;
  double mu2  = (decay[1]->dataPtr()->mass())/mW_;
  // scale
  scale_ = sqr(mW_);
  // now for the nlo loop correction
  double virt = CF_*aS_/Constants::pi;
  // now for the real correction
  double realFact=0.;
  for(int iemit=0;iemit<2;++iemit) {
    double phi  = UseRandom::rnd()*Constants::twopi;
    // set the emitter and the spectator
    double muj  = iemit==0 ? mu1 : mu2;
    double muk  = iemit==0 ? mu2 : mu1;
    double muj2 = sqr(muj);
    double muk2 = sqr(muk);
    // calculate y
    double yminus = 0.; 
    double yplus  = 1.-2.*muk*(1.-muk)/(1.-muj2-muk2);
    double y = yminus + UseRandom::rnd()*(yplus-yminus);
    double v = sqrt(sqr(2.*muk2 + (1.-muj2-muk2)*(1.-y))-4.*muk2)
      /(1.-muj2-muk2)/(1.-y);
    double zplus  = (1.+v)*(1.-muj2-muk2)*y/2./(muj2+(1.-muj2-muk2)*y);
    double zminus = (1.-v)*(1.-muj2-muk2)*y/2./(muj2+(1.-muj2-muk2)*y);
    double z = zminus + UseRandom::rnd()*(zplus-zminus);
    double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
    // calculate x1,x2,x3,xT
    double x2 = 1.-y*(1.-muj2-muk2)-muj2+muk2;
    double x1 = 1.+muj2-muk2-z*(x2-2.*muk2);
    // copy the particle objects over for calculateRealEmission
    vector<PPtr> hardProcess(3);
    hardProcess[0] = const_ptr_cast<PPtr>(&part);
    hardProcess[1] = decay[0];
    hardProcess[2] = decay[1];
    realFact += 0.25*jac*sqr(1.-muj2-muk2)/
      sqrt((1.-sqr(muj-muk))*(1.-sqr(muj+muk)))/Constants::twopi
      *2.*CF_*aS_*calculateRealEmission(x1, x2, hardProcess, phi, 
					muj, muk, iemit, true);
  }
  // the born + virtual + real
  output *= (1. + virt + realFact);
  return output;
}

void SMWDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(ix<6) quarkWeight_ [ix]=mode(ix)->maxWeight();
      else     leptonWeight_[ix-6]=mode(ix)->maxWeight();
    }
  }
}

void SMWDecayer::dataBaseOutput(ofstream & output,
				 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<quarkWeight_.size();++ix) {
    output << "newdef " << name() << ":QuarkMax " << ix << " "
	   << quarkWeight_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<leptonWeight_.size();++ix) {
    output << "newdef " << name() << ":LeptonMax " << ix << " "
	   << leptonWeight_[ix] << "\n";
  }
  // parameters for the PerturbativeDecayer base class
  PerturbativeDecayer::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}


void SMWDecayer::
initializeMECorrection(RealEmissionProcessPtr born, double & initial,
		       double & final) {
  // get the quark and antiquark
  ParticleVector qq; 
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix)
    qq.push_back(born->bornOutgoing()[ix]);
  // ensure quark first
  if(qq[0]->id()<0) swap(qq[0],qq[1]);
  // centre of mass energy
  d_Q_ = (qq[0]->momentum() + qq[1]->momentum()).m();
  // quark mass
  d_m_ = 0.5*(qq[0]->momentum().m()+qq[1]->momentum().m());
  // set the other parameters
  setRho(sqr(d_m_/d_Q_));
  setKtildeSymm();
  // otherwise can do it
  initial=1.;
  final  =1.;
}

bool SMWDecayer::softMatrixElementVeto(PPtr parent,
				       PPtr progenitor,
				       const bool & ,
				       const Energy & highestpT,
				       const vector<tcPDPtr> & ids,
				       const double & d_z,
				       const Energy & d_qt,
				       const Energy & ) {
  // check we should be applying the veto
  if(parent->id()!=progenitor->id()||
     ids[0]!=ids[1]||
     ids[2]->id()!=ParticleID::g) return false;
  // calculate pt
  Energy2 d_m2 = parent->momentum().m2();
  Energy2 pPerp2 = sqr(d_z*d_qt) - d_m2;
  if(pPerp2<ZERO) return true;
  Energy pPerp = (1.-d_z)*sqrt(pPerp2);
  // if not hardest so far don't apply veto
  if(pPerp<highestpT) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight and return
  return !UseRandom::rndbool(weight);
}


void SMWDecayer::setRho(double r) 
{ 
  d_rho_ = r;
  d_v_ = sqrt(1.-4.*d_rho_);
}

void SMWDecayer::setKtildeSymm() { 
  d_kt1_ = (1. + sqrt(1. - 4.*d_rho_))/2.;
  setKtilde2();
}

void SMWDecayer::setKtilde2() { 
   double num = d_rho_ * d_kt1_ + 0.25 * d_v_ *(1.+d_v_)*(1.+d_v_);
   double den = d_kt1_ - d_rho_;
   d_kt2_ = num/den;
}

double SMWDecayer::getZfromX(double x1, double x2) {
  double uval = u(x2);
  double num = x1 - (2. - x2)*uval;
  double den = sqrt(x2*x2 - 4.*d_rho_);
  return uval + num/den;
}

double SMWDecayer::getKfromX(double x1, double x2) {
   double zval = getZfromX(x1, x2);
   return (1.-x2)/(zval*(1.-zval));
}

double SMWDecayer::MEV(double x1, double x2) {
  // Vector part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    - 8.*d_rho_*(1.+2.*d_rho_);
  double den = (1.+2.*d_rho_)*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMWDecayer::MEA(double x1, double x2) {
  // Axial part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    + 2.*d_rho_*((5.-x1-x2)*(5.-x1-x2) - 19.0 + 4*d_rho_);
  double den = d_v_*d_v_*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMWDecayer::u(double x2) {
  return 0.5*(1. + d_rho_/(1.-x2+d_rho_));
}

void SMWDecayer::
getXXbar(double kti, double z, double &x, double &xbar) {
  double w = sqr(d_v_) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z);
  if (w < 0) {
    x = -1.; 
    xbar = -1;
  } else {
    x = (1. + sqr(d_v_)*(-1. + z) + sqr(kti*(-1. + z))*z*z*z 
	 + z*sqrt(w)
	 - kti*(-1. + z)*z*(2. + z*(-2 + sqrt(w))))/
      (1. - kti*(-1. + z)*z + sqrt(w));
    xbar = 1. + kti*(-1. + z)*z;
  }
}

double SMWDecayer::qWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the quark emission zone?
  if(k1 < d_kt1_) {
    rval = MEV(x, xbar)/PS(x, xbar);
    // is it also in the anti-quark emission zone?
    if(k2 < d_kt2_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMWDecayer::qbarWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the antiquark emission zone?
  if(k2 < d_kt2_) {
    rval = MEV(x, xbar)/PS(xbar, x);
    // is it also in the quark emission zone?
    if(k1 < d_kt1_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMWDecayer::qWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, x, xb);
  // if exceptionally out of phase space, leave this emission, as there 
  // is no good interpretation for the soft ME correction. 
  if (x < 0 || xb < 0) return 1.0; 
  return qWeight(x, xb); 
}

double SMWDecayer::qbarWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, xb, x);
  // see above in qWeightX. 
  if (x < 0 || xb < 0) return 1.0; 
  return qbarWeight(x, xb); 
}

double SMWDecayer::PS(double x, double xbar) {
  double u = 0.5*(1. + d_rho_ / (1.-xbar+d_rho_));
  double z = u + (x - (2.-xbar)*u)/sqrt(xbar*xbar - 4.*d_rho_);
  double brack = (1.+z*z)/(1.-z)- 2.*d_rho_/(1-xbar);
  // interesting: the splitting function without the subtraction
  // term. Actually gives a much worse approximation in the collinear
  // limit.  double brack = (1.+z*z)/(1.-z);
  double den = (1.-xbar)*sqrt(xbar*xbar - 4.*d_rho_);
  return brack/den;
}

double SMWDecayer::matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
				      const ParticleVector & decay3, MEOption,
				      ShowerInteraction inter) {
  // extract partons and LO momentas
  vector<cPDPtr> partons(1,inpart.dataPtr());
  vector<Lorentz5Momentum> lomom(1,inpart.momentum());
  for(unsigned int ix=0;ix<2;++ix) {
    partons.push_back(decay2[ix]->dataPtr());
    lomom.push_back(decay2[ix]->momentum());
  }
  vector<Lorentz5Momentum> realmom(1,inpart.momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==2) partons.push_back(decay3[ix]->dataPtr());
    realmom.push_back(decay3[ix]->momentum());
  }
  if(partons[0]->id()<0) {
    swap(partons[1],partons[2]);
    swap(lomom[1],lomom[2]);
    swap(realmom[1],realmom[2]);
  }
  scale_ = sqr(inpart.mass());
  double     lome = loME(partons,lomom);
  InvEnergy2 reme = realME(partons,realmom,inter);
  double ratio = reme/lome*sqr(inpart.mass())*4.*Constants::pi;
  if(inter==ShowerInteraction::QCD) ratio *= CF_;
  return ratio;
}

double SMWDecayer::meRatio(vector<cPDPtr> partons, 
			   vector<Lorentz5Momentum> momenta,
			   unsigned int iemitter, bool subtract) const {
  Lorentz5Momentum q = momenta[1]+momenta[2]+momenta[3];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[1].mass()+momenta[2].mass()))*
			(Q2-sqr(momenta[1].mass()-momenta[2].mass())));
  InvEnergy2 D[2];
  double lome(0.);
  for(unsigned int iemit=0;iemit<2;++iemit) {
    unsigned int ispect = iemit==0 ? 1 : 0;    
    Energy2 pipj = momenta[3      ] * momenta[1+iemit ];
    Energy2 pipk = momenta[3      ] * momenta[1+ispect];
    Energy2 pjpk = momenta[1+iemit] * momenta[1+ispect];
    double y = pipj/(pipj+pipk+pjpk); 
    double z = pipk/(     pipk+pjpk);
    Energy mij = sqrt(2.*pipj+sqr(momenta[1+iemit].mass()));
    Energy2 lamB = sqrt((Q2-sqr(mij+momenta[1+ispect].mass()))*
			(Q2-sqr(mij-momenta[1+ispect].mass())));
    Energy2 Qpk = q*momenta[1+ispect];
    Lorentz5Momentum pkt = 
      lambda/lamB*(momenta[1+ispect]-Qpk/Q2*q)
      +0.5/Q2*(Q2+sqr(momenta[1+ispect].mass())-sqr(momenta[1+ispect].mass()))*q;
    Lorentz5Momentum pijt = 
      q-pkt;
    double muj = momenta[1+iemit ].mass()/sqrt(Q2);
    double muk = momenta[1+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			 -vt/v*(2.-z+sqr(momenta[1+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(3);
    lomom[0] = momenta[0];
    if(iemit==0) {
      lomom[1] = pijt;
      lomom[2] = pkt ;
    }
    else {
      lomom[2] = pijt;
      lomom[1] = pkt ;
    }
    if(iemit==0) lome  = loME(partons,lomom);
  }
  InvEnergy2 ratio = realME(partons,momenta,ShowerInteraction::QCD)/lome*abs(D[iemitter])
    /(abs(D[0])+abs(D[1]));
  if(subtract)
    return Q2*(ratio-2.*D[iemitter]);
  else
    return Q2*ratio;
}

double SMWDecayer::loME(const vector<cPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>    vin;
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  VectorWaveFunction    win  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    win.reset(ix);
    vin.push_back(win);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total(0.);
  for(unsigned int inhel=0;inhel<3;++inhel) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	Complex diag1 = FFWVertex_->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	total += norm(diag1);
      }
    }
  }
  // return the answer
  return total;
}
 
InvEnergy2 SMWDecayer::realME(const vector<cPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta,
			      ShowerInteraction inter) const {
  // compute the spinors
  vector<VectorWaveFunction>     vin;
  vector<SpinorWaveFunction>     aout;
  vector<SpinorBarWaveFunction>  fout;
  vector<VectorWaveFunction>     gout;
  VectorWaveFunction    win  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  VectorWaveFunction    gluon(momenta[3],partons[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    gluon.reset(2*ix);
    gout.push_back(gluon);
  }
  for(unsigned int ix=0;ix<3;++ix){
    win.reset(ix);
    vin.push_back(win);
  }
  vector<Complex> diag(3,0.);

  double total(0.);

  AbstractFFVVertexPtr vertex = inter==ShowerInteraction::QCD ? FFGVertex_ : FFPVertex_;
  
  for(unsigned int inhel1=0;inhel1<3;++inhel1) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	  SpinorBarWaveFunction off1 =
	    vertex->evaluate(scale_,3,partons[1]->CC(),fout[outhel1],gout[outhel3]);
	  diag[0] = FFWVertex_->evaluate(scale_,aout[outhel2],off1,vin[inhel1]);
	  
	  SpinorWaveFunction off2 = 
	    vertex->evaluate(scale_,3,partons[2]->CC(),aout[outhel2],gout[outhel3]);
	  diag[1] = FFWVertex_->evaluate(scale_,off2,fout[outhel1],vin[inhel1]);

	  if(inter==ShowerInteraction::QED) {
	    VectorWaveFunction off3 =
	      WWWVertex_->evaluate(scale_,3,partons[0],vin[inhel1],gout[outhel3]);
	    diag[2] = FFWVertex_->evaluate(scale_,aout[outhel2],fout[outhel1],off3);
	  }
	  
	  // sum of diagrams
	  Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  // me2
	  total += norm(sum);
	}
      }
    }
  }
  // divide out the coupling
  total /= norm(vertex->norm());
  // double g = sqrt(2.)*abs(FFWVertex_->norm());
  // double xg = 2.*momenta[3].t()/momenta[0].mass();
  // double xe,mue2;
  // if(abs(partons[1]->id())==ParticleID::eminus) {
  //   xe = 2.*momenta[1].t()/momenta[0].mass();
  //   mue2 = sqr(momenta[1].mass()/momenta[0].mass());
  // }
  // else {
  //   xe = 2.*momenta[2].t()/momenta[0].mass();
  //   mue2 = sqr(momenta[2].mass()/momenta[0].mass());
  // }
  // double cg = -4. * g * g * (-pow(mue2, 3.) / 2. + (xg * xg / 4. + (xe / 2. + 1.) * xg + 5. / 2. * xe - 2.) * mue2 * mue2
  // 			     + (pow(xg, 3.) / 4. + (xe / 4. - 5. / 4.) * xg * xg + (-7. / 2. * xe + 3.) * xg - 3. * xe * xe
  // 				+ 11. / 2. * xe - 7. / 2.) * mue2 + (xg * xg / 2. + (xe - 2.) * xg + xe * xe - 2. * xe + 2.) * (-1. + xg + xe)) * (xe - mue2 - 1.) *
  //   pow(xg, -2.) * pow(-1. + xg + xe - mue2, -2.);
  
  // cerr << "real " << cg/total << "\n";
  // return the total
  return total*UnitRemoval::InvE2;
}

double SMWDecayer::calculateRealEmission(double x1, double x2, 
					 vector<PPtr> hardProcess,
					 double phi, double muj,
					 double muk, int iemit, 
					 bool subtract) const {
  // make partons data object for meRatio
  vector<cPDPtr> partons (3);
  for(int ix=0; ix<3; ++ix)
    partons[ix] = hardProcess[ix]->dataPtr();
  partons.push_back(gluon_);
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3)-0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1)-4.*sqr(muk)+4.*sqr(muj))
		       /(sqr(x2)-4.*sqr(muk))));
  // calculate the momenta
  Energy M = mW_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*sqr(muk),0.)),
			  0.5*M*x2,M*muk); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*sqr(muj),0.)),
			  0.5*M*x1,M*muj);
  Lorentz5Momentum pgluon(0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // boost and rotate momenta
  LorentzRotation eventFrame( ( hardProcess[1]->momentum() +
				hardProcess[2]->momentum() ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*hardProcess[iemit+1]->momentum();
  eventFrame.rotateZ( -spectator.phi()    );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  vector<Lorentz5Momentum> momenta(3);
  momenta[0]   = hardProcess[0]->momentum();
  if(iemit==0) {
    momenta[2] = eventFrame*pspect;
    momenta[1] = eventFrame*pemit ;
  }
  else {
    momenta[1] = eventFrame*pspect;
    momenta[2] = eventFrame*pemit ;
  }
  momenta.push_back(eventFrame*pgluon);
  // calculate the weight
  double realwgt(0.);
  if(1.-x1>1e-5 && 1.-x2>1e-5) 
    realwgt = meRatio(partons,momenta,iemit,subtract);
  return realwgt;
}
