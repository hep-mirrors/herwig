// -*- C++ -*-
//
// SMWDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWDecayer class.
//

#include "SMWDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Shower/Core/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "Herwig/Shower/Core/Base/Branching.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMWDecayer::EPS_=0.00000001;

SMWDecayer::SMWDecayer() 
  : quarkWeight_(6,0.), leptonWeight_(3,0.) {
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
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig StandardModel object in"
				  << "SMWDecayer::doinit()"
				  << Exception::runerror;
  FFWvertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFW());
  // make sure they are initialized
  FFWvertex_->init();
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
      if(!FFWvertex_->allowed(-ix,iy,ParticleID::Wminus))
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
    // if(!FFWvertex_->allowed(-ix,ix+1,ParticleID::Wminus))
    //   throw InitException() << "SMWDecayer::doinit() the W vertex" 
    // 			    << "cannot handle all the lepton modes" 
    // 			    << Exception::abortnow;
    extpart[1] = getParticleData(-ix);
    extpart[2] = getParticleData(ix+1);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    addMode(mode,leptonWeight_[(ix-11)/2],wgt);
  }
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
  os << FFWvertex_ << quarkWeight_ << leptonWeight_ << alpha_;
  
}

void SMWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFWvertex_ >> quarkWeight_ >> leptonWeight_ >> alpha_;
}

ClassDescription<SMWDecayer> SMWDecayer::initSMWDecayer;
// Definition of the static class description member.

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

  static Reference<SMWDecayer,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &SMWDecayer::alpha_, false, false, true, false, false);
}


// return the matrix element squared
double SMWDecayer::me2(const int, const Particle & inpart,
			const ParticleVector & decay,
			MEOption meopt) const {
  if(!ME()) 
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
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
  // compute the matrix element
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      for(unsigned int vhel=0;vhel<3;++vhel) {
	if(iferm>ianti) (*ME())(vhel,ia,ifm)=
	  FFWvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
	else            (*ME())(vhel,ifm,ia)=
	  FFWvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME()->contract(_rho)).real()*UnitRemoval::E2/scale;
  if(abs(decay[0]->id())<=6) output*=3.;
  if(decay[0]->hasColour())      decay[0]->antiColourNeighbour(decay[1]);
  else if(decay[1]->hasColour()) decay[1]->antiColourNeighbour(decay[0]);
  return output;
}

void SMWDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
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
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
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

RealEmissionProcessPtr SMWDecayer::
applyHardMatrixElementCorrection(RealEmissionProcessPtr born) {
  // get the quark and antiquark
  ParticleVector qq;
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix)
    qq.push_back(born->bornOutgoing()[ix]);
  if(!qq[0]->dataPtr()->coloured()) return RealEmissionProcessPtr();
  // ensure quark first
  bool order = qq[0]->id()<0;
  if(order) swap(qq[0],qq[1]);
  // get the momenta
  vector<Lorentz5Momentum> newfs = applyHard(qq);
  // return if no emission
  if(newfs.size()!=3) return RealEmissionProcessPtr();
  // perform final check to ensure energy greater than constituent mass
  for (int i=0; i<2; i++) {
    if (newfs[i].e() < qq[i]->data().constituentMass()) return RealEmissionProcessPtr();
  }
  if (newfs[2].e() < getParticleData(ParticleID::g)->constituentMass())
    return RealEmissionProcessPtr();
  // set masses
  for (int i=0; i<2; i++) newfs[i].setMass(qq[i]->mass());
  newfs[2].setMass(ZERO);
  // decide which particle emits
  bool firstEmits=
    newfs[2].vect().perp2(newfs[0].vect())<
    newfs[2].vect().perp2(newfs[1].vect());
  // create the new quark, antiquark and gluon
  PPtr newg = getParticleData(ParticleID::g)->produceParticle(newfs[2]);
  PPtr newq = qq[0]->dataPtr()->produceParticle(newfs[0]);
  PPtr newa = qq[1]->dataPtr()->produceParticle(newfs[1]);
  // create the output real emission process
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
    born->incoming().push_back(born->bornIncoming()[ix]);
  }
  if(!order) {
    born->outgoing().push_back(newq);
    born->outgoing().push_back(newa);
    born->outgoing().push_back(newg);
  }
  else {
    born->outgoing().push_back(newa);
    born->outgoing().push_back(newq);
    born->outgoing().push_back(newg);
    firstEmits = !firstEmits;
  }
  // make colour connections
  newg->colourNeighbour(newq);
  newa->colourNeighbour(newg);
  if(firstEmits) {
    born->emitter(1);
    born->spectator(2);
  }
  else {
    born->emitter(2);
    born->spectator(1);
  }
  born->emitted(3);
  born->interaction(ShowerInteraction::QCD);
  return born;
}

vector<Lorentz5Momentum> SMWDecayer::
applyHard(const ParticleVector &p) {
  double x, xbar;
  vector<Lorentz5Momentum> fs; 
  // return if no emission
  if (getHard(x, xbar) < UseRandom::rnd() || p.size() != 2) return fs; 
  // centre of mass energy
  Lorentz5Momentum pcm = p[0]->momentum() + p[1]->momentum(); 
  // momenta of quark,antiquark and gluon
  Lorentz5Momentum pq, pa, pg;
  if (p[0]->id() > 0) {
    pq = p[0]->momentum(); 
    pa = p[1]->momentum(); 
  } else {
    pa = p[0]->momentum(); 
    pq = p[1]->momentum(); 
  }
  // boost to boson rest frame
  Boost beta = (pcm.findBoostToCM()); 
  pq.boost(beta);    
  pa.boost(beta);
  // return if fails ?????
  double xg = 2.-x-xbar; 
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return fs;
  Axis u1, u2, u3;
  // moduli of momenta in units of Q and cos theta
  // stick to q direction?
  // p1 is the one that is kept, p2 is the other fermion, p3 the gluon.
  Energy e1, e2, e3; 
  Energy pp1, pp2, pp3;
  bool keepq = true; 
  if (UseRandom::rnd() > sqr(x)/(sqr(x)+sqr(xbar))) 
    keepq = false; 
  if (keepq) {
    pp1 = d_Q_*sqrt(sqr(x)-4.*d_rho_)/2.;
    pp2 = d_Q_*sqrt(sqr(xbar)-4.*d_rho_)/2.;
    e1 = d_Q_*x/2.; 
    e2 = d_Q_*xbar/2.; 
    u1 = pq.vect().unit();
  } else {
    pp2 = d_Q_*sqrt(sqr(x)-4.*d_rho_)/2.;
    pp1 = d_Q_*sqrt(sqr(xbar)-4.*d_rho_)/2.;
    e2 = d_Q_*x/2.; 
    e1 = d_Q_*xbar/2.; 
    u1 = pa.vect().unit();
  }
  pp3 = d_Q_*xg/2.;       
  e3 = pp3; 
  u2 = u1.orthogonal();
  u2 /= u2.mag();
  u3 = u1.cross(u2);
  u3 /= u3.mag();
  double ct2=-2., ct3=-2.;
  if (pp1 == ZERO || pp2 == ZERO || pp3 == ZERO) {
    bool touched = false;
    if (pp1 == ZERO) {
      ct2 = 1; 
      ct3 = -1; 
      touched = true;
    } 
    if (pp2 == ZERO || pp3 == ZERO) {
      ct2 = 1; 
      ct3 = 1; 
      touched = true;
    }
    if (!touched) 
      throw Exception() << "SMWDecayer::applyHard()"
			<< " did not set ct2/3" 
			<< Exception::abortnow;
  } else {
    ct3 = (sqr(pp1)+sqr(pp3)-sqr(pp2))/(2.*pp1*pp3);
    ct2 = (sqr(pp1)+sqr(pp2)-sqr(pp3))/(2.*pp1*pp2);
  }
  double phi = Constants::twopi*UseRandom::rnd();
  double cphi = cos(phi);
  double sphi = sin(phi); 
  double st2 = sqrt(1.-sqr(ct2));
  double st3 = sqrt(1.-sqr(ct3));
  ThreeVector<Energy> pv1, pv2, pv3; 
  pv1 = pp1*u1;
  pv2 = -ct2*pp2*u1 + st2*cphi*pp2*u2 + st2*sphi*pp2*u3;
  pv3 = -ct3*pp3*u1 - st3*cphi*pp3*u2 - st3*sphi*pp3*u3;
  if (keepq) {
    pq = Lorentz5Momentum(pv1, e1);
    pa = Lorentz5Momentum(pv2, e2);
  } else {
    pa = Lorentz5Momentum(pv1, e1);
    pq = Lorentz5Momentum(pv2, e2);
  }
  pg = Lorentz5Momentum(pv3, e3);
  pq.boost(-beta);
  pa.boost(-beta);
  pg.boost(-beta);
  fs.push_back(pq); 
  fs.push_back(pa); 
  fs.push_back(pg); 
  return fs;
}

double SMWDecayer::getHard(double &x1, double &x2) {
  double w = 0.0;
  double y1 = UseRandom::rnd(),y2 = UseRandom::rnd(); 
  // simply double MC efficiency 
  // -> weight has to be divided by two (Jacobian)
  if (y1 + y2 > 1) {
    y1 = 1.-y1; 
    y2 = 1.-y2;
  }
  bool inSoft = false; 
  if (y1 < 0.25) { 
    if (y2 < 0.25) {
      inSoft = true; 
      if (y1 < y2) {
	y1 = 0.25-y1;
	y2 = y1*(1.5 - 2.*y2);
      }	else {
	y2 = 0.25 - y2;
	y1 = y2*(1.5 - 2.*y1);
      }
    } else {
      if (y2 < y1 + 2.*sqr(y1)) return w;
    }
  } else {
    if (y2 < 0.25) {
      if (y1 < y2 + 2.*sqr(y2)) return w;
    }
  } 
  // inside PS?
  x1 = 1.-y1;
  x2 = 1.-y2;
  if(y1*y2*(1.-y1-y2) < d_rho_*sqr(y1+y2)) return w;
  double k1 = getKfromX(x1, x2);
  double k2 = getKfromX(x2, x1);
  // Is it in the quark emission zone?
  if (k1 < d_kt1_) return 0.0;
  // No...is it in the anti-quark emission zone?
  if (k2 < d_kt2_) return 0.0;  
  // Point is in dead zone: compute q qbar g weight
  w = MEV(x1, x2); 
  // for axial: 
  //  w = MEA(x1, x2); 
  // Reweight soft region
  if (inSoft) { 
    if (y1 < y2) w *= 2.*y1;
    else w *= 2.*y2;
  }
  // alpha and colour factors
  Energy2 pt2 = sqr(d_Q_)*(1.-x1)*(1.-x2);
  w *= 1./3./Constants::pi*alpha_->value(pt2); 
  return w; 
}

bool SMWDecayer::
softMatrixElementVeto(ShowerProgenitorPtr initial,ShowerParticlePtr parent,Branching br) {
  // check we should be applying the veto
  if(parent->id()!=initial->progenitor()->id()||
     br.ids[0]!=br.ids[1]||
     br.ids[2]->id()!=ParticleID::g) return false;
  // calculate pt
  double d_z = br.kinematics->z();
  Energy d_qt = br.kinematics->scale();
  Energy2 d_m2 = parent->momentum().m2();
  Energy2 pPerp2 = sqr(d_z*d_qt) - d_m2;
  if(pPerp2<ZERO) {
    parent->vetoEmission(br.type,br.kinematics->scale());
    return true;
  }
  Energy pPerp = (1.-d_z)*sqrt(pPerp2);
  // if not hardest so far don't apply veto
  if(pPerp<initial->highestpT()) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight
  bool veto = !UseRandom::rndbool(weight);
  // if vetoing reset the scale
  if(veto) parent->vetoEmission(br.type,br.kinematics->scale());
  // return the veto
  return veto;
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
