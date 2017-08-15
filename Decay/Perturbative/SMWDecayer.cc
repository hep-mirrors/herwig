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
#include "Herwig/Shower/Core/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "Herwig/Shower/Core/Base/Branching.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMWDecayer::EPS_=0.00000001;

SMWDecayer::SMWDecayer()
  : quarkWeight_(6,0.), leptonWeight_(3,0.), CF_(4./3.), pTmin_(1.*GeV), NLO_(false) {
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
  // make sure they are initialized
  FFGVertex_->init();
  FFWVertex_->init();
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
  os << FFWVertex_ << quarkWeight_ << leptonWeight_ << alpha_
  << FFGVertex_ << gluon_ << ounit( pTmin_, GeV ) << NLO_;  
}

void SMWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFWVertex_ >> quarkWeight_ >> leptonWeight_ >> alpha_
  >> FFGVertex_ >> gluon_ >> iunit( pTmin_, GeV ) >> NLO_;
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

  static Reference<SMWDecayer,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &SMWDecayer::alpha_, false, false, true, false, false);

  static Parameter<SMWDecayer, Energy> interfacePtMin
    ("minpT",
     "The pt cut on hardest emision generation",
     &SMWDecayer::pTmin_, GeV, 1.*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);

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
    realFact = 0.25*jac*sqr(1.-muj2-muk2)/
      sqrt((1.-sqr(muj-muk))*(1.-sqr(muj+muk)))/Constants::twopi
      *2.*CF_*aS_*calculateRealEmission(x1, x2, hardProcess, phi, 
  					muj, muk, iemit, true);
  }
  // the born + virtual + real
  output = output*(1. + virt + realFact);
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

RealEmissionProcessPtr SMWDecayer::
generateHardest(RealEmissionProcessPtr born) {
  assert(born->bornOutgoing().size()==2);
  // check coloured
  if(!born->bornOutgoing()[0]->dataPtr()->coloured()) return RealEmissionProcessPtr();
  // extract required info
  partons_.resize(2);
  quark_.resize(2);
  vector<PPtr> hardProcess;
  wboson_ = born->bornIncoming()[0];
  hardProcess.push_back(wboson_);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    partons_[ix] = born->bornOutgoing()[ix]->dataPtr();
    quark_[ix]   = born->bornOutgoing()[ix]->momentum();
    quark_[ix].setMass(partons_[ix]->mass());
    hardProcess.push_back(born->bornOutgoing()[ix]);
  }
  bool order = partons_[0]->id()<0;
  if(order) {
    swap(partons_[0]   ,partons_[1]   );
    swap(quark_[0]     ,quark_[1]     );
    swap(hardProcess[1],hardProcess[2]);
  }
  gauge_.setMass(0.*MeV);
  // Get the W boson mass.
  mw2_ = (quark_[0] + quark_[1]).m2();
  // Generate emission and set _quark[0,1] and _gauge to be the 
  // momenta of q, qbar and g after the hardest emission:
  if(!getEvent(hardProcess)) {
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
  // // assign the emitter based on evolution scales
  unsigned int iemitter   = quark_[0]*gauge_ > quark_[1]*gauge_ ? 2 : 1;
  unsigned int ispectator = iemitter==1                         ? 1 : 2;
  // create new partices and insert
  PPtr wboson = wboson_->dataPtr()->produceParticle(wboson_->momentum());
  born->incoming().push_back(wboson);
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

double SMWDecayer::meRatio(vector<cPDPtr> partons, 
					 vector<Lorentz5Momentum> momenta,
				 	 unsigned int iemitter, bool subtract) const {
  Lorentz5Momentum q = momenta[1]+momenta[2]+momenta[3];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[1].mass()+momenta[2].mass()))*
			(Q2-sqr(momenta[1].mass()-momenta[2].mass())));
  InvEnergy2 D[2];
  double lome[2];
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
    lome[iemit]  = loME(partons,lomom);
  }
  InvEnergy2 ratio = realME(partons,momenta)*abs(D[iemitter])
    /(abs(D[0]*lome[0])+abs(D[1]*lome[1]));
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
	Complex diag1 = FFWVertex()->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	total += norm(diag1);
      }
    }
  }
  // return the answer
  return total;
}
 
InvEnergy2 SMWDecayer::realME(const vector<cPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta) const {
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
  vector<Complex> diag(2,0.);

  double total(0.);
  for(unsigned int inhel1=0;inhel1<3;++inhel1) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	  SpinorBarWaveFunction off1 =
	    FFGVertex()->evaluate(scale_,3,partons[1],fout[outhel1],gout[outhel3]);
	  diag[0] = FFWVertex()->evaluate(scale_,aout[outhel2],off1,vin[inhel1]);

	  SpinorWaveFunction off2 = 
	    FFGVertex()->evaluate(scale_,3,partons[2],aout[outhel2],gout[outhel3]);
	  diag[1] = FFWVertex()->evaluate(scale_,off2,fout[outhel1],vin[inhel1]);

	  // sum of diagrams
	  Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  // me2
	  total += norm(sum);
	}
      }
    }
  }

  // divide out the coupling
  total /= norm(FFGVertex()->norm());
  // return the total
  return total*UnitRemoval::InvE2;
}
bool SMWDecayer::getEvent(vector<PPtr> hardProcess) {
  vector<Energy> particleMass;
  for(unsigned int ix=0;ix<hardProcess.size();++ix) {
    if(abs(hardProcess[ix]->id())==ParticleID::Wplus) {
      mW_ = hardProcess[ix]->mass();
    }
    else {
      particleMass.push_back(hardProcess[ix]->mass());
    }
  }
  if (particleMass.size()!=2)  {
    throw Exception()
      << "Number of outgoing particles is not equal to 2 in "
      << "SMWFermionPOWHEGDecayer::getEvent()" 
      << Exception::runerror;
  }
  // reduced mass
  mu1_ = particleMass[0]/mW_;
  mu2_ = particleMass[1]/mW_;
  // scale
  scale_ = sqr(mW_);
  // max pT
  Energy pTmax = 0.5*sqrt(mw2_);
  if(pTmax<pTmin_) return false;
  // Define over valued y_max & y_min according to the associated pt_min cut.
  double ymax  =  acosh(pTmax/pTmin_);
  double ymin  = -ymax;
  // pt of the emmission
  pT_ = pTmax;
  // prefactor
  double overEst = 4.;
  double prefactor = overEst*alphaS()->overestimateValue()*CF_*
    (ymax-ymin)/Constants::twopi;
  // loop to generate the pt and rapidity
  bool reject;  
  //arrays to hold the temporary  probabilities whilst the for loop progresses
  double probTemp[2][2]={{0.,0.},{0.,0.}};
  probTemp[0][0]=probTemp[0][1]=probTemp[1][0]=probTemp[1][1]=0.;
  double x1Solution[2][2] = {{0.,0.},{0.,0.}};
  double x2Solution[2][2] = {{0.,0.},{0.,0.}};
  double x3Solution[2]    = {0.,0.};
  Energy pT[2]            = {pTmax,pTmax};
  double yTemp[2]         = {0.,0.};
  double phi              = 0.;
  // do the competition
  for(int i=0; i<2; i++) {
    // set the emitter and the spectator
    double muj  = i==0 ? mu1_ : mu2_;
    double muk  = i==0 ? mu2_ : mu1_;
    double muj2 = sqr(muj);
    double muk2 = sqr(muk);
    do {
      // generation of phi
      phi = UseRandom::rnd() * Constants::twopi;
      // reject the emission
      reject = true; 
      // generate pt
      pT[i] *= pow(UseRandom::rnd(),1./prefactor);
      if(pT[i]<pTmin_) {
        pT[i] = -GeV;
        break;
      }
      // generate xT
      double xT2 = sqr(2./mW_*pT[i]);
      // generate y
      yTemp[i] = ymin + UseRandom::rnd()*(ymax-ymin);
      // generate x3 & x1 from pT & y
      double x1Plus  = 1-muk2+muj2;
      double x1Minus = 2.*muj;
      x3Solution[i]  = 2.*pT[i]*cosh(yTemp[i])/mW_;
      // prefactor
      double weightPrefactor = 0.5/sqrt((1.-sqr(muj-muk))*(1.-sqr(muj+muk)))/overEst;
      // calculate x1 & x2 solutions
      double discrim2 = (-sqr(x3Solution[i])+xT2)*
	(xT2*muk2+2.*x3Solution[i]-sqr(muj2)+2.*muk2+2.*muj2-sqr(x3Solution[i])-1.
	 +2.*muj2*muk2-sqr(muk2)-2.*muk2*x3Solution[i]-2.*muj2*x3Solution[i]);
      // check discrim2 is > 0
      if( discrim2 < ZERO) continue;
      double fact1 =2.*sqr(x3Solution[i])-4.*muk2-6.*x3Solution[i]+4.*muj2-xT2*x3Solution[i]
	+2.*xT2-2.*muj2*x3Solution[i]+2.*muk2*x3Solution[i]+4.;
      double fact2 = (4.-4.*x3Solution[i]+xT2);
      double discriminant = sqrt(discrim2);
      // two solns for x1
      x1Solution[i][0] = (fact1 + 2.*discriminant)/fact2;
      x1Solution[i][1] = (fact1 - 2.*discriminant)/fact2;
      bool found = false;
      for(unsigned int j=0;j<2;++j) {
	// calculate x2
	x2Solution[i][j] = 2.-x3Solution[i]-x1Solution[i][j];
	// set limits on x2
	double root = max(0.,sqr(x1Solution[i][j])-4.*muj2);
	root = sqrt(root);
	double x2Plus  = 1.+muk2-muj2
	  -0.5*(1.-x1Solution[i][j]+muj2-muk2)/(1.-x1Solution[i][j]+muj2)
	  *(x1Solution[i][j]-2.*muj2-root);
	double x2Minus = 1.+muk2-muj2
	  -0.5*(1.-x1Solution[i][j]+muj2-muk2)/(1.-x1Solution[i][j]+muj2)
	  *(x1Solution[i][j]-2.*muj2+root);

        if(x1Solution[i][j]>=x1Minus && x1Solution[i][j]<=x1Plus &&
	   x2Solution[i][j]>=x2Minus && x2Solution[i][j]<=x2Plus &&
           checkZMomenta(x1Solution[i][j], x2Solution[i][j], x3Solution[i], yTemp[i], pT[i], 
			 muj, muk)) {
          probTemp[i][j] = weightPrefactor*pT[i]*
            calculateJacobian(x1Solution[i][j], x2Solution[i][j], pT[i], muj, muk)*
	    calculateRealEmission(x1Solution[i][j], x2Solution[i][j], 
				  hardProcess, phi, muj, muk, i, false);
          found = true;
        }
        else {
          probTemp[i][j] = 0.;
        }
      }
      if(!found) continue;
      // alpha S piece
      double wgt = (probTemp[i][0]+probTemp[i][1])*alphaS()->ratio(sqr(pT[i]));
      // matrix element weight
      reject = UseRandom::rnd()>wgt;
    }
    while(reject);
  } // end of emitter for loop
  // no emission
  if(pT[0]<ZERO&&pT[1]<ZERO) return false;
  //pick the spectator and x1 x2 values
  double x1,x2,y;
  // particle 1 emits, particle 2 spectates
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
  // particle 2 emits, particle 1 spectates
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
  double muk = iemit == 0 ? mu2_ : mu1_;
  double muk2 = sqr(muk);
  double muj = iemit == 0 ? mu1_ : mu2_;
  double muj2 = sqr(muj);
  double xT2 = sqr(2./mW_*pT_);
  // Find the boost from the lab to the c.o.m with the spectator 
  // along the -z axis, and then invert it.
  LorentzRotation eventFrame( ( quark_[0] + quark_[1] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*quark_[ispect];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta() - Constants::pi );
  eventFrame.invert();
  // spectator
  quark_[ispect].setT( 0.5*x2*mW_ );
  quark_[ispect].setX( ZERO );
  quark_[ispect].setY( ZERO );
  quark_[ispect].setZ( -sqrt(0.25*mw2_*x2*x2-mw2_*muk2) );
  // gluon
  gauge_.setT( pT_*cosh(y)  );
  gauge_.setX( pT_*cos(phi) );
  gauge_.setY( pT_*sin(phi) );
  gauge_.setZ( pT_*sinh(y)  );
  gauge_.setMass(ZERO);
  // emitter
  quark_[iemit].setX( -pT_*cos(phi) );
  quark_[iemit].setY( -pT_*sin(phi) );
  quark_[iemit].setZ(  0.5*mW_*sqrt(sqr(x1)-xT2-4.*muj2) );
  if(sqrt(0.25*mw2_*x2*x2-mw2_*muk2)-pT_*sinh(y)<ZERO)
    quark_[iemit].setZ(-quark_[iemit].z());
  quark_[iemit].setT( 0.5*mW_*x1 );
  // boost constructed vectors into the event frame
  quark_[0] = eventFrame * quark_[0];
  quark_[1] = eventFrame * quark_[1];
  gauge_    = eventFrame * gauge_;
  // need to reset masses because for whatever reason the boost  
  // touches the mass component of the five-vector and can make  
  // zero mass objects acquire a floating point negative mass(!).
  gauge_.setMass( ZERO );
  quark_[iemit] .setMass(partons_[iemit ]->mass());
  quark_[ispect].setMass(partons_[ispect]->mass());

  return true;
}

InvEnergy SMWDecayer::calculateJacobian(double x1, double x2, Energy pT, 
						      double muj, double muk) const{
  double xPerp = abs(2.*pT/mW_);
  Energy jac = mW_/xPerp*fabs((x2*sqr(muj)+2.*sqr(muk)*x1
			       +sqr(muk)*x2-x1*x2-sqr(x2)+x2)/pow((sqr(x2)-4.*sqr(muk)),1.5));
  
  return 1./jac; //jacobian as defined is dptdy=jac*dx1dx2, therefore we have to divide by it
}

bool SMWDecayer::checkZMomenta(double x1, double x2, double x3, 
					     double y, Energy pT, double muj,
					     double muk) const {
  double xPerp2 = 4.*pT*pT/mW_/mW_;
  double root1 = sqrt(max(0.,sqr(x2)-4.*sqr(muk)));
  double root2 = sqrt(max(0.,sqr(x1)-xPerp2 - 4.*sqr(muj)));
  static double tolerance = 1e-6; 
  bool isMomentaReconstructed = false;  

  if(pT*sinh(y) > ZERO) {
    if(abs(-root1 + sqrt(sqr(x3)-xPerp2)  + root2) <= tolerance ||
       abs(-root1 + sqrt(sqr(x3)-xPerp2)  - root2)  <= tolerance)
      isMomentaReconstructed=true;
  }
  else if(pT*sinh(y) < ZERO){
    if(abs(-root1 - sqrt(sqr(x3)-xPerp2)  + root2) <= tolerance ||
       abs(-root1 - sqrt(sqr(x3)-xPerp2)  - root2)  <= tolerance)
	isMomentaReconstructed=true;
  }
  else 
    if(abs(-root1+ sqrt(sqr(x1)-xPerp2 - 4.*(muj))) <= tolerance)
      isMomentaReconstructed=true;
      
  return isMomentaReconstructed;
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
  // loop over the possible emitting partons
  double realwgt(0.);

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
  if(1.-x1>1e-5 && 1.-x2>1e-5) 
    realwgt += meRatio(partons,momenta,iemit,subtract);
  
  return realwgt;
}
