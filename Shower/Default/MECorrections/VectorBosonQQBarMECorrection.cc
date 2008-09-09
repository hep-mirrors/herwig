// -*- C++ -*-
//
// VectorBosonQQBarMECorrection.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorBosonQQBarMECorrection class.
//

#include "VectorBosonQQBarMECorrection.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"

using namespace Herwig;

const double VectorBosonQQBarMECorrection::EPS=0.00000001;

NoPIOClassDescription<VectorBosonQQBarMECorrection> VectorBosonQQBarMECorrection::initVectorBosonQQBarMECorrection;
// Definition of the static class description member.


IBPtr VectorBosonQQBarMECorrection::clone() const {
  return new_ptr(*this);
}

IBPtr VectorBosonQQBarMECorrection::fullclone() const {
  return new_ptr(*this);
}


void VectorBosonQQBarMECorrection::Init() {

  static ClassDocumentation<VectorBosonQQBarMECorrection> documentation
    ("The VectorBosonQQBarMECorrection class implements the matrix"
     " element correction for e+e- -> q qbar and Z/W-> q qbar");

}

bool VectorBosonQQBarMECorrection::canHandle(ShowerTreePtr tree,
					     double & initial, double & final,
					     EvolverPtr evolver) {
  // check radiation on
  if(!evolver->isFSRadiationON()) return false;
  // check 2 outgoing particles
  if(tree->outgoingLines().size()!=2) return false;
  // check 1 or 2 incoming particles
  if(tree->incomingLines().size()!=1&&tree->incomingLines().size()!=2) return false;
  // if one incoming particle must be W/Z or gamma
  int id[2];
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  cit=tree->incomingLines().begin();
  id[0]=cit->first->progenitor()->id();
  if(tree->incomingLines().size()==1) {
    if(abs(id[0])!=ParticleID::Wplus&&id[0]!=ParticleID::Z0&&
       abs(id[0])!=ParticleID::gamma) return false;
  }
  // check incoming leptons beams
  else {
    ++cit;
    id[1]=cit->first->progenitor()->id();
    if(id[0]!=-id[1]||abs(id[0])>15||abs(id[0])<11||abs(id[0])%2!=1) return false;
  }
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  // check outgoing quarks
  cjt=tree->outgoingLines().begin();
  id[0]=cjt->first->progenitor()->id();
  ++cjt;
  id[1]=cjt->first->progenitor()->id();
  if(id[0]!=-id[1]||abs(id[0])>6) return false;
  // otherwise can do it
  initial=1.;
  final  =1.;
  return true;
}

void VectorBosonQQBarMECorrection::
applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  // get the quark and antiquark
  ParticleVector qq; 
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    qq.push_back(cit->first->copy());
  // ensure quark first
  if(qq[0]->id()<0) swap(qq[0],qq[1]);
  // centre of mass energy
  d_Q = (qq[0]->momentum() + qq[1]->momentum()).m();
  // quark mass
  d_m = 0.5*(qq[0]->momentum().m()+qq[1]->momentum().m());
  // set the other parameters
  setRho(sqr(d_m/d_Q));
  setKtildeSymm();
  // get the momenta
  vector<Lorentz5Momentum> newfs = applyHard(qq);
  // return if no emission
  if(newfs.size()!=3) return;
  // perform final check to ensure energy greater than constituent mass
  for (int i=0; i<2; i++) {
    if (newfs[i].e() < qq[i]->data().constituentMass()) return;
  }
  if (newfs[2].e() < getParticleData(ParticleID::g)->constituentMass())
    return;
  // set masses
  for (int i=0; i<2; i++) newfs[i].setMass(qq[i]->mass());
  newfs[2].setMass(0.*MeV);
  // decide which particle emits
  bool firstEmits=
    newfs[2].vect().perp2(newfs[0].vect())<
    newfs[2].vect().perp2(newfs[1].vect());
  // create the new quark, antiquark and gluon
  PPtr newg = getParticleData(ParticleID::g)->produceParticle(newfs[2]);
  PPtr newq,newa;
  if(firstEmits) {
    newq = getParticleData(abs(qq[0]->id()))->produceParticle(newfs[0]);
    newa = new_ptr(Particle(*qq[1]));
    qq[1]->antiColourLine()->removeAntiColoured(newa);
    newa->set5Momentum(newfs[1]);
  }
  else {
    newq = new_ptr(Particle(*qq[0]));
    qq[0]->colourLine()->removeColoured(newq);
    newq->set5Momentum(newfs[0]);
    newa = getParticleData(-abs(qq[0]->id()))->produceParticle(newfs[1]);
  }
  // get the original colour line
  ColinePtr col;
  if(qq[0]->id()>0) col=qq[0]->colourLine();
  else              col=qq[0]->antiColourLine();
  // set the colour lines
  if(firstEmits) {
    col->addColoured(newq);
    col->addAntiColoured(newg);
    newa->colourNeighbour(newg);
  }
  else {
    col->addAntiColoured(newa);
    col->addColoured(newg);
    newq->antiColourNeighbour(newg);
  }
  // change the existing quark and antiquark
  PPtr orig;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(cit->first->progenitor()->id()==newq->id()) {
      // remove old particles from colour line
      col->removeColoured(cit->first->copy());
      col->removeColoured(cit->first->progenitor());
      // insert new particles
      cit->first->copy(newq);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newq,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(!firstEmits);
      if(firstEmits) orig=cit->first->original();
    }
    else {
      // remove old particles from colour line
      col->removeAntiColoured(cit->first->copy());
      col->removeColoured(cit->first->progenitor());
      // insert new particles
      cit->first->copy(newa);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newa,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(firstEmits);
      if(!firstEmits) orig=cit->first->original();
    }
  }
  // add the gluon
  ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
  ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
  gluon->perturbative(false);
  tree->outgoingLines().insert(make_pair(gluon,sg));
  tree->hardMatrixElementCorrection(true);
}

vector<Lorentz5Momentum> VectorBosonQQBarMECorrection::
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
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) return fs;
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
    pp1 = d_Q*sqrt(sqr(x)-4.*d_rho)/2.;
    pp2 = d_Q*sqrt(sqr(xbar)-4.*d_rho)/2.;
    e1 = d_Q*x/2.; 
    e2 = d_Q*xbar/2.; 
    u1 = pq.vect().unit();
  } else {
    pp2 = d_Q*sqrt(sqr(x)-4.*d_rho)/2.;
    pp1 = d_Q*sqrt(sqr(xbar)-4.*d_rho)/2.;
    e2 = d_Q*x/2.; 
    e1 = d_Q*xbar/2.; 
    u1 = pa.vect().unit();
  }
  pp3 = d_Q*xg/2.;       
  e3 = pp3; 
  u2 = u1.orthogonal();
  u2 /= u2.mag();
  u3 = u1.cross(u2);
  u3 /= u3.mag();
  double ct2=-2., ct3=-2.;
  if (pp1 == Energy() || pp2 == Energy() || pp3 == Energy()) {
    bool touched = false;
    if (pp1 == Energy()) {
      ct2 = 1; 
      ct3 = -1; 
      touched = true;
    } 
    if (pp2 == Energy() || pp3 == Energy()) {
      ct2 = 1; 
      ct3 = 1; 
      touched = true;
    }
    if (!touched) 
      throw Exception() << "VectorBosonQQBarMECorrection::applyHard()"
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
  Vector3<Energy> pv1, pv2, pv3; 
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

double VectorBosonQQBarMECorrection::getHard(double &x1, double &x2) {
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
  if(y1*y2*(1.-y1-y2) < d_rho*sqr(y1+y2)) return w;
  double k1 = getKfromX(x1, x2);
  double k2 = getKfromX(x2, x1);
  // Is it in the quark emission zone?
  if (k1 < d_kt1) return 0.0;
  // No...is it in the anti-quark emission zone?
  if (k2 < d_kt2) return 0.0;  
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
  // DGRELL  (this hard wired alpha_S still needs removing)
  w *= 1./3./Constants::pi*0.117997; 
  return w; 
}

bool VectorBosonQQBarMECorrection::
softMatrixElementVeto(ShowerProgenitorPtr initial,ShowerParticlePtr parent,Branching br)
{
  // check we should be applying the veto
  if(parent->id()!=initial->progenitor()->id()||
     br.ids[0]!=br.ids[1]||
     br.ids[2]!=ParticleID::g) return false;
  // calculate pt
  double d_z = br.kinematics->z();
  Energy d_qt = br.kinematics->scale();
  Energy2 d_m2 = parent->momentum().m2();
  Energy pPerp = (1.-d_z)*sqrt( sqr(d_z*d_qt) - d_m2);
  // if not hardest so far don't apply veto
  if(pPerp<initial->highestpT()) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight
  bool veto = !UseRandom::rndbool(weight);
  // if not vetoed reset max
  if(!veto) initial->highestpT(pPerp);
  // if vetoing reset the scale
  if(veto) parent->setEvolutionScale(br.kinematics->scale());
  // return the veto
  return veto;
}


void VectorBosonQQBarMECorrection::setRho(double r) 
{ 
  d_rho = r;
  d_v = sqrt(1.-4.*d_rho);
}

void VectorBosonQQBarMECorrection::setKtildeSymm() { 
  d_kt1 = (1. + sqrt(1. - 4.*d_rho))/2.;
  setKtilde2();
}

void VectorBosonQQBarMECorrection::setKtilde2() { 
   double num = d_rho * d_kt1 + 0.25 * d_v *(1.+d_v)*(1.+d_v);
   double den = d_kt1 - d_rho;
   d_kt2 = num/den;
}

double VectorBosonQQBarMECorrection::getZfromX(double x1, double x2) {
  double uval = u(x2);
  double num = x1 - (2. - x2)*uval;
  double den = sqrt(x2*x2 - 4.*d_rho);
  return uval + num/den;
}

double VectorBosonQQBarMECorrection::getKfromX(double x1, double x2) {
   double zval = getZfromX(x1, x2);
   return (1.-x2)/(zval*(1.-zval));
}

double VectorBosonQQBarMECorrection::MEV(double x1, double x2) {
  // Vector part
  double num = (x1+2.*d_rho)*(x1+2.*d_rho) + (x2+2.*d_rho)*(x2+2.*d_rho) 
    - 8.*d_rho*(1.+2.*d_rho);
  double den = (1.+2.*d_rho)*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho/((1.-x1)*(1.-x1)) 
	  - 2*d_rho/((1.-x2)*(1.-x2)))/d_v;
}

double VectorBosonQQBarMECorrection::MEA(double x1, double x2) {
  // Axial part
  double num = (x1+2.*d_rho)*(x1+2.*d_rho) + (x2+2.*d_rho)*(x2+2.*d_rho) 
    + 2.*d_rho*((5.-x1-x2)*(5.-x1-x2) - 19.0 + 4*d_rho);
  double den = d_v*d_v*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho/((1.-x1)*(1.-x1)) 
	  - 2*d_rho/((1.-x2)*(1.-x2)))/d_v;
}

double VectorBosonQQBarMECorrection::u(double x2) {
  return 0.5*(1. + d_rho/(1.-x2+d_rho));
}

void VectorBosonQQBarMECorrection::
getXXbar(double kti, double z, double &x, double &xbar) {
  x = (1. + sqr(d_v)*(-1. + z) + sqr(kti*(-1. + z))*z*z*z 
       + z*sqrt(sqr(d_v) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z))
       - kti*(-1. + z)*z*(2. + z*(-2 
       + sqrt(sqr(d_v)+ kti*(-1. + z)*z*(2. + kti*(-1. + z)*z))
				  )))
    /(1. - kti*(-1. + z)*z 
      + sqrt(sqr(d_v) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z)));
  xbar = 1. + kti*(-1. + z)*z;
}

double VectorBosonQQBarMECorrection::qWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the quark emission zone?
  if(k1 < d_kt1) {
    rval = MEV(x, xbar)/PS(x, xbar);
    // is it also in the anti-quark emission zone?
    if(k2 < d_kt2) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double VectorBosonQQBarMECorrection::qbarWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) return 0.0;
  double k1 = getKfromX(xbar, x);
  double k2 = getKfromX(x, xbar);
  // Is it in the antiquark emission zone?
  if(k2 < d_kt2) {
    rval = MEV(x, xbar)/PS(xbar, x);
    // is it also in the quark emission zone?
    if(k1 < d_kt1) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double VectorBosonQQBarMECorrection::qWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q), z, x, xb);
  return qWeight(x, xb); 
}

double VectorBosonQQBarMECorrection::qbarWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q), z, x, xb);
  return qbarWeight(x, xb); 
}

double VectorBosonQQBarMECorrection::PS(double x, double xbar) {
  double u = 0.5*(1. + d_rho / (1.-xbar+d_rho));
  double z = u + (x - (2.-xbar)*u)/sqrt(xbar*xbar - 4.*d_rho);
  double brack = (1.+z*z)/(1.-z)- 2.*d_rho/(1-xbar);
  // interesting: the splitting function without the subtraction
  // term. Actually gives a much worse approximation in the collinear
  // limit.  double brack = (1.+z*z)/(1.-z);
  double den = (1.-xbar)*sqrt(xbar*xbar - 4.*d_rho);
  return brack/den;
}
