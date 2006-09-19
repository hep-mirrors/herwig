// -*- C++ -*-
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

void VectorBosonQQBarMECorrection::Init() {

  static ClassDocumentation<VectorBosonQQBarMECorrection> documentation
    ("The VectorBosonQQBarMECorrection class implements the matrix"
     " element correction for e+e- -> q qbar and Z/W-> q qbar");

}

bool VectorBosonQQBarMECorrection::canHandle(ShowerTreePtr tree,
					     double & initial, double & final) {
   // check 4 external particles
   if(tree->incomingLines().size()!=2||tree->outgoingLines().size()!=2) return false;
   // check incoming leptons beams
   int id[2];
   map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
   cit=tree->incomingLines().begin();
   id[0]=cit->first->progenitor()->id();
   ++cit;
   id[1]=cit->first->progenitor()->id();
   if(id[0]!=-id[1]||abs(id[0])>15||abs(id[0])<11||abs(id[0])%2!=1) return false;
   // check outgoing quarks
   cit=tree->outgoingLines().begin();
   id[0]=cit->first->progenitor()->id();
   ++cit;
   id[1]=cit->first->progenitor()->id();
   if(id[0]!=-id[1]||abs(id[0])>6) return false;
   // otherwise can do it
   initial=1.;
   final  =1.;
   return true;
}

void VectorBosonQQBarMECorrection::
applyHardMatrixElementCorrection(ShowerTreePtr tree)
{
  // get the quark and antiquark
  ParticleVector qq; 
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    qq.push_back(cit->first->copy());
  PPtr temp;
  if(qq[0]->id()<0)
    {
      temp=qq[0];
      qq[0]=qq[1];
      qq[1]=temp;
    }
  // centre of mass energy
  d_Q = (qq[0]->momentum() + qq[1]->momentum()).m();
  // quark mass
  d_m = qq[0]->momentum().m();
  // set the other parameters
  setRho(sqr(d_m/d_Q));
  setKtildeSymm();
  // get the momenta
  vector<Lorentz5Momentum> newfs = applyHard(qq);
  // return if no emission
  if(newfs.size()!=3) return;
  // perform final check to ensure energy greater than constituent mass
  bool check = true; 
  for (int i=0; i<3; i++)
    if (newfs[i].e() < qq[i]->data().constituentMass()) check = false;
  // return if fails
  if (!check) return;
  // set masses
  for (int i=0; i<2; i++) newfs[i].setMass(qq[i]->mass());
  newfs[2].setMass(0.);
  // decide which particle emits
  bool firstEmits=
    newfs[2].vect().perp2(newfs[0].vect())<
    newfs[2].vect().perp2(newfs[1].vect());
  // create the new quark, antiquark and gluon
  PPtr newg = getParticleData(ParticleID::g)->produceParticle(newfs[2]);
  PPtr newq,newa;
  if(firstEmits)
    {
      newq = getParticleData(abs(qq[0]->id()))->produceParticle(newfs[0]);
      newa = new_ptr(Particle(*qq[1]));
      qq[1]->antiColourLine()->removeAntiColoured(newa);
      newa->set5Momentum(newfs[1]);
    }
  else
    {
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
  if(firstEmits)
    {
      col->addColoured(newq);
      col->addAntiColoured(newg);
      newa->colourNeighbour(newg);
    }
  else
    {
      col->addAntiColoured(newa);
      col->addColoured(newg);
      newq->antiColourNeighbour(newg);
    }
  // change the existing quark and antiquark
  PPtr orig;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    {
      if(cit->first->progenitor()->id()==newq->id())
	{
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
      else
	{
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
//   // safe 'largest pt so far'.  
//   Lorentz5Momentum ptot = newfs[0] + newfs[1] + newfs[2];
//   double x = 2*newfs[0]*ptot/sqr(d_Q);
//   double xb = 2*newfs[1]*ptot/sqr(d_Q);
//   Energy qt;
//   double z;
//   qt = d_Q*sqrt(getKfromX(x, xb));
//   z = getZfromX(x, xb);
//   particles[2].pT((1.-z)*sqrt(sqr(z*qt)-sqr(d_m)));
//   qt = d_Q*sqrt(getKfromX(xb, x));
//   z = getZfromX(xb, x);
//   particles[3].pT((1.-z)*sqrt(sqr(z*qt)-sqr(d_m)));
}

vector<Lorentz5Momentum> VectorBosonQQBarMECorrection::
applyHard(const ParticleVector &p) 
{
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
  Vector3 beta = (pcm.findBoostToCM()); 
  pq.boost(beta);    
  pa.boost(beta);
  // return if fails ?????
  double xg = 2.-x-xbar; 
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho*xg*xg) return fs;
  Vector3 u1, u2, u3;
  // moduli of momenta in units of Q and cos theta
  // stick to q direction?
  // p1 is the one that is kept, p2 is the other fermion, p3 the gluon.
  Energy e1, e2, e3; 
  double pp1, pp2, pp3;
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
  if (pp1 == 0 || pp2 == 0 || pp3 == 0) {
    bool touched = false;
    if (pp1 == 0) {
      ct2 = 1; 
      ct3 = -1; 
      touched = true;
    } 
    if (pp2 == 0 || pp3 == 0) {
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
  double phi = 2.*pi*UseRandom::rnd();
  double cphi = cos(phi);
  double sphi = sin(phi); 
  double st2 = sqrt(1.-sqr(ct2));
  double st3 = sqrt(1.-sqr(ct3));
  Vector3 pv1, pv2, pv3; 
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
  // alpha and colour factors (this hard wired alpha_S still needs removing)
  w *= 1./3./pi*0.117997; 
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
  if(pPerp<initial->pT()) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight
  bool veto = !UseRandom::rndbool(weight);
  // if not vetoed reset max
  if(!veto) initial->pT(pPerp);
  // if vetoing reset the scale
  if(veto) parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->scale());
  // return the veto
  return veto;
}
