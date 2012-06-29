// -*- C++ -*-
//
// MEee2gZ2qq.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2gZ2qq class.
//

#include "MEee2gZ2qq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
#include <numeric>

using namespace Herwig;

const double MEee2gZ2qq::EPS_=0.00000001;

void MEee2gZ2qq::doinit() {
  HwMEBase::doinit();
  massOption(vector<unsigned int>(2,massopt_));
  rescalingOption(3);
  if(minflav_>maxflav_)
    throw InitException() << "The minimum flavour " << minflav_  
			  << "must be lower the than maximum flavour " << maxflav_
			  << " in MEee2gZ2qq::doinit() " 
			  << Exception::runerror;
  // set the particle data objects
  Z0_    = getParticleData(ParticleID::Z0);
  gamma_ = getParticleData(ParticleID::gamma);
  gluon_ = getParticleData(ParticleID::g);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    FFZVertex_ = hwsm->vertexFFZ();
    FFPVertex_ = hwsm->vertexFFP();
    FFGVertex_ = hwsm->vertexFFG();
  }
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEee2gZ2qq::doinit() the Herwig++ version must be used" 
			     << Exception::runerror;
}

void MEee2gZ2qq::getDiagrams() const {
  // specific the diagrams
  tcPDPtr ep    = getParticleData(ParticleID::eplus);
  tcPDPtr em    = getParticleData(ParticleID::eminus);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0    = getParticleData(ParticleID::Z0);
  // setup the processes
  for ( int i =minflav_; i<=maxflav_; ++i ) {
    tcPDPtr qk = getParticleData(i);
    tcPDPtr qb = qk->CC();
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, gamma, 3, qk, 3, qb, -1)));
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, Z0   , 3, qk, 3, qb, -2)));
  }
}

Energy2 MEee2gZ2qq::scale() const {
  return sqr(getParticleData(ParticleID::Z0)->mass());
//   return sHat();
}

unsigned int MEee2gZ2qq::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2gZ2qq::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEee2gZ2qq::diagrams(const DiagramVector & diags) const {
  double lastCont(0.5),lastBW(0.5);
  if ( lastXCombPtr() ) {
    lastCont = meInfo()[0];
    lastBW = meInfo()[1];
  }
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(lastCont, i);
    else if ( diags[i]->id() == -2 ) sel.insert(lastBW, i);
  }
  return sel;
}

Selector<const ColourLines *>
MEee2gZ2qq::colourGeometries(tcDiagPtr ) const {
  static const ColourLines c("-5 4");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c);
  return sel;
}

void MEee2gZ2qq::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ << FFPVertex_ << FFGVertex_ 
     << Z0_ << gamma_ << gluon_ << minflav_ 
     << maxflav_ << massopt_ << alpha_ << ounit(pTmin_,GeV)
     << preFactor_;
}

void MEee2gZ2qq::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> FFGVertex_ 
     >> Z0_ >> gamma_ >> gluon_ >> minflav_ 
     >> maxflav_ >> massopt_ >> alpha_ >> iunit(pTmin_,GeV)
     >> preFactor_;
}

ClassDescription<MEee2gZ2qq> MEee2gZ2qq::initMEee2gZ2qq;
// Definition of the static class description member.
void MEee2gZ2qq::Init() {

  static ClassDocumentation<MEee2gZ2qq> documentation
    ("The MEee2gZ2qq class implements the matrix element for e+e- -> q qbar");

  static Parameter<MEee2gZ2qq,int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The PDG code of the quark with the lowest PDG code to produce.",
     &MEee2gZ2qq::minflav_, 1, 1, 6,
     false, false, Interface::limited);

  static Parameter<MEee2gZ2qq,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The PDG code of the quark with the highest PDG code to produce",
     &MEee2gZ2qq::maxflav_, 5, 1, 6,
     false, false, Interface::limited);

  static Switch<MEee2gZ2qq,unsigned int> interfaceTopMassOption
    ("TopMassOption",
     "Option for the treatment of the top quark mass",
     &MEee2gZ2qq::massopt_, 1, false, false);
  static SwitchOption interfaceTopMassOptionOnMassShell
    (interfaceTopMassOption,
     "OnMassShell",
     "The top is produced on its mass shell",
     1);
  static SwitchOption interfaceTopMassOption2
    (interfaceTopMassOption,
     "OffShell",
     "The top is generated off-shell using the mass and width generator.",
     2);

  static Reference<MEee2gZ2qq,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &MEee2gZ2qq::alpha_, false, false, true, false, false);

}

double MEee2gZ2qq::me2() const {
  return loME(mePartonData(),rescaledMomenta(),true); 
}

ProductionMatrixElement MEee2gZ2qq::HelicityME(vector<SpinorWaveFunction>    & fin,
					       vector<SpinorBarWaveFunction> & ain,
					       vector<SpinorBarWaveFunction> & fout,
					       vector<SpinorWaveFunction>    & aout,
					       double & me,
					       double & cont,
					       double & BW ) const {
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming antifermion (vbar spinor)
  // for the outgoing       
  // 0 outgoing fermion     (ubar spinor)
  // 1 outgoing antifermion (v    spinor)
  // me to be returned
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement gamma (PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement Zboson(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  // wavefunctions for the intermediate particles
  VectorWaveFunction interZ,interG;
  // temporary storage of the different diagrams
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double total[3]={0.,0.,0.};
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      interZ = FFZVertex_->evaluate(scale(),1,Z0_,fin[inhel1],ain[inhel2]);
      // intermediate photon
      interG = FFPVertex_->evaluate(scale(),1,gamma_,fin[inhel1],ain[inhel2]);
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {		
	  // first the Z exchange diagram
	  diag1 = FFZVertex_->evaluate(scale(),aout[outhel2],fout[outhel1],
				       interZ);
	  // then the photon exchange diagram
	  diag2 = FFPVertex_->evaluate(scale(),aout[outhel2],fout[outhel1],
				       interG);
	  // add up squares of individual terms
	  total[1] += norm(diag1);
	  Zboson(inhel1,inhel2,outhel1,outhel2) = diag1;
	  total[2] += norm(diag2);
	  gamma (inhel1,inhel2,outhel1,outhel2) = diag2;
	  // the full thing including interference
	  diag1 += diag2;
	  total[0] += norm(diag1);
	  output(inhel1,inhel2,outhel1,outhel2)=diag1;
	}
      }
    }
  }
  for(int ix=0;ix<3;++ix) total[ix] *= 0.25;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = 
      {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
       beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    total[0] = output.average(rho[0],rho[1]);
    total[1] = Zboson.average(rho[0],rho[1]);
    total[2] = gamma .average(rho[0],rho[1]);
  }
  // results
  for(int ix=0;ix<3;++ix) total[ix]*= 3.;
  cont = total[2];
  BW   = total[1];
  me   = total[0];
  return output;
}

void MEee2gZ2qq::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  if(hard[0]->id()<hard[1]->id()) swap(hard[0],hard[1]);
  if(hard[2]->id()<hard[3]->id()) swap(hard[2],hard[3]);
  vector<SpinorWaveFunction>    fin,aout;
  vector<SpinorBarWaveFunction> ain,fout;
  // get wave functions for off-shell momenta for later on
  SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
  SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
  SpinorBarWaveFunction(fout,hard[2],outgoing,true ,true);
  SpinorWaveFunction(   aout,hard[3],outgoing,true ,true);
  // now rescale the momenta and compute the matrix element with the
  // rescaled momenta for correlations
  vector<Lorentz5Momentum> momenta;
  cPDVector data;
  for(unsigned int ix=0;ix<4;++ix) {
    momenta.push_back(hard[ix]->momentum());
    data   .push_back(hard[ix]->dataPtr());
  }
  rescaleMomenta(momenta,data);
  SpinorWaveFunction    ein  (rescaledMomenta()[0],data[0],incoming);
  SpinorBarWaveFunction pin  (rescaledMomenta()[1],data[1],incoming);
  SpinorBarWaveFunction qkout(rescaledMomenta()[2],data[2],outgoing);
  SpinorWaveFunction    qbout(rescaledMomenta()[3],data[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ; fin [ix] = ein  ;
    pin.reset(ix)  ; ain [ix] = pin  ;
    qkout.reset(ix); fout[ix] = qkout;
    qbout.reset(ix); aout[ix] = qbout;
  }
  // calculate the matrix element
  double me,cont,BW;
  ProductionMatrixElement prodme=HelicityME(fin,ain,fout,aout,me,cont,BW);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(prodme);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix) {
    tSpinPtr spin = hard[ix]->spinInfo();
    if(ix<2) {
      tcPolarizedBeamPDPtr beam = 
	dynamic_ptr_cast<tcPolarizedBeamPDPtr>(hard[ix]->dataPtr());
      if(beam) spin->rhoMatrix() = beam->rhoMatrix();
    }
    spin->productionVertex(hardvertex);
  }
}

void MEee2gZ2qq::rebind(const TranslationMap & trans) {
  FFZVertex_ = trans.translate(FFZVertex_);
  FFPVertex_ = trans.translate(FFPVertex_);
  FFGVertex_ = trans.translate(FFGVertex_);
  Z0_        = trans.translate(Z0_);
  gamma_     = trans.translate(gamma_);
  gluon_     = trans.translate(gluon_);
  HwMEBase::rebind(trans);
}

IVector MEee2gZ2qq::getReferences() {
  IVector ret = HwMEBase::getReferences();
  ret.push_back(FFZVertex_);
  ret.push_back(FFPVertex_);
  ret.push_back(FFGVertex_);
  ret.push_back(Z0_       );
  ret.push_back(gamma_    );
  ret.push_back(gluon_    );
  return ret;
}

void MEee2gZ2qq::initializeMECorrection(ShowerTreePtr , double &  initial,
					double & final) {
  d_Q_ = sqrt(sHat());
  d_m_ = 0.5*(meMomenta()[2].mass()+meMomenta()[3].mass());
  // set the other parameters
  setRho(sqr(d_m_/d_Q_));
  setKtildeSymm();
  // otherwise can do it
  initial=1.;
  final  =1.;
}

void MEee2gZ2qq::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  // get the quark and antiquark
  ParticleVector qq; 
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    qq.push_back(cit->first->copy());
  // ensure quark first
  if(qq[0]->id()<0) swap(qq[0],qq[1]);
  // get the momenta
  vector<Lorentz5Momentum> newfs = applyHard(qq);
  // return if no emission
  if(newfs.size()!=3) return;
  // perform final check to ensure energy greater than constituent mass
  for (int i=0; i<2; i++) {
    if (newfs[i].e() < qq[i]->data().constituentMass()) return;
  }
  if (newfs[2].e() < gluon_->constituentMass())
    return;
  // set masses
  for (int i=0; i<2; i++) newfs[i].setMass(qq[i]->mass());
  newfs[2].setMass(ZERO);
  // decide which particle emits
  bool firstEmits=
    newfs[2].vect().perp2(newfs[0].vect())<
    newfs[2].vect().perp2(newfs[1].vect());
  // create the new quark, antiquark and gluon
  PPtr newg = gluon_->produceParticle(newfs[2]);
  PPtr newq,newa;
  if(firstEmits) {
    newq = qq[0]->dataPtr()->produceParticle(newfs[0]);
    newa = new_ptr(Particle(*qq[1]));
    qq[1]->antiColourLine()->removeAntiColoured(newa);
    newa->set5Momentum(newfs[1]);
  }
  else {
    newq = new_ptr(Particle(*qq[0]));
    qq[0]->colourLine()->removeColoured(newq);
    newq->set5Momentum(newfs[0]);
    newa = qq[1]->dataPtr()->produceParticle(newfs[1]);
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

vector<Lorentz5Momentum> MEee2gZ2qq::
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
      throw Exception() << "MEee2gZ2qq::applyHard()"
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

double MEee2gZ2qq::getHard(double &x1, double &x2) {
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

bool MEee2gZ2qq::softMatrixElementVeto(ShowerProgenitorPtr initial,
				       ShowerParticlePtr parent,Branching br) {
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

void MEee2gZ2qq::setRho(double r) { 
  d_rho_ = r;
  d_v_ = sqrt(1.-4.*d_rho_);
}

void MEee2gZ2qq::setKtildeSymm() { 
  d_kt1_ = (1. + sqrt(1. - 4.*d_rho_))/2.;
  setKtilde2();
}

void MEee2gZ2qq::setKtilde2() { 
   double num = d_rho_ * d_kt1_ + 0.25 * d_v_ *(1.+d_v_)*(1.+d_v_);
   double den = d_kt1_ - d_rho_;
   d_kt2_ = num/den;
}

double MEee2gZ2qq::getZfromX(double x1, double x2) {
  double uval = u(x2);
  double num = x1 - (2. - x2)*uval;
  double den = sqrt(x2*x2 - 4.*d_rho_);
  return uval + num/den;
}

double MEee2gZ2qq::getKfromX(double x1, double x2) {
   double zval = getZfromX(x1, x2);
   return (1.-x2)/(zval*(1.-zval));
}

double MEee2gZ2qq::MEV(double x1, double x2) {
  // Vector part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    - 8.*d_rho_*(1.+2.*d_rho_);
  double den = (1.+2.*d_rho_)*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double MEee2gZ2qq::MEA(double x1, double x2) {
  // Axial part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    + 2.*d_rho_*((5.-x1-x2)*(5.-x1-x2) - 19.0 + 4*d_rho_);
  double den = d_v_*d_v_*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double MEee2gZ2qq::u(double x2) {
  return 0.5*(1. + d_rho_/(1.-x2+d_rho_));
}

void MEee2gZ2qq::
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

double MEee2gZ2qq::qWeight(double x, double xbar) {
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

double MEee2gZ2qq::qbarWeight(double x, double xbar) {
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

double MEee2gZ2qq::qWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, x, xb);
  // if exceptionally out of phase space, leave this emission, as there 
  // is no good interpretation for the soft ME correction. 
  if (x < 0 || xb < 0) return 1.0; 
  return qWeight(x, xb); 
}

double MEee2gZ2qq::qbarWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, xb, x);
  // see above in qWeightX. 
  if (x < 0 || xb < 0) return 1.0; 
  return qbarWeight(x, xb); 
}

double MEee2gZ2qq::PS(double x, double xbar) {
  double u = 0.5*(1. + d_rho_ / (1.-xbar+d_rho_));
  double z = u + (x - (2.-xbar)*u)/sqrt(xbar*xbar - 4.*d_rho_);
  double brack = (1.+z*z)/(1.-z)- 2.*d_rho_/(1-xbar);
  // interesting: the splitting function without the subtraction
  // term. Actually gives a much worse approximation in the collinear
  // limit.  double brack = (1.+z*z)/(1.-z);
  double den = (1.-xbar)*sqrt(xbar*xbar - 4.*d_rho_);
  return brack/den;
}

HardTreePtr MEee2gZ2qq::generateHardest(ShowerTreePtr tree) {
  // get the momenta of the incoming and outgoing partons 
  // incoming
  ShowerProgenitorPtr 
    emProgenitor = tree->incomingLines().begin() ->first,
    epProgenitor = tree->incomingLines().rbegin()->first; 
  // outgoing
  ShowerProgenitorPtr 
    qkProgenitor = tree->outgoingLines().begin() ->first,
    qbProgenitor = tree->outgoingLines().rbegin()->first;
  // get the order right 
  if(emProgenitor->id()<0) swap(emProgenitor,epProgenitor);
  if(qkProgenitor->id()<0) swap(qkProgenitor,qbProgenitor);
  loMomenta_.resize(0);
  // extract the momenta 
  loMomenta_.push_back(emProgenitor->progenitor()->momentum());
  loMomenta_.push_back(epProgenitor->progenitor()->momentum());
  loMomenta_.push_back(qkProgenitor->progenitor()->momentum());
  loMomenta_.push_back(qbProgenitor->progenitor()->momentum());
  // and ParticleData objects
  partons_.resize(0);
  partons_.push_back(emProgenitor->progenitor()->dataPtr());
  partons_.push_back(epProgenitor->progenitor()->dataPtr());
  partons_.push_back(qkProgenitor->progenitor()->dataPtr());
  partons_.push_back(qbProgenitor->progenitor()->dataPtr());
  partons_.push_back(gluon_);
  // boost from lab to CMS frame with outgoing particles
  // along the z axis
  LorentzRotation eventFrame( ( loMomenta_[2] + loMomenta_[3] ).findBoostToCM() );
  Lorentz5Momentum spectator = eventFrame*loMomenta_[2];
  eventFrame.rotateZ( -spectator.phi() );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  // mass of the final-state system
  Energy2 M2 = (loMomenta_[2]+loMomenta_[3]).m2();
  Energy  M  = sqrt(M2);
  double mu1 = loMomenta_[2].mass()/M;
  double mu2 = loMomenta_[3].mass()/M;
  double mu12 = sqr(mu1), mu22 = sqr(mu2);
  double lambda = sqrt(1.+sqr(mu12)+sqr(mu22)-2.*mu12-2.*mu22-2.*mu12*mu22);
  // max pT
  Energy pTmax = 0.5*sqrt(M2)*
    (1.-sqr(loMomenta_[2].mass()+loMomenta_[3].mass())/M2);
  // max y
  if ( pTmax < pTmin_ ) 
    return HardTreePtr();
  double ymax = acosh(pTmax/pTmin_);
  // prefactor for the overestimate of the Sudakov
  double a = 4./3.*alpha_->overestimateValue()/Constants::twopi*
    2.*ymax*preFactor_;
  // variables for the emission
  Energy pT[2];
  double  y[2],phi[2],x3[2],x1[2][2],x2[2][2];
  double contrib[2][2];
  // storage of the real emmision momenta
  vector<Lorentz5Momentum> realMomenta[2][2]=
    {{vector<Lorentz5Momentum>(5),vector<Lorentz5Momentum>(5)},
     {vector<Lorentz5Momentum>(5),vector<Lorentz5Momentum>(5)}};
  for(unsigned int ix=0;ix<2;++ix)
    for(unsigned int iy=0;iy<2;++iy)
      for(unsigned int iz=0;iz<2;++iz)
	realMomenta[ix][iy][iz] = loMomenta_[iz];
  // generate the emission
  for(unsigned int ix=0;ix<2;++ix) {
    if(ix==1) {
      swap(mu1 ,mu2 );
      swap(mu12,mu22);
    }
    pT[ix] = pTmax;
    y [ix] = 0.;
    bool reject = true;
    do {
      // generate pT
      pT[ix] *= pow(UseRandom::rnd(),1./a);
      if(pT[ix]<pTmin_) {
        pT[ix] = -GeV;
        break;
      }
      // generate y
      y[ix] = -ymax+2.*UseRandom::rnd()*ymax;
      // generate phi
      phi[ix] = UseRandom::rnd()*Constants::twopi;
      // calculate x3 and check in allowed region
      x3[ix] = 2.*pT[ix]*cosh(y[ix])/M;
      if(x3[ix] < 0. || x3[ix] > 1. -sqr( mu1 + mu2 ) ) continue;
      // find the possible solutions for x1
      double xT2 = sqr(2./M*pT[ix]);
      double root = (-sqr(x3[ix])+xT2)*
	(xT2*mu22+2.*x3[ix]-sqr(mu12)+2.*mu22+2.*mu12-sqr(x3[ix])-1.
	 +2.*mu12*mu22-sqr(mu22)-2.*mu22*x3[ix]-2.*mu12*x3[ix]);
      double c1=2.*sqr(x3[ix])-4.*mu22-6.*x3[ix]+4.*mu12-xT2*x3[ix]
	+2.*xT2-2.*mu12*x3[ix]+2.*mu22*x3[ix]+4.;
      if(root<0.) continue;
      x1[ix][0] = 1./(4.-4.*x3[ix]+xT2)*(c1-2.*sqrt(root));
      x1[ix][1] = 1./(4.-4.*x3[ix]+xT2)*(c1+2.*sqrt(root));
      // change sign of y if 2nd particle emits
      if(ix==1) y[ix] *=-1.;
      // loop over the solutions
      for(unsigned int iy=0;iy<2;++iy) {
	contrib[ix][iy]=0.;
	// check x1 value allowed
	if(x1[ix][iy]<2.*mu1||x1[ix][iy]>1.+mu12-mu22) continue;
	// calculate x2 value and check allowed
	x2[ix][iy] = 2.-x3[ix]-x1[ix][iy];
	double root = max(0.,sqr(x1[ix][iy])-4.*mu12);
	root = sqrt(root);
	double x2min = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12+root);
	double x2max = 1.+mu22-mu12
	  -0.5*(1.-x1[ix][iy]+mu12-mu22)/(1.-x1[ix][iy]+mu12)*(x1[ix][iy]-2.*mu12-root);
	if(x2[ix][iy]<x2min||x2[ix][iy]>x2max) continue;
	// check the z components
	double z1 =  sqrt(sqr(x1[ix][iy])-4.*mu12-xT2);
	double z2 = -sqrt(sqr(x2[ix][iy])-4.*mu22);
	double z3 =  pT[ix]*sinh(y[ix])*2./M;
	if(ix==1) z3 *=-1.;
	if(abs(-z1+z2+z3)<1e-9) z1 *= -1.;
	if(abs(z1+z2+z3)>1e-5) continue;
	// construct the momenta
	realMomenta[ix][iy][4] =
	  Lorentz5Momentum(pT[ix]*cos(phi[ix]),pT[ix]*sin(phi[ix]),
			   pT[ix]*sinh(y[ix]) ,pT[ix]*cosh(y[ix]),ZERO);
	if(ix==0) {
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	  realMomenta[ix][iy][3] =
	    Lorentz5Momentum(ZERO,ZERO, z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	}
	else {
	  realMomenta[ix][iy][2] =
	    Lorentz5Momentum(ZERO,ZERO,-z2*0.5*M,x2[ix][iy]*0.5*M,M*mu2);
	  realMomenta[ix][iy][3] =
	    Lorentz5Momentum(-pT[ix]*cos(phi[ix]),-pT[ix]*sin(phi[ix]),
			     -z1*0.5*M,x1[ix][iy]*0.5*M,M*mu1);
	}
	// boost the momenta back to the lab
	for(unsigned int iz=2;iz<5;++iz)
	  realMomenta[ix][iy][iz] *= eventFrame;
	// jacobian and prefactors for the weight
	Energy J = M/sqrt(xT2)*abs(-x1[ix][iy]*x2[ix][iy]+2.*mu22*x1[ix][iy]
				   +x2[ix][iy]+x2[ix][iy]*mu12+mu22*x2[ix][iy]
				   -sqr(x2[ix][iy]))
	  /pow(sqr(x2[ix][iy])-4.*mu22,1.5);
	// prefactors etc
	contrib[ix][iy] = 0.5*pT[ix]/J/preFactor_/lambda;
	// matrix element piece
	contrib[ix][iy] *= alpha_->ratio(sqr(pT[ix]))*
	  meRatio(partons_,realMomenta[ix][iy],ix,false);
      }
      if(contrib[ix][0]+contrib[ix][1]>1.) {
	ostringstream s;
	s << "MEee2gZ2qq::generateHardest weight for channel " << ix
	  << "is " << contrib[ix][0]+contrib[ix][1] 
	  << " which is greater than 1";
	generator()->logWarning( Exception(s.str(), Exception::warning) );
      }
      reject =  UseRandom::rnd() > contrib[ix][0] + contrib[ix][1];
    }
    while (reject);
    if(pT[ix]<pTmin_)
      pT[ix] = -GeV;
  }
  if(pT[0]<ZERO && pT[1]<ZERO) {
    qkProgenitor->maximumpT(pTmin_);
    qbProgenitor->maximumpT(pTmin_);
    return HardTreePtr();
  }
  // now pick the emission with highest pT
  vector<Lorentz5Momentum> emmision;
  unsigned int iemit=0,ispect=0;
  Energy pTveto;
  if(pT[0]>pT[1]) {
    iemit  = 2;
    ispect = 3;
    pTveto = pT[0];
    if(UseRandom::rnd()<contrib[0][0]/(contrib[0][0]+contrib[0][1]))
      emmision = realMomenta[0][0];
    else
      emmision = realMomenta[0][1];
  }
  else {
    iemit  = 3;
    ispect = 2;
    pTveto = pT[1];
    if(UseRandom::rnd()<contrib[1][0]/(contrib[1][0]+contrib[1][1]))
      emmision = realMomenta[1][0];
    else
      emmision = realMomenta[1][1];
  }
  // Make the particles for the hard tree
  ShowerParticleVector hardParticles;
  for(unsigned int ix=0;ix<partons_.size();++ix) {
    hardParticles.push_back(new_ptr(ShowerParticle(partons_[ix],ix>=2)));
    hardParticles.back()->set5Momentum(emmision[ix]);
  }
  ShowerParticlePtr parent(new_ptr(ShowerParticle(partons_[iemit],true)));
  Lorentz5Momentum parentMomentum(emmision[iemit]+emmision[4]);
  parentMomentum.setMass(partons_[iemit]->mass());
  parent->set5Momentum(parentMomentum);
  // Create the vectors of HardBranchings to create the HardTree:
  vector<HardBranchingPtr> spaceBranchings,allBranchings;
  // Incoming boson:
  for(unsigned int ix=0;ix<2;++ix) {
    spaceBranchings.push_back(new_ptr(HardBranching(hardParticles[ix],SudakovPtr(),
						    HardBranchingPtr(),
						    HardBranching::Incoming)));
    allBranchings.push_back(spaceBranchings.back());
  }
  // Outgoing particles from hard emission:
  HardBranchingPtr spectatorBranch(new_ptr(HardBranching(hardParticles[ispect],
 							 SudakovPtr(),HardBranchingPtr(),
 							 HardBranching::Outgoing)));
  HardBranchingPtr emitterBranch(new_ptr(HardBranching(parent,SudakovPtr(),
						       HardBranchingPtr(),
						       HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(hardParticles[iemit], 
 						SudakovPtr(),HardBranchingPtr(),
 						HardBranching::Outgoing)));
  emitterBranch->addChild(new_ptr(HardBranching(hardParticles[4],
 						SudakovPtr(),HardBranchingPtr(),
 						HardBranching::Outgoing)));
  if(iemit==0) {
    allBranchings.push_back(emitterBranch);
    allBranchings.push_back(spectatorBranch);
  } 
  else {
    allBranchings.push_back( spectatorBranch );
    allBranchings.push_back( emitterBranch );
  }
  // Make the HardTree from the HardBranching vectors.
  HardTreePtr hardtree = new_ptr(HardTree(allBranchings,spaceBranchings,
					   ShowerInteraction::QCD));
  // Set the maximum pt for all other emissions
  qkProgenitor->maximumpT(pTveto);
  qbProgenitor->maximumpT(pTveto);
  // Connect the particles with the branchings in the HardTree
  hardtree->connect( qkProgenitor->progenitor(), allBranchings[2] );
  hardtree->connect( qbProgenitor->progenitor(), allBranchings[3] );
  // colour flow
  ColinePtr newline=new_ptr(ColourLine());
  for(set<HardBranchingPtr>::const_iterator cit=hardtree->branchings().begin();
      cit!=hardtree->branchings().end();++cit) {
    if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3)
      newline->addColoured((**cit).branchingParticle());
    else if((**cit).branchingParticle()->dataPtr()->iColour()==PDT::Colour3bar)
      newline->addAntiColoured((**cit).branchingParticle());
  }
  // Return the HardTree
  return hardtree;
}

double MEee2gZ2qq::meRatio(vector<cPDPtr> partons, 
			   vector<Lorentz5Momentum> momenta,
			   unsigned int iemitter,bool subtract) const {
  Lorentz5Momentum q = momenta[2]+momenta[3]+momenta[4];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[2].mass()+momenta[3].mass()))*
			(Q2-sqr(momenta[2].mass()-momenta[3].mass())));
  InvEnergy2 D[2];
  double lome[2];
  for(unsigned int iemit=0;iemit<2;++iemit) {
    unsigned int ispect = iemit==0 ? 1 : 0;
    Energy2 pipj = momenta[4      ] * momenta[2+iemit ];
    Energy2 pipk = momenta[4      ] * momenta[2+ispect];
    Energy2 pjpk = momenta[2+iemit] * momenta[2+ispect];
    double y = pipj/(pipj+pipk+pjpk);
    double z = pipk/(     pipk+pjpk);
    Energy mij = sqrt(2.*pipj+sqr(momenta[2+iemit].mass()));
    Energy2 lamB = sqrt((Q2-sqr(mij+momenta[2+ispect].mass()))*
			(Q2-sqr(mij-momenta[2+ispect].mass())));
    Energy2 Qpk = q*momenta[2+ispect];
    Lorentz5Momentum pkt = 
      lambda/lamB*(momenta[2+ispect]-Qpk/Q2*q)
      +0.5/Q2*(Q2+sqr(momenta[2+ispect].mass())-sqr(momenta[2+ispect].mass()))*q;
    Lorentz5Momentum pijt = 
      q-pkt;
    double muj = momenta[2+iemit ].mass()/sqrt(Q2);
    double muk = momenta[2+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			 -vt/v*(2.-z+sqr(momenta[2+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(4);
    lomom[0] = momenta[0];
    lomom[1] = momenta[1];
    if(iemit==0) {
      lomom[2] = pijt;
      lomom[3] = pkt ;
    }
    else {
      lomom[3] = pijt;
      lomom[2] = pkt ;
    }
    lome[iemit]  = loME(partons,lomom,false)/3.;
  }
  InvEnergy2 ratio = realME(partons,momenta)
    *abs(D[iemitter])/(abs(D[0]*lome[0])+abs(D[1]*lome[1]));
  if(subtract)
    return Q2*(ratio-2.*D[iemitter]);
  else
    return Q2*ratio;
}

double MEee2gZ2qq::loME(const vector<cPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta,
			bool first) const {
  // compute the spinors
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  SpinorWaveFunction    ein  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction pin  (momenta[1],partons[1],incoming);
  SpinorBarWaveFunction qkout(momenta[2],partons[2],outgoing);
  SpinorWaveFunction    qbout(momenta[3],partons[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ;
    fin.push_back( ein  );
    pin.reset(ix)  ;
    ain.push_back( pin  );
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  // compute the matrix element
  double me,lastCont,lastBW;
  HelicityME(fin,ain,fout,aout,me,lastCont,lastBW);
  // save the components
  if(first) {
    DVector save;
    save.push_back(lastCont);
    save.push_back(lastBW);
    meInfo(save);
  }
  // return the answer
  return me;
}

InvEnergy2 MEee2gZ2qq::realME(const vector<cPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  vector<VectorWaveFunction> gout;
  SpinorWaveFunction    ein  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction pin  (momenta[1],partons[1],incoming);
  SpinorBarWaveFunction qkout(momenta[2],partons[2],outgoing);
  SpinorWaveFunction    qbout(momenta[3],partons[3],outgoing);
  VectorWaveFunction    gluon(momenta[4],partons[4],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ;
    fin.push_back( ein  );
    pin.reset(ix)  ;
    ain.push_back( pin  );
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    gluon.reset(2*ix);
    gout.push_back(gluon);
  }
  vector<Complex> diag(4,0.);
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1);
  double total(0.);
  for(unsigned int inhel1=0;inhel1<2;++inhel1) {
    for(unsigned int inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      VectorWaveFunction interZ = 
	FFZVertex_->evaluate(scale(),1,Z0_,fin[inhel1],ain[inhel2]);
      // intermediate photon
      VectorWaveFunction interG = 
	FFPVertex_->evaluate(scale(),1,gamma_,fin[inhel1],ain[inhel2]);
      for(unsigned int outhel1=0;outhel1<2;++outhel1) {
	for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	  for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	    SpinorBarWaveFunction off1 =
	      FFGVertex_->evaluate(scale(),3,partons[2],fout[outhel1],gout[outhel3]);
	    diag[0] = FFZVertex_->evaluate(scale(),aout[outhel2],off1,interZ);
	    diag[1] = FFPVertex_->evaluate(scale(),aout[outhel2],off1,interG);
	    SpinorWaveFunction    off2 = 
	      FFGVertex_->evaluate(scale(),3,partons[3],aout[outhel2],gout[outhel3]);
	    diag[2] = FFZVertex_->evaluate(scale(),off2,fout[outhel1],interZ);
	    diag[3] = FFPVertex_->evaluate(scale(),off2,fout[outhel1],interG);
	    // sum of diagrams
	    Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    // matrix element
	    output(inhel1,inhel2,outhel1,outhel2,outhel3)=sum;
	    // me2
	    total += norm(sum);
	  }
	}
      }		
    }
  }
  // spin average
  total *= 0.25;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(partons[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(partons[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = 
      {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
       beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    total = output.average(rho[0],rho[1]);
  }
  // divide out the coupling
  total /= norm(FFGVertex_->norm());
  // return the total
  return total*UnitRemoval::InvE2;
}
