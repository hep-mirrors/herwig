// -*- C++ -*-
//
// MEee2gZ2qq.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
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
#include "Herwig++/Shower/Base/MECorrectionBase.h"

using namespace Herwig;

const double MEee2gZ2qq::EPS_=0.00000001;

void MEee2gZ2qq::doinit() {
  HwMEBase::doinit();
  massOption(vector<unsigned int>(2,_massopt));
  rescalingOption(3);
  if(_minflav>_maxflav)
    throw InitException() << "The minimum flavour " << _minflav  
			  << "must be lower the than maximum flavour " << _maxflav
			  << " in MEee2gZ2qq::doinit() " 
			  << Exception::runerror;
  // set the particle data objects
  _Z0=getParticleData(ThePEG::ParticleID::Z0);
  _gamma=getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _theFFZVertex = hwsm->vertexFFZ();
    _theFFPVertex = hwsm->vertexFFP();
  }
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEee2gZ2qq::doinit() the Herwig++ version must be used" 
			     << Exception::runerror;
}

void MEee2gZ2qq::getDiagrams() const {
  // specific the diagrams
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  // setup the processes
  for ( int i =_minflav; i<=_maxflav; ++i ) {
    tcPDPtr qk = getParticleData(i);
    tcPDPtr qb = qk->CC();
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, gamma, 3, qk, 3, qb, -1)));
    add(new_ptr((Tree2toNDiagram(2), em, ep, 1, Z0   , 3, qk, 3, qb, -2)));
  }
}

Energy2 MEee2gZ2qq::scale() const {
  return sHat();
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
  os << _theFFZVertex << _theFFPVertex << _Z0 << _gamma << _minflav 
     << _maxflav << _massopt << alpha_;
}

void MEee2gZ2qq::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theFFPVertex >> _Z0 >> _gamma >> _minflav 
     >> _maxflav >> _massopt >> alpha_;
}

ClassDescription<MEee2gZ2qq> MEee2gZ2qq::initMEee2gZ2qq;
// Definition of the static class description member.
void MEee2gZ2qq::Init() {

  static ClassDocumentation<MEee2gZ2qq> documentation
    ("The MEee2gZ2qq class implements the matrix element for e+e- -> q qbar");

  static Parameter<MEee2gZ2qq,int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The PDG code of the quark with the lowest PDG code to produce.",
     &MEee2gZ2qq::_minflav, 1, 1, 6,
     false, false, Interface::limited);

  static Parameter<MEee2gZ2qq,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The PDG code of the quark with the highest PDG code to produce",
     &MEee2gZ2qq::_maxflav, 5, 1, 6,
     false, false, Interface::limited);

  static Switch<MEee2gZ2qq,unsigned int> interfaceTopMassOption
    ("TopMassOption",
     "Option for the treatment of the top quark mass",
     &MEee2gZ2qq::_massopt, 1, false, false);
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
  int ie(0),ip(1),iqk(2),iqb(3);
  // get the order right
  if(mePartonData()[0]->id()!=11) swap(ie,ip);
  if(mePartonData()[2]->id()<0)   swap(iqk,iqb);
  // compute the spinors
  vector<SpinorWaveFunction> fin,aout;
  vector<SpinorBarWaveFunction>  ain,fout;
  SpinorWaveFunction    ein  (rescaledMomenta()[ie ],mePartonData()[ie ],incoming);
  SpinorBarWaveFunction pin  (rescaledMomenta()[ip ],mePartonData()[ip ],incoming);
  SpinorBarWaveFunction qkout(rescaledMomenta()[iqk],mePartonData()[iqk],outgoing);
  SpinorWaveFunction    qbout(rescaledMomenta()[iqb],mePartonData()[iqb],outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    ein.reset(ix)  ;fin.push_back( ein  );
    pin.reset(ix)  ;ain.push_back( pin  );
    qkout.reset(ix);fout.push_back(qkout);
    qbout.reset(ix);aout.push_back(qbout);
  }
  // compute the matrix element
  double me,lastCont,lastBW;
  HelicityME(fin,ain,fout,aout,me,lastCont,lastBW);
  // save the components
  DVector save;
  save.push_back(lastCont);
  save.push_back(lastBW);
  meInfo(save);
  // add the QCD K-factor
  int Nf = SM().Nf(scale());
  me *= (1.0 + alphaS()/Constants::pi 
	 + (1.986-0.115*Nf)*sqr(alphaS()/Constants::pi));
  // return the answer
  return me;
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
  // wavefunctions for the intermediate particles
  VectorWaveFunction interZ,interG;
  // temporary storage of the different diagrams
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1,outhel2;
  double total[3]={0.,0.};
  for(inhel1=0;inhel1<2;++inhel1) {
    for(inhel2=0;inhel2<2;++inhel2) {
      // intermediate Z
      interZ = _theFFZVertex->evaluate(sHat(),1,_Z0,fin[inhel1],ain[inhel2]);
      // intermediate photon
      interG = _theFFPVertex->evaluate(sHat(),1,_gamma,fin[inhel1],ain[inhel2]);
      for(outhel1=0;outhel1<2;++outhel1) {
	for(outhel2=0;outhel2<2;++outhel2) {		
	  // first the Z exchange diagram
	  diag1 = _theFFZVertex->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interZ);
	  // then the photon exchange diagram
	  diag2 = _theFFPVertex->evaluate(sHat(),aout[outhel2],fout[outhel1],
					  interG);
	  // add up squares of individual terms
	  total[1] += norm(diag1);
	  total[2] += norm(diag2);
	  // the full thing including interference
	  diag1 += diag2;
	  total[0] += norm(diag1);
	  output(inhel1,inhel2,outhel1,outhel2)=diag1;
	}
      }
    }
  }
  // results
  for(int ix=0;ix<3;++ix){total[ix]*=0.75;}
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
  for(unsigned int ix=0;ix<4;++ix)
    dynamic_ptr_cast<SpinfoPtr>(hard[ix]->spinInfo())->setProductionVertex(hardvertex);
}

void MEee2gZ2qq::rebind(const TranslationMap & trans)
  {
  // dummy = trans.translate(dummy);
  _theFFZVertex = trans.translate(_theFFZVertex);
  _theFFPVertex = trans.translate(_theFFPVertex);
  _Z0           = trans.translate(_Z0);
  _gamma        = trans.translate(_gamma);
  HwMEBase::rebind(trans);
}

IVector MEee2gZ2qq::getReferences() {
  IVector ret = HwMEBase::getReferences();
  ret.push_back(_theFFZVertex);
  ret.push_back(_theFFPVertex);
  ret.push_back(_Z0          );
  ret.push_back(_gamma       );
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
  if (newfs[2].e() < getParticleData(ParticleID::g)->constituentMass())
    return;
  // set masses
  for (int i=0; i<2; i++) newfs[i].setMass(qq[i]->mass());
  newfs[2].setMass(ZERO);
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
