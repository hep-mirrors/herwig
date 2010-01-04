// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISBase class.
//

#include "DISBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Utilities/Maths.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include <numeric>

using namespace Herwig;

DISBase::DISBase()  : _initial(6.), _final(3.),
		      _procprob(0.35),
		      _comptonint(0.), _bgfint(0.)
{}

void DISBase::persistentOutput(PersistentOStream & os) const {
  os << _comptonint << _bgfint << _procprob << _initial << _final;
}

void DISBase::persistentInput(PersistentIStream & is, int) {
  is >> _comptonint >> _bgfint >> _procprob  >> _initial >> _final;
}

AbstractClassDescription<DISBase> DISBase::initDISBase;
// Definition of the static class description member.

void DISBase::Init() {
  
  static ClassDocumentation<DISBase> documentation
    ("The DISBase class provides the base class for the "
     "implementation of DIS type processes including the "
     "hard corrections in either the old-fashioned matrix "
     "element correction of POWHEG approaches");

  static Parameter<DISBase,double> interfaceProcessProbability
    ("ProcessProbability",
     "The probabilty of the QCD compton process for the process selection",
     &DISBase::_procprob, 0.3, 0.0, 1.,
     false, false, Interface::limited);

  static Reference<DISBase,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &DISBase::_alpha, false, false, true, false, false);
  
}

void DISBase::doinit() {
  HwMEBase::doinit();
  // integrals of me over phase space
  double r5=sqrt(5.),darg((r5-1.)/(r5+1.)),ath(0.5*log((1.+1./r5)/(1.-1./r5)));
  _comptonint = 2.*(-21./20.-6./(5.*r5)*ath+sqr(Constants::pi)/3.
		    -2.*Math::ReLi2(1.-darg)-2.*Math::ReLi2(1.-1./darg));
  _bgfint = 121./9.-56./r5*ath;
}

void DISBase::initializeMECorrection(ShowerTreePtr, double & initial,
				     double & final) {
  initial = _initial;
  final   = _final;
}

void DISBase::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  static const double eps=1e-6;
  // find the incoming and outgoing quarks and leptons
  ShowerParticlePtr quark[2],lepton[2];
  PPtr hadron;
  tcBeamPtr beam;
  // incoming particles
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data())) {
      hadron = cit->first->original()->parents()[0];
      quark [0] = cit->first->progenitor();
      beam = cit->first->beam();
    }
    else if(LeptonMatcher::Check(cit->first->progenitor()->data())) {
      lepton[0] = cit->first->progenitor();
    }
  }
  tcPDFPtr pdf=beam->pdf();
  assert(beam&&pdf&&quark[0]&&lepton[0]);
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
    if(QuarkMatcher::Check(cit->first->progenitor()->data()))
      quark [1] = cit->first->progenitor();
    else if(LeptonMatcher::Check(cit->first->progenitor()->data()))
      lepton[1] = cit->first->progenitor();
  }
  assert(quark[1]&&lepton[1]);
  // extract the born variables
  Lorentz5Momentum q =lepton[0]->momentum()-lepton[1]->momentum();
  _q2 = -q.m2();
  double  xB = quark[0]->x();
  double  yB = 
    (                    q*quark[0]->momentum())/
    (lepton[0]->momentum()*quark[0]->momentum()); 
  _l = 2./yB-1.;
  // calculate the A coefficient for the correlations
  _acoeff = A(lepton[0]->dataPtr(),lepton[1]->dataPtr(),
	      quark [0]->dataPtr(),quark [1]->dataPtr(),_q2);
  vector<double> azicoeff;
  // select the type of process
  bool BGF = UseRandom::rnd()>_procprob;
  double xp,zp,wgt,x1,x2,x3,xperp;
  tcPDPtr gluon = getParticleData(ParticleID::g);
  // generate a QCD compton process
  if(!BGF) {
    wgt = generateComptonPoint(xp,zp);
    if(xp<eps) return;
    // common pieces
    Energy2 scale = _q2*((1.-xp)*(1-zp)*zp/xp+1.);
    wgt *= 2./3./Constants::pi*_alpha->value(scale)/_procprob;
    // PDF piece
    wgt *= pdf->xfx(beam,quark[0]->dataPtr(),scale,xB/xp)/
           pdf->xfx(beam,quark[0]->dataPtr(),_q2  ,xB);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = ComptonME(xp,x2,xperp,_acoeff,_l,true);
  }
  // generate a BGF process
  else {
    wgt = generateBGFPoint(xp,zp);
    if(xp<eps) return;
    // common pieces 
    Energy2 scale = _q2*((1.-xp)*(1-zp)*zp/xp+1);
    wgt *= 0.25/Constants::pi*_alpha->value(scale)/(1.-_procprob);
    // PDF piece
    wgt *= pdf->xfx(beam,gluon              ,scale,xB/xp)/
           pdf->xfx(beam,quark[0]->dataPtr(),_q2  ,xB);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = BGFME(xp,x2,x3,xperp,_acoeff,_l,true);
  }
  // compute the azimuthal average of the weight
  wgt *= (azicoeff[0]+0.5*azicoeff[2]);
  // decide whether or not to accept the weight
  if(UseRandom::rnd()>wgt) return;
  // if generate generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::applyHardMatrixElementCorrection() to"
				  << " generate phi" << Exception::eventerror;
  // compute the new incoming and outgoing momenta
  Energy Q(sqrt(_q2));
  Lorentz5Momentum p1 = Lorentz5Momentum( 0.5*Q*xperp*cos(phi), 0.5*Q*xperp*sin(phi),
					  -0.5*Q*x2,0.*GeV,0.*GeV);
  p1.rescaleEnergy();
  Lorentz5Momentum p2 = Lorentz5Momentum(-0.5*Q*xperp*cos(phi),-0.5*Q*xperp*sin(phi),
					 -0.5*Q*x3,0.*GeV,0.*GeV);
  p2.rescaleEnergy();
  Lorentz5Momentum pin(0.*GeV,0.*GeV,-0.5*x1*Q,-0.5*x1*Q,0.*GeV);
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleEnergy();
  Lorentz5Momentum pcmf = phadron+0.5/xB*q;
  pcmf.rescaleMass();
  LorentzRotation rot(-pcmf.boostVector());
  Lorentz5Momentum pbeam = rot*phadron;
  Axis axis(pbeam.vect().unit());
  double sinth(sqrt(1.-sqr(axis.z())));
  rot.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  Lorentz5Momentum pl    = rot*lepton[0]->momentum();
  rot.rotateZ(-atan2(pl.y(),pl.x()));
  // we need the Lorentz transform back to the lab
  rot.invert();
  // transform the momenta to lab frame
  pin *= rot;
  p1  *= rot;
  p2  *= rot;
  // test to ensure outgoing particles can be put on-shell
  if(!BGF) {
    if(p1.e()<quark[1]->dataPtr()->constituentMass()) return;
    if(p2.e()<gluon              ->constituentMass()) return;
  }
  else {
    if(p1.e()<quark[1]->dataPtr()      ->constituentMass()) return;
    if(p2.e()<quark[0]->dataPtr()->CC()->constituentMass()) return;
  }
  // create the new particles and add to ShowerTree
  bool isquark = quark[0]->colourLine();
  if(!BGF) {
    PPtr newin  = new_ptr(Particle(*quark[0]));
    newin->set5Momentum(pin);
    PPtr newg   = gluon              ->produceParticle(p2 );
    PPtr newout = quark[1]->dataPtr()->produceParticle(p1 ); 
    ColinePtr col=isquark ? 
      quark[0]->colourLine() : quark[0]->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine());
    // final-state emission
    if(xp>zp) {
      col->removeColoured(newout,!isquark);
      col->addColoured(newin,!isquark);
      col->addColoured(newg,!isquark);
      newline->addColoured(newg,isquark);
      newline->addColoured(newout,!isquark);
    }
    // initial-state emission
    else {
      col->removeColoured(newin ,!isquark);
      col->addColoured(newout,!isquark);
      col->addColoured(newg,isquark);
      newline->addColoured(newg,!isquark);
      newline->addColoured(newin,!isquark);
    }
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[0]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB/xp);
      cit->first->perturbative(xp>zp);
      if(xp<=zp) orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[1]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(xp<=zp);
      if(xp>zp) orig=cit->first->original();
    }
    assert(orig);
    // add the gluon
    ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sg));
    tree->hardMatrixElementCorrection(true);
  }
  else {
    PPtr newin   = gluon                    ->produceParticle(pin);
    PPtr newqbar = quark[0]->dataPtr()->CC()->produceParticle(p2 );
    PPtr newout  = quark[1]->dataPtr()      ->produceParticle(p1 );
    ColinePtr col=isquark ? quark[0]->colourLine() : quark[0]->antiColourLine();
    ColinePtr newline=new_ptr(ColourLine()); 
    col    ->addColoured(newin  ,!isquark);
    newline->addColoured(newin  , isquark);
    col    ->addColoured(newout ,!isquark);
    newline->addColoured(newqbar, isquark);
    PPtr orig;
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator 
	  cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[0]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newin);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
      cit->first->progenitor(sp);
      tree->incomingLines()[cit->first]=sp;
      sp->x(xB/xp);
      cit->first->perturbative(false);
      orig=cit->first->original();
    }
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	  cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit) {
      if(cit->first->progenitor()!=quark[1]) continue;
      // remove old particles from colour line
      col->removeColoured(cit->first->copy(),!isquark);
      col->removeColoured(cit->first->progenitor(),!isquark);
      // insert new particles
      cit->first->copy(newout);
      ShowerParticlePtr sp(new_ptr(ShowerParticle(*newout,1,true)));
      cit->first->progenitor(sp);
      tree->outgoingLines()[cit->first]=sp;
      cit->first->perturbative(true);
    }
    assert(orig);
    // add the (anti)quark
    ShowerParticlePtr sqbar=new_ptr(ShowerParticle(*newqbar,1,true));
    ShowerProgenitorPtr qbar=new_ptr(ShowerProgenitor(orig,newqbar,sqbar));
    qbar->perturbative(false);
    tree->outgoingLines().insert(make_pair(qbar,sqbar));
    tree->hardMatrixElementCorrection(true);
  }
}

bool DISBase::softMatrixElementVeto(ShowerProgenitorPtr initial,
				    ShowerParticlePtr parent, Branching br) {
  bool veto = !UseRandom::rndbool(parent->isFinalState() ? 1./_final : 1./_initial);
  // check if me correction should be applied
  long id[2]={initial->id(),parent->id()};
  if(id[0]!=id[1]||id[1]==ParticleID::g) return veto;
  // get the pT
  Energy pT=br.kinematics->pT();
  // check if hardest so far
  if(pT<initial->highestpT()) return veto;
  double kappa(sqr(br.kinematics->scale())/_q2),z(br.kinematics->z());
  double zk((1.-z)*kappa);
  // final-state
  double wgt(0.);
  if(parent->isFinalState()) {
    double zp=z,xp=1./(1.+z*zk);
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x2 = 1.-(1.-zp)/xp;
    vector<double> azicoeff = ComptonME(xp,x2,xperp,_acoeff,_l,false);
    wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.+sqr(z))/_final;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for FSR in DISMECorrection::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  else {
    double xp = 2.*z/(1.+zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double zp = 0.5* (1.-zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x1 = -1./xp, x2 = 1.-(1.-zp)/xp, x3 = 2.+x1-x2;
    // compton
    if(br.ids[0]!=ParticleID::g) {
      vector<double> azicoeff = ComptonME(xp,x2,xperp,_acoeff,_l,false);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp*(1.-z)/(1.-xp)/(1.+sqr(z))/
	(1.-zp+xp-2.*xp*(1.-zp));
    }
    // BGF
    else {
      vector<double> azicoeff = BGFME(xp,x2,x3,xperp,_acoeff,_l,true);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.-zp+xp-2.*xp*(1.-zp))/(sqr(z)+sqr(1.-z));
    }
    wgt /=_initial;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for ISR in DISMECorrection::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  // if not vetoed
  if(UseRandom::rndbool(wgt)) {
    initial->highestpT(pT);
    return false;
  }
  // otherwise
  parent->setEvolutionScale(br.kinematics->scale());
  return true;
}

double DISBase::generateComptonPoint(double &xp, double & zp) {
  static const double maxwgt = 1.;
  double wgt;
  do {
    xp  = UseRandom::rnd();
    double zpmin = xp, zpmax = 1./(1.+xp*(1.-xp));
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    if(UseRandom::rndbool()) swap(xp,zp);
    double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp,x2=1.-(1.-zp)/xp;
    wgt *= 2.*(1.+sqr(xp)*(sqr(x2)+1.5*xperp2))/(1.-xp)/(1.-zp);
    if(wgt>maxwgt) 
      generator()->logWarning( Exception("DISBase::generateComptonPoint "
					 "Weight greater than maximum", 
					 Exception::warning) );
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return _comptonint;
}

double DISBase::generateBGFPoint(double &xp, double & zp) {
  static const double maxwgt = 25.;
  double wgt;
  do {
    xp = UseRandom::rnd();
    double zpmax = 1./(1.+xp*(1.-xp)), zpmin = 1.-zpmax;
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    double x1 = -1./xp;
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.+x1-x2;
    double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
    wgt *= sqr(xp)/(1.-zp)*(sqr(x3)+sqr(x2)+3.*xperp2);
    if(wgt>maxwgt) 
      generator()->logWarning( Exception("DISBase::generateBGFPoint "
					 "Weight greater than maximum", 
					 Exception::warning) );
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return _bgfint;
//   static const double maxwgt = 2.,npow=0.34,ac=1.0;
//   double wgt;
//   do {
//     double rho = UseRandom::rnd();
//     xp = 1.-pow(rho,1./(1.-npow));
//     wgt = (sqr(xp)+ac+sqr(1.-xp));
//     if(wgt>1.+ac) cerr << "testing violates BGF maxA " << wgt << "\n";
//   }
//   while(wgt<UseRandom::rnd()*(1.+ac));
//   double xpwgt = -((6.-5.*npow+sqr(npow))*ac-3.*npow+sqr(npow)+4) 
//     /(sqr(npow)*(npow-6.)+11.*npow-6.);
//   xpwgt *= pow(1.-xp,npow)/wgt;
//   double xp2(sqr(xp)),lxp(log(xp)),xp4(sqr(xp2)),lxp1(log(1.-xp));
//   double zpwgt = (2.*xp4*(lxp+lxp1-3.)+4.*xp*xp2*(3.-lxp-lxp1)
// 		  +xp2*(-13.+lxp+lxp1)+xp*(+7.+lxp+lxp1)-lxp-lxp1-1.)/(1.+xp-xp2);
//   do {
//     double zpmax = 1./(1.+xp*(1.-xp)), zpmin = 1.-zpmax;
//     zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
//     wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
//     double x1 = -1./xp;
//     double x2 = 1.-(1.-zp)/xp;
//     double x3 = 2.+x1-x2;
//     double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
//     wgt *= sqr(xp)/(1.-zp)*(sqr(x3)+sqr(x2)+3.*xperp2);
//     if(wgt>maxwgt*zpwgt) cerr << "testing violates BGF maxB " << wgt/xpwgt << "\n";
//   }
//   while(wgt<UseRandom::rnd()*maxwgt);
//   return zpwgt*xpwgt;
}

vector<double> DISBase::ComptonME(double xp, double x2, double xperp,
				  double A, double l, bool norm) {
  vector<double> output(3,0.);
  double cos2 =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2 = xperp/sqrt(sqr(x2)+sqr(xperp));
  double root = sqrt(sqr(l)-1.);
  output[0] = sqr(cos2)-A*cos2*l+sqr(l);
  output[1] = A*cos2*root*sin2-2.*l*root*sin2;
  output[2] = sqr(root)*sqr(sin2);
  double lo(1+A*l+sqr(l));
  double denom = norm ? 1.+sqr(xp)*(sqr(x2)+1.5*sqr(xperp)) : 1.;
  double fact  = sqr(xp)*(sqr(x2)+sqr(xperp))/lo;
  for(unsigned int ix=0;ix<output.size();++ix) 
    output[ix] = ((ix==0 ? 1. : 0.) +fact*output[ix])/denom;
  return output;
}

vector<double> DISBase::BGFME(double xp, double x2, double x3, 
			      double xperp, double A, double l,
			      bool norm) {
  vector<double> output(3,0.);
  double cos2  =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2  = xperp/sqrt(sqr(x2)+sqr(xperp));
  double fact2 = sqr(xp)*(sqr(x2)+sqr(xperp));
  double cos3  =   x3 /sqrt(sqr(x3)+sqr(xperp));
  double sin3  = xperp/sqrt(sqr(x3)+sqr(xperp));
  double fact3 = sqr(xp)*(sqr(x3)+sqr(xperp));
  double root = sqrt(sqr(l)-1.);
  output[0] = fact3*(sqr(cos3)-A*cos3*l+sqr(l))+
    fact2*(sqr(cos2)-A*cos2*l+sqr(l));
  output[1] =-fact3*(A*cos3*root*sin3-2.*l*root*sin3)+
    fact2*(A*cos2*root*sin2-2.*l*root*sin2);
  output[2] = fact3*(sqr(root)*sqr(sin3))+
    fact2*(sqr(root)*sqr(sin2));
  double lo(1+A*l+sqr(l));
  double denom = norm ? sqr(xp)*(sqr(x3)+sqr(x2)+3.*sqr(xperp))*lo : lo;
  for(unsigned int ix=0;ix<output.size();++ix) output[ix] /= denom;
  return output;
}
