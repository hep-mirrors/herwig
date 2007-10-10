// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CSVectorBosonQQbarHardGenerator class.
//

#include "CSVectorBosonQQbarHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace Herwig;

void CSVectorBosonQQbarHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << ounit(_ptmin,GeV);
}

void CSVectorBosonQQbarHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> iunit(_ptmin,GeV);
}

ClassDescription<CSVectorBosonQQbarHardGenerator> 
CSVectorBosonQQbarHardGenerator::initCSVectorBosonQQbarHardGenerator;
// Definition of the static class description member.

void CSVectorBosonQQbarHardGenerator::Init() {

  static ClassDocumentation<CSVectorBosonQQbarHardGenerator> documentation
    ("There is no documentation for the CSVectorBosonQQbarHardGenerator class");

  static Reference<CSVectorBosonQQbarHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &CSVectorBosonQQbarHardGenerator::_alphaS, false, false, true, false, false);

  static Parameter<CSVectorBosonQQbarHardGenerator,Energy> interfacepTMin
    ("pTMin",
     "The minimum pt of the hardest emission",
     &CSVectorBosonQQbarHardGenerator::_ptmin, GeV, 0.5*GeV, 0.2*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

bool CSVectorBosonQQbarHardGenerator::canHandle(ShowerTreePtr tree) {
  if(tree->incomingLines().size()!=1) return false;    
  if((tree->incomingLines().begin()->first->id()==22)&&
  (tree->incomingLines().begin()->first->progenitor()->id()==23)) return false;
  map<ShowerProgenitorPtr,tShowerParticlePtr> outgoing=tree->outgoingLines();
  if(outgoing.size()!=2) return false;
  if(abs(outgoing.begin()->first->progenitor()->id())>6)  return false;
  if(outgoing.begin()->first->progenitor()->id()!=
     -1*outgoing.rbegin()->first->progenitor()->id())     return false;
  return true;
}

NasonTreePtr CSVectorBosonQQbarHardGenerator::generateHardest(ShowerTreePtr tree) {
  cerr << *generator()->currentEvent();
  // Get the progenitors: Q and Qbar.
  vector<tcPDPtr> partons(2);
  ShowerProgenitorPtr 
    QProgenitor   =tree->outgoingLines().begin()->first,
    QbarProgenitor=tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap(QProgenitor   ,QbarProgenitor);
  partons[0]=QProgenitor->progenitor()->dataPtr();
  partons[1]=QbarProgenitor->progenitor()->dataPtr();
  // incoming particles
  ShowerProgenitorPtr incoming=tree->incomingLines().begin()->first;
  tPPtr fin;
  for(unsigned int ix=0;ix<incoming->original()->parents().size();++ix) {
    if(incoming->original()->parents()[ix]->id()>0) 
      fin=incoming->original()->parents()[ix];
  }
  assert(fin);
  Lorentz5Momentum pin=fin->momentum(); 
  cerr << "testing " << *fin << "\n";
  cerr << "testing " << *incoming->original() << "\n";
  // generate the variables for the emission
  Lorentz5Momentum pboson=
    QProgenitor   ->progenitor()->momentum()+
    QbarProgenitor->progenitor()->momentum();
  pin.boost(-pboson.boostVector());
  cerr << "testing pboson " << pboson/GeV << "\n";
  cerr << "testing pin    " << pin/GeV << "\n";
  Axis dout(QProgenitor   ->progenitor()->momentum().vect().unit());
  Axis din (pin.vect().unit());
  dout.rotateUz(din);
  double theta = dout.theta();
  double phi   = dout.phi();

  pboson.rescaleMass();
  Energy2 s=pboson.mass2();
  cerr << "testing " << s/GeV2 << "\n";
  generate(s,theta,phi);

  cerr << "testing in generate hardest\n";
  return NasonTreePtr();
}

void CSVectorBosonQQbarHardGenerator::generate(Energy2 s,double theta,
					       double phiout) {
  double c = 2.*_alphaS->overestimateValue()/3./Constants::pi;
  Energy2 pt2(s);
  double y,x1,x2,z,cosd,sind,phi,cos1,wgt;
  double cbeam = cos(theta), sbeam = sin(theta);
  double C1,C2;
  do {
    // pt2 and y from overestimate of the function
    pt2 = s*exp(-sqrt(sqr(log(pt2/s))-log(UseRandom::rnd())/c));
    y = exp(log(pt2/s)*(1.-0.5*UseRandom::rnd()));
    // random phi
    phi = 2.*Constants::pi*UseRandom::rnd();
    // calculate z
    z = (y-pt2/s)/y/(1.-y);
    // compute R
    x2 = 1.-y;
    x1 = x2*z+y;
    cosd = -(2.*x1+2.*x2-2.-x1*x2)/x1/x2;
    sind = sin(acos(cosd));
    cos1 = -sbeam*sind*cos(phi)+cosd*cbeam;
    double R = (sqr(x1)*(1.+sqr(cos1 ))+
		sqr(x2)*(1.+sqr(cbeam)))/(1.-x1)/(1.-x2);
    R /= 1.+sqr(cbeam);
    cerr << "testing cos " << cos1 << " " << cbeam << "\n";

    cerr << "testing x" << x1 << " " << x2 << "\n";
    cerr << "testing y" << y  << " " << z << "\n";

    C1 = 1./(1.-x2)*(2./(2.-x1-x2)-(1.+x1))+(1.-x1)/x2;
    C2 = 1./(1.-x1)*(2./(2.-x1-x2)-(1.+x2))+(1.-x2)/x1;

    wgt = 0.25*_alphaS->ratio(pt2)*R*y*(1.-z*(1.-y))*C1/(C1+C2);
    cerr << "testing pieces " 
	 << _alphaS->ratio(pt2) << " " << R << " " 
	 << 1./(0.25*y*(1.-z*(1.-y))) << "\n";
    if(wgt>1.) cerr << "testing weight violates max " << wgt << "\n";
  }
  while(pt2>sqr(_ptmin)&&UseRandom::rnd()>wgt);
  unsigned int iemit = UseRandom::rnd()<C1/(C1+C2) ? 1 : 2;
  Energy roots = 0.5*sqrt(s);
  Lorentz5Momentum p1(x1*roots*sind*cos(phi),
		      x1*roots*sind*sin(phi),x1*roots*cosd,x1*roots);
  Lorentz5Momentum p2(0.*GeV,0.*GeV,x2*roots,x2*roots);
  Lorentz5Momentum p3(Lorentz5Momentum(2.*roots)-p1-p2);
  cerr << "testing momenta " << p1/GeV << " " << p1.m()/GeV << "\n";
  cerr << "testing momenta " << p2/GeV << " " << p2.m()/GeV << "\n";
  cerr << "testing momenta " << p3/GeV << " " << p3.m()/GeV << "\n";
  cerr << "testing iemit" << iemit << "\n";
  LorentzRotation rot;
  // 2 is spectator
  if(iemit==1) {
    rot.setRotateY(Constants::pi+theta);
    rot.rotateZ(phiout);
  }
  // 1 is spectator
  else {
    rot.setRotateY(theta);
    rot.rotateZ(phiout);
  }
  p1.transform(rot);
  p2.transform(rot);
  p3.transform(rot);
  cerr << "testing momenta " << p1/GeV << " " << p1.m()/GeV << "\n";
  cerr << "testing momenta " << p2/GeV << " " << p2.m()/GeV << "\n";
  cerr << "testing momenta " << p3/GeV << " " << p3.m()/GeV << "\n";
}
