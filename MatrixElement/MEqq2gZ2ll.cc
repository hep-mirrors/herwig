// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEqq2gZ2ll class.
//

#include "MEqq2gZ2ll.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Timer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Herwig;

MEqq2gZ2ll::MEqq2gZ2ll()
  : coefs(20), mZ2(0.0*GeV2), GZ2(0.0*GeV2), lastCont(0.0), lastBW(0.0) {}

MEqq2gZ2ll::MEqq2gZ2ll(const MEqq2gZ2ll & x)
  : ME2to2QCD(x), coefs(x.coefs), mZ2(x.mZ2), GZ2(x.GZ2),
    lastCont(x.lastCont), lastBW(x.lastBW) {}

MEqq2gZ2ll::~MEqq2gZ2ll() {}

unsigned int MEqq2gZ2ll::orderInAlphaS() const {
  return 0;
}

unsigned int MEqq2gZ2ll::orderInAlphaEW() const {
  return 2;
}

void MEqq2gZ2ll::getDiagrams() const {
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  for(int i = 1; i <= maxFlavour(); ++i) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    for(long lid = ParticleID::eminus; lid <= ParticleID::tauminus; lid+=2) {
      tcPDPtr l = getParticleData(lid);
      tcPDPtr lb = l->CC();
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gamma, 3, l, 3, lb, -1)));
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0, 3, l, 3, lb, -2)));
    }
  }
}

Energy2 MEqq2gZ2ll::scale() const {
  return sHat();
}

double MEqq2gZ2ll::me2() const {

  Energy2 m2 = meMomenta()[2].mass2();
  //  Energy2 p1p3 = 0.5*(m2 - tHat());
  //  Energy2 p1p2 = 0.5*sHat();
  Energy2 p1p3 = meMomenta()[0].dot(meMomenta()[2]);
  Energy2 p1p2 = meMomenta()[0].dot(meMomenta()[1]);
  Energy4 pt2 = sqr(p1p3);
  Energy4 pts = p1p3*p1p2;
  Energy4 ps2 = sqr(p1p2);
  Energy4 psm = p1p2*m2;

  int up = abs(mePartonData()[0]->id() + 1)%2;

  lastCont =
    (coefs[0 + up]*(pt2 - pts) + coefs[2 + up]*(ps2 + psm))/sqr(sHat());
  double intr = 0.25*(coefs[4 + up]*pt2 + coefs[6 + up]*pts +
		 coefs[8 + up]*ps2 + coefs[10 + up]*psm)*
    (sHat() - mZ2)/(sHat()*(sqr(sHat() - mZ2) + mZ2*GZ2));
  lastBW = 0.25*(coefs[12 + up]*pt2 + coefs[14 + up]*pts +
	    coefs[16 + up]*ps2 + coefs[18 + up]*psm)/
    (sqr(sHat() - mZ2) + mZ2*GZ2);

  double alphaS = SM().alphaS(scale());
  //int Nf = SM().Nf(scale());
  int Nf = 3;
  DVector save;
  meInfo(save << lastCont << lastBW);
  return (lastCont + intr + lastBW)*sqr(SM().alphaEM(scale()))*
    (1.0 + alphaS/pi + (1.986-0.115*Nf)*sqr(alphaS/pi));
}

Selector<MEqq2gZ2ll::DiagramIndex>
MEqq2gZ2ll::diagrams(const DiagramVector & diags) const {
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
MEqq2gZ2ll::colourGeometries(tcDiagPtr diag) const {

  static ColourLines c("1 -2");

  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c);
  return sel;
}

IBPtr MEqq2gZ2ll::clone() const {
  return new_ptr(*this);
}

IBPtr MEqq2gZ2ll::fullclone() const {
  return new_ptr(*this);
}

void MEqq2gZ2ll::doinit() throw(InitException) {
  double C = sqr(4.0*Constants::pi)/3.0;
  double SW2 = SM().sin2ThetaW();
  double SW4 = sqr(SW2);
  double SW6 = SW2*SW4;
  double SW8 = SW2*SW6;
  double CW2 = 1.0 - SW2;
  coefs[0] =  16.0*C;
  coefs[1] =  64.0*C;
  coefs[2] =   8.0*C;
  coefs[3] =  32.0*C;
  C /= (CW2*SW2);
  coefs[4] =  4.0*(32.0*SW4 - 32.0*SW2 +  6.0)*C;
  coefs[5] =  8.0*(64.0*SW4 - 40.0*SW2 +  6.0)*C;
  coefs[6] = -4.0*(32.0*SW4 - 32.0*SW2 + 12.0)*C;
  coefs[7] = -8.0*(64.0*SW4 - 40.0*SW2 + 12.0)*C;
  coefs[8] =  4.0*(16.0*SW4 - 16.0*SW2 +  6.0)*C;
  coefs[9] =  8.0*(32.0*SW4 - 20.0*SW2 +  6.0)*C;
  coefs[10] = 4.0*(16.0*SW4 - 16.0*SW2 +  3.0)*C;
  coefs[11] = 8.0*(32.0*SW4 - 20.0*SW2 +  3.0)*C;
  C /= (CW2*SW2);
  coefs[12] =  ( 64.0*SW8 - 128.0*SW6 + 128.0*SW4 -  48.0*SW2 +  9.0)*C;
  coefs[13] =  (256.0*SW8 - 320.0*SW6 + 200.0*SW4 -  60.0*SW2 +  9.0)*C;
  coefs[14] = -( 64.0*SW8 - 128.0*SW6 + 176.0*SW4 -  96.0*SW2 + 18.0)*C;
  coefs[15] = -(256.0*SW8 - 320.0*SW6 + 296.0*SW4 - 120.0*SW2 + 18.0)*C;
  coefs[16] =  ( 32.0*SW8 -  64.0*SW6 +  88.0*SW4 -  48.0*SW2 +  9.0)*C;
  coefs[17] =  (128.0*SW8 - 160.0*SW6 + 148.0*SW4 -  60.0*SW2 +  9.0)*C;
  coefs[18] =  ( 32.0*SW8 -  64.0*SW6 +  28.0*SW4 -   6.0*SW2)*C;
  coefs[19] =  (128.0*SW8 - 160.0*SW6 +  64.0*SW4 -  12.0*SW2)*C;

  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  mZ2 = sqr(Z0->mass());
  GZ2 = sqr(Z0->width());

  ME2to2QCD::doinit();
}

void MEqq2gZ2ll::persistentOutput(PersistentOStream & os) const {
  os << coefs << ounit(mZ2, GeV2) << ounit(GZ2, GeV2) << lastCont << lastBW;
}

void MEqq2gZ2ll::persistentInput(PersistentIStream & is, int) {
  is >> coefs >> iunit(mZ2, GeV2) >> iunit(GZ2, GeV2) >> lastCont >> lastBW;
}

ClassDescription<MEqq2gZ2ll> MEqq2gZ2ll::initMEqq2gZ2ll;

void MEqq2gZ2ll::Init() {

  static ClassDocumentation<MEqq2gZ2ll> documentation
    ("The \\classname{MEqq2gZ2ll} class implements the full"
     "$q + \\bar{q} \\rightarrow \\gamma/Z^0 \\rightarrow \\ell^+ + \\ell^-$ "
     "matrix element including the interference terms.");

}

