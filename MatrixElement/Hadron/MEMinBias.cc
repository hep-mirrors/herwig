// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEMinBias class.
//

#include "MEMinBias.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
//#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Interface/Parameter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

void MEMinBias::getDiagrams() const {
  int maxflav(2);
  // Pomeron data
  tcPDPtr pom = getParticleData(990);

  for ( int i = 1; i <= maxflav; ++i ) {
    for( int j=1; j <= i; ++j){
      tcPDPtr q1 = getParticleData(i);
      tcPDPtr q1b = q1->CC();
      tcPDPtr q2 = getParticleData(j);
      tcPDPtr q2b = q2->CC();

      // For each flavour we add:
      //qq -> qq
      add(new_ptr((Tree2toNDiagram(3), q1, pom, q2, 1, q1, 2, q2, -1)));
      //qqb -> qqb
      add(new_ptr((Tree2toNDiagram(3), q1, pom, q2b, 1, q1, 2, q2b, -2)));
      //qbqb -> qbqb
      add(new_ptr((Tree2toNDiagram(3), q1b, pom, q2b, 1, q1b, 2, q2b, -3)));
    }
  }
}

Energy2 MEMinBias::scale() const {
  return sqr(10*GeV);
}

int MEMinBias::nDim() const {
  return 0;
}

void MEMinBias::setKinematics() {
  HwMEBase::setKinematics(); // Always call the base class method first.
}

bool MEMinBias::generateKinematics(const double *) {
  // generate the masses of the particles
  for ( int i = 2, N = meMomenta().size(); i < N; ++i ) {
    meMomenta()[i] = Lorentz5Momentum(mePartonData()[i]->generateMass());
  }

  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } catch ( ImpossibleKinematics ) {
    return false;
  }

  Energy pt = ZERO;
  meMomenta()[2].setVect(Momentum3( pt,  pt, q));
  meMomenta()[3].setVect(Momentum3(-pt, -pt, -q));

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  jacobian(1.0);
  return true;
}

double MEMinBias::me2() const {
  //tuned so it gives the correct normalization for xmin = 0.11
  return csNorm_*(sqr(generator()->maximumCMEnergy())/GeV2);
}

CrossSection MEMinBias::dSigHatDR() const {
  return me2()*jacobian()/sHat()*sqr(hbarc);
}

unsigned int MEMinBias::orderInAlphaS() const {
  return 2;
}

unsigned int MEMinBias::orderInAlphaEW() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
MEMinBias::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(1.0, i);

  return sel;
}

Selector<const ColourLines *>
MEMinBias::colourGeometries(tcDiagPtr diag) const {

  static ColourLines qq("1 4, 3 5");
  static ColourLines qqb("1 4, -3 -5");
  static ColourLines qbqb("-1 -4, -3 -5");

  Selector<const ColourLines *> sel;
  
  switch(diag->id()){
  case -1:
    sel.insert(1.0, &qq);
    break;
  case -2:
    sel.insert(1.0, &qqb);
    break;
  case -3:
    sel.insert(1.0, &qbqb);
    break;
  }
  return sel;
}


IBPtr MEMinBias::clone() const {
  return new_ptr(*this);
}

IBPtr MEMinBias::fullclone() const {
  return new_ptr(*this);
}


ClassDescription<MEMinBias> MEMinBias::initMEMinBias;
// Definition of the static class description member.

void MEMinBias::persistentOutput(PersistentOStream & os) const {
  os << csNorm_;
}

void MEMinBias::persistentInput(PersistentIStream & is, int) {
  is >> csNorm_;
}

void MEMinBias::Init() {

  static ClassDocumentation<MEMinBias> documentation
    ("There is no documentation for the MEMinBias class");
     
  static Parameter<MEMinBias,double> interfacecsNorm
    ("csNorm",
     "Normalization of the min-bias cross section.",
     &MEMinBias::csNorm_, 
     1.0, 0.0, 100.0, 
     false, false, Interface::limited);
}

