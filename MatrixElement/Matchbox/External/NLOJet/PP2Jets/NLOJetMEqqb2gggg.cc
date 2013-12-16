// -*- C++ -*-
//
// NLOJetMEqqb2gggg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEqqb2gggg class.
//

#include "NLOJetMEqqb2gggg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEqqb2gggg::NLOJetMEqqb2gggg() {}

NLOJetMEqqb2gggg::~NLOJetMEqqb2gggg() {}

IBPtr NLOJetMEqqb2gggg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEqqb2gggg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEqqb2gggg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEqqb2gggg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEqqb2gggg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEqqb2gggg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEqqb2gggg("Herwig::NLOJetMEqqb2gggg", "HwMatchboxNLOJet.so");

void NLOJetMEqqb2gggg::Init() {

  static ClassDocumentation<NLOJetMEqqb2gggg> documentation
    ("NLOJetMEqqb2gggg");

}


void NLOJetMEqqb2gggg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr q = quark[xi1]; 
    tcPDPtr qb = antiquark[xi1]; 
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 3, g, 4, g, 4, g, 5, g, 5, g, -1)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 3, g, 4, g, 5, g, 4, g, 5, g, -2)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 3, g, 4, g, 5, g, 5, g, 4, g, -3)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 1, g, 3, g, 5, g, 5, g, -4)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 1, g, 5, g, 3, g, 5, g, -5)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 1, g, 5, g, 5, g, 3, g, -6)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 3, g, 1, g, 5, g, 5, g, -7)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 5, g, 1, g, 3, g, 5, g, -8)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 5, g, 1, g, 5, g, 3, g, -9)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 3, g, 5, g, 1, g, 5, g, -10)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 5, g, 3, g, 1, g, 5, g, -11)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 5, g, 5, g, 1, g, 3, g, -12)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 3, g, 5, g, 5, g, 1, g, -13)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 5, g, 3, g, 5, g, 1, g, -14)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 2, g, 5, g, 5, g, 3, g, 1, g, -15)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 5, g, 5, g, 3, g, 4, g, -16)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 5, g, 5, g, 4, g, 3, g, -17)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 5, g, 3, g, 5, g, 4, g, -18)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 5, g, 4, g, 5, g, 3, g, -19)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 5, g, 3, g, 4, g, 5, g, -20)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 5, g, 4, g, 3, g, 5, g, -21)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 3, g, 5, g, 5, g, 4, g, -22)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 4, g, 5, g, 5, g, 3, g, -23)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 3, g, 5, g, 4, g, 5, g, -24)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 4, g, 5, g, 3, g, 5, g, -25)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 3, g, 4, g, 5, g, 5, g, -26)));
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, g, 4, g, 3, g, 5, g, 5, g, -27)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 1, g, 4, g, 2, g, 3, g, -28)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 1, g, 4, g, 3, g, 2, g, -29)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 1, g, 2, g, 4, g, 3, g, -30)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 1, g, 3, g, 4, g, 2, g, -31)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 1, g, 2, g, 3, g, 4, g, -32)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 1, g, 3, g, 2, g, 4, g, -33)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 1, g, 5, g, 5, g, 4, g, -34)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 1, g, 5, g, 5, g, 2, g, -35)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 1, g, 5, g, 4, g, 5, g, -36)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 1, g, 5, g, 2, g, 5, g, -37)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 1, g, 4, g, 5, g, 5, g, -38)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 1, g, 2, g, 5, g, 5, g, -39)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 4, g, 1, g, 2, g, 3, g, -40)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 4, g, 1, g, 3, g, 2, g, -41)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 2, g, 1, g, 4, g, 3, g, -42)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 3, g, 1, g, 4, g, 2, g, -43)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 2, g, 1, g, 3, g, 4, g, -44)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 3, g, 1, g, 2, g, 4, g, -45)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 5, g, 1, g, 5, g, 4, g, -46)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 5, g, 1, g, 5, g, 2, g, -47)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 5, g, 1, g, 4, g, 5, g, -48)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 5, g, 1, g, 2, g, 5, g, -49)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 4, g, 1, g, 5, g, 5, g, -50)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 2, g, 1, g, 5, g, 5, g, -51)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 4, g, 2, g, 1, g, 3, g, -52)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 4, g, 3, g, 1, g, 2, g, -53)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 2, g, 4, g, 1, g, 3, g, -54)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 3, g, 4, g, 1, g, 2, g, -55)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 2, g, 3, g, 1, g, 4, g, -56)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 3, g, 2, g, 1, g, 4, g, -57)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 5, g, 5, g, 1, g, 4, g, -58)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 5, g, 5, g, 1, g, 2, g, -59)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 5, g, 4, g, 1, g, 5, g, -60)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 5, g, 2, g, 1, g, 5, g, -61)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 4, g, 5, g, 1, g, 5, g, -62)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 2, g, 5, g, 1, g, 5, g, -63)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 4, g, 2, g, 3, g, 1, g, -64)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 4, g, 3, g, 2, g, 1, g, -65)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 2, g, 4, g, 3, g, 1, g, -66)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 3, g, 4, g, 2, g, 1, g, -67)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 2, g, 3, g, 4, g, 1, g, -68)));
    addSafe(new_ptr((Tree2toNDiagram(5), q, q, q, q, qb, 3, g, 2, g, 4, g, 1, g, -69)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 5, g, 5, g, 4, g, 1, g, -70)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 5, g, 5, g, 2, g, 1, g, -71)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 5, g, 4, g, 5, g, 1, g, -72)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 5, g, 2, g, 5, g, 1, g, -73)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, g, 4, g, 5, g, 5, g, 1, g, -74)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 3, g, 2, g, 5, g, 5, g, 1, g, -75)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 2, g, 5, g, 5, g, 4, g, -76)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 3, g, 5, g, 5, g, 2, g, -77)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 2, g, 5, g, 4, g, 5, g, -78)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 3, g, 5, g, 2, g, 5, g, -79)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 2, g, 4, g, 5, g, 5, g, -80)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 3, g, 2, g, 5, g, 5, g, -81)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 5, g, 2, g, 5, g, 4, g, -82)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 5, g, 3, g, 5, g, 2, g, -83)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 5, g, 2, g, 4, g, 5, g, -84)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 5, g, 3, g, 2, g, 5, g, -85)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 4, g, 2, g, 5, g, 5, g, -86)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 2, g, 3, g, 5, g, 5, g, -87)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 5, g, 5, g, 2, g, 4, g, -88)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 5, g, 5, g, 3, g, 2, g, -89)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 5, g, 4, g, 2, g, 5, g, -90)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 5, g, 2, g, 3, g, 5, g, -91)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 4, g, 5, g, 2, g, 5, g, -92)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 2, g, 5, g, 3, g, 5, g, -93)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 5, g, 5, g, 4, g, 2, g, -94)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 5, g, 5, g, 2, g, 3, g, -95)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 5, g, 4, g, 5, g, 2, g, -96)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 5, g, 2, g, 5, g, 3, g, -97)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, g, 4, g, 5, g, 5, g, 2, g, -98)));
    addSafe(new_ptr((Tree2toNDiagram(4), q, q, q, qb, 1, g, 2, g, 5, g, 5, g, 3, g, -99)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 2, g, 4, g, 4, g, 5, g, 5, g, -100)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 1, g, 4, g, 4, g, 5, g, 5, g, -101)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 2, g, 4, g, 5, g, 4, g, 5, g, -102)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 1, g, 4, g, 5, g, 4, g, 5, g, -103)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 2, g, 4, g, 5, g, 5, g, 4, g, -104)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 1, g, 4, g, 5, g, 5, g, 4, g, -105)));
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEqqb2gggg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv1at1[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv1at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv1at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv1at4[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv1at5[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv1at6[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv1at7[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv2at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv2at1[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv2at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv2at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv2at4[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv2at5[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv2at6[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv2at7[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv3at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv3at1[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv3at2[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv3at3[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv3at4[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv3at5[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv3at6[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv3at7[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv4at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv4at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv5at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv5at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv6at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv6at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv7at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv7at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv8at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv8at1[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv9at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv9at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv10at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv10at1[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv11at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv11at1[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv12at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv12at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv13at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv13at1[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv14at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv14at1[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv15at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv15at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv16at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv16at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv16at2[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv16at3[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv16at4[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv16at5[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv16at6[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv16at7[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv17at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv17at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv17at2[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv17at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv17at4[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv17at5[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv17at6[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv17at7[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv18at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv18at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv18at2[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv18at3[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv18at4[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv18at5[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv18at6[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv18at7[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv19at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv19at1[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv19at2[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv19at3[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv19at4[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv19at5[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv19at6[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv19at7[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv20at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv20at1[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv20at2[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv20at3[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv20at4[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv20at5[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv20at6[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv20at7[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv21at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv21at1[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv21at2[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv21at3[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv21at4[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv21at5[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv21at6[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv21at7[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv22at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv22at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv22at2[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv22at3[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv22at4[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv22at5[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv22at6[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv22at7[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv23at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv23at1[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv23at2[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv23at3[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv23at4[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv23at5[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv23at6[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv23at7[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv24at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv24at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv24at2[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv24at3[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv24at4[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv24at5[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv24at6[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv24at7[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv25at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv25at1[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv25at2[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv25at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv25at4[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv25at5[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv25at6[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv25at7[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv26at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv26at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv26at2[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv26at3[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv26at4[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv26at5[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv26at6[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv26at7[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv27at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv27at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv27at2[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv27at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv27at4[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv27at5[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv27at6[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv27at7[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv28at0[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv29at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv30at0[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv31at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv32at0[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv33at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv34at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv34at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv34at2[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv34at3[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv35at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv35at1[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv36at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv36at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv36at2[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv36at3[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv37at0[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv37at1[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv38at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv38at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv38at2[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv38at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv39at0[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv39at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv40at0[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv41at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv42at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv43at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv44at0[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv45at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv46at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv46at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv46at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv46at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv47at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv47at1[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv48at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv48at1[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv48at2[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv48at3[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv49at0[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv49at1[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv50at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv50at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv50at2[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv50at3[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv51at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv51at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv52at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv53at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv54at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv55at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv56at0[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv57at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv58at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv58at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv58at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv58at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv59at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv59at1[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv60at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv60at1[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv60at2[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv60at3[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv61at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv61at1[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv62at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv62at1[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv62at2[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv62at3[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv63at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv63at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv64at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv65at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv66at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv67at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv68at0[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv69at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv70at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv70at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv70at2[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv70at3[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv71at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv71at1[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv72at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv72at1[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv72at2[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv72at3[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv73at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv73at1[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv74at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv74at1[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv74at2[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv74at3[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv75at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv75at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv76at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv76at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv76at2[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv76at3[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv77at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv77at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv78at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv78at1[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv78at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv78at3[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv79at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv79at1[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv80at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv80at1[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv80at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv80at3[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv81at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv81at1[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv82at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv82at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv82at2[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv82at3[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv83at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv83at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv84at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv84at1[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv84at2[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv84at3[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv85at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv85at1[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv86at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv86at1[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv86at2[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv86at3[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv87at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv87at1[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv88at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv88at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv88at2[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv88at3[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv89at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv89at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv90at0[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv90at1[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv90at2[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv90at3[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv91at0[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv91at1[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv92at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv92at1[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv92at2[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv92at3[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv93at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv93at1[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv94at0[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv94at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv94at2[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv94at3[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv95at0[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv95at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv96at0[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv96at1[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv96at2[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv96at3[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv97at0[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv97at1[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv98at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv98at1[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv98at2[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv98at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv99at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv99at1[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv100at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv100at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv100at2[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv100at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv101at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv101at1[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv101at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv101at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv102at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv102at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv102at2[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv102at3[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv103at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv103at1[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv103at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv103at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv104at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv104at1[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv104at2[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv104at3[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv105at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv105at1[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv105at2[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv105at3[]  = { 0, 4, 1, 3, 2, -1, -999};

  Selector <MEBase::DiagramIndex> sel;
  for(MEBase::DiagramIndex i=0; i < diags.size(); ++i){
    if( diags[i]->id() == -1 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at2, sizeof(cv1at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at3, sizeof(cv1at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at4, sizeof(cv1at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at5, sizeof(cv1at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at6, sizeof(cv1at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at7, sizeof(cv1at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at2, sizeof(cv2at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at3, sizeof(cv2at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at4, sizeof(cv2at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at5, sizeof(cv2at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at6, sizeof(cv2at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at7, sizeof(cv2at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at2, sizeof(cv3at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at3, sizeof(cv3at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at4, sizeof(cv3at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at5, sizeof(cv3at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at6, sizeof(cv3at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at7, sizeof(cv3at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -4 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -5 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -6 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -7 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at1, sizeof(cv7at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -8 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -9 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -10 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -11 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -12 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at1, sizeof(cv12at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -13 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at1, sizeof(cv13at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -14 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at1, sizeof(cv14at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -15 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -16 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv16at0, sizeof(cv16at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at1, sizeof(cv16at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at2, sizeof(cv16at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at3, sizeof(cv16at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at4, sizeof(cv16at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at5, sizeof(cv16at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at6, sizeof(cv16at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at7, sizeof(cv16at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -17 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv17at0, sizeof(cv17at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at1, sizeof(cv17at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at2, sizeof(cv17at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at3, sizeof(cv17at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at4, sizeof(cv17at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at5, sizeof(cv17at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at6, sizeof(cv17at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at7, sizeof(cv17at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -18 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv18at0, sizeof(cv18at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at1, sizeof(cv18at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at2, sizeof(cv18at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at3, sizeof(cv18at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at4, sizeof(cv18at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at5, sizeof(cv18at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at6, sizeof(cv18at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at7, sizeof(cv18at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -19 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at1, sizeof(cv19at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at2, sizeof(cv19at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at3, sizeof(cv19at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at4, sizeof(cv19at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at5, sizeof(cv19at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at6, sizeof(cv19at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at7, sizeof(cv19at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -20 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at1, sizeof(cv20at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at2, sizeof(cv20at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at3, sizeof(cv20at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at4, sizeof(cv20at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at5, sizeof(cv20at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at6, sizeof(cv20at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at7, sizeof(cv20at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -21 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at1, sizeof(cv21at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at2, sizeof(cv21at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at3, sizeof(cv21at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at4, sizeof(cv21at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at5, sizeof(cv21at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at6, sizeof(cv21at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at7, sizeof(cv21at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -22 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at1, sizeof(cv22at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at2, sizeof(cv22at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at3, sizeof(cv22at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at4, sizeof(cv22at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at5, sizeof(cv22at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at6, sizeof(cv22at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at7, sizeof(cv22at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -23 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at1, sizeof(cv23at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at2, sizeof(cv23at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at3, sizeof(cv23at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at4, sizeof(cv23at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at5, sizeof(cv23at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at6, sizeof(cv23at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at7, sizeof(cv23at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -24 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at1, sizeof(cv24at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at2, sizeof(cv24at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at3, sizeof(cv24at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at4, sizeof(cv24at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at5, sizeof(cv24at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at6, sizeof(cv24at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at7, sizeof(cv24at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -25 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at1, sizeof(cv25at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at2, sizeof(cv25at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at3, sizeof(cv25at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at4, sizeof(cv25at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at5, sizeof(cv25at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at6, sizeof(cv25at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at7, sizeof(cv25at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -26 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at1, sizeof(cv26at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at2, sizeof(cv26at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at3, sizeof(cv26at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at4, sizeof(cv26at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at5, sizeof(cv26at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at6, sizeof(cv26at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at7, sizeof(cv26at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -27 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at1, sizeof(cv27at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at2, sizeof(cv27at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at3, sizeof(cv27at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at4, sizeof(cv27at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at5, sizeof(cv27at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at6, sizeof(cv27at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at7, sizeof(cv27at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -28 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -29 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -30 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -31 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -32 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -33 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -34 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv34at0, sizeof(cv34at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at1, sizeof(cv34at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at2, sizeof(cv34at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at3, sizeof(cv34at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -35 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv35at0, sizeof(cv35at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at1, sizeof(cv35at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -36 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at1, sizeof(cv36at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at2, sizeof(cv36at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at3, sizeof(cv36at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -37 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at1, sizeof(cv37at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -38 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at1, sizeof(cv38at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at2, sizeof(cv38at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at3, sizeof(cv38at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -39 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at1, sizeof(cv39at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -40 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -41 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -42 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -43 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv43at0, sizeof(cv43at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -44 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv44at0, sizeof(cv44at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -45 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv45at0, sizeof(cv45at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -46 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv46at0, sizeof(cv46at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at1, sizeof(cv46at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at2, sizeof(cv46at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at3, sizeof(cv46at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -47 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at1, sizeof(cv47at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -48 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at1, sizeof(cv48at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at2, sizeof(cv48at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at3, sizeof(cv48at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -49 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at1, sizeof(cv49at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -50 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv50at0, sizeof(cv50at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at1, sizeof(cv50at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at2, sizeof(cv50at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at3, sizeof(cv50at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -51 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv51at0, sizeof(cv51at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at1, sizeof(cv51at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -52 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv52at0, sizeof(cv52at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -53 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -54 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -55 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -56 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -57 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -58 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv58at0, sizeof(cv58at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at1, sizeof(cv58at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at2, sizeof(cv58at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at3, sizeof(cv58at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -59 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv59at0, sizeof(cv59at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at1, sizeof(cv59at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -60 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv60at0, sizeof(cv60at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at1, sizeof(cv60at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at2, sizeof(cv60at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at3, sizeof(cv60at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -61 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv61at0, sizeof(cv61at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at1, sizeof(cv61at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -62 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv62at0, sizeof(cv62at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at1, sizeof(cv62at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at2, sizeof(cv62at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at3, sizeof(cv62at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -63 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv63at0, sizeof(cv63at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at1, sizeof(cv63at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -64 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv64at0, sizeof(cv64at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -65 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv65at0, sizeof(cv65at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -66 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv66at0, sizeof(cv66at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -67 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv67at0, sizeof(cv67at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -68 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -69 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -70 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv70at0, sizeof(cv70at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at1, sizeof(cv70at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at2, sizeof(cv70at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at3, sizeof(cv70at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -71 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv71at0, sizeof(cv71at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at1, sizeof(cv71at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -72 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv72at0, sizeof(cv72at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at1, sizeof(cv72at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at2, sizeof(cv72at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at3, sizeof(cv72at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -73 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv73at0, sizeof(cv73at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at1, sizeof(cv73at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -74 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv74at0, sizeof(cv74at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at1, sizeof(cv74at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at2, sizeof(cv74at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at3, sizeof(cv74at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -75 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv75at0, sizeof(cv75at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at1, sizeof(cv75at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -76 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv76at0, sizeof(cv76at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at1, sizeof(cv76at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at2, sizeof(cv76at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at3, sizeof(cv76at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -77 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv77at0, sizeof(cv77at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at1, sizeof(cv77at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -78 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv78at0, sizeof(cv78at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at1, sizeof(cv78at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at2, sizeof(cv78at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at3, sizeof(cv78at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -79 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv79at0, sizeof(cv79at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at1, sizeof(cv79at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -80 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv80at0, sizeof(cv80at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at1, sizeof(cv80at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at2, sizeof(cv80at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at3, sizeof(cv80at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -81 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv81at0, sizeof(cv81at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at1, sizeof(cv81at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -82 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv82at0, sizeof(cv82at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at1, sizeof(cv82at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at2, sizeof(cv82at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at3, sizeof(cv82at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -83 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv83at0, sizeof(cv83at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at1, sizeof(cv83at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -84 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv84at0, sizeof(cv84at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at1, sizeof(cv84at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at2, sizeof(cv84at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at3, sizeof(cv84at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -85 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv85at0, sizeof(cv85at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at1, sizeof(cv85at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -86 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv86at0, sizeof(cv86at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at1, sizeof(cv86at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at2, sizeof(cv86at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at3, sizeof(cv86at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -87 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv87at0, sizeof(cv87at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at1, sizeof(cv87at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -88 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv88at0, sizeof(cv88at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at1, sizeof(cv88at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at2, sizeof(cv88at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at3, sizeof(cv88at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -89 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv89at0, sizeof(cv89at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at1, sizeof(cv89at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -90 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv90at0, sizeof(cv90at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at1, sizeof(cv90at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at2, sizeof(cv90at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at3, sizeof(cv90at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -91 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv91at0, sizeof(cv91at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at1, sizeof(cv91at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -92 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv92at0, sizeof(cv92at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at1, sizeof(cv92at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at2, sizeof(cv92at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at3, sizeof(cv92at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -93 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv93at0, sizeof(cv93at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at1, sizeof(cv93at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -94 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv94at0, sizeof(cv94at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at1, sizeof(cv94at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at2, sizeof(cv94at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at3, sizeof(cv94at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -95 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv95at0, sizeof(cv95at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at1, sizeof(cv95at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -96 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv96at0, sizeof(cv96at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at1, sizeof(cv96at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at2, sizeof(cv96at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at3, sizeof(cv96at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -97 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv97at0, sizeof(cv97at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at1, sizeof(cv97at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -98 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv98at0, sizeof(cv98at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at1, sizeof(cv98at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at2, sizeof(cv98at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at3, sizeof(cv98at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -99 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv99at0, sizeof(cv99at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at1, sizeof(cv99at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -100 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv100at0, sizeof(cv100at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at1, sizeof(cv100at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at2, sizeof(cv100at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at3, sizeof(cv100at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -101 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv101at0, sizeof(cv101at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at1, sizeof(cv101at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at2, sizeof(cv101at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at3, sizeof(cv101at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -102 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv102at0, sizeof(cv102at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at1, sizeof(cv102at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at2, sizeof(cv102at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at3, sizeof(cv102at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -103 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv103at0, sizeof(cv103at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at1, sizeof(cv103at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at2, sizeof(cv103at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at3, sizeof(cv103at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -104 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv104at0, sizeof(cv104at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at1, sizeof(cv104at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at2, sizeof(cv104at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at3, sizeof(cv104at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -105 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv105at0, sizeof(cv105at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at1, sizeof(cv105at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at2, sizeof(cv105at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at3, sizeof(cv105at3)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEqqb2gggg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEqqb2gggg

  static const ColourLines diag1[8] = { 
    ColourLines("1 3 5 9, -2 -3 -4 -6, 6 -7, 7 4 -5 -8, 8 -9"), 
    ColourLines("1 3 5 9, -2 -3 -4 -7, 6 4 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 3 5 8, -2 -3 -4 -6, 6 -7, 7 4 -5 -9, -8 9"), 
    ColourLines("1 3 5 8, -2 -3 -4 -7, 6 4 -5 -9, -6 7, -8 9"), 
    ColourLines("1 3 4 7, -2 -3 -5 -8, 6 -7, -6 -4 5 9, 8 -9"), 
    ColourLines("1 3 4 6, -2 -3 -5 -8, -6 7, -7 -4 5 9, 8 -9"), 
    ColourLines("1 3 4 7, -2 -3 -5 -9, 6 -7, -6 -4 5 8, -8 9"), 
    ColourLines("1 3 4 6, -2 -3 -5 -9, -6 7, -7 -4 5 8, -8 9")
  }; 
  static const ColourLines diag2[8] = { 
    ColourLines("1 3 5 9, -2 -3 -4 -6, 6 -8, 7 -9, -7 -5 4 8"), 
    ColourLines("1 3 5 9, -2 -3 -4 -8, 6 4 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 3 5 7, -2 -3 -4 -6, 6 -8, -7 9, 8 4 -5 -9"), 
    ColourLines("1 3 5 7, -2 -3 -4 -8, 6 4 -5 -9, -6 8, -7 9"), 
    ColourLines("1 3 4 8, -2 -3 -5 -7, 6 -8, -6 -4 5 9, 7 -9"), 
    ColourLines("1 3 4 6, -2 -3 -5 -7, -6 8, 7 -9, -8 -4 5 9"), 
    ColourLines("1 3 4 8, -2 -3 -5 -9, 6 -8, -6 -4 5 7, -7 9"), 
    ColourLines("1 3 4 6, -2 -3 -5 -9, -6 8, 7 5 -4 -8, -7 9")
  }; 
  static const ColourLines diag3[8] = { 
    ColourLines("1 3 5 8, -2 -3 -4 -6, 6 -9, 7 -8, -7 -5 4 9"), 
    ColourLines("1 3 5 8, -2 -3 -4 -9, 6 4 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 3 5 7, -2 -3 -4 -6, 6 -9, -7 8, -8 -5 4 9"), 
    ColourLines("1 3 5 7, -2 -3 -4 -9, 6 4 -5 -8, -6 9, -7 8"), 
    ColourLines("1 3 4 9, -2 -3 -5 -7, 6 -9, -6 -4 5 8, 7 -8"), 
    ColourLines("1 3 4 6, -2 -3 -5 -7, -6 9, 7 -8, 8 5 -4 -9"), 
    ColourLines("1 3 4 9, -2 -3 -5 -8, 6 -9, -6 -4 5 7, -7 8"), 
    ColourLines("1 3 4 6, -2 -3 -5 -8, -6 9, 7 5 -4 -9, -7 8")
  }; 
  static const ColourLines diag4[2] = { 
    ColourLines("1 6, -4 -7, -6 2 5 9, 7 3 -5 -8, 8 -9"), 
    ColourLines("1 6, -4 -7, -6 2 5 8, 7 3 -5 -9, -8 9")
  }; 
  static const ColourLines diag5[2] = { 
    ColourLines("1 6, -4 -8, -6 2 5 9, 7 -9, -7 -5 3 8"), 
    ColourLines("1 6, -4 -8, -6 2 5 7, -7 9, 8 3 -5 -9")
  }; 
  static const ColourLines diag6[2] = { 
    ColourLines("1 6, -4 -9, -6 2 5 8, 7 -8, -7 -5 3 9"), 
    ColourLines("1 6, -4 -9, -6 2 5 7, -7 8, -8 -5 3 9")
  }; 
  static const ColourLines diag7[2] = { 
    ColourLines("1 7, -4 -6, 6 3 -5 -8, -7 2 5 9, 8 -9"), 
    ColourLines("1 7, -4 -6, 6 3 -5 -9, -7 2 5 8, -8 9")
  }; 
  static const ColourLines diag8[2] = { 
    ColourLines("1 7, -4 -8, 6 -9, -6 -5 3 8, -7 2 5 9"), 
    ColourLines("1 7, -4 -8, 6 5 2 -7, -6 9, 8 3 -5 -9")
  }; 
  static const ColourLines diag9[2] = { 
    ColourLines("1 7, -4 -9, 6 -8, -6 -5 3 9, -7 2 5 8"), 
    ColourLines("1 7, -4 -9, 6 5 2 -7, -6 8, -8 -5 3 9")
  }; 
  static const ColourLines diag10[2] = { 
    ColourLines("1 8, -4 -6, 6 3 -5 -7, 7 -9, -8 2 5 9"), 
    ColourLines("1 8, -4 -6, 6 3 -5 -9, 7 5 2 -8, -7 9")
  }; 
  static const ColourLines diag11[2] = { 
    ColourLines("1 8, -4 -7, 6 -9, -6 -5 3 7, -8 2 5 9"), 
    ColourLines("1 8, -4 -7, 6 5 2 -8, -6 9, 7 3 -5 -9")
  }; 
  static const ColourLines diag12[2] = { 
    ColourLines("1 8, -4 -9, 6 -7, -6 -5 3 9, 7 5 2 -8"), 
    ColourLines("1 8, -4 -9, 6 5 2 -8, -6 7, -7 -5 3 9")
  }; 
  static const ColourLines diag13[2] = { 
    ColourLines("1 9, -4 -6, 6 3 -5 -7, 7 -8, 8 5 2 -9"), 
    ColourLines("1 9, -4 -6, 6 3 -5 -8, 7 5 2 -9, -7 8")
  }; 
  static const ColourLines diag14[2] = { 
    ColourLines("1 9, -4 -7, 6 -8, -6 -5 3 7, 8 5 2 -9"), 
    ColourLines("1 9, -4 -7, 6 5 2 -9, -6 8, 7 3 -5 -8")
  }; 
  static const ColourLines diag15[2] = { 
    ColourLines("1 9, -4 -8, 6 -7, -6 -5 3 8, 7 5 2 -9"), 
    ColourLines("1 9, -4 -8, 6 5 2 -9, -6 7, -7 -5 3 8")
  }; 
  static const ColourLines diag16[8] = { 
    ColourLines("1 3 8, -2 -3 -4 -9, 6 -7, -6 -5 9, 7 5 4 -8"), 
    ColourLines("1 3 8, -2 -3 -4 -9, 6 5 4 -8, -6 7, -7 -5 9"), 
    ColourLines("1 3 4 5 7, -2 -3 -8, 6 -7, -6 -5 9, 8 -4 -9"), 
    ColourLines("1 3 4 5 6, -2 -3 -8, -6 7, -7 -5 9, 8 -4 -9"), 
    ColourLines("1 3 8, -2 -3 -4 -5 -6, 6 -7, 7 5 -9, -8 4 9"), 
    ColourLines("1 3 8, -2 -3 -4 -5 -7, 6 5 -9, -6 7, -8 4 9"), 
    ColourLines("1 3 4 9, -2 -3 -8, 6 -7, -6 -5 -4 8, 7 5 -9"), 
    ColourLines("1 3 4 9, -2 -3 -8, 6 5 -9, -6 7, -7 -5 -4 8")
  }; 
  static const ColourLines diag17[8] = { 
    ColourLines("1 3 9, -2 -3 -4 -8, 6 -7, -6 -5 8, 7 5 4 -9"), 
    ColourLines("1 3 9, -2 -3 -4 -8, 6 5 4 -9, -6 7, -7 -5 8"), 
    ColourLines("1 3 4 5 7, -2 -3 -9, 6 -7, -6 -5 8, -8 -4 9"), 
    ColourLines("1 3 4 5 6, -2 -3 -9, -6 7, -7 -5 8, -8 -4 9"), 
    ColourLines("1 3 9, -2 -3 -4 -5 -6, 6 -7, 7 5 -8, 8 4 -9"), 
    ColourLines("1 3 9, -2 -3 -4 -5 -7, 6 5 -8, -6 7, 8 4 -9"), 
    ColourLines("1 3 4 8, -2 -3 -9, 6 -7, -6 -5 -4 9, 7 5 -8"), 
    ColourLines("1 3 4 8, -2 -3 -9, 6 5 -8, -6 7, -7 -5 -4 9")
  }; 
  static const ColourLines diag18[8] = { 
    ColourLines("1 3 7, -2 -3 -4 -9, 6 -8, -6 -5 9, -7 4 5 8"), 
    ColourLines("1 3 7, -2 -3 -4 -9, 6 5 4 -7, -6 8, -8 -5 9"), 
    ColourLines("1 3 4 5 8, -2 -3 -7, 6 -8, -6 -5 9, 7 -4 -9"), 
    ColourLines("1 3 4 5 6, -2 -3 -7, -6 8, 7 -4 -9, -8 -5 9"), 
    ColourLines("1 3 7, -2 -3 -4 -5 -6, 6 -8, -7 4 9, 8 5 -9"), 
    ColourLines("1 3 7, -2 -3 -4 -5 -8, 6 5 -9, -6 8, -7 4 9"), 
    ColourLines("1 3 4 9, -2 -3 -7, 6 -8, -6 -5 -4 7, 8 5 -9"), 
    ColourLines("1 3 4 9, -2 -3 -7, 6 5 -9, -6 8, 7 -4 -5 -8")
  }; 
  static const ColourLines diag19[8] = { 
    ColourLines("1 3 9, -2 -3 -4 -7, 6 -8, -6 -5 7, 8 5 4 -9"), 
    ColourLines("1 3 9, -2 -3 -4 -7, 6 5 4 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 3 4 5 8, -2 -3 -9, 6 -8, -6 -5 7, -7 -4 9"), 
    ColourLines("1 3 4 5 6, -2 -3 -9, -6 8, 7 -5 -8, -7 -4 9"), 
    ColourLines("1 3 9, -2 -3 -4 -5 -6, 6 -8, 7 4 -9, -7 5 8"), 
    ColourLines("1 3 9, -2 -3 -4 -5 -8, 6 5 -7, -6 8, 7 4 -9"), 
    ColourLines("1 3 4 7, -2 -3 -9, 6 -8, -6 -5 -4 9, -7 5 8"), 
    ColourLines("1 3 4 7, -2 -3 -9, 6 5 -7, -6 8, -8 -5 -4 9")
  }; 
  static const ColourLines diag20[8] = { 
    ColourLines("1 3 7, -2 -3 -4 -8, 6 -9, -6 -5 8, -7 4 5 9"), 
    ColourLines("1 3 7, -2 -3 -4 -8, 6 5 4 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 3 4 5 9, -2 -3 -7, 6 -9, -6 -5 8, 7 -4 -8"), 
    ColourLines("1 3 4 5 6, -2 -3 -7, -6 9, 7 -4 -8, 8 -5 -9"), 
    ColourLines("1 3 7, -2 -3 -4 -5 -6, 6 -9, -7 4 8, -8 5 9"), 
    ColourLines("1 3 7, -2 -3 -4 -5 -9, 6 5 -8, -6 9, -7 4 8"), 
    ColourLines("1 3 4 8, -2 -3 -7, 6 -9, -6 -5 -4 7, -8 5 9"), 
    ColourLines("1 3 4 8, -2 -3 -7, 6 5 -8, -6 9, 7 -4 -5 -9")
  }; 
  static const ColourLines diag21[8] = { 
    ColourLines("1 3 8, -2 -3 -4 -7, 6 -9, -6 -5 7, -8 4 5 9"), 
    ColourLines("1 3 8, -2 -3 -4 -7, 6 5 4 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 3 4 5 9, -2 -3 -8, 6 -9, -6 -5 7, -7 -4 8"), 
    ColourLines("1 3 4 5 6, -2 -3 -8, -6 9, 7 -5 -9, -7 -4 8"), 
    ColourLines("1 3 8, -2 -3 -4 -5 -6, 6 -9, 7 4 -8, -7 5 9"), 
    ColourLines("1 3 8, -2 -3 -4 -5 -9, 6 5 -7, -6 9, 7 4 -8"), 
    ColourLines("1 3 4 7, -2 -3 -8, 6 -9, -6 -5 -4 8, -7 5 9"), 
    ColourLines("1 3 4 7, -2 -3 -8, 6 5 -7, -6 9, 8 -4 -5 -9")
  }; 
  static const ColourLines diag22[8] = { 
    ColourLines("1 3 6, -2 -3 -4 -9, -6 4 5 8, 7 -8, -7 -5 9"), 
    ColourLines("1 3 6, -2 -3 -4 -9, -6 4 5 7, -7 8, -8 -5 9"), 
    ColourLines("1 3 4 5 8, -2 -3 -6, 6 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 3 4 5 7, -2 -3 -6, 6 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 3 6, -2 -3 -4 -5 -7, -6 4 9, 7 -8, 8 5 -9"), 
    ColourLines("1 3 6, -2 -3 -4 -5 -8, -6 4 9, 7 5 -9, -7 8"), 
    ColourLines("1 3 4 9, -2 -3 -6, 6 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 3 4 9, -2 -3 -6, 6 -4 -5 -8, 7 5 -9, -7 8")
  }; 
  static const ColourLines diag23[8] = { 
    ColourLines("1 3 9, -2 -3 -4 -6, 6 -5 -7, 7 -8, 8 5 4 -9"), 
    ColourLines("1 3 9, -2 -3 -4 -6, 6 -5 -8, 7 5 4 -9, -7 8"), 
    ColourLines("1 3 4 5 8, -2 -3 -9, 6 -5 -7, -6 -4 9, 7 -8"), 
    ColourLines("1 3 4 5 7, -2 -3 -9, 6 -5 -8, -6 -4 9, -7 8"), 
    ColourLines("1 3 9, -2 -3 -4 -5 -7, 6 4 -9, -6 5 8, 7 -8"), 
    ColourLines("1 3 9, -2 -3 -4 -5 -8, 6 4 -9, -6 5 7, -7 8"), 
    ColourLines("1 3 4 6, -2 -3 -9, -6 5 8, 7 -8, -7 -5 -4 9"), 
    ColourLines("1 3 4 6, -2 -3 -9, -6 5 7, -7 8, -8 -5 -4 9")
  }; 
  static const ColourLines diag24[8] = { 
    ColourLines("1 3 6, -2 -3 -4 -8, -6 4 5 9, 7 -9, -7 -5 8"), 
    ColourLines("1 3 6, -2 -3 -4 -8, -6 4 5 7, -7 9, 8 -5 -9"), 
    ColourLines("1 3 4 5 9, -2 -3 -6, 6 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 3 4 5 7, -2 -3 -6, 6 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 3 6, -2 -3 -4 -5 -7, -6 4 8, 7 -9, -8 5 9"), 
    ColourLines("1 3 6, -2 -3 -4 -5 -9, -6 4 8, 7 5 -8, -7 9"), 
    ColourLines("1 3 4 8, -2 -3 -6, 6 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 3 4 8, -2 -3 -6, 6 -4 -5 -9, 7 5 -8, -7 9")
  }; 
  static const ColourLines diag25[8] = { 
    ColourLines("1 3 8, -2 -3 -4 -6, 6 -5 -7, 7 -9, -8 4 5 9"), 
    ColourLines("1 3 8, -2 -3 -4 -6, 6 -5 -9, 7 5 4 -8, -7 9"), 
    ColourLines("1 3 4 5 9, -2 -3 -8, 6 -5 -7, -6 -4 8, 7 -9"), 
    ColourLines("1 3 4 5 7, -2 -3 -8, 6 -5 -9, -6 -4 8, -7 9"), 
    ColourLines("1 3 8, -2 -3 -4 -5 -7, 6 4 -8, -6 5 9, 7 -9"), 
    ColourLines("1 3 8, -2 -3 -4 -5 -9, 6 4 -8, -6 5 7, -7 9"), 
    ColourLines("1 3 4 6, -2 -3 -8, -6 5 9, 7 -9, -7 -5 -4 8"), 
    ColourLines("1 3 4 6, -2 -3 -8, -6 5 7, -7 9, 8 -4 -5 -9")
  }; 
  static const ColourLines diag26[8] = { 
    ColourLines("1 3 6, -2 -3 -4 -7, -6 4 5 9, 7 -5 -8, 8 -9"), 
    ColourLines("1 3 6, -2 -3 -4 -7, -6 4 5 8, 7 -5 -9, -8 9"), 
    ColourLines("1 3 4 5 9, -2 -3 -6, 6 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 3 4 5 8, -2 -3 -6, 6 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 3 6, -2 -3 -4 -5 -8, -6 4 7, -7 5 9, 8 -9"), 
    ColourLines("1 3 6, -2 -3 -4 -5 -9, -6 4 7, -7 5 8, -8 9"), 
    ColourLines("1 3 4 7, -2 -3 -6, 6 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 3 4 7, -2 -3 -6, 6 -4 -5 -9, -7 5 8, -8 9")
  }; 
  static const ColourLines diag27[8] = { 
    ColourLines("1 3 7, -2 -3 -4 -6, 6 -5 -8, -7 4 5 9, 8 -9"), 
    ColourLines("1 3 7, -2 -3 -4 -6, 6 -5 -9, -7 4 5 8, -8 9"), 
    ColourLines("1 3 4 5 9, -2 -3 -7, 6 -5 -8, -6 -4 7, 8 -9"), 
    ColourLines("1 3 4 5 8, -2 -3 -7, 6 -5 -9, -6 -4 7, -8 9"), 
    ColourLines("1 3 7, -2 -3 -4 -5 -8, 6 4 -7, -6 5 9, 8 -9"), 
    ColourLines("1 3 7, -2 -3 -4 -5 -9, 6 4 -7, -6 5 8, -8 9"), 
    ColourLines("1 3 4 6, -2 -3 -7, -6 5 9, 7 -4 -5 -8, 8 -9"), 
    ColourLines("1 3 4 6, -2 -3 -7, -6 5 8, 7 -4 -5 -9, -8 9")
  }; 
  static const ColourLines diag28[1] = { 
    ColourLines("1 6, -5 -7, -6 2 8, 7 4 -9, -8 3 9")
  }; 
  static const ColourLines diag29[1] = { 
    ColourLines("1 6, -5 -7, -6 2 9, 7 4 -8, 8 3 -9")
  }; 
  static const ColourLines diag30[1] = { 
    ColourLines("1 6, -5 -8, -6 2 7, -7 3 9, 8 4 -9")
  }; 
  static const ColourLines diag31[1] = { 
    ColourLines("1 6, -5 -8, -6 2 9, 7 3 -9, -7 4 8")
  }; 
  static const ColourLines diag32[1] = { 
    ColourLines("1 6, -5 -9, -6 2 7, -7 3 8, -8 4 9")
  }; 
  static const ColourLines diag33[1] = { 
    ColourLines("1 6, -5 -9, -6 2 8, 7 3 -8, -7 4 9")
  }; 
  static const ColourLines diag34[4] = { 
    ColourLines("1 6, -3 -4 -9, -6 2 4 5 8, 7 -8, -7 -5 9"), 
    ColourLines("1 6, -3 -4 -9, -6 2 4 5 7, -7 8, -8 -5 9"), 
    ColourLines("1 6, -3 -4 -5 -7, -6 2 4 9, 7 -8, 8 5 -9"), 
    ColourLines("1 6, -3 -4 -5 -8, -6 2 4 9, 7 5 -9, -7 8")
  }; 
  static const ColourLines diag35[2] = { 
    ColourLines("1 6, -4 -5 -7, -6 2 9, 7 -8, 8 5 3 -9"), 
    ColourLines("1 6, -4 -5 -8, -6 2 9, 7 5 3 -9, -7 8")
  }; 
  static const ColourLines diag36[4] = { 
    ColourLines("1 6, -3 -4 -8, -6 2 4 5 9, 7 -9, -7 -5 8"), 
    ColourLines("1 6, -3 -4 -8, -6 2 4 5 7, -7 9, 8 -5 -9"), 
    ColourLines("1 6, -3 -4 -5 -7, -6 2 4 8, 7 -9, -8 5 9"), 
    ColourLines("1 6, -3 -4 -5 -9, -6 2 4 8, 7 5 -8, -7 9")
  }; 
  static const ColourLines diag37[2] = { 
    ColourLines("1 6, -4 -5 -7, -6 2 8, 7 -9, -8 3 5 9"), 
    ColourLines("1 6, -4 -5 -9, -6 2 8, 7 5 3 -8, -7 9")
  }; 
  static const ColourLines diag38[4] = { 
    ColourLines("1 6, -3 -4 -7, -6 2 4 5 9, 7 -5 -8, 8 -9"), 
    ColourLines("1 6, -3 -4 -7, -6 2 4 5 8, 7 -5 -9, -8 9"), 
    ColourLines("1 6, -3 -4 -5 -8, -6 2 4 7, -7 5 9, 8 -9"), 
    ColourLines("1 6, -3 -4 -5 -9, -6 2 4 7, -7 5 8, -8 9")
  }; 
  static const ColourLines diag39[2] = { 
    ColourLines("1 6, -4 -5 -8, -6 2 7, -7 3 5 9, 8 -9"), 
    ColourLines("1 6, -4 -5 -9, -6 2 7, -7 3 5 8, -8 9")
  }; 
  static const ColourLines diag40[1] = { 
    ColourLines("1 7, -5 -6, 6 4 -9, -7 2 8, -8 3 9")
  }; 
  static const ColourLines diag41[1] = { 
    ColourLines("1 7, -5 -6, 6 4 -8, -7 2 9, 8 3 -9")
  }; 
  static const ColourLines diag42[1] = { 
    ColourLines("1 7, -5 -8, 6 2 -7, -6 3 9, 8 4 -9")
  }; 
  static const ColourLines diag43[1] = { 
    ColourLines("1 7, -5 -8, 6 3 -9, -6 4 8, -7 2 9")
  }; 
  static const ColourLines diag44[1] = { 
    ColourLines("1 7, -5 -9, 6 2 -7, -6 3 8, -8 4 9")
  }; 
  static const ColourLines diag45[1] = { 
    ColourLines("1 7, -5 -9, 6 3 -8, -6 4 9, -7 2 8")
  }; 
  static const ColourLines diag46[4] = { 
    ColourLines("1 7, -3 -4 -9, 6 -8, -6 -5 9, -7 2 4 5 8"), 
    ColourLines("1 7, -3 -4 -9, 6 5 4 2 -7, -6 8, -8 -5 9"), 
    ColourLines("1 7, -3 -4 -5 -6, 6 -8, -7 2 4 9, 8 5 -9"), 
    ColourLines("1 7, -3 -4 -5 -8, 6 5 -9, -6 8, -7 2 4 9")
  }; 
  static const ColourLines diag47[2] = { 
    ColourLines("1 7, -4 -5 -6, 6 -8, -7 2 9, 8 5 3 -9"), 
    ColourLines("1 7, -4 -5 -8, 6 5 3 -9, -6 8, -7 2 9")
  }; 
  static const ColourLines diag48[4] = { 
    ColourLines("1 7, -3 -4 -8, 6 -9, -6 -5 8, -7 2 4 5 9"), 
    ColourLines("1 7, -3 -4 -8, 6 5 4 2 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 7, -3 -4 -5 -6, 6 -9, -7 2 4 8, -8 5 9"), 
    ColourLines("1 7, -3 -4 -5 -9, 6 5 -8, -6 9, -7 2 4 8")
  }; 
  static const ColourLines diag49[2] = { 
    ColourLines("1 7, -4 -5 -6, 6 -9, -7 2 8, -8 3 5 9"), 
    ColourLines("1 7, -4 -5 -9, 6 5 3 -8, -6 9, -7 2 8")
  }; 
  static const ColourLines diag50[4] = { 
    ColourLines("1 7, -3 -4 -6, 6 -5 -8, -7 2 4 5 9, 8 -9"), 
    ColourLines("1 7, -3 -4 -6, 6 -5 -9, -7 2 4 5 8, -8 9"), 
    ColourLines("1 7, -3 -4 -5 -8, 6 4 2 -7, -6 5 9, 8 -9"), 
    ColourLines("1 7, -3 -4 -5 -9, 6 4 2 -7, -6 5 8, -8 9")
  }; 
  static const ColourLines diag51[2] = { 
    ColourLines("1 7, -4 -5 -8, 6 2 -7, -6 3 5 9, 8 -9"), 
    ColourLines("1 7, -4 -5 -9, 6 2 -7, -6 3 5 8, -8 9")
  }; 
  static const ColourLines diag52[1] = { 
    ColourLines("1 8, -5 -6, 6 4 -9, 7 2 -8, -7 3 9")
  }; 
  static const ColourLines diag53[1] = { 
    ColourLines("1 8, -5 -6, 6 4 -7, 7 3 -9, -8 2 9")
  }; 
  static const ColourLines diag54[1] = { 
    ColourLines("1 8, -5 -7, 6 2 -8, -6 3 9, 7 4 -9")
  }; 
  static const ColourLines diag55[1] = { 
    ColourLines("1 8, -5 -7, 6 3 -9, -6 4 7, -8 2 9")
  }; 
  static const ColourLines diag56[1] = { 
    ColourLines("1 8, -5 -9, 6 2 -8, -6 3 7, -7 4 9")
  }; 
  static const ColourLines diag57[1] = { 
    ColourLines("1 8, -5 -9, 6 3 -7, -6 4 9, 7 2 -8")
  }; 
  static const ColourLines diag58[4] = { 
    ColourLines("1 8, -3 -4 -9, 6 -7, -6 -5 9, 7 5 4 2 -8"), 
    ColourLines("1 8, -3 -4 -9, 6 5 4 2 -8, -6 7, -7 -5 9"), 
    ColourLines("1 8, -3 -4 -5 -6, 6 -7, 7 5 -9, -8 2 4 9"), 
    ColourLines("1 8, -3 -4 -5 -7, 6 5 -9, -6 7, -8 2 4 9")
  }; 
  static const ColourLines diag59[2] = { 
    ColourLines("1 8, -4 -5 -6, 6 -7, 7 5 3 -9, -8 2 9"), 
    ColourLines("1 8, -4 -5 -7, 6 5 3 -9, -6 7, -8 2 9")
  }; 
  static const ColourLines diag60[4] = { 
    ColourLines("1 8, -3 -4 -7, 6 -9, -6 -5 7, -8 2 4 5 9"), 
    ColourLines("1 8, -3 -4 -7, 6 5 4 2 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 8, -3 -4 -5 -6, 6 -9, 7 4 2 -8, -7 5 9"), 
    ColourLines("1 8, -3 -4 -5 -9, 6 5 -7, -6 9, 7 4 2 -8")
  }; 
  static const ColourLines diag61[2] = { 
    ColourLines("1 8, -4 -5 -6, 6 -9, 7 2 -8, -7 3 5 9"), 
    ColourLines("1 8, -4 -5 -9, 6 5 3 -7, -6 9, 7 2 -8")
  }; 
  static const ColourLines diag62[4] = { 
    ColourLines("1 8, -3 -4 -6, 6 -5 -7, 7 -9, -8 2 4 5 9"), 
    ColourLines("1 8, -3 -4 -6, 6 -5 -9, 7 5 4 2 -8, -7 9"), 
    ColourLines("1 8, -3 -4 -5 -7, 6 4 2 -8, -6 5 9, 7 -9"), 
    ColourLines("1 8, -3 -4 -5 -9, 6 4 2 -8, -6 5 7, -7 9")
  }; 
  static const ColourLines diag63[2] = { 
    ColourLines("1 8, -4 -5 -7, 6 2 -8, -6 3 5 9, 7 -9"), 
    ColourLines("1 8, -4 -5 -9, 6 2 -8, -6 3 5 7, -7 9")
  }; 
  static const ColourLines diag64[1] = { 
    ColourLines("1 9, -5 -6, 6 4 -8, 7 2 -9, -7 3 8")
  }; 
  static const ColourLines diag65[1] = { 
    ColourLines("1 9, -5 -6, 6 4 -7, 7 3 -8, 8 2 -9")
  }; 
  static const ColourLines diag66[1] = { 
    ColourLines("1 9, -5 -7, 6 2 -9, -6 3 8, 7 4 -8")
  }; 
  static const ColourLines diag67[1] = { 
    ColourLines("1 9, -5 -7, 6 3 -8, -6 4 7, 8 2 -9")
  }; 
  static const ColourLines diag68[1] = { 
    ColourLines("1 9, -5 -8, 6 2 -9, -6 3 7, -7 4 8")
  }; 
  static const ColourLines diag69[1] = { 
    ColourLines("1 9, -5 -8, 6 3 -7, -6 4 8, 7 2 -9")
  }; 
  static const ColourLines diag70[4] = { 
    ColourLines("1 9, -3 -4 -8, 6 -7, -6 -5 8, 7 5 4 2 -9"), 
    ColourLines("1 9, -3 -4 -8, 6 5 4 2 -9, -6 7, -7 -5 8"), 
    ColourLines("1 9, -3 -4 -5 -6, 6 -7, 7 5 -8, 8 4 2 -9"), 
    ColourLines("1 9, -3 -4 -5 -7, 6 5 -8, -6 7, 8 4 2 -9")
  }; 
  static const ColourLines diag71[2] = { 
    ColourLines("1 9, -4 -5 -6, 6 -7, 7 5 3 -8, 8 2 -9"), 
    ColourLines("1 9, -4 -5 -7, 6 5 3 -8, -6 7, 8 2 -9")
  }; 
  static const ColourLines diag72[4] = { 
    ColourLines("1 9, -3 -4 -7, 6 -8, -6 -5 7, 8 5 4 2 -9"), 
    ColourLines("1 9, -3 -4 -7, 6 5 4 2 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 9, -3 -4 -5 -6, 6 -8, 7 4 2 -9, -7 5 8"), 
    ColourLines("1 9, -3 -4 -5 -8, 6 5 -7, -6 8, 7 4 2 -9")
  }; 
  static const ColourLines diag73[2] = { 
    ColourLines("1 9, -4 -5 -6, 6 -8, 7 2 -9, -7 3 5 8"), 
    ColourLines("1 9, -4 -5 -8, 6 5 3 -7, -6 8, 7 2 -9")
  }; 
  static const ColourLines diag74[4] = { 
    ColourLines("1 9, -3 -4 -6, 6 -5 -7, 7 -8, 8 5 4 2 -9"), 
    ColourLines("1 9, -3 -4 -6, 6 -5 -8, 7 5 4 2 -9, -7 8"), 
    ColourLines("1 9, -3 -4 -5 -7, 6 4 2 -9, -6 5 8, 7 -8"), 
    ColourLines("1 9, -3 -4 -5 -8, 6 4 2 -9, -6 5 7, -7 8")
  }; 
  static const ColourLines diag75[2] = { 
    ColourLines("1 9, -4 -5 -7, 6 2 -9, -6 3 5 8, 7 -8"), 
    ColourLines("1 9, -4 -5 -8, 6 2 -9, -6 3 5 7, -7 8")
  }; 
  static const ColourLines diag76[4] = { 
    ColourLines("1 4 5 8, -3 -6, 6 2 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 4 5 7, -3 -6, 6 2 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 4 9, -3 -6, 6 2 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 4 9, -3 -6, 6 2 -4 -5 -8, 7 5 -9, -7 8")
  }; 
  static const ColourLines diag77[2] = { 
    ColourLines("1 5 8, -4 -6, 6 3 -9, 7 -8, -7 -5 2 9"), 
    ColourLines("1 5 7, -4 -6, 6 3 -9, -7 8, -8 -5 2 9")
  }; 
  static const ColourLines diag78[4] = { 
    ColourLines("1 4 5 9, -3 -6, 6 2 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 4 5 7, -3 -6, 6 2 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 4 8, -3 -6, 6 2 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 4 8, -3 -6, 6 2 -4 -5 -9, 7 5 -8, -7 9")
  }; 
  static const ColourLines diag79[2] = { 
    ColourLines("1 5 9, -4 -6, 6 3 -8, 7 -9, -7 -5 2 8"), 
    ColourLines("1 5 7, -4 -6, 6 3 -8, -7 9, 8 2 -5 -9")
  }; 
  static const ColourLines diag80[4] = { 
    ColourLines("1 4 5 9, -3 -6, 6 2 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 4 5 8, -3 -6, 6 2 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 4 7, -3 -6, 6 2 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 4 7, -3 -6, 6 2 -4 -5 -9, -7 5 8, -8 9")
  }; 
  static const ColourLines diag81[2] = { 
    ColourLines("1 5 9, -4 -6, 6 3 -7, 7 2 -5 -8, 8 -9"), 
    ColourLines("1 5 8, -4 -6, 6 3 -7, 7 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag82[4] = { 
    ColourLines("1 4 5 8, -3 -7, 6 -8, -6 -5 9, 7 2 -4 -9"), 
    ColourLines("1 4 5 6, -3 -7, -6 8, 7 2 -4 -9, -8 -5 9"), 
    ColourLines("1 4 9, -3 -7, 6 -8, -6 -5 -4 2 7, 8 5 -9"), 
    ColourLines("1 4 9, -3 -7, 6 5 -9, -6 8, 7 2 -4 -5 -8")
  }; 
  static const ColourLines diag83[2] = { 
    ColourLines("1 5 8, -4 -7, 6 -8, -6 -5 2 9, 7 3 -9"), 
    ColourLines("1 5 6, -4 -7, -6 8, 7 3 -9, -8 -5 2 9")
  }; 
  static const ColourLines diag84[4] = { 
    ColourLines("1 4 5 9, -3 -7, 6 -9, -6 -5 8, 7 2 -4 -8"), 
    ColourLines("1 4 5 6, -3 -7, -6 9, 7 2 -4 -8, 8 -5 -9"), 
    ColourLines("1 4 8, -3 -7, 6 -9, -6 -5 -4 2 7, -8 5 9"), 
    ColourLines("1 4 8, -3 -7, 6 5 -8, -6 9, 7 2 -4 -5 -9")
  }; 
  static const ColourLines diag85[2] = { 
    ColourLines("1 5 9, -4 -7, 6 -9, -6 -5 2 8, 7 3 -8"), 
    ColourLines("1 5 6, -4 -7, -6 9, 7 3 -8, 8 2 -5 -9")
  }; 
  static const ColourLines diag86[4] = { 
    ColourLines("1 4 5 9, -3 -7, 6 -5 -8, -6 -4 2 7, 8 -9"), 
    ColourLines("1 4 5 8, -3 -7, 6 -5 -9, -6 -4 2 7, -8 9"), 
    ColourLines("1 4 6, -3 -7, -6 5 9, 7 2 -4 -5 -8, 8 -9"), 
    ColourLines("1 4 6, -3 -7, -6 5 8, 7 2 -4 -5 -9, -8 9")
  }; 
  static const ColourLines diag87[2] = { 
    ColourLines("1 5 9, -4 -7, 6 2 -5 -8, -6 3 7, 8 -9"), 
    ColourLines("1 5 8, -4 -7, 6 2 -5 -9, -6 3 7, -8 9")
  }; 
  static const ColourLines diag88[4] = { 
    ColourLines("1 4 5 7, -3 -8, 6 -7, -6 -5 9, 8 2 -4 -9"), 
    ColourLines("1 4 5 6, -3 -8, -6 7, -7 -5 9, 8 2 -4 -9"), 
    ColourLines("1 4 9, -3 -8, 6 -7, -6 -5 -4 2 8, 7 5 -9"), 
    ColourLines("1 4 9, -3 -8, 6 5 -9, -6 7, -7 -5 -4 2 8")
  }; 
  static const ColourLines diag89[2] = { 
    ColourLines("1 5 7, -4 -8, 6 -7, -6 -5 2 9, 8 3 -9"), 
    ColourLines("1 5 6, -4 -8, -6 7, -7 -5 2 9, 8 3 -9")
  }; 
  static const ColourLines diag90[4] = { 
    ColourLines("1 4 5 9, -3 -8, 6 -9, -6 -5 7, -7 -4 2 8"), 
    ColourLines("1 4 5 6, -3 -8, -6 9, 7 -5 -9, -7 -4 2 8"), 
    ColourLines("1 4 7, -3 -8, 6 -9, -6 -5 -4 2 8, -7 5 9"), 
    ColourLines("1 4 7, -3 -8, 6 5 -7, -6 9, 8 2 -4 -5 -9")
  }; 
  static const ColourLines diag91[2] = { 
    ColourLines("1 5 9, -4 -8, 6 -9, -6 -5 2 7, -7 3 8"), 
    ColourLines("1 5 6, -4 -8, -6 9, 7 2 -5 -9, -7 3 8")
  }; 
  static const ColourLines diag92[4] = { 
    ColourLines("1 4 5 9, -3 -8, 6 -5 -7, -6 -4 2 8, 7 -9"), 
    ColourLines("1 4 5 7, -3 -8, 6 -5 -9, -6 -4 2 8, -7 9"), 
    ColourLines("1 4 6, -3 -8, -6 5 9, 7 -9, -7 -5 -4 2 8"), 
    ColourLines("1 4 6, -3 -8, -6 5 7, -7 9, 8 2 -4 -5 -9")
  }; 
  static const ColourLines diag93[2] = { 
    ColourLines("1 5 9, -4 -8, 6 2 -5 -7, -6 3 8, 7 -9"), 
    ColourLines("1 5 7, -4 -8, 6 2 -5 -9, -6 3 8, -7 9")
  }; 
  static const ColourLines diag94[4] = { 
    ColourLines("1 4 5 7, -3 -9, 6 -7, -6 -5 8, -8 -4 2 9"), 
    ColourLines("1 4 5 6, -3 -9, -6 7, -7 -5 8, -8 -4 2 9"), 
    ColourLines("1 4 8, -3 -9, 6 -7, -6 -5 -4 2 9, 7 5 -8"), 
    ColourLines("1 4 8, -3 -9, 6 5 -8, -6 7, -7 -5 -4 2 9")
  }; 
  static const ColourLines diag95[2] = { 
    ColourLines("1 5 7, -4 -9, 6 -7, -6 -5 2 8, -8 3 9"), 
    ColourLines("1 5 6, -4 -9, -6 7, -7 -5 2 8, -8 3 9")
  }; 
  static const ColourLines diag96[4] = { 
    ColourLines("1 4 5 8, -3 -9, 6 -8, -6 -5 7, -7 -4 2 9"), 
    ColourLines("1 4 5 6, -3 -9, -6 8, 7 -5 -8, -7 -4 2 9"), 
    ColourLines("1 4 7, -3 -9, 6 -8, -6 -5 -4 2 9, -7 5 8"), 
    ColourLines("1 4 7, -3 -9, 6 5 -7, -6 8, -8 -5 -4 2 9")
  }; 
  static const ColourLines diag97[2] = { 
    ColourLines("1 5 8, -4 -9, 6 -8, -6 -5 2 7, -7 3 9"), 
    ColourLines("1 5 6, -4 -9, -6 8, 7 2 -5 -8, -7 3 9")
  }; 
  static const ColourLines diag98[4] = { 
    ColourLines("1 4 5 8, -3 -9, 6 -5 -7, -6 -4 2 9, 7 -8"), 
    ColourLines("1 4 5 7, -3 -9, 6 -5 -8, -6 -4 2 9, -7 8"), 
    ColourLines("1 4 6, -3 -9, -6 5 8, 7 -8, -7 -5 -4 2 9"), 
    ColourLines("1 4 6, -3 -9, -6 5 7, -7 8, -8 -5 -4 2 9")
  }; 
  static const ColourLines diag99[2] = { 
    ColourLines("1 5 8, -4 -9, 6 2 -5 -7, -6 3 9, 7 -8"), 
    ColourLines("1 5 7, -4 -9, 6 2 -5 -8, -6 3 9, -7 8")
  }; 
  static const ColourLines diag100[4] = { 
    ColourLines("1 4 7, -3 -5 -8, 6 -7, -6 -4 2 5 9, 8 -9"), 
    ColourLines("1 4 6, -3 -5 -8, -6 7, -7 -4 2 5 9, 8 -9"), 
    ColourLines("1 4 7, -3 -5 -9, 6 -7, -6 -4 2 5 8, -8 9"), 
    ColourLines("1 4 6, -3 -5 -9, -6 7, -7 -4 2 5 8, -8 9")
  }; 
  static const ColourLines diag101[4] = { 
    ColourLines("1 5 9, -3 -4 -6, 6 -7, 7 4 2 -5 -8, 8 -9"), 
    ColourLines("1 5 9, -3 -4 -7, 6 4 2 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 5 8, -3 -4 -6, 6 -7, 7 4 2 -5 -9, -8 9"), 
    ColourLines("1 5 8, -3 -4 -7, 6 4 2 -5 -9, -6 7, -8 9")
  }; 
  static const ColourLines diag102[4] = { 
    ColourLines("1 4 8, -3 -5 -7, 6 -8, -6 -4 2 5 9, 7 -9"), 
    ColourLines("1 4 6, -3 -5 -7, -6 8, 7 -9, -8 -4 2 5 9"), 
    ColourLines("1 4 8, -3 -5 -9, 6 -8, -6 -4 2 5 7, -7 9"), 
    ColourLines("1 4 6, -3 -5 -9, -6 8, 7 5 2 -4 -8, -7 9")
  }; 
  static const ColourLines diag103[4] = { 
    ColourLines("1 5 9, -3 -4 -6, 6 -8, 7 -9, -7 -5 2 4 8"), 
    ColourLines("1 5 9, -3 -4 -8, 6 4 2 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 5 7, -3 -4 -6, 6 -8, -7 9, 8 4 2 -5 -9"), 
    ColourLines("1 5 7, -3 -4 -8, 6 4 2 -5 -9, -6 8, -7 9")
  }; 
  static const ColourLines diag104[4] = { 
    ColourLines("1 4 9, -3 -5 -7, 6 -9, -6 -4 2 5 8, 7 -8"), 
    ColourLines("1 4 6, -3 -5 -7, -6 9, 7 -8, 8 5 2 -4 -9"), 
    ColourLines("1 4 9, -3 -5 -8, 6 -9, -6 -4 2 5 7, -7 8"), 
    ColourLines("1 4 6, -3 -5 -8, -6 9, 7 5 2 -4 -9, -7 8")
  }; 
  static const ColourLines diag105[4] = { 
    ColourLines("1 5 8, -3 -4 -6, 6 -9, 7 -8, -7 -5 2 4 9"), 
    ColourLines("1 5 8, -3 -4 -9, 6 4 2 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 5 7, -3 -4 -6, 6 -9, -7 8, -8 -5 2 4 9"), 
    ColourLines("1 5 7, -3 -4 -9, 6 4 2 -5 -8, -6 9, -7 8")
  }; 

  static const int cv1at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv1at1[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv1at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv1at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv1at4[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv1at5[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv1at6[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv1at7[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv2at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv2at1[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv2at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv2at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv2at4[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv2at5[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv2at6[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv2at7[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv3at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv3at1[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv3at2[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv3at3[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv3at4[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv3at5[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv3at6[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv3at7[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv4at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv4at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv5at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv5at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv6at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv6at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv7at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv7at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv8at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv8at1[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv9at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv9at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv10at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv10at1[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv11at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv11at1[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv12at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv12at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv13at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv13at1[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv14at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv14at1[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv15at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv15at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv16at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv16at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv16at2[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv16at3[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv16at4[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv16at5[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv16at6[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv16at7[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv17at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv17at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv17at2[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv17at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv17at4[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv17at5[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv17at6[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv17at7[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv18at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv18at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv18at2[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv18at3[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv18at4[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv18at5[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv18at6[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv18at7[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv19at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv19at1[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv19at2[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv19at3[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv19at4[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv19at5[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv19at6[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv19at7[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv20at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv20at1[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv20at2[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv20at3[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv20at4[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv20at5[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv20at6[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv20at7[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv21at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv21at1[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv21at2[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv21at3[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv21at4[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv21at5[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv21at6[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv21at7[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv22at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv22at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv22at2[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv22at3[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv22at4[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv22at5[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv22at6[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv22at7[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv23at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv23at1[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv23at2[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv23at3[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv23at4[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv23at5[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv23at6[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv23at7[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv24at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv24at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv24at2[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv24at3[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv24at4[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv24at5[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv24at6[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv24at7[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv25at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv25at1[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv25at2[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv25at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv25at4[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv25at5[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv25at6[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv25at7[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv26at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv26at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv26at2[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv26at3[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv26at4[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv26at5[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv26at6[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv26at7[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv27at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv27at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv27at2[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv27at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv27at4[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv27at5[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv27at6[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv27at7[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv28at0[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv29at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv30at0[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv31at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv32at0[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv33at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv34at0[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv34at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv34at2[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv34at3[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv35at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv35at1[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv36at0[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv36at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv36at2[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv36at3[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv37at0[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv37at1[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv38at0[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv38at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv38at2[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv38at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv39at0[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv39at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv40at0[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv41at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv42at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv43at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv44at0[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv45at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv46at0[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv46at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv46at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv46at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv47at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv47at1[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv48at0[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv48at1[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv48at2[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv48at3[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv49at0[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv49at1[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv50at0[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv50at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv50at2[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv50at3[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv51at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv51at1[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv52at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv53at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv54at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv55at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv56at0[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv57at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv58at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv58at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv58at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv58at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv59at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv59at1[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv60at0[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv60at1[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv60at2[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv60at3[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv61at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv61at1[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv62at0[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv62at1[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv62at2[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv62at3[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv63at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv63at1[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv64at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv65at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv66at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv67at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv68at0[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv69at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv70at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv70at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv70at2[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv70at3[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv71at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv71at1[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv72at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv72at1[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv72at2[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv72at3[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv73at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv73at1[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv74at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv74at1[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv74at2[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv74at3[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv75at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv75at1[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv76at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv76at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv76at2[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv76at3[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv77at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv77at1[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv78at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv78at1[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv78at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv78at3[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv79at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv79at1[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv80at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv80at1[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv80at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv80at3[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv81at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv81at1[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv82at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv82at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv82at2[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv82at3[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv83at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv83at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv84at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv84at1[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv84at2[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv84at3[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv85at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv85at1[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv86at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv86at1[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv86at2[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv86at3[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv87at0[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv87at1[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv88at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv88at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv88at2[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv88at3[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv89at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv89at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv90at0[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv90at1[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv90at2[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv90at3[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv91at0[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv91at1[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv92at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv92at1[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv92at2[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv92at3[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv93at0[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv93at1[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv94at0[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv94at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv94at2[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv94at3[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv95at0[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv95at1[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv96at0[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv96at1[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv96at2[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv96at3[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv97at0[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv97at1[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv98at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv98at1[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv98at2[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv98at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv99at0[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv99at1[]  = { 0, 4, 1, 3, 2, -1, -999};
  static const int cv100at0[]  = { 0, 3, 4, 1, 2, -1, -999};
  static const int cv100at1[]  = { 0, 3, 4, 2, 1, -1, -999};
  static const int cv100at2[]  = { 0, 4, 3, 1, 2, -1, -999};
  static const int cv100at3[]  = { 0, 4, 3, 2, 1, -1, -999};
  static const int cv101at0[]  = { 0, 1, 2, 3, 4, -1, -999};
  static const int cv101at1[]  = { 0, 2, 1, 3, 4, -1, -999};
  static const int cv101at2[]  = { 0, 1, 2, 4, 3, -1, -999};
  static const int cv101at3[]  = { 0, 2, 1, 4, 3, -1, -999};
  static const int cv102at0[]  = { 0, 2, 4, 1, 3, -1, -999};
  static const int cv102at1[]  = { 0, 2, 4, 3, 1, -1, -999};
  static const int cv102at2[]  = { 0, 4, 2, 1, 3, -1, -999};
  static const int cv102at3[]  = { 0, 4, 2, 3, 1, -1, -999};
  static const int cv103at0[]  = { 0, 1, 3, 2, 4, -1, -999};
  static const int cv103at1[]  = { 0, 3, 1, 2, 4, -1, -999};
  static const int cv103at2[]  = { 0, 1, 3, 4, 2, -1, -999};
  static const int cv103at3[]  = { 0, 3, 1, 4, 2, -1, -999};
  static const int cv104at0[]  = { 0, 2, 3, 1, 4, -1, -999};
  static const int cv104at1[]  = { 0, 2, 3, 4, 1, -1, -999};
  static const int cv104at2[]  = { 0, 3, 2, 1, 4, -1, -999};
  static const int cv104at3[]  = { 0, 3, 2, 4, 1, -1, -999};
  static const int cv105at0[]  = { 0, 1, 4, 2, 3, -1, -999};
  static const int cv105at1[]  = { 0, 4, 1, 2, 3, -1, -999};
  static const int cv105at2[]  = { 0, 1, 4, 3, 2, -1, -999};
  static const int cv105at3[]  = { 0, 4, 1, 3, 2, -1, -999};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int)),  &(diag1[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at2, sizeof(cv1at2)/sizeof(int)),  &(diag1[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at3, sizeof(cv1at3)/sizeof(int)),  &(diag1[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at4, sizeof(cv1at4)/sizeof(int)),  &(diag1[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at5, sizeof(cv1at5)/sizeof(int)),  &(diag1[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at6, sizeof(cv1at6)/sizeof(int)),  &(diag1[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at7, sizeof(cv1at7)/sizeof(int)),  &(diag1[7]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int)),  &(diag2[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at2, sizeof(cv2at2)/sizeof(int)),  &(diag2[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at3, sizeof(cv2at3)/sizeof(int)),  &(diag2[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at4, sizeof(cv2at4)/sizeof(int)),  &(diag2[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at5, sizeof(cv2at5)/sizeof(int)),  &(diag2[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at6, sizeof(cv2at6)/sizeof(int)),  &(diag2[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at7, sizeof(cv2at7)/sizeof(int)),  &(diag2[7]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int)),  &(diag3[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at2, sizeof(cv3at2)/sizeof(int)),  &(diag3[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at3, sizeof(cv3at3)/sizeof(int)),  &(diag3[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at4, sizeof(cv3at4)/sizeof(int)),  &(diag3[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at5, sizeof(cv3at5)/sizeof(int)),  &(diag3[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at6, sizeof(cv3at6)/sizeof(int)),  &(diag3[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at7, sizeof(cv3at7)/sizeof(int)),  &(diag3[7]) );
  } 
  else if( diag->id() == -4 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int)),  &(diag4[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int)),  &(diag4[1]) );
  } 
  else if( diag->id() == -5 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int)),  &(diag5[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int)),  &(diag5[1]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int)),  &(diag6[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int)),  &(diag6[1]) );
  } 
  else if( diag->id() == -7 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int)),  &(diag7[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at1, sizeof(cv7at1)/sizeof(int)),  &(diag7[1]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int)),  &(diag8[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int)),  &(diag8[1]) );
  } 
  else if( diag->id() == -9 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int)),  &(diag9[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int)),  &(diag9[1]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)),  &(diag10[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int)),  &(diag10[1]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)),  &(diag11[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int)),  &(diag11[1]) );
  } 
  else if( diag->id() == -12 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int)),  &(diag12[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at1, sizeof(cv12at1)/sizeof(int)),  &(diag12[1]) );
  } 
  else if( diag->id() == -13 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int)),  &(diag13[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at1, sizeof(cv13at1)/sizeof(int)),  &(diag13[1]) );
  } 
  else if( diag->id() == -14 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int)),  &(diag14[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at1, sizeof(cv14at1)/sizeof(int)),  &(diag14[1]) );
  } 
  else if( diag->id() == -15 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int)),  &(diag15[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int)),  &(diag15[1]) );
  } 
  else if( diag->id() == -16 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at0, sizeof(cv16at0)/sizeof(int)),  &(diag16[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at1, sizeof(cv16at1)/sizeof(int)),  &(diag16[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at2, sizeof(cv16at2)/sizeof(int)),  &(diag16[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at3, sizeof(cv16at3)/sizeof(int)),  &(diag16[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at4, sizeof(cv16at4)/sizeof(int)),  &(diag16[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at5, sizeof(cv16at5)/sizeof(int)),  &(diag16[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at6, sizeof(cv16at6)/sizeof(int)),  &(diag16[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at7, sizeof(cv16at7)/sizeof(int)),  &(diag16[7]) );
  } 
  else if( diag->id() == -17 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at0, sizeof(cv17at0)/sizeof(int)),  &(diag17[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at1, sizeof(cv17at1)/sizeof(int)),  &(diag17[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at2, sizeof(cv17at2)/sizeof(int)),  &(diag17[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at3, sizeof(cv17at3)/sizeof(int)),  &(diag17[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at4, sizeof(cv17at4)/sizeof(int)),  &(diag17[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at5, sizeof(cv17at5)/sizeof(int)),  &(diag17[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at6, sizeof(cv17at6)/sizeof(int)),  &(diag17[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at7, sizeof(cv17at7)/sizeof(int)),  &(diag17[7]) );
  } 
  else if( diag->id() == -18 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at0, sizeof(cv18at0)/sizeof(int)),  &(diag18[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at1, sizeof(cv18at1)/sizeof(int)),  &(diag18[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at2, sizeof(cv18at2)/sizeof(int)),  &(diag18[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at3, sizeof(cv18at3)/sizeof(int)),  &(diag18[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at4, sizeof(cv18at4)/sizeof(int)),  &(diag18[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at5, sizeof(cv18at5)/sizeof(int)),  &(diag18[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at6, sizeof(cv18at6)/sizeof(int)),  &(diag18[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at7, sizeof(cv18at7)/sizeof(int)),  &(diag18[7]) );
  } 
  else if( diag->id() == -19 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int)),  &(diag19[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at1, sizeof(cv19at1)/sizeof(int)),  &(diag19[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at2, sizeof(cv19at2)/sizeof(int)),  &(diag19[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at3, sizeof(cv19at3)/sizeof(int)),  &(diag19[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at4, sizeof(cv19at4)/sizeof(int)),  &(diag19[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at5, sizeof(cv19at5)/sizeof(int)),  &(diag19[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at6, sizeof(cv19at6)/sizeof(int)),  &(diag19[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at7, sizeof(cv19at7)/sizeof(int)),  &(diag19[7]) );
  } 
  else if( diag->id() == -20 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int)),  &(diag20[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at1, sizeof(cv20at1)/sizeof(int)),  &(diag20[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at2, sizeof(cv20at2)/sizeof(int)),  &(diag20[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at3, sizeof(cv20at3)/sizeof(int)),  &(diag20[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at4, sizeof(cv20at4)/sizeof(int)),  &(diag20[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at5, sizeof(cv20at5)/sizeof(int)),  &(diag20[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at6, sizeof(cv20at6)/sizeof(int)),  &(diag20[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at7, sizeof(cv20at7)/sizeof(int)),  &(diag20[7]) );
  } 
  else if( diag->id() == -21 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int)),  &(diag21[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at1, sizeof(cv21at1)/sizeof(int)),  &(diag21[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at2, sizeof(cv21at2)/sizeof(int)),  &(diag21[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at3, sizeof(cv21at3)/sizeof(int)),  &(diag21[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at4, sizeof(cv21at4)/sizeof(int)),  &(diag21[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at5, sizeof(cv21at5)/sizeof(int)),  &(diag21[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at6, sizeof(cv21at6)/sizeof(int)),  &(diag21[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at7, sizeof(cv21at7)/sizeof(int)),  &(diag21[7]) );
  } 
  else if( diag->id() == -22 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int)),  &(diag22[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at1, sizeof(cv22at1)/sizeof(int)),  &(diag22[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at2, sizeof(cv22at2)/sizeof(int)),  &(diag22[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at3, sizeof(cv22at3)/sizeof(int)),  &(diag22[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at4, sizeof(cv22at4)/sizeof(int)),  &(diag22[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at5, sizeof(cv22at5)/sizeof(int)),  &(diag22[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at6, sizeof(cv22at6)/sizeof(int)),  &(diag22[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at7, sizeof(cv22at7)/sizeof(int)),  &(diag22[7]) );
  } 
  else if( diag->id() == -23 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int)),  &(diag23[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at1, sizeof(cv23at1)/sizeof(int)),  &(diag23[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at2, sizeof(cv23at2)/sizeof(int)),  &(diag23[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at3, sizeof(cv23at3)/sizeof(int)),  &(diag23[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at4, sizeof(cv23at4)/sizeof(int)),  &(diag23[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at5, sizeof(cv23at5)/sizeof(int)),  &(diag23[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at6, sizeof(cv23at6)/sizeof(int)),  &(diag23[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at7, sizeof(cv23at7)/sizeof(int)),  &(diag23[7]) );
  } 
  else if( diag->id() == -24 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int)),  &(diag24[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at1, sizeof(cv24at1)/sizeof(int)),  &(diag24[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at2, sizeof(cv24at2)/sizeof(int)),  &(diag24[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at3, sizeof(cv24at3)/sizeof(int)),  &(diag24[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at4, sizeof(cv24at4)/sizeof(int)),  &(diag24[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at5, sizeof(cv24at5)/sizeof(int)),  &(diag24[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at6, sizeof(cv24at6)/sizeof(int)),  &(diag24[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at7, sizeof(cv24at7)/sizeof(int)),  &(diag24[7]) );
  } 
  else if( diag->id() == -25 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int)),  &(diag25[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at1, sizeof(cv25at1)/sizeof(int)),  &(diag25[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at2, sizeof(cv25at2)/sizeof(int)),  &(diag25[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at3, sizeof(cv25at3)/sizeof(int)),  &(diag25[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at4, sizeof(cv25at4)/sizeof(int)),  &(diag25[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at5, sizeof(cv25at5)/sizeof(int)),  &(diag25[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at6, sizeof(cv25at6)/sizeof(int)),  &(diag25[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at7, sizeof(cv25at7)/sizeof(int)),  &(diag25[7]) );
  } 
  else if( diag->id() == -26 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int)),  &(diag26[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at1, sizeof(cv26at1)/sizeof(int)),  &(diag26[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at2, sizeof(cv26at2)/sizeof(int)),  &(diag26[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at3, sizeof(cv26at3)/sizeof(int)),  &(diag26[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at4, sizeof(cv26at4)/sizeof(int)),  &(diag26[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at5, sizeof(cv26at5)/sizeof(int)),  &(diag26[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at6, sizeof(cv26at6)/sizeof(int)),  &(diag26[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at7, sizeof(cv26at7)/sizeof(int)),  &(diag26[7]) );
  } 
  else if( diag->id() == -27 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int)),  &(diag27[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at1, sizeof(cv27at1)/sizeof(int)),  &(diag27[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at2, sizeof(cv27at2)/sizeof(int)),  &(diag27[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at3, sizeof(cv27at3)/sizeof(int)),  &(diag27[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at4, sizeof(cv27at4)/sizeof(int)),  &(diag27[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at5, sizeof(cv27at5)/sizeof(int)),  &(diag27[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at6, sizeof(cv27at6)/sizeof(int)),  &(diag27[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at7, sizeof(cv27at7)/sizeof(int)),  &(diag27[7]) );
  } 
  else if( diag->id() == -28 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)),  &(diag28[0]) );
  } 
  else if( diag->id() == -29 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int)),  &(diag29[0]) );
  } 
  else if( diag->id() == -30 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int)),  &(diag30[0]) );
  } 
  else if( diag->id() == -31 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int)),  &(diag31[0]) );
  } 
  else if( diag->id() == -32 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int)),  &(diag32[0]) );
  } 
  else if( diag->id() == -33 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int)),  &(diag33[0]) );
  } 
  else if( diag->id() == -34 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at0, sizeof(cv34at0)/sizeof(int)),  &(diag34[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at1, sizeof(cv34at1)/sizeof(int)),  &(diag34[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at2, sizeof(cv34at2)/sizeof(int)),  &(diag34[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at3, sizeof(cv34at3)/sizeof(int)),  &(diag34[3]) );
  } 
  else if( diag->id() == -35 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at0, sizeof(cv35at0)/sizeof(int)),  &(diag35[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at1, sizeof(cv35at1)/sizeof(int)),  &(diag35[1]) );
  } 
  else if( diag->id() == -36 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int)),  &(diag36[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at1, sizeof(cv36at1)/sizeof(int)),  &(diag36[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at2, sizeof(cv36at2)/sizeof(int)),  &(diag36[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at3, sizeof(cv36at3)/sizeof(int)),  &(diag36[3]) );
  } 
  else if( diag->id() == -37 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int)),  &(diag37[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at1, sizeof(cv37at1)/sizeof(int)),  &(diag37[1]) );
  } 
  else if( diag->id() == -38 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int)),  &(diag38[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at1, sizeof(cv38at1)/sizeof(int)),  &(diag38[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at2, sizeof(cv38at2)/sizeof(int)),  &(diag38[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at3, sizeof(cv38at3)/sizeof(int)),  &(diag38[3]) );
  } 
  else if( diag->id() == -39 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int)),  &(diag39[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at1, sizeof(cv39at1)/sizeof(int)),  &(diag39[1]) );
  } 
  else if( diag->id() == -40 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int)),  &(diag40[0]) );
  } 
  else if( diag->id() == -41 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int)),  &(diag41[0]) );
  } 
  else if( diag->id() == -42 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int)),  &(diag42[0]) );
  } 
  else if( diag->id() == -43 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at0, sizeof(cv43at0)/sizeof(int)),  &(diag43[0]) );
  } 
  else if( diag->id() == -44 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at0, sizeof(cv44at0)/sizeof(int)),  &(diag44[0]) );
  } 
  else if( diag->id() == -45 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at0, sizeof(cv45at0)/sizeof(int)),  &(diag45[0]) );
  } 
  else if( diag->id() == -46 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at0, sizeof(cv46at0)/sizeof(int)),  &(diag46[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at1, sizeof(cv46at1)/sizeof(int)),  &(diag46[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at2, sizeof(cv46at2)/sizeof(int)),  &(diag46[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at3, sizeof(cv46at3)/sizeof(int)),  &(diag46[3]) );
  } 
  else if( diag->id() == -47 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int)),  &(diag47[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at1, sizeof(cv47at1)/sizeof(int)),  &(diag47[1]) );
  } 
  else if( diag->id() == -48 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int)),  &(diag48[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at1, sizeof(cv48at1)/sizeof(int)),  &(diag48[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at2, sizeof(cv48at2)/sizeof(int)),  &(diag48[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at3, sizeof(cv48at3)/sizeof(int)),  &(diag48[3]) );
  } 
  else if( diag->id() == -49 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int)),  &(diag49[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at1, sizeof(cv49at1)/sizeof(int)),  &(diag49[1]) );
  } 
  else if( diag->id() == -50 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at0, sizeof(cv50at0)/sizeof(int)),  &(diag50[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at1, sizeof(cv50at1)/sizeof(int)),  &(diag50[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at2, sizeof(cv50at2)/sizeof(int)),  &(diag50[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at3, sizeof(cv50at3)/sizeof(int)),  &(diag50[3]) );
  } 
  else if( diag->id() == -51 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at0, sizeof(cv51at0)/sizeof(int)),  &(diag51[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at1, sizeof(cv51at1)/sizeof(int)),  &(diag51[1]) );
  } 
  else if( diag->id() == -52 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at0, sizeof(cv52at0)/sizeof(int)),  &(diag52[0]) );
  } 
  else if( diag->id() == -53 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int)),  &(diag53[0]) );
  } 
  else if( diag->id() == -54 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int)),  &(diag54[0]) );
  } 
  else if( diag->id() == -55 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int)),  &(diag55[0]) );
  } 
  else if( diag->id() == -56 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int)),  &(diag56[0]) );
  } 
  else if( diag->id() == -57 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int)),  &(diag57[0]) );
  } 
  else if( diag->id() == -58 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at0, sizeof(cv58at0)/sizeof(int)),  &(diag58[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at1, sizeof(cv58at1)/sizeof(int)),  &(diag58[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at2, sizeof(cv58at2)/sizeof(int)),  &(diag58[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at3, sizeof(cv58at3)/sizeof(int)),  &(diag58[3]) );
  } 
  else if( diag->id() == -59 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at0, sizeof(cv59at0)/sizeof(int)),  &(diag59[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at1, sizeof(cv59at1)/sizeof(int)),  &(diag59[1]) );
  } 
  else if( diag->id() == -60 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at0, sizeof(cv60at0)/sizeof(int)),  &(diag60[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at1, sizeof(cv60at1)/sizeof(int)),  &(diag60[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at2, sizeof(cv60at2)/sizeof(int)),  &(diag60[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at3, sizeof(cv60at3)/sizeof(int)),  &(diag60[3]) );
  } 
  else if( diag->id() == -61 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at0, sizeof(cv61at0)/sizeof(int)),  &(diag61[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at1, sizeof(cv61at1)/sizeof(int)),  &(diag61[1]) );
  } 
  else if( diag->id() == -62 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at0, sizeof(cv62at0)/sizeof(int)),  &(diag62[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at1, sizeof(cv62at1)/sizeof(int)),  &(diag62[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at2, sizeof(cv62at2)/sizeof(int)),  &(diag62[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at3, sizeof(cv62at3)/sizeof(int)),  &(diag62[3]) );
  } 
  else if( diag->id() == -63 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at0, sizeof(cv63at0)/sizeof(int)),  &(diag63[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at1, sizeof(cv63at1)/sizeof(int)),  &(diag63[1]) );
  } 
  else if( diag->id() == -64 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at0, sizeof(cv64at0)/sizeof(int)),  &(diag64[0]) );
  } 
  else if( diag->id() == -65 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at0, sizeof(cv65at0)/sizeof(int)),  &(diag65[0]) );
  } 
  else if( diag->id() == -66 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at0, sizeof(cv66at0)/sizeof(int)),  &(diag66[0]) );
  } 
  else if( diag->id() == -67 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at0, sizeof(cv67at0)/sizeof(int)),  &(diag67[0]) );
  } 
  else if( diag->id() == -68 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int)),  &(diag68[0]) );
  } 
  else if( diag->id() == -69 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int)),  &(diag69[0]) );
  } 
  else if( diag->id() == -70 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at0, sizeof(cv70at0)/sizeof(int)),  &(diag70[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at1, sizeof(cv70at1)/sizeof(int)),  &(diag70[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at2, sizeof(cv70at2)/sizeof(int)),  &(diag70[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at3, sizeof(cv70at3)/sizeof(int)),  &(diag70[3]) );
  } 
  else if( diag->id() == -71 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at0, sizeof(cv71at0)/sizeof(int)),  &(diag71[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at1, sizeof(cv71at1)/sizeof(int)),  &(diag71[1]) );
  } 
  else if( diag->id() == -72 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at0, sizeof(cv72at0)/sizeof(int)),  &(diag72[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at1, sizeof(cv72at1)/sizeof(int)),  &(diag72[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at2, sizeof(cv72at2)/sizeof(int)),  &(diag72[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at3, sizeof(cv72at3)/sizeof(int)),  &(diag72[3]) );
  } 
  else if( diag->id() == -73 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at0, sizeof(cv73at0)/sizeof(int)),  &(diag73[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at1, sizeof(cv73at1)/sizeof(int)),  &(diag73[1]) );
  } 
  else if( diag->id() == -74 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at0, sizeof(cv74at0)/sizeof(int)),  &(diag74[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at1, sizeof(cv74at1)/sizeof(int)),  &(diag74[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at2, sizeof(cv74at2)/sizeof(int)),  &(diag74[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at3, sizeof(cv74at3)/sizeof(int)),  &(diag74[3]) );
  } 
  else if( diag->id() == -75 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at0, sizeof(cv75at0)/sizeof(int)),  &(diag75[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at1, sizeof(cv75at1)/sizeof(int)),  &(diag75[1]) );
  } 
  else if( diag->id() == -76 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at0, sizeof(cv76at0)/sizeof(int)),  &(diag76[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at1, sizeof(cv76at1)/sizeof(int)),  &(diag76[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at2, sizeof(cv76at2)/sizeof(int)),  &(diag76[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at3, sizeof(cv76at3)/sizeof(int)),  &(diag76[3]) );
  } 
  else if( diag->id() == -77 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at0, sizeof(cv77at0)/sizeof(int)),  &(diag77[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at1, sizeof(cv77at1)/sizeof(int)),  &(diag77[1]) );
  } 
  else if( diag->id() == -78 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at0, sizeof(cv78at0)/sizeof(int)),  &(diag78[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at1, sizeof(cv78at1)/sizeof(int)),  &(diag78[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at2, sizeof(cv78at2)/sizeof(int)),  &(diag78[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at3, sizeof(cv78at3)/sizeof(int)),  &(diag78[3]) );
  } 
  else if( diag->id() == -79 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at0, sizeof(cv79at0)/sizeof(int)),  &(diag79[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at1, sizeof(cv79at1)/sizeof(int)),  &(diag79[1]) );
  } 
  else if( diag->id() == -80 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at0, sizeof(cv80at0)/sizeof(int)),  &(diag80[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at1, sizeof(cv80at1)/sizeof(int)),  &(diag80[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at2, sizeof(cv80at2)/sizeof(int)),  &(diag80[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at3, sizeof(cv80at3)/sizeof(int)),  &(diag80[3]) );
  } 
  else if( diag->id() == -81 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at0, sizeof(cv81at0)/sizeof(int)),  &(diag81[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at1, sizeof(cv81at1)/sizeof(int)),  &(diag81[1]) );
  } 
  else if( diag->id() == -82 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at0, sizeof(cv82at0)/sizeof(int)),  &(diag82[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at1, sizeof(cv82at1)/sizeof(int)),  &(diag82[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at2, sizeof(cv82at2)/sizeof(int)),  &(diag82[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at3, sizeof(cv82at3)/sizeof(int)),  &(diag82[3]) );
  } 
  else if( diag->id() == -83 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at0, sizeof(cv83at0)/sizeof(int)),  &(diag83[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at1, sizeof(cv83at1)/sizeof(int)),  &(diag83[1]) );
  } 
  else if( diag->id() == -84 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at0, sizeof(cv84at0)/sizeof(int)),  &(diag84[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at1, sizeof(cv84at1)/sizeof(int)),  &(diag84[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at2, sizeof(cv84at2)/sizeof(int)),  &(diag84[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at3, sizeof(cv84at3)/sizeof(int)),  &(diag84[3]) );
  } 
  else if( diag->id() == -85 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at0, sizeof(cv85at0)/sizeof(int)),  &(diag85[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at1, sizeof(cv85at1)/sizeof(int)),  &(diag85[1]) );
  } 
  else if( diag->id() == -86 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at0, sizeof(cv86at0)/sizeof(int)),  &(diag86[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at1, sizeof(cv86at1)/sizeof(int)),  &(diag86[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at2, sizeof(cv86at2)/sizeof(int)),  &(diag86[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at3, sizeof(cv86at3)/sizeof(int)),  &(diag86[3]) );
  } 
  else if( diag->id() == -87 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at0, sizeof(cv87at0)/sizeof(int)),  &(diag87[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at1, sizeof(cv87at1)/sizeof(int)),  &(diag87[1]) );
  } 
  else if( diag->id() == -88 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at0, sizeof(cv88at0)/sizeof(int)),  &(diag88[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at1, sizeof(cv88at1)/sizeof(int)),  &(diag88[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at2, sizeof(cv88at2)/sizeof(int)),  &(diag88[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at3, sizeof(cv88at3)/sizeof(int)),  &(diag88[3]) );
  } 
  else if( diag->id() == -89 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at0, sizeof(cv89at0)/sizeof(int)),  &(diag89[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at1, sizeof(cv89at1)/sizeof(int)),  &(diag89[1]) );
  } 
  else if( diag->id() == -90 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at0, sizeof(cv90at0)/sizeof(int)),  &(diag90[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at1, sizeof(cv90at1)/sizeof(int)),  &(diag90[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at2, sizeof(cv90at2)/sizeof(int)),  &(diag90[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at3, sizeof(cv90at3)/sizeof(int)),  &(diag90[3]) );
  } 
  else if( diag->id() == -91 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at0, sizeof(cv91at0)/sizeof(int)),  &(diag91[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at1, sizeof(cv91at1)/sizeof(int)),  &(diag91[1]) );
  } 
  else if( diag->id() == -92 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at0, sizeof(cv92at0)/sizeof(int)),  &(diag92[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at1, sizeof(cv92at1)/sizeof(int)),  &(diag92[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at2, sizeof(cv92at2)/sizeof(int)),  &(diag92[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at3, sizeof(cv92at3)/sizeof(int)),  &(diag92[3]) );
  } 
  else if( diag->id() == -93 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at0, sizeof(cv93at0)/sizeof(int)),  &(diag93[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at1, sizeof(cv93at1)/sizeof(int)),  &(diag93[1]) );
  } 
  else if( diag->id() == -94 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at0, sizeof(cv94at0)/sizeof(int)),  &(diag94[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at1, sizeof(cv94at1)/sizeof(int)),  &(diag94[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at2, sizeof(cv94at2)/sizeof(int)),  &(diag94[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at3, sizeof(cv94at3)/sizeof(int)),  &(diag94[3]) );
  } 
  else if( diag->id() == -95 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at0, sizeof(cv95at0)/sizeof(int)),  &(diag95[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at1, sizeof(cv95at1)/sizeof(int)),  &(diag95[1]) );
  } 
  else if( diag->id() == -96 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at0, sizeof(cv96at0)/sizeof(int)),  &(diag96[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at1, sizeof(cv96at1)/sizeof(int)),  &(diag96[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at2, sizeof(cv96at2)/sizeof(int)),  &(diag96[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at3, sizeof(cv96at3)/sizeof(int)),  &(diag96[3]) );
  } 
  else if( diag->id() == -97 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at0, sizeof(cv97at0)/sizeof(int)),  &(diag97[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at1, sizeof(cv97at1)/sizeof(int)),  &(diag97[1]) );
  } 
  else if( diag->id() == -98 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at0, sizeof(cv98at0)/sizeof(int)),  &(diag98[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at1, sizeof(cv98at1)/sizeof(int)),  &(diag98[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at2, sizeof(cv98at2)/sizeof(int)),  &(diag98[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at3, sizeof(cv98at3)/sizeof(int)),  &(diag98[3]) );
  } 
  else if( diag->id() == -99 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at0, sizeof(cv99at0)/sizeof(int)),  &(diag99[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at1, sizeof(cv99at1)/sizeof(int)),  &(diag99[1]) );
  } 
  else if( diag->id() == -100 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at0, sizeof(cv100at0)/sizeof(int)),  &(diag100[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at1, sizeof(cv100at1)/sizeof(int)),  &(diag100[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at2, sizeof(cv100at2)/sizeof(int)),  &(diag100[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at3, sizeof(cv100at3)/sizeof(int)),  &(diag100[3]) );
  } 
  else if( diag->id() == -101 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at0, sizeof(cv101at0)/sizeof(int)),  &(diag101[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at1, sizeof(cv101at1)/sizeof(int)),  &(diag101[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at2, sizeof(cv101at2)/sizeof(int)),  &(diag101[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at3, sizeof(cv101at3)/sizeof(int)),  &(diag101[3]) );
  } 
  else if( diag->id() == -102 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at0, sizeof(cv102at0)/sizeof(int)),  &(diag102[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at1, sizeof(cv102at1)/sizeof(int)),  &(diag102[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at2, sizeof(cv102at2)/sizeof(int)),  &(diag102[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at3, sizeof(cv102at3)/sizeof(int)),  &(diag102[3]) );
  } 
  else if( diag->id() == -103 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at0, sizeof(cv103at0)/sizeof(int)),  &(diag103[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at1, sizeof(cv103at1)/sizeof(int)),  &(diag103[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at2, sizeof(cv103at2)/sizeof(int)),  &(diag103[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at3, sizeof(cv103at3)/sizeof(int)),  &(diag103[3]) );
  } 
  else if( diag->id() == -104 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at0, sizeof(cv104at0)/sizeof(int)),  &(diag104[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at1, sizeof(cv104at1)/sizeof(int)),  &(diag104[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at2, sizeof(cv104at2)/sizeof(int)),  &(diag104[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at3, sizeof(cv104at3)/sizeof(int)),  &(diag104[3]) );
  } 
  else if( diag->id() == -105 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at0, sizeof(cv105at0)/sizeof(int)),  &(diag105[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at1, sizeof(cv105at1)/sizeof(int)),  &(diag105[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at2, sizeof(cv105at2)/sizeof(int)),  &(diag105[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at3, sizeof(cv105at3)/sizeof(int)),  &(diag105[3]) );
  } 
  return sel;
}
