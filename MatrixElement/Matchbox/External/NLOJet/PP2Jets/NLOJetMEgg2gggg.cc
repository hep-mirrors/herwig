// -*- C++ -*-
//
// NLOJetMEgg2gggg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEgg2gggg class.
//

#include "NLOJetMEgg2gggg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEgg2gggg::NLOJetMEgg2gggg() {}

NLOJetMEgg2gggg::~NLOJetMEgg2gggg() {}

IBPtr NLOJetMEgg2gggg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEgg2gggg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEgg2gggg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEgg2gggg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEgg2gggg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEgg2gggg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEgg2gggg("Herwig::NLOJetMEgg2gggg", "HwMatchboxNLOJet.so");

void NLOJetMEgg2gggg::Init() {

  static ClassDocumentation<NLOJetMEgg2gggg> documentation
    ("NLOJetMEgg2gggg");

}


void NLOJetMEgg2gggg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 3, g, 4, g, 4, g, 5, g, 5, g, -1)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 3, g, 4, g, 5, g, 4, g, 5, g, -2)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 3, g, 4, g, 5, g, 5, g, 4, g, -3)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 1, g, 3, g, 5, g, 5, g, -4)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 1, g, 5, g, 3, g, 5, g, -5)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 1, g, 5, g, 5, g, 3, g, -6)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 3, g, 1, g, 5, g, 5, g, -7)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, g, 1, g, 3, g, 5, g, -8)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, g, 1, g, 5, g, 3, g, -9)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 3, g, 5, g, 1, g, 5, g, -10)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, g, 3, g, 1, g, 5, g, -11)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, g, 5, g, 1, g, 3, g, -12)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 3, g, 5, g, 5, g, 1, g, -13)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, g, 3, g, 5, g, 1, g, -14)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, g, 5, g, 3, g, 1, g, -15)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, g, 5, g, 3, g, 4, g, -16)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, g, 5, g, 4, g, 3, g, -17)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, g, 3, g, 5, g, 4, g, -18)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, g, 4, g, 5, g, 3, g, -19)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, g, 3, g, 4, g, 5, g, -20)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, g, 4, g, 3, g, 5, g, -21)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 3, g, 5, g, 5, g, 4, g, -22)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 4, g, 5, g, 5, g, 3, g, -23)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 3, g, 5, g, 4, g, 5, g, -24)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 4, g, 5, g, 3, g, 5, g, -25)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 3, g, 4, g, 5, g, 5, g, -26)));
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 4, g, 3, g, 5, g, 5, g, -27)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 1, g, 4, g, 2, g, 3, g, -28)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 1, g, 4, g, 3, g, 2, g, -29)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 1, g, 2, g, 4, g, 3, g, -30)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 1, g, 3, g, 4, g, 2, g, -31)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 1, g, 2, g, 3, g, 4, g, -32)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 1, g, 3, g, 2, g, 4, g, -33)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 1, g, 5, g, 5, g, 4, g, -34)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 1, g, 5, g, 5, g, 2, g, -35)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 1, g, 5, g, 4, g, 5, g, -36)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 1, g, 5, g, 2, g, 5, g, -37)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 1, g, 4, g, 5, g, 5, g, -38)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 1, g, 2, g, 5, g, 5, g, -39)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 4, g, 1, g, 2, g, 3, g, -40)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 4, g, 1, g, 3, g, 2, g, -41)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 2, g, 1, g, 4, g, 3, g, -42)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 3, g, 1, g, 4, g, 2, g, -43)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 2, g, 1, g, 3, g, 4, g, -44)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 3, g, 1, g, 2, g, 4, g, -45)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, g, 1, g, 5, g, 4, g, -46)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, g, 1, g, 5, g, 2, g, -47)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, g, 1, g, 4, g, 5, g, -48)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, g, 1, g, 2, g, 5, g, -49)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 4, g, 1, g, 5, g, 5, g, -50)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 2, g, 1, g, 5, g, 5, g, -51)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 4, g, 2, g, 1, g, 3, g, -52)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 4, g, 3, g, 1, g, 2, g, -53)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 2, g, 4, g, 1, g, 3, g, -54)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 3, g, 4, g, 1, g, 2, g, -55)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 2, g, 3, g, 1, g, 4, g, -56)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 3, g, 2, g, 1, g, 4, g, -57)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, g, 5, g, 1, g, 4, g, -58)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, g, 5, g, 1, g, 2, g, -59)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, g, 4, g, 1, g, 5, g, -60)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, g, 2, g, 1, g, 5, g, -61)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 4, g, 5, g, 1, g, 5, g, -62)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 2, g, 5, g, 1, g, 5, g, -63)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 4, g, 2, g, 3, g, 1, g, -64)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 4, g, 3, g, 2, g, 1, g, -65)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 2, g, 4, g, 3, g, 1, g, -66)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 3, g, 4, g, 2, g, 1, g, -67)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 2, g, 3, g, 4, g, 1, g, -68)));
  addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, g, g, 3, g, 2, g, 4, g, 1, g, -69)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, g, 5, g, 4, g, 1, g, -70)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, g, 5, g, 2, g, 1, g, -71)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, g, 4, g, 5, g, 1, g, -72)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, g, 2, g, 5, g, 1, g, -73)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 4, g, 5, g, 5, g, 1, g, -74)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 2, g, 5, g, 5, g, 1, g, -75)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 2, g, 5, g, 5, g, 4, g, -76)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 3, g, 5, g, 5, g, 2, g, -77)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 2, g, 5, g, 4, g, 5, g, -78)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 3, g, 5, g, 2, g, 5, g, -79)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 2, g, 4, g, 5, g, 5, g, -80)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 3, g, 2, g, 5, g, 5, g, -81)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, g, 2, g, 5, g, 4, g, -82)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, g, 3, g, 5, g, 2, g, -83)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, g, 2, g, 4, g, 5, g, -84)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, g, 3, g, 2, g, 5, g, -85)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 4, g, 2, g, 5, g, 5, g, -86)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 2, g, 3, g, 5, g, 5, g, -87)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, g, 5, g, 2, g, 4, g, -88)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, g, 5, g, 3, g, 2, g, -89)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, g, 4, g, 2, g, 5, g, -90)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, g, 2, g, 3, g, 5, g, -91)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 4, g, 5, g, 2, g, 5, g, -92)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 2, g, 5, g, 3, g, 5, g, -93)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, g, 5, g, 4, g, 2, g, -94)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, g, 5, g, 2, g, 3, g, -95)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, g, 4, g, 5, g, 2, g, -96)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, g, 2, g, 5, g, 3, g, -97)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 4, g, 5, g, 5, g, 2, g, -98)));
  addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 2, g, 5, g, 5, g, 3, g, -99)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 2, g, 4, g, 4, g, 5, g, 5, g, -100)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 1, g, 4, g, 4, g, 5, g, 5, g, -101)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 2, g, 4, g, 5, g, 4, g, 5, g, -102)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 1, g, 4, g, 5, g, 4, g, 5, g, -103)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 2, g, 4, g, 5, g, 5, g, 4, g, -104)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 1, g, 4, g, 5, g, 5, g, 4, g, -105)));
}


Selector <MEBase::DiagramIndex> NLOJetMEgg2gggg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv1at1[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv1at2[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv1at3[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv1at4[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv1at5[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv1at6[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv1at7[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv1at8[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv1at9[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv1at10[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv1at11[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv1at12[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv1at13[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv1at14[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv1at15[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv2at0[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv2at1[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv2at2[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv2at3[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv2at4[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv2at5[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv2at6[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv2at7[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv2at8[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv2at9[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv2at10[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv2at11[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv2at12[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv2at13[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv2at14[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv2at15[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv3at0[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv3at1[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv3at2[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv3at3[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv3at4[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv3at5[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv3at6[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv3at7[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv3at8[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv3at9[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv3at10[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv3at11[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv3at12[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv3at13[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv3at14[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv3at15[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv4at0[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv4at1[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv4at2[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv4at3[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv4at4[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv4at5[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv4at6[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv4at7[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv4at8[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv4at9[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv4at10[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv4at11[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv4at12[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv4at13[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv4at14[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv4at15[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv5at0[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv5at1[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv5at2[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv5at3[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv5at4[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv5at5[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv5at6[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv5at7[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv5at8[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv5at9[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv5at10[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv5at11[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv5at12[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv5at13[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv5at14[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv5at15[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv6at0[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv6at1[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv6at2[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv6at3[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv6at4[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv6at5[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv6at6[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv6at7[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv6at8[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv6at9[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv6at10[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv6at11[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv6at12[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv6at13[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv6at14[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv6at15[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv7at0[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv7at1[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv7at2[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv7at3[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv7at4[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv7at5[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv7at6[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv7at7[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv7at8[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv7at9[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv7at10[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv7at11[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv7at12[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv7at13[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv7at14[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv7at15[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv8at0[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv8at1[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv8at2[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv8at3[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv8at4[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv8at5[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv8at6[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv8at7[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv8at8[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv8at9[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv8at10[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv8at11[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv8at12[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv8at13[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv8at14[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv8at15[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv9at0[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv9at1[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv9at2[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv9at3[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv9at4[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv9at5[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv9at6[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv9at7[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv9at8[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv9at9[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv9at10[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv9at11[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv9at12[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv9at13[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv9at14[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv9at15[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv10at0[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv10at1[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv10at2[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv10at3[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv10at4[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv10at5[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv10at6[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv10at7[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv10at8[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv10at9[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv10at10[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv10at11[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv10at12[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv10at13[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv10at14[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv10at15[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv11at0[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv11at1[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv11at2[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv11at3[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv11at4[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv11at5[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv11at6[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv11at7[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv11at8[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv11at9[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv11at10[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv11at11[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv11at12[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv11at13[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv11at14[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv11at15[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv12at0[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv12at1[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv12at2[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv12at3[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv12at4[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv12at5[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv12at6[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv12at7[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv12at8[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv12at9[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv12at10[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv12at11[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv12at12[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv12at13[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv12at14[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv12at15[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv13at0[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv13at1[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv13at2[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv13at3[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv13at4[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv13at5[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv13at6[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv13at7[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv13at8[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv13at9[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv13at10[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv13at11[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv13at12[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv13at13[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv13at14[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv13at15[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv14at0[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv14at1[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv14at2[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv14at3[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv14at4[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv14at5[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv14at6[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv14at7[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv14at8[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv14at9[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv14at10[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv14at11[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv14at12[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv14at13[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv14at14[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv14at15[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv15at0[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv15at1[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv15at2[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv15at3[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv15at4[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv15at5[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv15at6[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv15at7[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv15at8[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv15at9[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv15at10[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv15at11[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv15at12[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv15at13[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv15at14[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv15at15[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv16at0[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv16at1[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv16at2[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv16at3[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv16at4[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv16at5[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv16at6[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv16at7[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv16at8[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv16at9[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv16at10[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv16at11[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv16at12[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv16at13[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv16at14[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv16at15[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv17at0[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv17at1[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv17at2[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv17at3[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv17at4[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv17at5[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv17at6[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv17at7[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv17at8[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv17at9[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv17at10[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv17at11[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv17at12[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv17at13[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv17at14[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv17at15[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv18at0[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv18at1[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv18at2[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv18at3[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv18at4[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv18at5[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv18at6[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv18at7[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv18at8[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv18at9[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv18at10[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv18at11[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv18at12[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv18at13[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv18at14[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv18at15[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv19at0[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv19at1[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv19at2[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv19at3[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv19at4[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv19at5[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv19at6[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv19at7[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv19at8[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv19at9[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv19at10[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv19at11[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv19at12[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv19at13[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv19at14[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv19at15[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv20at0[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv20at1[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv20at2[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv20at3[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv20at4[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv20at5[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv20at6[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv20at7[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv20at8[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv20at9[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv20at10[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv20at11[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv20at12[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv20at13[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv20at14[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv20at15[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv21at0[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv21at1[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv21at2[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv21at3[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv21at4[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv21at5[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv21at6[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv21at7[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv21at8[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv21at9[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv21at10[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv21at11[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv21at12[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv21at13[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv21at14[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv21at15[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv22at0[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv22at1[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv22at2[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv22at3[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv22at4[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv22at5[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv22at6[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv22at7[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv22at8[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv22at9[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv22at10[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv22at11[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv22at12[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv22at13[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv22at14[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv22at15[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv23at0[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv23at1[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv23at2[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv23at3[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv23at4[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv23at5[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv23at6[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv23at7[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv23at8[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv23at9[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv23at10[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv23at11[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv23at12[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv23at13[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv23at14[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv23at15[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv24at0[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv24at1[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv24at2[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv24at3[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv24at4[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv24at5[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv24at6[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv24at7[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv24at8[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv24at9[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv24at10[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv24at11[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv24at12[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv24at13[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv24at14[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv24at15[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv25at0[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv25at1[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv25at2[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv25at3[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv25at4[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv25at5[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv25at6[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv25at7[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv25at8[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv25at9[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv25at10[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv25at11[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv25at12[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv25at13[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv25at14[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv25at15[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv26at0[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv26at1[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv26at2[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv26at3[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv26at4[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv26at5[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv26at6[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv26at7[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv26at8[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv26at9[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv26at10[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv26at11[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv26at12[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv26at13[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv26at14[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv26at15[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv27at0[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv27at1[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv27at2[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv27at3[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv27at4[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv27at5[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv27at6[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv27at7[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv27at8[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv27at9[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv27at10[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv27at11[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv27at12[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv27at13[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv27at14[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv27at15[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv28at0[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv28at1[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv28at2[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv28at3[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv28at4[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv28at5[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv28at6[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv28at7[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv28at8[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv28at9[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv28at10[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv28at11[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv28at12[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv28at13[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv28at14[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv28at15[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv29at0[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv29at1[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv29at2[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv29at3[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv29at4[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv29at5[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv29at6[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv29at7[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv29at8[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv29at9[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv29at10[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv29at11[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv29at12[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv29at13[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv29at14[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv29at15[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv30at0[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv30at1[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv30at2[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv30at3[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv30at4[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv30at5[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv30at6[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv30at7[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv30at8[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv30at9[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv30at10[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv30at11[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv30at12[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv30at13[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv30at14[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv30at15[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv31at0[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv31at1[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv31at2[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv31at3[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv31at4[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv31at5[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv31at6[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv31at7[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv31at8[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv31at9[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv31at10[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv31at11[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv31at12[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv31at13[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv31at14[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv31at15[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv32at0[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv32at1[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv32at2[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv32at3[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv32at4[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv32at5[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv32at6[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv32at7[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv32at8[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv32at9[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv32at10[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv32at11[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv32at12[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv32at13[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv32at14[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv32at15[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv33at0[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv33at1[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv33at2[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv33at3[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv33at4[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv33at5[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv33at6[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv33at7[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv33at8[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv33at9[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv33at10[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv33at11[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv33at12[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv33at13[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv33at14[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv33at15[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv34at0[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv34at1[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv34at2[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv34at3[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv34at4[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv34at5[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv34at6[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv34at7[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv34at8[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv34at9[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv34at10[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv34at11[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv34at12[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv34at13[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv34at14[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv34at15[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv35at0[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv35at1[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv35at2[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv35at3[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv35at4[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv35at5[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv35at6[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv35at7[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv35at8[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv35at9[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv35at10[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv35at11[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv35at12[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv35at13[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv35at14[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv35at15[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv36at0[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv36at1[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv36at2[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv36at3[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv36at4[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv36at5[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv36at6[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv36at7[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv36at8[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv36at9[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv36at10[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv36at11[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv36at12[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv36at13[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv36at14[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv36at15[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv37at0[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv37at1[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv37at2[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv37at3[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv37at4[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv37at5[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv37at6[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv37at7[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv37at8[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv37at9[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv37at10[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv37at11[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv37at12[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv37at13[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv37at14[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv37at15[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv38at0[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv38at1[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv38at2[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv38at3[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv38at4[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv38at5[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv38at6[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv38at7[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv38at8[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv38at9[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv38at10[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv38at11[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv38at12[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv38at13[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv38at14[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv38at15[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv39at0[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv39at1[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv39at2[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv39at3[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv39at4[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv39at5[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv39at6[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv39at7[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv39at8[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv39at9[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv39at10[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv39at11[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv39at12[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv39at13[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv39at14[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv39at15[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv40at0[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv40at1[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv40at2[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv40at3[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv40at4[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv40at5[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv40at6[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv40at7[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv40at8[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv40at9[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv40at10[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv40at11[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv40at12[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv40at13[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv40at14[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv40at15[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv41at0[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv41at1[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv41at2[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv41at3[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv41at4[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv41at5[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv41at6[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv41at7[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv41at8[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv41at9[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv41at10[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv41at11[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv41at12[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv41at13[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv41at14[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv41at15[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv42at0[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv42at1[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv42at2[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv42at3[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv42at4[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv42at5[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv42at6[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv42at7[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv42at8[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv42at9[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv42at10[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv42at11[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv42at12[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv42at13[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv42at14[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv42at15[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv43at0[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv43at1[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv43at2[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv43at3[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv43at4[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv43at5[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv43at6[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv43at7[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv43at8[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv43at9[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv43at10[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv43at11[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv43at12[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv43at13[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv43at14[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv43at15[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv44at0[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv44at1[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv44at2[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv44at3[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv44at4[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv44at5[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv44at6[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv44at7[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv44at8[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv44at9[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv44at10[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv44at11[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv44at12[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv44at13[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv44at14[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv44at15[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv45at0[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv45at1[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv45at2[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv45at3[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv45at4[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv45at5[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv45at6[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv45at7[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv45at8[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv45at9[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv45at10[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv45at11[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv45at12[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv45at13[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv45at14[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv45at15[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv46at0[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv46at1[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv46at2[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv46at3[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv46at4[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv46at5[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv46at6[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv46at7[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv46at8[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv46at9[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv46at10[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv46at11[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv46at12[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv46at13[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv46at14[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv46at15[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv47at0[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv47at1[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv47at2[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv47at3[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv47at4[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv47at5[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv47at6[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv47at7[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv47at8[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv47at9[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv47at10[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv47at11[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv47at12[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv47at13[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv47at14[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv47at15[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv48at0[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv48at1[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv48at2[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv48at3[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv48at4[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv48at5[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv48at6[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv48at7[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv48at8[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv48at9[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv48at10[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv48at11[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv48at12[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv48at13[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv48at14[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv48at15[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv49at0[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv49at1[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv49at2[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv49at3[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv49at4[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv49at5[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv49at6[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv49at7[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv49at8[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv49at9[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv49at10[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv49at11[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv49at12[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv49at13[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv49at14[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv49at15[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv50at0[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv50at1[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv50at2[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv50at3[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv50at4[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv50at5[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv50at6[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv50at7[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv50at8[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv50at9[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv50at10[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv50at11[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv50at12[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv50at13[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv50at14[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv50at15[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv51at0[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv51at1[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv51at2[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv51at3[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv51at4[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv51at5[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv51at6[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv51at7[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv51at8[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv51at9[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv51at10[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv51at11[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv51at12[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv51at13[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv51at14[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv51at15[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv52at0[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv52at1[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv52at2[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv52at3[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv52at4[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv52at5[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv52at6[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv52at7[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv52at8[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv52at9[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv52at10[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv52at11[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv52at12[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv52at13[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv52at14[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv52at15[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv53at0[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv53at1[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv53at2[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv53at3[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv53at4[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv53at5[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv53at6[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv53at7[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv53at8[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv53at9[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv53at10[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv53at11[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv53at12[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv53at13[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv53at14[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv53at15[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv54at0[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv54at1[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv54at2[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv54at3[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv54at4[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv54at5[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv54at6[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv54at7[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv54at8[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv54at9[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv54at10[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv54at11[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv54at12[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv54at13[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv54at14[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv54at15[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv55at0[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv55at1[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv55at2[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv55at3[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv55at4[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv55at5[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv55at6[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv55at7[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv55at8[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv55at9[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv55at10[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv55at11[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv55at12[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv55at13[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv55at14[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv55at15[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv56at0[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv56at1[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv56at2[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv56at3[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv56at4[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv56at5[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv56at6[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv56at7[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv56at8[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv56at9[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv56at10[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv56at11[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv56at12[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv56at13[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv56at14[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv56at15[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv57at0[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv57at1[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv57at2[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv57at3[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv57at4[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv57at5[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv57at6[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv57at7[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv57at8[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv57at9[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv57at10[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv57at11[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv57at12[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv57at13[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv57at14[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv57at15[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv58at0[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv58at1[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv58at2[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv58at3[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv58at4[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv58at5[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv58at6[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv58at7[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv58at8[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv58at9[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv58at10[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv58at11[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv58at12[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv58at13[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv58at14[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv58at15[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv59at0[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv59at1[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv59at2[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv59at3[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv59at4[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv59at5[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv59at6[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv59at7[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv59at8[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv59at9[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv59at10[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv59at11[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv59at12[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv59at13[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv59at14[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv59at15[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv60at0[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv60at1[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv60at2[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv60at3[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv60at4[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv60at5[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv60at6[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv60at7[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv60at8[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv60at9[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv60at10[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv60at11[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv60at12[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv60at13[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv60at14[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv60at15[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv61at0[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv61at1[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv61at2[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv61at3[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv61at4[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv61at5[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv61at6[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv61at7[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv61at8[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv61at9[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv61at10[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv61at11[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv61at12[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv61at13[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv61at14[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv61at15[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv62at0[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv62at1[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv62at2[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv62at3[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv62at4[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv62at5[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv62at6[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv62at7[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv62at8[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv62at9[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv62at10[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv62at11[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv62at12[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv62at13[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv62at14[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv62at15[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv63at0[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv63at1[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv63at2[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv63at3[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv63at4[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv63at5[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv63at6[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv63at7[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv63at8[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv63at9[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv63at10[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv63at11[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv63at12[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv63at13[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv63at14[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv63at15[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv64at0[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv64at1[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv64at2[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv64at3[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv64at4[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv64at5[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv64at6[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv64at7[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv64at8[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv64at9[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv64at10[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv64at11[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv64at12[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv64at13[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv64at14[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv64at15[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv65at0[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv65at1[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv65at2[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv65at3[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv65at4[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv65at5[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv65at6[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv65at7[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv65at8[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv65at9[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv65at10[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv65at11[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv65at12[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv65at13[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv65at14[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv65at15[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv66at0[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv66at1[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv66at2[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv66at3[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv66at4[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv66at5[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv66at6[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv66at7[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv66at8[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv66at9[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv66at10[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv66at11[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv66at12[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv66at13[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv66at14[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv66at15[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv67at0[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv67at1[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv67at2[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv67at3[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv67at4[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv67at5[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv67at6[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv67at7[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv67at8[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv67at9[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv67at10[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv67at11[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv67at12[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv67at13[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv67at14[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv67at15[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv68at0[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv68at1[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv68at2[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv68at3[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv68at4[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv68at5[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv68at6[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv68at7[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv68at8[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv68at9[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv68at10[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv68at11[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv68at12[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv68at13[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv68at14[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv68at15[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv69at0[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv69at1[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv69at2[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv69at3[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv69at4[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv69at5[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv69at6[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv69at7[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv69at8[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv69at9[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv69at10[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv69at11[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv69at12[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv69at13[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv69at14[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv69at15[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv70at0[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv70at1[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv70at2[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv70at3[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv70at4[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv70at5[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv70at6[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv70at7[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv70at8[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv70at9[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv70at10[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv70at11[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv70at12[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv70at13[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv70at14[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv70at15[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv71at0[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv71at1[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv71at2[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv71at3[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv71at4[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv71at5[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv71at6[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv71at7[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv71at8[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv71at9[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv71at10[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv71at11[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv71at12[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv71at13[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv71at14[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv71at15[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv72at0[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv72at1[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv72at2[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv72at3[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv72at4[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv72at5[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv72at6[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv72at7[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv72at8[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv72at9[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv72at10[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv72at11[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv72at12[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv72at13[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv72at14[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv72at15[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv73at0[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv73at1[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv73at2[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv73at3[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv73at4[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv73at5[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv73at6[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv73at7[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv73at8[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv73at9[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv73at10[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv73at11[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv73at12[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv73at13[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv73at14[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv73at15[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv74at0[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv74at1[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv74at2[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv74at3[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv74at4[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv74at5[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv74at6[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv74at7[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv74at8[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv74at9[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv74at10[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv74at11[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv74at12[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv74at13[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv74at14[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv74at15[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv75at0[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv75at1[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv75at2[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv75at3[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv75at4[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv75at5[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv75at6[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv75at7[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv75at8[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv75at9[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv75at10[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv75at11[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv75at12[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv75at13[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv75at14[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv75at15[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv76at0[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv76at1[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv76at2[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv76at3[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv76at4[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv76at5[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv76at6[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv76at7[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv76at8[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv76at9[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv76at10[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv76at11[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv76at12[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv76at13[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv76at14[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv76at15[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv77at0[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv77at1[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv77at2[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv77at3[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv77at4[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv77at5[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv77at6[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv77at7[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv77at8[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv77at9[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv77at10[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv77at11[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv77at12[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv77at13[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv77at14[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv77at15[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv78at0[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv78at1[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv78at2[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv78at3[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv78at4[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv78at5[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv78at6[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv78at7[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv78at8[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv78at9[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv78at10[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv78at11[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv78at12[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv78at13[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv78at14[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv78at15[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv79at0[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv79at1[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv79at2[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv79at3[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv79at4[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv79at5[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv79at6[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv79at7[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv79at8[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv79at9[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv79at10[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv79at11[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv79at12[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv79at13[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv79at14[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv79at15[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv80at0[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv80at1[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv80at2[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv80at3[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv80at4[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv80at5[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv80at6[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv80at7[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv80at8[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv80at9[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv80at10[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv80at11[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv80at12[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv80at13[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv80at14[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv80at15[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv81at0[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv81at1[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv81at2[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv81at3[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv81at4[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv81at5[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv81at6[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv81at7[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv81at8[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv81at9[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv81at10[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv81at11[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv81at12[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv81at13[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv81at14[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv81at15[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv82at0[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv82at1[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv82at2[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv82at3[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv82at4[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv82at5[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv82at6[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv82at7[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv82at8[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv82at9[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv82at10[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv82at11[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv82at12[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv82at13[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv82at14[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv82at15[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv83at0[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv83at1[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv83at2[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv83at3[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv83at4[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv83at5[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv83at6[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv83at7[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv83at8[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv83at9[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv83at10[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv83at11[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv83at12[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv83at13[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv83at14[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv83at15[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv84at0[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv84at1[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv84at2[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv84at3[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv84at4[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv84at5[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv84at6[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv84at7[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv84at8[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv84at9[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv84at10[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv84at11[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv84at12[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv84at13[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv84at14[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv84at15[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv85at0[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv85at1[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv85at2[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv85at3[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv85at4[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv85at5[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv85at6[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv85at7[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv85at8[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv85at9[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv85at10[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv85at11[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv85at12[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv85at13[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv85at14[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv85at15[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv86at0[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv86at1[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv86at2[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv86at3[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv86at4[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv86at5[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv86at6[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv86at7[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv86at8[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv86at9[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv86at10[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv86at11[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv86at12[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv86at13[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv86at14[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv86at15[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv87at0[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv87at1[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv87at2[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv87at3[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv87at4[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv87at5[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv87at6[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv87at7[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv87at8[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv87at9[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv87at10[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv87at11[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv87at12[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv87at13[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv87at14[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv87at15[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv88at0[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv88at1[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv88at2[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv88at3[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv88at4[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv88at5[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv88at6[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv88at7[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv88at8[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv88at9[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv88at10[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv88at11[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv88at12[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv88at13[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv88at14[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv88at15[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv89at0[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv89at1[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv89at2[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv89at3[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv89at4[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv89at5[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv89at6[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv89at7[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv89at8[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv89at9[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv89at10[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv89at11[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv89at12[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv89at13[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv89at14[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv89at15[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv90at0[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv90at1[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv90at2[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv90at3[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv90at4[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv90at5[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv90at6[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv90at7[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv90at8[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv90at9[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv90at10[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv90at11[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv90at12[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv90at13[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv90at14[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv90at15[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv91at0[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv91at1[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv91at2[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv91at3[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv91at4[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv91at5[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv91at6[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv91at7[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv91at8[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv91at9[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv91at10[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv91at11[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv91at12[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv91at13[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv91at14[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv91at15[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv92at0[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv92at1[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv92at2[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv92at3[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv92at4[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv92at5[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv92at6[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv92at7[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv92at8[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv92at9[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv92at10[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv92at11[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv92at12[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv92at13[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv92at14[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv92at15[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv93at0[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv93at1[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv93at2[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv93at3[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv93at4[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv93at5[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv93at6[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv93at7[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv93at8[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv93at9[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv93at10[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv93at11[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv93at12[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv93at13[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv93at14[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv93at15[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv94at0[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv94at1[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv94at2[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv94at3[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv94at4[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv94at5[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv94at6[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv94at7[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv94at8[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv94at9[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv94at10[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv94at11[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv94at12[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv94at13[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv94at14[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv94at15[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv95at0[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv95at1[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv95at2[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv95at3[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv95at4[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv95at5[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv95at6[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv95at7[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv95at8[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv95at9[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv95at10[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv95at11[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv95at12[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv95at13[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv95at14[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv95at15[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv96at0[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv96at1[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv96at2[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv96at3[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv96at4[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv96at5[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv96at6[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv96at7[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv96at8[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv96at9[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv96at10[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv96at11[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv96at12[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv96at13[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv96at14[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv96at15[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv97at0[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv97at1[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv97at2[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv97at3[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv97at4[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv97at5[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv97at6[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv97at7[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv97at8[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv97at9[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv97at10[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv97at11[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv97at12[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv97at13[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv97at14[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv97at15[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv98at0[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv98at1[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv98at2[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv98at3[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv98at4[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv98at5[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv98at6[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv98at7[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv98at8[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv98at9[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv98at10[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv98at11[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv98at12[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv98at13[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv98at14[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv98at15[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv99at0[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv99at1[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv99at2[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv99at3[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv99at4[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv99at5[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv99at6[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv99at7[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv99at8[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv99at9[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv99at10[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv99at11[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv99at12[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv99at13[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv99at14[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv99at15[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv100at0[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv100at1[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv100at2[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv100at3[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv100at4[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv100at5[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv100at6[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv100at7[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv100at8[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv100at9[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv100at10[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv100at11[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv100at12[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv100at13[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv100at14[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv100at15[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv101at0[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv101at1[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv101at2[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv101at3[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv101at4[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv101at5[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv101at6[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv101at7[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv101at8[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv101at9[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv101at10[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv101at11[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv101at12[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv101at13[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv101at14[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv101at15[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv102at0[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv102at1[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv102at2[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv102at3[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv102at4[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv102at5[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv102at6[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv102at7[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv102at8[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv102at9[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv102at10[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv102at11[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv102at12[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv102at13[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv102at14[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv102at15[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv103at0[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv103at1[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv103at2[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv103at3[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv103at4[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv103at5[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv103at6[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv103at7[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv103at8[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv103at9[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv103at10[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv103at11[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv103at12[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv103at13[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv103at14[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv103at15[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv104at0[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv104at1[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv104at2[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv104at3[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv104at4[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv104at5[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv104at6[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv104at7[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv104at8[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv104at9[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv104at10[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv104at11[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv104at12[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv104at13[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv104at14[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv104at15[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv105at0[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv105at1[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv105at2[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv105at3[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv105at4[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv105at5[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv105at6[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv105at7[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv105at8[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv105at9[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv105at10[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv105at11[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv105at12[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv105at13[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv105at14[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv105at15[]  = { -1, 2, 3, 0, 1, 4, -998};

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
        + nloJetAmplitude()->colourOrdered2(cv1at7, sizeof(cv1at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at8, sizeof(cv1at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at9, sizeof(cv1at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at10, sizeof(cv1at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at11, sizeof(cv1at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at12, sizeof(cv1at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at13, sizeof(cv1at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at14, sizeof(cv1at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at15, sizeof(cv1at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at2, sizeof(cv2at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at3, sizeof(cv2at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at4, sizeof(cv2at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at5, sizeof(cv2at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at6, sizeof(cv2at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at7, sizeof(cv2at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at8, sizeof(cv2at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at9, sizeof(cv2at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at10, sizeof(cv2at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at11, sizeof(cv2at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at12, sizeof(cv2at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at13, sizeof(cv2at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at14, sizeof(cv2at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at15, sizeof(cv2at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at2, sizeof(cv3at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at3, sizeof(cv3at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at4, sizeof(cv3at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at5, sizeof(cv3at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at6, sizeof(cv3at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at7, sizeof(cv3at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at8, sizeof(cv3at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at9, sizeof(cv3at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at10, sizeof(cv3at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at11, sizeof(cv3at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at12, sizeof(cv3at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at13, sizeof(cv3at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at14, sizeof(cv3at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at15, sizeof(cv3at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -4 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at2, sizeof(cv4at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at3, sizeof(cv4at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at4, sizeof(cv4at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at5, sizeof(cv4at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at6, sizeof(cv4at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at7, sizeof(cv4at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at8, sizeof(cv4at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at9, sizeof(cv4at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at10, sizeof(cv4at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at11, sizeof(cv4at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at12, sizeof(cv4at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at13, sizeof(cv4at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at14, sizeof(cv4at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at15, sizeof(cv4at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -5 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at2, sizeof(cv5at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at3, sizeof(cv5at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at4, sizeof(cv5at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at5, sizeof(cv5at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at6, sizeof(cv5at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at7, sizeof(cv5at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at8, sizeof(cv5at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at9, sizeof(cv5at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at10, sizeof(cv5at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at11, sizeof(cv5at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at12, sizeof(cv5at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at13, sizeof(cv5at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at14, sizeof(cv5at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at15, sizeof(cv5at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -6 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at2, sizeof(cv6at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at3, sizeof(cv6at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at4, sizeof(cv6at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at5, sizeof(cv6at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at6, sizeof(cv6at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at7, sizeof(cv6at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at8, sizeof(cv6at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at9, sizeof(cv6at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at10, sizeof(cv6at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at11, sizeof(cv6at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at12, sizeof(cv6at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at13, sizeof(cv6at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at14, sizeof(cv6at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at15, sizeof(cv6at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -7 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at1, sizeof(cv7at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at2, sizeof(cv7at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at3, sizeof(cv7at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at4, sizeof(cv7at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at5, sizeof(cv7at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at6, sizeof(cv7at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at7, sizeof(cv7at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at8, sizeof(cv7at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at9, sizeof(cv7at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at10, sizeof(cv7at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at11, sizeof(cv7at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at12, sizeof(cv7at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at13, sizeof(cv7at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at14, sizeof(cv7at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv7at15, sizeof(cv7at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -8 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at2, sizeof(cv8at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at3, sizeof(cv8at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at4, sizeof(cv8at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at5, sizeof(cv8at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at6, sizeof(cv8at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at7, sizeof(cv8at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at8, sizeof(cv8at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at9, sizeof(cv8at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at10, sizeof(cv8at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at11, sizeof(cv8at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at12, sizeof(cv8at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at13, sizeof(cv8at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at14, sizeof(cv8at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at15, sizeof(cv8at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -9 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at2, sizeof(cv9at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at3, sizeof(cv9at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at4, sizeof(cv9at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at5, sizeof(cv9at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at6, sizeof(cv9at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at7, sizeof(cv9at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at8, sizeof(cv9at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at9, sizeof(cv9at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at10, sizeof(cv9at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at11, sizeof(cv9at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at12, sizeof(cv9at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at13, sizeof(cv9at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at14, sizeof(cv9at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at15, sizeof(cv9at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -10 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at2, sizeof(cv10at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at3, sizeof(cv10at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at4, sizeof(cv10at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at5, sizeof(cv10at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at6, sizeof(cv10at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at7, sizeof(cv10at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at8, sizeof(cv10at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at9, sizeof(cv10at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at10, sizeof(cv10at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at11, sizeof(cv10at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at12, sizeof(cv10at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at13, sizeof(cv10at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at14, sizeof(cv10at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at15, sizeof(cv10at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -11 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at2, sizeof(cv11at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at3, sizeof(cv11at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at4, sizeof(cv11at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at5, sizeof(cv11at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at6, sizeof(cv11at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at7, sizeof(cv11at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at8, sizeof(cv11at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at9, sizeof(cv11at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at10, sizeof(cv11at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at11, sizeof(cv11at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at12, sizeof(cv11at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at13, sizeof(cv11at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at14, sizeof(cv11at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at15, sizeof(cv11at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -12 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at1, sizeof(cv12at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at2, sizeof(cv12at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at3, sizeof(cv12at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at4, sizeof(cv12at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at5, sizeof(cv12at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at6, sizeof(cv12at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at7, sizeof(cv12at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at8, sizeof(cv12at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at9, sizeof(cv12at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at10, sizeof(cv12at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at11, sizeof(cv12at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at12, sizeof(cv12at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at13, sizeof(cv12at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at14, sizeof(cv12at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at15, sizeof(cv12at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -13 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at1, sizeof(cv13at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at2, sizeof(cv13at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at3, sizeof(cv13at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at4, sizeof(cv13at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at5, sizeof(cv13at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at6, sizeof(cv13at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at7, sizeof(cv13at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at8, sizeof(cv13at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at9, sizeof(cv13at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at10, sizeof(cv13at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at11, sizeof(cv13at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at12, sizeof(cv13at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at13, sizeof(cv13at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at14, sizeof(cv13at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv13at15, sizeof(cv13at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -14 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at1, sizeof(cv14at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at2, sizeof(cv14at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at3, sizeof(cv14at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at4, sizeof(cv14at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at5, sizeof(cv14at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at6, sizeof(cv14at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at7, sizeof(cv14at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at8, sizeof(cv14at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at9, sizeof(cv14at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at10, sizeof(cv14at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at11, sizeof(cv14at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at12, sizeof(cv14at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at13, sizeof(cv14at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at14, sizeof(cv14at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv14at15, sizeof(cv14at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -15 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at2, sizeof(cv15at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at3, sizeof(cv15at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at4, sizeof(cv15at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at5, sizeof(cv15at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at6, sizeof(cv15at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at7, sizeof(cv15at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at8, sizeof(cv15at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at9, sizeof(cv15at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at10, sizeof(cv15at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at11, sizeof(cv15at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at12, sizeof(cv15at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at13, sizeof(cv15at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at14, sizeof(cv15at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at15, sizeof(cv15at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -16 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv16at0, sizeof(cv16at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at1, sizeof(cv16at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at2, sizeof(cv16at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at3, sizeof(cv16at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at4, sizeof(cv16at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at5, sizeof(cv16at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at6, sizeof(cv16at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at7, sizeof(cv16at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at8, sizeof(cv16at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at9, sizeof(cv16at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at10, sizeof(cv16at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at11, sizeof(cv16at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at12, sizeof(cv16at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at13, sizeof(cv16at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at14, sizeof(cv16at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv16at15, sizeof(cv16at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -17 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv17at0, sizeof(cv17at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at1, sizeof(cv17at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at2, sizeof(cv17at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at3, sizeof(cv17at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at4, sizeof(cv17at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at5, sizeof(cv17at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at6, sizeof(cv17at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at7, sizeof(cv17at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at8, sizeof(cv17at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at9, sizeof(cv17at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at10, sizeof(cv17at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at11, sizeof(cv17at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at12, sizeof(cv17at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at13, sizeof(cv17at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at14, sizeof(cv17at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at15, sizeof(cv17at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -18 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv18at0, sizeof(cv18at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at1, sizeof(cv18at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at2, sizeof(cv18at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at3, sizeof(cv18at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at4, sizeof(cv18at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at5, sizeof(cv18at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at6, sizeof(cv18at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at7, sizeof(cv18at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at8, sizeof(cv18at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at9, sizeof(cv18at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at10, sizeof(cv18at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at11, sizeof(cv18at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at12, sizeof(cv18at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at13, sizeof(cv18at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at14, sizeof(cv18at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at15, sizeof(cv18at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -19 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at1, sizeof(cv19at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at2, sizeof(cv19at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at3, sizeof(cv19at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at4, sizeof(cv19at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at5, sizeof(cv19at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at6, sizeof(cv19at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at7, sizeof(cv19at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at8, sizeof(cv19at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at9, sizeof(cv19at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at10, sizeof(cv19at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at11, sizeof(cv19at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at12, sizeof(cv19at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at13, sizeof(cv19at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at14, sizeof(cv19at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at15, sizeof(cv19at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -20 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at1, sizeof(cv20at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at2, sizeof(cv20at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at3, sizeof(cv20at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at4, sizeof(cv20at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at5, sizeof(cv20at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at6, sizeof(cv20at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at7, sizeof(cv20at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at8, sizeof(cv20at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at9, sizeof(cv20at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at10, sizeof(cv20at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at11, sizeof(cv20at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at12, sizeof(cv20at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at13, sizeof(cv20at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at14, sizeof(cv20at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at15, sizeof(cv20at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -21 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at1, sizeof(cv21at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at2, sizeof(cv21at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at3, sizeof(cv21at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at4, sizeof(cv21at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at5, sizeof(cv21at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at6, sizeof(cv21at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at7, sizeof(cv21at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at8, sizeof(cv21at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at9, sizeof(cv21at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at10, sizeof(cv21at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at11, sizeof(cv21at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at12, sizeof(cv21at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at13, sizeof(cv21at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at14, sizeof(cv21at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at15, sizeof(cv21at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -22 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at1, sizeof(cv22at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at2, sizeof(cv22at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at3, sizeof(cv22at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at4, sizeof(cv22at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at5, sizeof(cv22at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at6, sizeof(cv22at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at7, sizeof(cv22at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at8, sizeof(cv22at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at9, sizeof(cv22at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at10, sizeof(cv22at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at11, sizeof(cv22at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at12, sizeof(cv22at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at13, sizeof(cv22at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at14, sizeof(cv22at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at15, sizeof(cv22at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -23 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at1, sizeof(cv23at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at2, sizeof(cv23at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at3, sizeof(cv23at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at4, sizeof(cv23at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at5, sizeof(cv23at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at6, sizeof(cv23at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at7, sizeof(cv23at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at8, sizeof(cv23at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at9, sizeof(cv23at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at10, sizeof(cv23at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at11, sizeof(cv23at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at12, sizeof(cv23at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at13, sizeof(cv23at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at14, sizeof(cv23at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at15, sizeof(cv23at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -24 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at1, sizeof(cv24at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at2, sizeof(cv24at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at3, sizeof(cv24at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at4, sizeof(cv24at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at5, sizeof(cv24at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at6, sizeof(cv24at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at7, sizeof(cv24at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at8, sizeof(cv24at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at9, sizeof(cv24at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at10, sizeof(cv24at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at11, sizeof(cv24at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at12, sizeof(cv24at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at13, sizeof(cv24at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at14, sizeof(cv24at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at15, sizeof(cv24at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -25 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at1, sizeof(cv25at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at2, sizeof(cv25at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at3, sizeof(cv25at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at4, sizeof(cv25at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at5, sizeof(cv25at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at6, sizeof(cv25at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at7, sizeof(cv25at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at8, sizeof(cv25at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at9, sizeof(cv25at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at10, sizeof(cv25at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at11, sizeof(cv25at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at12, sizeof(cv25at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at13, sizeof(cv25at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at14, sizeof(cv25at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at15, sizeof(cv25at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -26 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at1, sizeof(cv26at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at2, sizeof(cv26at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at3, sizeof(cv26at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at4, sizeof(cv26at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at5, sizeof(cv26at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at6, sizeof(cv26at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at7, sizeof(cv26at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at8, sizeof(cv26at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at9, sizeof(cv26at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at10, sizeof(cv26at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at11, sizeof(cv26at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at12, sizeof(cv26at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at13, sizeof(cv26at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at14, sizeof(cv26at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at15, sizeof(cv26at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -27 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at1, sizeof(cv27at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at2, sizeof(cv27at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at3, sizeof(cv27at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at4, sizeof(cv27at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at5, sizeof(cv27at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at6, sizeof(cv27at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at7, sizeof(cv27at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at8, sizeof(cv27at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at9, sizeof(cv27at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at10, sizeof(cv27at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at11, sizeof(cv27at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at12, sizeof(cv27at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at13, sizeof(cv27at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at14, sizeof(cv27at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at15, sizeof(cv27at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -28 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at1, sizeof(cv28at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at2, sizeof(cv28at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at3, sizeof(cv28at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at4, sizeof(cv28at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at5, sizeof(cv28at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at6, sizeof(cv28at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at7, sizeof(cv28at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at8, sizeof(cv28at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at9, sizeof(cv28at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at10, sizeof(cv28at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at11, sizeof(cv28at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at12, sizeof(cv28at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at13, sizeof(cv28at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at14, sizeof(cv28at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at15, sizeof(cv28at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -29 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at1, sizeof(cv29at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at2, sizeof(cv29at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at3, sizeof(cv29at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at4, sizeof(cv29at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at5, sizeof(cv29at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at6, sizeof(cv29at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at7, sizeof(cv29at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at8, sizeof(cv29at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at9, sizeof(cv29at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at10, sizeof(cv29at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at11, sizeof(cv29at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at12, sizeof(cv29at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at13, sizeof(cv29at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at14, sizeof(cv29at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at15, sizeof(cv29at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -30 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at1, sizeof(cv30at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at2, sizeof(cv30at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at3, sizeof(cv30at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at4, sizeof(cv30at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at5, sizeof(cv30at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at6, sizeof(cv30at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at7, sizeof(cv30at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at8, sizeof(cv30at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at9, sizeof(cv30at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at10, sizeof(cv30at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at11, sizeof(cv30at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at12, sizeof(cv30at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at13, sizeof(cv30at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at14, sizeof(cv30at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at15, sizeof(cv30at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -31 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at1, sizeof(cv31at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at2, sizeof(cv31at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at3, sizeof(cv31at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at4, sizeof(cv31at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at5, sizeof(cv31at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at6, sizeof(cv31at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at7, sizeof(cv31at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at8, sizeof(cv31at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at9, sizeof(cv31at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at10, sizeof(cv31at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at11, sizeof(cv31at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at12, sizeof(cv31at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at13, sizeof(cv31at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at14, sizeof(cv31at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at15, sizeof(cv31at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -32 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at1, sizeof(cv32at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at2, sizeof(cv32at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at3, sizeof(cv32at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at4, sizeof(cv32at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at5, sizeof(cv32at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at6, sizeof(cv32at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at7, sizeof(cv32at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at8, sizeof(cv32at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at9, sizeof(cv32at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at10, sizeof(cv32at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at11, sizeof(cv32at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at12, sizeof(cv32at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at13, sizeof(cv32at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at14, sizeof(cv32at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at15, sizeof(cv32at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -33 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at1, sizeof(cv33at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at2, sizeof(cv33at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at3, sizeof(cv33at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at4, sizeof(cv33at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at5, sizeof(cv33at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at6, sizeof(cv33at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at7, sizeof(cv33at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at8, sizeof(cv33at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at9, sizeof(cv33at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at10, sizeof(cv33at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at11, sizeof(cv33at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at12, sizeof(cv33at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at13, sizeof(cv33at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at14, sizeof(cv33at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at15, sizeof(cv33at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -34 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv34at0, sizeof(cv34at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at1, sizeof(cv34at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at2, sizeof(cv34at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at3, sizeof(cv34at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at4, sizeof(cv34at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at5, sizeof(cv34at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at6, sizeof(cv34at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at7, sizeof(cv34at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at8, sizeof(cv34at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at9, sizeof(cv34at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at10, sizeof(cv34at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at11, sizeof(cv34at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at12, sizeof(cv34at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at13, sizeof(cv34at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at14, sizeof(cv34at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv34at15, sizeof(cv34at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -35 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv35at0, sizeof(cv35at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at1, sizeof(cv35at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at2, sizeof(cv35at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at3, sizeof(cv35at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at4, sizeof(cv35at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at5, sizeof(cv35at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at6, sizeof(cv35at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at7, sizeof(cv35at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at8, sizeof(cv35at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at9, sizeof(cv35at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at10, sizeof(cv35at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at11, sizeof(cv35at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at12, sizeof(cv35at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at13, sizeof(cv35at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at14, sizeof(cv35at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at15, sizeof(cv35at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -36 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at1, sizeof(cv36at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at2, sizeof(cv36at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at3, sizeof(cv36at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at4, sizeof(cv36at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at5, sizeof(cv36at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at6, sizeof(cv36at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at7, sizeof(cv36at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at8, sizeof(cv36at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at9, sizeof(cv36at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at10, sizeof(cv36at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at11, sizeof(cv36at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at12, sizeof(cv36at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at13, sizeof(cv36at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at14, sizeof(cv36at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv36at15, sizeof(cv36at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -37 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at1, sizeof(cv37at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at2, sizeof(cv37at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at3, sizeof(cv37at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at4, sizeof(cv37at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at5, sizeof(cv37at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at6, sizeof(cv37at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at7, sizeof(cv37at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at8, sizeof(cv37at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at9, sizeof(cv37at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at10, sizeof(cv37at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at11, sizeof(cv37at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at12, sizeof(cv37at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at13, sizeof(cv37at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at14, sizeof(cv37at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv37at15, sizeof(cv37at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -38 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at1, sizeof(cv38at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at2, sizeof(cv38at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at3, sizeof(cv38at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at4, sizeof(cv38at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at5, sizeof(cv38at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at6, sizeof(cv38at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at7, sizeof(cv38at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at8, sizeof(cv38at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at9, sizeof(cv38at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at10, sizeof(cv38at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at11, sizeof(cv38at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at12, sizeof(cv38at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at13, sizeof(cv38at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at14, sizeof(cv38at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at15, sizeof(cv38at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -39 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at1, sizeof(cv39at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at2, sizeof(cv39at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at3, sizeof(cv39at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at4, sizeof(cv39at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at5, sizeof(cv39at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at6, sizeof(cv39at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at7, sizeof(cv39at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at8, sizeof(cv39at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at9, sizeof(cv39at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at10, sizeof(cv39at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at11, sizeof(cv39at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at12, sizeof(cv39at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at13, sizeof(cv39at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at14, sizeof(cv39at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at15, sizeof(cv39at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -40 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at1, sizeof(cv40at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at2, sizeof(cv40at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at3, sizeof(cv40at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at4, sizeof(cv40at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at5, sizeof(cv40at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at6, sizeof(cv40at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at7, sizeof(cv40at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at8, sizeof(cv40at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at9, sizeof(cv40at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at10, sizeof(cv40at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at11, sizeof(cv40at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at12, sizeof(cv40at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at13, sizeof(cv40at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at14, sizeof(cv40at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv40at15, sizeof(cv40at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -41 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at1, sizeof(cv41at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at2, sizeof(cv41at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at3, sizeof(cv41at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at4, sizeof(cv41at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at5, sizeof(cv41at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at6, sizeof(cv41at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at7, sizeof(cv41at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at8, sizeof(cv41at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at9, sizeof(cv41at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at10, sizeof(cv41at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at11, sizeof(cv41at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at12, sizeof(cv41at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at13, sizeof(cv41at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at14, sizeof(cv41at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at15, sizeof(cv41at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -42 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at1, sizeof(cv42at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at2, sizeof(cv42at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at3, sizeof(cv42at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at4, sizeof(cv42at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at5, sizeof(cv42at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at6, sizeof(cv42at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at7, sizeof(cv42at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at8, sizeof(cv42at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at9, sizeof(cv42at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at10, sizeof(cv42at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at11, sizeof(cv42at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at12, sizeof(cv42at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at13, sizeof(cv42at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at14, sizeof(cv42at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at15, sizeof(cv42at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -43 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv43at0, sizeof(cv43at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at1, sizeof(cv43at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at2, sizeof(cv43at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at3, sizeof(cv43at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at4, sizeof(cv43at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at5, sizeof(cv43at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at6, sizeof(cv43at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at7, sizeof(cv43at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at8, sizeof(cv43at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at9, sizeof(cv43at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at10, sizeof(cv43at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at11, sizeof(cv43at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at12, sizeof(cv43at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at13, sizeof(cv43at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at14, sizeof(cv43at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at15, sizeof(cv43at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -44 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv44at0, sizeof(cv44at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at1, sizeof(cv44at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at2, sizeof(cv44at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at3, sizeof(cv44at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at4, sizeof(cv44at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at5, sizeof(cv44at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at6, sizeof(cv44at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at7, sizeof(cv44at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at8, sizeof(cv44at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at9, sizeof(cv44at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at10, sizeof(cv44at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at11, sizeof(cv44at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at12, sizeof(cv44at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at13, sizeof(cv44at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at14, sizeof(cv44at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at15, sizeof(cv44at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -45 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv45at0, sizeof(cv45at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at1, sizeof(cv45at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at2, sizeof(cv45at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at3, sizeof(cv45at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at4, sizeof(cv45at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at5, sizeof(cv45at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at6, sizeof(cv45at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at7, sizeof(cv45at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at8, sizeof(cv45at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at9, sizeof(cv45at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at10, sizeof(cv45at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at11, sizeof(cv45at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at12, sizeof(cv45at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at13, sizeof(cv45at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at14, sizeof(cv45at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at15, sizeof(cv45at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -46 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv46at0, sizeof(cv46at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at1, sizeof(cv46at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at2, sizeof(cv46at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at3, sizeof(cv46at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at4, sizeof(cv46at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at5, sizeof(cv46at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at6, sizeof(cv46at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at7, sizeof(cv46at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at8, sizeof(cv46at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at9, sizeof(cv46at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at10, sizeof(cv46at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at11, sizeof(cv46at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at12, sizeof(cv46at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at13, sizeof(cv46at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at14, sizeof(cv46at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv46at15, sizeof(cv46at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -47 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at1, sizeof(cv47at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at2, sizeof(cv47at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at3, sizeof(cv47at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at4, sizeof(cv47at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at5, sizeof(cv47at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at6, sizeof(cv47at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at7, sizeof(cv47at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at8, sizeof(cv47at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at9, sizeof(cv47at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at10, sizeof(cv47at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at11, sizeof(cv47at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at12, sizeof(cv47at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at13, sizeof(cv47at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at14, sizeof(cv47at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv47at15, sizeof(cv47at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -48 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at1, sizeof(cv48at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at2, sizeof(cv48at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at3, sizeof(cv48at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at4, sizeof(cv48at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at5, sizeof(cv48at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at6, sizeof(cv48at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at7, sizeof(cv48at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at8, sizeof(cv48at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at9, sizeof(cv48at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at10, sizeof(cv48at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at11, sizeof(cv48at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at12, sizeof(cv48at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at13, sizeof(cv48at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at14, sizeof(cv48at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at15, sizeof(cv48at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -49 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at1, sizeof(cv49at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at2, sizeof(cv49at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at3, sizeof(cv49at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at4, sizeof(cv49at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at5, sizeof(cv49at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at6, sizeof(cv49at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at7, sizeof(cv49at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at8, sizeof(cv49at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at9, sizeof(cv49at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at10, sizeof(cv49at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at11, sizeof(cv49at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at12, sizeof(cv49at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at13, sizeof(cv49at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at14, sizeof(cv49at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at15, sizeof(cv49at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -50 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv50at0, sizeof(cv50at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at1, sizeof(cv50at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at2, sizeof(cv50at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at3, sizeof(cv50at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at4, sizeof(cv50at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at5, sizeof(cv50at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at6, sizeof(cv50at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at7, sizeof(cv50at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at8, sizeof(cv50at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at9, sizeof(cv50at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at10, sizeof(cv50at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at11, sizeof(cv50at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at12, sizeof(cv50at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at13, sizeof(cv50at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at14, sizeof(cv50at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at15, sizeof(cv50at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -51 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv51at0, sizeof(cv51at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at1, sizeof(cv51at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at2, sizeof(cv51at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at3, sizeof(cv51at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at4, sizeof(cv51at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at5, sizeof(cv51at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at6, sizeof(cv51at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at7, sizeof(cv51at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at8, sizeof(cv51at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at9, sizeof(cv51at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at10, sizeof(cv51at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at11, sizeof(cv51at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at12, sizeof(cv51at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at13, sizeof(cv51at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at14, sizeof(cv51at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at15, sizeof(cv51at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -52 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv52at0, sizeof(cv52at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at1, sizeof(cv52at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at2, sizeof(cv52at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at3, sizeof(cv52at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at4, sizeof(cv52at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at5, sizeof(cv52at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at6, sizeof(cv52at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at7, sizeof(cv52at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at8, sizeof(cv52at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at9, sizeof(cv52at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at10, sizeof(cv52at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at11, sizeof(cv52at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at12, sizeof(cv52at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at13, sizeof(cv52at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at14, sizeof(cv52at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at15, sizeof(cv52at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -53 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at1, sizeof(cv53at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at2, sizeof(cv53at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at3, sizeof(cv53at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at4, sizeof(cv53at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at5, sizeof(cv53at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at6, sizeof(cv53at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at7, sizeof(cv53at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at8, sizeof(cv53at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at9, sizeof(cv53at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at10, sizeof(cv53at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at11, sizeof(cv53at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at12, sizeof(cv53at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at13, sizeof(cv53at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at14, sizeof(cv53at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at15, sizeof(cv53at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -54 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at1, sizeof(cv54at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at2, sizeof(cv54at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at3, sizeof(cv54at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at4, sizeof(cv54at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at5, sizeof(cv54at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at6, sizeof(cv54at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at7, sizeof(cv54at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at8, sizeof(cv54at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at9, sizeof(cv54at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at10, sizeof(cv54at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at11, sizeof(cv54at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at12, sizeof(cv54at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at13, sizeof(cv54at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at14, sizeof(cv54at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at15, sizeof(cv54at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -55 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at1, sizeof(cv55at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at2, sizeof(cv55at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at3, sizeof(cv55at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at4, sizeof(cv55at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at5, sizeof(cv55at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at6, sizeof(cv55at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at7, sizeof(cv55at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at8, sizeof(cv55at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at9, sizeof(cv55at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at10, sizeof(cv55at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at11, sizeof(cv55at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at12, sizeof(cv55at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at13, sizeof(cv55at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at14, sizeof(cv55at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at15, sizeof(cv55at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -56 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at1, sizeof(cv56at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at2, sizeof(cv56at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at3, sizeof(cv56at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at4, sizeof(cv56at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at5, sizeof(cv56at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at6, sizeof(cv56at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at7, sizeof(cv56at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at8, sizeof(cv56at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at9, sizeof(cv56at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at10, sizeof(cv56at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at11, sizeof(cv56at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at12, sizeof(cv56at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at13, sizeof(cv56at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at14, sizeof(cv56at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at15, sizeof(cv56at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -57 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at1, sizeof(cv57at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at2, sizeof(cv57at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at3, sizeof(cv57at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at4, sizeof(cv57at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at5, sizeof(cv57at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at6, sizeof(cv57at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at7, sizeof(cv57at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at8, sizeof(cv57at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at9, sizeof(cv57at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at10, sizeof(cv57at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at11, sizeof(cv57at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at12, sizeof(cv57at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at13, sizeof(cv57at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at14, sizeof(cv57at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at15, sizeof(cv57at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -58 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv58at0, sizeof(cv58at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at1, sizeof(cv58at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at2, sizeof(cv58at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at3, sizeof(cv58at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at4, sizeof(cv58at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at5, sizeof(cv58at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at6, sizeof(cv58at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at7, sizeof(cv58at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at8, sizeof(cv58at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at9, sizeof(cv58at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at10, sizeof(cv58at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at11, sizeof(cv58at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at12, sizeof(cv58at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at13, sizeof(cv58at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at14, sizeof(cv58at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at15, sizeof(cv58at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -59 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv59at0, sizeof(cv59at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at1, sizeof(cv59at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at2, sizeof(cv59at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at3, sizeof(cv59at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at4, sizeof(cv59at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at5, sizeof(cv59at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at6, sizeof(cv59at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at7, sizeof(cv59at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at8, sizeof(cv59at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at9, sizeof(cv59at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at10, sizeof(cv59at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at11, sizeof(cv59at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at12, sizeof(cv59at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at13, sizeof(cv59at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at14, sizeof(cv59at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at15, sizeof(cv59at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -60 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv60at0, sizeof(cv60at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at1, sizeof(cv60at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at2, sizeof(cv60at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at3, sizeof(cv60at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at4, sizeof(cv60at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at5, sizeof(cv60at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at6, sizeof(cv60at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at7, sizeof(cv60at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at8, sizeof(cv60at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at9, sizeof(cv60at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at10, sizeof(cv60at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at11, sizeof(cv60at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at12, sizeof(cv60at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at13, sizeof(cv60at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at14, sizeof(cv60at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at15, sizeof(cv60at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -61 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv61at0, sizeof(cv61at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at1, sizeof(cv61at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at2, sizeof(cv61at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at3, sizeof(cv61at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at4, sizeof(cv61at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at5, sizeof(cv61at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at6, sizeof(cv61at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at7, sizeof(cv61at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at8, sizeof(cv61at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at9, sizeof(cv61at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at10, sizeof(cv61at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at11, sizeof(cv61at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at12, sizeof(cv61at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at13, sizeof(cv61at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at14, sizeof(cv61at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv61at15, sizeof(cv61at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -62 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv62at0, sizeof(cv62at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at1, sizeof(cv62at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at2, sizeof(cv62at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at3, sizeof(cv62at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at4, sizeof(cv62at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at5, sizeof(cv62at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at6, sizeof(cv62at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at7, sizeof(cv62at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at8, sizeof(cv62at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at9, sizeof(cv62at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at10, sizeof(cv62at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at11, sizeof(cv62at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at12, sizeof(cv62at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at13, sizeof(cv62at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at14, sizeof(cv62at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at15, sizeof(cv62at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -63 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv63at0, sizeof(cv63at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at1, sizeof(cv63at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at2, sizeof(cv63at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at3, sizeof(cv63at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at4, sizeof(cv63at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at5, sizeof(cv63at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at6, sizeof(cv63at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at7, sizeof(cv63at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at8, sizeof(cv63at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at9, sizeof(cv63at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at10, sizeof(cv63at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at11, sizeof(cv63at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at12, sizeof(cv63at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at13, sizeof(cv63at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at14, sizeof(cv63at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv63at15, sizeof(cv63at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -64 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv64at0, sizeof(cv64at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at1, sizeof(cv64at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at2, sizeof(cv64at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at3, sizeof(cv64at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at4, sizeof(cv64at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at5, sizeof(cv64at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at6, sizeof(cv64at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at7, sizeof(cv64at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at8, sizeof(cv64at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at9, sizeof(cv64at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at10, sizeof(cv64at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at11, sizeof(cv64at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at12, sizeof(cv64at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at13, sizeof(cv64at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at14, sizeof(cv64at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at15, sizeof(cv64at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -65 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv65at0, sizeof(cv65at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at1, sizeof(cv65at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at2, sizeof(cv65at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at3, sizeof(cv65at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at4, sizeof(cv65at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at5, sizeof(cv65at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at6, sizeof(cv65at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at7, sizeof(cv65at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at8, sizeof(cv65at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at9, sizeof(cv65at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at10, sizeof(cv65at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at11, sizeof(cv65at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at12, sizeof(cv65at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at13, sizeof(cv65at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at14, sizeof(cv65at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at15, sizeof(cv65at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -66 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv66at0, sizeof(cv66at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at1, sizeof(cv66at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at2, sizeof(cv66at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at3, sizeof(cv66at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at4, sizeof(cv66at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at5, sizeof(cv66at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at6, sizeof(cv66at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at7, sizeof(cv66at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at8, sizeof(cv66at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at9, sizeof(cv66at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at10, sizeof(cv66at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at11, sizeof(cv66at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at12, sizeof(cv66at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at13, sizeof(cv66at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at14, sizeof(cv66at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at15, sizeof(cv66at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -67 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv67at0, sizeof(cv67at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at1, sizeof(cv67at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at2, sizeof(cv67at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at3, sizeof(cv67at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at4, sizeof(cv67at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at5, sizeof(cv67at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at6, sizeof(cv67at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at7, sizeof(cv67at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at8, sizeof(cv67at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at9, sizeof(cv67at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at10, sizeof(cv67at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at11, sizeof(cv67at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at12, sizeof(cv67at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at13, sizeof(cv67at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at14, sizeof(cv67at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at15, sizeof(cv67at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -68 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at1, sizeof(cv68at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at2, sizeof(cv68at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at3, sizeof(cv68at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at4, sizeof(cv68at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at5, sizeof(cv68at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at6, sizeof(cv68at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at7, sizeof(cv68at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at8, sizeof(cv68at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at9, sizeof(cv68at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at10, sizeof(cv68at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at11, sizeof(cv68at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at12, sizeof(cv68at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at13, sizeof(cv68at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at14, sizeof(cv68at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at15, sizeof(cv68at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -69 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at1, sizeof(cv69at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at2, sizeof(cv69at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at3, sizeof(cv69at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at4, sizeof(cv69at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at5, sizeof(cv69at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at6, sizeof(cv69at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at7, sizeof(cv69at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at8, sizeof(cv69at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at9, sizeof(cv69at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at10, sizeof(cv69at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at11, sizeof(cv69at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at12, sizeof(cv69at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at13, sizeof(cv69at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at14, sizeof(cv69at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at15, sizeof(cv69at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -70 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv70at0, sizeof(cv70at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at1, sizeof(cv70at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at2, sizeof(cv70at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at3, sizeof(cv70at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at4, sizeof(cv70at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at5, sizeof(cv70at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at6, sizeof(cv70at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at7, sizeof(cv70at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at8, sizeof(cv70at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at9, sizeof(cv70at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at10, sizeof(cv70at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at11, sizeof(cv70at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at12, sizeof(cv70at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at13, sizeof(cv70at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at14, sizeof(cv70at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at15, sizeof(cv70at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -71 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv71at0, sizeof(cv71at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at1, sizeof(cv71at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at2, sizeof(cv71at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at3, sizeof(cv71at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at4, sizeof(cv71at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at5, sizeof(cv71at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at6, sizeof(cv71at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at7, sizeof(cv71at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at8, sizeof(cv71at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at9, sizeof(cv71at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at10, sizeof(cv71at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at11, sizeof(cv71at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at12, sizeof(cv71at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at13, sizeof(cv71at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at14, sizeof(cv71at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at15, sizeof(cv71at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -72 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv72at0, sizeof(cv72at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at1, sizeof(cv72at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at2, sizeof(cv72at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at3, sizeof(cv72at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at4, sizeof(cv72at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at5, sizeof(cv72at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at6, sizeof(cv72at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at7, sizeof(cv72at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at8, sizeof(cv72at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at9, sizeof(cv72at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at10, sizeof(cv72at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at11, sizeof(cv72at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at12, sizeof(cv72at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at13, sizeof(cv72at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at14, sizeof(cv72at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv72at15, sizeof(cv72at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -73 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv73at0, sizeof(cv73at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at1, sizeof(cv73at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at2, sizeof(cv73at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at3, sizeof(cv73at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at4, sizeof(cv73at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at5, sizeof(cv73at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at6, sizeof(cv73at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at7, sizeof(cv73at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at8, sizeof(cv73at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at9, sizeof(cv73at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at10, sizeof(cv73at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at11, sizeof(cv73at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at12, sizeof(cv73at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at13, sizeof(cv73at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at14, sizeof(cv73at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv73at15, sizeof(cv73at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -74 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv74at0, sizeof(cv74at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at1, sizeof(cv74at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at2, sizeof(cv74at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at3, sizeof(cv74at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at4, sizeof(cv74at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at5, sizeof(cv74at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at6, sizeof(cv74at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at7, sizeof(cv74at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at8, sizeof(cv74at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at9, sizeof(cv74at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at10, sizeof(cv74at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at11, sizeof(cv74at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at12, sizeof(cv74at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at13, sizeof(cv74at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at14, sizeof(cv74at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv74at15, sizeof(cv74at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -75 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv75at0, sizeof(cv75at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at1, sizeof(cv75at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at2, sizeof(cv75at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at3, sizeof(cv75at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at4, sizeof(cv75at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at5, sizeof(cv75at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at6, sizeof(cv75at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at7, sizeof(cv75at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at8, sizeof(cv75at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at9, sizeof(cv75at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at10, sizeof(cv75at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at11, sizeof(cv75at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at12, sizeof(cv75at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at13, sizeof(cv75at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at14, sizeof(cv75at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv75at15, sizeof(cv75at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -76 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv76at0, sizeof(cv76at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at1, sizeof(cv76at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at2, sizeof(cv76at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at3, sizeof(cv76at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at4, sizeof(cv76at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at5, sizeof(cv76at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at6, sizeof(cv76at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at7, sizeof(cv76at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at8, sizeof(cv76at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at9, sizeof(cv76at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at10, sizeof(cv76at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at11, sizeof(cv76at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at12, sizeof(cv76at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at13, sizeof(cv76at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at14, sizeof(cv76at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv76at15, sizeof(cv76at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -77 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv77at0, sizeof(cv77at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at1, sizeof(cv77at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at2, sizeof(cv77at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at3, sizeof(cv77at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at4, sizeof(cv77at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at5, sizeof(cv77at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at6, sizeof(cv77at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at7, sizeof(cv77at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at8, sizeof(cv77at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at9, sizeof(cv77at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at10, sizeof(cv77at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at11, sizeof(cv77at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at12, sizeof(cv77at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at13, sizeof(cv77at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at14, sizeof(cv77at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv77at15, sizeof(cv77at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -78 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv78at0, sizeof(cv78at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at1, sizeof(cv78at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at2, sizeof(cv78at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at3, sizeof(cv78at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at4, sizeof(cv78at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at5, sizeof(cv78at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at6, sizeof(cv78at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at7, sizeof(cv78at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at8, sizeof(cv78at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at9, sizeof(cv78at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at10, sizeof(cv78at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at11, sizeof(cv78at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at12, sizeof(cv78at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at13, sizeof(cv78at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at14, sizeof(cv78at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv78at15, sizeof(cv78at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -79 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv79at0, sizeof(cv79at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at1, sizeof(cv79at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at2, sizeof(cv79at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at3, sizeof(cv79at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at4, sizeof(cv79at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at5, sizeof(cv79at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at6, sizeof(cv79at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at7, sizeof(cv79at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at8, sizeof(cv79at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at9, sizeof(cv79at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at10, sizeof(cv79at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at11, sizeof(cv79at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at12, sizeof(cv79at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at13, sizeof(cv79at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at14, sizeof(cv79at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv79at15, sizeof(cv79at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -80 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv80at0, sizeof(cv80at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at1, sizeof(cv80at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at2, sizeof(cv80at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at3, sizeof(cv80at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at4, sizeof(cv80at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at5, sizeof(cv80at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at6, sizeof(cv80at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at7, sizeof(cv80at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at8, sizeof(cv80at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at9, sizeof(cv80at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at10, sizeof(cv80at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at11, sizeof(cv80at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at12, sizeof(cv80at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at13, sizeof(cv80at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at14, sizeof(cv80at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at15, sizeof(cv80at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -81 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv81at0, sizeof(cv81at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at1, sizeof(cv81at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at2, sizeof(cv81at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at3, sizeof(cv81at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at4, sizeof(cv81at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at5, sizeof(cv81at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at6, sizeof(cv81at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at7, sizeof(cv81at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at8, sizeof(cv81at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at9, sizeof(cv81at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at10, sizeof(cv81at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at11, sizeof(cv81at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at12, sizeof(cv81at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at13, sizeof(cv81at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at14, sizeof(cv81at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at15, sizeof(cv81at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -82 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv82at0, sizeof(cv82at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at1, sizeof(cv82at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at2, sizeof(cv82at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at3, sizeof(cv82at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at4, sizeof(cv82at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at5, sizeof(cv82at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at6, sizeof(cv82at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at7, sizeof(cv82at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at8, sizeof(cv82at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at9, sizeof(cv82at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at10, sizeof(cv82at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at11, sizeof(cv82at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at12, sizeof(cv82at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at13, sizeof(cv82at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at14, sizeof(cv82at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv82at15, sizeof(cv82at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -83 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv83at0, sizeof(cv83at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at1, sizeof(cv83at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at2, sizeof(cv83at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at3, sizeof(cv83at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at4, sizeof(cv83at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at5, sizeof(cv83at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at6, sizeof(cv83at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at7, sizeof(cv83at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at8, sizeof(cv83at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at9, sizeof(cv83at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at10, sizeof(cv83at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at11, sizeof(cv83at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at12, sizeof(cv83at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at13, sizeof(cv83at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at14, sizeof(cv83at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv83at15, sizeof(cv83at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -84 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv84at0, sizeof(cv84at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at1, sizeof(cv84at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at2, sizeof(cv84at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at3, sizeof(cv84at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at4, sizeof(cv84at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at5, sizeof(cv84at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at6, sizeof(cv84at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at7, sizeof(cv84at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at8, sizeof(cv84at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at9, sizeof(cv84at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at10, sizeof(cv84at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at11, sizeof(cv84at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at12, sizeof(cv84at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at13, sizeof(cv84at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at14, sizeof(cv84at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv84at15, sizeof(cv84at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -85 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv85at0, sizeof(cv85at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at1, sizeof(cv85at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at2, sizeof(cv85at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at3, sizeof(cv85at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at4, sizeof(cv85at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at5, sizeof(cv85at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at6, sizeof(cv85at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at7, sizeof(cv85at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at8, sizeof(cv85at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at9, sizeof(cv85at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at10, sizeof(cv85at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at11, sizeof(cv85at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at12, sizeof(cv85at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at13, sizeof(cv85at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at14, sizeof(cv85at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv85at15, sizeof(cv85at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -86 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv86at0, sizeof(cv86at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at1, sizeof(cv86at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at2, sizeof(cv86at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at3, sizeof(cv86at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at4, sizeof(cv86at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at5, sizeof(cv86at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at6, sizeof(cv86at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at7, sizeof(cv86at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at8, sizeof(cv86at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at9, sizeof(cv86at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at10, sizeof(cv86at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at11, sizeof(cv86at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at12, sizeof(cv86at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at13, sizeof(cv86at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at14, sizeof(cv86at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at15, sizeof(cv86at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -87 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv87at0, sizeof(cv87at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at1, sizeof(cv87at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at2, sizeof(cv87at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at3, sizeof(cv87at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at4, sizeof(cv87at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at5, sizeof(cv87at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at6, sizeof(cv87at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at7, sizeof(cv87at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at8, sizeof(cv87at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at9, sizeof(cv87at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at10, sizeof(cv87at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at11, sizeof(cv87at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at12, sizeof(cv87at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at13, sizeof(cv87at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at14, sizeof(cv87at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at15, sizeof(cv87at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -88 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv88at0, sizeof(cv88at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at1, sizeof(cv88at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at2, sizeof(cv88at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at3, sizeof(cv88at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at4, sizeof(cv88at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at5, sizeof(cv88at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at6, sizeof(cv88at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at7, sizeof(cv88at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at8, sizeof(cv88at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at9, sizeof(cv88at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at10, sizeof(cv88at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at11, sizeof(cv88at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at12, sizeof(cv88at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at13, sizeof(cv88at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at14, sizeof(cv88at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at15, sizeof(cv88at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -89 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv89at0, sizeof(cv89at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at1, sizeof(cv89at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at2, sizeof(cv89at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at3, sizeof(cv89at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at4, sizeof(cv89at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at5, sizeof(cv89at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at6, sizeof(cv89at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at7, sizeof(cv89at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at8, sizeof(cv89at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at9, sizeof(cv89at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at10, sizeof(cv89at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at11, sizeof(cv89at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at12, sizeof(cv89at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at13, sizeof(cv89at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at14, sizeof(cv89at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at15, sizeof(cv89at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -90 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv90at0, sizeof(cv90at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at1, sizeof(cv90at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at2, sizeof(cv90at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at3, sizeof(cv90at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at4, sizeof(cv90at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at5, sizeof(cv90at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at6, sizeof(cv90at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at7, sizeof(cv90at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at8, sizeof(cv90at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at9, sizeof(cv90at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at10, sizeof(cv90at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at11, sizeof(cv90at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at12, sizeof(cv90at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at13, sizeof(cv90at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at14, sizeof(cv90at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv90at15, sizeof(cv90at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -91 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv91at0, sizeof(cv91at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at1, sizeof(cv91at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at2, sizeof(cv91at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at3, sizeof(cv91at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at4, sizeof(cv91at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at5, sizeof(cv91at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at6, sizeof(cv91at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at7, sizeof(cv91at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at8, sizeof(cv91at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at9, sizeof(cv91at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at10, sizeof(cv91at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at11, sizeof(cv91at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at12, sizeof(cv91at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at13, sizeof(cv91at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at14, sizeof(cv91at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv91at15, sizeof(cv91at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -92 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv92at0, sizeof(cv92at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at1, sizeof(cv92at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at2, sizeof(cv92at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at3, sizeof(cv92at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at4, sizeof(cv92at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at5, sizeof(cv92at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at6, sizeof(cv92at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at7, sizeof(cv92at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at8, sizeof(cv92at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at9, sizeof(cv92at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at10, sizeof(cv92at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at11, sizeof(cv92at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at12, sizeof(cv92at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at13, sizeof(cv92at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at14, sizeof(cv92at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv92at15, sizeof(cv92at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -93 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv93at0, sizeof(cv93at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at1, sizeof(cv93at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at2, sizeof(cv93at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at3, sizeof(cv93at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at4, sizeof(cv93at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at5, sizeof(cv93at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at6, sizeof(cv93at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at7, sizeof(cv93at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at8, sizeof(cv93at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at9, sizeof(cv93at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at10, sizeof(cv93at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at11, sizeof(cv93at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at12, sizeof(cv93at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at13, sizeof(cv93at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at14, sizeof(cv93at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv93at15, sizeof(cv93at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -94 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv94at0, sizeof(cv94at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at1, sizeof(cv94at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at2, sizeof(cv94at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at3, sizeof(cv94at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at4, sizeof(cv94at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at5, sizeof(cv94at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at6, sizeof(cv94at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at7, sizeof(cv94at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at8, sizeof(cv94at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at9, sizeof(cv94at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at10, sizeof(cv94at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at11, sizeof(cv94at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at12, sizeof(cv94at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at13, sizeof(cv94at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at14, sizeof(cv94at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at15, sizeof(cv94at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -95 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv95at0, sizeof(cv95at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at1, sizeof(cv95at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at2, sizeof(cv95at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at3, sizeof(cv95at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at4, sizeof(cv95at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at5, sizeof(cv95at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at6, sizeof(cv95at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at7, sizeof(cv95at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at8, sizeof(cv95at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at9, sizeof(cv95at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at10, sizeof(cv95at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at11, sizeof(cv95at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at12, sizeof(cv95at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at13, sizeof(cv95at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at14, sizeof(cv95at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at15, sizeof(cv95at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -96 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv96at0, sizeof(cv96at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at1, sizeof(cv96at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at2, sizeof(cv96at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at3, sizeof(cv96at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at4, sizeof(cv96at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at5, sizeof(cv96at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at6, sizeof(cv96at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at7, sizeof(cv96at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at8, sizeof(cv96at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at9, sizeof(cv96at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at10, sizeof(cv96at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at11, sizeof(cv96at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at12, sizeof(cv96at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at13, sizeof(cv96at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at14, sizeof(cv96at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv96at15, sizeof(cv96at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -97 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv97at0, sizeof(cv97at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at1, sizeof(cv97at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at2, sizeof(cv97at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at3, sizeof(cv97at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at4, sizeof(cv97at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at5, sizeof(cv97at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at6, sizeof(cv97at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at7, sizeof(cv97at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at8, sizeof(cv97at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at9, sizeof(cv97at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at10, sizeof(cv97at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at11, sizeof(cv97at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at12, sizeof(cv97at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at13, sizeof(cv97at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at14, sizeof(cv97at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv97at15, sizeof(cv97at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -98 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv98at0, sizeof(cv98at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at1, sizeof(cv98at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at2, sizeof(cv98at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at3, sizeof(cv98at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at4, sizeof(cv98at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at5, sizeof(cv98at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at6, sizeof(cv98at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at7, sizeof(cv98at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at8, sizeof(cv98at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at9, sizeof(cv98at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at10, sizeof(cv98at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at11, sizeof(cv98at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at12, sizeof(cv98at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at13, sizeof(cv98at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at14, sizeof(cv98at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv98at15, sizeof(cv98at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -99 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv99at0, sizeof(cv99at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at1, sizeof(cv99at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at2, sizeof(cv99at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at3, sizeof(cv99at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at4, sizeof(cv99at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at5, sizeof(cv99at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at6, sizeof(cv99at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at7, sizeof(cv99at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at8, sizeof(cv99at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at9, sizeof(cv99at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at10, sizeof(cv99at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at11, sizeof(cv99at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at12, sizeof(cv99at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at13, sizeof(cv99at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at14, sizeof(cv99at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv99at15, sizeof(cv99at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -100 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv100at0, sizeof(cv100at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at1, sizeof(cv100at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at2, sizeof(cv100at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at3, sizeof(cv100at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at4, sizeof(cv100at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at5, sizeof(cv100at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at6, sizeof(cv100at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at7, sizeof(cv100at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at8, sizeof(cv100at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at9, sizeof(cv100at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at10, sizeof(cv100at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at11, sizeof(cv100at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at12, sizeof(cv100at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at13, sizeof(cv100at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at14, sizeof(cv100at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at15, sizeof(cv100at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -101 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv101at0, sizeof(cv101at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at1, sizeof(cv101at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at2, sizeof(cv101at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at3, sizeof(cv101at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at4, sizeof(cv101at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at5, sizeof(cv101at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at6, sizeof(cv101at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at7, sizeof(cv101at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at8, sizeof(cv101at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at9, sizeof(cv101at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at10, sizeof(cv101at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at11, sizeof(cv101at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at12, sizeof(cv101at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at13, sizeof(cv101at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at14, sizeof(cv101at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at15, sizeof(cv101at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -102 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv102at0, sizeof(cv102at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at1, sizeof(cv102at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at2, sizeof(cv102at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at3, sizeof(cv102at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at4, sizeof(cv102at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at5, sizeof(cv102at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at6, sizeof(cv102at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at7, sizeof(cv102at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at8, sizeof(cv102at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at9, sizeof(cv102at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at10, sizeof(cv102at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at11, sizeof(cv102at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at12, sizeof(cv102at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at13, sizeof(cv102at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at14, sizeof(cv102at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv102at15, sizeof(cv102at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -103 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv103at0, sizeof(cv103at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at1, sizeof(cv103at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at2, sizeof(cv103at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at3, sizeof(cv103at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at4, sizeof(cv103at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at5, sizeof(cv103at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at6, sizeof(cv103at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at7, sizeof(cv103at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at8, sizeof(cv103at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at9, sizeof(cv103at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at10, sizeof(cv103at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at11, sizeof(cv103at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at12, sizeof(cv103at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at13, sizeof(cv103at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at14, sizeof(cv103at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv103at15, sizeof(cv103at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -104 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv104at0, sizeof(cv104at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at1, sizeof(cv104at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at2, sizeof(cv104at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at3, sizeof(cv104at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at4, sizeof(cv104at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at5, sizeof(cv104at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at6, sizeof(cv104at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at7, sizeof(cv104at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at8, sizeof(cv104at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at9, sizeof(cv104at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at10, sizeof(cv104at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at11, sizeof(cv104at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at12, sizeof(cv104at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at13, sizeof(cv104at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at14, sizeof(cv104at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv104at15, sizeof(cv104at15)/sizeof(int)), i);   
    else if( diags[i]->id() == -105 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv105at0, sizeof(cv105at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at1, sizeof(cv105at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at2, sizeof(cv105at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at3, sizeof(cv105at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at4, sizeof(cv105at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at5, sizeof(cv105at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at6, sizeof(cv105at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at7, sizeof(cv105at7)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at8, sizeof(cv105at8)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at9, sizeof(cv105at9)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at10, sizeof(cv105at10)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at11, sizeof(cv105at11)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at12, sizeof(cv105at12)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at13, sizeof(cv105at13)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at14, sizeof(cv105at14)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv105at15, sizeof(cv105at15)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEgg2gggg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEgg2gggg

  static const ColourLines diag1[16] = { 
    ColourLines("1 3 5 9, -1 2, -2 -3 -4 -6, 6 -7, 7 4 -5 -8, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 5 9, 6 -7, 7 4 -5 -8, 8 -9"), 
    ColourLines("1 3 5 9, -1 2, -2 -3 -4 -7, 6 4 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 5 9, 6 4 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 3 5 8, -1 2, -2 -3 -4 -6, 6 -7, 7 4 -5 -9, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 5 8, 6 -7, 7 4 -5 -9, -8 9"), 
    ColourLines("1 3 5 8, -1 2, -2 -3 -4 -7, 6 4 -5 -9, -6 7, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 5 8, 6 4 -5 -9, -6 7, -8 9"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -5 -8, 6 -7, -6 -4 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -5 -8, 2 3 4 7, 6 -7, -6 -4 5 9, 8 -9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -8, -6 7, -7 -4 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -5 -8, 2 3 4 6, -6 7, -7 -4 5 9, 8 -9"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -5 -9, 6 -7, -6 -4 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -5 -9, 2 3 4 7, 6 -7, -6 -4 5 8, -8 9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -9, -6 7, -7 -4 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -5 -9, 2 3 4 6, -6 7, -7 -4 5 8, -8 9")
  }; 
  static const ColourLines diag2[16] = { 
    ColourLines("1 3 5 9, -1 2, -2 -3 -4 -6, 6 -8, 7 -9, -7 -5 4 8"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 5 9, 6 -8, 7 -9, -7 -5 4 8"), 
    ColourLines("1 3 5 9, -1 2, -2 -3 -4 -8, 6 4 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 5 9, 6 4 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 3 5 7, -1 2, -2 -3 -4 -6, 6 -8, -7 9, 8 4 -5 -9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 5 7, 6 -8, -7 9, 8 4 -5 -9"), 
    ColourLines("1 3 5 7, -1 2, -2 -3 -4 -8, 6 4 -5 -9, -6 8, -7 9"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 5 7, 6 4 -5 -9, -6 8, -7 9"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -5 -7, 6 -8, -6 -4 5 9, 7 -9"), 
    ColourLines("1 -2, -1 -3 -5 -7, 2 3 4 8, 6 -8, -6 -4 5 9, 7 -9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -7, -6 8, 7 -9, -8 -4 5 9"), 
    ColourLines("1 -2, -1 -3 -5 -7, 2 3 4 6, -6 8, 7 -9, -8 -4 5 9"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -5 -9, 6 -8, -6 -4 5 7, -7 9"), 
    ColourLines("1 -2, -1 -3 -5 -9, 2 3 4 8, 6 -8, -6 -4 5 7, -7 9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -9, -6 8, 7 5 -4 -8, -7 9"), 
    ColourLines("1 -2, -1 -3 -5 -9, 2 3 4 6, -6 8, 7 5 -4 -8, -7 9")
  }; 
  static const ColourLines diag3[16] = { 
    ColourLines("1 3 5 8, -1 2, -2 -3 -4 -6, 6 -9, 7 -8, -7 -5 4 9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 5 8, 6 -9, 7 -8, -7 -5 4 9"), 
    ColourLines("1 3 5 8, -1 2, -2 -3 -4 -9, 6 4 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 5 8, 6 4 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 3 5 7, -1 2, -2 -3 -4 -6, 6 -9, -7 8, -8 -5 4 9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 5 7, 6 -9, -7 8, -8 -5 4 9"), 
    ColourLines("1 3 5 7, -1 2, -2 -3 -4 -9, 6 4 -5 -8, -6 9, -7 8"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 5 7, 6 4 -5 -8, -6 9, -7 8"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -5 -7, 6 -9, -6 -4 5 8, 7 -8"), 
    ColourLines("1 -2, -1 -3 -5 -7, 2 3 4 9, 6 -9, -6 -4 5 8, 7 -8"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -7, -6 9, 7 -8, 8 5 -4 -9"), 
    ColourLines("1 -2, -1 -3 -5 -7, 2 3 4 6, -6 9, 7 -8, 8 5 -4 -9"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -5 -8, 6 -9, -6 -4 5 7, -7 8"), 
    ColourLines("1 -2, -1 -3 -5 -8, 2 3 4 9, 6 -9, -6 -4 5 7, -7 8"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -8, -6 9, 7 5 -4 -9, -7 8"), 
    ColourLines("1 -2, -1 -3 -5 -8, 2 3 4 6, -6 9, 7 5 -4 -9, -7 8")
  }; 
  static const ColourLines diag4[16] = { 
    ColourLines("1 2 5 9, -1 -6, 4 -3 -2 6, -4 -7, 7 3 -5 -8, 8 -9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -7, -6 2 5 9, 7 3 -5 -8, 8 -9"), 
    ColourLines("1 2 5 9, -1 -6, 4 7, -4 3 -5 -8, 6 -2 -3 -7, 8 -9"), 
    ColourLines("1 6, -1 -2 -3 -7, 4 7, -4 3 -5 -8, -6 2 5 9, 8 -9"), 
    ColourLines("1 2 5 8, -1 -6, 4 -3 -2 6, -4 -7, 7 3 -5 -9, -8 9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -7, -6 2 5 8, 7 3 -5 -9, -8 9"), 
    ColourLines("1 2 5 8, -1 -6, 4 7, -4 3 -5 -9, 6 -2 -3 -7, -8 9"), 
    ColourLines("1 6, -1 -2 -3 -7, 4 7, -4 3 -5 -9, -6 2 5 8, -8 9"), 
    ColourLines("1 2 3 7, -1 -6, 4 -3 5 9, -4 -7, 6 -2 -5 -8, 8 -9"), 
    ColourLines("1 6, -1 -2 -5 -8, 4 -3 5 9, -4 -7, -6 2 3 7, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 7, 6 -2 -5 -8, -7 -3 5 9, 8 -9"), 
    ColourLines("1 6, -1 -2 -5 -8, 4 7, -4 3 2 -6, -7 -3 5 9, 8 -9"), 
    ColourLines("1 2 3 7, -1 -6, 4 -3 5 8, -4 -7, 6 -2 -5 -9, -8 9"), 
    ColourLines("1 6, -1 -2 -5 -9, 4 -3 5 8, -4 -7, -6 2 3 7, -8 9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 7, 6 -2 -5 -9, -7 -3 5 8, -8 9"), 
    ColourLines("1 6, -1 -2 -5 -9, 4 7, -4 3 2 -6, -7 -3 5 8, -8 9")
  }; 
  static const ColourLines diag5[16] = { 
    ColourLines("1 2 5 9, -1 -6, 4 -3 -2 6, -4 -8, 7 -9, -7 -5 3 8"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -8, -6 2 5 9, 7 -9, -7 -5 3 8"), 
    ColourLines("1 2 5 9, -1 -6, 4 8, -4 3 -5 -7, 6 -2 -3 -8, 7 -9"), 
    ColourLines("1 6, -1 -2 -3 -8, 4 8, -4 3 -5 -7, -6 2 5 9, 7 -9"), 
    ColourLines("1 2 5 7, -1 -6, 4 -3 -2 6, -4 -8, -7 9, 8 3 -5 -9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -8, -6 2 5 7, -7 9, 8 3 -5 -9"), 
    ColourLines("1 2 5 7, -1 -6, 4 8, -4 3 -5 -9, 6 -2 -3 -8, -7 9"), 
    ColourLines("1 6, -1 -2 -3 -8, 4 8, -4 3 -5 -9, -6 2 5 7, -7 9"), 
    ColourLines("1 2 3 8, -1 -6, 4 -3 5 9, -4 -8, 6 -2 -5 -7, 7 -9"), 
    ColourLines("1 6, -1 -2 -5 -7, 4 -3 5 9, -4 -8, -6 2 3 8, 7 -9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 8, 6 -2 -5 -7, 7 -9, -8 -3 5 9"), 
    ColourLines("1 6, -1 -2 -5 -7, 4 8, -4 3 2 -6, 7 -9, -8 -3 5 9"), 
    ColourLines("1 2 3 8, -1 -6, 4 -3 5 7, -4 -8, 6 -2 -5 -9, -7 9"), 
    ColourLines("1 6, -1 -2 -5 -9, 4 -3 5 7, -4 -8, -6 2 3 8, -7 9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 8, 6 -2 -5 -9, 7 5 -3 -8, -7 9"), 
    ColourLines("1 6, -1 -2 -5 -9, 4 8, -4 3 2 -6, 7 5 -3 -8, -7 9")
  }; 
  static const ColourLines diag6[16] = { 
    ColourLines("1 2 5 8, -1 -6, 4 -3 -2 6, -4 -9, 7 -8, -7 -5 3 9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -9, -6 2 5 8, 7 -8, -7 -5 3 9"), 
    ColourLines("1 2 5 8, -1 -6, 4 9, -4 3 -5 -7, 6 -2 -3 -9, 7 -8"), 
    ColourLines("1 6, -1 -2 -3 -9, 4 9, -4 3 -5 -7, -6 2 5 8, 7 -8"), 
    ColourLines("1 2 5 7, -1 -6, 4 -3 -2 6, -4 -9, -7 8, -8 -5 3 9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -9, -6 2 5 7, -7 8, -8 -5 3 9"), 
    ColourLines("1 2 5 7, -1 -6, 4 9, -4 3 -5 -8, 6 -2 -3 -9, -7 8"), 
    ColourLines("1 6, -1 -2 -3 -9, 4 9, -4 3 -5 -8, -6 2 5 7, -7 8"), 
    ColourLines("1 2 3 9, -1 -6, 4 -3 5 8, -4 -9, 6 -2 -5 -7, 7 -8"), 
    ColourLines("1 6, -1 -2 -5 -7, 4 -3 5 8, -4 -9, -6 2 3 9, 7 -8"), 
    ColourLines("1 2 3 -4, -1 -6, 4 9, 6 -2 -5 -7, 7 -8, 8 5 -3 -9"), 
    ColourLines("1 6, -1 -2 -5 -7, 4 9, -4 3 2 -6, 7 -8, 8 5 -3 -9"), 
    ColourLines("1 2 3 9, -1 -6, 4 -3 5 7, -4 -9, 6 -2 -5 -8, -7 8"), 
    ColourLines("1 6, -1 -2 -5 -8, 4 -3 5 7, -4 -9, -6 2 3 9, -7 8"), 
    ColourLines("1 2 3 -4, -1 -6, 4 9, 6 -2 -5 -8, 7 5 -3 -9, -7 8"), 
    ColourLines("1 6, -1 -2 -5 -8, 4 9, -4 3 2 -6, 7 5 -3 -9, -7 8")
  }; 
  static const ColourLines diag7[16] = { 
    ColourLines("1 2 5 9, -1 -7, 4 -3 -2 7, -4 -6, 6 3 -5 -8, 8 -9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -6, 6 3 -5 -8, -7 2 5 9, 8 -9"), 
    ColourLines("1 2 5 9, -1 -7, 4 6, -4 3 -5 -8, -6 -3 -2 7, 8 -9"), 
    ColourLines("1 7, -1 -2 -3 -6, 4 6, -4 3 -5 -8, -7 2 5 9, 8 -9"), 
    ColourLines("1 2 5 8, -1 -7, 4 -3 -2 7, -4 -6, 6 3 -5 -9, -8 9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -6, 6 3 -5 -9, -7 2 5 8, -8 9"), 
    ColourLines("1 2 5 8, -1 -7, 4 6, -4 3 -5 -9, -6 -3 -2 7, -8 9"), 
    ColourLines("1 7, -1 -2 -3 -6, 4 6, -4 3 -5 -9, -7 2 5 8, -8 9"), 
    ColourLines("1 2 3 6, -1 -7, 4 -3 5 9, -4 -6, 7 -2 -5 -8, 8 -9"), 
    ColourLines("1 7, -1 -2 -5 -8, 4 -3 5 9, -4 -6, 6 3 2 -7, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 6, -6 -3 5 9, 7 -2 -5 -8, 8 -9"), 
    ColourLines("1 7, -1 -2 -5 -8, 4 6, -4 3 2 -7, -6 -3 5 9, 8 -9"), 
    ColourLines("1 2 3 6, -1 -7, 4 -3 5 8, -4 -6, 7 -2 -5 -9, -8 9"), 
    ColourLines("1 7, -1 -2 -5 -9, 4 -3 5 8, -4 -6, 6 3 2 -7, -8 9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 6, -6 -3 5 8, 7 -2 -5 -9, -8 9"), 
    ColourLines("1 7, -1 -2 -5 -9, 4 6, -4 3 2 -7, -6 -3 5 8, -8 9")
  }; 
  static const ColourLines diag8[16] = { 
    ColourLines("1 2 5 9, -1 -7, 4 -3 -2 7, -4 -8, 6 -9, -6 -5 3 8"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -8, 6 -9, -6 -5 3 8, -7 2 5 9"), 
    ColourLines("1 2 5 9, -1 -7, 4 8, -4 3 -5 -6, 6 -9, 7 -2 -3 -8"), 
    ColourLines("1 7, -1 -2 -3 -8, 4 8, -4 3 -5 -6, 6 -9, -7 2 5 9"), 
    ColourLines("1 2 5 6, -1 -7, 4 -3 -2 7, -4 -8, -6 9, 8 3 -5 -9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -8, 6 5 2 -7, -6 9, 8 3 -5 -9"), 
    ColourLines("1 2 5 6, -1 -7, 4 8, -4 3 -5 -9, -6 9, 7 -2 -3 -8"), 
    ColourLines("1 7, -1 -2 -3 -8, 4 8, -4 3 -5 -9, 6 5 2 -7, -6 9"), 
    ColourLines("1 2 3 8, -1 -7, 4 -3 5 9, -4 -8, 6 -9, -6 -5 -2 7"), 
    ColourLines("1 7, -1 -2 -5 -6, 4 -3 5 9, -4 -8, 6 -9, -7 2 3 8"), 
    ColourLines("1 2 3 -4, -1 -7, 4 8, 6 -9, -6 -5 -2 7, -8 -3 5 9"), 
    ColourLines("1 7, -1 -2 -5 -6, 4 8, -4 3 2 -7, 6 -9, -8 -3 5 9"), 
    ColourLines("1 2 3 8, -1 -7, 4 -3 5 6, -4 -8, -6 9, 7 -2 -5 -9"), 
    ColourLines("1 7, -1 -2 -5 -9, 4 -3 5 6, -4 -8, -6 9, -7 2 3 8"), 
    ColourLines("1 2 3 -4, -1 -7, 4 8, 6 5 -3 -8, -6 9, 7 -2 -5 -9"), 
    ColourLines("1 7, -1 -2 -5 -9, 4 8, -4 3 2 -7, 6 5 -3 -8, -6 9")
  }; 
  static const ColourLines diag9[16] = { 
    ColourLines("1 2 5 8, -1 -7, 4 -3 -2 7, -4 -9, 6 -8, -6 -5 3 9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -9, 6 -8, -6 -5 3 9, -7 2 5 8"), 
    ColourLines("1 2 5 8, -1 -7, 4 9, -4 3 -5 -6, 6 -8, 7 -2 -3 -9"), 
    ColourLines("1 7, -1 -2 -3 -9, 4 9, -4 3 -5 -6, 6 -8, -7 2 5 8"), 
    ColourLines("1 2 5 6, -1 -7, 4 -3 -2 7, -4 -9, -6 8, -8 -5 3 9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -9, 6 5 2 -7, -6 8, -8 -5 3 9"), 
    ColourLines("1 2 5 6, -1 -7, 4 9, -4 3 -5 -8, -6 8, 7 -2 -3 -9"), 
    ColourLines("1 7, -1 -2 -3 -9, 4 9, -4 3 -5 -8, 6 5 2 -7, -6 8"), 
    ColourLines("1 2 3 9, -1 -7, 4 -3 5 8, -4 -9, 6 -8, -6 -5 -2 7"), 
    ColourLines("1 7, -1 -2 -5 -6, 4 -3 5 8, -4 -9, 6 -8, -7 2 3 9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 9, 6 -8, -6 -5 -2 7, 8 5 -3 -9"), 
    ColourLines("1 7, -1 -2 -5 -6, 4 9, -4 3 2 -7, 6 -8, 8 5 -3 -9"), 
    ColourLines("1 2 3 9, -1 -7, 4 -3 5 6, -4 -9, -6 8, 7 -2 -5 -8"), 
    ColourLines("1 7, -1 -2 -5 -8, 4 -3 5 6, -4 -9, -6 8, -7 2 3 9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 9, 6 5 -3 -9, -6 8, 7 -2 -5 -8"), 
    ColourLines("1 7, -1 -2 -5 -8, 4 9, -4 3 2 -7, 6 5 -3 -9, -6 8")
  }; 
  static const ColourLines diag10[16] = { 
    ColourLines("1 2 5 9, -1 -8, 4 -3 -2 8, -4 -6, 6 3 -5 -7, 7 -9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -6, 6 3 -5 -7, 7 -9, -8 2 5 9"), 
    ColourLines("1 2 5 9, -1 -8, 4 6, -4 3 -5 -7, -6 -3 -2 8, 7 -9"), 
    ColourLines("1 8, -1 -2 -3 -6, 4 6, -4 3 -5 -7, 7 -9, -8 2 5 9"), 
    ColourLines("1 2 5 7, -1 -8, 4 -3 -2 8, -4 -6, 6 3 -5 -9, -7 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -6, 6 3 -5 -9, 7 5 2 -8, -7 9"), 
    ColourLines("1 2 5 7, -1 -8, 4 6, -4 3 -5 -9, -6 -3 -2 8, -7 9"), 
    ColourLines("1 8, -1 -2 -3 -6, 4 6, -4 3 -5 -9, 7 5 2 -8, -7 9"), 
    ColourLines("1 2 3 6, -1 -8, 4 -3 5 9, -4 -6, 7 -9, -7 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -7, 4 -3 5 9, -4 -6, 6 3 2 -8, 7 -9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 6, -6 -3 5 9, 7 -9, -7 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -7, 4 6, -4 3 2 -8, -6 -3 5 9, 7 -9"), 
    ColourLines("1 2 3 6, -1 -8, 4 -3 5 7, -4 -6, -7 9, 8 -2 -5 -9"), 
    ColourLines("1 8, -1 -2 -5 -9, 4 -3 5 7, -4 -6, 6 3 2 -8, -7 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 6, -6 -3 5 7, -7 9, 8 -2 -5 -9"), 
    ColourLines("1 8, -1 -2 -5 -9, 4 6, -4 3 2 -8, -6 -3 5 7, -7 9")
  }; 
  static const ColourLines diag11[16] = { 
    ColourLines("1 2 5 9, -1 -8, 4 -3 -2 8, -4 -7, 6 -9, -6 -5 3 7"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -7, 6 -9, -6 -5 3 7, -8 2 5 9"), 
    ColourLines("1 2 5 9, -1 -8, 4 7, -4 3 -5 -6, 6 -9, -7 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -7, 4 7, -4 3 -5 -6, 6 -9, -8 2 5 9"), 
    ColourLines("1 2 5 6, -1 -8, 4 -3 -2 8, -4 -7, -6 9, 7 3 -5 -9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -7, 6 5 2 -8, -6 9, 7 3 -5 -9"), 
    ColourLines("1 2 5 6, -1 -8, 4 7, -4 3 -5 -9, -6 9, -7 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -7, 4 7, -4 3 -5 -9, 6 5 2 -8, -6 9"), 
    ColourLines("1 2 3 7, -1 -8, 4 -3 5 9, -4 -7, 6 -9, -6 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -6, 4 -3 5 9, -4 -7, 6 -9, 7 3 2 -8"), 
    ColourLines("1 2 3 -4, -1 -8, 4 7, 6 -9, -6 -5 -2 8, -7 -3 5 9"), 
    ColourLines("1 8, -1 -2 -5 -6, 4 7, -4 3 2 -8, 6 -9, -7 -3 5 9"), 
    ColourLines("1 2 3 7, -1 -8, 4 -3 5 6, -4 -7, -6 9, 8 -2 -5 -9"), 
    ColourLines("1 8, -1 -2 -5 -9, 4 -3 5 6, -4 -7, -6 9, 7 3 2 -8"), 
    ColourLines("1 2 3 -4, -1 -8, 4 7, 6 5 -3 -7, -6 9, 8 -2 -5 -9"), 
    ColourLines("1 8, -1 -2 -5 -9, 4 7, -4 3 2 -8, 6 5 -3 -7, -6 9")
  }; 
  static const ColourLines diag12[16] = { 
    ColourLines("1 2 5 7, -1 -8, 4 -3 -2 8, -4 -9, 6 -7, -6 -5 3 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -9, 6 -7, -6 -5 3 9, 7 5 2 -8"), 
    ColourLines("1 2 5 7, -1 -8, 4 9, -4 3 -5 -6, 6 -7, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 4 9, -4 3 -5 -6, 6 -7, 7 5 2 -8"), 
    ColourLines("1 2 5 6, -1 -8, 4 -3 -2 8, -4 -9, -6 7, -7 -5 3 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -9, 6 5 2 -8, -6 7, -7 -5 3 9"), 
    ColourLines("1 2 5 6, -1 -8, 4 9, -4 3 -5 -7, -6 7, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 4 9, -4 3 -5 -7, 6 5 2 -8, -6 7"), 
    ColourLines("1 2 3 9, -1 -8, 4 -3 5 7, -4 -9, 6 -7, -6 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -6, 4 -3 5 7, -4 -9, 6 -7, -8 2 3 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 9, 6 -7, -6 -5 -2 8, 7 5 -3 -9"), 
    ColourLines("1 8, -1 -2 -5 -6, 4 9, -4 3 2 -8, 6 -7, 7 5 -3 -9"), 
    ColourLines("1 2 3 9, -1 -8, 4 -3 5 6, -4 -9, -6 7, -7 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -7, 4 -3 5 6, -4 -9, -6 7, -8 2 3 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 9, 6 5 -3 -9, -6 7, -7 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -7, 4 9, -4 3 2 -8, 6 5 -3 -9, -6 7")
  }; 
  static const ColourLines diag13[16] = { 
    ColourLines("1 2 5 8, -1 -9, 4 -3 -2 9, -4 -6, 6 3 -5 -7, 7 -8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -6, 6 3 -5 -7, 7 -8, 8 5 2 -9"), 
    ColourLines("1 2 5 8, -1 -9, 4 6, -4 3 -5 -7, -6 -3 -2 9, 7 -8"), 
    ColourLines("1 9, -1 -2 -3 -6, 4 6, -4 3 -5 -7, 7 -8, 8 5 2 -9"), 
    ColourLines("1 2 5 7, -1 -9, 4 -3 -2 9, -4 -6, 6 3 -5 -8, -7 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -6, 6 3 -5 -8, 7 5 2 -9, -7 8"), 
    ColourLines("1 2 5 7, -1 -9, 4 6, -4 3 -5 -8, -6 -3 -2 9, -7 8"), 
    ColourLines("1 9, -1 -2 -3 -6, 4 6, -4 3 -5 -8, 7 5 2 -9, -7 8"), 
    ColourLines("1 2 3 6, -1 -9, 4 -3 5 8, -4 -6, 7 -8, -7 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -7, 4 -3 5 8, -4 -6, 6 3 2 -9, 7 -8"), 
    ColourLines("1 2 3 -4, -1 -9, 4 6, -6 -3 5 8, 7 -8, -7 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -7, 4 6, -4 3 2 -9, -6 -3 5 8, 7 -8"), 
    ColourLines("1 2 3 6, -1 -9, 4 -3 5 7, -4 -6, -7 8, -8 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -8, 4 -3 5 7, -4 -6, 6 3 2 -9, -7 8"), 
    ColourLines("1 2 3 -4, -1 -9, 4 6, -6 -3 5 7, -7 8, -8 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -8, 4 6, -4 3 2 -9, -6 -3 5 7, -7 8")
  }; 
  static const ColourLines diag14[16] = { 
    ColourLines("1 2 5 8, -1 -9, 4 -3 -2 9, -4 -7, 6 -8, -6 -5 3 7"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -7, 6 -8, -6 -5 3 7, 8 5 2 -9"), 
    ColourLines("1 2 5 8, -1 -9, 4 7, -4 3 -5 -6, 6 -8, -7 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -7, 4 7, -4 3 -5 -6, 6 -8, 8 5 2 -9"), 
    ColourLines("1 2 5 6, -1 -9, 4 -3 -2 9, -4 -7, -6 8, 7 3 -5 -8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -7, 6 5 2 -9, -6 8, 7 3 -5 -8"), 
    ColourLines("1 2 5 6, -1 -9, 4 7, -4 3 -5 -8, -6 8, -7 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -7, 4 7, -4 3 -5 -8, 6 5 2 -9, -6 8"), 
    ColourLines("1 2 3 7, -1 -9, 4 -3 5 8, -4 -7, 6 -8, -6 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -6, 4 -3 5 8, -4 -7, 6 -8, 7 3 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 7, 6 -8, -6 -5 -2 9, -7 -3 5 8"), 
    ColourLines("1 9, -1 -2 -5 -6, 4 7, -4 3 2 -9, 6 -8, -7 -3 5 8"), 
    ColourLines("1 2 3 7, -1 -9, 4 -3 5 6, -4 -7, -6 8, -8 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -8, 4 -3 5 6, -4 -7, -6 8, 7 3 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 7, 6 5 -3 -7, -6 8, -8 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -8, 4 7, -4 3 2 -9, 6 5 -3 -7, -6 8")
  }; 
  static const ColourLines diag15[16] = { 
    ColourLines("1 2 5 7, -1 -9, 4 -3 -2 9, -4 -8, 6 -7, -6 -5 3 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -8, 6 -7, -6 -5 3 8, 7 5 2 -9"), 
    ColourLines("1 2 5 7, -1 -9, 4 8, -4 3 -5 -6, 6 -7, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 4 8, -4 3 -5 -6, 6 -7, 7 5 2 -9"), 
    ColourLines("1 2 5 6, -1 -9, 4 -3 -2 9, -4 -8, -6 7, -7 -5 3 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -8, 6 5 2 -9, -6 7, -7 -5 3 8"), 
    ColourLines("1 2 5 6, -1 -9, 4 8, -4 3 -5 -7, -6 7, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 4 8, -4 3 -5 -7, 6 5 2 -9, -6 7"), 
    ColourLines("1 2 3 8, -1 -9, 4 -3 5 7, -4 -8, 6 -7, -6 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -6, 4 -3 5 7, -4 -8, 6 -7, 8 3 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 8, 6 -7, -6 -5 -2 9, 7 5 -3 -8"), 
    ColourLines("1 9, -1 -2 -5 -6, 4 8, -4 3 2 -9, 6 -7, 7 5 -3 -8"), 
    ColourLines("1 2 3 8, -1 -9, 4 -3 5 6, -4 -8, -6 7, -7 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -7, 4 -3 5 6, -4 -8, -6 7, 8 3 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 8, 6 5 -3 -8, -6 7, -7 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -7, 4 8, -4 3 2 -9, 6 5 -3 -8, -6 7")
  }; 
  static const ColourLines diag16[16] = { 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -9, 6 -7, -6 -5 9, 7 5 4 -8"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 8, 6 -7, -6 -5 9, 7 5 4 -8"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -9, 6 5 4 -8, -6 7, -7 -5 9"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 8, 6 5 4 -8, -6 7, -7 -5 9"), 
    ColourLines("1 3 4 5 7, -1 2, -2 -3 -8, 6 -7, -6 -5 9, 8 -4 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 7, 6 -7, -6 -5 9, 8 -4 -9"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -8, -6 7, -7 -5 9, 8 -4 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 6, -6 7, -7 -5 9, 8 -4 -9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -6, 6 -7, 7 5 -9, -8 4 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -6, 2 3 8, 6 -7, 7 5 -9, -8 4 9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -7, 6 5 -9, -6 7, -8 4 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 8, 6 5 -9, -6 7, -8 4 9"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -8, 6 -7, -6 -5 -4 8, 7 5 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 9, 6 -7, -6 -5 -4 8, 7 5 -9"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -8, 6 5 -9, -6 7, -7 -5 -4 8"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 9, 6 5 -9, -6 7, -7 -5 -4 8")
  }; 
  static const ColourLines diag17[16] = { 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -8, 6 -7, -6 -5 8, 7 5 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 9, 6 -7, -6 -5 8, 7 5 4 -9"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -8, 6 5 4 -9, -6 7, -7 -5 8"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 9, 6 5 4 -9, -6 7, -7 -5 8"), 
    ColourLines("1 3 4 5 7, -1 2, -2 -3 -9, 6 -7, -6 -5 8, -8 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 7, 6 -7, -6 -5 8, -8 -4 9"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -9, -6 7, -7 -5 8, -8 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 6, -6 7, -7 -5 8, -8 -4 9"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -6, 6 -7, 7 5 -8, 8 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -6, 2 3 9, 6 -7, 7 5 -8, 8 4 -9"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -7, 6 5 -8, -6 7, 8 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 9, 6 5 -8, -6 7, 8 4 -9"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -9, 6 -7, -6 -5 -4 9, 7 5 -8"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 8, 6 -7, -6 -5 -4 9, 7 5 -8"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -9, 6 5 -8, -6 7, -7 -5 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 8, 6 5 -8, -6 7, -7 -5 -4 9")
  }; 
  static const ColourLines diag18[16] = { 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -9, 6 -8, -6 -5 9, -7 4 5 8"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 7, 6 -8, -6 -5 9, -7 4 5 8"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -9, 6 5 4 -7, -6 8, -8 -5 9"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 7, 6 5 4 -7, -6 8, -8 -5 9"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -7, 6 -8, -6 -5 9, 7 -4 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 8, 6 -8, -6 -5 9, 7 -4 -9"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -7, -6 8, 7 -4 -9, -8 -5 9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 6, -6 8, 7 -4 -9, -8 -5 9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -5 -6, 6 -8, -7 4 9, 8 5 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -6, 2 3 7, 6 -8, -7 4 9, 8 5 -9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -5 -8, 6 5 -9, -6 8, -7 4 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 7, 6 5 -9, -6 8, -7 4 9"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -7, 6 -8, -6 -5 -4 7, 8 5 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 9, 6 -8, -6 -5 -4 7, 8 5 -9"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -7, 6 5 -9, -6 8, 7 -4 -5 -8"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 9, 6 5 -9, -6 8, 7 -4 -5 -8")
  }; 
  static const ColourLines diag19[16] = { 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -7, 6 -8, -6 -5 7, 8 5 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 9, 6 -8, -6 -5 7, 8 5 4 -9"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -7, 6 5 4 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 9, 6 5 4 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -9, 6 -8, -6 -5 7, -7 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 8, 6 -8, -6 -5 7, -7 -4 9"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -9, -6 8, 7 -5 -8, -7 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 6, -6 8, 7 -5 -8, -7 -4 9"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -6, 6 -8, 7 4 -9, -7 5 8"), 
    ColourLines("1 -2, -1 -3 -4 -5 -6, 2 3 9, 6 -8, 7 4 -9, -7 5 8"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -8, 6 5 -7, -6 8, 7 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 9, 6 5 -7, -6 8, 7 4 -9"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -9, 6 -8, -6 -5 -4 9, -7 5 8"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 7, 6 -8, -6 -5 -4 9, -7 5 8"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -9, 6 5 -7, -6 8, -8 -5 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 7, 6 5 -7, -6 8, -8 -5 -4 9")
  }; 
  static const ColourLines diag20[16] = { 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -8, 6 -9, -6 -5 8, -7 4 5 9"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 7, 6 -9, -6 -5 8, -7 4 5 9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -8, 6 5 4 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 7, 6 5 4 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -7, 6 -9, -6 -5 8, 7 -4 -8"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 9, 6 -9, -6 -5 8, 7 -4 -8"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -7, -6 9, 7 -4 -8, 8 -5 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 6, -6 9, 7 -4 -8, 8 -5 -9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -5 -6, 6 -9, -7 4 8, -8 5 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -6, 2 3 7, 6 -9, -7 4 8, -8 5 9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -5 -9, 6 5 -8, -6 9, -7 4 8"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 7, 6 5 -8, -6 9, -7 4 8"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -7, 6 -9, -6 -5 -4 7, -8 5 9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 8, 6 -9, -6 -5 -4 7, -8 5 9"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -7, 6 5 -8, -6 9, 7 -4 -5 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 8, 6 5 -8, -6 9, 7 -4 -5 -9")
  }; 
  static const ColourLines diag21[16] = { 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -7, 6 -9, -6 -5 7, -8 4 5 9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 8, 6 -9, -6 -5 7, -8 4 5 9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -7, 6 5 4 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 8, 6 5 4 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -8, 6 -9, -6 -5 7, -7 -4 8"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 9, 6 -9, -6 -5 7, -7 -4 8"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -8, -6 9, 7 -5 -9, -7 -4 8"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 6, -6 9, 7 -5 -9, -7 -4 8"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -6, 6 -9, 7 4 -8, -7 5 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -6, 2 3 8, 6 -9, 7 4 -8, -7 5 9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -9, 6 5 -7, -6 9, 7 4 -8"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 8, 6 5 -7, -6 9, 7 4 -8"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -8, 6 -9, -6 -5 -4 8, -7 5 9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 7, 6 -9, -6 -5 -4 8, -7 5 9"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -8, 6 5 -7, -6 9, 8 -4 -5 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 7, 6 5 -7, -6 9, 8 -4 -5 -9")
  }; 
  static const ColourLines diag22[16] = { 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -9, -6 4 5 8, 7 -8, -7 -5 9"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 6, -6 4 5 8, 7 -8, -7 -5 9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -9, -6 4 5 7, -7 8, -8 -5 9"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 6, -6 4 5 7, -7 8, -8 -5 9"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -6, 6 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 5 8, 6 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 3 4 5 7, -1 2, -2 -3 -6, 6 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 5 7, 6 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -7, -6 4 9, 7 -8, 8 5 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 6, -6 4 9, 7 -8, 8 5 -9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -8, -6 4 9, 7 5 -9, -7 8"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 6, -6 4 9, 7 5 -9, -7 8"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -6, 6 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 9, 6 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -6, 6 -4 -5 -8, 7 5 -9, -7 8"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 9, 6 -4 -5 -8, 7 5 -9, -7 8")
  }; 
  static const ColourLines diag23[16] = { 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -6, 6 -5 -7, 7 -8, 8 5 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 9, 6 -5 -7, 7 -8, 8 5 4 -9"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -6, 6 -5 -8, 7 5 4 -9, -7 8"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 9, 6 -5 -8, 7 5 4 -9, -7 8"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -9, 6 -5 -7, -6 -4 9, 7 -8"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 8, 6 -5 -7, -6 -4 9, 7 -8"), 
    ColourLines("1 3 4 5 7, -1 2, -2 -3 -9, 6 -5 -8, -6 -4 9, -7 8"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 7, 6 -5 -8, -6 -4 9, -7 8"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -7, 6 4 -9, -6 5 8, 7 -8"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 9, 6 4 -9, -6 5 8, 7 -8"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -8, 6 4 -9, -6 5 7, -7 8"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 9, 6 4 -9, -6 5 7, -7 8"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -9, -6 5 8, 7 -8, -7 -5 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 6, -6 5 8, 7 -8, -7 -5 -4 9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -9, -6 5 7, -7 8, -8 -5 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 6, -6 5 7, -7 8, -8 -5 -4 9")
  }; 
  static const ColourLines diag24[16] = { 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -8, -6 4 5 9, 7 -9, -7 -5 8"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 6, -6 4 5 9, 7 -9, -7 -5 8"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -8, -6 4 5 7, -7 9, 8 -5 -9"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 6, -6 4 5 7, -7 9, 8 -5 -9"), 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -6, 6 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 5 9, 6 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 3 4 5 7, -1 2, -2 -3 -6, 6 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 5 7, 6 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -7, -6 4 8, 7 -9, -8 5 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 6, -6 4 8, 7 -9, -8 5 9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -9, -6 4 8, 7 5 -8, -7 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 6, -6 4 8, 7 5 -8, -7 9"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -6, 6 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 8, 6 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -6, 6 -4 -5 -9, 7 5 -8, -7 9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 8, 6 -4 -5 -9, 7 5 -8, -7 9")
  }; 
  static const ColourLines diag25[16] = { 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -6, 6 -5 -7, 7 -9, -8 4 5 9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 8, 6 -5 -7, 7 -9, -8 4 5 9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -6, 6 -5 -9, 7 5 4 -8, -7 9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 8, 6 -5 -9, 7 5 4 -8, -7 9"), 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -8, 6 -5 -7, -6 -4 8, 7 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 9, 6 -5 -7, -6 -4 8, 7 -9"), 
    ColourLines("1 3 4 5 7, -1 2, -2 -3 -8, 6 -5 -9, -6 -4 8, -7 9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 7, 6 -5 -9, -6 -4 8, -7 9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -7, 6 4 -8, -6 5 9, 7 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 8, 6 4 -8, -6 5 9, 7 -9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -9, 6 4 -8, -6 5 7, -7 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 8, 6 4 -8, -6 5 7, -7 9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -8, -6 5 9, 7 -9, -7 -5 -4 8"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 6, -6 5 9, 7 -9, -7 -5 -4 8"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -8, -6 5 7, -7 9, 8 -4 -5 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 6, -6 5 7, -7 9, 8 -4 -5 -9")
  }; 
  static const ColourLines diag26[16] = { 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -7, -6 4 5 9, 7 -5 -8, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 6, -6 4 5 9, 7 -5 -8, 8 -9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -7, -6 4 5 8, 7 -5 -9, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 6, -6 4 5 8, 7 -5 -9, -8 9"), 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -6, 6 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 5 9, 6 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -6, 6 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 5 8, 6 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -8, -6 4 7, -7 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 6, -6 4 7, -7 5 9, 8 -9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -9, -6 4 7, -7 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 6, -6 4 7, -7 5 8, -8 9"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -6, 6 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 7, 6 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 3 4 7, -1 2, -2 -3 -6, 6 -4 -5 -9, -7 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -6, 2 3 4 7, 6 -4 -5 -9, -7 5 8, -8 9")
  }; 
  static const ColourLines diag27[16] = { 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -6, 6 -5 -8, -7 4 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 7, 6 -5 -8, -7 4 5 9, 8 -9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -6, 6 -5 -9, -7 4 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -6, 2 3 7, 6 -5 -9, -7 4 5 8, -8 9"), 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -7, 6 -5 -8, -6 -4 7, 8 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 9, 6 -5 -8, -6 -4 7, 8 -9"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -7, 6 -5 -9, -6 -4 7, -8 9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 8, 6 -5 -9, -6 -4 7, -8 9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -5 -8, 6 4 -7, -6 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 7, 6 4 -7, -6 5 9, 8 -9"), 
    ColourLines("1 3 7, -1 2, -2 -3 -4 -5 -9, 6 4 -7, -6 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 7, 6 4 -7, -6 5 8, -8 9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -7, -6 5 9, 7 -4 -5 -8, 8 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 6, -6 5 9, 7 -4 -5 -8, 8 -9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -7, -6 5 8, 7 -4 -5 -9, -8 9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 6, -6 5 8, 7 -4 -5 -9, -8 9")
  }; 
  static const ColourLines diag28[16] = { 
    ColourLines("1 2 8, -1 -6, 5 -4 9, -5 -7, 6 -2 -3 -9, 7 4 3 -8"), 
    ColourLines("1 6, -1 -2 -3 -9, 5 -4 9, -5 -7, -6 2 8, 7 4 3 -8"), 
    ColourLines("1 2 8, -1 -6, 5 7, -5 4 3 -8, 6 -2 -3 -9, -7 -4 9"), 
    ColourLines("1 6, -1 -2 -3 -9, 5 7, -5 4 3 -8, -6 2 8, -7 -4 9"), 
    ColourLines("1 2 3 4 7, -1 -6, 5 -4 9, -5 -7, 6 -2 -8, 8 -3 -9"), 
    ColourLines("1 6, -1 -2 -8, 5 -4 9, -5 -7, -6 2 3 4 7, 8 -3 -9"), 
    ColourLines("1 2 3 4 -5, -1 -6, 5 7, 6 -2 -8, -7 -4 9, 8 -3 -9"), 
    ColourLines("1 6, -1 -2 -8, 5 7, -5 4 3 2 -6, -7 -4 9, 8 -3 -9"), 
    ColourLines("1 2 8, -1 -6, 5 -4 -3 -2 6, -5 -7, 7 4 -9, -8 3 9"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -7, -6 2 8, 7 4 -9, -8 3 9"), 
    ColourLines("1 2 8, -1 -6, 5 7, -5 4 -9, 6 -2 -3 -4 -7, -8 3 9"), 
    ColourLines("1 6, -1 -2 -3 -4 -7, 5 7, -5 4 -9, -6 2 8, -8 3 9"), 
    ColourLines("1 2 3 9, -1 -6, 5 -4 -3 8, -5 -7, 6 -2 -8, 7 4 -9"), 
    ColourLines("1 6, -1 -2 -8, 5 -4 -3 8, -5 -7, -6 2 3 9, 7 4 -9"), 
    ColourLines("1 2 3 9, -1 -6, 5 7, -5 4 -9, 6 -2 -8, -7 -4 -3 8"), 
    ColourLines("1 6, -1 -2 -8, 5 7, -5 4 -9, -6 2 3 9, -7 -4 -3 8")
  }; 
  static const ColourLines diag29[16] = { 
    ColourLines("1 2 9, -1 -6, 5 -4 8, -5 -7, 6 -2 -3 -8, 7 4 3 -9"), 
    ColourLines("1 6, -1 -2 -3 -8, 5 -4 8, -5 -7, -6 2 9, 7 4 3 -9"), 
    ColourLines("1 2 9, -1 -6, 5 7, -5 4 3 -9, 6 -2 -3 -8, -7 -4 8"), 
    ColourLines("1 6, -1 -2 -3 -8, 5 7, -5 4 3 -9, -6 2 9, -7 -4 8"), 
    ColourLines("1 2 3 4 7, -1 -6, 5 -4 8, -5 -7, 6 -2 -9, -8 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 5 -4 8, -5 -7, -6 2 3 4 7, -8 -3 9"), 
    ColourLines("1 2 3 4 -5, -1 -6, 5 7, 6 -2 -9, -7 -4 8, -8 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 5 7, -5 4 3 2 -6, -7 -4 8, -8 -3 9"), 
    ColourLines("1 2 9, -1 -6, 5 -4 -3 -2 6, -5 -7, 7 4 -8, 8 3 -9"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -7, -6 2 9, 7 4 -8, 8 3 -9"), 
    ColourLines("1 2 9, -1 -6, 5 7, -5 4 -8, 6 -2 -3 -4 -7, 8 3 -9"), 
    ColourLines("1 6, -1 -2 -3 -4 -7, 5 7, -5 4 -8, -6 2 9, 8 3 -9"), 
    ColourLines("1 2 3 8, -1 -6, 5 -4 -3 9, -5 -7, 6 -2 -9, 7 4 -8"), 
    ColourLines("1 6, -1 -2 -9, 5 -4 -3 9, -5 -7, -6 2 3 8, 7 4 -8"), 
    ColourLines("1 2 3 8, -1 -6, 5 7, -5 4 -8, 6 -2 -9, -7 -4 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 5 7, -5 4 -8, -6 2 3 8, -7 -4 -3 9")
  }; 
  static const ColourLines diag30[16] = { 
    ColourLines("1 2 7, -1 -6, 5 -4 9, -5 -8, 6 -2 -3 -9, -7 3 4 8"), 
    ColourLines("1 6, -1 -2 -3 -9, 5 -4 9, -5 -8, -6 2 7, -7 3 4 8"), 
    ColourLines("1 2 7, -1 -6, 5 8, -5 4 3 -7, 6 -2 -3 -9, -8 -4 9"), 
    ColourLines("1 6, -1 -2 -3 -9, 5 8, -5 4 3 -7, -6 2 7, -8 -4 9"), 
    ColourLines("1 2 3 4 8, -1 -6, 5 -4 9, -5 -8, 6 -2 -7, 7 -3 -9"), 
    ColourLines("1 6, -1 -2 -7, 5 -4 9, -5 -8, -6 2 3 4 8, 7 -3 -9"), 
    ColourLines("1 2 3 4 -5, -1 -6, 5 8, 6 -2 -7, 7 -3 -9, -8 -4 9"), 
    ColourLines("1 6, -1 -2 -7, 5 8, -5 4 3 2 -6, 7 -3 -9, -8 -4 9"), 
    ColourLines("1 2 7, -1 -6, 5 -4 -3 -2 6, -5 -8, -7 3 9, 8 4 -9"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -8, -6 2 7, -7 3 9, 8 4 -9"), 
    ColourLines("1 2 7, -1 -6, 5 8, -5 4 -9, 6 -2 -3 -4 -8, -7 3 9"), 
    ColourLines("1 6, -1 -2 -3 -4 -8, 5 8, -5 4 -9, -6 2 7, -7 3 9"), 
    ColourLines("1 2 3 9, -1 -6, 5 -4 -3 7, -5 -8, 6 -2 -7, 8 4 -9"), 
    ColourLines("1 6, -1 -2 -7, 5 -4 -3 7, -5 -8, -6 2 3 9, 8 4 -9"), 
    ColourLines("1 2 3 9, -1 -6, 5 8, -5 4 -9, 6 -2 -7, 7 -3 -4 -8"), 
    ColourLines("1 6, -1 -2 -7, 5 8, -5 4 -9, -6 2 3 9, 7 -3 -4 -8")
  }; 
  static const ColourLines diag31[16] = { 
    ColourLines("1 2 9, -1 -6, 5 -4 7, -5 -8, 6 -2 -3 -7, 8 4 3 -9"), 
    ColourLines("1 6, -1 -2 -3 -7, 5 -4 7, -5 -8, -6 2 9, 8 4 3 -9"), 
    ColourLines("1 2 9, -1 -6, 5 8, -5 4 3 -9, 6 -2 -3 -7, 7 -4 -8"), 
    ColourLines("1 6, -1 -2 -3 -7, 5 8, -5 4 3 -9, -6 2 9, 7 -4 -8"), 
    ColourLines("1 2 3 4 8, -1 -6, 5 -4 7, -5 -8, 6 -2 -9, -7 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 5 -4 7, -5 -8, -6 2 3 4 8, -7 -3 9"), 
    ColourLines("1 2 3 4 -5, -1 -6, 5 8, 6 -2 -9, 7 -4 -8, -7 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 5 8, -5 4 3 2 -6, 7 -4 -8, -7 -3 9"), 
    ColourLines("1 2 9, -1 -6, 5 -4 -3 -2 6, -5 -8, 7 3 -9, -7 4 8"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -8, -6 2 9, 7 3 -9, -7 4 8"), 
    ColourLines("1 2 9, -1 -6, 5 8, -5 4 -7, 6 -2 -3 -4 -8, 7 3 -9"), 
    ColourLines("1 6, -1 -2 -3 -4 -8, 5 8, -5 4 -7, -6 2 9, 7 3 -9"), 
    ColourLines("1 2 3 7, -1 -6, 5 -4 -3 9, -5 -8, 6 -2 -9, -7 4 8"), 
    ColourLines("1 6, -1 -2 -9, 5 -4 -3 9, -5 -8, -6 2 3 7, -7 4 8"), 
    ColourLines("1 2 3 7, -1 -6, 5 8, -5 4 -7, 6 -2 -9, -8 -4 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 5 8, -5 4 -7, -6 2 3 7, -8 -4 -3 9")
  }; 
  static const ColourLines diag32[16] = { 
    ColourLines("1 2 7, -1 -6, 5 -4 8, -5 -9, 6 -2 -3 -8, -7 3 4 9"), 
    ColourLines("1 6, -1 -2 -3 -8, 5 -4 8, -5 -9, -6 2 7, -7 3 4 9"), 
    ColourLines("1 2 7, -1 -6, 5 9, -5 4 3 -7, 6 -2 -3 -8, 8 -4 -9"), 
    ColourLines("1 6, -1 -2 -3 -8, 5 9, -5 4 3 -7, -6 2 7, 8 -4 -9"), 
    ColourLines("1 2 3 4 9, -1 -6, 5 -4 8, -5 -9, 6 -2 -7, 7 -3 -8"), 
    ColourLines("1 6, -1 -2 -7, 5 -4 8, -5 -9, -6 2 3 4 9, 7 -3 -8"), 
    ColourLines("1 2 3 4 -5, -1 -6, 5 9, 6 -2 -7, 7 -3 -8, 8 -4 -9"), 
    ColourLines("1 6, -1 -2 -7, 5 9, -5 4 3 2 -6, 7 -3 -8, 8 -4 -9"), 
    ColourLines("1 2 7, -1 -6, 5 -4 -3 -2 6, -5 -9, -7 3 8, -8 4 9"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -9, -6 2 7, -7 3 8, -8 4 9"), 
    ColourLines("1 2 7, -1 -6, 5 9, -5 4 -8, 6 -2 -3 -4 -9, -7 3 8"), 
    ColourLines("1 6, -1 -2 -3 -4 -9, 5 9, -5 4 -8, -6 2 7, -7 3 8"), 
    ColourLines("1 2 3 8, -1 -6, 5 -4 -3 7, -5 -9, 6 -2 -7, -8 4 9"), 
    ColourLines("1 6, -1 -2 -7, 5 -4 -3 7, -5 -9, -6 2 3 8, -8 4 9"), 
    ColourLines("1 2 3 8, -1 -6, 5 9, -5 4 -8, 6 -2 -7, 7 -3 -4 -9"), 
    ColourLines("1 6, -1 -2 -7, 5 9, -5 4 -8, -6 2 3 8, 7 -3 -4 -9")
  }; 
  static const ColourLines diag33[16] = { 
    ColourLines("1 2 8, -1 -6, 5 -4 7, -5 -9, 6 -2 -3 -7, -8 3 4 9"), 
    ColourLines("1 6, -1 -2 -3 -7, 5 -4 7, -5 -9, -6 2 8, -8 3 4 9"), 
    ColourLines("1 2 8, -1 -6, 5 9, -5 4 3 -8, 6 -2 -3 -7, 7 -4 -9"), 
    ColourLines("1 6, -1 -2 -3 -7, 5 9, -5 4 3 -8, -6 2 8, 7 -4 -9"), 
    ColourLines("1 2 3 4 9, -1 -6, 5 -4 7, -5 -9, 6 -2 -8, -7 -3 8"), 
    ColourLines("1 6, -1 -2 -8, 5 -4 7, -5 -9, -6 2 3 4 9, -7 -3 8"), 
    ColourLines("1 2 3 4 -5, -1 -6, 5 9, 6 -2 -8, 7 -4 -9, -7 -3 8"), 
    ColourLines("1 6, -1 -2 -8, 5 9, -5 4 3 2 -6, 7 -4 -9, -7 -3 8"), 
    ColourLines("1 2 8, -1 -6, 5 -4 -3 -2 6, -5 -9, 7 3 -8, -7 4 9"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -9, -6 2 8, 7 3 -8, -7 4 9"), 
    ColourLines("1 2 8, -1 -6, 5 9, -5 4 -7, 6 -2 -3 -4 -9, 7 3 -8"), 
    ColourLines("1 6, -1 -2 -3 -4 -9, 5 9, -5 4 -7, -6 2 8, 7 3 -8"), 
    ColourLines("1 2 3 7, -1 -6, 5 -4 -3 8, -5 -9, 6 -2 -8, -7 4 9"), 
    ColourLines("1 6, -1 -2 -8, 5 -4 -3 8, -5 -9, -6 2 3 7, -7 4 9"), 
    ColourLines("1 2 3 7, -1 -6, 5 9, -5 4 -7, 6 -2 -8, 8 -3 -4 -9"), 
    ColourLines("1 6, -1 -2 -8, 5 9, -5 4 -7, -6 2 3 7, 8 -3 -4 -9")
  }; 
  static const ColourLines diag34[16] = { 
    ColourLines("1 2 -3, -1 -6, 3 4 5 8, 6 -2 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 6, -1 -2 -4 -9, 3 4 5 8, -3 2 -6, 7 -8, -7 -5 9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 5 7, 6 -2 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 6, -1 -2 -4 -9, 3 4 5 7, -3 2 -6, -7 8, -8 -5 9"), 
    ColourLines("1 2 4 5 8, -1 -6, 3 -2 6, -3 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -9, -6 2 4 5 8, 7 -8, -7 -5 9"), 
    ColourLines("1 2 4 5 7, -1 -6, 3 -2 6, -3 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -9, -6 2 4 5 7, -7 8, -8 -5 9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 9, 6 -2 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 6, -1 -2 -4 -5 -7, 3 4 9, -3 2 -6, 7 -8, 8 5 -9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 9, 6 -2 -4 -5 -8, 7 5 -9, -7 8"), 
    ColourLines("1 6, -1 -2 -4 -5 -8, 3 4 9, -3 2 -6, 7 5 -9, -7 8"), 
    ColourLines("1 2 4 9, -1 -6, 3 -2 6, -3 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -7, -6 2 4 9, 7 -8, 8 5 -9"), 
    ColourLines("1 2 4 9, -1 -6, 3 -2 6, -3 -4 -5 -8, 7 5 -9, -7 8"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -8, -6 2 4 9, 7 5 -9, -7 8")
  }; 
  static const ColourLines diag35[16] = { 
    ColourLines("1 2 9, -1 -6, 4 -3 -2 6, -4 -5 -7, 7 -8, 8 5 3 -9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -7, -6 2 9, 7 -8, 8 5 3 -9"), 
    ColourLines("1 2 9, -1 -6, 4 -3 -2 6, -4 -5 -8, 7 5 3 -9, -7 8"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -8, -6 2 9, 7 5 3 -9, -7 8"), 
    ColourLines("1 2 3 5 8, -1 -6, 4 -3 9, -4 -5 -7, 6 -2 -9, 7 -8"), 
    ColourLines("1 6, -1 -2 -9, 4 -3 9, -4 -5 -7, -6 2 3 5 8, 7 -8"), 
    ColourLines("1 2 3 5 7, -1 -6, 4 -3 9, -4 -5 -8, 6 -2 -9, -7 8"), 
    ColourLines("1 6, -1 -2 -9, 4 -3 9, -4 -5 -8, -6 2 3 5 7, -7 8"), 
    ColourLines("1 2 9, -1 -6, 4 5 8, -4 3 -9, 6 -2 -3 -5 -7, 7 -8"), 
    ColourLines("1 6, -1 -2 -3 -5 -7, 4 5 8, -4 3 -9, -6 2 9, 7 -8"), 
    ColourLines("1 2 9, -1 -6, 4 5 7, -4 3 -9, 6 -2 -3 -5 -8, -7 8"), 
    ColourLines("1 6, -1 -2 -3 -5 -8, 4 5 7, -4 3 -9, -6 2 9, -7 8"), 
    ColourLines("1 2 3 -4, -1 -6, 4 5 8, 6 -2 -9, 7 -8, -7 -5 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 4 5 8, -4 3 2 -6, 7 -8, -7 -5 -3 9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 5 7, 6 -2 -9, -7 8, -8 -5 -3 9"), 
    ColourLines("1 6, -1 -2 -9, 4 5 7, -4 3 2 -6, -7 8, -8 -5 -3 9")
  }; 
  static const ColourLines diag36[16] = { 
    ColourLines("1 2 -3, -1 -6, 3 4 5 9, 6 -2 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 6, -1 -2 -4 -8, 3 4 5 9, -3 2 -6, 7 -9, -7 -5 8"), 
    ColourLines("1 2 -3, -1 -6, 3 4 5 7, 6 -2 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 6, -1 -2 -4 -8, 3 4 5 7, -3 2 -6, -7 9, 8 -5 -9"), 
    ColourLines("1 2 4 5 9, -1 -6, 3 -2 6, -3 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -8, -6 2 4 5 9, 7 -9, -7 -5 8"), 
    ColourLines("1 2 4 5 7, -1 -6, 3 -2 6, -3 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -8, -6 2 4 5 7, -7 9, 8 -5 -9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 8, 6 -2 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 6, -1 -2 -4 -5 -7, 3 4 8, -3 2 -6, 7 -9, -8 5 9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 8, 6 -2 -4 -5 -9, 7 5 -8, -7 9"), 
    ColourLines("1 6, -1 -2 -4 -5 -9, 3 4 8, -3 2 -6, 7 5 -8, -7 9"), 
    ColourLines("1 2 4 8, -1 -6, 3 -2 6, -3 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -7, -6 2 4 8, 7 -9, -8 5 9"), 
    ColourLines("1 2 4 8, -1 -6, 3 -2 6, -3 -4 -5 -9, 7 5 -8, -7 9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -9, -6 2 4 8, 7 5 -8, -7 9")
  }; 
  static const ColourLines diag37[16] = { 
    ColourLines("1 2 8, -1 -6, 4 -3 -2 6, -4 -5 -7, 7 -9, -8 3 5 9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -7, -6 2 8, 7 -9, -8 3 5 9"), 
    ColourLines("1 2 8, -1 -6, 4 -3 -2 6, -4 -5 -9, 7 5 3 -8, -7 9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -9, -6 2 8, 7 5 3 -8, -7 9"), 
    ColourLines("1 2 3 5 9, -1 -6, 4 -3 8, -4 -5 -7, 6 -2 -8, 7 -9"), 
    ColourLines("1 6, -1 -2 -8, 4 -3 8, -4 -5 -7, -6 2 3 5 9, 7 -9"), 
    ColourLines("1 2 3 5 7, -1 -6, 4 -3 8, -4 -5 -9, 6 -2 -8, -7 9"), 
    ColourLines("1 6, -1 -2 -8, 4 -3 8, -4 -5 -9, -6 2 3 5 7, -7 9"), 
    ColourLines("1 2 8, -1 -6, 4 5 9, -4 3 -8, 6 -2 -3 -5 -7, 7 -9"), 
    ColourLines("1 6, -1 -2 -3 -5 -7, 4 5 9, -4 3 -8, -6 2 8, 7 -9"), 
    ColourLines("1 2 8, -1 -6, 4 5 7, -4 3 -8, 6 -2 -3 -5 -9, -7 9"), 
    ColourLines("1 6, -1 -2 -3 -5 -9, 4 5 7, -4 3 -8, -6 2 8, -7 9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 5 9, 6 -2 -8, 7 -9, -7 -5 -3 8"), 
    ColourLines("1 6, -1 -2 -8, 4 5 9, -4 3 2 -6, 7 -9, -7 -5 -3 8"), 
    ColourLines("1 2 3 -4, -1 -6, 4 5 7, 6 -2 -8, -7 9, 8 -3 -5 -9"), 
    ColourLines("1 6, -1 -2 -8, 4 5 7, -4 3 2 -6, -7 9, 8 -3 -5 -9")
  }; 
  static const ColourLines diag38[16] = { 
    ColourLines("1 2 -3, -1 -6, 3 4 5 9, 6 -2 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 6, -1 -2 -4 -7, 3 4 5 9, -3 2 -6, 7 -5 -8, 8 -9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 5 8, 6 -2 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 6, -1 -2 -4 -7, 3 4 5 8, -3 2 -6, 7 -5 -9, -8 9"), 
    ColourLines("1 2 4 5 9, -1 -6, 3 -2 6, -3 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -7, -6 2 4 5 9, 7 -5 -8, 8 -9"), 
    ColourLines("1 2 4 5 8, -1 -6, 3 -2 6, -3 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -7, -6 2 4 5 8, 7 -5 -9, -8 9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 7, 6 -2 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 6, -1 -2 -4 -5 -8, 3 4 7, -3 2 -6, -7 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -6, 3 4 7, 6 -2 -4 -5 -9, -7 5 8, -8 9"), 
    ColourLines("1 6, -1 -2 -4 -5 -9, 3 4 7, -3 2 -6, -7 5 8, -8 9"), 
    ColourLines("1 2 4 7, -1 -6, 3 -2 6, -3 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -8, -6 2 4 7, -7 5 9, 8 -9"), 
    ColourLines("1 2 4 7, -1 -6, 3 -2 6, -3 -4 -5 -9, -7 5 8, -8 9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -9, -6 2 4 7, -7 5 8, -8 9")
  }; 
  static const ColourLines diag39[16] = { 
    ColourLines("1 2 7, -1 -6, 4 -3 -2 6, -4 -5 -8, -7 3 5 9, 8 -9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -8, -6 2 7, -7 3 5 9, 8 -9"), 
    ColourLines("1 2 7, -1 -6, 4 -3 -2 6, -4 -5 -9, -7 3 5 8, -8 9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -9, -6 2 7, -7 3 5 8, -8 9"), 
    ColourLines("1 2 3 5 9, -1 -6, 4 -3 7, -4 -5 -8, 6 -2 -7, 8 -9"), 
    ColourLines("1 6, -1 -2 -7, 4 -3 7, -4 -5 -8, -6 2 3 5 9, 8 -9"), 
    ColourLines("1 2 3 5 8, -1 -6, 4 -3 7, -4 -5 -9, 6 -2 -7, -8 9"), 
    ColourLines("1 6, -1 -2 -7, 4 -3 7, -4 -5 -9, -6 2 3 5 8, -8 9"), 
    ColourLines("1 2 7, -1 -6, 4 5 9, -4 3 -7, 6 -2 -3 -5 -8, 8 -9"), 
    ColourLines("1 6, -1 -2 -3 -5 -8, 4 5 9, -4 3 -7, -6 2 7, 8 -9"), 
    ColourLines("1 2 7, -1 -6, 4 5 8, -4 3 -7, 6 -2 -3 -5 -9, -8 9"), 
    ColourLines("1 6, -1 -2 -3 -5 -9, 4 5 8, -4 3 -7, -6 2 7, -8 9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 5 9, 6 -2 -7, 7 -3 -5 -8, 8 -9"), 
    ColourLines("1 6, -1 -2 -7, 4 5 9, -4 3 2 -6, 7 -3 -5 -8, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -6, 4 5 8, 6 -2 -7, 7 -3 -5 -9, -8 9"), 
    ColourLines("1 6, -1 -2 -7, 4 5 8, -4 3 2 -6, 7 -3 -5 -9, -8 9")
  }; 
  static const ColourLines diag40[16] = { 
    ColourLines("1 2 8, -1 -7, 5 -4 9, -5 -6, 6 4 3 -8, 7 -2 -3 -9"), 
    ColourLines("1 7, -1 -2 -3 -9, 5 -4 9, -5 -6, 6 4 3 -8, -7 2 8"), 
    ColourLines("1 2 8, -1 -7, 5 6, -5 4 3 -8, -6 -4 9, 7 -2 -3 -9"), 
    ColourLines("1 7, -1 -2 -3 -9, 5 6, -5 4 3 -8, -6 -4 9, -7 2 8"), 
    ColourLines("1 2 3 4 6, -1 -7, 5 -4 9, -5 -6, 7 -2 -8, 8 -3 -9"), 
    ColourLines("1 7, -1 -2 -8, 5 -4 9, -5 -6, 6 4 3 2 -7, 8 -3 -9"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 6, -6 -4 9, 7 -2 -8, 8 -3 -9"), 
    ColourLines("1 7, -1 -2 -8, 5 6, -5 4 3 2 -7, -6 -4 9, 8 -3 -9"), 
    ColourLines("1 2 8, -1 -7, 5 -4 -3 -2 7, -5 -6, 6 4 -9, -8 3 9"), 
    ColourLines("1 7, -1 -2 -3 -4 5, -5 -6, 6 4 -9, -7 2 8, -8 3 9"), 
    ColourLines("1 2 8, -1 -7, 5 6, -5 4 -9, -6 -4 -3 -2 7, -8 3 9"), 
    ColourLines("1 7, -1 -2 -3 -4 -6, 5 6, -5 4 -9, -7 2 8, -8 3 9"), 
    ColourLines("1 2 3 9, -1 -7, 5 -4 -3 8, -5 -6, 6 4 -9, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 5 -4 -3 8, -5 -6, 6 4 -9, -7 2 3 9"), 
    ColourLines("1 2 3 9, -1 -7, 5 6, -5 4 -9, -6 -4 -3 8, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 5 6, -5 4 -9, -6 -4 -3 8, -7 2 3 9")
  }; 
  static const ColourLines diag41[16] = { 
    ColourLines("1 2 9, -1 -7, 5 -4 8, -5 -6, 6 4 3 -9, 7 -2 -3 -8"), 
    ColourLines("1 7, -1 -2 -3 -8, 5 -4 8, -5 -6, 6 4 3 -9, -7 2 9"), 
    ColourLines("1 2 9, -1 -7, 5 6, -5 4 3 -9, -6 -4 8, 7 -2 -3 -8"), 
    ColourLines("1 7, -1 -2 -3 -8, 5 6, -5 4 3 -9, -6 -4 8, -7 2 9"), 
    ColourLines("1 2 3 4 6, -1 -7, 5 -4 8, -5 -6, 7 -2 -9, -8 -3 9"), 
    ColourLines("1 7, -1 -2 -9, 5 -4 8, -5 -6, 6 4 3 2 -7, -8 -3 9"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 6, -6 -4 8, 7 -2 -9, -8 -3 9"), 
    ColourLines("1 7, -1 -2 -9, 5 6, -5 4 3 2 -7, -6 -4 8, -8 -3 9"), 
    ColourLines("1 2 9, -1 -7, 5 -4 -3 -2 7, -5 -6, 6 4 -8, 8 3 -9"), 
    ColourLines("1 7, -1 -2 -3 -4 5, -5 -6, 6 4 -8, -7 2 9, 8 3 -9"), 
    ColourLines("1 2 9, -1 -7, 5 6, -5 4 -8, -6 -4 -3 -2 7, 8 3 -9"), 
    ColourLines("1 7, -1 -2 -3 -4 -6, 5 6, -5 4 -8, -7 2 9, 8 3 -9"), 
    ColourLines("1 2 3 8, -1 -7, 5 -4 -3 9, -5 -6, 6 4 -8, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 5 -4 -3 9, -5 -6, 6 4 -8, -7 2 3 8"), 
    ColourLines("1 2 3 8, -1 -7, 5 6, -5 4 -8, -6 -4 -3 9, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 5 6, -5 4 -8, -6 -4 -3 9, -7 2 3 8")
  }; 
  static const ColourLines diag42[16] = { 
    ColourLines("1 2 6, -1 -7, 5 -4 9, -5 -8, -6 3 4 8, 7 -2 -3 -9"), 
    ColourLines("1 7, -1 -2 -3 -9, 5 -4 9, -5 -8, 6 2 -7, -6 3 4 8"), 
    ColourLines("1 2 6, -1 -7, 5 8, -5 4 3 -6, 7 -2 -3 -9, -8 -4 9"), 
    ColourLines("1 7, -1 -2 -3 -9, 5 8, -5 4 3 -6, 6 2 -7, -8 -4 9"), 
    ColourLines("1 2 3 4 8, -1 -7, 5 -4 9, -5 -8, 6 -3 -9, -6 -2 7"), 
    ColourLines("1 7, -1 -2 -6, 5 -4 9, -5 -8, 6 -3 -9, -7 2 3 4 8"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 8, 6 -3 -9, -6 -2 7, -8 -4 9"), 
    ColourLines("1 7, -1 -2 -6, 5 8, -5 4 3 2 -7, 6 -3 -9, -8 -4 9"), 
    ColourLines("1 2 6, -1 -7, 5 -4 -3 -2 7, -5 -8, -6 3 9, 8 4 -9"), 
    ColourLines("1 7, -1 -2 -3 -4 5, -5 -8, 6 2 -7, -6 3 9, 8 4 -9"), 
    ColourLines("1 2 6, -1 -7, 5 8, -5 4 -9, -6 3 9, 7 -2 -3 -4 -8"), 
    ColourLines("1 7, -1 -2 -3 -4 -8, 5 8, -5 4 -9, 6 2 -7, -6 3 9"), 
    ColourLines("1 2 3 9, -1 -7, 5 -4 -3 6, -5 -8, -6 -2 7, 8 4 -9"), 
    ColourLines("1 7, -1 -2 -6, 5 -4 -3 6, -5 -8, -7 2 3 9, 8 4 -9"), 
    ColourLines("1 2 3 9, -1 -7, 5 8, -5 4 -9, 6 -3 -4 -8, -6 -2 7"), 
    ColourLines("1 7, -1 -2 -6, 5 8, -5 4 -9, 6 -3 -4 -8, -7 2 3 9")
  }; 
  static const ColourLines diag43[16] = { 
    ColourLines("1 2 9, -1 -7, 5 -4 6, -5 -8, -6 -3 -2 7, 8 4 3 -9"), 
    ColourLines("1 7, -1 -2 -3 -6, 5 -4 6, -5 -8, -7 2 9, 8 4 3 -9"), 
    ColourLines("1 2 9, -1 -7, 5 8, -5 4 3 -9, 6 -4 -8, -6 -3 -2 7"), 
    ColourLines("1 7, -1 -2 -3 -6, 5 8, -5 4 3 -9, 6 -4 -8, -7 2 9"), 
    ColourLines("1 2 3 4 8, -1 -7, 5 -4 6, -5 -8, -6 -3 9, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 5 -4 6, -5 -8, -6 -3 9, -7 2 3 4 8"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 8, 6 -4 -8, -6 -3 9, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 5 8, -5 4 3 2 -7, 6 -4 -8, -6 -3 9"), 
    ColourLines("1 2 9, -1 -7, 5 -4 -3 -2 7, -5 -8, 6 3 -9, -6 4 8"), 
    ColourLines("1 7, -1 -2 -3 -4 5, -5 -8, 6 3 -9, -6 4 8, -7 2 9"), 
    ColourLines("1 2 9, -1 -7, 5 8, -5 4 -6, 6 3 -9, 7 -2 -3 -4 -8"), 
    ColourLines("1 7, -1 -2 -3 -4 -8, 5 8, -5 4 -6, 6 3 -9, -7 2 9"), 
    ColourLines("1 2 3 6, -1 -7, 5 -4 -3 9, -5 -8, -6 4 8, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 5 -4 -3 9, -5 -8, 6 3 2 -7, -6 4 8"), 
    ColourLines("1 2 3 6, -1 -7, 5 8, -5 4 -6, 7 -2 -9, -8 -4 -3 9"), 
    ColourLines("1 7, -1 -2 -9, 5 8, -5 4 -6, 6 3 2 -7, -8 -4 -3 9")
  }; 
  static const ColourLines diag44[16] = { 
    ColourLines("1 2 6, -1 -7, 5 -4 8, -5 -9, -6 3 4 9, 7 -2 -3 -8"), 
    ColourLines("1 7, -1 -2 -3 -8, 5 -4 8, -5 -9, 6 2 -7, -6 3 4 9"), 
    ColourLines("1 2 6, -1 -7, 5 9, -5 4 3 -6, 7 -2 -3 -8, 8 -4 -9"), 
    ColourLines("1 7, -1 -2 -3 -8, 5 9, -5 4 3 -6, 6 2 -7, 8 -4 -9"), 
    ColourLines("1 2 3 4 9, -1 -7, 5 -4 8, -5 -9, 6 -3 -8, -6 -2 7"), 
    ColourLines("1 7, -1 -2 -6, 5 -4 8, -5 -9, 6 -3 -8, -7 2 3 4 9"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 9, 6 -3 -8, -6 -2 7, 8 -4 -9"), 
    ColourLines("1 7, -1 -2 -6, 5 9, -5 4 3 2 -7, 6 -3 -8, 8 -4 -9"), 
    ColourLines("1 2 6, -1 -7, 5 -4 -3 -2 7, -5 -9, -6 3 8, -8 4 9"), 
    ColourLines("1 7, -1 -2 -3 -4 5, -5 -9, 6 2 -7, -6 3 8, -8 4 9"), 
    ColourLines("1 2 6, -1 -7, 5 9, -5 4 -8, -6 3 8, 7 -2 -3 -4 -9"), 
    ColourLines("1 7, -1 -2 -3 -4 -9, 5 9, -5 4 -8, 6 2 -7, -6 3 8"), 
    ColourLines("1 2 3 8, -1 -7, 5 -4 -3 6, -5 -9, -6 -2 7, -8 4 9"), 
    ColourLines("1 7, -1 -2 -6, 5 -4 -3 6, -5 -9, -7 2 3 8, -8 4 9"), 
    ColourLines("1 2 3 8, -1 -7, 5 9, -5 4 -8, 6 -3 -4 -9, -6 -2 7"), 
    ColourLines("1 7, -1 -2 -6, 5 9, -5 4 -8, 6 -3 -4 -9, -7 2 3 8")
  }; 
  static const ColourLines diag45[16] = { 
    ColourLines("1 2 8, -1 -7, 5 -4 6, -5 -9, -6 -3 -2 7, -8 3 4 9"), 
    ColourLines("1 7, -1 -2 -3 -6, 5 -4 6, -5 -9, -7 2 8, -8 3 4 9"), 
    ColourLines("1 2 8, -1 -7, 5 9, -5 4 3 -8, 6 -4 -9, -6 -3 -2 7"), 
    ColourLines("1 7, -1 -2 -3 -6, 5 9, -5 4 3 -8, 6 -4 -9, -7 2 8"), 
    ColourLines("1 2 3 4 9, -1 -7, 5 -4 6, -5 -9, -6 -3 8, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 5 -4 6, -5 -9, -6 -3 8, -7 2 3 4 9"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 9, 6 -4 -9, -6 -3 8, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 5 9, -5 4 3 2 -7, 6 -4 -9, -6 -3 8"), 
    ColourLines("1 2 8, -1 -7, 5 -4 -3 -2 7, -5 -9, 6 3 -8, -6 4 9"), 
    ColourLines("1 7, -1 -2 -3 -4 5, -5 -9, 6 3 -8, -6 4 9, -7 2 8"), 
    ColourLines("1 2 8, -1 -7, 5 9, -5 4 -6, 6 3 -8, 7 -2 -3 -4 -9"), 
    ColourLines("1 7, -1 -2 -3 -4 -9, 5 9, -5 4 -6, 6 3 -8, -7 2 8"), 
    ColourLines("1 2 3 6, -1 -7, 5 -4 -3 8, -5 -9, -6 4 9, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 5 -4 -3 8, -5 -9, 6 3 2 -7, -6 4 9"), 
    ColourLines("1 2 3 6, -1 -7, 5 9, -5 4 -6, 7 -2 -8, 8 -3 -4 -9"), 
    ColourLines("1 7, -1 -2 -8, 5 9, -5 4 -6, 6 3 2 -7, 8 -3 -4 -9")
  }; 
  static const ColourLines diag46[16] = { 
    ColourLines("1 2 -3, -1 -7, 3 4 5 8, 6 -8, -6 -5 9, 7 -2 -4 -9"), 
    ColourLines("1 7, -1 -2 -4 -9, 3 4 5 8, -3 2 -7, 6 -8, -6 -5 9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 5 6, -6 8, 7 -2 -4 -9, -8 -5 9"), 
    ColourLines("1 7, -1 -2 -4 -9, 3 4 5 6, -3 2 -7, -6 8, -8 -5 9"), 
    ColourLines("1 2 4 5 8, -1 -7, 3 -2 7, -3 -4 -9, 6 -8, -6 -5 9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -9, 6 -8, -6 -5 9, -7 2 4 5 8"), 
    ColourLines("1 2 4 5 6, -1 -7, 3 -2 7, -3 -4 -9, -6 8, -8 -5 9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -9, 6 5 4 2 -7, -6 8, -8 -5 9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 9, 6 -8, -6 -5 -4 -2 7, 8 5 -9"), 
    ColourLines("1 7, -1 -2 -4 -5 -6, 3 4 9, -3 2 -7, 6 -8, 8 5 -9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 9, 6 5 -9, -6 8, 7 -2 -4 -5 -8"), 
    ColourLines("1 7, -1 -2 -4 -5 -8, 3 4 9, -3 2 -7, 6 5 -9, -6 8"), 
    ColourLines("1 2 4 9, -1 -7, 3 -2 7, -3 -4 -5 -6, 6 -8, 8 5 -9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -5 -6, 6 -8, -7 2 4 9, 8 5 -9"), 
    ColourLines("1 2 4 9, -1 -7, 3 -2 7, -3 -4 -5 -8, 6 5 -9, -6 8"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -5 -8, 6 5 -9, -6 8, -7 2 4 9")
  }; 
  static const ColourLines diag47[16] = { 
    ColourLines("1 2 9, -1 -7, 4 -3 -2 7, -4 -5 -6, 6 -8, 8 5 3 -9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -5 -6, 6 -8, -7 2 9, 8 5 3 -9"), 
    ColourLines("1 2 9, -1 -7, 4 -3 -2 7, -4 -5 -8, 6 5 3 -9, -6 8"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -5 -8, 6 5 3 -9, -6 8, -7 2 9"), 
    ColourLines("1 2 3 5 8, -1 -7, 4 -3 9, -4 -5 -6, 6 -8, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 4 -3 9, -4 -5 -6, 6 -8, -7 2 3 5 8"), 
    ColourLines("1 2 3 5 6, -1 -7, 4 -3 9, -4 -5 -8, -6 8, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 4 -3 9, -4 -5 -8, 6 5 3 2 -7, -6 8"), 
    ColourLines("1 2 9, -1 -7, 4 5 8, -4 3 -9, 6 -8, -6 -5 -3 -2 7"), 
    ColourLines("1 7, -1 -2 -3 -5 -6, 4 5 8, -4 3 -9, 6 -8, -7 2 9"), 
    ColourLines("1 2 9, -1 -7, 4 5 6, -4 3 -9, -6 8, 7 -2 -3 -5 -8"), 
    ColourLines("1 7, -1 -2 -3 -5 -8, 4 5 6, -4 3 -9, -6 8, -7 2 9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 8, 6 -8, -6 -5 -3 9, 7 -2 -9"), 
    ColourLines("1 7, -1 -2 -9, 4 5 8, -4 3 2 -7, 6 -8, -6 -5 -3 9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 6, -6 8, 7 -2 -9, -8 -5 -3 9"), 
    ColourLines("1 7, -1 -2 -9, 4 5 6, -4 3 2 -7, -6 8, -8 -5 -3 9")
  }; 
  static const ColourLines diag48[16] = { 
    ColourLines("1 2 -3, -1 -7, 3 4 5 9, 6 -9, -6 -5 8, 7 -2 -4 -8"), 
    ColourLines("1 7, -1 -2 -4 -8, 3 4 5 9, -3 2 -7, 6 -9, -6 -5 8"), 
    ColourLines("1 2 -3, -1 -7, 3 4 5 6, -6 9, 7 -2 -4 -8, 8 -5 -9"), 
    ColourLines("1 7, -1 -2 -4 -8, 3 4 5 6, -3 2 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 2 4 5 9, -1 -7, 3 -2 7, -3 -4 -8, 6 -9, -6 -5 8"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -8, 6 -9, -6 -5 8, -7 2 4 5 9"), 
    ColourLines("1 2 4 5 6, -1 -7, 3 -2 7, -3 -4 -8, -6 9, 8 -5 -9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -8, 6 5 4 2 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 8, 6 -9, -6 -5 -4 -2 7, -8 5 9"), 
    ColourLines("1 7, -1 -2 -4 -5 -6, 3 4 8, -3 2 -7, 6 -9, -8 5 9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 8, 6 5 -8, -6 9, 7 -2 -4 -5 -9"), 
    ColourLines("1 7, -1 -2 -4 -5 -9, 3 4 8, -3 2 -7, 6 5 -8, -6 9"), 
    ColourLines("1 2 4 8, -1 -7, 3 -2 7, -3 -4 -5 -6, 6 -9, -8 5 9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -5 -6, 6 -9, -7 2 4 8, -8 5 9"), 
    ColourLines("1 2 4 8, -1 -7, 3 -2 7, -3 -4 -5 -9, 6 5 -8, -6 9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -5 -9, 6 5 -8, -6 9, -7 2 4 8")
  }; 
  static const ColourLines diag49[16] = { 
    ColourLines("1 2 8, -1 -7, 4 -3 -2 7, -4 -5 -6, 6 -9, -8 3 5 9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -5 -6, 6 -9, -7 2 8, -8 3 5 9"), 
    ColourLines("1 2 8, -1 -7, 4 -3 -2 7, -4 -5 -9, 6 5 3 -8, -6 9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -5 -9, 6 5 3 -8, -6 9, -7 2 8"), 
    ColourLines("1 2 3 5 9, -1 -7, 4 -3 8, -4 -5 -6, 6 -9, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 4 -3 8, -4 -5 -6, 6 -9, -7 2 3 5 9"), 
    ColourLines("1 2 3 5 6, -1 -7, 4 -3 8, -4 -5 -9, -6 9, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 4 -3 8, -4 -5 -9, 6 5 3 2 -7, -6 9"), 
    ColourLines("1 2 8, -1 -7, 4 5 9, -4 3 -8, 6 -9, -6 -5 -3 -2 7"), 
    ColourLines("1 7, -1 -2 -3 -5 -6, 4 5 9, -4 3 -8, 6 -9, -7 2 8"), 
    ColourLines("1 2 8, -1 -7, 4 5 6, -4 3 -8, -6 9, 7 -2 -3 -5 -9"), 
    ColourLines("1 7, -1 -2 -3 -5 -9, 4 5 6, -4 3 -8, -6 9, -7 2 8"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 9, 6 -9, -6 -5 -3 8, 7 -2 -8"), 
    ColourLines("1 7, -1 -2 -8, 4 5 9, -4 3 2 -7, 6 -9, -6 -5 -3 8"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 6, -6 9, 7 -2 -8, 8 -3 -5 -9"), 
    ColourLines("1 7, -1 -2 -8, 4 5 6, -4 3 2 -7, -6 9, 8 -3 -5 -9")
  }; 
  static const ColourLines diag50[16] = { 
    ColourLines("1 2 -3, -1 -7, 3 4 5 9, 6 -5 -8, -6 -4 -2 7, 8 -9"), 
    ColourLines("1 7, -1 -2 -4 -6, 3 4 5 9, -3 2 -7, 6 -5 -8, 8 -9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 5 8, 6 -5 -9, -6 -4 -2 7, -8 9"), 
    ColourLines("1 7, -1 -2 -4 -6, 3 4 5 8, -3 2 -7, 6 -5 -9, -8 9"), 
    ColourLines("1 2 4 5 9, -1 -7, 3 -2 7, -3 -4 -6, 6 -5 -8, 8 -9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -6, 6 -5 -8, -7 2 4 5 9, 8 -9"), 
    ColourLines("1 2 4 5 8, -1 -7, 3 -2 7, -3 -4 -6, 6 -5 -9, -8 9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -6, 6 -5 -9, -7 2 4 5 8, -8 9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 6, -6 5 9, 7 -2 -4 -5 -8, 8 -9"), 
    ColourLines("1 7, -1 -2 -4 -5 -8, 3 4 6, -3 2 -7, -6 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 6, -6 5 8, 7 -2 -4 -5 -9, -8 9"), 
    ColourLines("1 7, -1 -2 -4 -5 -9, 3 4 6, -3 2 -7, -6 5 8, -8 9"), 
    ColourLines("1 2 4 6, -1 -7, 3 -2 7, -3 -4 -5 -8, -6 5 9, 8 -9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -5 -8, 6 4 2 -7, -6 5 9, 8 -9"), 
    ColourLines("1 2 4 6, -1 -7, 3 -2 7, -3 -4 -5 -9, -6 5 8, -8 9"), 
    ColourLines("1 7, -1 -2 3, -3 -4 -5 -9, 6 4 2 -7, -6 5 8, -8 9")
  }; 
  static const ColourLines diag51[16] = { 
    ColourLines("1 2 6, -1 -7, 4 -3 -2 7, -4 -5 -8, -6 3 5 9, 8 -9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -5 -8, 6 2 -7, -6 3 5 9, 8 -9"), 
    ColourLines("1 2 6, -1 -7, 4 -3 -2 7, -4 -5 -9, -6 3 5 8, -8 9"), 
    ColourLines("1 7, -1 -2 -3 4, -4 -5 -9, 6 2 -7, -6 3 5 8, -8 9"), 
    ColourLines("1 2 3 5 9, -1 -7, 4 -3 6, -4 -5 -8, -6 -2 7, 8 -9"), 
    ColourLines("1 7, -1 -2 -6, 4 -3 6, -4 -5 -8, -7 2 3 5 9, 8 -9"), 
    ColourLines("1 2 3 5 8, -1 -7, 4 -3 6, -4 -5 -9, -6 -2 7, -8 9"), 
    ColourLines("1 7, -1 -2 -6, 4 -3 6, -4 -5 -9, -7 2 3 5 8, -8 9"), 
    ColourLines("1 2 6, -1 -7, 4 5 9, -4 3 -6, 7 -2 -3 -5 -8, 8 -9"), 
    ColourLines("1 7, -1 -2 -3 -5 -8, 4 5 9, -4 3 -6, 6 2 -7, 8 -9"), 
    ColourLines("1 2 6, -1 -7, 4 5 8, -4 3 -6, 7 -2 -3 -5 -9, -8 9"), 
    ColourLines("1 7, -1 -2 -3 -5 -9, 4 5 8, -4 3 -6, 6 2 -7, -8 9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 9, 6 -3 -5 -8, -6 -2 7, 8 -9"), 
    ColourLines("1 7, -1 -2 -6, 4 5 9, -4 3 2 -7, 6 -3 -5 -8, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 8, 6 -3 -5 -9, -6 -2 7, -8 9"), 
    ColourLines("1 7, -1 -2 -6, 4 5 8, -4 3 2 -7, 6 -3 -5 -9, -8 9")
  }; 
  static const ColourLines diag52[16] = { 
    ColourLines("1 2 7, -1 -8, 5 -4 9, -5 -6, 6 4 3 -7, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 5 -4 9, -5 -6, 6 4 3 -7, 7 2 -8"), 
    ColourLines("1 2 7, -1 -8, 5 6, -5 4 3 -7, -6 -4 9, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 5 6, -5 4 3 -7, -6 -4 9, 7 2 -8"), 
    ColourLines("1 2 3 4 6, -1 -8, 5 -4 9, -5 -6, 7 -3 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 -4 9, -5 -6, 6 4 3 2 -8, 7 -3 -9"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 6, -6 -4 9, 7 -3 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 6, -5 4 3 2 -8, -6 -4 9, 7 -3 -9"), 
    ColourLines("1 2 7, -1 -8, 5 -4 -3 -2 8, -5 -6, 6 4 -9, -7 3 9"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -6, 6 4 -9, 7 2 -8, -7 3 9"), 
    ColourLines("1 2 7, -1 -8, 5 6, -5 4 -9, -6 -4 -3 -2 8, -7 3 9"), 
    ColourLines("1 8, -1 -2 -3 -4 -6, 5 6, -5 4 -9, 7 2 -8, -7 3 9"), 
    ColourLines("1 2 3 9, -1 -8, 5 -4 -3 7, -5 -6, 6 4 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 -4 -3 7, -5 -6, 6 4 -9, -8 2 3 9"), 
    ColourLines("1 2 3 9, -1 -8, 5 6, -5 4 -9, -6 -4 -3 7, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 6, -5 4 -9, -6 -4 -3 7, -8 2 3 9")
  }; 
  static const ColourLines diag53[16] = { 
    ColourLines("1 2 9, -1 -8, 5 -4 7, -5 -6, 6 4 3 -9, -7 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -7, 5 -4 7, -5 -6, 6 4 3 -9, -8 2 9"), 
    ColourLines("1 2 9, -1 -8, 5 6, -5 4 3 -9, -6 -4 7, -7 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -7, 5 6, -5 4 3 -9, -6 -4 7, -8 2 9"), 
    ColourLines("1 2 3 4 6, -1 -8, 5 -4 7, -5 -6, -7 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 -4 7, -5 -6, 6 4 3 2 -8, -7 -3 9"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 6, -6 -4 7, -7 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 6, -5 4 3 2 -8, -6 -4 7, -7 -3 9"), 
    ColourLines("1 2 9, -1 -8, 5 -4 -3 -2 8, -5 -6, 6 4 -7, 7 3 -9"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -6, 6 4 -7, 7 3 -9, -8 2 9"), 
    ColourLines("1 2 9, -1 -8, 5 6, -5 4 -7, -6 -4 -3 -2 8, 7 3 -9"), 
    ColourLines("1 8, -1 -2 -3 -4 -6, 5 6, -5 4 -7, 7 3 -9, -8 2 9"), 
    ColourLines("1 2 3 7, -1 -8, 5 -4 -3 9, -5 -6, 6 4 -7, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 -4 -3 9, -5 -6, 6 4 -7, 7 3 2 -8"), 
    ColourLines("1 2 3 7, -1 -8, 5 6, -5 4 -7, -6 -4 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 6, -5 4 -7, -6 -4 -3 9, 7 3 2 -8")
  }; 
  static const ColourLines diag54[16] = { 
    ColourLines("1 2 6, -1 -8, 5 -4 9, -5 -7, -6 3 4 7, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 5 -4 9, -5 -7, 6 2 -8, -6 3 4 7"), 
    ColourLines("1 2 6, -1 -8, 5 7, -5 4 3 -6, -7 -4 9, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 5 7, -5 4 3 -6, 6 2 -8, -7 -4 9"), 
    ColourLines("1 2 3 4 7, -1 -8, 5 -4 9, -5 -7, 6 -3 -9, -6 -2 8"), 
    ColourLines("1 8, -1 -2 -6, 5 -4 9, -5 -7, 6 -3 -9, 7 4 3 2 -8"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 7, 6 -3 -9, -6 -2 8, -7 -4 9"), 
    ColourLines("1 8, -1 -2 -6, 5 7, -5 4 3 2 -8, 6 -3 -9, -7 -4 9"), 
    ColourLines("1 2 6, -1 -8, 5 -4 -3 -2 8, -5 -7, -6 3 9, 7 4 -9"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -7, 6 2 -8, -6 3 9, 7 4 -9"), 
    ColourLines("1 2 6, -1 -8, 5 7, -5 4 -9, -6 3 9, -7 -4 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -4 -7, 5 7, -5 4 -9, 6 2 -8, -6 3 9"), 
    ColourLines("1 2 3 9, -1 -8, 5 -4 -3 6, -5 -7, -6 -2 8, 7 4 -9"), 
    ColourLines("1 8, -1 -2 -6, 5 -4 -3 6, -5 -7, 7 4 -9, -8 2 3 9"), 
    ColourLines("1 2 3 9, -1 -8, 5 7, -5 4 -9, 6 -3 -4 -7, -6 -2 8"), 
    ColourLines("1 8, -1 -2 -6, 5 7, -5 4 -9, 6 -3 -4 -7, -8 2 3 9")
  }; 
  static const ColourLines diag55[16] = { 
    ColourLines("1 2 9, -1 -8, 5 -4 6, -5 -7, -6 -3 -2 8, 7 4 3 -9"), 
    ColourLines("1 8, -1 -2 -3 -6, 5 -4 6, -5 -7, 7 4 3 -9, -8 2 9"), 
    ColourLines("1 2 9, -1 -8, 5 7, -5 4 3 -9, 6 -4 -7, -6 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -6, 5 7, -5 4 3 -9, 6 -4 -7, -8 2 9"), 
    ColourLines("1 2 3 4 7, -1 -8, 5 -4 6, -5 -7, -6 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 -4 6, -5 -7, -6 -3 9, 7 4 3 2 -8"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 7, 6 -4 -7, -6 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 7, -5 4 3 2 -8, 6 -4 -7, -6 -3 9"), 
    ColourLines("1 2 9, -1 -8, 5 -4 -3 -2 8, -5 -7, 6 3 -9, -6 4 7"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -7, 6 3 -9, -6 4 7, -8 2 9"), 
    ColourLines("1 2 9, -1 -8, 5 7, -5 4 -6, 6 3 -9, -7 -4 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -4 -7, 5 7, -5 4 -6, 6 3 -9, -8 2 9"), 
    ColourLines("1 2 3 6, -1 -8, 5 -4 -3 9, -5 -7, -6 4 7, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 -4 -3 9, -5 -7, 6 3 2 -8, -6 4 7"), 
    ColourLines("1 2 3 6, -1 -8, 5 7, -5 4 -6, -7 -4 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 7, -5 4 -6, 6 3 2 -8, -7 -4 -3 9")
  }; 
  static const ColourLines diag56[16] = { 
    ColourLines("1 2 6, -1 -8, 5 -4 7, -5 -9, -6 3 4 9, -7 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -7, 5 -4 7, -5 -9, 6 2 -8, -6 3 4 9"), 
    ColourLines("1 2 6, -1 -8, 5 9, -5 4 3 -6, 7 -4 -9, -7 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -7, 5 9, -5 4 3 -6, 6 2 -8, 7 -4 -9"), 
    ColourLines("1 2 3 4 9, -1 -8, 5 -4 7, -5 -9, 6 -3 -7, -6 -2 8"), 
    ColourLines("1 8, -1 -2 -6, 5 -4 7, -5 -9, 6 -3 -7, -8 2 3 4 9"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 9, 6 -3 -7, -6 -2 8, 7 -4 -9"), 
    ColourLines("1 8, -1 -2 -6, 5 9, -5 4 3 2 -8, 6 -3 -7, 7 -4 -9"), 
    ColourLines("1 2 6, -1 -8, 5 -4 -3 -2 8, -5 -9, -6 3 7, -7 4 9"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -9, 6 2 -8, -6 3 7, -7 4 9"), 
    ColourLines("1 2 6, -1 -8, 5 9, -5 4 -7, -6 3 7, 8 -2 -3 -4 -9"), 
    ColourLines("1 8, -1 -2 -3 -4 -9, 5 9, -5 4 -7, 6 2 -8, -6 3 7"), 
    ColourLines("1 2 3 7, -1 -8, 5 -4 -3 6, -5 -9, -6 -2 8, -7 4 9"), 
    ColourLines("1 8, -1 -2 -6, 5 -4 -3 6, -5 -9, 7 3 2 -8, -7 4 9"), 
    ColourLines("1 2 3 7, -1 -8, 5 9, -5 4 -7, 6 -3 -4 -9, -6 -2 8"), 
    ColourLines("1 8, -1 -2 -6, 5 9, -5 4 -7, 6 -3 -4 -9, 7 3 2 -8")
  }; 
  static const ColourLines diag57[16] = { 
    ColourLines("1 2 7, -1 -8, 5 -4 6, -5 -9, -6 -3 -2 8, -7 3 4 9"), 
    ColourLines("1 8, -1 -2 -3 -6, 5 -4 6, -5 -9, 7 2 -8, -7 3 4 9"), 
    ColourLines("1 2 7, -1 -8, 5 9, -5 4 3 -7, 6 -4 -9, -6 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -6, 5 9, -5 4 3 -7, 6 -4 -9, 7 2 -8"), 
    ColourLines("1 2 3 4 9, -1 -8, 5 -4 6, -5 -9, -6 -3 7, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 -4 6, -5 -9, -6 -3 7, -8 2 3 4 9"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 9, 6 -4 -9, -6 -3 7, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 9, -5 4 3 2 -8, 6 -4 -9, -6 -3 7"), 
    ColourLines("1 2 7, -1 -8, 5 -4 -3 -2 8, -5 -9, 6 3 -7, -6 4 9"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -9, 6 3 -7, -6 4 9, 7 2 -8"), 
    ColourLines("1 2 7, -1 -8, 5 9, -5 4 -6, 6 3 -7, 8 -2 -3 -4 -9"), 
    ColourLines("1 8, -1 -2 -3 -4 -9, 5 9, -5 4 -6, 6 3 -7, 7 2 -8"), 
    ColourLines("1 2 3 6, -1 -8, 5 -4 -3 7, -5 -9, -6 4 9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 -4 -3 7, -5 -9, 6 3 2 -8, -6 4 9"), 
    ColourLines("1 2 3 6, -1 -8, 5 9, -5 4 -6, 7 -3 -4 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 9, -5 4 -6, 6 3 2 -8, 7 -3 -4 -9")
  }; 
  static const ColourLines diag58[16] = { 
    ColourLines("1 2 -3, -1 -8, 3 4 5 7, 6 -7, -6 -5 9, 8 -2 -4 -9"), 
    ColourLines("1 8, -1 -2 -4 -9, 3 4 5 7, -3 2 -8, 6 -7, -6 -5 9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 5 6, -6 7, -7 -5 9, 8 -2 -4 -9"), 
    ColourLines("1 8, -1 -2 -4 -9, 3 4 5 6, -3 2 -8, -6 7, -7 -5 9"), 
    ColourLines("1 2 4 5 7, -1 -8, 3 -2 8, -3 -4 -9, 6 -7, -6 -5 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -9, 6 -7, -6 -5 9, 7 5 4 2 -8"), 
    ColourLines("1 2 4 5 6, -1 -8, 3 -2 8, -3 -4 -9, -6 7, -7 -5 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -9, 6 5 4 2 -8, -6 7, -7 -5 9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 9, 6 -7, -6 -5 -4 -2 8, 7 5 -9"), 
    ColourLines("1 8, -1 -2 -4 -5 -6, 3 4 9, -3 2 -8, 6 -7, 7 5 -9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 9, 6 5 -9, -6 7, -7 -5 -4 -2 8"), 
    ColourLines("1 8, -1 -2 -4 -5 -7, 3 4 9, -3 2 -8, 6 5 -9, -6 7"), 
    ColourLines("1 2 4 9, -1 -8, 3 -2 8, -3 -4 -5 -6, 6 -7, 7 5 -9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -6, 6 -7, 7 5 -9, -8 2 4 9"), 
    ColourLines("1 2 4 9, -1 -8, 3 -2 8, -3 -4 -5 -7, 6 5 -9, -6 7"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -7, 6 5 -9, -6 7, -8 2 4 9")
  }; 
  static const ColourLines diag59[16] = { 
    ColourLines("1 2 9, -1 -8, 4 -3 -2 8, -4 -5 -6, 6 -7, 7 5 3 -9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -6, 6 -7, 7 5 3 -9, -8 2 9"), 
    ColourLines("1 2 9, -1 -8, 4 -3 -2 8, -4 -5 -7, 6 5 3 -9, -6 7"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -7, 6 5 3 -9, -6 7, -8 2 9"), 
    ColourLines("1 2 3 5 7, -1 -8, 4 -3 9, -4 -5 -6, 6 -7, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 4 -3 9, -4 -5 -6, 6 -7, 7 5 3 2 -8"), 
    ColourLines("1 2 3 5 6, -1 -8, 4 -3 9, -4 -5 -7, -6 7, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 4 -3 9, -4 -5 -7, 6 5 3 2 -8, -6 7"), 
    ColourLines("1 2 9, -1 -8, 4 5 7, -4 3 -9, 6 -7, -6 -5 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -5 -6, 4 5 7, -4 3 -9, 6 -7, -8 2 9"), 
    ColourLines("1 2 9, -1 -8, 4 5 6, -4 3 -9, -6 7, -7 -5 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -5 -7, 4 5 6, -4 3 -9, -6 7, -8 2 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 5 7, 6 -7, -6 -5 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 4 5 7, -4 3 2 -8, 6 -7, -6 -5 -3 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 5 6, -6 7, -7 -5 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 4 5 6, -4 3 2 -8, -6 7, -7 -5 -3 9")
  }; 
  static const ColourLines diag60[16] = { 
    ColourLines("1 2 -3, -1 -8, 3 4 5 9, 6 -9, -6 -5 7, -7 -4 -2 8"), 
    ColourLines("1 8, -1 -2 -4 -7, 3 4 5 9, -3 2 -8, 6 -9, -6 -5 7"), 
    ColourLines("1 2 -3, -1 -8, 3 4 5 6, -6 9, 7 -5 -9, -7 -4 -2 8"), 
    ColourLines("1 8, -1 -2 -4 -7, 3 4 5 6, -3 2 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 2 4 5 9, -1 -8, 3 -2 8, -3 -4 -7, 6 -9, -6 -5 7"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -7, 6 -9, -6 -5 7, -8 2 4 5 9"), 
    ColourLines("1 2 4 5 6, -1 -8, 3 -2 8, -3 -4 -7, -6 9, 7 -5 -9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -7, 6 5 4 2 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 7, 6 -9, -6 -5 -4 -2 8, -7 5 9"), 
    ColourLines("1 8, -1 -2 -4 -5 -6, 3 4 7, -3 2 -8, 6 -9, -7 5 9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 7, 6 5 -7, -6 9, 8 -2 -4 -5 -9"), 
    ColourLines("1 8, -1 -2 -4 -5 -9, 3 4 7, -3 2 -8, 6 5 -7, -6 9"), 
    ColourLines("1 2 4 7, -1 -8, 3 -2 8, -3 -4 -5 -6, 6 -9, -7 5 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -6, 6 -9, 7 4 2 -8, -7 5 9"), 
    ColourLines("1 2 4 7, -1 -8, 3 -2 8, -3 -4 -5 -9, 6 5 -7, -6 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -9, 6 5 -7, -6 9, 7 4 2 -8")
  }; 
  static const ColourLines diag61[16] = { 
    ColourLines("1 2 7, -1 -8, 4 -3 -2 8, -4 -5 -6, 6 -9, -7 3 5 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -6, 6 -9, 7 2 -8, -7 3 5 9"), 
    ColourLines("1 2 7, -1 -8, 4 -3 -2 8, -4 -5 -9, 6 5 3 -7, -6 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -9, 6 5 3 -7, -6 9, 7 2 -8"), 
    ColourLines("1 2 3 5 9, -1 -8, 4 -3 7, -4 -5 -6, 6 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 4 -3 7, -4 -5 -6, 6 -9, -8 2 3 5 9"), 
    ColourLines("1 2 3 5 6, -1 -8, 4 -3 7, -4 -5 -9, -6 9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 4 -3 7, -4 -5 -9, 6 5 3 2 -8, -6 9"), 
    ColourLines("1 2 7, -1 -8, 4 5 9, -4 3 -7, 6 -9, -6 -5 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -5 -6, 4 5 9, -4 3 -7, 6 -9, 7 2 -8"), 
    ColourLines("1 2 7, -1 -8, 4 5 6, -4 3 -7, -6 9, 8 -2 -3 -5 -9"), 
    ColourLines("1 8, -1 -2 -3 -5 -9, 4 5 6, -4 3 -7, -6 9, 7 2 -8"), 
    ColourLines("1 2 3 -4, -1 -8, 4 5 9, 6 -9, -6 -5 -3 7, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 4 5 9, -4 3 2 -8, 6 -9, -6 -5 -3 7"), 
    ColourLines("1 2 3 -4, -1 -8, 4 5 6, -6 9, 7 -3 -5 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 4 5 6, -4 3 2 -8, -6 9, 7 -3 -5 -9")
  }; 
  static const ColourLines diag62[16] = { 
    ColourLines("1 2 -3, -1 -8, 3 4 5 9, 6 -5 -7, -6 -4 -2 8, 7 -9"), 
    ColourLines("1 8, -1 -2 -4 -6, 3 4 5 9, -3 2 -8, 6 -5 -7, 7 -9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 5 7, 6 -5 -9, -6 -4 -2 8, -7 9"), 
    ColourLines("1 8, -1 -2 -4 -6, 3 4 5 7, -3 2 -8, 6 -5 -9, -7 9"), 
    ColourLines("1 2 4 5 9, -1 -8, 3 -2 8, -3 -4 -6, 6 -5 -7, 7 -9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -6, 6 -5 -7, 7 -9, -8 2 4 5 9"), 
    ColourLines("1 2 4 5 7, -1 -8, 3 -2 8, -3 -4 -6, 6 -5 -9, -7 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -6, 6 -5 -9, 7 5 4 2 -8, -7 9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 6, -6 5 9, 7 -9, -7 -5 -4 -2 8"), 
    ColourLines("1 8, -1 -2 -4 -5 -7, 3 4 6, -3 2 -8, -6 5 9, 7 -9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 6, -6 5 7, -7 9, 8 -2 -4 -5 -9"), 
    ColourLines("1 8, -1 -2 -4 -5 -9, 3 4 6, -3 2 -8, -6 5 7, -7 9"), 
    ColourLines("1 2 4 6, -1 -8, 3 -2 8, -3 -4 -5 -7, -6 5 9, 7 -9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -7, 6 4 2 -8, -6 5 9, 7 -9"), 
    ColourLines("1 2 4 6, -1 -8, 3 -2 8, -3 -4 -5 -9, -6 5 7, -7 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -9, 6 4 2 -8, -6 5 7, -7 9")
  }; 
  static const ColourLines diag63[16] = { 
    ColourLines("1 2 6, -1 -8, 4 -3 -2 8, -4 -5 -7, -6 3 5 9, 7 -9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -7, 6 2 -8, -6 3 5 9, 7 -9"), 
    ColourLines("1 2 6, -1 -8, 4 -3 -2 8, -4 -5 -9, -6 3 5 7, -7 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -9, 6 2 -8, -6 3 5 7, -7 9"), 
    ColourLines("1 2 3 5 9, -1 -8, 4 -3 6, -4 -5 -7, -6 -2 8, 7 -9"), 
    ColourLines("1 8, -1 -2 -6, 4 -3 6, -4 -5 -7, 7 -9, -8 2 3 5 9"), 
    ColourLines("1 2 3 5 7, -1 -8, 4 -3 6, -4 -5 -9, -6 -2 8, -7 9"), 
    ColourLines("1 8, -1 -2 -6, 4 -3 6, -4 -5 -9, 7 5 3 2 -8, -7 9"), 
    ColourLines("1 2 6, -1 -8, 4 5 9, -4 3 -6, 7 -9, -7 -5 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -5 -7, 4 5 9, -4 3 -6, 6 2 -8, 7 -9"), 
    ColourLines("1 2 6, -1 -8, 4 5 7, -4 3 -6, -7 9, 8 -2 -3 -5 -9"), 
    ColourLines("1 8, -1 -2 -3 -5 -9, 4 5 7, -4 3 -6, 6 2 -8, -7 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 5 9, 6 -3 -5 -7, -6 -2 8, 7 -9"), 
    ColourLines("1 8, -1 -2 -6, 4 5 9, -4 3 2 -8, 6 -3 -5 -7, 7 -9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 5 7, 6 -3 -5 -9, -6 -2 8, -7 9"), 
    ColourLines("1 8, -1 -2 -6, 4 5 7, -4 3 2 -8, 6 -3 -5 -9, -7 9")
  }; 
  static const ColourLines diag64[16] = { 
    ColourLines("1 2 7, -1 -9, 5 -4 8, -5 -6, 6 4 3 -7, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 5 -4 8, -5 -6, 6 4 3 -7, 7 2 -9"), 
    ColourLines("1 2 7, -1 -9, 5 6, -5 4 3 -7, -6 -4 8, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 5 6, -5 4 3 -7, -6 -4 8, 7 2 -9"), 
    ColourLines("1 2 3 4 6, -1 -9, 5 -4 8, -5 -6, 7 -3 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 -4 8, -5 -6, 6 4 3 2 -9, 7 -3 -8"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 6, -6 -4 8, 7 -3 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 6, -5 4 3 2 -9, -6 -4 8, 7 -3 -8"), 
    ColourLines("1 2 7, -1 -9, 5 -4 -3 -2 9, -5 -6, 6 4 -8, -7 3 8"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -6, 6 4 -8, 7 2 -9, -7 3 8"), 
    ColourLines("1 2 7, -1 -9, 5 6, -5 4 -8, -6 -4 -3 -2 9, -7 3 8"), 
    ColourLines("1 9, -1 -2 -3 -4 -6, 5 6, -5 4 -8, 7 2 -9, -7 3 8"), 
    ColourLines("1 2 3 8, -1 -9, 5 -4 -3 7, -5 -6, 6 4 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 -4 -3 7, -5 -6, 6 4 -8, 8 3 2 -9"), 
    ColourLines("1 2 3 8, -1 -9, 5 6, -5 4 -8, -6 -4 -3 7, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 6, -5 4 -8, -6 -4 -3 7, 8 3 2 -9")
  }; 
  static const ColourLines diag65[16] = { 
    ColourLines("1 2 8, -1 -9, 5 -4 7, -5 -6, 6 4 3 -8, -7 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -7, 5 -4 7, -5 -6, 6 4 3 -8, 8 2 -9"), 
    ColourLines("1 2 8, -1 -9, 5 6, -5 4 3 -8, -6 -4 7, -7 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -7, 5 6, -5 4 3 -8, -6 -4 7, 8 2 -9"), 
    ColourLines("1 2 3 4 6, -1 -9, 5 -4 7, -5 -6, -7 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 -4 7, -5 -6, 6 4 3 2 -9, -7 -3 8"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 6, -6 -4 7, -7 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 6, -5 4 3 2 -9, -6 -4 7, -7 -3 8"), 
    ColourLines("1 2 8, -1 -9, 5 -4 -3 -2 9, -5 -6, 6 4 -7, 7 3 -8"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -6, 6 4 -7, 7 3 -8, 8 2 -9"), 
    ColourLines("1 2 8, -1 -9, 5 6, -5 4 -7, -6 -4 -3 -2 9, 7 3 -8"), 
    ColourLines("1 9, -1 -2 -3 -4 -6, 5 6, -5 4 -7, 7 3 -8, 8 2 -9"), 
    ColourLines("1 2 3 7, -1 -9, 5 -4 -3 8, -5 -6, 6 4 -7, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 -4 -3 8, -5 -6, 6 4 -7, 7 3 2 -9"), 
    ColourLines("1 2 3 7, -1 -9, 5 6, -5 4 -7, -6 -4 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 6, -5 4 -7, -6 -4 -3 8, 7 3 2 -9")
  }; 
  static const ColourLines diag66[16] = { 
    ColourLines("1 2 6, -1 -9, 5 -4 8, -5 -7, -6 3 4 7, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 5 -4 8, -5 -7, 6 2 -9, -6 3 4 7"), 
    ColourLines("1 2 6, -1 -9, 5 7, -5 4 3 -6, -7 -4 8, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 5 7, -5 4 3 -6, 6 2 -9, -7 -4 8"), 
    ColourLines("1 2 3 4 7, -1 -9, 5 -4 8, -5 -7, 6 -3 -8, -6 -2 9"), 
    ColourLines("1 9, -1 -2 -6, 5 -4 8, -5 -7, 6 -3 -8, 7 4 3 2 -9"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 7, 6 -3 -8, -6 -2 9, -7 -4 8"), 
    ColourLines("1 9, -1 -2 -6, 5 7, -5 4 3 2 -9, 6 -3 -8, -7 -4 8"), 
    ColourLines("1 2 6, -1 -9, 5 -4 -3 -2 9, -5 -7, -6 3 8, 7 4 -8"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -7, 6 2 -9, -6 3 8, 7 4 -8"), 
    ColourLines("1 2 6, -1 -9, 5 7, -5 4 -8, -6 3 8, -7 -4 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -4 -7, 5 7, -5 4 -8, 6 2 -9, -6 3 8"), 
    ColourLines("1 2 3 8, -1 -9, 5 -4 -3 6, -5 -7, -6 -2 9, 7 4 -8"), 
    ColourLines("1 9, -1 -2 -6, 5 -4 -3 6, -5 -7, 7 4 -8, 8 3 2 -9"), 
    ColourLines("1 2 3 8, -1 -9, 5 7, -5 4 -8, 6 -3 -4 -7, -6 -2 9"), 
    ColourLines("1 9, -1 -2 -6, 5 7, -5 4 -8, 6 -3 -4 -7, 8 3 2 -9")
  }; 
  static const ColourLines diag67[16] = { 
    ColourLines("1 2 8, -1 -9, 5 -4 6, -5 -7, -6 -3 -2 9, 7 4 3 -8"), 
    ColourLines("1 9, -1 -2 -3 -6, 5 -4 6, -5 -7, 7 4 3 -8, 8 2 -9"), 
    ColourLines("1 2 8, -1 -9, 5 7, -5 4 3 -8, 6 -4 -7, -6 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -6, 5 7, -5 4 3 -8, 6 -4 -7, 8 2 -9"), 
    ColourLines("1 2 3 4 7, -1 -9, 5 -4 6, -5 -7, -6 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 -4 6, -5 -7, -6 -3 8, 7 4 3 2 -9"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 7, 6 -4 -7, -6 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 7, -5 4 3 2 -9, 6 -4 -7, -6 -3 8"), 
    ColourLines("1 2 8, -1 -9, 5 -4 -3 -2 9, -5 -7, 6 3 -8, -6 4 7"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -7, 6 3 -8, -6 4 7, 8 2 -9"), 
    ColourLines("1 2 8, -1 -9, 5 7, -5 4 -6, 6 3 -8, -7 -4 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -4 -7, 5 7, -5 4 -6, 6 3 -8, 8 2 -9"), 
    ColourLines("1 2 3 6, -1 -9, 5 -4 -3 8, -5 -7, -6 4 7, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 -4 -3 8, -5 -7, 6 3 2 -9, -6 4 7"), 
    ColourLines("1 2 3 6, -1 -9, 5 7, -5 4 -6, -7 -4 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 7, -5 4 -6, 6 3 2 -9, -7 -4 -3 8")
  }; 
  static const ColourLines diag68[16] = { 
    ColourLines("1 2 6, -1 -9, 5 -4 7, -5 -8, -6 3 4 8, -7 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -7, 5 -4 7, -5 -8, 6 2 -9, -6 3 4 8"), 
    ColourLines("1 2 6, -1 -9, 5 8, -5 4 3 -6, 7 -4 -8, -7 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -7, 5 8, -5 4 3 -6, 6 2 -9, 7 -4 -8"), 
    ColourLines("1 2 3 4 8, -1 -9, 5 -4 7, -5 -8, 6 -3 -7, -6 -2 9"), 
    ColourLines("1 9, -1 -2 -6, 5 -4 7, -5 -8, 6 -3 -7, 8 4 3 2 -9"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 8, 6 -3 -7, -6 -2 9, 7 -4 -8"), 
    ColourLines("1 9, -1 -2 -6, 5 8, -5 4 3 2 -9, 6 -3 -7, 7 -4 -8"), 
    ColourLines("1 2 6, -1 -9, 5 -4 -3 -2 9, -5 -8, -6 3 7, -7 4 8"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -8, 6 2 -9, -6 3 7, -7 4 8"), 
    ColourLines("1 2 6, -1 -9, 5 8, -5 4 -7, -6 3 7, -8 -4 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -4 -8, 5 8, -5 4 -7, 6 2 -9, -6 3 7"), 
    ColourLines("1 2 3 7, -1 -9, 5 -4 -3 6, -5 -8, -6 -2 9, -7 4 8"), 
    ColourLines("1 9, -1 -2 -6, 5 -4 -3 6, -5 -8, 7 3 2 -9, -7 4 8"), 
    ColourLines("1 2 3 7, -1 -9, 5 8, -5 4 -7, 6 -3 -4 -8, -6 -2 9"), 
    ColourLines("1 9, -1 -2 -6, 5 8, -5 4 -7, 6 -3 -4 -8, 7 3 2 -9")
  }; 
  static const ColourLines diag69[16] = { 
    ColourLines("1 2 7, -1 -9, 5 -4 6, -5 -8, -6 -3 -2 9, -7 3 4 8"), 
    ColourLines("1 9, -1 -2 -3 -6, 5 -4 6, -5 -8, 7 2 -9, -7 3 4 8"), 
    ColourLines("1 2 7, -1 -9, 5 8, -5 4 3 -7, 6 -4 -8, -6 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -6, 5 8, -5 4 3 -7, 6 -4 -8, 7 2 -9"), 
    ColourLines("1 2 3 4 8, -1 -9, 5 -4 6, -5 -8, -6 -3 7, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 -4 6, -5 -8, -6 -3 7, 8 4 3 2 -9"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 8, 6 -4 -8, -6 -3 7, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 8, -5 4 3 2 -9, 6 -4 -8, -6 -3 7"), 
    ColourLines("1 2 7, -1 -9, 5 -4 -3 -2 9, -5 -8, 6 3 -7, -6 4 8"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -8, 6 3 -7, -6 4 8, 7 2 -9"), 
    ColourLines("1 2 7, -1 -9, 5 8, -5 4 -6, 6 3 -7, -8 -4 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -4 -8, 5 8, -5 4 -6, 6 3 -7, 7 2 -9"), 
    ColourLines("1 2 3 6, -1 -9, 5 -4 -3 7, -5 -8, -6 4 8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 -4 -3 7, -5 -8, 6 3 2 -9, -6 4 8"), 
    ColourLines("1 2 3 6, -1 -9, 5 8, -5 4 -6, 7 -3 -4 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 8, -5 4 -6, 6 3 2 -9, 7 -3 -4 -8")
  }; 
  static const ColourLines diag70[16] = { 
    ColourLines("1 2 -3, -1 -9, 3 4 5 7, 6 -7, -6 -5 8, -8 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -8, 3 4 5 7, -3 2 -9, 6 -7, -6 -5 8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 5 6, -6 7, -7 -5 8, -8 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -8, 3 4 5 6, -3 2 -9, -6 7, -7 -5 8"), 
    ColourLines("1 2 4 5 7, -1 -9, 3 -2 9, -3 -4 -8, 6 -7, -6 -5 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -8, 6 -7, -6 -5 8, 7 5 4 2 -9"), 
    ColourLines("1 2 4 5 6, -1 -9, 3 -2 9, -3 -4 -8, -6 7, -7 -5 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -8, 6 5 4 2 -9, -6 7, -7 -5 8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 8, 6 -7, -6 -5 -4 -2 9, 7 5 -8"), 
    ColourLines("1 9, -1 -2 -4 -5 -6, 3 4 8, -3 2 -9, 6 -7, 7 5 -8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 8, 6 5 -8, -6 7, -7 -5 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -5 -7, 3 4 8, -3 2 -9, 6 5 -8, -6 7"), 
    ColourLines("1 2 4 8, -1 -9, 3 -2 9, -3 -4 -5 -6, 6 -7, 7 5 -8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -6, 6 -7, 7 5 -8, 8 4 2 -9"), 
    ColourLines("1 2 4 8, -1 -9, 3 -2 9, -3 -4 -5 -7, 6 5 -8, -6 7"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -7, 6 5 -8, -6 7, 8 4 2 -9")
  }; 
  static const ColourLines diag71[16] = { 
    ColourLines("1 2 8, -1 -9, 4 -3 -2 9, -4 -5 -6, 6 -7, 7 5 3 -8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -6, 6 -7, 7 5 3 -8, 8 2 -9"), 
    ColourLines("1 2 8, -1 -9, 4 -3 -2 9, -4 -5 -7, 6 5 3 -8, -6 7"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -7, 6 5 3 -8, -6 7, 8 2 -9"), 
    ColourLines("1 2 3 5 7, -1 -9, 4 -3 8, -4 -5 -6, 6 -7, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 4 -3 8, -4 -5 -6, 6 -7, 7 5 3 2 -9"), 
    ColourLines("1 2 3 5 6, -1 -9, 4 -3 8, -4 -5 -7, -6 7, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 4 -3 8, -4 -5 -7, 6 5 3 2 -9, -6 7"), 
    ColourLines("1 2 8, -1 -9, 4 5 7, -4 3 -8, 6 -7, -6 -5 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -5 -6, 4 5 7, -4 3 -8, 6 -7, 8 2 -9"), 
    ColourLines("1 2 8, -1 -9, 4 5 6, -4 3 -8, -6 7, -7 -5 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -5 -7, 4 5 6, -4 3 -8, -6 7, 8 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 5 7, 6 -7, -6 -5 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 4 5 7, -4 3 2 -9, 6 -7, -6 -5 -3 8"), 
    ColourLines("1 2 3 -4, -1 -9, 4 5 6, -6 7, -7 -5 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 4 5 6, -4 3 2 -9, -6 7, -7 -5 -3 8")
  }; 
  static const ColourLines diag72[16] = { 
    ColourLines("1 2 -3, -1 -9, 3 4 5 8, 6 -8, -6 -5 7, -7 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -7, 3 4 5 8, -3 2 -9, 6 -8, -6 -5 7"), 
    ColourLines("1 2 -3, -1 -9, 3 4 5 6, -6 8, 7 -5 -8, -7 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -7, 3 4 5 6, -3 2 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 2 4 5 8, -1 -9, 3 -2 9, -3 -4 -7, 6 -8, -6 -5 7"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -7, 6 -8, -6 -5 7, 8 5 4 2 -9"), 
    ColourLines("1 2 4 5 6, -1 -9, 3 -2 9, -3 -4 -7, -6 8, 7 -5 -8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -7, 6 5 4 2 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 7, 6 -8, -6 -5 -4 -2 9, -7 5 8"), 
    ColourLines("1 9, -1 -2 -4 -5 -6, 3 4 7, -3 2 -9, 6 -8, -7 5 8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 7, 6 5 -7, -6 8, -8 -5 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -5 -8, 3 4 7, -3 2 -9, 6 5 -7, -6 8"), 
    ColourLines("1 2 4 7, -1 -9, 3 -2 9, -3 -4 -5 -6, 6 -8, -7 5 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -6, 6 -8, 7 4 2 -9, -7 5 8"), 
    ColourLines("1 2 4 7, -1 -9, 3 -2 9, -3 -4 -5 -8, 6 5 -7, -6 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -8, 6 5 -7, -6 8, 7 4 2 -9")
  }; 
  static const ColourLines diag73[16] = { 
    ColourLines("1 2 7, -1 -9, 4 -3 -2 9, -4 -5 -6, 6 -8, -7 3 5 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -6, 6 -8, 7 2 -9, -7 3 5 8"), 
    ColourLines("1 2 7, -1 -9, 4 -3 -2 9, -4 -5 -8, 6 5 3 -7, -6 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -8, 6 5 3 -7, -6 8, 7 2 -9"), 
    ColourLines("1 2 3 5 8, -1 -9, 4 -3 7, -4 -5 -6, 6 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 4 -3 7, -4 -5 -6, 6 -8, 8 5 3 2 -9"), 
    ColourLines("1 2 3 5 6, -1 -9, 4 -3 7, -4 -5 -8, -6 8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 4 -3 7, -4 -5 -8, 6 5 3 2 -9, -6 8"), 
    ColourLines("1 2 7, -1 -9, 4 5 8, -4 3 -7, 6 -8, -6 -5 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -5 -6, 4 5 8, -4 3 -7, 6 -8, 7 2 -9"), 
    ColourLines("1 2 7, -1 -9, 4 5 6, -4 3 -7, -6 8, -8 -5 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -5 -8, 4 5 6, -4 3 -7, -6 8, 7 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 5 8, 6 -8, -6 -5 -3 7, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 4 5 8, -4 3 2 -9, 6 -8, -6 -5 -3 7"), 
    ColourLines("1 2 3 -4, -1 -9, 4 5 6, -6 8, 7 -3 -5 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 4 5 6, -4 3 2 -9, -6 8, 7 -3 -5 -8")
  }; 
  static const ColourLines diag74[16] = { 
    ColourLines("1 2 -3, -1 -9, 3 4 5 8, 6 -5 -7, -6 -4 -2 9, 7 -8"), 
    ColourLines("1 9, -1 -2 -4 -6, 3 4 5 8, -3 2 -9, 6 -5 -7, 7 -8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 5 7, 6 -5 -8, -6 -4 -2 9, -7 8"), 
    ColourLines("1 9, -1 -2 -4 -6, 3 4 5 7, -3 2 -9, 6 -5 -8, -7 8"), 
    ColourLines("1 2 4 5 8, -1 -9, 3 -2 9, -3 -4 -6, 6 -5 -7, 7 -8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -6, 6 -5 -7, 7 -8, 8 5 4 2 -9"), 
    ColourLines("1 2 4 5 7, -1 -9, 3 -2 9, -3 -4 -6, 6 -5 -8, -7 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -6, 6 -5 -8, 7 5 4 2 -9, -7 8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 6, -6 5 8, 7 -8, -7 -5 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -5 -7, 3 4 6, -3 2 -9, -6 5 8, 7 -8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 6, -6 5 7, -7 8, -8 -5 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -5 -8, 3 4 6, -3 2 -9, -6 5 7, -7 8"), 
    ColourLines("1 2 4 6, -1 -9, 3 -2 9, -3 -4 -5 -7, -6 5 8, 7 -8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -7, 6 4 2 -9, -6 5 8, 7 -8"), 
    ColourLines("1 2 4 6, -1 -9, 3 -2 9, -3 -4 -5 -8, -6 5 7, -7 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -8, 6 4 2 -9, -6 5 7, -7 8")
  }; 
  static const ColourLines diag75[16] = { 
    ColourLines("1 2 6, -1 -9, 4 -3 -2 9, -4 -5 -7, -6 3 5 8, 7 -8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -7, 6 2 -9, -6 3 5 8, 7 -8"), 
    ColourLines("1 2 6, -1 -9, 4 -3 -2 9, -4 -5 -8, -6 3 5 7, -7 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -8, 6 2 -9, -6 3 5 7, -7 8"), 
    ColourLines("1 2 3 5 8, -1 -9, 4 -3 6, -4 -5 -7, -6 -2 9, 7 -8"), 
    ColourLines("1 9, -1 -2 -6, 4 -3 6, -4 -5 -7, 7 -8, 8 5 3 2 -9"), 
    ColourLines("1 2 3 5 7, -1 -9, 4 -3 6, -4 -5 -8, -6 -2 9, -7 8"), 
    ColourLines("1 9, -1 -2 -6, 4 -3 6, -4 -5 -8, 7 5 3 2 -9, -7 8"), 
    ColourLines("1 2 6, -1 -9, 4 5 8, -4 3 -6, 7 -8, -7 -5 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -5 -7, 4 5 8, -4 3 -6, 6 2 -9, 7 -8"), 
    ColourLines("1 2 6, -1 -9, 4 5 7, -4 3 -6, -7 8, -8 -5 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -5 -8, 4 5 7, -4 3 -6, 6 2 -9, -7 8"), 
    ColourLines("1 2 3 -4, -1 -9, 4 5 8, 6 -3 -5 -7, -6 -2 9, 7 -8"), 
    ColourLines("1 9, -1 -2 -6, 4 5 8, -4 3 2 -9, 6 -3 -5 -7, 7 -8"), 
    ColourLines("1 2 3 -4, -1 -9, 4 5 7, 6 -3 -5 -8, -6 -2 9, -7 8"), 
    ColourLines("1 9, -1 -2 -6, 4 5 7, -4 3 2 -9, 6 -3 -5 -8, -7 8")
  }; 
  static const ColourLines diag76[16] = { 
    ColourLines("1 4 5 8, -1 -2 3, -3 -6, 6 2 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 4 5 8, -1 -2 -6, 3 6, -3 2 -4 -9, 7 -8, -7 -5 9"), 
    ColourLines("1 4 5 7, -1 -2 3, -3 -6, 6 2 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 4 5 7, -1 -2 -6, 3 6, -3 2 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 2 6, -1 -4 -9, 3 -2 4 5 8, -3 -6, 7 -8, -7 -5 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 6, -6 -2 4 5 8, 7 -8, -7 -5 9"), 
    ColourLines("1 2 6, -1 -4 -9, 3 -2 4 5 7, -3 -6, -7 8, -8 -5 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 6, -6 -2 4 5 7, -7 8, -8 -5 9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -6, 6 2 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 4 9, -1 -2 -6, 3 6, -3 2 -4 -5 -7, 7 -8, 8 5 -9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -6, 6 2 -4 -5 -8, 7 5 -9, -7 8"), 
    ColourLines("1 4 9, -1 -2 -6, 3 6, -3 2 -4 -5 -8, 7 5 -9, -7 8"), 
    ColourLines("1 2 6, -1 -4 -5 -7, 3 -2 4 9, -3 -6, 7 -8, 8 5 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 6, -6 -2 4 9, 7 -8, 8 5 -9"), 
    ColourLines("1 2 6, -1 -4 -5 -8, 3 -2 4 9, -3 -6, 7 5 -9, -7 8"), 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 6, -6 -2 4 9, 7 5 -9, -7 8")
  }; 
  static const ColourLines diag77[16] = { 
    ColourLines("1 2 3 6, -1 -5 -7, 4 -3 9, -4 -6, 7 -8, 8 5 -2 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 6, -6 -3 9, 7 -8, 8 5 -2 -9"), 
    ColourLines("1 2 3 6, -1 -5 -8, 4 -3 9, -4 -6, 7 5 -2 -9, -7 8"), 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 6, -6 -3 9, 7 5 -2 -9, -7 8"), 
    ColourLines("1 2 9, -1 -5 -7, 4 -3 -2 5 8, -4 -6, 6 3 -9, 7 -8"), 
    ColourLines("1 2 9, -1 -5 -7, 4 6, -4 3 -9, -6 -3 -2 5 8, 7 -8"), 
    ColourLines("1 2 9, -1 -5 -8, 4 -3 -2 5 7, -4 -6, 6 3 -9, -7 8"), 
    ColourLines("1 2 9, -1 -5 -8, 4 6, -4 3 -9, -6 -3 -2 5 7, -7 8"), 
    ColourLines("1 5 8, -1 -2 -9, 4 -3 9, -4 -6, 6 3 2 -5 -7, 7 -8"), 
    ColourLines("1 5 8, -1 -2 -9, 4 6, -4 3 2 -5 -7, -6 -3 9, 7 -8"), 
    ColourLines("1 5 7, -1 -2 -9, 4 -3 9, -4 -6, 6 3 2 -5 -8, -7 8"), 
    ColourLines("1 5 7, -1 -2 -9, 4 6, -4 3 2 -5 -8, -6 -3 9, -7 8"), 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -6, 6 3 -9, 7 -8, -7 -5 2 9"), 
    ColourLines("1 5 8, -1 -2 -3 -6, 4 6, -4 3 -9, 7 -8, -7 -5 2 9"), 
    ColourLines("1 5 7, -1 -2 -3 4, -4 -6, 6 3 -9, -7 8, -8 -5 2 9"), 
    ColourLines("1 5 7, -1 -2 -3 -6, 4 6, -4 3 -9, -7 8, -8 -5 2 9")
  }; 
  static const ColourLines diag78[16] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -6, 6 2 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 4 5 9, -1 -2 -6, 3 6, -3 2 -4 -8, 7 -9, -7 -5 8"), 
    ColourLines("1 4 5 7, -1 -2 3, -3 -6, 6 2 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 4 5 7, -1 -2 -6, 3 6, -3 2 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 2 6, -1 -4 -8, 3 -2 4 5 9, -3 -6, 7 -9, -7 -5 8"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 6, -6 -2 4 5 9, 7 -9, -7 -5 8"), 
    ColourLines("1 2 6, -1 -4 -8, 3 -2 4 5 7, -3 -6, -7 9, 8 -5 -9"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 6, -6 -2 4 5 7, -7 9, 8 -5 -9"), 
    ColourLines("1 4 8, -1 -2 3, -3 -6, 6 2 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 4 8, -1 -2 -6, 3 6, -3 2 -4 -5 -7, 7 -9, -8 5 9"), 
    ColourLines("1 4 8, -1 -2 3, -3 -6, 6 2 -4 -5 -9, 7 5 -8, -7 9"), 
    ColourLines("1 4 8, -1 -2 -6, 3 6, -3 2 -4 -5 -9, 7 5 -8, -7 9"), 
    ColourLines("1 2 6, -1 -4 -5 -7, 3 -2 4 8, -3 -6, 7 -9, -8 5 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 6, -6 -2 4 8, 7 -9, -8 5 9"), 
    ColourLines("1 2 6, -1 -4 -5 -9, 3 -2 4 8, -3 -6, 7 5 -8, -7 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 6, -6 -2 4 8, 7 5 -8, -7 9")
  }; 
  static const ColourLines diag79[16] = { 
    ColourLines("1 2 3 6, -1 -5 -7, 4 -3 8, -4 -6, 7 -9, -8 -2 5 9"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 6, -6 -3 8, 7 -9, -8 -2 5 9"), 
    ColourLines("1 2 3 6, -1 -5 -9, 4 -3 8, -4 -6, 7 5 -2 -8, -7 9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 6, -6 -3 8, 7 5 -2 -8, -7 9"), 
    ColourLines("1 2 8, -1 -5 -7, 4 -3 -2 5 9, -4 -6, 6 3 -8, 7 -9"), 
    ColourLines("1 2 8, -1 -5 -7, 4 6, -4 3 -8, -6 -3 -2 5 9, 7 -9"), 
    ColourLines("1 2 8, -1 -5 -9, 4 -3 -2 5 7, -4 -6, 6 3 -8, -7 9"), 
    ColourLines("1 2 8, -1 -5 -9, 4 6, -4 3 -8, -6 -3 -2 5 7, -7 9"), 
    ColourLines("1 5 9, -1 -2 -8, 4 -3 8, -4 -6, 6 3 2 -5 -7, 7 -9"), 
    ColourLines("1 5 9, -1 -2 -8, 4 6, -4 3 2 -5 -7, -6 -3 8, 7 -9"), 
    ColourLines("1 5 7, -1 -2 -8, 4 -3 8, -4 -6, 6 3 2 -5 -9, -7 9"), 
    ColourLines("1 5 7, -1 -2 -8, 4 6, -4 3 2 -5 -9, -6 -3 8, -7 9"), 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -6, 6 3 -8, 7 -9, -7 -5 2 8"), 
    ColourLines("1 5 9, -1 -2 -3 -6, 4 6, -4 3 -8, 7 -9, -7 -5 2 8"), 
    ColourLines("1 5 7, -1 -2 -3 4, -4 -6, 6 3 -8, -7 9, 8 2 -5 -9"), 
    ColourLines("1 5 7, -1 -2 -3 -6, 4 6, -4 3 -8, -7 9, 8 2 -5 -9")
  }; 
  static const ColourLines diag80[16] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -6, 6 2 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 4 5 9, -1 -2 -6, 3 6, -3 2 -4 -7, 7 -5 -8, 8 -9"), 
    ColourLines("1 4 5 8, -1 -2 3, -3 -6, 6 2 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 4 5 8, -1 -2 -6, 3 6, -3 2 -4 -7, 7 -5 -9, -8 9"), 
    ColourLines("1 2 6, -1 -4 -7, 3 -2 4 5 9, -3 -6, 7 -5 -8, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 6, -6 -2 4 5 9, 7 -5 -8, 8 -9"), 
    ColourLines("1 2 6, -1 -4 -7, 3 -2 4 5 8, -3 -6, 7 -5 -9, -8 9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 6, -6 -2 4 5 8, 7 -5 -9, -8 9"), 
    ColourLines("1 4 7, -1 -2 3, -3 -6, 6 2 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 4 7, -1 -2 -6, 3 6, -3 2 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 4 7, -1 -2 3, -3 -6, 6 2 -4 -5 -9, -7 5 8, -8 9"), 
    ColourLines("1 4 7, -1 -2 -6, 3 6, -3 2 -4 -5 -9, -7 5 8, -8 9"), 
    ColourLines("1 2 6, -1 -4 -5 -8, 3 -2 4 7, -3 -6, -7 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 6, -6 -2 4 7, -7 5 9, 8 -9"), 
    ColourLines("1 2 6, -1 -4 -5 -9, 3 -2 4 7, -3 -6, -7 5 8, -8 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 6, -6 -2 4 7, -7 5 8, -8 9")
  }; 
  static const ColourLines diag81[16] = { 
    ColourLines("1 2 3 6, -1 -5 -8, 4 -3 7, -4 -6, -7 -2 5 9, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 6, -6 -3 7, -7 -2 5 9, 8 -9"), 
    ColourLines("1 2 3 6, -1 -5 -9, 4 -3 7, -4 -6, -7 -2 5 8, -8 9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 6, -6 -3 7, -7 -2 5 8, -8 9"), 
    ColourLines("1 2 7, -1 -5 -8, 4 -3 -2 5 9, -4 -6, 6 3 -7, 8 -9"), 
    ColourLines("1 2 7, -1 -5 -8, 4 6, -4 3 -7, -6 -3 -2 5 9, 8 -9"), 
    ColourLines("1 2 7, -1 -5 -9, 4 -3 -2 5 8, -4 -6, 6 3 -7, -8 9"), 
    ColourLines("1 2 7, -1 -5 -9, 4 6, -4 3 -7, -6 -3 -2 5 8, -8 9"), 
    ColourLines("1 5 9, -1 -2 -7, 4 -3 7, -4 -6, 6 3 2 -5 -8, 8 -9"), 
    ColourLines("1 5 9, -1 -2 -7, 4 6, -4 3 2 -5 -8, -6 -3 7, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -7, 4 -3 7, -4 -6, 6 3 2 -5 -9, -8 9"), 
    ColourLines("1 5 8, -1 -2 -7, 4 6, -4 3 2 -5 -9, -6 -3 7, -8 9"), 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -6, 6 3 -7, 7 2 -5 -8, 8 -9"), 
    ColourLines("1 5 9, -1 -2 -3 -6, 4 6, -4 3 -7, 7 2 -5 -8, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -6, 6 3 -7, 7 2 -5 -9, -8 9"), 
    ColourLines("1 5 8, -1 -2 -3 -6, 4 6, -4 3 -7, 7 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag82[16] = { 
    ColourLines("1 4 5 8, -1 -2 3, -3 -7, 6 -8, -6 -5 9, 7 2 -4 -9"), 
    ColourLines("1 4 5 8, -1 -2 -7, 3 7, -3 2 -4 -9, 6 -8, -6 -5 9"), 
    ColourLines("1 4 5 6, -1 -2 3, -3 -7, -6 8, 7 2 -4 -9, -8 -5 9"), 
    ColourLines("1 4 5 6, -1 -2 -7, 3 7, -3 2 -4 -9, -6 8, -8 -5 9"), 
    ColourLines("1 2 7, -1 -4 -9, 3 -2 4 5 8, -3 -7, 6 -8, -6 -5 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 7, 6 -8, -6 -5 9, -7 -2 4 5 8"), 
    ColourLines("1 2 7, -1 -4 -9, 3 -2 4 5 6, -3 -7, -6 8, -8 -5 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 7, 6 5 4 -2 -7, -6 8, -8 -5 9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -7, 6 -8, -6 -5 -4 2 7, 8 5 -9"), 
    ColourLines("1 4 9, -1 -2 -7, 3 7, -3 2 -4 -5 -6, 6 -8, 8 5 -9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -7, 6 5 -9, -6 8, 7 2 -4 -5 -8"), 
    ColourLines("1 4 9, -1 -2 -7, 3 7, -3 2 -4 -5 -8, 6 5 -9, -6 8"), 
    ColourLines("1 2 7, -1 -4 -5 -6, 3 -2 4 9, -3 -7, 6 -8, 8 5 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -6, 3 7, 6 -8, -7 -2 4 9, 8 5 -9"), 
    ColourLines("1 2 7, -1 -4 -5 -8, 3 -2 4 9, -3 -7, 6 5 -9, -6 8"), 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 7, 6 5 -9, -6 8, -7 -2 4 9")
  }; 
  static const ColourLines diag83[16] = { 
    ColourLines("1 2 3 7, -1 -5 -6, 4 -3 9, -4 -7, 6 -8, 8 5 -2 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -6, 4 7, 6 -8, -7 -3 9, 8 5 -2 -9"), 
    ColourLines("1 2 3 7, -1 -5 -8, 4 -3 9, -4 -7, 6 5 -2 -9, -6 8"), 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 7, 6 5 -2 -9, -6 8, -7 -3 9"), 
    ColourLines("1 2 9, -1 -5 -6, 4 -3 -2 5 8, -4 -7, 6 -8, 7 3 -9"), 
    ColourLines("1 2 9, -1 -5 -6, 4 7, -4 3 -9, 6 -8, -7 -3 -2 5 8"), 
    ColourLines("1 2 9, -1 -5 -8, 4 -3 -2 5 6, -4 -7, -6 8, 7 3 -9"), 
    ColourLines("1 2 9, -1 -5 -8, 4 7, -4 3 -9, 6 5 -2 -3 -7, -6 8"), 
    ColourLines("1 5 8, -1 -2 -9, 4 -3 9, -4 -7, 6 -8, -6 -5 2 3 7"), 
    ColourLines("1 5 8, -1 -2 -9, 4 7, -4 3 2 -5 -6, 6 -8, -7 -3 9"), 
    ColourLines("1 5 6, -1 -2 -9, 4 -3 9, -4 -7, -6 8, 7 3 2 -5 -8"), 
    ColourLines("1 5 6, -1 -2 -9, 4 7, -4 3 2 -5 -8, -6 8, -7 -3 9"), 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -7, 6 -8, -6 -5 2 9, 7 3 -9"), 
    ColourLines("1 5 8, -1 -2 -3 -7, 4 7, -4 3 -9, 6 -8, -6 -5 2 9"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -7, -6 8, 7 3 -9, -8 -5 2 9"), 
    ColourLines("1 5 6, -1 -2 -3 -7, 4 7, -4 3 -9, -6 8, -8 -5 2 9")
  }; 
  static const ColourLines diag84[16] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -7, 6 -9, -6 -5 8, 7 2 -4 -8"), 
    ColourLines("1 4 5 9, -1 -2 -7, 3 7, -3 2 -4 -8, 6 -9, -6 -5 8"), 
    ColourLines("1 4 5 6, -1 -2 3, -3 -7, -6 9, 7 2 -4 -8, 8 -5 -9"), 
    ColourLines("1 4 5 6, -1 -2 -7, 3 7, -3 2 -4 -8, -6 9, 8 -5 -9"), 
    ColourLines("1 2 7, -1 -4 -8, 3 -2 4 5 9, -3 -7, 6 -9, -6 -5 8"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 7, 6 -9, -6 -5 8, -7 -2 4 5 9"), 
    ColourLines("1 2 7, -1 -4 -8, 3 -2 4 5 6, -3 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 7, 6 5 4 -2 -7, -6 9, 8 -5 -9"), 
    ColourLines("1 4 8, -1 -2 3, -3 -7, 6 -9, -6 -5 -4 2 7, -8 5 9"), 
    ColourLines("1 4 8, -1 -2 -7, 3 7, -3 2 -4 -5 -6, 6 -9, -8 5 9"), 
    ColourLines("1 4 8, -1 -2 3, -3 -7, 6 5 -8, -6 9, 7 2 -4 -5 -9"), 
    ColourLines("1 4 8, -1 -2 -7, 3 7, -3 2 -4 -5 -9, 6 5 -8, -6 9"), 
    ColourLines("1 2 7, -1 -4 -5 -6, 3 -2 4 8, -3 -7, 6 -9, -8 5 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -6, 3 7, 6 -9, -7 -2 4 8, -8 5 9"), 
    ColourLines("1 2 7, -1 -4 -5 -9, 3 -2 4 8, -3 -7, 6 5 -8, -6 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 7, 6 5 -8, -6 9, -7 -2 4 8")
  }; 
  static const ColourLines diag85[16] = { 
    ColourLines("1 2 3 7, -1 -5 -6, 4 -3 8, -4 -7, 6 -9, -8 -2 5 9"), 
    ColourLines("1 2 3 -4, -1 -5 -6, 4 7, 6 -9, -7 -3 8, -8 -2 5 9"), 
    ColourLines("1 2 3 7, -1 -5 -9, 4 -3 8, -4 -7, 6 5 -2 -8, -6 9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 7, 6 5 -2 -8, -6 9, -7 -3 8"), 
    ColourLines("1 2 8, -1 -5 -6, 4 -3 -2 5 9, -4 -7, 6 -9, 7 3 -8"), 
    ColourLines("1 2 8, -1 -5 -6, 4 7, -4 3 -8, 6 -9, -7 -3 -2 5 9"), 
    ColourLines("1 2 8, -1 -5 -9, 4 -3 -2 5 6, -4 -7, -6 9, 7 3 -8"), 
    ColourLines("1 2 8, -1 -5 -9, 4 7, -4 3 -8, 6 5 -2 -3 -7, -6 9"), 
    ColourLines("1 5 9, -1 -2 -8, 4 -3 8, -4 -7, 6 -9, -6 -5 2 3 7"), 
    ColourLines("1 5 9, -1 -2 -8, 4 7, -4 3 2 -5 -6, 6 -9, -7 -3 8"), 
    ColourLines("1 5 6, -1 -2 -8, 4 -3 8, -4 -7, -6 9, 7 3 2 -5 -9"), 
    ColourLines("1 5 6, -1 -2 -8, 4 7, -4 3 2 -5 -9, -6 9, -7 -3 8"), 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -7, 6 -9, -6 -5 2 8, 7 3 -8"), 
    ColourLines("1 5 9, -1 -2 -3 -7, 4 7, -4 3 -8, 6 -9, -6 -5 2 8"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -7, -6 9, 7 3 -8, 8 2 -5 -9"), 
    ColourLines("1 5 6, -1 -2 -3 -7, 4 7, -4 3 -8, -6 9, 8 2 -5 -9")
  }; 
  static const ColourLines diag86[16] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -7, 6 -5 -8, -6 -4 2 7, 8 -9"), 
    ColourLines("1 4 5 9, -1 -2 -7, 3 7, -3 2 -4 -6, 6 -5 -8, 8 -9"), 
    ColourLines("1 4 5 8, -1 -2 3, -3 -7, 6 -5 -9, -6 -4 2 7, -8 9"), 
    ColourLines("1 4 5 8, -1 -2 -7, 3 7, -3 2 -4 -6, 6 -5 -9, -8 9"), 
    ColourLines("1 2 7, -1 -4 -6, 3 -2 4 5 9, -3 -7, 6 -5 -8, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 7, 6 -5 -8, -7 -2 4 5 9, 8 -9"), 
    ColourLines("1 2 7, -1 -4 -6, 3 -2 4 5 8, -3 -7, 6 -5 -9, -8 9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 7, 6 -5 -9, -7 -2 4 5 8, -8 9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -7, -6 5 9, 7 2 -4 -5 -8, 8 -9"), 
    ColourLines("1 4 6, -1 -2 -7, 3 7, -3 2 -4 -5 -8, -6 5 9, 8 -9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -7, -6 5 8, 7 2 -4 -5 -9, -8 9"), 
    ColourLines("1 4 6, -1 -2 -7, 3 7, -3 2 -4 -5 -9, -6 5 8, -8 9"), 
    ColourLines("1 2 7, -1 -4 -5 -8, 3 -2 4 6, -3 -7, -6 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 7, 6 4 -2 -7, -6 5 9, 8 -9"), 
    ColourLines("1 2 7, -1 -4 -5 -9, 3 -2 4 6, -3 -7, -6 5 8, -8 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 7, 6 4 -2 -7, -6 5 8, -8 9")
  }; 
  static const ColourLines diag87[16] = { 
    ColourLines("1 2 3 7, -1 -5 -8, 4 -3 6, -4 -7, -6 -2 5 9, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 7, 6 -3 -7, -6 -2 5 9, 8 -9"), 
    ColourLines("1 2 3 7, -1 -5 -9, 4 -3 6, -4 -7, -6 -2 5 8, -8 9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 7, 6 -3 -7, -6 -2 5 8, -8 9"), 
    ColourLines("1 2 6, -1 -5 -8, 4 -3 -2 5 9, -4 -7, -6 3 7, 8 -9"), 
    ColourLines("1 2 6, -1 -5 -8, 4 7, -4 3 -6, -7 -3 -2 5 9, 8 -9"), 
    ColourLines("1 2 6, -1 -5 -9, 4 -3 -2 5 8, -4 -7, -6 3 7, -8 9"), 
    ColourLines("1 2 6, -1 -5 -9, 4 7, -4 3 -6, -7 -3 -2 5 8, -8 9"), 
    ColourLines("1 5 9, -1 -2 -6, 4 -3 6, -4 -7, 7 3 2 -5 -8, 8 -9"), 
    ColourLines("1 5 9, -1 -2 -6, 4 7, -4 3 2 -5 -8, 6 -3 -7, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -6, 4 -3 6, -4 -7, 7 3 2 -5 -9, -8 9"), 
    ColourLines("1 5 8, -1 -2 -6, 4 7, -4 3 2 -5 -9, 6 -3 -7, -8 9"), 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -7, 6 2 -5 -8, -6 3 7, 8 -9"), 
    ColourLines("1 5 9, -1 -2 -3 -7, 4 7, -4 3 -6, 6 2 -5 -8, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -7, 6 2 -5 -9, -6 3 7, -8 9"), 
    ColourLines("1 5 8, -1 -2 -3 -7, 4 7, -4 3 -6, 6 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag88[16] = { 
    ColourLines("1 4 5 7, -1 -2 3, -3 -8, 6 -7, -6 -5 9, 8 2 -4 -9"), 
    ColourLines("1 4 5 7, -1 -2 -8, 3 8, -3 2 -4 -9, 6 -7, -6 -5 9"), 
    ColourLines("1 4 5 6, -1 -2 3, -3 -8, -6 7, -7 -5 9, 8 2 -4 -9"), 
    ColourLines("1 4 5 6, -1 -2 -8, 3 8, -3 2 -4 -9, -6 7, -7 -5 9"), 
    ColourLines("1 2 8, -1 -4 -9, 3 -2 4 5 7, -3 -8, 6 -7, -6 -5 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 8, 6 -7, -6 -5 9, 7 5 4 -2 -8"), 
    ColourLines("1 2 8, -1 -4 -9, 3 -2 4 5 6, -3 -8, -6 7, -7 -5 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 8, 6 5 4 -2 -8, -6 7, -7 -5 9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -8, 6 -7, -6 -5 -4 2 8, 7 5 -9"), 
    ColourLines("1 4 9, -1 -2 -8, 3 8, -3 2 -4 -5 -6, 6 -7, 7 5 -9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -8, 6 5 -9, -6 7, -7 -5 -4 2 8"), 
    ColourLines("1 4 9, -1 -2 -8, 3 8, -3 2 -4 -5 -7, 6 5 -9, -6 7"), 
    ColourLines("1 2 8, -1 -4 -5 -6, 3 -2 4 9, -3 -8, 6 -7, 7 5 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -6, 3 8, 6 -7, 7 5 -9, -8 -2 4 9"), 
    ColourLines("1 2 8, -1 -4 -5 -7, 3 -2 4 9, -3 -8, 6 5 -9, -6 7"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 8, 6 5 -9, -6 7, -8 -2 4 9")
  }; 
  static const ColourLines diag89[16] = { 
    ColourLines("1 2 3 8, -1 -5 -6, 4 -3 9, -4 -8, 6 -7, 7 5 -2 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -6, 4 8, 6 -7, 7 5 -2 -9, -8 -3 9"), 
    ColourLines("1 2 3 8, -1 -5 -7, 4 -3 9, -4 -8, 6 5 -2 -9, -6 7"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 8, 6 5 -2 -9, -6 7, -8 -3 9"), 
    ColourLines("1 2 9, -1 -5 -6, 4 -3 -2 5 7, -4 -8, 6 -7, 8 3 -9"), 
    ColourLines("1 2 9, -1 -5 -6, 4 8, -4 3 -9, 6 -7, 7 5 -2 -3 -8"), 
    ColourLines("1 2 9, -1 -5 -7, 4 -3 -2 5 6, -4 -8, -6 7, 8 3 -9"), 
    ColourLines("1 2 9, -1 -5 -7, 4 8, -4 3 -9, 6 5 -2 -3 -8, -6 7"), 
    ColourLines("1 5 7, -1 -2 -9, 4 -3 9, -4 -8, 6 -7, -6 -5 2 3 8"), 
    ColourLines("1 5 7, -1 -2 -9, 4 8, -4 3 2 -5 -6, 6 -7, -8 -3 9"), 
    ColourLines("1 5 6, -1 -2 -9, 4 -3 9, -4 -8, -6 7, -7 -5 2 3 8"), 
    ColourLines("1 5 6, -1 -2 -9, 4 8, -4 3 2 -5 -7, -6 7, -8 -3 9"), 
    ColourLines("1 5 7, -1 -2 -3 4, -4 -8, 6 -7, -6 -5 2 9, 8 3 -9"), 
    ColourLines("1 5 7, -1 -2 -3 -8, 4 8, -4 3 -9, 6 -7, -6 -5 2 9"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -8, -6 7, -7 -5 2 9, 8 3 -9"), 
    ColourLines("1 5 6, -1 -2 -3 -8, 4 8, -4 3 -9, -6 7, -7 -5 2 9")
  }; 
  static const ColourLines diag90[16] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -8, 6 -9, -6 -5 7, -7 -4 2 8"), 
    ColourLines("1 4 5 9, -1 -2 -8, 3 8, -3 2 -4 -7, 6 -9, -6 -5 7"), 
    ColourLines("1 4 5 6, -1 -2 3, -3 -8, -6 9, 7 -5 -9, -7 -4 2 8"), 
    ColourLines("1 4 5 6, -1 -2 -8, 3 8, -3 2 -4 -7, -6 9, 7 -5 -9"), 
    ColourLines("1 2 8, -1 -4 -7, 3 -2 4 5 9, -3 -8, 6 -9, -6 -5 7"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 8, 6 -9, -6 -5 7, -8 -2 4 5 9"), 
    ColourLines("1 2 8, -1 -4 -7, 3 -2 4 5 6, -3 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 8, 6 5 4 -2 -8, -6 9, 7 -5 -9"), 
    ColourLines("1 4 7, -1 -2 3, -3 -8, 6 -9, -6 -5 -4 2 8, -7 5 9"), 
    ColourLines("1 4 7, -1 -2 -8, 3 8, -3 2 -4 -5 -6, 6 -9, -7 5 9"), 
    ColourLines("1 4 7, -1 -2 3, -3 -8, 6 5 -7, -6 9, 8 2 -4 -5 -9"), 
    ColourLines("1 4 7, -1 -2 -8, 3 8, -3 2 -4 -5 -9, 6 5 -7, -6 9"), 
    ColourLines("1 2 8, -1 -4 -5 -6, 3 -2 4 7, -3 -8, 6 -9, -7 5 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -6, 3 8, 6 -9, 7 4 -2 -8, -7 5 9"), 
    ColourLines("1 2 8, -1 -4 -5 -9, 3 -2 4 7, -3 -8, 6 5 -7, -6 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 8, 6 5 -7, -6 9, 7 4 -2 -8")
  }; 
  static const ColourLines diag91[16] = { 
    ColourLines("1 2 3 8, -1 -5 -6, 4 -3 7, -4 -8, 6 -9, -7 -2 5 9"), 
    ColourLines("1 2 3 -4, -1 -5 -6, 4 8, 6 -9, 7 -3 -8, -7 -2 5 9"), 
    ColourLines("1 2 3 8, -1 -5 -9, 4 -3 7, -4 -8, 6 5 -2 -7, -6 9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 8, 6 5 -2 -7, -6 9, 7 -3 -8"), 
    ColourLines("1 2 7, -1 -5 -6, 4 -3 -2 5 9, -4 -8, 6 -9, -7 3 8"), 
    ColourLines("1 2 7, -1 -5 -6, 4 8, -4 3 -7, 6 -9, -8 -3 -2 5 9"), 
    ColourLines("1 2 7, -1 -5 -9, 4 -3 -2 5 6, -4 -8, -6 9, -7 3 8"), 
    ColourLines("1 2 7, -1 -5 -9, 4 8, -4 3 -7, 6 5 -2 -3 -8, -6 9"), 
    ColourLines("1 5 9, -1 -2 -7, 4 -3 7, -4 -8, 6 -9, -6 -5 2 3 8"), 
    ColourLines("1 5 9, -1 -2 -7, 4 8, -4 3 2 -5 -6, 6 -9, 7 -3 -8"), 
    ColourLines("1 5 6, -1 -2 -7, 4 -3 7, -4 -8, -6 9, 8 3 2 -5 -9"), 
    ColourLines("1 5 6, -1 -2 -7, 4 8, -4 3 2 -5 -9, -6 9, 7 -3 -8"), 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -8, 6 -9, -6 -5 2 7, -7 3 8"), 
    ColourLines("1 5 9, -1 -2 -3 -8, 4 8, -4 3 -7, 6 -9, -6 -5 2 7"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -8, -6 9, 7 2 -5 -9, -7 3 8"), 
    ColourLines("1 5 6, -1 -2 -3 -8, 4 8, -4 3 -7, -6 9, 7 2 -5 -9")
  }; 
  static const ColourLines diag92[16] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -8, 6 -5 -7, -6 -4 2 8, 7 -9"), 
    ColourLines("1 4 5 9, -1 -2 -8, 3 8, -3 2 -4 -6, 6 -5 -7, 7 -9"), 
    ColourLines("1 4 5 7, -1 -2 3, -3 -8, 6 -5 -9, -6 -4 2 8, -7 9"), 
    ColourLines("1 4 5 7, -1 -2 -8, 3 8, -3 2 -4 -6, 6 -5 -9, -7 9"), 
    ColourLines("1 2 8, -1 -4 -6, 3 -2 4 5 9, -3 -8, 6 -5 -7, 7 -9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 8, 6 -5 -7, 7 -9, -8 -2 4 5 9"), 
    ColourLines("1 2 8, -1 -4 -6, 3 -2 4 5 7, -3 -8, 6 -5 -9, -7 9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 8, 6 -5 -9, 7 5 4 -2 -8, -7 9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -8, -6 5 9, 7 -9, -7 -5 -4 2 8"), 
    ColourLines("1 4 6, -1 -2 -8, 3 8, -3 2 -4 -5 -7, -6 5 9, 7 -9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -8, -6 5 7, -7 9, 8 2 -4 -5 -9"), 
    ColourLines("1 4 6, -1 -2 -8, 3 8, -3 2 -4 -5 -9, -6 5 7, -7 9"), 
    ColourLines("1 2 8, -1 -4 -5 -7, 3 -2 4 6, -3 -8, -6 5 9, 7 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 8, 6 4 -2 -8, -6 5 9, 7 -9"), 
    ColourLines("1 2 8, -1 -4 -5 -9, 3 -2 4 6, -3 -8, -6 5 7, -7 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 8, 6 4 -2 -8, -6 5 7, -7 9")
  }; 
  static const ColourLines diag93[16] = { 
    ColourLines("1 2 3 8, -1 -5 -7, 4 -3 6, -4 -8, -6 -2 5 9, 7 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 8, 6 -3 -8, -6 -2 5 9, 7 -9"), 
    ColourLines("1 2 3 8, -1 -5 -9, 4 -3 6, -4 -8, -6 -2 5 7, -7 9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 8, 6 -3 -8, -6 -2 5 7, -7 9"), 
    ColourLines("1 2 6, -1 -5 -7, 4 -3 -2 5 9, -4 -8, -6 3 8, 7 -9"), 
    ColourLines("1 2 6, -1 -5 -7, 4 8, -4 3 -6, 7 -9, -8 -3 -2 5 9"), 
    ColourLines("1 2 6, -1 -5 -9, 4 -3 -2 5 7, -4 -8, -6 3 8, -7 9"), 
    ColourLines("1 2 6, -1 -5 -9, 4 8, -4 3 -6, 7 5 -2 -3 -8, -7 9"), 
    ColourLines("1 5 9, -1 -2 -6, 4 -3 6, -4 -8, 7 -9, -7 -5 2 3 8"), 
    ColourLines("1 5 9, -1 -2 -6, 4 8, -4 3 2 -5 -7, 6 -3 -8, 7 -9"), 
    ColourLines("1 5 7, -1 -2 -6, 4 -3 6, -4 -8, -7 9, 8 3 2 -5 -9"), 
    ColourLines("1 5 7, -1 -2 -6, 4 8, -4 3 2 -5 -9, 6 -3 -8, -7 9"), 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -8, 6 2 -5 -7, -6 3 8, 7 -9"), 
    ColourLines("1 5 9, -1 -2 -3 -8, 4 8, -4 3 -6, 6 2 -5 -7, 7 -9"), 
    ColourLines("1 5 7, -1 -2 -3 4, -4 -8, 6 2 -5 -9, -6 3 8, -7 9"), 
    ColourLines("1 5 7, -1 -2 -3 -8, 4 8, -4 3 -6, 6 2 -5 -9, -7 9")
  }; 
  static const ColourLines diag94[16] = { 
    ColourLines("1 4 5 7, -1 -2 3, -3 -9, 6 -7, -6 -5 8, -8 -4 2 9"), 
    ColourLines("1 4 5 7, -1 -2 -9, 3 9, -3 2 -4 -8, 6 -7, -6 -5 8"), 
    ColourLines("1 4 5 6, -1 -2 3, -3 -9, -6 7, -7 -5 8, -8 -4 2 9"), 
    ColourLines("1 4 5 6, -1 -2 -9, 3 9, -3 2 -4 -8, -6 7, -7 -5 8"), 
    ColourLines("1 2 9, -1 -4 -8, 3 -2 4 5 7, -3 -9, 6 -7, -6 -5 8"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 9, 6 -7, -6 -5 8, 7 5 4 -2 -9"), 
    ColourLines("1 2 9, -1 -4 -8, 3 -2 4 5 6, -3 -9, -6 7, -7 -5 8"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 9, 6 5 4 -2 -9, -6 7, -7 -5 8"), 
    ColourLines("1 4 8, -1 -2 3, -3 -9, 6 -7, -6 -5 -4 2 9, 7 5 -8"), 
    ColourLines("1 4 8, -1 -2 -9, 3 9, -3 2 -4 -5 -6, 6 -7, 7 5 -8"), 
    ColourLines("1 4 8, -1 -2 3, -3 -9, 6 5 -8, -6 7, -7 -5 -4 2 9"), 
    ColourLines("1 4 8, -1 -2 -9, 3 9, -3 2 -4 -5 -7, 6 5 -8, -6 7"), 
    ColourLines("1 2 9, -1 -4 -5 -6, 3 -2 4 8, -3 -9, 6 -7, 7 5 -8"), 
    ColourLines("1 2 -3, -1 -4 -5 -6, 3 9, 6 -7, 7 5 -8, 8 4 -2 -9"), 
    ColourLines("1 2 9, -1 -4 -5 -7, 3 -2 4 8, -3 -9, 6 5 -8, -6 7"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 9, 6 5 -8, -6 7, 8 4 -2 -9")
  }; 
  static const ColourLines diag95[16] = { 
    ColourLines("1 2 3 9, -1 -5 -6, 4 -3 8, -4 -9, 6 -7, 7 5 -2 -8"), 
    ColourLines("1 2 3 -4, -1 -5 -6, 4 9, 6 -7, 7 5 -2 -8, 8 -3 -9"), 
    ColourLines("1 2 3 9, -1 -5 -7, 4 -3 8, -4 -9, 6 5 -2 -8, -6 7"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 9, 6 5 -2 -8, -6 7, 8 -3 -9"), 
    ColourLines("1 2 8, -1 -5 -6, 4 -3 -2 5 7, -4 -9, 6 -7, -8 3 9"), 
    ColourLines("1 2 8, -1 -5 -6, 4 9, -4 3 -8, 6 -7, 7 5 -2 -3 -9"), 
    ColourLines("1 2 8, -1 -5 -7, 4 -3 -2 5 6, -4 -9, -6 7, -8 3 9"), 
    ColourLines("1 2 8, -1 -5 -7, 4 9, -4 3 -8, 6 5 -2 -3 -9, -6 7"), 
    ColourLines("1 5 7, -1 -2 -8, 4 -3 8, -4 -9, 6 -7, -6 -5 2 3 9"), 
    ColourLines("1 5 7, -1 -2 -8, 4 9, -4 3 2 -5 -6, 6 -7, 8 -3 -9"), 
    ColourLines("1 5 6, -1 -2 -8, 4 -3 8, -4 -9, -6 7, -7 -5 2 3 9"), 
    ColourLines("1 5 6, -1 -2 -8, 4 9, -4 3 2 -5 -7, -6 7, 8 -3 -9"), 
    ColourLines("1 5 7, -1 -2 -3 4, -4 -9, 6 -7, -6 -5 2 8, -8 3 9"), 
    ColourLines("1 5 7, -1 -2 -3 -9, 4 9, -4 3 -8, 6 -7, -6 -5 2 8"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -9, -6 7, -7 -5 2 8, -8 3 9"), 
    ColourLines("1 5 6, -1 -2 -3 -9, 4 9, -4 3 -8, -6 7, -7 -5 2 8")
  }; 
  static const ColourLines diag96[16] = { 
    ColourLines("1 4 5 8, -1 -2 3, -3 -9, 6 -8, -6 -5 7, -7 -4 2 9"), 
    ColourLines("1 4 5 8, -1 -2 -9, 3 9, -3 2 -4 -7, 6 -8, -6 -5 7"), 
    ColourLines("1 4 5 6, -1 -2 3, -3 -9, -6 8, 7 -5 -8, -7 -4 2 9"), 
    ColourLines("1 4 5 6, -1 -2 -9, 3 9, -3 2 -4 -7, -6 8, 7 -5 -8"), 
    ColourLines("1 2 9, -1 -4 -7, 3 -2 4 5 8, -3 -9, 6 -8, -6 -5 7"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 9, 6 -8, -6 -5 7, 8 5 4 -2 -9"), 
    ColourLines("1 2 9, -1 -4 -7, 3 -2 4 5 6, -3 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 9, 6 5 4 -2 -9, -6 8, 7 -5 -8"), 
    ColourLines("1 4 7, -1 -2 3, -3 -9, 6 -8, -6 -5 -4 2 9, -7 5 8"), 
    ColourLines("1 4 7, -1 -2 -9, 3 9, -3 2 -4 -5 -6, 6 -8, -7 5 8"), 
    ColourLines("1 4 7, -1 -2 3, -3 -9, 6 5 -7, -6 8, -8 -5 -4 2 9"), 
    ColourLines("1 4 7, -1 -2 -9, 3 9, -3 2 -4 -5 -8, 6 5 -7, -6 8"), 
    ColourLines("1 2 9, -1 -4 -5 -6, 3 -2 4 7, -3 -9, 6 -8, -7 5 8"), 
    ColourLines("1 2 -3, -1 -4 -5 -6, 3 9, 6 -8, 7 4 -2 -9, -7 5 8"), 
    ColourLines("1 2 9, -1 -4 -5 -8, 3 -2 4 7, -3 -9, 6 5 -7, -6 8"), 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 9, 6 5 -7, -6 8, 7 4 -2 -9")
  }; 
  static const ColourLines diag97[16] = { 
    ColourLines("1 2 3 9, -1 -5 -6, 4 -3 7, -4 -9, 6 -8, -7 -2 5 8"), 
    ColourLines("1 2 3 -4, -1 -5 -6, 4 9, 6 -8, 7 -3 -9, -7 -2 5 8"), 
    ColourLines("1 2 3 9, -1 -5 -8, 4 -3 7, -4 -9, 6 5 -2 -7, -6 8"), 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 9, 6 5 -2 -7, -6 8, 7 -3 -9"), 
    ColourLines("1 2 7, -1 -5 -6, 4 -3 -2 5 8, -4 -9, 6 -8, -7 3 9"), 
    ColourLines("1 2 7, -1 -5 -6, 4 9, -4 3 -7, 6 -8, 8 5 -2 -3 -9"), 
    ColourLines("1 2 7, -1 -5 -8, 4 -3 -2 5 6, -4 -9, -6 8, -7 3 9"), 
    ColourLines("1 2 7, -1 -5 -8, 4 9, -4 3 -7, 6 5 -2 -3 -9, -6 8"), 
    ColourLines("1 5 8, -1 -2 -7, 4 -3 7, -4 -9, 6 -8, -6 -5 2 3 9"), 
    ColourLines("1 5 8, -1 -2 -7, 4 9, -4 3 2 -5 -6, 6 -8, 7 -3 -9"), 
    ColourLines("1 5 6, -1 -2 -7, 4 -3 7, -4 -9, -6 8, -8 -5 2 3 9"), 
    ColourLines("1 5 6, -1 -2 -7, 4 9, -4 3 2 -5 -8, -6 8, 7 -3 -9"), 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -9, 6 -8, -6 -5 2 7, -7 3 9"), 
    ColourLines("1 5 8, -1 -2 -3 -9, 4 9, -4 3 -7, 6 -8, -6 -5 2 7"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -9, -6 8, 7 2 -5 -8, -7 3 9"), 
    ColourLines("1 5 6, -1 -2 -3 -9, 4 9, -4 3 -7, -6 8, 7 2 -5 -8")
  }; 
  static const ColourLines diag98[16] = { 
    ColourLines("1 4 5 8, -1 -2 3, -3 -9, 6 -5 -7, -6 -4 2 9, 7 -8"), 
    ColourLines("1 4 5 8, -1 -2 -9, 3 9, -3 2 -4 -6, 6 -5 -7, 7 -8"), 
    ColourLines("1 4 5 7, -1 -2 3, -3 -9, 6 -5 -8, -6 -4 2 9, -7 8"), 
    ColourLines("1 4 5 7, -1 -2 -9, 3 9, -3 2 -4 -6, 6 -5 -8, -7 8"), 
    ColourLines("1 2 9, -1 -4 -6, 3 -2 4 5 8, -3 -9, 6 -5 -7, 7 -8"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 9, 6 -5 -7, 7 -8, 8 5 4 -2 -9"), 
    ColourLines("1 2 9, -1 -4 -6, 3 -2 4 5 7, -3 -9, 6 -5 -8, -7 8"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 9, 6 -5 -8, 7 5 4 -2 -9, -7 8"), 
    ColourLines("1 4 6, -1 -2 3, -3 -9, -6 5 8, 7 -8, -7 -5 -4 2 9"), 
    ColourLines("1 4 6, -1 -2 -9, 3 9, -3 2 -4 -5 -7, -6 5 8, 7 -8"), 
    ColourLines("1 4 6, -1 -2 3, -3 -9, -6 5 7, -7 8, -8 -5 -4 2 9"), 
    ColourLines("1 4 6, -1 -2 -9, 3 9, -3 2 -4 -5 -8, -6 5 7, -7 8"), 
    ColourLines("1 2 9, -1 -4 -5 -7, 3 -2 4 6, -3 -9, -6 5 8, 7 -8"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 9, 6 4 -2 -9, -6 5 8, 7 -8"), 
    ColourLines("1 2 9, -1 -4 -5 -8, 3 -2 4 6, -3 -9, -6 5 7, -7 8"), 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 9, 6 4 -2 -9, -6 5 7, -7 8")
  }; 
  static const ColourLines diag99[16] = { 
    ColourLines("1 2 3 9, -1 -5 -7, 4 -3 6, -4 -9, -6 -2 5 8, 7 -8"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 9, 6 -3 -9, -6 -2 5 8, 7 -8"), 
    ColourLines("1 2 3 9, -1 -5 -8, 4 -3 6, -4 -9, -6 -2 5 7, -7 8"), 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 9, 6 -3 -9, -6 -2 5 7, -7 8"), 
    ColourLines("1 2 6, -1 -5 -7, 4 -3 -2 5 8, -4 -9, -6 3 9, 7 -8"), 
    ColourLines("1 2 6, -1 -5 -7, 4 9, -4 3 -6, 7 -8, 8 5 -2 -3 -9"), 
    ColourLines("1 2 6, -1 -5 -8, 4 -3 -2 5 7, -4 -9, -6 3 9, -7 8"), 
    ColourLines("1 2 6, -1 -5 -8, 4 9, -4 3 -6, 7 5 -2 -3 -9, -7 8"), 
    ColourLines("1 5 8, -1 -2 -6, 4 -3 6, -4 -9, 7 -8, -7 -5 2 3 9"), 
    ColourLines("1 5 8, -1 -2 -6, 4 9, -4 3 2 -5 -7, 6 -3 -9, 7 -8"), 
    ColourLines("1 5 7, -1 -2 -6, 4 -3 6, -4 -9, -7 8, -8 -5 2 3 9"), 
    ColourLines("1 5 7, -1 -2 -6, 4 9, -4 3 2 -5 -8, 6 -3 -9, -7 8"), 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -9, 6 2 -5 -7, -6 3 9, 7 -8"), 
    ColourLines("1 5 8, -1 -2 -3 -9, 4 9, -4 3 -6, 6 2 -5 -7, 7 -8"), 
    ColourLines("1 5 7, -1 -2 -3 4, -4 -9, 6 2 -5 -8, -6 3 9, -7 8"), 
    ColourLines("1 5 7, -1 -2 -3 -9, 4 9, -4 3 -6, 6 2 -5 -8, -7 8")
  }; 
  static const ColourLines diag100[16] = { 
    ColourLines("1 2 5 9, -1 -4 -6, 3 -2 4 7, -3 -5 -8, 6 -7, 8 -9"), 
    ColourLines("1 2 5 9, -1 -4 -7, 3 -2 4 6, -3 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 2 5 8, -1 -4 -6, 3 -2 4 7, -3 -5 -9, 6 -7, -8 9"), 
    ColourLines("1 2 5 8, -1 -4 -7, 3 -2 4 6, -3 -5 -9, -6 7, -8 9"), 
    ColourLines("1 4 7, -1 -2 3, -3 -5 -8, 6 -7, -6 -4 2 5 9, 8 -9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -8, -6 7, -7 -4 2 5 9, 8 -9"), 
    ColourLines("1 4 7, -1 -2 3, -3 -5 -9, 6 -7, -6 -4 2 5 8, -8 9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -9, -6 7, -7 -4 2 5 8, -8 9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 5 9, 6 -7, 7 4 -2 -5 -8, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 5 9, 6 4 -2 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 5 8, 6 -7, 7 4 -2 -5 -9, -8 9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 5 8, 6 4 -2 -5 -9, -6 7, -8 9"), 
    ColourLines("1 4 7, -1 -2 -5 -8, 3 5 9, -3 2 -4 -6, 6 -7, 8 -9"), 
    ColourLines("1 4 6, -1 -2 -5 -8, 3 5 9, -3 2 -4 -7, -6 7, 8 -9"), 
    ColourLines("1 4 7, -1 -2 -5 -9, 3 5 8, -3 2 -4 -6, 6 -7, -8 9"), 
    ColourLines("1 4 6, -1 -2 -5 -9, 3 5 8, -3 2 -4 -7, -6 7, -8 9")
  }; 
  static const ColourLines diag101[16] = { 
    ColourLines("1 2 4 7, -1 -5 -8, 3 -2 5 9, -3 -4 -6, 6 -7, 8 -9"), 
    ColourLines("1 2 4 6, -1 -5 -8, 3 -2 5 9, -3 -4 -7, -6 7, 8 -9"), 
    ColourLines("1 2 4 7, -1 -5 -9, 3 -2 5 8, -3 -4 -6, 6 -7, -8 9"), 
    ColourLines("1 2 4 6, -1 -5 -9, 3 -2 5 8, -3 -4 -7, -6 7, -8 9"), 
    ColourLines("1 2 -3, -1 -5 -8, 3 4 7, 6 -7, -6 -4 -2 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -5 -8, 3 4 6, -6 7, -7 -4 -2 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -5 -9, 3 4 7, 6 -7, -6 -4 -2 5 8, -8 9"), 
    ColourLines("1 2 -3, -1 -5 -9, 3 4 6, -6 7, -7 -4 -2 5 8, -8 9"), 
    ColourLines("1 5 9, -1 -2 3, -3 -4 -6, 6 -7, 7 4 2 -5 -8, 8 -9"), 
    ColourLines("1 5 9, -1 -2 3, -3 -4 -7, 6 4 2 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 5 8, -1 -2 3, -3 -4 -6, 6 -7, 7 4 2 -5 -9, -8 9"), 
    ColourLines("1 5 8, -1 -2 3, -3 -4 -7, 6 4 2 -5 -9, -6 7, -8 9"), 
    ColourLines("1 5 9, -1 -2 -4 -6, 3 4 7, -3 2 -5 -8, 6 -7, 8 -9"), 
    ColourLines("1 5 9, -1 -2 -4 -7, 3 4 6, -3 2 -5 -8, -6 7, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -4 -6, 3 4 7, -3 2 -5 -9, 6 -7, -8 9"), 
    ColourLines("1 5 8, -1 -2 -4 -7, 3 4 6, -3 2 -5 -9, -6 7, -8 9")
  }; 
  static const ColourLines diag102[16] = { 
    ColourLines("1 2 5 9, -1 -4 -6, 3 -2 4 8, -3 -5 -7, 6 -8, 7 -9"), 
    ColourLines("1 2 5 9, -1 -4 -8, 3 -2 4 6, -3 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 2 5 7, -1 -4 -6, 3 -2 4 8, -3 -5 -9, 6 -8, -7 9"), 
    ColourLines("1 2 5 7, -1 -4 -8, 3 -2 4 6, -3 -5 -9, -6 8, -7 9"), 
    ColourLines("1 4 8, -1 -2 3, -3 -5 -7, 6 -8, -6 -4 2 5 9, 7 -9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -7, -6 8, 7 -9, -8 -4 2 5 9"), 
    ColourLines("1 4 8, -1 -2 3, -3 -5 -9, 6 -8, -6 -4 2 5 7, -7 9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -9, -6 8, 7 5 2 -4 -8, -7 9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 5 9, 6 -8, 7 -9, -7 -5 -2 4 8"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 5 9, 6 4 -2 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 5 7, 6 -8, -7 9, 8 4 -2 -5 -9"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 5 7, 6 4 -2 -5 -9, -6 8, -7 9"), 
    ColourLines("1 4 8, -1 -2 -5 -7, 3 5 9, -3 2 -4 -6, 6 -8, 7 -9"), 
    ColourLines("1 4 6, -1 -2 -5 -7, 3 5 9, -3 2 -4 -8, -6 8, 7 -9"), 
    ColourLines("1 4 8, -1 -2 -5 -9, 3 5 7, -3 2 -4 -6, 6 -8, -7 9"), 
    ColourLines("1 4 6, -1 -2 -5 -9, 3 5 7, -3 2 -4 -8, -6 8, -7 9")
  }; 
  static const ColourLines diag103[16] = { 
    ColourLines("1 2 4 8, -1 -5 -7, 3 -2 5 9, -3 -4 -6, 6 -8, 7 -9"), 
    ColourLines("1 2 4 6, -1 -5 -7, 3 -2 5 9, -3 -4 -8, -6 8, 7 -9"), 
    ColourLines("1 2 4 8, -1 -5 -9, 3 -2 5 7, -3 -4 -6, 6 -8, -7 9"), 
    ColourLines("1 2 4 6, -1 -5 -9, 3 -2 5 7, -3 -4 -8, -6 8, -7 9"), 
    ColourLines("1 2 -3, -1 -5 -7, 3 4 8, 6 -8, -6 -4 -2 5 9, 7 -9"), 
    ColourLines("1 2 -3, -1 -5 -7, 3 4 6, -6 8, 7 -9, -8 -4 -2 5 9"), 
    ColourLines("1 2 -3, -1 -5 -9, 3 4 8, 6 -8, -6 -4 -2 5 7, -7 9"), 
    ColourLines("1 2 -3, -1 -5 -9, 3 4 6, -6 8, 7 5 -2 -4 -8, -7 9"), 
    ColourLines("1 5 9, -1 -2 3, -3 -4 -6, 6 -8, 7 -9, -7 -5 2 4 8"), 
    ColourLines("1 5 9, -1 -2 3, -3 -4 -8, 6 4 2 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 5 7, -1 -2 3, -3 -4 -6, 6 -8, -7 9, 8 4 2 -5 -9"), 
    ColourLines("1 5 7, -1 -2 3, -3 -4 -8, 6 4 2 -5 -9, -6 8, -7 9"), 
    ColourLines("1 5 9, -1 -2 -4 -6, 3 4 8, -3 2 -5 -7, 6 -8, 7 -9"), 
    ColourLines("1 5 9, -1 -2 -4 -8, 3 4 6, -3 2 -5 -7, -6 8, 7 -9"), 
    ColourLines("1 5 7, -1 -2 -4 -6, 3 4 8, -3 2 -5 -9, 6 -8, -7 9"), 
    ColourLines("1 5 7, -1 -2 -4 -8, 3 4 6, -3 2 -5 -9, -6 8, -7 9")
  }; 
  static const ColourLines diag104[16] = { 
    ColourLines("1 2 5 8, -1 -4 -6, 3 -2 4 9, -3 -5 -7, 6 -9, 7 -8"), 
    ColourLines("1 2 5 8, -1 -4 -9, 3 -2 4 6, -3 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 2 5 7, -1 -4 -6, 3 -2 4 9, -3 -5 -8, 6 -9, -7 8"), 
    ColourLines("1 2 5 7, -1 -4 -9, 3 -2 4 6, -3 -5 -8, -6 9, -7 8"), 
    ColourLines("1 4 9, -1 -2 3, -3 -5 -7, 6 -9, -6 -4 2 5 8, 7 -8"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -7, -6 9, 7 -8, 8 5 2 -4 -9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -5 -8, 6 -9, -6 -4 2 5 7, -7 8"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -8, -6 9, 7 5 2 -4 -9, -7 8"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 5 8, 6 -9, 7 -8, -7 -5 -2 4 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 5 8, 6 4 -2 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 2 -3, -1 -4 -6, 3 5 7, 6 -9, -7 8, -8 -5 -2 4 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 5 7, 6 4 -2 -5 -8, -6 9, -7 8"), 
    ColourLines("1 4 9, -1 -2 -5 -7, 3 5 8, -3 2 -4 -6, 6 -9, 7 -8"), 
    ColourLines("1 4 6, -1 -2 -5 -7, 3 5 8, -3 2 -4 -9, -6 9, 7 -8"), 
    ColourLines("1 4 9, -1 -2 -5 -8, 3 5 7, -3 2 -4 -6, 6 -9, -7 8"), 
    ColourLines("1 4 6, -1 -2 -5 -8, 3 5 7, -3 2 -4 -9, -6 9, -7 8")
  }; 
  static const ColourLines diag105[16] = { 
    ColourLines("1 2 4 9, -1 -5 -7, 3 -2 5 8, -3 -4 -6, 6 -9, 7 -8"), 
    ColourLines("1 2 4 6, -1 -5 -7, 3 -2 5 8, -3 -4 -9, -6 9, 7 -8"), 
    ColourLines("1 2 4 9, -1 -5 -8, 3 -2 5 7, -3 -4 -6, 6 -9, -7 8"), 
    ColourLines("1 2 4 6, -1 -5 -8, 3 -2 5 7, -3 -4 -9, -6 9, -7 8"), 
    ColourLines("1 2 -3, -1 -5 -7, 3 4 9, 6 -9, -6 -4 -2 5 8, 7 -8"), 
    ColourLines("1 2 -3, -1 -5 -7, 3 4 6, -6 9, 7 -8, 8 5 -2 -4 -9"), 
    ColourLines("1 2 -3, -1 -5 -8, 3 4 9, 6 -9, -6 -4 -2 5 7, -7 8"), 
    ColourLines("1 2 -3, -1 -5 -8, 3 4 6, -6 9, 7 5 -2 -4 -9, -7 8"), 
    ColourLines("1 5 8, -1 -2 3, -3 -4 -6, 6 -9, 7 -8, -7 -5 2 4 9"), 
    ColourLines("1 5 8, -1 -2 3, -3 -4 -9, 6 4 2 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 5 7, -1 -2 3, -3 -4 -6, 6 -9, -7 8, -8 -5 2 4 9"), 
    ColourLines("1 5 7, -1 -2 3, -3 -4 -9, 6 4 2 -5 -8, -6 9, -7 8"), 
    ColourLines("1 5 8, -1 -2 -4 -6, 3 4 9, -3 2 -5 -7, 6 -9, 7 -8"), 
    ColourLines("1 5 8, -1 -2 -4 -9, 3 4 6, -3 2 -5 -7, -6 9, 7 -8"), 
    ColourLines("1 5 7, -1 -2 -4 -6, 3 4 9, -3 2 -5 -8, 6 -9, -7 8"), 
    ColourLines("1 5 7, -1 -2 -4 -9, 3 4 6, -3 2 -5 -8, -6 9, -7 8")
  }; 

  static const int cv1at0[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv1at1[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv1at2[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv1at3[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv1at4[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv1at5[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv1at6[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv1at7[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv1at8[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv1at9[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv1at10[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv1at11[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv1at12[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv1at13[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv1at14[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv1at15[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv2at0[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv2at1[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv2at2[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv2at3[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv2at4[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv2at5[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv2at6[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv2at7[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv2at8[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv2at9[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv2at10[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv2at11[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv2at12[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv2at13[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv2at14[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv2at15[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv3at0[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv3at1[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv3at2[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv3at3[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv3at4[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv3at5[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv3at6[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv3at7[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv3at8[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv3at9[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv3at10[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv3at11[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv3at12[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv3at13[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv3at14[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv3at15[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv4at0[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv4at1[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv4at2[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv4at3[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv4at4[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv4at5[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv4at6[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv4at7[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv4at8[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv4at9[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv4at10[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv4at11[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv4at12[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv4at13[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv4at14[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv4at15[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv5at0[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv5at1[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv5at2[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv5at3[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv5at4[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv5at5[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv5at6[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv5at7[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv5at8[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv5at9[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv5at10[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv5at11[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv5at12[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv5at13[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv5at14[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv5at15[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv6at0[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv6at1[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv6at2[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv6at3[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv6at4[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv6at5[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv6at6[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv6at7[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv6at8[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv6at9[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv6at10[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv6at11[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv6at12[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv6at13[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv6at14[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv6at15[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv7at0[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv7at1[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv7at2[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv7at3[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv7at4[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv7at5[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv7at6[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv7at7[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv7at8[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv7at9[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv7at10[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv7at11[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv7at12[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv7at13[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv7at14[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv7at15[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv8at0[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv8at1[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv8at2[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv8at3[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv8at4[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv8at5[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv8at6[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv8at7[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv8at8[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv8at9[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv8at10[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv8at11[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv8at12[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv8at13[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv8at14[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv8at15[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv9at0[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv9at1[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv9at2[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv9at3[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv9at4[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv9at5[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv9at6[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv9at7[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv9at8[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv9at9[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv9at10[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv9at11[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv9at12[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv9at13[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv9at14[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv9at15[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv10at0[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv10at1[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv10at2[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv10at3[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv10at4[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv10at5[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv10at6[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv10at7[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv10at8[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv10at9[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv10at10[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv10at11[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv10at12[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv10at13[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv10at14[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv10at15[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv11at0[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv11at1[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv11at2[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv11at3[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv11at4[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv11at5[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv11at6[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv11at7[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv11at8[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv11at9[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv11at10[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv11at11[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv11at12[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv11at13[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv11at14[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv11at15[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv12at0[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv12at1[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv12at2[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv12at3[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv12at4[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv12at5[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv12at6[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv12at7[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv12at8[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv12at9[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv12at10[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv12at11[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv12at12[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv12at13[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv12at14[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv12at15[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv13at0[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv13at1[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv13at2[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv13at3[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv13at4[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv13at5[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv13at6[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv13at7[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv13at8[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv13at9[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv13at10[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv13at11[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv13at12[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv13at13[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv13at14[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv13at15[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv14at0[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv14at1[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv14at2[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv14at3[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv14at4[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv14at5[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv14at6[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv14at7[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv14at8[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv14at9[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv14at10[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv14at11[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv14at12[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv14at13[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv14at14[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv14at15[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv15at0[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv15at1[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv15at2[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv15at3[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv15at4[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv15at5[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv15at6[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv15at7[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv15at8[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv15at9[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv15at10[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv15at11[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv15at12[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv15at13[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv15at14[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv15at15[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv16at0[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv16at1[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv16at2[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv16at3[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv16at4[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv16at5[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv16at6[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv16at7[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv16at8[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv16at9[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv16at10[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv16at11[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv16at12[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv16at13[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv16at14[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv16at15[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv17at0[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv17at1[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv17at2[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv17at3[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv17at4[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv17at5[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv17at6[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv17at7[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv17at8[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv17at9[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv17at10[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv17at11[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv17at12[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv17at13[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv17at14[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv17at15[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv18at0[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv18at1[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv18at2[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv18at3[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv18at4[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv18at5[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv18at6[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv18at7[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv18at8[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv18at9[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv18at10[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv18at11[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv18at12[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv18at13[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv18at14[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv18at15[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv19at0[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv19at1[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv19at2[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv19at3[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv19at4[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv19at5[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv19at6[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv19at7[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv19at8[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv19at9[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv19at10[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv19at11[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv19at12[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv19at13[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv19at14[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv19at15[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv20at0[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv20at1[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv20at2[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv20at3[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv20at4[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv20at5[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv20at6[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv20at7[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv20at8[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv20at9[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv20at10[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv20at11[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv20at12[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv20at13[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv20at14[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv20at15[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv21at0[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv21at1[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv21at2[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv21at3[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv21at4[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv21at5[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv21at6[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv21at7[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv21at8[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv21at9[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv21at10[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv21at11[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv21at12[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv21at13[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv21at14[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv21at15[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv22at0[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv22at1[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv22at2[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv22at3[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv22at4[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv22at5[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv22at6[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv22at7[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv22at8[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv22at9[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv22at10[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv22at11[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv22at12[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv22at13[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv22at14[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv22at15[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv23at0[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv23at1[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv23at2[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv23at3[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv23at4[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv23at5[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv23at6[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv23at7[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv23at8[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv23at9[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv23at10[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv23at11[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv23at12[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv23at13[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv23at14[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv23at15[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv24at0[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv24at1[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv24at2[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv24at3[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv24at4[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv24at5[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv24at6[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv24at7[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv24at8[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv24at9[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv24at10[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv24at11[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv24at12[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv24at13[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv24at14[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv24at15[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv25at0[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv25at1[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv25at2[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv25at3[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv25at4[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv25at5[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv25at6[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv25at7[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv25at8[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv25at9[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv25at10[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv25at11[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv25at12[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv25at13[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv25at14[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv25at15[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv26at0[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv26at1[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv26at2[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv26at3[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv26at4[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv26at5[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv26at6[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv26at7[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv26at8[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv26at9[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv26at10[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv26at11[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv26at12[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv26at13[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv26at14[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv26at15[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv27at0[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv27at1[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv27at2[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv27at3[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv27at4[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv27at5[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv27at6[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv27at7[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv27at8[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv27at9[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv27at10[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv27at11[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv27at12[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv27at13[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv27at14[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv27at15[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv28at0[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv28at1[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv28at2[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv28at3[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv28at4[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv28at5[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv28at6[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv28at7[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv28at8[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv28at9[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv28at10[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv28at11[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv28at12[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv28at13[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv28at14[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv28at15[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv29at0[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv29at1[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv29at2[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv29at3[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv29at4[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv29at5[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv29at6[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv29at7[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv29at8[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv29at9[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv29at10[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv29at11[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv29at12[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv29at13[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv29at14[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv29at15[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv30at0[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv30at1[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv30at2[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv30at3[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv30at4[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv30at5[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv30at6[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv30at7[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv30at8[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv30at9[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv30at10[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv30at11[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv30at12[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv30at13[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv30at14[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv30at15[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv31at0[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv31at1[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv31at2[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv31at3[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv31at4[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv31at5[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv31at6[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv31at7[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv31at8[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv31at9[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv31at10[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv31at11[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv31at12[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv31at13[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv31at14[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv31at15[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv32at0[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv32at1[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv32at2[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv32at3[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv32at4[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv32at5[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv32at6[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv32at7[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv32at8[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv32at9[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv32at10[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv32at11[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv32at12[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv32at13[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv32at14[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv32at15[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv33at0[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv33at1[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv33at2[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv33at3[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv33at4[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv33at5[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv33at6[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv33at7[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv33at8[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv33at9[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv33at10[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv33at11[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv33at12[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv33at13[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv33at14[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv33at15[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv34at0[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv34at1[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv34at2[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv34at3[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv34at4[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv34at5[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv34at6[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv34at7[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv34at8[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv34at9[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv34at10[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv34at11[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv34at12[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv34at13[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv34at14[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv34at15[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv35at0[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv35at1[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv35at2[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv35at3[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv35at4[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv35at5[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv35at6[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv35at7[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv35at8[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv35at9[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv35at10[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv35at11[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv35at12[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv35at13[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv35at14[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv35at15[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv36at0[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv36at1[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv36at2[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv36at3[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv36at4[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv36at5[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv36at6[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv36at7[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv36at8[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv36at9[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv36at10[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv36at11[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv36at12[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv36at13[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv36at14[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv36at15[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv37at0[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv37at1[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv37at2[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv37at3[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv37at4[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv37at5[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv37at6[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv37at7[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv37at8[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv37at9[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv37at10[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv37at11[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv37at12[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv37at13[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv37at14[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv37at15[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv38at0[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv38at1[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv38at2[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv38at3[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv38at4[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv38at5[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv38at6[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv38at7[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv38at8[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv38at9[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv38at10[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv38at11[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv38at12[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv38at13[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv38at14[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv38at15[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv39at0[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv39at1[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv39at2[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv39at3[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv39at4[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv39at5[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv39at6[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv39at7[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv39at8[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv39at9[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv39at10[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv39at11[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv39at12[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv39at13[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv39at14[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv39at15[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv40at0[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv40at1[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv40at2[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv40at3[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv40at4[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv40at5[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv40at6[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv40at7[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv40at8[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv40at9[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv40at10[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv40at11[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv40at12[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv40at13[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv40at14[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv40at15[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv41at0[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv41at1[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv41at2[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv41at3[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv41at4[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv41at5[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv41at6[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv41at7[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv41at8[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv41at9[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv41at10[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv41at11[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv41at12[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv41at13[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv41at14[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv41at15[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv42at0[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv42at1[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv42at2[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv42at3[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv42at4[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv42at5[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv42at6[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv42at7[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv42at8[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv42at9[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv42at10[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv42at11[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv42at12[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv42at13[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv42at14[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv42at15[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv43at0[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv43at1[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv43at2[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv43at3[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv43at4[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv43at5[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv43at6[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv43at7[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv43at8[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv43at9[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv43at10[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv43at11[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv43at12[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv43at13[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv43at14[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv43at15[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv44at0[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv44at1[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv44at2[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv44at3[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv44at4[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv44at5[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv44at6[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv44at7[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv44at8[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv44at9[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv44at10[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv44at11[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv44at12[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv44at13[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv44at14[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv44at15[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv45at0[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv45at1[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv45at2[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv45at3[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv45at4[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv45at5[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv45at6[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv45at7[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv45at8[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv45at9[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv45at10[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv45at11[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv45at12[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv45at13[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv45at14[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv45at15[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv46at0[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv46at1[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv46at2[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv46at3[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv46at4[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv46at5[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv46at6[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv46at7[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv46at8[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv46at9[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv46at10[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv46at11[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv46at12[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv46at13[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv46at14[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv46at15[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv47at0[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv47at1[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv47at2[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv47at3[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv47at4[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv47at5[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv47at6[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv47at7[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv47at8[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv47at9[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv47at10[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv47at11[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv47at12[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv47at13[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv47at14[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv47at15[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv48at0[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv48at1[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv48at2[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv48at3[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv48at4[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv48at5[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv48at6[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv48at7[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv48at8[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv48at9[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv48at10[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv48at11[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv48at12[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv48at13[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv48at14[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv48at15[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv49at0[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv49at1[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv49at2[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv49at3[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv49at4[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv49at5[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv49at6[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv49at7[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv49at8[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv49at9[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv49at10[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv49at11[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv49at12[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv49at13[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv49at14[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv49at15[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv50at0[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv50at1[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv50at2[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv50at3[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv50at4[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv50at5[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv50at6[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv50at7[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv50at8[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv50at9[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv50at10[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv50at11[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv50at12[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv50at13[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv50at14[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv50at15[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv51at0[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv51at1[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv51at2[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv51at3[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv51at4[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv51at5[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv51at6[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv51at7[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv51at8[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv51at9[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv51at10[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv51at11[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv51at12[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv51at13[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv51at14[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv51at15[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv52at0[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv52at1[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv52at2[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv52at3[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv52at4[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv52at5[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv52at6[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv52at7[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv52at8[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv52at9[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv52at10[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv52at11[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv52at12[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv52at13[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv52at14[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv52at15[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv53at0[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv53at1[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv53at2[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv53at3[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv53at4[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv53at5[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv53at6[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv53at7[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv53at8[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv53at9[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv53at10[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv53at11[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv53at12[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv53at13[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv53at14[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv53at15[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv54at0[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv54at1[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv54at2[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv54at3[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv54at4[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv54at5[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv54at6[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv54at7[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv54at8[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv54at9[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv54at10[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv54at11[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv54at12[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv54at13[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv54at14[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv54at15[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv55at0[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv55at1[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv55at2[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv55at3[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv55at4[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv55at5[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv55at6[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv55at7[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv55at8[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv55at9[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv55at10[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv55at11[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv55at12[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv55at13[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv55at14[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv55at15[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv56at0[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv56at1[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv56at2[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv56at3[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv56at4[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv56at5[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv56at6[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv56at7[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv56at8[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv56at9[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv56at10[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv56at11[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv56at12[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv56at13[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv56at14[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv56at15[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv57at0[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv57at1[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv57at2[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv57at3[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv57at4[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv57at5[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv57at6[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv57at7[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv57at8[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv57at9[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv57at10[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv57at11[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv57at12[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv57at13[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv57at14[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv57at15[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv58at0[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv58at1[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv58at2[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv58at3[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv58at4[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv58at5[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv58at6[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv58at7[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv58at8[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv58at9[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv58at10[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv58at11[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv58at12[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv58at13[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv58at14[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv58at15[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv59at0[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv59at1[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv59at2[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv59at3[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv59at4[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv59at5[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv59at6[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv59at7[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv59at8[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv59at9[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv59at10[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv59at11[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv59at12[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv59at13[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv59at14[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv59at15[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv60at0[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv60at1[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv60at2[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv60at3[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv60at4[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv60at5[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv60at6[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv60at7[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv60at8[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv60at9[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv60at10[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv60at11[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv60at12[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv60at13[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv60at14[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv60at15[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv61at0[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv61at1[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv61at2[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv61at3[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv61at4[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv61at5[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv61at6[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv61at7[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv61at8[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv61at9[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv61at10[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv61at11[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv61at12[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv61at13[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv61at14[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv61at15[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv62at0[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv62at1[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv62at2[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv62at3[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv62at4[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv62at5[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv62at6[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv62at7[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv62at8[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv62at9[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv62at10[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv62at11[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv62at12[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv62at13[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv62at14[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv62at15[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv63at0[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv63at1[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv63at2[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv63at3[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv63at4[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv63at5[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv63at6[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv63at7[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv63at8[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv63at9[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv63at10[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv63at11[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv63at12[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv63at13[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv63at14[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv63at15[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv64at0[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv64at1[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv64at2[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv64at3[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv64at4[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv64at5[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv64at6[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv64at7[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv64at8[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv64at9[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv64at10[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv64at11[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv64at12[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv64at13[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv64at14[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv64at15[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv65at0[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv65at1[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv65at2[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv65at3[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv65at4[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv65at5[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv65at6[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv65at7[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv65at8[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv65at9[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv65at10[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv65at11[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv65at12[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv65at13[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv65at14[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv65at15[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv66at0[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv66at1[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv66at2[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv66at3[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv66at4[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv66at5[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv66at6[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv66at7[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv66at8[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv66at9[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv66at10[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv66at11[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv66at12[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv66at13[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv66at14[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv66at15[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv67at0[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv67at1[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv67at2[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv67at3[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv67at4[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv67at5[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv67at6[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv67at7[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv67at8[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv67at9[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv67at10[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv67at11[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv67at12[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv67at13[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv67at14[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv67at15[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv68at0[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv68at1[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv68at2[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv68at3[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv68at4[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv68at5[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv68at6[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv68at7[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv68at8[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv68at9[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv68at10[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv68at11[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv68at12[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv68at13[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv68at14[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv68at15[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv69at0[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv69at1[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv69at2[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv69at3[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv69at4[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv69at5[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv69at6[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv69at7[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv69at8[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv69at9[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv69at10[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv69at11[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv69at12[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv69at13[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv69at14[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv69at15[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv70at0[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv70at1[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv70at2[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv70at3[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv70at4[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv70at5[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv70at6[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv70at7[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv70at8[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv70at9[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv70at10[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv70at11[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv70at12[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv70at13[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv70at14[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv70at15[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv71at0[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv71at1[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv71at2[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv71at3[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv71at4[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv71at5[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv71at6[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv71at7[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv71at8[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv71at9[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv71at10[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv71at11[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv71at12[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv71at13[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv71at14[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv71at15[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv72at0[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv72at1[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv72at2[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv72at3[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv72at4[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv72at5[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv72at6[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv72at7[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv72at8[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv72at9[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv72at10[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv72at11[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv72at12[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv72at13[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv72at14[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv72at15[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv73at0[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv73at1[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv73at2[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv73at3[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv73at4[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv73at5[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv73at6[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv73at7[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv73at8[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv73at9[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv73at10[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv73at11[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv73at12[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv73at13[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv73at14[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv73at15[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv74at0[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv74at1[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv74at2[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv74at3[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv74at4[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv74at5[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv74at6[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv74at7[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv74at8[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv74at9[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv74at10[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv74at11[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv74at12[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv74at13[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv74at14[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv74at15[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv75at0[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv75at1[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv75at2[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv75at3[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv75at4[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv75at5[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv75at6[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv75at7[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv75at8[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv75at9[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv75at10[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv75at11[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv75at12[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv75at13[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv75at14[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv75at15[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv76at0[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv76at1[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv76at2[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv76at3[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv76at4[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv76at5[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv76at6[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv76at7[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv76at8[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv76at9[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv76at10[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv76at11[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv76at12[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv76at13[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv76at14[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv76at15[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv77at0[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv77at1[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv77at2[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv77at3[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv77at4[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv77at5[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv77at6[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv77at7[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv77at8[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv77at9[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv77at10[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv77at11[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv77at12[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv77at13[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv77at14[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv77at15[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv78at0[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv78at1[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv78at2[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv78at3[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv78at4[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv78at5[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv78at6[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv78at7[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv78at8[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv78at9[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv78at10[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv78at11[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv78at12[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv78at13[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv78at14[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv78at15[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv79at0[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv79at1[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv79at2[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv79at3[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv79at4[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv79at5[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv79at6[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv79at7[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv79at8[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv79at9[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv79at10[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv79at11[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv79at12[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv79at13[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv79at14[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv79at15[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv80at0[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv80at1[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv80at2[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv80at3[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv80at4[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv80at5[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv80at6[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv80at7[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv80at8[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv80at9[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv80at10[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv80at11[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv80at12[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv80at13[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv80at14[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv80at15[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv81at0[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv81at1[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv81at2[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv81at3[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv81at4[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv81at5[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv81at6[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv81at7[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv81at8[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv81at9[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv81at10[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv81at11[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv81at12[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv81at13[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv81at14[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv81at15[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv82at0[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv82at1[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv82at2[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv82at3[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv82at4[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv82at5[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv82at6[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv82at7[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv82at8[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv82at9[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv82at10[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv82at11[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv82at12[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv82at13[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv82at14[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv82at15[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv83at0[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv83at1[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv83at2[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv83at3[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv83at4[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv83at5[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv83at6[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv83at7[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv83at8[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv83at9[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv83at10[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv83at11[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv83at12[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv83at13[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv83at14[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv83at15[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv84at0[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv84at1[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv84at2[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv84at3[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv84at4[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv84at5[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv84at6[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv84at7[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv84at8[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv84at9[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv84at10[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv84at11[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv84at12[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv84at13[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv84at14[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv84at15[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv85at0[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv85at1[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv85at2[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv85at3[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv85at4[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv85at5[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv85at6[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv85at7[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv85at8[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv85at9[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv85at10[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv85at11[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv85at12[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv85at13[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv85at14[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv85at15[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv86at0[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv86at1[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv86at2[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv86at3[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv86at4[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv86at5[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv86at6[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv86at7[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv86at8[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv86at9[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv86at10[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv86at11[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv86at12[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv86at13[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv86at14[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv86at15[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv87at0[]  = { -1, 2, 0, 1, 4, 3, -998};
  static const int cv87at1[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv87at2[]  = { -1, 2, 0, 1, 3, 4, -998};
  static const int cv87at3[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv87at4[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv87at5[]  = { -1, 1, 0, 2, 4, 3, -998};
  static const int cv87at6[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv87at7[]  = { -1, 1, 0, 2, 3, 4, -998};
  static const int cv87at8[]  = { -1, 4, 3, 2, 0, 1, -998};
  static const int cv87at9[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv87at10[]  = { -1, 3, 4, 2, 0, 1, -998};
  static const int cv87at11[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv87at12[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv87at13[]  = { -1, 4, 3, 1, 0, 2, -998};
  static const int cv87at14[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv87at15[]  = { -1, 3, 4, 1, 0, 2, -998};
  static const int cv88at0[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv88at1[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv88at2[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv88at3[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv88at4[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv88at5[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv88at6[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv88at7[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv88at8[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv88at9[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv88at10[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv88at11[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv88at12[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv88at13[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv88at14[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv88at15[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv89at0[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv89at1[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv89at2[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv89at3[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv89at4[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv89at5[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv89at6[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv89at7[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv89at8[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv89at9[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv89at10[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv89at11[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv89at12[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv89at13[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv89at14[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv89at15[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv90at0[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv90at1[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv90at2[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv90at3[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv90at4[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv90at5[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv90at6[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv90at7[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv90at8[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv90at9[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv90at10[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv90at11[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv90at12[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv90at13[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv90at14[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv90at15[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv91at0[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv91at1[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv91at2[]  = { -1, 3, 0, 2, 1, 4, -998};
  static const int cv91at3[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv91at4[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv91at5[]  = { -1, 2, 0, 3, 4, 1, -998};
  static const int cv91at6[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv91at7[]  = { -1, 2, 0, 3, 1, 4, -998};
  static const int cv91at8[]  = { -1, 4, 1, 3, 0, 2, -998};
  static const int cv91at9[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv91at10[]  = { -1, 1, 4, 3, 0, 2, -998};
  static const int cv91at11[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv91at12[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv91at13[]  = { -1, 4, 1, 2, 0, 3, -998};
  static const int cv91at14[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv91at15[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv92at0[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv92at1[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv92at2[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv92at3[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv92at4[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv92at5[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv92at6[]  = { -1, 3, 0, 2, 4, 1, -998};
  static const int cv92at7[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv92at8[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv92at9[]  = { -1, 1, 4, 2, 0, 3, -998};
  static const int cv92at10[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv92at11[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv92at12[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv92at13[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv92at14[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv92at15[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv93at0[]  = { -1, 3, 0, 1, 4, 2, -998};
  static const int cv93at1[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv93at2[]  = { -1, 3, 0, 1, 2, 4, -998};
  static const int cv93at3[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv93at4[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv93at5[]  = { -1, 1, 0, 3, 4, 2, -998};
  static const int cv93at6[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv93at7[]  = { -1, 1, 0, 3, 2, 4, -998};
  static const int cv93at8[]  = { -1, 4, 2, 3, 0, 1, -998};
  static const int cv93at9[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv93at10[]  = { -1, 2, 4, 3, 0, 1, -998};
  static const int cv93at11[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv93at12[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv93at13[]  = { -1, 4, 2, 1, 0, 3, -998};
  static const int cv93at14[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv93at15[]  = { -1, 2, 4, 1, 0, 3, -998};
  static const int cv94at0[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv94at1[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv94at2[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv94at3[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv94at4[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv94at5[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv94at6[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv94at7[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv94at8[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv94at9[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv94at10[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv94at11[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv94at12[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv94at13[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv94at14[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv94at15[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv95at0[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv95at1[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv95at2[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv95at3[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv95at4[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv95at5[]  = { -1, 3, 0, 4, 2, 1, -998};
  static const int cv95at6[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv95at7[]  = { -1, 3, 0, 4, 1, 2, -998};
  static const int cv95at8[]  = { -1, 2, 1, 4, 0, 3, -998};
  static const int cv95at9[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv95at10[]  = { -1, 1, 2, 4, 0, 3, -998};
  static const int cv95at11[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv95at12[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv95at13[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv95at14[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv95at15[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv96at0[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv96at1[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv96at2[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv96at3[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv96at4[]  = { -1, 4, 0, 3, 1, 2, -998};
  static const int cv96at5[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv96at6[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv96at7[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv96at8[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv96at9[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv96at10[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv96at11[]  = { -1, 2, 1, 3, 0, 4, -998};
  static const int cv96at12[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv96at13[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv96at14[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv96at15[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv97at0[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv97at1[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv97at2[]  = { -1, 4, 0, 2, 1, 3, -998};
  static const int cv97at3[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv97at4[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv97at5[]  = { -1, 2, 0, 4, 3, 1, -998};
  static const int cv97at6[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv97at7[]  = { -1, 2, 0, 4, 1, 3, -998};
  static const int cv97at8[]  = { -1, 3, 1, 4, 0, 2, -998};
  static const int cv97at9[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv97at10[]  = { -1, 1, 3, 4, 0, 2, -998};
  static const int cv97at11[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv97at12[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv97at13[]  = { -1, 3, 1, 2, 0, 4, -998};
  static const int cv97at14[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv97at15[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv98at0[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv98at1[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv98at2[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv98at3[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv98at4[]  = { -1, 4, 0, 3, 2, 1, -998};
  static const int cv98at5[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv98at6[]  = { -1, 4, 0, 2, 3, 1, -998};
  static const int cv98at7[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv98at8[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv98at9[]  = { -1, 1, 3, 2, 0, 4, -998};
  static const int cv98at10[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv98at11[]  = { -1, 1, 2, 3, 0, 4, -998};
  static const int cv98at12[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv98at13[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv98at14[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv98at15[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv99at0[]  = { -1, 4, 0, 1, 3, 2, -998};
  static const int cv99at1[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv99at2[]  = { -1, 4, 0, 1, 2, 3, -998};
  static const int cv99at3[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv99at4[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv99at5[]  = { -1, 1, 0, 4, 3, 2, -998};
  static const int cv99at6[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv99at7[]  = { -1, 1, 0, 4, 2, 3, -998};
  static const int cv99at8[]  = { -1, 3, 2, 4, 0, 1, -998};
  static const int cv99at9[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv99at10[]  = { -1, 2, 3, 4, 0, 1, -998};
  static const int cv99at11[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv99at12[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv99at13[]  = { -1, 3, 2, 1, 0, 4, -998};
  static const int cv99at14[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv99at15[]  = { -1, 2, 3, 1, 0, 4, -998};
  static const int cv100at0[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv100at1[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv100at2[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv100at3[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv100at4[]  = { -1, 2, 1, 4, 3, 0, -998};
  static const int cv100at5[]  = { -1, 1, 2, 4, 3, 0, -998};
  static const int cv100at6[]  = { -1, 2, 1, 3, 4, 0, -998};
  static const int cv100at7[]  = { -1, 1, 2, 3, 4, 0, -998};
  static const int cv100at8[]  = { -1, 0, 4, 3, 2, 1, -998};
  static const int cv100at9[]  = { -1, 0, 4, 3, 1, 2, -998};
  static const int cv100at10[]  = { -1, 0, 3, 4, 2, 1, -998};
  static const int cv100at11[]  = { -1, 0, 3, 4, 1, 2, -998};
  static const int cv100at12[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv100at13[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv100at14[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv100at15[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv101at0[]  = { -1, 2, 1, 0, 4, 3, -998};
  static const int cv101at1[]  = { -1, 1, 2, 0, 4, 3, -998};
  static const int cv101at2[]  = { -1, 2, 1, 0, 3, 4, -998};
  static const int cv101at3[]  = { -1, 1, 2, 0, 3, 4, -998};
  static const int cv101at4[]  = { -1, 0, 2, 1, 4, 3, -998};
  static const int cv101at5[]  = { -1, 0, 1, 2, 4, 3, -998};
  static const int cv101at6[]  = { -1, 0, 2, 1, 3, 4, -998};
  static const int cv101at7[]  = { -1, 0, 1, 2, 3, 4, -998};
  static const int cv101at8[]  = { -1, 4, 3, 2, 1, 0, -998};
  static const int cv101at9[]  = { -1, 4, 3, 1, 2, 0, -998};
  static const int cv101at10[]  = { -1, 3, 4, 2, 1, 0, -998};
  static const int cv101at11[]  = { -1, 3, 4, 1, 2, 0, -998};
  static const int cv101at12[]  = { -1, 4, 3, 0, 2, 1, -998};
  static const int cv101at13[]  = { -1, 4, 3, 0, 1, 2, -998};
  static const int cv101at14[]  = { -1, 3, 4, 0, 2, 1, -998};
  static const int cv101at15[]  = { -1, 3, 4, 0, 1, 2, -998};
  static const int cv102at0[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv102at1[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv102at2[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv102at3[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv102at4[]  = { -1, 3, 1, 4, 2, 0, -998};
  static const int cv102at5[]  = { -1, 1, 3, 4, 2, 0, -998};
  static const int cv102at6[]  = { -1, 3, 1, 2, 4, 0, -998};
  static const int cv102at7[]  = { -1, 1, 3, 2, 4, 0, -998};
  static const int cv102at8[]  = { -1, 0, 4, 2, 3, 1, -998};
  static const int cv102at9[]  = { -1, 0, 4, 2, 1, 3, -998};
  static const int cv102at10[]  = { -1, 0, 2, 4, 3, 1, -998};
  static const int cv102at11[]  = { -1, 0, 2, 4, 1, 3, -998};
  static const int cv102at12[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv102at13[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv102at14[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv102at15[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv103at0[]  = { -1, 3, 1, 0, 4, 2, -998};
  static const int cv103at1[]  = { -1, 1, 3, 0, 4, 2, -998};
  static const int cv103at2[]  = { -1, 3, 1, 0, 2, 4, -998};
  static const int cv103at3[]  = { -1, 1, 3, 0, 2, 4, -998};
  static const int cv103at4[]  = { -1, 0, 3, 1, 4, 2, -998};
  static const int cv103at5[]  = { -1, 0, 1, 3, 4, 2, -998};
  static const int cv103at6[]  = { -1, 0, 3, 1, 2, 4, -998};
  static const int cv103at7[]  = { -1, 0, 1, 3, 2, 4, -998};
  static const int cv103at8[]  = { -1, 4, 2, 3, 1, 0, -998};
  static const int cv103at9[]  = { -1, 4, 2, 1, 3, 0, -998};
  static const int cv103at10[]  = { -1, 2, 4, 3, 1, 0, -998};
  static const int cv103at11[]  = { -1, 2, 4, 1, 3, 0, -998};
  static const int cv103at12[]  = { -1, 4, 2, 0, 3, 1, -998};
  static const int cv103at13[]  = { -1, 4, 2, 0, 1, 3, -998};
  static const int cv103at14[]  = { -1, 2, 4, 0, 3, 1, -998};
  static const int cv103at15[]  = { -1, 2, 4, 0, 1, 3, -998};
  static const int cv104at0[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv104at1[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv104at2[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv104at3[]  = { -1, 2, 3, 0, 1, 4, -998};
  static const int cv104at4[]  = { -1, 4, 1, 3, 2, 0, -998};
  static const int cv104at5[]  = { -1, 1, 4, 3, 2, 0, -998};
  static const int cv104at6[]  = { -1, 4, 1, 2, 3, 0, -998};
  static const int cv104at7[]  = { -1, 1, 4, 2, 3, 0, -998};
  static const int cv104at8[]  = { -1, 0, 3, 2, 4, 1, -998};
  static const int cv104at9[]  = { -1, 0, 3, 2, 1, 4, -998};
  static const int cv104at10[]  = { -1, 0, 2, 3, 4, 1, -998};
  static const int cv104at11[]  = { -1, 0, 2, 3, 1, 4, -998};
  static const int cv104at12[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv104at13[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv104at14[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv104at15[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv105at0[]  = { -1, 4, 1, 0, 3, 2, -998};
  static const int cv105at1[]  = { -1, 1, 4, 0, 3, 2, -998};
  static const int cv105at2[]  = { -1, 4, 1, 0, 2, 3, -998};
  static const int cv105at3[]  = { -1, 1, 4, 0, 2, 3, -998};
  static const int cv105at4[]  = { -1, 0, 4, 1, 3, 2, -998};
  static const int cv105at5[]  = { -1, 0, 1, 4, 3, 2, -998};
  static const int cv105at6[]  = { -1, 0, 4, 1, 2, 3, -998};
  static const int cv105at7[]  = { -1, 0, 1, 4, 2, 3, -998};
  static const int cv105at8[]  = { -1, 3, 2, 4, 1, 0, -998};
  static const int cv105at9[]  = { -1, 3, 2, 1, 4, 0, -998};
  static const int cv105at10[]  = { -1, 2, 3, 4, 1, 0, -998};
  static const int cv105at11[]  = { -1, 2, 3, 1, 4, 0, -998};
  static const int cv105at12[]  = { -1, 3, 2, 0, 4, 1, -998};
  static const int cv105at13[]  = { -1, 3, 2, 0, 1, 4, -998};
  static const int cv105at14[]  = { -1, 2, 3, 0, 4, 1, -998};
  static const int cv105at15[]  = { -1, 2, 3, 0, 1, 4, -998};

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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at8, sizeof(cv1at8)/sizeof(int)),  &(diag1[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at9, sizeof(cv1at9)/sizeof(int)),  &(diag1[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at10, sizeof(cv1at10)/sizeof(int)),  &(diag1[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at11, sizeof(cv1at11)/sizeof(int)),  &(diag1[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at12, sizeof(cv1at12)/sizeof(int)),  &(diag1[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at13, sizeof(cv1at13)/sizeof(int)),  &(diag1[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at14, sizeof(cv1at14)/sizeof(int)),  &(diag1[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at15, sizeof(cv1at15)/sizeof(int)),  &(diag1[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at8, sizeof(cv2at8)/sizeof(int)),  &(diag2[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at9, sizeof(cv2at9)/sizeof(int)),  &(diag2[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at10, sizeof(cv2at10)/sizeof(int)),  &(diag2[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at11, sizeof(cv2at11)/sizeof(int)),  &(diag2[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at12, sizeof(cv2at12)/sizeof(int)),  &(diag2[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at13, sizeof(cv2at13)/sizeof(int)),  &(diag2[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at14, sizeof(cv2at14)/sizeof(int)),  &(diag2[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at15, sizeof(cv2at15)/sizeof(int)),  &(diag2[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at8, sizeof(cv3at8)/sizeof(int)),  &(diag3[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at9, sizeof(cv3at9)/sizeof(int)),  &(diag3[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at10, sizeof(cv3at10)/sizeof(int)),  &(diag3[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at11, sizeof(cv3at11)/sizeof(int)),  &(diag3[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at12, sizeof(cv3at12)/sizeof(int)),  &(diag3[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at13, sizeof(cv3at13)/sizeof(int)),  &(diag3[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at14, sizeof(cv3at14)/sizeof(int)),  &(diag3[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at15, sizeof(cv3at15)/sizeof(int)),  &(diag3[15]) );
  } 
  else if( diag->id() == -4 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int)),  &(diag4[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int)),  &(diag4[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at2, sizeof(cv4at2)/sizeof(int)),  &(diag4[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at3, sizeof(cv4at3)/sizeof(int)),  &(diag4[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at4, sizeof(cv4at4)/sizeof(int)),  &(diag4[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at5, sizeof(cv4at5)/sizeof(int)),  &(diag4[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at6, sizeof(cv4at6)/sizeof(int)),  &(diag4[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at7, sizeof(cv4at7)/sizeof(int)),  &(diag4[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at8, sizeof(cv4at8)/sizeof(int)),  &(diag4[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at9, sizeof(cv4at9)/sizeof(int)),  &(diag4[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at10, sizeof(cv4at10)/sizeof(int)),  &(diag4[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at11, sizeof(cv4at11)/sizeof(int)),  &(diag4[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at12, sizeof(cv4at12)/sizeof(int)),  &(diag4[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at13, sizeof(cv4at13)/sizeof(int)),  &(diag4[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at14, sizeof(cv4at14)/sizeof(int)),  &(diag4[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at15, sizeof(cv4at15)/sizeof(int)),  &(diag4[15]) );
  } 
  else if( diag->id() == -5 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int)),  &(diag5[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int)),  &(diag5[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at2, sizeof(cv5at2)/sizeof(int)),  &(diag5[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at3, sizeof(cv5at3)/sizeof(int)),  &(diag5[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at4, sizeof(cv5at4)/sizeof(int)),  &(diag5[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at5, sizeof(cv5at5)/sizeof(int)),  &(diag5[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at6, sizeof(cv5at6)/sizeof(int)),  &(diag5[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at7, sizeof(cv5at7)/sizeof(int)),  &(diag5[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at8, sizeof(cv5at8)/sizeof(int)),  &(diag5[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at9, sizeof(cv5at9)/sizeof(int)),  &(diag5[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at10, sizeof(cv5at10)/sizeof(int)),  &(diag5[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at11, sizeof(cv5at11)/sizeof(int)),  &(diag5[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at12, sizeof(cv5at12)/sizeof(int)),  &(diag5[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at13, sizeof(cv5at13)/sizeof(int)),  &(diag5[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at14, sizeof(cv5at14)/sizeof(int)),  &(diag5[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at15, sizeof(cv5at15)/sizeof(int)),  &(diag5[15]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int)),  &(diag6[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int)),  &(diag6[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at2, sizeof(cv6at2)/sizeof(int)),  &(diag6[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at3, sizeof(cv6at3)/sizeof(int)),  &(diag6[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at4, sizeof(cv6at4)/sizeof(int)),  &(diag6[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at5, sizeof(cv6at5)/sizeof(int)),  &(diag6[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at6, sizeof(cv6at6)/sizeof(int)),  &(diag6[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at7, sizeof(cv6at7)/sizeof(int)),  &(diag6[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at8, sizeof(cv6at8)/sizeof(int)),  &(diag6[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at9, sizeof(cv6at9)/sizeof(int)),  &(diag6[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at10, sizeof(cv6at10)/sizeof(int)),  &(diag6[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at11, sizeof(cv6at11)/sizeof(int)),  &(diag6[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at12, sizeof(cv6at12)/sizeof(int)),  &(diag6[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at13, sizeof(cv6at13)/sizeof(int)),  &(diag6[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at14, sizeof(cv6at14)/sizeof(int)),  &(diag6[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at15, sizeof(cv6at15)/sizeof(int)),  &(diag6[15]) );
  } 
  else if( diag->id() == -7 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int)),  &(diag7[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at1, sizeof(cv7at1)/sizeof(int)),  &(diag7[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at2, sizeof(cv7at2)/sizeof(int)),  &(diag7[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at3, sizeof(cv7at3)/sizeof(int)),  &(diag7[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at4, sizeof(cv7at4)/sizeof(int)),  &(diag7[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at5, sizeof(cv7at5)/sizeof(int)),  &(diag7[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at6, sizeof(cv7at6)/sizeof(int)),  &(diag7[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at7, sizeof(cv7at7)/sizeof(int)),  &(diag7[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at8, sizeof(cv7at8)/sizeof(int)),  &(diag7[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at9, sizeof(cv7at9)/sizeof(int)),  &(diag7[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at10, sizeof(cv7at10)/sizeof(int)),  &(diag7[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at11, sizeof(cv7at11)/sizeof(int)),  &(diag7[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at12, sizeof(cv7at12)/sizeof(int)),  &(diag7[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at13, sizeof(cv7at13)/sizeof(int)),  &(diag7[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at14, sizeof(cv7at14)/sizeof(int)),  &(diag7[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at15, sizeof(cv7at15)/sizeof(int)),  &(diag7[15]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int)),  &(diag8[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int)),  &(diag8[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at2, sizeof(cv8at2)/sizeof(int)),  &(diag8[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at3, sizeof(cv8at3)/sizeof(int)),  &(diag8[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at4, sizeof(cv8at4)/sizeof(int)),  &(diag8[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at5, sizeof(cv8at5)/sizeof(int)),  &(diag8[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at6, sizeof(cv8at6)/sizeof(int)),  &(diag8[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at7, sizeof(cv8at7)/sizeof(int)),  &(diag8[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at8, sizeof(cv8at8)/sizeof(int)),  &(diag8[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at9, sizeof(cv8at9)/sizeof(int)),  &(diag8[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at10, sizeof(cv8at10)/sizeof(int)),  &(diag8[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at11, sizeof(cv8at11)/sizeof(int)),  &(diag8[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at12, sizeof(cv8at12)/sizeof(int)),  &(diag8[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at13, sizeof(cv8at13)/sizeof(int)),  &(diag8[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at14, sizeof(cv8at14)/sizeof(int)),  &(diag8[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at15, sizeof(cv8at15)/sizeof(int)),  &(diag8[15]) );
  } 
  else if( diag->id() == -9 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int)),  &(diag9[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int)),  &(diag9[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at2, sizeof(cv9at2)/sizeof(int)),  &(diag9[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at3, sizeof(cv9at3)/sizeof(int)),  &(diag9[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at4, sizeof(cv9at4)/sizeof(int)),  &(diag9[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at5, sizeof(cv9at5)/sizeof(int)),  &(diag9[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at6, sizeof(cv9at6)/sizeof(int)),  &(diag9[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at7, sizeof(cv9at7)/sizeof(int)),  &(diag9[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at8, sizeof(cv9at8)/sizeof(int)),  &(diag9[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at9, sizeof(cv9at9)/sizeof(int)),  &(diag9[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at10, sizeof(cv9at10)/sizeof(int)),  &(diag9[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at11, sizeof(cv9at11)/sizeof(int)),  &(diag9[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at12, sizeof(cv9at12)/sizeof(int)),  &(diag9[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at13, sizeof(cv9at13)/sizeof(int)),  &(diag9[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at14, sizeof(cv9at14)/sizeof(int)),  &(diag9[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at15, sizeof(cv9at15)/sizeof(int)),  &(diag9[15]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)),  &(diag10[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int)),  &(diag10[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at2, sizeof(cv10at2)/sizeof(int)),  &(diag10[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at3, sizeof(cv10at3)/sizeof(int)),  &(diag10[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at4, sizeof(cv10at4)/sizeof(int)),  &(diag10[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at5, sizeof(cv10at5)/sizeof(int)),  &(diag10[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at6, sizeof(cv10at6)/sizeof(int)),  &(diag10[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at7, sizeof(cv10at7)/sizeof(int)),  &(diag10[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at8, sizeof(cv10at8)/sizeof(int)),  &(diag10[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at9, sizeof(cv10at9)/sizeof(int)),  &(diag10[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at10, sizeof(cv10at10)/sizeof(int)),  &(diag10[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at11, sizeof(cv10at11)/sizeof(int)),  &(diag10[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at12, sizeof(cv10at12)/sizeof(int)),  &(diag10[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at13, sizeof(cv10at13)/sizeof(int)),  &(diag10[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at14, sizeof(cv10at14)/sizeof(int)),  &(diag10[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at15, sizeof(cv10at15)/sizeof(int)),  &(diag10[15]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)),  &(diag11[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int)),  &(diag11[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at2, sizeof(cv11at2)/sizeof(int)),  &(diag11[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at3, sizeof(cv11at3)/sizeof(int)),  &(diag11[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at4, sizeof(cv11at4)/sizeof(int)),  &(diag11[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at5, sizeof(cv11at5)/sizeof(int)),  &(diag11[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at6, sizeof(cv11at6)/sizeof(int)),  &(diag11[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at7, sizeof(cv11at7)/sizeof(int)),  &(diag11[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at8, sizeof(cv11at8)/sizeof(int)),  &(diag11[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at9, sizeof(cv11at9)/sizeof(int)),  &(diag11[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at10, sizeof(cv11at10)/sizeof(int)),  &(diag11[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at11, sizeof(cv11at11)/sizeof(int)),  &(diag11[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at12, sizeof(cv11at12)/sizeof(int)),  &(diag11[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at13, sizeof(cv11at13)/sizeof(int)),  &(diag11[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at14, sizeof(cv11at14)/sizeof(int)),  &(diag11[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at15, sizeof(cv11at15)/sizeof(int)),  &(diag11[15]) );
  } 
  else if( diag->id() == -12 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int)),  &(diag12[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at1, sizeof(cv12at1)/sizeof(int)),  &(diag12[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at2, sizeof(cv12at2)/sizeof(int)),  &(diag12[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at3, sizeof(cv12at3)/sizeof(int)),  &(diag12[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at4, sizeof(cv12at4)/sizeof(int)),  &(diag12[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at5, sizeof(cv12at5)/sizeof(int)),  &(diag12[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at6, sizeof(cv12at6)/sizeof(int)),  &(diag12[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at7, sizeof(cv12at7)/sizeof(int)),  &(diag12[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at8, sizeof(cv12at8)/sizeof(int)),  &(diag12[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at9, sizeof(cv12at9)/sizeof(int)),  &(diag12[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at10, sizeof(cv12at10)/sizeof(int)),  &(diag12[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at11, sizeof(cv12at11)/sizeof(int)),  &(diag12[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at12, sizeof(cv12at12)/sizeof(int)),  &(diag12[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at13, sizeof(cv12at13)/sizeof(int)),  &(diag12[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at14, sizeof(cv12at14)/sizeof(int)),  &(diag12[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at15, sizeof(cv12at15)/sizeof(int)),  &(diag12[15]) );
  } 
  else if( diag->id() == -13 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int)),  &(diag13[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at1, sizeof(cv13at1)/sizeof(int)),  &(diag13[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at2, sizeof(cv13at2)/sizeof(int)),  &(diag13[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at3, sizeof(cv13at3)/sizeof(int)),  &(diag13[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at4, sizeof(cv13at4)/sizeof(int)),  &(diag13[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at5, sizeof(cv13at5)/sizeof(int)),  &(diag13[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at6, sizeof(cv13at6)/sizeof(int)),  &(diag13[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at7, sizeof(cv13at7)/sizeof(int)),  &(diag13[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at8, sizeof(cv13at8)/sizeof(int)),  &(diag13[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at9, sizeof(cv13at9)/sizeof(int)),  &(diag13[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at10, sizeof(cv13at10)/sizeof(int)),  &(diag13[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at11, sizeof(cv13at11)/sizeof(int)),  &(diag13[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at12, sizeof(cv13at12)/sizeof(int)),  &(diag13[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at13, sizeof(cv13at13)/sizeof(int)),  &(diag13[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at14, sizeof(cv13at14)/sizeof(int)),  &(diag13[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at15, sizeof(cv13at15)/sizeof(int)),  &(diag13[15]) );
  } 
  else if( diag->id() == -14 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int)),  &(diag14[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at1, sizeof(cv14at1)/sizeof(int)),  &(diag14[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at2, sizeof(cv14at2)/sizeof(int)),  &(diag14[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at3, sizeof(cv14at3)/sizeof(int)),  &(diag14[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at4, sizeof(cv14at4)/sizeof(int)),  &(diag14[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at5, sizeof(cv14at5)/sizeof(int)),  &(diag14[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at6, sizeof(cv14at6)/sizeof(int)),  &(diag14[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at7, sizeof(cv14at7)/sizeof(int)),  &(diag14[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at8, sizeof(cv14at8)/sizeof(int)),  &(diag14[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at9, sizeof(cv14at9)/sizeof(int)),  &(diag14[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at10, sizeof(cv14at10)/sizeof(int)),  &(diag14[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at11, sizeof(cv14at11)/sizeof(int)),  &(diag14[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at12, sizeof(cv14at12)/sizeof(int)),  &(diag14[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at13, sizeof(cv14at13)/sizeof(int)),  &(diag14[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at14, sizeof(cv14at14)/sizeof(int)),  &(diag14[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at15, sizeof(cv14at15)/sizeof(int)),  &(diag14[15]) );
  } 
  else if( diag->id() == -15 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int)),  &(diag15[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int)),  &(diag15[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at2, sizeof(cv15at2)/sizeof(int)),  &(diag15[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at3, sizeof(cv15at3)/sizeof(int)),  &(diag15[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at4, sizeof(cv15at4)/sizeof(int)),  &(diag15[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at5, sizeof(cv15at5)/sizeof(int)),  &(diag15[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at6, sizeof(cv15at6)/sizeof(int)),  &(diag15[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at7, sizeof(cv15at7)/sizeof(int)),  &(diag15[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at8, sizeof(cv15at8)/sizeof(int)),  &(diag15[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at9, sizeof(cv15at9)/sizeof(int)),  &(diag15[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at10, sizeof(cv15at10)/sizeof(int)),  &(diag15[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at11, sizeof(cv15at11)/sizeof(int)),  &(diag15[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at12, sizeof(cv15at12)/sizeof(int)),  &(diag15[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at13, sizeof(cv15at13)/sizeof(int)),  &(diag15[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at14, sizeof(cv15at14)/sizeof(int)),  &(diag15[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at15, sizeof(cv15at15)/sizeof(int)),  &(diag15[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at8, sizeof(cv16at8)/sizeof(int)),  &(diag16[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at9, sizeof(cv16at9)/sizeof(int)),  &(diag16[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at10, sizeof(cv16at10)/sizeof(int)),  &(diag16[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at11, sizeof(cv16at11)/sizeof(int)),  &(diag16[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at12, sizeof(cv16at12)/sizeof(int)),  &(diag16[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at13, sizeof(cv16at13)/sizeof(int)),  &(diag16[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at14, sizeof(cv16at14)/sizeof(int)),  &(diag16[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at15, sizeof(cv16at15)/sizeof(int)),  &(diag16[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at8, sizeof(cv17at8)/sizeof(int)),  &(diag17[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at9, sizeof(cv17at9)/sizeof(int)),  &(diag17[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at10, sizeof(cv17at10)/sizeof(int)),  &(diag17[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at11, sizeof(cv17at11)/sizeof(int)),  &(diag17[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at12, sizeof(cv17at12)/sizeof(int)),  &(diag17[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at13, sizeof(cv17at13)/sizeof(int)),  &(diag17[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at14, sizeof(cv17at14)/sizeof(int)),  &(diag17[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at15, sizeof(cv17at15)/sizeof(int)),  &(diag17[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at8, sizeof(cv18at8)/sizeof(int)),  &(diag18[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at9, sizeof(cv18at9)/sizeof(int)),  &(diag18[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at10, sizeof(cv18at10)/sizeof(int)),  &(diag18[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at11, sizeof(cv18at11)/sizeof(int)),  &(diag18[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at12, sizeof(cv18at12)/sizeof(int)),  &(diag18[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at13, sizeof(cv18at13)/sizeof(int)),  &(diag18[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at14, sizeof(cv18at14)/sizeof(int)),  &(diag18[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at15, sizeof(cv18at15)/sizeof(int)),  &(diag18[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at8, sizeof(cv19at8)/sizeof(int)),  &(diag19[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at9, sizeof(cv19at9)/sizeof(int)),  &(diag19[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at10, sizeof(cv19at10)/sizeof(int)),  &(diag19[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at11, sizeof(cv19at11)/sizeof(int)),  &(diag19[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at12, sizeof(cv19at12)/sizeof(int)),  &(diag19[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at13, sizeof(cv19at13)/sizeof(int)),  &(diag19[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at14, sizeof(cv19at14)/sizeof(int)),  &(diag19[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at15, sizeof(cv19at15)/sizeof(int)),  &(diag19[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at8, sizeof(cv20at8)/sizeof(int)),  &(diag20[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at9, sizeof(cv20at9)/sizeof(int)),  &(diag20[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at10, sizeof(cv20at10)/sizeof(int)),  &(diag20[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at11, sizeof(cv20at11)/sizeof(int)),  &(diag20[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at12, sizeof(cv20at12)/sizeof(int)),  &(diag20[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at13, sizeof(cv20at13)/sizeof(int)),  &(diag20[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at14, sizeof(cv20at14)/sizeof(int)),  &(diag20[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at15, sizeof(cv20at15)/sizeof(int)),  &(diag20[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at8, sizeof(cv21at8)/sizeof(int)),  &(diag21[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at9, sizeof(cv21at9)/sizeof(int)),  &(diag21[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at10, sizeof(cv21at10)/sizeof(int)),  &(diag21[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at11, sizeof(cv21at11)/sizeof(int)),  &(diag21[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at12, sizeof(cv21at12)/sizeof(int)),  &(diag21[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at13, sizeof(cv21at13)/sizeof(int)),  &(diag21[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at14, sizeof(cv21at14)/sizeof(int)),  &(diag21[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at15, sizeof(cv21at15)/sizeof(int)),  &(diag21[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at8, sizeof(cv22at8)/sizeof(int)),  &(diag22[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at9, sizeof(cv22at9)/sizeof(int)),  &(diag22[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at10, sizeof(cv22at10)/sizeof(int)),  &(diag22[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at11, sizeof(cv22at11)/sizeof(int)),  &(diag22[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at12, sizeof(cv22at12)/sizeof(int)),  &(diag22[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at13, sizeof(cv22at13)/sizeof(int)),  &(diag22[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at14, sizeof(cv22at14)/sizeof(int)),  &(diag22[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at15, sizeof(cv22at15)/sizeof(int)),  &(diag22[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at8, sizeof(cv23at8)/sizeof(int)),  &(diag23[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at9, sizeof(cv23at9)/sizeof(int)),  &(diag23[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at10, sizeof(cv23at10)/sizeof(int)),  &(diag23[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at11, sizeof(cv23at11)/sizeof(int)),  &(diag23[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at12, sizeof(cv23at12)/sizeof(int)),  &(diag23[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at13, sizeof(cv23at13)/sizeof(int)),  &(diag23[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at14, sizeof(cv23at14)/sizeof(int)),  &(diag23[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at15, sizeof(cv23at15)/sizeof(int)),  &(diag23[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at8, sizeof(cv24at8)/sizeof(int)),  &(diag24[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at9, sizeof(cv24at9)/sizeof(int)),  &(diag24[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at10, sizeof(cv24at10)/sizeof(int)),  &(diag24[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at11, sizeof(cv24at11)/sizeof(int)),  &(diag24[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at12, sizeof(cv24at12)/sizeof(int)),  &(diag24[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at13, sizeof(cv24at13)/sizeof(int)),  &(diag24[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at14, sizeof(cv24at14)/sizeof(int)),  &(diag24[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at15, sizeof(cv24at15)/sizeof(int)),  &(diag24[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at8, sizeof(cv25at8)/sizeof(int)),  &(diag25[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at9, sizeof(cv25at9)/sizeof(int)),  &(diag25[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at10, sizeof(cv25at10)/sizeof(int)),  &(diag25[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at11, sizeof(cv25at11)/sizeof(int)),  &(diag25[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at12, sizeof(cv25at12)/sizeof(int)),  &(diag25[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at13, sizeof(cv25at13)/sizeof(int)),  &(diag25[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at14, sizeof(cv25at14)/sizeof(int)),  &(diag25[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at15, sizeof(cv25at15)/sizeof(int)),  &(diag25[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at8, sizeof(cv26at8)/sizeof(int)),  &(diag26[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at9, sizeof(cv26at9)/sizeof(int)),  &(diag26[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at10, sizeof(cv26at10)/sizeof(int)),  &(diag26[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at11, sizeof(cv26at11)/sizeof(int)),  &(diag26[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at12, sizeof(cv26at12)/sizeof(int)),  &(diag26[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at13, sizeof(cv26at13)/sizeof(int)),  &(diag26[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at14, sizeof(cv26at14)/sizeof(int)),  &(diag26[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at15, sizeof(cv26at15)/sizeof(int)),  &(diag26[15]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at8, sizeof(cv27at8)/sizeof(int)),  &(diag27[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at9, sizeof(cv27at9)/sizeof(int)),  &(diag27[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at10, sizeof(cv27at10)/sizeof(int)),  &(diag27[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at11, sizeof(cv27at11)/sizeof(int)),  &(diag27[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at12, sizeof(cv27at12)/sizeof(int)),  &(diag27[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at13, sizeof(cv27at13)/sizeof(int)),  &(diag27[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at14, sizeof(cv27at14)/sizeof(int)),  &(diag27[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at15, sizeof(cv27at15)/sizeof(int)),  &(diag27[15]) );
  } 
  else if( diag->id() == -28 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)),  &(diag28[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at1, sizeof(cv28at1)/sizeof(int)),  &(diag28[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at2, sizeof(cv28at2)/sizeof(int)),  &(diag28[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at3, sizeof(cv28at3)/sizeof(int)),  &(diag28[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at4, sizeof(cv28at4)/sizeof(int)),  &(diag28[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at5, sizeof(cv28at5)/sizeof(int)),  &(diag28[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at6, sizeof(cv28at6)/sizeof(int)),  &(diag28[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at7, sizeof(cv28at7)/sizeof(int)),  &(diag28[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at8, sizeof(cv28at8)/sizeof(int)),  &(diag28[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at9, sizeof(cv28at9)/sizeof(int)),  &(diag28[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at10, sizeof(cv28at10)/sizeof(int)),  &(diag28[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at11, sizeof(cv28at11)/sizeof(int)),  &(diag28[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at12, sizeof(cv28at12)/sizeof(int)),  &(diag28[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at13, sizeof(cv28at13)/sizeof(int)),  &(diag28[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at14, sizeof(cv28at14)/sizeof(int)),  &(diag28[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at15, sizeof(cv28at15)/sizeof(int)),  &(diag28[15]) );
  } 
  else if( diag->id() == -29 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int)),  &(diag29[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at1, sizeof(cv29at1)/sizeof(int)),  &(diag29[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at2, sizeof(cv29at2)/sizeof(int)),  &(diag29[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at3, sizeof(cv29at3)/sizeof(int)),  &(diag29[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at4, sizeof(cv29at4)/sizeof(int)),  &(diag29[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at5, sizeof(cv29at5)/sizeof(int)),  &(diag29[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at6, sizeof(cv29at6)/sizeof(int)),  &(diag29[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at7, sizeof(cv29at7)/sizeof(int)),  &(diag29[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at8, sizeof(cv29at8)/sizeof(int)),  &(diag29[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at9, sizeof(cv29at9)/sizeof(int)),  &(diag29[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at10, sizeof(cv29at10)/sizeof(int)),  &(diag29[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at11, sizeof(cv29at11)/sizeof(int)),  &(diag29[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at12, sizeof(cv29at12)/sizeof(int)),  &(diag29[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at13, sizeof(cv29at13)/sizeof(int)),  &(diag29[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at14, sizeof(cv29at14)/sizeof(int)),  &(diag29[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at15, sizeof(cv29at15)/sizeof(int)),  &(diag29[15]) );
  } 
  else if( diag->id() == -30 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int)),  &(diag30[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at1, sizeof(cv30at1)/sizeof(int)),  &(diag30[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at2, sizeof(cv30at2)/sizeof(int)),  &(diag30[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at3, sizeof(cv30at3)/sizeof(int)),  &(diag30[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at4, sizeof(cv30at4)/sizeof(int)),  &(diag30[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at5, sizeof(cv30at5)/sizeof(int)),  &(diag30[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at6, sizeof(cv30at6)/sizeof(int)),  &(diag30[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at7, sizeof(cv30at7)/sizeof(int)),  &(diag30[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at8, sizeof(cv30at8)/sizeof(int)),  &(diag30[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at9, sizeof(cv30at9)/sizeof(int)),  &(diag30[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at10, sizeof(cv30at10)/sizeof(int)),  &(diag30[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at11, sizeof(cv30at11)/sizeof(int)),  &(diag30[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at12, sizeof(cv30at12)/sizeof(int)),  &(diag30[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at13, sizeof(cv30at13)/sizeof(int)),  &(diag30[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at14, sizeof(cv30at14)/sizeof(int)),  &(diag30[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at15, sizeof(cv30at15)/sizeof(int)),  &(diag30[15]) );
  } 
  else if( diag->id() == -31 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int)),  &(diag31[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at1, sizeof(cv31at1)/sizeof(int)),  &(diag31[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at2, sizeof(cv31at2)/sizeof(int)),  &(diag31[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at3, sizeof(cv31at3)/sizeof(int)),  &(diag31[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at4, sizeof(cv31at4)/sizeof(int)),  &(diag31[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at5, sizeof(cv31at5)/sizeof(int)),  &(diag31[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at6, sizeof(cv31at6)/sizeof(int)),  &(diag31[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at7, sizeof(cv31at7)/sizeof(int)),  &(diag31[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at8, sizeof(cv31at8)/sizeof(int)),  &(diag31[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at9, sizeof(cv31at9)/sizeof(int)),  &(diag31[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at10, sizeof(cv31at10)/sizeof(int)),  &(diag31[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at11, sizeof(cv31at11)/sizeof(int)),  &(diag31[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at12, sizeof(cv31at12)/sizeof(int)),  &(diag31[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at13, sizeof(cv31at13)/sizeof(int)),  &(diag31[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at14, sizeof(cv31at14)/sizeof(int)),  &(diag31[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at15, sizeof(cv31at15)/sizeof(int)),  &(diag31[15]) );
  } 
  else if( diag->id() == -32 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int)),  &(diag32[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at1, sizeof(cv32at1)/sizeof(int)),  &(diag32[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at2, sizeof(cv32at2)/sizeof(int)),  &(diag32[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at3, sizeof(cv32at3)/sizeof(int)),  &(diag32[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at4, sizeof(cv32at4)/sizeof(int)),  &(diag32[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at5, sizeof(cv32at5)/sizeof(int)),  &(diag32[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at6, sizeof(cv32at6)/sizeof(int)),  &(diag32[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at7, sizeof(cv32at7)/sizeof(int)),  &(diag32[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at8, sizeof(cv32at8)/sizeof(int)),  &(diag32[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at9, sizeof(cv32at9)/sizeof(int)),  &(diag32[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at10, sizeof(cv32at10)/sizeof(int)),  &(diag32[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at11, sizeof(cv32at11)/sizeof(int)),  &(diag32[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at12, sizeof(cv32at12)/sizeof(int)),  &(diag32[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at13, sizeof(cv32at13)/sizeof(int)),  &(diag32[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at14, sizeof(cv32at14)/sizeof(int)),  &(diag32[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at15, sizeof(cv32at15)/sizeof(int)),  &(diag32[15]) );
  } 
  else if( diag->id() == -33 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int)),  &(diag33[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at1, sizeof(cv33at1)/sizeof(int)),  &(diag33[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at2, sizeof(cv33at2)/sizeof(int)),  &(diag33[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at3, sizeof(cv33at3)/sizeof(int)),  &(diag33[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at4, sizeof(cv33at4)/sizeof(int)),  &(diag33[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at5, sizeof(cv33at5)/sizeof(int)),  &(diag33[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at6, sizeof(cv33at6)/sizeof(int)),  &(diag33[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at7, sizeof(cv33at7)/sizeof(int)),  &(diag33[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at8, sizeof(cv33at8)/sizeof(int)),  &(diag33[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at9, sizeof(cv33at9)/sizeof(int)),  &(diag33[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at10, sizeof(cv33at10)/sizeof(int)),  &(diag33[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at11, sizeof(cv33at11)/sizeof(int)),  &(diag33[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at12, sizeof(cv33at12)/sizeof(int)),  &(diag33[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at13, sizeof(cv33at13)/sizeof(int)),  &(diag33[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at14, sizeof(cv33at14)/sizeof(int)),  &(diag33[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at15, sizeof(cv33at15)/sizeof(int)),  &(diag33[15]) );
  } 
  else if( diag->id() == -34 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at0, sizeof(cv34at0)/sizeof(int)),  &(diag34[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at1, sizeof(cv34at1)/sizeof(int)),  &(diag34[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at2, sizeof(cv34at2)/sizeof(int)),  &(diag34[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at3, sizeof(cv34at3)/sizeof(int)),  &(diag34[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at4, sizeof(cv34at4)/sizeof(int)),  &(diag34[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at5, sizeof(cv34at5)/sizeof(int)),  &(diag34[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at6, sizeof(cv34at6)/sizeof(int)),  &(diag34[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at7, sizeof(cv34at7)/sizeof(int)),  &(diag34[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at8, sizeof(cv34at8)/sizeof(int)),  &(diag34[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at9, sizeof(cv34at9)/sizeof(int)),  &(diag34[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at10, sizeof(cv34at10)/sizeof(int)),  &(diag34[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at11, sizeof(cv34at11)/sizeof(int)),  &(diag34[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at12, sizeof(cv34at12)/sizeof(int)),  &(diag34[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at13, sizeof(cv34at13)/sizeof(int)),  &(diag34[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at14, sizeof(cv34at14)/sizeof(int)),  &(diag34[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at15, sizeof(cv34at15)/sizeof(int)),  &(diag34[15]) );
  } 
  else if( diag->id() == -35 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at0, sizeof(cv35at0)/sizeof(int)),  &(diag35[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at1, sizeof(cv35at1)/sizeof(int)),  &(diag35[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at2, sizeof(cv35at2)/sizeof(int)),  &(diag35[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at3, sizeof(cv35at3)/sizeof(int)),  &(diag35[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at4, sizeof(cv35at4)/sizeof(int)),  &(diag35[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at5, sizeof(cv35at5)/sizeof(int)),  &(diag35[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at6, sizeof(cv35at6)/sizeof(int)),  &(diag35[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at7, sizeof(cv35at7)/sizeof(int)),  &(diag35[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at8, sizeof(cv35at8)/sizeof(int)),  &(diag35[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at9, sizeof(cv35at9)/sizeof(int)),  &(diag35[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at10, sizeof(cv35at10)/sizeof(int)),  &(diag35[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at11, sizeof(cv35at11)/sizeof(int)),  &(diag35[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at12, sizeof(cv35at12)/sizeof(int)),  &(diag35[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at13, sizeof(cv35at13)/sizeof(int)),  &(diag35[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at14, sizeof(cv35at14)/sizeof(int)),  &(diag35[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at15, sizeof(cv35at15)/sizeof(int)),  &(diag35[15]) );
  } 
  else if( diag->id() == -36 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int)),  &(diag36[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at1, sizeof(cv36at1)/sizeof(int)),  &(diag36[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at2, sizeof(cv36at2)/sizeof(int)),  &(diag36[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at3, sizeof(cv36at3)/sizeof(int)),  &(diag36[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at4, sizeof(cv36at4)/sizeof(int)),  &(diag36[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at5, sizeof(cv36at5)/sizeof(int)),  &(diag36[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at6, sizeof(cv36at6)/sizeof(int)),  &(diag36[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at7, sizeof(cv36at7)/sizeof(int)),  &(diag36[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at8, sizeof(cv36at8)/sizeof(int)),  &(diag36[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at9, sizeof(cv36at9)/sizeof(int)),  &(diag36[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at10, sizeof(cv36at10)/sizeof(int)),  &(diag36[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at11, sizeof(cv36at11)/sizeof(int)),  &(diag36[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at12, sizeof(cv36at12)/sizeof(int)),  &(diag36[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at13, sizeof(cv36at13)/sizeof(int)),  &(diag36[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at14, sizeof(cv36at14)/sizeof(int)),  &(diag36[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at15, sizeof(cv36at15)/sizeof(int)),  &(diag36[15]) );
  } 
  else if( diag->id() == -37 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int)),  &(diag37[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at1, sizeof(cv37at1)/sizeof(int)),  &(diag37[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at2, sizeof(cv37at2)/sizeof(int)),  &(diag37[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at3, sizeof(cv37at3)/sizeof(int)),  &(diag37[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at4, sizeof(cv37at4)/sizeof(int)),  &(diag37[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at5, sizeof(cv37at5)/sizeof(int)),  &(diag37[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at6, sizeof(cv37at6)/sizeof(int)),  &(diag37[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at7, sizeof(cv37at7)/sizeof(int)),  &(diag37[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at8, sizeof(cv37at8)/sizeof(int)),  &(diag37[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at9, sizeof(cv37at9)/sizeof(int)),  &(diag37[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at10, sizeof(cv37at10)/sizeof(int)),  &(diag37[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at11, sizeof(cv37at11)/sizeof(int)),  &(diag37[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at12, sizeof(cv37at12)/sizeof(int)),  &(diag37[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at13, sizeof(cv37at13)/sizeof(int)),  &(diag37[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at14, sizeof(cv37at14)/sizeof(int)),  &(diag37[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at15, sizeof(cv37at15)/sizeof(int)),  &(diag37[15]) );
  } 
  else if( diag->id() == -38 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int)),  &(diag38[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at1, sizeof(cv38at1)/sizeof(int)),  &(diag38[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at2, sizeof(cv38at2)/sizeof(int)),  &(diag38[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at3, sizeof(cv38at3)/sizeof(int)),  &(diag38[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at4, sizeof(cv38at4)/sizeof(int)),  &(diag38[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at5, sizeof(cv38at5)/sizeof(int)),  &(diag38[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at6, sizeof(cv38at6)/sizeof(int)),  &(diag38[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at7, sizeof(cv38at7)/sizeof(int)),  &(diag38[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at8, sizeof(cv38at8)/sizeof(int)),  &(diag38[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at9, sizeof(cv38at9)/sizeof(int)),  &(diag38[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at10, sizeof(cv38at10)/sizeof(int)),  &(diag38[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at11, sizeof(cv38at11)/sizeof(int)),  &(diag38[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at12, sizeof(cv38at12)/sizeof(int)),  &(diag38[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at13, sizeof(cv38at13)/sizeof(int)),  &(diag38[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at14, sizeof(cv38at14)/sizeof(int)),  &(diag38[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at15, sizeof(cv38at15)/sizeof(int)),  &(diag38[15]) );
  } 
  else if( diag->id() == -39 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int)),  &(diag39[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at1, sizeof(cv39at1)/sizeof(int)),  &(diag39[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at2, sizeof(cv39at2)/sizeof(int)),  &(diag39[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at3, sizeof(cv39at3)/sizeof(int)),  &(diag39[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at4, sizeof(cv39at4)/sizeof(int)),  &(diag39[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at5, sizeof(cv39at5)/sizeof(int)),  &(diag39[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at6, sizeof(cv39at6)/sizeof(int)),  &(diag39[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at7, sizeof(cv39at7)/sizeof(int)),  &(diag39[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at8, sizeof(cv39at8)/sizeof(int)),  &(diag39[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at9, sizeof(cv39at9)/sizeof(int)),  &(diag39[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at10, sizeof(cv39at10)/sizeof(int)),  &(diag39[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at11, sizeof(cv39at11)/sizeof(int)),  &(diag39[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at12, sizeof(cv39at12)/sizeof(int)),  &(diag39[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at13, sizeof(cv39at13)/sizeof(int)),  &(diag39[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at14, sizeof(cv39at14)/sizeof(int)),  &(diag39[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at15, sizeof(cv39at15)/sizeof(int)),  &(diag39[15]) );
  } 
  else if( diag->id() == -40 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int)),  &(diag40[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at1, sizeof(cv40at1)/sizeof(int)),  &(diag40[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at2, sizeof(cv40at2)/sizeof(int)),  &(diag40[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at3, sizeof(cv40at3)/sizeof(int)),  &(diag40[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at4, sizeof(cv40at4)/sizeof(int)),  &(diag40[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at5, sizeof(cv40at5)/sizeof(int)),  &(diag40[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at6, sizeof(cv40at6)/sizeof(int)),  &(diag40[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at7, sizeof(cv40at7)/sizeof(int)),  &(diag40[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at8, sizeof(cv40at8)/sizeof(int)),  &(diag40[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at9, sizeof(cv40at9)/sizeof(int)),  &(diag40[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at10, sizeof(cv40at10)/sizeof(int)),  &(diag40[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at11, sizeof(cv40at11)/sizeof(int)),  &(diag40[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at12, sizeof(cv40at12)/sizeof(int)),  &(diag40[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at13, sizeof(cv40at13)/sizeof(int)),  &(diag40[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at14, sizeof(cv40at14)/sizeof(int)),  &(diag40[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at15, sizeof(cv40at15)/sizeof(int)),  &(diag40[15]) );
  } 
  else if( diag->id() == -41 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int)),  &(diag41[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at1, sizeof(cv41at1)/sizeof(int)),  &(diag41[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at2, sizeof(cv41at2)/sizeof(int)),  &(diag41[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at3, sizeof(cv41at3)/sizeof(int)),  &(diag41[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at4, sizeof(cv41at4)/sizeof(int)),  &(diag41[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at5, sizeof(cv41at5)/sizeof(int)),  &(diag41[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at6, sizeof(cv41at6)/sizeof(int)),  &(diag41[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at7, sizeof(cv41at7)/sizeof(int)),  &(diag41[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at8, sizeof(cv41at8)/sizeof(int)),  &(diag41[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at9, sizeof(cv41at9)/sizeof(int)),  &(diag41[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at10, sizeof(cv41at10)/sizeof(int)),  &(diag41[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at11, sizeof(cv41at11)/sizeof(int)),  &(diag41[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at12, sizeof(cv41at12)/sizeof(int)),  &(diag41[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at13, sizeof(cv41at13)/sizeof(int)),  &(diag41[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at14, sizeof(cv41at14)/sizeof(int)),  &(diag41[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at15, sizeof(cv41at15)/sizeof(int)),  &(diag41[15]) );
  } 
  else if( diag->id() == -42 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int)),  &(diag42[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at1, sizeof(cv42at1)/sizeof(int)),  &(diag42[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at2, sizeof(cv42at2)/sizeof(int)),  &(diag42[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at3, sizeof(cv42at3)/sizeof(int)),  &(diag42[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at4, sizeof(cv42at4)/sizeof(int)),  &(diag42[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at5, sizeof(cv42at5)/sizeof(int)),  &(diag42[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at6, sizeof(cv42at6)/sizeof(int)),  &(diag42[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at7, sizeof(cv42at7)/sizeof(int)),  &(diag42[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at8, sizeof(cv42at8)/sizeof(int)),  &(diag42[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at9, sizeof(cv42at9)/sizeof(int)),  &(diag42[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at10, sizeof(cv42at10)/sizeof(int)),  &(diag42[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at11, sizeof(cv42at11)/sizeof(int)),  &(diag42[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at12, sizeof(cv42at12)/sizeof(int)),  &(diag42[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at13, sizeof(cv42at13)/sizeof(int)),  &(diag42[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at14, sizeof(cv42at14)/sizeof(int)),  &(diag42[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at15, sizeof(cv42at15)/sizeof(int)),  &(diag42[15]) );
  } 
  else if( diag->id() == -43 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at0, sizeof(cv43at0)/sizeof(int)),  &(diag43[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at1, sizeof(cv43at1)/sizeof(int)),  &(diag43[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at2, sizeof(cv43at2)/sizeof(int)),  &(diag43[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at3, sizeof(cv43at3)/sizeof(int)),  &(diag43[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at4, sizeof(cv43at4)/sizeof(int)),  &(diag43[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at5, sizeof(cv43at5)/sizeof(int)),  &(diag43[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at6, sizeof(cv43at6)/sizeof(int)),  &(diag43[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at7, sizeof(cv43at7)/sizeof(int)),  &(diag43[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at8, sizeof(cv43at8)/sizeof(int)),  &(diag43[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at9, sizeof(cv43at9)/sizeof(int)),  &(diag43[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at10, sizeof(cv43at10)/sizeof(int)),  &(diag43[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at11, sizeof(cv43at11)/sizeof(int)),  &(diag43[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at12, sizeof(cv43at12)/sizeof(int)),  &(diag43[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at13, sizeof(cv43at13)/sizeof(int)),  &(diag43[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at14, sizeof(cv43at14)/sizeof(int)),  &(diag43[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at15, sizeof(cv43at15)/sizeof(int)),  &(diag43[15]) );
  } 
  else if( diag->id() == -44 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at0, sizeof(cv44at0)/sizeof(int)),  &(diag44[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at1, sizeof(cv44at1)/sizeof(int)),  &(diag44[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at2, sizeof(cv44at2)/sizeof(int)),  &(diag44[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at3, sizeof(cv44at3)/sizeof(int)),  &(diag44[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at4, sizeof(cv44at4)/sizeof(int)),  &(diag44[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at5, sizeof(cv44at5)/sizeof(int)),  &(diag44[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at6, sizeof(cv44at6)/sizeof(int)),  &(diag44[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at7, sizeof(cv44at7)/sizeof(int)),  &(diag44[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at8, sizeof(cv44at8)/sizeof(int)),  &(diag44[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at9, sizeof(cv44at9)/sizeof(int)),  &(diag44[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at10, sizeof(cv44at10)/sizeof(int)),  &(diag44[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at11, sizeof(cv44at11)/sizeof(int)),  &(diag44[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at12, sizeof(cv44at12)/sizeof(int)),  &(diag44[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at13, sizeof(cv44at13)/sizeof(int)),  &(diag44[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at14, sizeof(cv44at14)/sizeof(int)),  &(diag44[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at15, sizeof(cv44at15)/sizeof(int)),  &(diag44[15]) );
  } 
  else if( diag->id() == -45 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at0, sizeof(cv45at0)/sizeof(int)),  &(diag45[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at1, sizeof(cv45at1)/sizeof(int)),  &(diag45[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at2, sizeof(cv45at2)/sizeof(int)),  &(diag45[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at3, sizeof(cv45at3)/sizeof(int)),  &(diag45[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at4, sizeof(cv45at4)/sizeof(int)),  &(diag45[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at5, sizeof(cv45at5)/sizeof(int)),  &(diag45[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at6, sizeof(cv45at6)/sizeof(int)),  &(diag45[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at7, sizeof(cv45at7)/sizeof(int)),  &(diag45[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at8, sizeof(cv45at8)/sizeof(int)),  &(diag45[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at9, sizeof(cv45at9)/sizeof(int)),  &(diag45[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at10, sizeof(cv45at10)/sizeof(int)),  &(diag45[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at11, sizeof(cv45at11)/sizeof(int)),  &(diag45[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at12, sizeof(cv45at12)/sizeof(int)),  &(diag45[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at13, sizeof(cv45at13)/sizeof(int)),  &(diag45[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at14, sizeof(cv45at14)/sizeof(int)),  &(diag45[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at15, sizeof(cv45at15)/sizeof(int)),  &(diag45[15]) );
  } 
  else if( diag->id() == -46 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at0, sizeof(cv46at0)/sizeof(int)),  &(diag46[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at1, sizeof(cv46at1)/sizeof(int)),  &(diag46[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at2, sizeof(cv46at2)/sizeof(int)),  &(diag46[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at3, sizeof(cv46at3)/sizeof(int)),  &(diag46[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at4, sizeof(cv46at4)/sizeof(int)),  &(diag46[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at5, sizeof(cv46at5)/sizeof(int)),  &(diag46[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at6, sizeof(cv46at6)/sizeof(int)),  &(diag46[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at7, sizeof(cv46at7)/sizeof(int)),  &(diag46[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at8, sizeof(cv46at8)/sizeof(int)),  &(diag46[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at9, sizeof(cv46at9)/sizeof(int)),  &(diag46[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at10, sizeof(cv46at10)/sizeof(int)),  &(diag46[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at11, sizeof(cv46at11)/sizeof(int)),  &(diag46[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at12, sizeof(cv46at12)/sizeof(int)),  &(diag46[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at13, sizeof(cv46at13)/sizeof(int)),  &(diag46[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at14, sizeof(cv46at14)/sizeof(int)),  &(diag46[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at15, sizeof(cv46at15)/sizeof(int)),  &(diag46[15]) );
  } 
  else if( diag->id() == -47 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int)),  &(diag47[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at1, sizeof(cv47at1)/sizeof(int)),  &(diag47[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at2, sizeof(cv47at2)/sizeof(int)),  &(diag47[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at3, sizeof(cv47at3)/sizeof(int)),  &(diag47[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at4, sizeof(cv47at4)/sizeof(int)),  &(diag47[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at5, sizeof(cv47at5)/sizeof(int)),  &(diag47[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at6, sizeof(cv47at6)/sizeof(int)),  &(diag47[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at7, sizeof(cv47at7)/sizeof(int)),  &(diag47[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at8, sizeof(cv47at8)/sizeof(int)),  &(diag47[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at9, sizeof(cv47at9)/sizeof(int)),  &(diag47[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at10, sizeof(cv47at10)/sizeof(int)),  &(diag47[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at11, sizeof(cv47at11)/sizeof(int)),  &(diag47[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at12, sizeof(cv47at12)/sizeof(int)),  &(diag47[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at13, sizeof(cv47at13)/sizeof(int)),  &(diag47[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at14, sizeof(cv47at14)/sizeof(int)),  &(diag47[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at15, sizeof(cv47at15)/sizeof(int)),  &(diag47[15]) );
  } 
  else if( diag->id() == -48 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int)),  &(diag48[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at1, sizeof(cv48at1)/sizeof(int)),  &(diag48[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at2, sizeof(cv48at2)/sizeof(int)),  &(diag48[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at3, sizeof(cv48at3)/sizeof(int)),  &(diag48[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at4, sizeof(cv48at4)/sizeof(int)),  &(diag48[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at5, sizeof(cv48at5)/sizeof(int)),  &(diag48[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at6, sizeof(cv48at6)/sizeof(int)),  &(diag48[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at7, sizeof(cv48at7)/sizeof(int)),  &(diag48[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at8, sizeof(cv48at8)/sizeof(int)),  &(diag48[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at9, sizeof(cv48at9)/sizeof(int)),  &(diag48[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at10, sizeof(cv48at10)/sizeof(int)),  &(diag48[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at11, sizeof(cv48at11)/sizeof(int)),  &(diag48[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at12, sizeof(cv48at12)/sizeof(int)),  &(diag48[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at13, sizeof(cv48at13)/sizeof(int)),  &(diag48[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at14, sizeof(cv48at14)/sizeof(int)),  &(diag48[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at15, sizeof(cv48at15)/sizeof(int)),  &(diag48[15]) );
  } 
  else if( diag->id() == -49 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int)),  &(diag49[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at1, sizeof(cv49at1)/sizeof(int)),  &(diag49[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at2, sizeof(cv49at2)/sizeof(int)),  &(diag49[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at3, sizeof(cv49at3)/sizeof(int)),  &(diag49[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at4, sizeof(cv49at4)/sizeof(int)),  &(diag49[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at5, sizeof(cv49at5)/sizeof(int)),  &(diag49[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at6, sizeof(cv49at6)/sizeof(int)),  &(diag49[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at7, sizeof(cv49at7)/sizeof(int)),  &(diag49[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at8, sizeof(cv49at8)/sizeof(int)),  &(diag49[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at9, sizeof(cv49at9)/sizeof(int)),  &(diag49[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at10, sizeof(cv49at10)/sizeof(int)),  &(diag49[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at11, sizeof(cv49at11)/sizeof(int)),  &(diag49[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at12, sizeof(cv49at12)/sizeof(int)),  &(diag49[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at13, sizeof(cv49at13)/sizeof(int)),  &(diag49[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at14, sizeof(cv49at14)/sizeof(int)),  &(diag49[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at15, sizeof(cv49at15)/sizeof(int)),  &(diag49[15]) );
  } 
  else if( diag->id() == -50 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at0, sizeof(cv50at0)/sizeof(int)),  &(diag50[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at1, sizeof(cv50at1)/sizeof(int)),  &(diag50[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at2, sizeof(cv50at2)/sizeof(int)),  &(diag50[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at3, sizeof(cv50at3)/sizeof(int)),  &(diag50[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at4, sizeof(cv50at4)/sizeof(int)),  &(diag50[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at5, sizeof(cv50at5)/sizeof(int)),  &(diag50[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at6, sizeof(cv50at6)/sizeof(int)),  &(diag50[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at7, sizeof(cv50at7)/sizeof(int)),  &(diag50[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at8, sizeof(cv50at8)/sizeof(int)),  &(diag50[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at9, sizeof(cv50at9)/sizeof(int)),  &(diag50[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at10, sizeof(cv50at10)/sizeof(int)),  &(diag50[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at11, sizeof(cv50at11)/sizeof(int)),  &(diag50[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at12, sizeof(cv50at12)/sizeof(int)),  &(diag50[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at13, sizeof(cv50at13)/sizeof(int)),  &(diag50[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at14, sizeof(cv50at14)/sizeof(int)),  &(diag50[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at15, sizeof(cv50at15)/sizeof(int)),  &(diag50[15]) );
  } 
  else if( diag->id() == -51 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at0, sizeof(cv51at0)/sizeof(int)),  &(diag51[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at1, sizeof(cv51at1)/sizeof(int)),  &(diag51[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at2, sizeof(cv51at2)/sizeof(int)),  &(diag51[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at3, sizeof(cv51at3)/sizeof(int)),  &(diag51[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at4, sizeof(cv51at4)/sizeof(int)),  &(diag51[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at5, sizeof(cv51at5)/sizeof(int)),  &(diag51[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at6, sizeof(cv51at6)/sizeof(int)),  &(diag51[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at7, sizeof(cv51at7)/sizeof(int)),  &(diag51[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at8, sizeof(cv51at8)/sizeof(int)),  &(diag51[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at9, sizeof(cv51at9)/sizeof(int)),  &(diag51[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at10, sizeof(cv51at10)/sizeof(int)),  &(diag51[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at11, sizeof(cv51at11)/sizeof(int)),  &(diag51[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at12, sizeof(cv51at12)/sizeof(int)),  &(diag51[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at13, sizeof(cv51at13)/sizeof(int)),  &(diag51[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at14, sizeof(cv51at14)/sizeof(int)),  &(diag51[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at15, sizeof(cv51at15)/sizeof(int)),  &(diag51[15]) );
  } 
  else if( diag->id() == -52 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at0, sizeof(cv52at0)/sizeof(int)),  &(diag52[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at1, sizeof(cv52at1)/sizeof(int)),  &(diag52[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at2, sizeof(cv52at2)/sizeof(int)),  &(diag52[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at3, sizeof(cv52at3)/sizeof(int)),  &(diag52[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at4, sizeof(cv52at4)/sizeof(int)),  &(diag52[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at5, sizeof(cv52at5)/sizeof(int)),  &(diag52[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at6, sizeof(cv52at6)/sizeof(int)),  &(diag52[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at7, sizeof(cv52at7)/sizeof(int)),  &(diag52[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at8, sizeof(cv52at8)/sizeof(int)),  &(diag52[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at9, sizeof(cv52at9)/sizeof(int)),  &(diag52[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at10, sizeof(cv52at10)/sizeof(int)),  &(diag52[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at11, sizeof(cv52at11)/sizeof(int)),  &(diag52[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at12, sizeof(cv52at12)/sizeof(int)),  &(diag52[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at13, sizeof(cv52at13)/sizeof(int)),  &(diag52[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at14, sizeof(cv52at14)/sizeof(int)),  &(diag52[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at15, sizeof(cv52at15)/sizeof(int)),  &(diag52[15]) );
  } 
  else if( diag->id() == -53 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int)),  &(diag53[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at1, sizeof(cv53at1)/sizeof(int)),  &(diag53[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at2, sizeof(cv53at2)/sizeof(int)),  &(diag53[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at3, sizeof(cv53at3)/sizeof(int)),  &(diag53[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at4, sizeof(cv53at4)/sizeof(int)),  &(diag53[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at5, sizeof(cv53at5)/sizeof(int)),  &(diag53[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at6, sizeof(cv53at6)/sizeof(int)),  &(diag53[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at7, sizeof(cv53at7)/sizeof(int)),  &(diag53[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at8, sizeof(cv53at8)/sizeof(int)),  &(diag53[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at9, sizeof(cv53at9)/sizeof(int)),  &(diag53[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at10, sizeof(cv53at10)/sizeof(int)),  &(diag53[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at11, sizeof(cv53at11)/sizeof(int)),  &(diag53[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at12, sizeof(cv53at12)/sizeof(int)),  &(diag53[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at13, sizeof(cv53at13)/sizeof(int)),  &(diag53[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at14, sizeof(cv53at14)/sizeof(int)),  &(diag53[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at15, sizeof(cv53at15)/sizeof(int)),  &(diag53[15]) );
  } 
  else if( diag->id() == -54 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int)),  &(diag54[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at1, sizeof(cv54at1)/sizeof(int)),  &(diag54[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at2, sizeof(cv54at2)/sizeof(int)),  &(diag54[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at3, sizeof(cv54at3)/sizeof(int)),  &(diag54[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at4, sizeof(cv54at4)/sizeof(int)),  &(diag54[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at5, sizeof(cv54at5)/sizeof(int)),  &(diag54[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at6, sizeof(cv54at6)/sizeof(int)),  &(diag54[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at7, sizeof(cv54at7)/sizeof(int)),  &(diag54[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at8, sizeof(cv54at8)/sizeof(int)),  &(diag54[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at9, sizeof(cv54at9)/sizeof(int)),  &(diag54[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at10, sizeof(cv54at10)/sizeof(int)),  &(diag54[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at11, sizeof(cv54at11)/sizeof(int)),  &(diag54[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at12, sizeof(cv54at12)/sizeof(int)),  &(diag54[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at13, sizeof(cv54at13)/sizeof(int)),  &(diag54[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at14, sizeof(cv54at14)/sizeof(int)),  &(diag54[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at15, sizeof(cv54at15)/sizeof(int)),  &(diag54[15]) );
  } 
  else if( diag->id() == -55 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int)),  &(diag55[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at1, sizeof(cv55at1)/sizeof(int)),  &(diag55[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at2, sizeof(cv55at2)/sizeof(int)),  &(diag55[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at3, sizeof(cv55at3)/sizeof(int)),  &(diag55[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at4, sizeof(cv55at4)/sizeof(int)),  &(diag55[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at5, sizeof(cv55at5)/sizeof(int)),  &(diag55[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at6, sizeof(cv55at6)/sizeof(int)),  &(diag55[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at7, sizeof(cv55at7)/sizeof(int)),  &(diag55[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at8, sizeof(cv55at8)/sizeof(int)),  &(diag55[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at9, sizeof(cv55at9)/sizeof(int)),  &(diag55[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at10, sizeof(cv55at10)/sizeof(int)),  &(diag55[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at11, sizeof(cv55at11)/sizeof(int)),  &(diag55[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at12, sizeof(cv55at12)/sizeof(int)),  &(diag55[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at13, sizeof(cv55at13)/sizeof(int)),  &(diag55[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at14, sizeof(cv55at14)/sizeof(int)),  &(diag55[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at15, sizeof(cv55at15)/sizeof(int)),  &(diag55[15]) );
  } 
  else if( diag->id() == -56 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int)),  &(diag56[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at1, sizeof(cv56at1)/sizeof(int)),  &(diag56[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at2, sizeof(cv56at2)/sizeof(int)),  &(diag56[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at3, sizeof(cv56at3)/sizeof(int)),  &(diag56[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at4, sizeof(cv56at4)/sizeof(int)),  &(diag56[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at5, sizeof(cv56at5)/sizeof(int)),  &(diag56[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at6, sizeof(cv56at6)/sizeof(int)),  &(diag56[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at7, sizeof(cv56at7)/sizeof(int)),  &(diag56[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at8, sizeof(cv56at8)/sizeof(int)),  &(diag56[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at9, sizeof(cv56at9)/sizeof(int)),  &(diag56[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at10, sizeof(cv56at10)/sizeof(int)),  &(diag56[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at11, sizeof(cv56at11)/sizeof(int)),  &(diag56[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at12, sizeof(cv56at12)/sizeof(int)),  &(diag56[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at13, sizeof(cv56at13)/sizeof(int)),  &(diag56[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at14, sizeof(cv56at14)/sizeof(int)),  &(diag56[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at15, sizeof(cv56at15)/sizeof(int)),  &(diag56[15]) );
  } 
  else if( diag->id() == -57 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int)),  &(diag57[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at1, sizeof(cv57at1)/sizeof(int)),  &(diag57[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at2, sizeof(cv57at2)/sizeof(int)),  &(diag57[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at3, sizeof(cv57at3)/sizeof(int)),  &(diag57[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at4, sizeof(cv57at4)/sizeof(int)),  &(diag57[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at5, sizeof(cv57at5)/sizeof(int)),  &(diag57[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at6, sizeof(cv57at6)/sizeof(int)),  &(diag57[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at7, sizeof(cv57at7)/sizeof(int)),  &(diag57[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at8, sizeof(cv57at8)/sizeof(int)),  &(diag57[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at9, sizeof(cv57at9)/sizeof(int)),  &(diag57[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at10, sizeof(cv57at10)/sizeof(int)),  &(diag57[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at11, sizeof(cv57at11)/sizeof(int)),  &(diag57[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at12, sizeof(cv57at12)/sizeof(int)),  &(diag57[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at13, sizeof(cv57at13)/sizeof(int)),  &(diag57[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at14, sizeof(cv57at14)/sizeof(int)),  &(diag57[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at15, sizeof(cv57at15)/sizeof(int)),  &(diag57[15]) );
  } 
  else if( diag->id() == -58 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at0, sizeof(cv58at0)/sizeof(int)),  &(diag58[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at1, sizeof(cv58at1)/sizeof(int)),  &(diag58[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at2, sizeof(cv58at2)/sizeof(int)),  &(diag58[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at3, sizeof(cv58at3)/sizeof(int)),  &(diag58[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at4, sizeof(cv58at4)/sizeof(int)),  &(diag58[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at5, sizeof(cv58at5)/sizeof(int)),  &(diag58[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at6, sizeof(cv58at6)/sizeof(int)),  &(diag58[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at7, sizeof(cv58at7)/sizeof(int)),  &(diag58[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at8, sizeof(cv58at8)/sizeof(int)),  &(diag58[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at9, sizeof(cv58at9)/sizeof(int)),  &(diag58[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at10, sizeof(cv58at10)/sizeof(int)),  &(diag58[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at11, sizeof(cv58at11)/sizeof(int)),  &(diag58[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at12, sizeof(cv58at12)/sizeof(int)),  &(diag58[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at13, sizeof(cv58at13)/sizeof(int)),  &(diag58[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at14, sizeof(cv58at14)/sizeof(int)),  &(diag58[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at15, sizeof(cv58at15)/sizeof(int)),  &(diag58[15]) );
  } 
  else if( diag->id() == -59 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at0, sizeof(cv59at0)/sizeof(int)),  &(diag59[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at1, sizeof(cv59at1)/sizeof(int)),  &(diag59[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at2, sizeof(cv59at2)/sizeof(int)),  &(diag59[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at3, sizeof(cv59at3)/sizeof(int)),  &(diag59[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at4, sizeof(cv59at4)/sizeof(int)),  &(diag59[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at5, sizeof(cv59at5)/sizeof(int)),  &(diag59[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at6, sizeof(cv59at6)/sizeof(int)),  &(diag59[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at7, sizeof(cv59at7)/sizeof(int)),  &(diag59[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at8, sizeof(cv59at8)/sizeof(int)),  &(diag59[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at9, sizeof(cv59at9)/sizeof(int)),  &(diag59[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at10, sizeof(cv59at10)/sizeof(int)),  &(diag59[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at11, sizeof(cv59at11)/sizeof(int)),  &(diag59[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at12, sizeof(cv59at12)/sizeof(int)),  &(diag59[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at13, sizeof(cv59at13)/sizeof(int)),  &(diag59[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at14, sizeof(cv59at14)/sizeof(int)),  &(diag59[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at15, sizeof(cv59at15)/sizeof(int)),  &(diag59[15]) );
  } 
  else if( diag->id() == -60 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at0, sizeof(cv60at0)/sizeof(int)),  &(diag60[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at1, sizeof(cv60at1)/sizeof(int)),  &(diag60[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at2, sizeof(cv60at2)/sizeof(int)),  &(diag60[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at3, sizeof(cv60at3)/sizeof(int)),  &(diag60[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at4, sizeof(cv60at4)/sizeof(int)),  &(diag60[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at5, sizeof(cv60at5)/sizeof(int)),  &(diag60[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at6, sizeof(cv60at6)/sizeof(int)),  &(diag60[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at7, sizeof(cv60at7)/sizeof(int)),  &(diag60[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at8, sizeof(cv60at8)/sizeof(int)),  &(diag60[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at9, sizeof(cv60at9)/sizeof(int)),  &(diag60[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at10, sizeof(cv60at10)/sizeof(int)),  &(diag60[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at11, sizeof(cv60at11)/sizeof(int)),  &(diag60[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at12, sizeof(cv60at12)/sizeof(int)),  &(diag60[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at13, sizeof(cv60at13)/sizeof(int)),  &(diag60[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at14, sizeof(cv60at14)/sizeof(int)),  &(diag60[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at15, sizeof(cv60at15)/sizeof(int)),  &(diag60[15]) );
  } 
  else if( diag->id() == -61 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at0, sizeof(cv61at0)/sizeof(int)),  &(diag61[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at1, sizeof(cv61at1)/sizeof(int)),  &(diag61[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at2, sizeof(cv61at2)/sizeof(int)),  &(diag61[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at3, sizeof(cv61at3)/sizeof(int)),  &(diag61[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at4, sizeof(cv61at4)/sizeof(int)),  &(diag61[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at5, sizeof(cv61at5)/sizeof(int)),  &(diag61[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at6, sizeof(cv61at6)/sizeof(int)),  &(diag61[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at7, sizeof(cv61at7)/sizeof(int)),  &(diag61[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at8, sizeof(cv61at8)/sizeof(int)),  &(diag61[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at9, sizeof(cv61at9)/sizeof(int)),  &(diag61[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at10, sizeof(cv61at10)/sizeof(int)),  &(diag61[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at11, sizeof(cv61at11)/sizeof(int)),  &(diag61[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at12, sizeof(cv61at12)/sizeof(int)),  &(diag61[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at13, sizeof(cv61at13)/sizeof(int)),  &(diag61[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at14, sizeof(cv61at14)/sizeof(int)),  &(diag61[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at15, sizeof(cv61at15)/sizeof(int)),  &(diag61[15]) );
  } 
  else if( diag->id() == -62 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at0, sizeof(cv62at0)/sizeof(int)),  &(diag62[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at1, sizeof(cv62at1)/sizeof(int)),  &(diag62[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at2, sizeof(cv62at2)/sizeof(int)),  &(diag62[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at3, sizeof(cv62at3)/sizeof(int)),  &(diag62[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at4, sizeof(cv62at4)/sizeof(int)),  &(diag62[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at5, sizeof(cv62at5)/sizeof(int)),  &(diag62[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at6, sizeof(cv62at6)/sizeof(int)),  &(diag62[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at7, sizeof(cv62at7)/sizeof(int)),  &(diag62[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at8, sizeof(cv62at8)/sizeof(int)),  &(diag62[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at9, sizeof(cv62at9)/sizeof(int)),  &(diag62[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at10, sizeof(cv62at10)/sizeof(int)),  &(diag62[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at11, sizeof(cv62at11)/sizeof(int)),  &(diag62[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at12, sizeof(cv62at12)/sizeof(int)),  &(diag62[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at13, sizeof(cv62at13)/sizeof(int)),  &(diag62[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at14, sizeof(cv62at14)/sizeof(int)),  &(diag62[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at15, sizeof(cv62at15)/sizeof(int)),  &(diag62[15]) );
  } 
  else if( diag->id() == -63 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at0, sizeof(cv63at0)/sizeof(int)),  &(diag63[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at1, sizeof(cv63at1)/sizeof(int)),  &(diag63[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at2, sizeof(cv63at2)/sizeof(int)),  &(diag63[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at3, sizeof(cv63at3)/sizeof(int)),  &(diag63[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at4, sizeof(cv63at4)/sizeof(int)),  &(diag63[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at5, sizeof(cv63at5)/sizeof(int)),  &(diag63[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at6, sizeof(cv63at6)/sizeof(int)),  &(diag63[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at7, sizeof(cv63at7)/sizeof(int)),  &(diag63[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at8, sizeof(cv63at8)/sizeof(int)),  &(diag63[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at9, sizeof(cv63at9)/sizeof(int)),  &(diag63[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at10, sizeof(cv63at10)/sizeof(int)),  &(diag63[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at11, sizeof(cv63at11)/sizeof(int)),  &(diag63[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at12, sizeof(cv63at12)/sizeof(int)),  &(diag63[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at13, sizeof(cv63at13)/sizeof(int)),  &(diag63[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at14, sizeof(cv63at14)/sizeof(int)),  &(diag63[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at15, sizeof(cv63at15)/sizeof(int)),  &(diag63[15]) );
  } 
  else if( diag->id() == -64 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at0, sizeof(cv64at0)/sizeof(int)),  &(diag64[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at1, sizeof(cv64at1)/sizeof(int)),  &(diag64[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at2, sizeof(cv64at2)/sizeof(int)),  &(diag64[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at3, sizeof(cv64at3)/sizeof(int)),  &(diag64[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at4, sizeof(cv64at4)/sizeof(int)),  &(diag64[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at5, sizeof(cv64at5)/sizeof(int)),  &(diag64[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at6, sizeof(cv64at6)/sizeof(int)),  &(diag64[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at7, sizeof(cv64at7)/sizeof(int)),  &(diag64[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at8, sizeof(cv64at8)/sizeof(int)),  &(diag64[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at9, sizeof(cv64at9)/sizeof(int)),  &(diag64[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at10, sizeof(cv64at10)/sizeof(int)),  &(diag64[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at11, sizeof(cv64at11)/sizeof(int)),  &(diag64[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at12, sizeof(cv64at12)/sizeof(int)),  &(diag64[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at13, sizeof(cv64at13)/sizeof(int)),  &(diag64[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at14, sizeof(cv64at14)/sizeof(int)),  &(diag64[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at15, sizeof(cv64at15)/sizeof(int)),  &(diag64[15]) );
  } 
  else if( diag->id() == -65 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at0, sizeof(cv65at0)/sizeof(int)),  &(diag65[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at1, sizeof(cv65at1)/sizeof(int)),  &(diag65[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at2, sizeof(cv65at2)/sizeof(int)),  &(diag65[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at3, sizeof(cv65at3)/sizeof(int)),  &(diag65[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at4, sizeof(cv65at4)/sizeof(int)),  &(diag65[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at5, sizeof(cv65at5)/sizeof(int)),  &(diag65[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at6, sizeof(cv65at6)/sizeof(int)),  &(diag65[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at7, sizeof(cv65at7)/sizeof(int)),  &(diag65[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at8, sizeof(cv65at8)/sizeof(int)),  &(diag65[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at9, sizeof(cv65at9)/sizeof(int)),  &(diag65[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at10, sizeof(cv65at10)/sizeof(int)),  &(diag65[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at11, sizeof(cv65at11)/sizeof(int)),  &(diag65[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at12, sizeof(cv65at12)/sizeof(int)),  &(diag65[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at13, sizeof(cv65at13)/sizeof(int)),  &(diag65[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at14, sizeof(cv65at14)/sizeof(int)),  &(diag65[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at15, sizeof(cv65at15)/sizeof(int)),  &(diag65[15]) );
  } 
  else if( diag->id() == -66 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at0, sizeof(cv66at0)/sizeof(int)),  &(diag66[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at1, sizeof(cv66at1)/sizeof(int)),  &(diag66[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at2, sizeof(cv66at2)/sizeof(int)),  &(diag66[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at3, sizeof(cv66at3)/sizeof(int)),  &(diag66[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at4, sizeof(cv66at4)/sizeof(int)),  &(diag66[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at5, sizeof(cv66at5)/sizeof(int)),  &(diag66[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at6, sizeof(cv66at6)/sizeof(int)),  &(diag66[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at7, sizeof(cv66at7)/sizeof(int)),  &(diag66[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at8, sizeof(cv66at8)/sizeof(int)),  &(diag66[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at9, sizeof(cv66at9)/sizeof(int)),  &(diag66[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at10, sizeof(cv66at10)/sizeof(int)),  &(diag66[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at11, sizeof(cv66at11)/sizeof(int)),  &(diag66[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at12, sizeof(cv66at12)/sizeof(int)),  &(diag66[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at13, sizeof(cv66at13)/sizeof(int)),  &(diag66[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at14, sizeof(cv66at14)/sizeof(int)),  &(diag66[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at15, sizeof(cv66at15)/sizeof(int)),  &(diag66[15]) );
  } 
  else if( diag->id() == -67 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at0, sizeof(cv67at0)/sizeof(int)),  &(diag67[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at1, sizeof(cv67at1)/sizeof(int)),  &(diag67[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at2, sizeof(cv67at2)/sizeof(int)),  &(diag67[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at3, sizeof(cv67at3)/sizeof(int)),  &(diag67[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at4, sizeof(cv67at4)/sizeof(int)),  &(diag67[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at5, sizeof(cv67at5)/sizeof(int)),  &(diag67[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at6, sizeof(cv67at6)/sizeof(int)),  &(diag67[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at7, sizeof(cv67at7)/sizeof(int)),  &(diag67[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at8, sizeof(cv67at8)/sizeof(int)),  &(diag67[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at9, sizeof(cv67at9)/sizeof(int)),  &(diag67[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at10, sizeof(cv67at10)/sizeof(int)),  &(diag67[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at11, sizeof(cv67at11)/sizeof(int)),  &(diag67[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at12, sizeof(cv67at12)/sizeof(int)),  &(diag67[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at13, sizeof(cv67at13)/sizeof(int)),  &(diag67[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at14, sizeof(cv67at14)/sizeof(int)),  &(diag67[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at15, sizeof(cv67at15)/sizeof(int)),  &(diag67[15]) );
  } 
  else if( diag->id() == -68 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int)),  &(diag68[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at1, sizeof(cv68at1)/sizeof(int)),  &(diag68[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at2, sizeof(cv68at2)/sizeof(int)),  &(diag68[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at3, sizeof(cv68at3)/sizeof(int)),  &(diag68[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at4, sizeof(cv68at4)/sizeof(int)),  &(diag68[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at5, sizeof(cv68at5)/sizeof(int)),  &(diag68[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at6, sizeof(cv68at6)/sizeof(int)),  &(diag68[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at7, sizeof(cv68at7)/sizeof(int)),  &(diag68[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at8, sizeof(cv68at8)/sizeof(int)),  &(diag68[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at9, sizeof(cv68at9)/sizeof(int)),  &(diag68[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at10, sizeof(cv68at10)/sizeof(int)),  &(diag68[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at11, sizeof(cv68at11)/sizeof(int)),  &(diag68[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at12, sizeof(cv68at12)/sizeof(int)),  &(diag68[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at13, sizeof(cv68at13)/sizeof(int)),  &(diag68[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at14, sizeof(cv68at14)/sizeof(int)),  &(diag68[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at15, sizeof(cv68at15)/sizeof(int)),  &(diag68[15]) );
  } 
  else if( diag->id() == -69 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int)),  &(diag69[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at1, sizeof(cv69at1)/sizeof(int)),  &(diag69[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at2, sizeof(cv69at2)/sizeof(int)),  &(diag69[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at3, sizeof(cv69at3)/sizeof(int)),  &(diag69[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at4, sizeof(cv69at4)/sizeof(int)),  &(diag69[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at5, sizeof(cv69at5)/sizeof(int)),  &(diag69[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at6, sizeof(cv69at6)/sizeof(int)),  &(diag69[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at7, sizeof(cv69at7)/sizeof(int)),  &(diag69[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at8, sizeof(cv69at8)/sizeof(int)),  &(diag69[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at9, sizeof(cv69at9)/sizeof(int)),  &(diag69[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at10, sizeof(cv69at10)/sizeof(int)),  &(diag69[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at11, sizeof(cv69at11)/sizeof(int)),  &(diag69[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at12, sizeof(cv69at12)/sizeof(int)),  &(diag69[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at13, sizeof(cv69at13)/sizeof(int)),  &(diag69[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at14, sizeof(cv69at14)/sizeof(int)),  &(diag69[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at15, sizeof(cv69at15)/sizeof(int)),  &(diag69[15]) );
  } 
  else if( diag->id() == -70 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at0, sizeof(cv70at0)/sizeof(int)),  &(diag70[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at1, sizeof(cv70at1)/sizeof(int)),  &(diag70[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at2, sizeof(cv70at2)/sizeof(int)),  &(diag70[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at3, sizeof(cv70at3)/sizeof(int)),  &(diag70[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at4, sizeof(cv70at4)/sizeof(int)),  &(diag70[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at5, sizeof(cv70at5)/sizeof(int)),  &(diag70[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at6, sizeof(cv70at6)/sizeof(int)),  &(diag70[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at7, sizeof(cv70at7)/sizeof(int)),  &(diag70[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at8, sizeof(cv70at8)/sizeof(int)),  &(diag70[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at9, sizeof(cv70at9)/sizeof(int)),  &(diag70[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at10, sizeof(cv70at10)/sizeof(int)),  &(diag70[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at11, sizeof(cv70at11)/sizeof(int)),  &(diag70[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at12, sizeof(cv70at12)/sizeof(int)),  &(diag70[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at13, sizeof(cv70at13)/sizeof(int)),  &(diag70[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at14, sizeof(cv70at14)/sizeof(int)),  &(diag70[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at15, sizeof(cv70at15)/sizeof(int)),  &(diag70[15]) );
  } 
  else if( diag->id() == -71 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at0, sizeof(cv71at0)/sizeof(int)),  &(diag71[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at1, sizeof(cv71at1)/sizeof(int)),  &(diag71[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at2, sizeof(cv71at2)/sizeof(int)),  &(diag71[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at3, sizeof(cv71at3)/sizeof(int)),  &(diag71[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at4, sizeof(cv71at4)/sizeof(int)),  &(diag71[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at5, sizeof(cv71at5)/sizeof(int)),  &(diag71[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at6, sizeof(cv71at6)/sizeof(int)),  &(diag71[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at7, sizeof(cv71at7)/sizeof(int)),  &(diag71[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at8, sizeof(cv71at8)/sizeof(int)),  &(diag71[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at9, sizeof(cv71at9)/sizeof(int)),  &(diag71[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at10, sizeof(cv71at10)/sizeof(int)),  &(diag71[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at11, sizeof(cv71at11)/sizeof(int)),  &(diag71[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at12, sizeof(cv71at12)/sizeof(int)),  &(diag71[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at13, sizeof(cv71at13)/sizeof(int)),  &(diag71[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at14, sizeof(cv71at14)/sizeof(int)),  &(diag71[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv71at15, sizeof(cv71at15)/sizeof(int)),  &(diag71[15]) );
  } 
  else if( diag->id() == -72 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at0, sizeof(cv72at0)/sizeof(int)),  &(diag72[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at1, sizeof(cv72at1)/sizeof(int)),  &(diag72[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at2, sizeof(cv72at2)/sizeof(int)),  &(diag72[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at3, sizeof(cv72at3)/sizeof(int)),  &(diag72[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at4, sizeof(cv72at4)/sizeof(int)),  &(diag72[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at5, sizeof(cv72at5)/sizeof(int)),  &(diag72[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at6, sizeof(cv72at6)/sizeof(int)),  &(diag72[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at7, sizeof(cv72at7)/sizeof(int)),  &(diag72[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at8, sizeof(cv72at8)/sizeof(int)),  &(diag72[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at9, sizeof(cv72at9)/sizeof(int)),  &(diag72[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at10, sizeof(cv72at10)/sizeof(int)),  &(diag72[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at11, sizeof(cv72at11)/sizeof(int)),  &(diag72[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at12, sizeof(cv72at12)/sizeof(int)),  &(diag72[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at13, sizeof(cv72at13)/sizeof(int)),  &(diag72[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at14, sizeof(cv72at14)/sizeof(int)),  &(diag72[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv72at15, sizeof(cv72at15)/sizeof(int)),  &(diag72[15]) );
  } 
  else if( diag->id() == -73 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at0, sizeof(cv73at0)/sizeof(int)),  &(diag73[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at1, sizeof(cv73at1)/sizeof(int)),  &(diag73[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at2, sizeof(cv73at2)/sizeof(int)),  &(diag73[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at3, sizeof(cv73at3)/sizeof(int)),  &(diag73[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at4, sizeof(cv73at4)/sizeof(int)),  &(diag73[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at5, sizeof(cv73at5)/sizeof(int)),  &(diag73[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at6, sizeof(cv73at6)/sizeof(int)),  &(diag73[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at7, sizeof(cv73at7)/sizeof(int)),  &(diag73[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at8, sizeof(cv73at8)/sizeof(int)),  &(diag73[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at9, sizeof(cv73at9)/sizeof(int)),  &(diag73[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at10, sizeof(cv73at10)/sizeof(int)),  &(diag73[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at11, sizeof(cv73at11)/sizeof(int)),  &(diag73[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at12, sizeof(cv73at12)/sizeof(int)),  &(diag73[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at13, sizeof(cv73at13)/sizeof(int)),  &(diag73[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at14, sizeof(cv73at14)/sizeof(int)),  &(diag73[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv73at15, sizeof(cv73at15)/sizeof(int)),  &(diag73[15]) );
  } 
  else if( diag->id() == -74 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at0, sizeof(cv74at0)/sizeof(int)),  &(diag74[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at1, sizeof(cv74at1)/sizeof(int)),  &(diag74[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at2, sizeof(cv74at2)/sizeof(int)),  &(diag74[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at3, sizeof(cv74at3)/sizeof(int)),  &(diag74[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at4, sizeof(cv74at4)/sizeof(int)),  &(diag74[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at5, sizeof(cv74at5)/sizeof(int)),  &(diag74[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at6, sizeof(cv74at6)/sizeof(int)),  &(diag74[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at7, sizeof(cv74at7)/sizeof(int)),  &(diag74[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at8, sizeof(cv74at8)/sizeof(int)),  &(diag74[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at9, sizeof(cv74at9)/sizeof(int)),  &(diag74[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at10, sizeof(cv74at10)/sizeof(int)),  &(diag74[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at11, sizeof(cv74at11)/sizeof(int)),  &(diag74[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at12, sizeof(cv74at12)/sizeof(int)),  &(diag74[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at13, sizeof(cv74at13)/sizeof(int)),  &(diag74[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at14, sizeof(cv74at14)/sizeof(int)),  &(diag74[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv74at15, sizeof(cv74at15)/sizeof(int)),  &(diag74[15]) );
  } 
  else if( diag->id() == -75 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at0, sizeof(cv75at0)/sizeof(int)),  &(diag75[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at1, sizeof(cv75at1)/sizeof(int)),  &(diag75[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at2, sizeof(cv75at2)/sizeof(int)),  &(diag75[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at3, sizeof(cv75at3)/sizeof(int)),  &(diag75[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at4, sizeof(cv75at4)/sizeof(int)),  &(diag75[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at5, sizeof(cv75at5)/sizeof(int)),  &(diag75[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at6, sizeof(cv75at6)/sizeof(int)),  &(diag75[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at7, sizeof(cv75at7)/sizeof(int)),  &(diag75[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at8, sizeof(cv75at8)/sizeof(int)),  &(diag75[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at9, sizeof(cv75at9)/sizeof(int)),  &(diag75[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at10, sizeof(cv75at10)/sizeof(int)),  &(diag75[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at11, sizeof(cv75at11)/sizeof(int)),  &(diag75[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at12, sizeof(cv75at12)/sizeof(int)),  &(diag75[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at13, sizeof(cv75at13)/sizeof(int)),  &(diag75[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at14, sizeof(cv75at14)/sizeof(int)),  &(diag75[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv75at15, sizeof(cv75at15)/sizeof(int)),  &(diag75[15]) );
  } 
  else if( diag->id() == -76 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at0, sizeof(cv76at0)/sizeof(int)),  &(diag76[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at1, sizeof(cv76at1)/sizeof(int)),  &(diag76[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at2, sizeof(cv76at2)/sizeof(int)),  &(diag76[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at3, sizeof(cv76at3)/sizeof(int)),  &(diag76[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at4, sizeof(cv76at4)/sizeof(int)),  &(diag76[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at5, sizeof(cv76at5)/sizeof(int)),  &(diag76[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at6, sizeof(cv76at6)/sizeof(int)),  &(diag76[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at7, sizeof(cv76at7)/sizeof(int)),  &(diag76[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at8, sizeof(cv76at8)/sizeof(int)),  &(diag76[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at9, sizeof(cv76at9)/sizeof(int)),  &(diag76[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at10, sizeof(cv76at10)/sizeof(int)),  &(diag76[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at11, sizeof(cv76at11)/sizeof(int)),  &(diag76[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at12, sizeof(cv76at12)/sizeof(int)),  &(diag76[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at13, sizeof(cv76at13)/sizeof(int)),  &(diag76[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at14, sizeof(cv76at14)/sizeof(int)),  &(diag76[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv76at15, sizeof(cv76at15)/sizeof(int)),  &(diag76[15]) );
  } 
  else if( diag->id() == -77 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at0, sizeof(cv77at0)/sizeof(int)),  &(diag77[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at1, sizeof(cv77at1)/sizeof(int)),  &(diag77[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at2, sizeof(cv77at2)/sizeof(int)),  &(diag77[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at3, sizeof(cv77at3)/sizeof(int)),  &(diag77[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at4, sizeof(cv77at4)/sizeof(int)),  &(diag77[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at5, sizeof(cv77at5)/sizeof(int)),  &(diag77[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at6, sizeof(cv77at6)/sizeof(int)),  &(diag77[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at7, sizeof(cv77at7)/sizeof(int)),  &(diag77[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at8, sizeof(cv77at8)/sizeof(int)),  &(diag77[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at9, sizeof(cv77at9)/sizeof(int)),  &(diag77[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at10, sizeof(cv77at10)/sizeof(int)),  &(diag77[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at11, sizeof(cv77at11)/sizeof(int)),  &(diag77[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at12, sizeof(cv77at12)/sizeof(int)),  &(diag77[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at13, sizeof(cv77at13)/sizeof(int)),  &(diag77[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at14, sizeof(cv77at14)/sizeof(int)),  &(diag77[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at15, sizeof(cv77at15)/sizeof(int)),  &(diag77[15]) );
  } 
  else if( diag->id() == -78 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at0, sizeof(cv78at0)/sizeof(int)),  &(diag78[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at1, sizeof(cv78at1)/sizeof(int)),  &(diag78[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at2, sizeof(cv78at2)/sizeof(int)),  &(diag78[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at3, sizeof(cv78at3)/sizeof(int)),  &(diag78[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at4, sizeof(cv78at4)/sizeof(int)),  &(diag78[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at5, sizeof(cv78at5)/sizeof(int)),  &(diag78[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at6, sizeof(cv78at6)/sizeof(int)),  &(diag78[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at7, sizeof(cv78at7)/sizeof(int)),  &(diag78[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at8, sizeof(cv78at8)/sizeof(int)),  &(diag78[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at9, sizeof(cv78at9)/sizeof(int)),  &(diag78[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at10, sizeof(cv78at10)/sizeof(int)),  &(diag78[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at11, sizeof(cv78at11)/sizeof(int)),  &(diag78[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at12, sizeof(cv78at12)/sizeof(int)),  &(diag78[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at13, sizeof(cv78at13)/sizeof(int)),  &(diag78[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at14, sizeof(cv78at14)/sizeof(int)),  &(diag78[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at15, sizeof(cv78at15)/sizeof(int)),  &(diag78[15]) );
  } 
  else if( diag->id() == -79 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at0, sizeof(cv79at0)/sizeof(int)),  &(diag79[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at1, sizeof(cv79at1)/sizeof(int)),  &(diag79[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at2, sizeof(cv79at2)/sizeof(int)),  &(diag79[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at3, sizeof(cv79at3)/sizeof(int)),  &(diag79[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at4, sizeof(cv79at4)/sizeof(int)),  &(diag79[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at5, sizeof(cv79at5)/sizeof(int)),  &(diag79[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at6, sizeof(cv79at6)/sizeof(int)),  &(diag79[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at7, sizeof(cv79at7)/sizeof(int)),  &(diag79[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at8, sizeof(cv79at8)/sizeof(int)),  &(diag79[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at9, sizeof(cv79at9)/sizeof(int)),  &(diag79[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at10, sizeof(cv79at10)/sizeof(int)),  &(diag79[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at11, sizeof(cv79at11)/sizeof(int)),  &(diag79[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at12, sizeof(cv79at12)/sizeof(int)),  &(diag79[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at13, sizeof(cv79at13)/sizeof(int)),  &(diag79[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at14, sizeof(cv79at14)/sizeof(int)),  &(diag79[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at15, sizeof(cv79at15)/sizeof(int)),  &(diag79[15]) );
  } 
  else if( diag->id() == -80 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at0, sizeof(cv80at0)/sizeof(int)),  &(diag80[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at1, sizeof(cv80at1)/sizeof(int)),  &(diag80[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at2, sizeof(cv80at2)/sizeof(int)),  &(diag80[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at3, sizeof(cv80at3)/sizeof(int)),  &(diag80[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at4, sizeof(cv80at4)/sizeof(int)),  &(diag80[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at5, sizeof(cv80at5)/sizeof(int)),  &(diag80[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at6, sizeof(cv80at6)/sizeof(int)),  &(diag80[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at7, sizeof(cv80at7)/sizeof(int)),  &(diag80[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at8, sizeof(cv80at8)/sizeof(int)),  &(diag80[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at9, sizeof(cv80at9)/sizeof(int)),  &(diag80[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at10, sizeof(cv80at10)/sizeof(int)),  &(diag80[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at11, sizeof(cv80at11)/sizeof(int)),  &(diag80[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at12, sizeof(cv80at12)/sizeof(int)),  &(diag80[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at13, sizeof(cv80at13)/sizeof(int)),  &(diag80[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at14, sizeof(cv80at14)/sizeof(int)),  &(diag80[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at15, sizeof(cv80at15)/sizeof(int)),  &(diag80[15]) );
  } 
  else if( diag->id() == -81 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at0, sizeof(cv81at0)/sizeof(int)),  &(diag81[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at1, sizeof(cv81at1)/sizeof(int)),  &(diag81[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at2, sizeof(cv81at2)/sizeof(int)),  &(diag81[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at3, sizeof(cv81at3)/sizeof(int)),  &(diag81[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at4, sizeof(cv81at4)/sizeof(int)),  &(diag81[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at5, sizeof(cv81at5)/sizeof(int)),  &(diag81[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at6, sizeof(cv81at6)/sizeof(int)),  &(diag81[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at7, sizeof(cv81at7)/sizeof(int)),  &(diag81[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at8, sizeof(cv81at8)/sizeof(int)),  &(diag81[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at9, sizeof(cv81at9)/sizeof(int)),  &(diag81[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at10, sizeof(cv81at10)/sizeof(int)),  &(diag81[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at11, sizeof(cv81at11)/sizeof(int)),  &(diag81[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at12, sizeof(cv81at12)/sizeof(int)),  &(diag81[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at13, sizeof(cv81at13)/sizeof(int)),  &(diag81[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at14, sizeof(cv81at14)/sizeof(int)),  &(diag81[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at15, sizeof(cv81at15)/sizeof(int)),  &(diag81[15]) );
  } 
  else if( diag->id() == -82 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at0, sizeof(cv82at0)/sizeof(int)),  &(diag82[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at1, sizeof(cv82at1)/sizeof(int)),  &(diag82[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at2, sizeof(cv82at2)/sizeof(int)),  &(diag82[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at3, sizeof(cv82at3)/sizeof(int)),  &(diag82[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at4, sizeof(cv82at4)/sizeof(int)),  &(diag82[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at5, sizeof(cv82at5)/sizeof(int)),  &(diag82[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at6, sizeof(cv82at6)/sizeof(int)),  &(diag82[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at7, sizeof(cv82at7)/sizeof(int)),  &(diag82[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at8, sizeof(cv82at8)/sizeof(int)),  &(diag82[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at9, sizeof(cv82at9)/sizeof(int)),  &(diag82[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at10, sizeof(cv82at10)/sizeof(int)),  &(diag82[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at11, sizeof(cv82at11)/sizeof(int)),  &(diag82[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at12, sizeof(cv82at12)/sizeof(int)),  &(diag82[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at13, sizeof(cv82at13)/sizeof(int)),  &(diag82[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at14, sizeof(cv82at14)/sizeof(int)),  &(diag82[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at15, sizeof(cv82at15)/sizeof(int)),  &(diag82[15]) );
  } 
  else if( diag->id() == -83 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at0, sizeof(cv83at0)/sizeof(int)),  &(diag83[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at1, sizeof(cv83at1)/sizeof(int)),  &(diag83[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at2, sizeof(cv83at2)/sizeof(int)),  &(diag83[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at3, sizeof(cv83at3)/sizeof(int)),  &(diag83[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at4, sizeof(cv83at4)/sizeof(int)),  &(diag83[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at5, sizeof(cv83at5)/sizeof(int)),  &(diag83[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at6, sizeof(cv83at6)/sizeof(int)),  &(diag83[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at7, sizeof(cv83at7)/sizeof(int)),  &(diag83[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at8, sizeof(cv83at8)/sizeof(int)),  &(diag83[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at9, sizeof(cv83at9)/sizeof(int)),  &(diag83[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at10, sizeof(cv83at10)/sizeof(int)),  &(diag83[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at11, sizeof(cv83at11)/sizeof(int)),  &(diag83[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at12, sizeof(cv83at12)/sizeof(int)),  &(diag83[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at13, sizeof(cv83at13)/sizeof(int)),  &(diag83[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at14, sizeof(cv83at14)/sizeof(int)),  &(diag83[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at15, sizeof(cv83at15)/sizeof(int)),  &(diag83[15]) );
  } 
  else if( diag->id() == -84 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at0, sizeof(cv84at0)/sizeof(int)),  &(diag84[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at1, sizeof(cv84at1)/sizeof(int)),  &(diag84[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at2, sizeof(cv84at2)/sizeof(int)),  &(diag84[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at3, sizeof(cv84at3)/sizeof(int)),  &(diag84[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at4, sizeof(cv84at4)/sizeof(int)),  &(diag84[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at5, sizeof(cv84at5)/sizeof(int)),  &(diag84[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at6, sizeof(cv84at6)/sizeof(int)),  &(diag84[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at7, sizeof(cv84at7)/sizeof(int)),  &(diag84[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at8, sizeof(cv84at8)/sizeof(int)),  &(diag84[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at9, sizeof(cv84at9)/sizeof(int)),  &(diag84[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at10, sizeof(cv84at10)/sizeof(int)),  &(diag84[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at11, sizeof(cv84at11)/sizeof(int)),  &(diag84[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at12, sizeof(cv84at12)/sizeof(int)),  &(diag84[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at13, sizeof(cv84at13)/sizeof(int)),  &(diag84[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at14, sizeof(cv84at14)/sizeof(int)),  &(diag84[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at15, sizeof(cv84at15)/sizeof(int)),  &(diag84[15]) );
  } 
  else if( diag->id() == -85 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at0, sizeof(cv85at0)/sizeof(int)),  &(diag85[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at1, sizeof(cv85at1)/sizeof(int)),  &(diag85[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at2, sizeof(cv85at2)/sizeof(int)),  &(diag85[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at3, sizeof(cv85at3)/sizeof(int)),  &(diag85[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at4, sizeof(cv85at4)/sizeof(int)),  &(diag85[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at5, sizeof(cv85at5)/sizeof(int)),  &(diag85[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at6, sizeof(cv85at6)/sizeof(int)),  &(diag85[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at7, sizeof(cv85at7)/sizeof(int)),  &(diag85[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at8, sizeof(cv85at8)/sizeof(int)),  &(diag85[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at9, sizeof(cv85at9)/sizeof(int)),  &(diag85[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at10, sizeof(cv85at10)/sizeof(int)),  &(diag85[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at11, sizeof(cv85at11)/sizeof(int)),  &(diag85[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at12, sizeof(cv85at12)/sizeof(int)),  &(diag85[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at13, sizeof(cv85at13)/sizeof(int)),  &(diag85[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at14, sizeof(cv85at14)/sizeof(int)),  &(diag85[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at15, sizeof(cv85at15)/sizeof(int)),  &(diag85[15]) );
  } 
  else if( diag->id() == -86 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at0, sizeof(cv86at0)/sizeof(int)),  &(diag86[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at1, sizeof(cv86at1)/sizeof(int)),  &(diag86[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at2, sizeof(cv86at2)/sizeof(int)),  &(diag86[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at3, sizeof(cv86at3)/sizeof(int)),  &(diag86[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at4, sizeof(cv86at4)/sizeof(int)),  &(diag86[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at5, sizeof(cv86at5)/sizeof(int)),  &(diag86[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at6, sizeof(cv86at6)/sizeof(int)),  &(diag86[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at7, sizeof(cv86at7)/sizeof(int)),  &(diag86[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at8, sizeof(cv86at8)/sizeof(int)),  &(diag86[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at9, sizeof(cv86at9)/sizeof(int)),  &(diag86[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at10, sizeof(cv86at10)/sizeof(int)),  &(diag86[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at11, sizeof(cv86at11)/sizeof(int)),  &(diag86[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at12, sizeof(cv86at12)/sizeof(int)),  &(diag86[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at13, sizeof(cv86at13)/sizeof(int)),  &(diag86[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at14, sizeof(cv86at14)/sizeof(int)),  &(diag86[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at15, sizeof(cv86at15)/sizeof(int)),  &(diag86[15]) );
  } 
  else if( diag->id() == -87 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at0, sizeof(cv87at0)/sizeof(int)),  &(diag87[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at1, sizeof(cv87at1)/sizeof(int)),  &(diag87[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at2, sizeof(cv87at2)/sizeof(int)),  &(diag87[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at3, sizeof(cv87at3)/sizeof(int)),  &(diag87[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at4, sizeof(cv87at4)/sizeof(int)),  &(diag87[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at5, sizeof(cv87at5)/sizeof(int)),  &(diag87[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at6, sizeof(cv87at6)/sizeof(int)),  &(diag87[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at7, sizeof(cv87at7)/sizeof(int)),  &(diag87[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at8, sizeof(cv87at8)/sizeof(int)),  &(diag87[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at9, sizeof(cv87at9)/sizeof(int)),  &(diag87[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at10, sizeof(cv87at10)/sizeof(int)),  &(diag87[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at11, sizeof(cv87at11)/sizeof(int)),  &(diag87[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at12, sizeof(cv87at12)/sizeof(int)),  &(diag87[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at13, sizeof(cv87at13)/sizeof(int)),  &(diag87[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at14, sizeof(cv87at14)/sizeof(int)),  &(diag87[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at15, sizeof(cv87at15)/sizeof(int)),  &(diag87[15]) );
  } 
  else if( diag->id() == -88 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at0, sizeof(cv88at0)/sizeof(int)),  &(diag88[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at1, sizeof(cv88at1)/sizeof(int)),  &(diag88[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at2, sizeof(cv88at2)/sizeof(int)),  &(diag88[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at3, sizeof(cv88at3)/sizeof(int)),  &(diag88[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at4, sizeof(cv88at4)/sizeof(int)),  &(diag88[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at5, sizeof(cv88at5)/sizeof(int)),  &(diag88[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at6, sizeof(cv88at6)/sizeof(int)),  &(diag88[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at7, sizeof(cv88at7)/sizeof(int)),  &(diag88[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at8, sizeof(cv88at8)/sizeof(int)),  &(diag88[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at9, sizeof(cv88at9)/sizeof(int)),  &(diag88[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at10, sizeof(cv88at10)/sizeof(int)),  &(diag88[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at11, sizeof(cv88at11)/sizeof(int)),  &(diag88[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at12, sizeof(cv88at12)/sizeof(int)),  &(diag88[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at13, sizeof(cv88at13)/sizeof(int)),  &(diag88[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at14, sizeof(cv88at14)/sizeof(int)),  &(diag88[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv88at15, sizeof(cv88at15)/sizeof(int)),  &(diag88[15]) );
  } 
  else if( diag->id() == -89 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at0, sizeof(cv89at0)/sizeof(int)),  &(diag89[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at1, sizeof(cv89at1)/sizeof(int)),  &(diag89[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at2, sizeof(cv89at2)/sizeof(int)),  &(diag89[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at3, sizeof(cv89at3)/sizeof(int)),  &(diag89[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at4, sizeof(cv89at4)/sizeof(int)),  &(diag89[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at5, sizeof(cv89at5)/sizeof(int)),  &(diag89[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at6, sizeof(cv89at6)/sizeof(int)),  &(diag89[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at7, sizeof(cv89at7)/sizeof(int)),  &(diag89[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at8, sizeof(cv89at8)/sizeof(int)),  &(diag89[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at9, sizeof(cv89at9)/sizeof(int)),  &(diag89[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at10, sizeof(cv89at10)/sizeof(int)),  &(diag89[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at11, sizeof(cv89at11)/sizeof(int)),  &(diag89[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at12, sizeof(cv89at12)/sizeof(int)),  &(diag89[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at13, sizeof(cv89at13)/sizeof(int)),  &(diag89[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at14, sizeof(cv89at14)/sizeof(int)),  &(diag89[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv89at15, sizeof(cv89at15)/sizeof(int)),  &(diag89[15]) );
  } 
  else if( diag->id() == -90 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at0, sizeof(cv90at0)/sizeof(int)),  &(diag90[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at1, sizeof(cv90at1)/sizeof(int)),  &(diag90[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at2, sizeof(cv90at2)/sizeof(int)),  &(diag90[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at3, sizeof(cv90at3)/sizeof(int)),  &(diag90[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at4, sizeof(cv90at4)/sizeof(int)),  &(diag90[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at5, sizeof(cv90at5)/sizeof(int)),  &(diag90[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at6, sizeof(cv90at6)/sizeof(int)),  &(diag90[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at7, sizeof(cv90at7)/sizeof(int)),  &(diag90[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at8, sizeof(cv90at8)/sizeof(int)),  &(diag90[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at9, sizeof(cv90at9)/sizeof(int)),  &(diag90[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at10, sizeof(cv90at10)/sizeof(int)),  &(diag90[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at11, sizeof(cv90at11)/sizeof(int)),  &(diag90[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at12, sizeof(cv90at12)/sizeof(int)),  &(diag90[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at13, sizeof(cv90at13)/sizeof(int)),  &(diag90[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at14, sizeof(cv90at14)/sizeof(int)),  &(diag90[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv90at15, sizeof(cv90at15)/sizeof(int)),  &(diag90[15]) );
  } 
  else if( diag->id() == -91 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at0, sizeof(cv91at0)/sizeof(int)),  &(diag91[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at1, sizeof(cv91at1)/sizeof(int)),  &(diag91[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at2, sizeof(cv91at2)/sizeof(int)),  &(diag91[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at3, sizeof(cv91at3)/sizeof(int)),  &(diag91[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at4, sizeof(cv91at4)/sizeof(int)),  &(diag91[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at5, sizeof(cv91at5)/sizeof(int)),  &(diag91[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at6, sizeof(cv91at6)/sizeof(int)),  &(diag91[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at7, sizeof(cv91at7)/sizeof(int)),  &(diag91[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at8, sizeof(cv91at8)/sizeof(int)),  &(diag91[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at9, sizeof(cv91at9)/sizeof(int)),  &(diag91[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at10, sizeof(cv91at10)/sizeof(int)),  &(diag91[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at11, sizeof(cv91at11)/sizeof(int)),  &(diag91[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at12, sizeof(cv91at12)/sizeof(int)),  &(diag91[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at13, sizeof(cv91at13)/sizeof(int)),  &(diag91[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at14, sizeof(cv91at14)/sizeof(int)),  &(diag91[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv91at15, sizeof(cv91at15)/sizeof(int)),  &(diag91[15]) );
  } 
  else if( diag->id() == -92 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at0, sizeof(cv92at0)/sizeof(int)),  &(diag92[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at1, sizeof(cv92at1)/sizeof(int)),  &(diag92[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at2, sizeof(cv92at2)/sizeof(int)),  &(diag92[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at3, sizeof(cv92at3)/sizeof(int)),  &(diag92[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at4, sizeof(cv92at4)/sizeof(int)),  &(diag92[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at5, sizeof(cv92at5)/sizeof(int)),  &(diag92[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at6, sizeof(cv92at6)/sizeof(int)),  &(diag92[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at7, sizeof(cv92at7)/sizeof(int)),  &(diag92[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at8, sizeof(cv92at8)/sizeof(int)),  &(diag92[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at9, sizeof(cv92at9)/sizeof(int)),  &(diag92[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at10, sizeof(cv92at10)/sizeof(int)),  &(diag92[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at11, sizeof(cv92at11)/sizeof(int)),  &(diag92[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at12, sizeof(cv92at12)/sizeof(int)),  &(diag92[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at13, sizeof(cv92at13)/sizeof(int)),  &(diag92[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at14, sizeof(cv92at14)/sizeof(int)),  &(diag92[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv92at15, sizeof(cv92at15)/sizeof(int)),  &(diag92[15]) );
  } 
  else if( diag->id() == -93 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at0, sizeof(cv93at0)/sizeof(int)),  &(diag93[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at1, sizeof(cv93at1)/sizeof(int)),  &(diag93[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at2, sizeof(cv93at2)/sizeof(int)),  &(diag93[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at3, sizeof(cv93at3)/sizeof(int)),  &(diag93[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at4, sizeof(cv93at4)/sizeof(int)),  &(diag93[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at5, sizeof(cv93at5)/sizeof(int)),  &(diag93[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at6, sizeof(cv93at6)/sizeof(int)),  &(diag93[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at7, sizeof(cv93at7)/sizeof(int)),  &(diag93[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at8, sizeof(cv93at8)/sizeof(int)),  &(diag93[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at9, sizeof(cv93at9)/sizeof(int)),  &(diag93[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at10, sizeof(cv93at10)/sizeof(int)),  &(diag93[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at11, sizeof(cv93at11)/sizeof(int)),  &(diag93[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at12, sizeof(cv93at12)/sizeof(int)),  &(diag93[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at13, sizeof(cv93at13)/sizeof(int)),  &(diag93[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at14, sizeof(cv93at14)/sizeof(int)),  &(diag93[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv93at15, sizeof(cv93at15)/sizeof(int)),  &(diag93[15]) );
  } 
  else if( diag->id() == -94 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at0, sizeof(cv94at0)/sizeof(int)),  &(diag94[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at1, sizeof(cv94at1)/sizeof(int)),  &(diag94[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at2, sizeof(cv94at2)/sizeof(int)),  &(diag94[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at3, sizeof(cv94at3)/sizeof(int)),  &(diag94[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at4, sizeof(cv94at4)/sizeof(int)),  &(diag94[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at5, sizeof(cv94at5)/sizeof(int)),  &(diag94[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at6, sizeof(cv94at6)/sizeof(int)),  &(diag94[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at7, sizeof(cv94at7)/sizeof(int)),  &(diag94[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at8, sizeof(cv94at8)/sizeof(int)),  &(diag94[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at9, sizeof(cv94at9)/sizeof(int)),  &(diag94[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at10, sizeof(cv94at10)/sizeof(int)),  &(diag94[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at11, sizeof(cv94at11)/sizeof(int)),  &(diag94[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at12, sizeof(cv94at12)/sizeof(int)),  &(diag94[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at13, sizeof(cv94at13)/sizeof(int)),  &(diag94[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at14, sizeof(cv94at14)/sizeof(int)),  &(diag94[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at15, sizeof(cv94at15)/sizeof(int)),  &(diag94[15]) );
  } 
  else if( diag->id() == -95 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at0, sizeof(cv95at0)/sizeof(int)),  &(diag95[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at1, sizeof(cv95at1)/sizeof(int)),  &(diag95[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at2, sizeof(cv95at2)/sizeof(int)),  &(diag95[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at3, sizeof(cv95at3)/sizeof(int)),  &(diag95[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at4, sizeof(cv95at4)/sizeof(int)),  &(diag95[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at5, sizeof(cv95at5)/sizeof(int)),  &(diag95[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at6, sizeof(cv95at6)/sizeof(int)),  &(diag95[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at7, sizeof(cv95at7)/sizeof(int)),  &(diag95[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at8, sizeof(cv95at8)/sizeof(int)),  &(diag95[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at9, sizeof(cv95at9)/sizeof(int)),  &(diag95[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at10, sizeof(cv95at10)/sizeof(int)),  &(diag95[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at11, sizeof(cv95at11)/sizeof(int)),  &(diag95[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at12, sizeof(cv95at12)/sizeof(int)),  &(diag95[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at13, sizeof(cv95at13)/sizeof(int)),  &(diag95[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at14, sizeof(cv95at14)/sizeof(int)),  &(diag95[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv95at15, sizeof(cv95at15)/sizeof(int)),  &(diag95[15]) );
  } 
  else if( diag->id() == -96 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at0, sizeof(cv96at0)/sizeof(int)),  &(diag96[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at1, sizeof(cv96at1)/sizeof(int)),  &(diag96[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at2, sizeof(cv96at2)/sizeof(int)),  &(diag96[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at3, sizeof(cv96at3)/sizeof(int)),  &(diag96[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at4, sizeof(cv96at4)/sizeof(int)),  &(diag96[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at5, sizeof(cv96at5)/sizeof(int)),  &(diag96[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at6, sizeof(cv96at6)/sizeof(int)),  &(diag96[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at7, sizeof(cv96at7)/sizeof(int)),  &(diag96[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at8, sizeof(cv96at8)/sizeof(int)),  &(diag96[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at9, sizeof(cv96at9)/sizeof(int)),  &(diag96[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at10, sizeof(cv96at10)/sizeof(int)),  &(diag96[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at11, sizeof(cv96at11)/sizeof(int)),  &(diag96[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at12, sizeof(cv96at12)/sizeof(int)),  &(diag96[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at13, sizeof(cv96at13)/sizeof(int)),  &(diag96[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at14, sizeof(cv96at14)/sizeof(int)),  &(diag96[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv96at15, sizeof(cv96at15)/sizeof(int)),  &(diag96[15]) );
  } 
  else if( diag->id() == -97 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at0, sizeof(cv97at0)/sizeof(int)),  &(diag97[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at1, sizeof(cv97at1)/sizeof(int)),  &(diag97[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at2, sizeof(cv97at2)/sizeof(int)),  &(diag97[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at3, sizeof(cv97at3)/sizeof(int)),  &(diag97[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at4, sizeof(cv97at4)/sizeof(int)),  &(diag97[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at5, sizeof(cv97at5)/sizeof(int)),  &(diag97[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at6, sizeof(cv97at6)/sizeof(int)),  &(diag97[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at7, sizeof(cv97at7)/sizeof(int)),  &(diag97[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at8, sizeof(cv97at8)/sizeof(int)),  &(diag97[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at9, sizeof(cv97at9)/sizeof(int)),  &(diag97[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at10, sizeof(cv97at10)/sizeof(int)),  &(diag97[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at11, sizeof(cv97at11)/sizeof(int)),  &(diag97[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at12, sizeof(cv97at12)/sizeof(int)),  &(diag97[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at13, sizeof(cv97at13)/sizeof(int)),  &(diag97[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at14, sizeof(cv97at14)/sizeof(int)),  &(diag97[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv97at15, sizeof(cv97at15)/sizeof(int)),  &(diag97[15]) );
  } 
  else if( diag->id() == -98 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at0, sizeof(cv98at0)/sizeof(int)),  &(diag98[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at1, sizeof(cv98at1)/sizeof(int)),  &(diag98[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at2, sizeof(cv98at2)/sizeof(int)),  &(diag98[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at3, sizeof(cv98at3)/sizeof(int)),  &(diag98[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at4, sizeof(cv98at4)/sizeof(int)),  &(diag98[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at5, sizeof(cv98at5)/sizeof(int)),  &(diag98[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at6, sizeof(cv98at6)/sizeof(int)),  &(diag98[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at7, sizeof(cv98at7)/sizeof(int)),  &(diag98[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at8, sizeof(cv98at8)/sizeof(int)),  &(diag98[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at9, sizeof(cv98at9)/sizeof(int)),  &(diag98[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at10, sizeof(cv98at10)/sizeof(int)),  &(diag98[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at11, sizeof(cv98at11)/sizeof(int)),  &(diag98[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at12, sizeof(cv98at12)/sizeof(int)),  &(diag98[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at13, sizeof(cv98at13)/sizeof(int)),  &(diag98[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at14, sizeof(cv98at14)/sizeof(int)),  &(diag98[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv98at15, sizeof(cv98at15)/sizeof(int)),  &(diag98[15]) );
  } 
  else if( diag->id() == -99 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at0, sizeof(cv99at0)/sizeof(int)),  &(diag99[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at1, sizeof(cv99at1)/sizeof(int)),  &(diag99[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at2, sizeof(cv99at2)/sizeof(int)),  &(diag99[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at3, sizeof(cv99at3)/sizeof(int)),  &(diag99[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at4, sizeof(cv99at4)/sizeof(int)),  &(diag99[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at5, sizeof(cv99at5)/sizeof(int)),  &(diag99[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at6, sizeof(cv99at6)/sizeof(int)),  &(diag99[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at7, sizeof(cv99at7)/sizeof(int)),  &(diag99[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at8, sizeof(cv99at8)/sizeof(int)),  &(diag99[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at9, sizeof(cv99at9)/sizeof(int)),  &(diag99[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at10, sizeof(cv99at10)/sizeof(int)),  &(diag99[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at11, sizeof(cv99at11)/sizeof(int)),  &(diag99[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at12, sizeof(cv99at12)/sizeof(int)),  &(diag99[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at13, sizeof(cv99at13)/sizeof(int)),  &(diag99[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at14, sizeof(cv99at14)/sizeof(int)),  &(diag99[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv99at15, sizeof(cv99at15)/sizeof(int)),  &(diag99[15]) );
  } 
  else if( diag->id() == -100 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at0, sizeof(cv100at0)/sizeof(int)),  &(diag100[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at1, sizeof(cv100at1)/sizeof(int)),  &(diag100[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at2, sizeof(cv100at2)/sizeof(int)),  &(diag100[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at3, sizeof(cv100at3)/sizeof(int)),  &(diag100[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at4, sizeof(cv100at4)/sizeof(int)),  &(diag100[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at5, sizeof(cv100at5)/sizeof(int)),  &(diag100[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at6, sizeof(cv100at6)/sizeof(int)),  &(diag100[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at7, sizeof(cv100at7)/sizeof(int)),  &(diag100[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at8, sizeof(cv100at8)/sizeof(int)),  &(diag100[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at9, sizeof(cv100at9)/sizeof(int)),  &(diag100[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at10, sizeof(cv100at10)/sizeof(int)),  &(diag100[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at11, sizeof(cv100at11)/sizeof(int)),  &(diag100[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at12, sizeof(cv100at12)/sizeof(int)),  &(diag100[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at13, sizeof(cv100at13)/sizeof(int)),  &(diag100[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at14, sizeof(cv100at14)/sizeof(int)),  &(diag100[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at15, sizeof(cv100at15)/sizeof(int)),  &(diag100[15]) );
  } 
  else if( diag->id() == -101 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at0, sizeof(cv101at0)/sizeof(int)),  &(diag101[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at1, sizeof(cv101at1)/sizeof(int)),  &(diag101[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at2, sizeof(cv101at2)/sizeof(int)),  &(diag101[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at3, sizeof(cv101at3)/sizeof(int)),  &(diag101[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at4, sizeof(cv101at4)/sizeof(int)),  &(diag101[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at5, sizeof(cv101at5)/sizeof(int)),  &(diag101[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at6, sizeof(cv101at6)/sizeof(int)),  &(diag101[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at7, sizeof(cv101at7)/sizeof(int)),  &(diag101[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at8, sizeof(cv101at8)/sizeof(int)),  &(diag101[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at9, sizeof(cv101at9)/sizeof(int)),  &(diag101[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at10, sizeof(cv101at10)/sizeof(int)),  &(diag101[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at11, sizeof(cv101at11)/sizeof(int)),  &(diag101[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at12, sizeof(cv101at12)/sizeof(int)),  &(diag101[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at13, sizeof(cv101at13)/sizeof(int)),  &(diag101[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at14, sizeof(cv101at14)/sizeof(int)),  &(diag101[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv101at15, sizeof(cv101at15)/sizeof(int)),  &(diag101[15]) );
  } 
  else if( diag->id() == -102 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at0, sizeof(cv102at0)/sizeof(int)),  &(diag102[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at1, sizeof(cv102at1)/sizeof(int)),  &(diag102[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at2, sizeof(cv102at2)/sizeof(int)),  &(diag102[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at3, sizeof(cv102at3)/sizeof(int)),  &(diag102[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at4, sizeof(cv102at4)/sizeof(int)),  &(diag102[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at5, sizeof(cv102at5)/sizeof(int)),  &(diag102[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at6, sizeof(cv102at6)/sizeof(int)),  &(diag102[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at7, sizeof(cv102at7)/sizeof(int)),  &(diag102[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at8, sizeof(cv102at8)/sizeof(int)),  &(diag102[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at9, sizeof(cv102at9)/sizeof(int)),  &(diag102[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at10, sizeof(cv102at10)/sizeof(int)),  &(diag102[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at11, sizeof(cv102at11)/sizeof(int)),  &(diag102[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at12, sizeof(cv102at12)/sizeof(int)),  &(diag102[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at13, sizeof(cv102at13)/sizeof(int)),  &(diag102[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at14, sizeof(cv102at14)/sizeof(int)),  &(diag102[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at15, sizeof(cv102at15)/sizeof(int)),  &(diag102[15]) );
  } 
  else if( diag->id() == -103 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at0, sizeof(cv103at0)/sizeof(int)),  &(diag103[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at1, sizeof(cv103at1)/sizeof(int)),  &(diag103[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at2, sizeof(cv103at2)/sizeof(int)),  &(diag103[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at3, sizeof(cv103at3)/sizeof(int)),  &(diag103[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at4, sizeof(cv103at4)/sizeof(int)),  &(diag103[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at5, sizeof(cv103at5)/sizeof(int)),  &(diag103[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at6, sizeof(cv103at6)/sizeof(int)),  &(diag103[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at7, sizeof(cv103at7)/sizeof(int)),  &(diag103[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at8, sizeof(cv103at8)/sizeof(int)),  &(diag103[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at9, sizeof(cv103at9)/sizeof(int)),  &(diag103[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at10, sizeof(cv103at10)/sizeof(int)),  &(diag103[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at11, sizeof(cv103at11)/sizeof(int)),  &(diag103[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at12, sizeof(cv103at12)/sizeof(int)),  &(diag103[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at13, sizeof(cv103at13)/sizeof(int)),  &(diag103[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at14, sizeof(cv103at14)/sizeof(int)),  &(diag103[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at15, sizeof(cv103at15)/sizeof(int)),  &(diag103[15]) );
  } 
  else if( diag->id() == -104 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at0, sizeof(cv104at0)/sizeof(int)),  &(diag104[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at1, sizeof(cv104at1)/sizeof(int)),  &(diag104[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at2, sizeof(cv104at2)/sizeof(int)),  &(diag104[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at3, sizeof(cv104at3)/sizeof(int)),  &(diag104[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at4, sizeof(cv104at4)/sizeof(int)),  &(diag104[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at5, sizeof(cv104at5)/sizeof(int)),  &(diag104[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at6, sizeof(cv104at6)/sizeof(int)),  &(diag104[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at7, sizeof(cv104at7)/sizeof(int)),  &(diag104[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at8, sizeof(cv104at8)/sizeof(int)),  &(diag104[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at9, sizeof(cv104at9)/sizeof(int)),  &(diag104[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at10, sizeof(cv104at10)/sizeof(int)),  &(diag104[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at11, sizeof(cv104at11)/sizeof(int)),  &(diag104[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at12, sizeof(cv104at12)/sizeof(int)),  &(diag104[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at13, sizeof(cv104at13)/sizeof(int)),  &(diag104[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at14, sizeof(cv104at14)/sizeof(int)),  &(diag104[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at15, sizeof(cv104at15)/sizeof(int)),  &(diag104[15]) );
  } 
  else if( diag->id() == -105 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at0, sizeof(cv105at0)/sizeof(int)),  &(diag105[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at1, sizeof(cv105at1)/sizeof(int)),  &(diag105[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at2, sizeof(cv105at2)/sizeof(int)),  &(diag105[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at3, sizeof(cv105at3)/sizeof(int)),  &(diag105[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at4, sizeof(cv105at4)/sizeof(int)),  &(diag105[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at5, sizeof(cv105at5)/sizeof(int)),  &(diag105[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at6, sizeof(cv105at6)/sizeof(int)),  &(diag105[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at7, sizeof(cv105at7)/sizeof(int)),  &(diag105[7]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at8, sizeof(cv105at8)/sizeof(int)),  &(diag105[8]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at9, sizeof(cv105at9)/sizeof(int)),  &(diag105[9]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at10, sizeof(cv105at10)/sizeof(int)),  &(diag105[10]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at11, sizeof(cv105at11)/sizeof(int)),  &(diag105[11]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at12, sizeof(cv105at12)/sizeof(int)),  &(diag105[12]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at13, sizeof(cv105at13)/sizeof(int)),  &(diag105[13]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at14, sizeof(cv105at14)/sizeof(int)),  &(diag105[14]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at15, sizeof(cv105at15)/sizeof(int)),  &(diag105[15]) );
  } 
  return sel;
}
