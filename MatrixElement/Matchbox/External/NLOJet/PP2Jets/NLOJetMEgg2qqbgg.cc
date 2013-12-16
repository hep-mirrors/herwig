// -*- C++ -*-
//
// NLOJetMEgg2qqbgg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEgg2qqbgg class.
//

#include "NLOJetMEgg2qqbgg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEgg2qqbgg::NLOJetMEgg2qqbgg() {}

NLOJetMEgg2qqbgg::~NLOJetMEgg2qqbgg() {}

IBPtr NLOJetMEgg2qqbgg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEgg2qqbgg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEgg2qqbgg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEgg2qqbgg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEgg2qqbgg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEgg2qqbgg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEgg2qqbgg("Herwig::NLOJetMEgg2qqbgg", "HwMatchboxNLOJet.so");

void NLOJetMEgg2qqbgg::Init() {

  static ClassDocumentation<NLOJetMEgg2qqbgg> documentation
    ("NLOJetMEgg2qqbgg");

}


void NLOJetMEgg2qqbgg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr q = quark[xi1]; 
    tcPDPtr qb = antiquark[xi1]; 
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 3, g, 4, q, 4, qb, 5, g, 5, g, -1)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, q, 3, qb, 4, q, 5, qb, 4, g, 5, g, -2)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, q, 3, qb, 4, q, 5, qb, 5, g, 4, g, -3)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, qb, g, 2, g, 1, q, 3, qb, 5, g, 5, g, -4)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, g, g, 2, qb, 1, q, 5, qb, 3, g, 5, g, -5)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, g, g, 2, qb, 1, q, 5, qb, 5, g, 3, g, -6)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, q, g, 2, g, 3, q, 1, qb, 5, g, 5, g, -7)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, g, g, 2, q, 5, q, 1, qb, 3, g, 5, g, -8)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, g, g, 2, q, 5, q, 1, qb, 5, g, 3, g, -9)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, q, g, 2, qb, 3, q, 5, qb, 1, g, 5, g, -10)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, qb, g, 2, q, 5, q, 3, qb, 1, g, 5, g, -11)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, q, 5, qb, 1, g, 3, g, -12)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, q, g, 2, qb, 3, q, 5, qb, 5, g, 1, g, -13)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, qb, g, 2, q, 5, q, 3, qb, 5, g, 1, g, -14)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 2, g, 5, q, 5, qb, 3, g, 1, g, -15)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, q, 5, qb, 3, g, 4, g, -16)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, g, 5, q, 5, qb, 4, g, 3, g, -17)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, q, 4, q, 5, q, 3, qb, 5, g, 4, g, -18)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, q, 5, q, 4, qb, 5, g, 3, g, -19)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, q, 4, q, 5, q, 3, qb, 4, g, 5, g, -20)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, q, 5, q, 4, qb, 3, g, 5, g, -21)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, qb, 4, qb, 3, q, 5, qb, 5, g, 4, g, -22)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, qb, 4, q, 5, qb, 5, g, 3, g, -23)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, qb, 4, qb, 3, q, 5, qb, 4, g, 5, g, -24)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 4, qb, 4, q, 5, qb, 3, g, 5, g, -25)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, qb, 4, g, 3, q, 4, qb, 5, g, 5, g, -26)));
    addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, q, 4, g, 4, q, 3, qb, 5, g, 5, g, -27)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, qb, qb, qb, g, 1, q, 4, qb, 2, g, 3, g, -28)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, qb, qb, qb, g, 1, q, 4, qb, 3, g, 2, g, -29)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, qb, g, g, g, 1, q, 2, qb, 4, g, 3, g, -30)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, qb, qb, g, g, 1, q, 3, qb, 4, g, 2, g, -31)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, qb, g, g, g, 1, q, 2, qb, 3, g, 4, g, -32)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, qb, qb, g, g, 1, q, 3, qb, 2, g, 4, g, -33)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 2, qb, 4, qb, 1, q, 5, qb, 5, g, 4, g, -34)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, qb, g, 3, qb, 1, q, 5, qb, 5, g, 2, g, -35)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 2, qb, 4, qb, 1, q, 5, qb, 4, g, 5, g, -36)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, qb, g, 3, qb, 1, q, 5, qb, 2, g, 5, g, -37)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 2, qb, 4, g, 1, q, 4, qb, 5, g, 5, g, -38)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, g, g, 3, g, 1, q, 2, qb, 5, g, 5, g, -39)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, q, q, q, g, 4, q, 1, qb, 2, g, 3, g, -40)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, q, q, q, g, 4, q, 1, qb, 3, g, 2, g, -41)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, q, g, g, g, 2, q, 1, qb, 4, g, 3, g, -42)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, q, q, g, g, 3, q, 1, qb, 4, g, 2, g, -43)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, q, g, g, g, 2, q, 1, qb, 3, g, 4, g, -44)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, q, q, g, g, 3, q, 1, qb, 2, g, 4, g, -45)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 2, q, 4, q, 5, q, 1, qb, 5, g, 4, g, -46)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, q, g, 3, q, 5, q, 1, qb, 5, g, 2, g, -47)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 2, q, 4, q, 5, q, 1, qb, 4, g, 5, g, -48)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, q, g, 3, q, 5, q, 1, qb, 2, g, 5, g, -49)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 2, q, 4, g, 4, q, 1, qb, 5, g, 5, g, -50)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, g, g, 3, g, 2, q, 1, qb, 5, g, 5, g, -51)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, q, q, g, 4, q, 2, qb, 1, g, 3, g, -52)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, q, g, 4, q, 3, qb, 1, g, 2, g, -53)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, qb, qb, g, 2, q, 4, qb, 1, g, 3, g, -54)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, qb, g, 3, q, 4, qb, 1, g, 2, g, -55)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, qb, g, g, 2, q, 3, qb, 1, g, 4, g, -56)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, q, g, g, 3, q, 2, qb, 1, g, 4, g, -57)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, q, 5, qb, 1, g, 4, g, -58)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, q, 5, qb, 1, g, 2, g, -59)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, q, 5, q, 4, qb, 1, g, 5, g, -60)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, q, g, 3, q, 5, q, 2, qb, 1, g, 5, g, -61)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, qb, 4, q, 5, qb, 1, g, 5, g, -62)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, qb, g, 3, qb, 2, q, 5, qb, 1, g, 5, g, -63)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, q, q, g, 4, q, 2, qb, 3, g, 1, g, -64)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, q, g, 4, q, 3, qb, 2, g, 1, g, -65)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, qb, qb, g, 2, q, 4, qb, 3, g, 1, g, -66)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, g, qb, g, 3, q, 4, qb, 2, g, 1, g, -67)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, qb, g, g, 2, q, 3, qb, 4, g, 1, g, -68)));
    addSafe(new_ptr((Tree2toNDiagram(5), g, g, q, g, g, 3, q, 2, qb, 4, g, 1, g, -69)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, g, 5, q, 5, qb, 4, g, 1, g, -70)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 3, g, 5, q, 5, qb, 2, g, 1, g, -71)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, q, 5, q, 4, qb, 5, g, 1, g, -72)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, q, g, 3, q, 5, q, 2, qb, 5, g, 1, g, -73)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 4, qb, 4, q, 5, qb, 5, g, 1, g, -74)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, qb, g, 3, qb, 2, q, 5, qb, 5, g, 1, g, -75)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 1, qb, 4, qb, 2, q, 5, qb, 5, g, 4, g, -76)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, q, g, 1, qb, 3, q, 5, qb, 5, g, 2, g, -77)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 1, qb, 4, qb, 2, q, 5, qb, 4, g, 5, g, -78)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, q, g, 1, qb, 3, q, 5, qb, 2, g, 5, g, -79)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 1, qb, 4, g, 2, q, 4, qb, 5, g, 5, g, -80)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, q, g, 1, g, 3, q, 2, qb, 5, g, 5, g, -81)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 1, q, 4, q, 5, q, 2, qb, 5, g, 4, g, -82)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, qb, g, 1, q, 5, q, 3, qb, 5, g, 2, g, -83)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 1, q, 4, q, 5, q, 2, qb, 4, g, 5, g, -84)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, qb, g, 1, q, 5, q, 3, qb, 2, g, 5, g, -85)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 1, q, 4, g, 4, q, 2, qb, 5, g, 5, g, -86)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, qb, g, 1, g, 2, q, 3, qb, 5, g, 5, g, -87)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, q, 5, qb, 2, g, 4, g, -88)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, q, 5, qb, 3, g, 2, g, -89)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, q, 5, q, 4, qb, 2, g, 5, g, -90)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, g, g, 1, q, 5, q, 2, qb, 3, g, 5, g, -91)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, qb, 4, q, 5, qb, 2, g, 5, g, -92)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, g, g, 1, qb, 2, q, 5, qb, 3, g, 5, g, -93)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, g, 5, q, 5, qb, 4, g, 2, g, -94)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, g, g, g, 1, g, 5, q, 5, qb, 2, g, 3, g, -95)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, q, 5, q, 4, qb, 5, g, 2, g, -96)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, qb, g, g, 1, q, 5, q, 2, qb, 5, g, 3, g, -97)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 4, qb, 4, q, 5, qb, 5, g, 2, g, -98)));
    addSafe(new_ptr((Tree2toNDiagram(4), g, q, g, g, 1, qb, 2, q, 5, qb, 5, g, 3, g, -99)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 2, g, 4, q, 4, qb, 5, g, 5, g, -100)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 1, g, 4, q, 4, qb, 5, g, 5, g, -101)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 1, q, 2, qb, 4, q, 5, qb, 4, g, 5, g, -102)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 2, q, 1, qb, 4, q, 5, qb, 4, g, 5, g, -103)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, qb, g, 1, q, 2, qb, 4, q, 5, qb, 5, g, 4, g, -104)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, g, 2, q, 1, qb, 4, q, 5, qb, 5, g, 4, g, -105)));
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEgg2qqbgg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv1at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv1at2[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv1at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv1at4[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv1at5[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv1at6[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv1at7[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv2at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv2at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv3at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv3at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv4at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv4at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv5at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv5at1[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv6at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv6at1[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv7at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv7at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv8at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv8at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv9at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv9at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv10at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv10at1[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv11at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv11at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv12at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv12at1[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv12at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv12at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv12at4[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv12at5[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv12at6[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv12at7[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv13at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv13at1[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv14at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv14at1[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv15at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv15at1[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv15at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv15at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv15at4[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv15at5[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv15at6[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv15at7[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv16at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv16at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv16at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv16at3[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv16at4[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv16at5[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv16at6[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv16at7[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv17at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv17at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv17at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv17at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv17at4[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv17at5[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv17at6[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv17at7[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv18at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv18at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv19at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv19at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv19at2[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv19at3[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv20at0[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv20at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv21at0[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv21at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv21at2[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv21at3[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv22at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv22at1[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv23at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv23at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv23at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv23at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv24at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv24at1[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv25at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv25at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv25at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv25at3[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv26at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv26at1[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv26at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv26at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv27at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv27at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv27at2[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv27at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv28at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv29at0[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv30at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv30at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv30at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv30at3[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv31at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv31at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv32at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv32at1[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv32at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv32at3[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv33at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv33at1[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv34at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv35at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv36at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv37at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv38at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv38at1[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv39at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv39at1[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv39at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv39at3[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv40at0[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv41at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv42at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv42at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv42at2[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv42at3[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv43at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv43at1[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv44at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv44at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv44at2[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv44at3[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv45at0[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv45at1[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv46at0[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv47at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv48at0[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv49at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv50at0[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv50at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv51at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv51at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv51at2[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv51at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv52at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv52at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv53at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv53at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv53at2[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv53at3[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv54at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv54at1[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv55at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv55at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv55at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv55at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv56at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv56at1[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv56at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv56at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv57at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv57at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv57at2[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv57at3[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv58at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv58at1[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv58at2[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv58at3[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv58at4[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv58at5[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv58at6[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv58at7[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv59at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv59at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv59at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv59at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv59at4[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv59at5[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv59at6[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv59at7[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv60at0[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv60at1[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv60at2[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv60at3[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv61at0[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv61at1[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv62at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv62at1[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv62at2[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv62at3[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv63at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv63at1[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv64at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv64at1[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv65at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv65at1[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv65at2[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv65at3[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv66at0[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv66at1[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv67at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv67at1[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv67at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv67at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv68at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv68at1[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv68at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv68at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv69at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv69at1[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv69at2[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv69at3[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv70at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv70at1[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv70at2[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv70at3[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv70at4[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv70at5[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv70at6[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv70at7[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv71at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv71at1[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv71at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv71at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv71at4[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv71at5[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv71at6[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv71at7[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv72at0[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv72at1[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv72at2[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv72at3[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv73at0[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv73at1[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv74at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv74at1[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv74at2[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv74at3[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv75at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv75at1[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv76at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv77at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv78at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv79at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv80at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv80at1[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv81at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv81at1[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv81at2[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv81at3[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv82at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv83at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv84at0[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv85at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv86at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv86at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv87at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv87at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv87at2[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv87at3[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv88at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv88at1[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv88at2[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv88at3[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv88at4[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv88at5[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv88at6[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv88at7[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv89at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv89at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv89at2[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv89at3[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv89at4[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv89at5[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv89at6[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv89at7[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv90at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv90at1[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv90at2[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv90at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv91at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv91at1[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv92at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv92at1[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv92at2[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv92at3[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv93at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv93at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv94at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv94at1[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv94at2[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv94at3[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv94at4[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv94at5[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv94at6[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv94at7[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv95at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv95at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv95at2[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv95at3[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv95at4[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv95at5[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv95at6[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv95at7[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv96at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv96at1[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv96at2[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv96at3[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv97at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv97at1[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv98at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv98at1[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv98at2[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv98at3[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv99at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv99at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv100at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv100at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv100at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv100at3[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv100at4[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv100at5[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv100at6[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv100at7[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv101at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv101at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv101at2[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv101at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv101at4[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv101at5[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv101at6[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv101at7[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv102at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv103at0[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv104at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv105at0[]  = { 1, 4, 0, -1, 3, 2, -999};

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
        + nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv12at1, sizeof(cv12at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at2, sizeof(cv12at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at3, sizeof(cv12at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at4, sizeof(cv12at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at5, sizeof(cv12at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at6, sizeof(cv12at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv12at7, sizeof(cv12at7)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at2, sizeof(cv15at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at3, sizeof(cv15at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at4, sizeof(cv15at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at5, sizeof(cv15at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at6, sizeof(cv15at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at7, sizeof(cv15at7)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv18at1, sizeof(cv18at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -19 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at1, sizeof(cv19at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at2, sizeof(cv19at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv19at3, sizeof(cv19at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -20 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at1, sizeof(cv20at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -21 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at1, sizeof(cv21at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at2, sizeof(cv21at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv21at3, sizeof(cv21at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -22 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at1, sizeof(cv22at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -23 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at1, sizeof(cv23at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at2, sizeof(cv23at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv23at3, sizeof(cv23at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -24 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv24at1, sizeof(cv24at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -25 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at1, sizeof(cv25at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at2, sizeof(cv25at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at3, sizeof(cv25at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -26 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at1, sizeof(cv26at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at2, sizeof(cv26at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv26at3, sizeof(cv26at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -27 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at1, sizeof(cv27at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at2, sizeof(cv27at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv27at3, sizeof(cv27at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -28 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -29 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -30 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at1, sizeof(cv30at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at2, sizeof(cv30at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at3, sizeof(cv30at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -31 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv31at1, sizeof(cv31at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -32 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at1, sizeof(cv32at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at2, sizeof(cv32at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at3, sizeof(cv32at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -33 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at1, sizeof(cv33at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -34 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv34at0, sizeof(cv34at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -35 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv35at0, sizeof(cv35at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -36 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -37 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -38 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv38at1, sizeof(cv38at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -39 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at1, sizeof(cv39at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at2, sizeof(cv39at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv39at3, sizeof(cv39at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -40 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -41 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -42 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at1, sizeof(cv42at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at2, sizeof(cv42at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv42at3, sizeof(cv42at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -43 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv43at0, sizeof(cv43at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv43at1, sizeof(cv43at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -44 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv44at0, sizeof(cv44at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at1, sizeof(cv44at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at2, sizeof(cv44at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv44at3, sizeof(cv44at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -45 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv45at0, sizeof(cv45at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv45at1, sizeof(cv45at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -46 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv46at0, sizeof(cv46at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -47 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -48 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -49 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -50 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv50at0, sizeof(cv50at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv50at1, sizeof(cv50at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -51 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv51at0, sizeof(cv51at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at1, sizeof(cv51at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at2, sizeof(cv51at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv51at3, sizeof(cv51at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -52 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv52at0, sizeof(cv52at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at1, sizeof(cv52at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -53 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at1, sizeof(cv53at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at2, sizeof(cv53at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at3, sizeof(cv53at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -54 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv54at1, sizeof(cv54at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -55 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at1, sizeof(cv55at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at2, sizeof(cv55at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at3, sizeof(cv55at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -56 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at1, sizeof(cv56at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at2, sizeof(cv56at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv56at3, sizeof(cv56at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -57 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at1, sizeof(cv57at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at2, sizeof(cv57at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv57at3, sizeof(cv57at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -58 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv58at0, sizeof(cv58at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at1, sizeof(cv58at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at2, sizeof(cv58at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at3, sizeof(cv58at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at4, sizeof(cv58at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at5, sizeof(cv58at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at6, sizeof(cv58at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv58at7, sizeof(cv58at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -59 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv59at0, sizeof(cv59at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at1, sizeof(cv59at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at2, sizeof(cv59at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at3, sizeof(cv59at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at4, sizeof(cv59at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at5, sizeof(cv59at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at6, sizeof(cv59at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv59at7, sizeof(cv59at7)/sizeof(int)), i);   
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
          nloJetAmplitude()->colourOrdered2(cv64at0, sizeof(cv64at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv64at1, sizeof(cv64at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -65 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv65at0, sizeof(cv65at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at1, sizeof(cv65at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at2, sizeof(cv65at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv65at3, sizeof(cv65at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -66 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv66at0, sizeof(cv66at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv66at1, sizeof(cv66at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -67 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv67at0, sizeof(cv67at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at1, sizeof(cv67at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at2, sizeof(cv67at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at3, sizeof(cv67at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -68 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at1, sizeof(cv68at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at2, sizeof(cv68at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at3, sizeof(cv68at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -69 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at1, sizeof(cv69at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at2, sizeof(cv69at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv69at3, sizeof(cv69at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -70 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv70at0, sizeof(cv70at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at1, sizeof(cv70at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at2, sizeof(cv70at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at3, sizeof(cv70at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at4, sizeof(cv70at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at5, sizeof(cv70at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at6, sizeof(cv70at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at7, sizeof(cv70at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -71 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv71at0, sizeof(cv71at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at1, sizeof(cv71at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at2, sizeof(cv71at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at3, sizeof(cv71at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at4, sizeof(cv71at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at5, sizeof(cv71at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at6, sizeof(cv71at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv71at7, sizeof(cv71at7)/sizeof(int)), i);   
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
          nloJetAmplitude()->colourOrdered2(cv76at0, sizeof(cv76at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -77 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv77at0, sizeof(cv77at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -78 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv78at0, sizeof(cv78at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -79 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv79at0, sizeof(cv79at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -80 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv80at0, sizeof(cv80at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv80at1, sizeof(cv80at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -81 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv81at0, sizeof(cv81at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at1, sizeof(cv81at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at2, sizeof(cv81at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv81at3, sizeof(cv81at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -82 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv82at0, sizeof(cv82at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -83 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv83at0, sizeof(cv83at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -84 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv84at0, sizeof(cv84at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -85 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv85at0, sizeof(cv85at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -86 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv86at0, sizeof(cv86at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv86at1, sizeof(cv86at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -87 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv87at0, sizeof(cv87at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at1, sizeof(cv87at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at2, sizeof(cv87at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv87at3, sizeof(cv87at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -88 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv88at0, sizeof(cv88at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at1, sizeof(cv88at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at2, sizeof(cv88at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at3, sizeof(cv88at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at4, sizeof(cv88at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at5, sizeof(cv88at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at6, sizeof(cv88at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv88at7, sizeof(cv88at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -89 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv89at0, sizeof(cv89at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at1, sizeof(cv89at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at2, sizeof(cv89at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at3, sizeof(cv89at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at4, sizeof(cv89at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at5, sizeof(cv89at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at6, sizeof(cv89at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv89at7, sizeof(cv89at7)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv94at3, sizeof(cv94at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at4, sizeof(cv94at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at5, sizeof(cv94at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at6, sizeof(cv94at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv94at7, sizeof(cv94at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -95 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv95at0, sizeof(cv95at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at1, sizeof(cv95at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at2, sizeof(cv95at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at3, sizeof(cv95at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at4, sizeof(cv95at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at5, sizeof(cv95at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at6, sizeof(cv95at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv95at7, sizeof(cv95at7)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv100at3, sizeof(cv100at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at4, sizeof(cv100at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at5, sizeof(cv100at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at6, sizeof(cv100at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv100at7, sizeof(cv100at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -101 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv101at0, sizeof(cv101at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at1, sizeof(cv101at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at2, sizeof(cv101at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at3, sizeof(cv101at3)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at4, sizeof(cv101at4)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at5, sizeof(cv101at5)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at6, sizeof(cv101at6)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv101at7, sizeof(cv101at7)/sizeof(int)), i);   
    else if( diags[i]->id() == -102 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv102at0, sizeof(cv102at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -103 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv103at0, sizeof(cv103at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -104 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv104at0, sizeof(cv104at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -105 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv105at0, sizeof(cv105at0)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEgg2qqbgg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEgg2qqbgg

  static const ColourLines diag1[8] = { 
    ColourLines("1 3 5 9, -1 2, -2 -3 -4 -7, 6 4 -5 -8, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 5 9, 6 4 -5 -8, 8 -9"), 
    ColourLines("1 3 5 8, -1 2, -2 -3 -4 -7, 6 4 -5 -9, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 5 8, 6 4 -5 -9, -8 9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -8, -7 -4 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -5 -8, 2 3 4 6, -7 -4 5 9, 8 -9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -5 -9, -7 -4 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -5 -9, 2 3 4 6, -7 -4 5 8, -8 9")
  }; 
  static const ColourLines diag2[2] = { 
    ColourLines("1 3 4 8, -1 2, -2 -3 -5 -9, 6 -8, -7 9"), 
    ColourLines("1 -2, -1 -3 -5 -9, 2 3 4 8, 6 -8, -7 9")
  }; 
  static const ColourLines diag3[2] = { 
    ColourLines("1 3 4 9, -1 2, -2 -3 -5 -8, 6 -9, -7 8"), 
    ColourLines("1 -2, -1 -3 -5 -8, 2 3 4 9, 6 -9, -7 8")
  }; 
  static const ColourLines diag4[2] = { 
    ColourLines("1 6, -1 -2 -5 -8, 4 -3 5 9, -4 -7, 8 -9"), 
    ColourLines("1 6, -1 -2 -5 -9, 4 -3 5 8, -4 -7, -8 9")
  }; 
  static const ColourLines diag5[2] = { 
    ColourLines("1 6, -1 -2 -3 4, -4 -8, -7 9, 8 3 -5 -9"), 
    ColourLines("1 6, -1 -2 -3 -8, 4 8, -4 3 -5 -9, -7 9")
  }; 
  static const ColourLines diag6[2] = { 
    ColourLines("1 6, -1 -2 -3 4, -4 -9, -7 8, -8 -5 3 9"), 
    ColourLines("1 6, -1 -2 -3 -9, 4 9, -4 3 -5 -8, -7 8")
  }; 
  static const ColourLines diag7[2] = { 
    ColourLines("1 2 5 9, -1 -7, 4 6, -4 3 -5 -8, 8 -9"), 
    ColourLines("1 2 5 8, -1 -7, 4 6, -4 3 -5 -9, -8 9")
  }; 
  static const ColourLines diag8[2] = { 
    ColourLines("1 2 3 8, -1 -7, 4 -3 5 9, -4 -8, 6 -9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 8, 6 -9, -8 -3 5 9")
  }; 
  static const ColourLines diag9[2] = { 
    ColourLines("1 2 3 9, -1 -7, 4 -3 5 8, -4 -9, 6 -8"), 
    ColourLines("1 2 3 -4, -1 -7, 4 9, 6 -8, 8 5 -3 -9")
  }; 
  static const ColourLines diag10[2] = { 
    ColourLines("1 2 3 -4, -1 -8, 4 6, -7 9, 8 -2 -5 -9"), 
    ColourLines("1 8, -1 -2 -5 -9, 4 6, -4 3 2 -8, -7 9")
  }; 
  static const ColourLines diag11[2] = { 
    ColourLines("1 2 5 9, -1 -8, 4 -3 -2 8, -4 -7, 6 -9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -7, 6 -9, -8 2 5 9")
  }; 
  static const ColourLines diag12[8] = { 
    ColourLines("1 2 5 6, -1 -8, 4 -3 -2 8, -4 -9, -7 -5 3 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -9, 6 5 2 -8, -7 -5 3 9"), 
    ColourLines("1 2 5 6, -1 -8, 4 9, -4 3 -5 -7, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 4 9, -4 3 -5 -7, 6 5 2 -8"), 
    ColourLines("1 2 3 9, -1 -8, 4 -3 5 6, -4 -9, -7 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -7, 4 -3 5 6, -4 -9, -8 2 3 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 9, 6 5 -3 -9, -7 -5 -2 8"), 
    ColourLines("1 8, -1 -2 -5 -7, 4 9, -4 3 2 -8, 6 5 -3 -9")
  }; 
  static const ColourLines diag13[2] = { 
    ColourLines("1 2 3 -4, -1 -9, 4 6, -7 8, -8 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -8, 4 6, -4 3 2 -9, -7 8")
  }; 
  static const ColourLines diag14[2] = { 
    ColourLines("1 2 5 8, -1 -9, 4 -3 -2 9, -4 -7, 6 -8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -7, 6 -8, 8 5 2 -9")
  }; 
  static const ColourLines diag15[8] = { 
    ColourLines("1 2 5 6, -1 -9, 4 -3 -2 9, -4 -8, -7 -5 3 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -8, 6 5 2 -9, -7 -5 3 8"), 
    ColourLines("1 2 5 6, -1 -9, 4 8, -4 3 -5 -7, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 4 8, -4 3 -5 -7, 6 5 2 -9"), 
    ColourLines("1 2 3 8, -1 -9, 4 -3 5 6, -4 -8, -7 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -7, 4 -3 5 6, -4 -8, 8 3 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 8, 6 5 -3 -8, -7 -5 -2 9"), 
    ColourLines("1 9, -1 -2 -5 -7, 4 8, -4 3 2 -9, 6 5 -3 -8")
  }; 
  static const ColourLines diag16[8] = { 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -9, 6 5 4 -8, -7 -5 9"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 8, 6 5 4 -8, -7 -5 9"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -8, -7 -5 9, 8 -4 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 6, -7 -5 9, 8 -4 -9"), 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -7, 6 5 -9, -8 4 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 8, 6 5 -9, -8 4 9"), 
    ColourLines("1 3 4 9, -1 2, -2 -3 -8, 6 5 -9, -7 -5 -4 8"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 9, 6 5 -9, -7 -5 -4 8")
  }; 
  static const ColourLines diag17[8] = { 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -8, 6 5 4 -9, -7 -5 8"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 9, 6 5 4 -9, -7 -5 8"), 
    ColourLines("1 3 4 5 6, -1 2, -2 -3 -9, -7 -5 8, -8 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 6, -7 -5 8, -8 -4 9"), 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -7, 6 5 -8, 8 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -7, 2 3 9, 6 5 -8, 8 4 -9"), 
    ColourLines("1 3 4 8, -1 2, -2 -3 -9, 6 5 -8, -7 -5 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 8, 6 5 -8, -7 -5 -4 9")
  }; 
  static const ColourLines diag18[2] = { 
    ColourLines("1 3 4 9, -1 2, -2 -3 -7, 6 -8, 8 5 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 9, 6 -8, 8 5 -9")
  }; 
  static const ColourLines diag19[4] = { 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -7, 6 -8, 8 5 4 -9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 9, 6 -8, 8 5 4 -9"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -9, 6 -8, -7 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 5 8, 6 -8, -7 -4 9")
  }; 
  static const ColourLines diag20[2] = { 
    ColourLines("1 3 4 8, -1 2, -2 -3 -7, 6 -9, -8 5 9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 8, 6 -9, -8 5 9")
  }; 
  static const ColourLines diag21[4] = { 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -7, 6 -9, -8 4 5 9"), 
    ColourLines("1 -2, -1 -3 -4 -7, 2 3 8, 6 -9, -8 4 5 9"), 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -8, 6 -9, -7 -4 8"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 5 9, 6 -9, -7 -4 8")
  }; 
  static const ColourLines diag22[2] = { 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -9, -7 8, -8 -5 9"), 
    ColourLines("1 -2, -1 -3 -4 -9, 2 3 6, -7 8, -8 -5 9")
  }; 
  static const ColourLines diag23[4] = { 
    ColourLines("1 3 9, -1 2, -2 -3 -4 -5 -8, 6 4 -9, -7 8"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 9, 6 4 -9, -7 8"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -9, -7 8, -8 -5 -4 9"), 
    ColourLines("1 -2, -1 -3 -9, 2 3 4 6, -7 8, -8 -5 -4 9")
  }; 
  static const ColourLines diag24[2] = { 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -8, -7 9, 8 -5 -9"), 
    ColourLines("1 -2, -1 -3 -4 -8, 2 3 6, -7 9, 8 -5 -9")
  }; 
  static const ColourLines diag25[4] = { 
    ColourLines("1 3 8, -1 2, -2 -3 -4 -5 -9, 6 4 -8, -7 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 8, 6 4 -8, -7 9"), 
    ColourLines("1 3 4 6, -1 2, -2 -3 -8, -7 9, 8 -4 -5 -9"), 
    ColourLines("1 -2, -1 -3 -8, 2 3 4 6, -7 9, 8 -4 -5 -9")
  }; 
  static const ColourLines diag26[4] = { 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -8, 2 3 6, -7 5 9, 8 -9"), 
    ColourLines("1 3 6, -1 2, -2 -3 -4 -5 -9, -7 5 8, -8 9"), 
    ColourLines("1 -2, -1 -3 -4 -5 -9, 2 3 6, -7 5 8, -8 9")
  }; 
  static const ColourLines diag27[4] = { 
    ColourLines("1 3 4 5 9, -1 2, -2 -3 -7, 6 -5 -8, 8 -9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 9, 6 -5 -8, 8 -9"), 
    ColourLines("1 3 4 5 8, -1 2, -2 -3 -7, 6 -5 -9, -8 9"), 
    ColourLines("1 -2, -1 -3 -7, 2 3 4 5 8, 6 -5 -9, -8 9")
  }; 
  static const ColourLines diag28[1] = { 
    ColourLines("1 6, -1 -2 -8, 5 -4 9, -5 -7, 8 -3 -9")
  }; 
  static const ColourLines diag29[1] = { 
    ColourLines("1 6, -1 -2 -9, 5 -4 8, -5 -7, -8 -3 9")
  }; 
  static const ColourLines diag30[4] = { 
    ColourLines("1 6, -1 -2 -3 -9, 5 -4 9, -5 -8, -7 3 4 8"), 
    ColourLines("1 6, -1 -2 -3 -9, 5 8, -5 4 3 -7, -8 -4 9"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -8, -7 3 9, 8 4 -9"), 
    ColourLines("1 6, -1 -2 -3 -4 -8, 5 8, -5 4 -9, -7 3 9")
  }; 
  static const ColourLines diag31[2] = { 
    ColourLines("1 6, -1 -2 -9, 5 -4 -3 9, -5 -8, -7 4 8"), 
    ColourLines("1 6, -1 -2 -9, 5 8, -5 4 -7, -8 -4 -3 9")
  }; 
  static const ColourLines diag32[4] = { 
    ColourLines("1 6, -1 -2 -3 -8, 5 -4 8, -5 -9, -7 3 4 9"), 
    ColourLines("1 6, -1 -2 -3 -8, 5 9, -5 4 3 -7, 8 -4 -9"), 
    ColourLines("1 6, -1 -2 -3 -4 5, -5 -9, -7 3 8, -8 4 9"), 
    ColourLines("1 6, -1 -2 -3 -4 -9, 5 9, -5 4 -8, -7 3 8")
  }; 
  static const ColourLines diag33[2] = { 
    ColourLines("1 6, -1 -2 -8, 5 -4 -3 8, -5 -9, -7 4 9"), 
    ColourLines("1 6, -1 -2 -8, 5 9, -5 4 -7, 8 -3 -4 -9")
  }; 
  static const ColourLines diag34[1] = { 
    ColourLines("1 6, -1 -2 3, -3 -4 -9, -7 8, -8 -5 9")
  }; 
  static const ColourLines diag35[1] = { 
    ColourLines("1 6, -1 -2 -9, 4 -3 9, -4 -5 -8, -7 8")
  }; 
  static const ColourLines diag36[1] = { 
    ColourLines("1 6, -1 -2 3, -3 -4 -8, -7 9, 8 -5 -9")
  }; 
  static const ColourLines diag37[1] = { 
    ColourLines("1 6, -1 -2 -8, 4 -3 8, -4 -5 -9, -7 9")
  }; 
  static const ColourLines diag38[2] = { 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -8, -7 5 9, 8 -9"), 
    ColourLines("1 6, -1 -2 3, -3 -4 -5 -9, -7 5 8, -8 9")
  }; 
  static const ColourLines diag39[4] = { 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -8, -7 3 5 9, 8 -9"), 
    ColourLines("1 6, -1 -2 -3 4, -4 -5 -9, -7 3 5 8, -8 9"), 
    ColourLines("1 6, -1 -2 -3 -5 -8, 4 5 9, -4 3 -7, 8 -9"), 
    ColourLines("1 6, -1 -2 -3 -5 -9, 4 5 8, -4 3 -7, -8 9")
  }; 
  static const ColourLines diag40[1] = { 
    ColourLines("1 2 8, -1 -7, 5 6, -5 4 -9, -8 3 9")
  }; 
  static const ColourLines diag41[1] = { 
    ColourLines("1 2 9, -1 -7, 5 6, -5 4 -8, 8 3 -9")
  }; 
  static const ColourLines diag42[4] = { 
    ColourLines("1 2 3 4 8, -1 -7, 5 -4 9, -5 -8, 6 -3 -9"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 8, 6 -3 -9, -8 -4 9"), 
    ColourLines("1 2 3 9, -1 -7, 5 -4 -3 6, -5 -8, 8 4 -9"), 
    ColourLines("1 2 3 9, -1 -7, 5 8, -5 4 -9, 6 -3 -4 -8")
  }; 
  static const ColourLines diag43[2] = { 
    ColourLines("1 2 9, -1 -7, 5 -4 6, -5 -8, 8 4 3 -9"), 
    ColourLines("1 2 9, -1 -7, 5 8, -5 4 3 -9, 6 -4 -8")
  }; 
  static const ColourLines diag44[4] = { 
    ColourLines("1 2 3 4 9, -1 -7, 5 -4 8, -5 -9, 6 -3 -8"), 
    ColourLines("1 2 3 4 -5, -1 -7, 5 9, 6 -3 -8, 8 -4 -9"), 
    ColourLines("1 2 3 8, -1 -7, 5 -4 -3 6, -5 -9, -8 4 9"), 
    ColourLines("1 2 3 8, -1 -7, 5 9, -5 4 -8, 6 -3 -4 -9")
  }; 
  static const ColourLines diag45[2] = { 
    ColourLines("1 2 8, -1 -7, 5 -4 6, -5 -9, -8 3 4 9"), 
    ColourLines("1 2 8, -1 -7, 5 9, -5 4 3 -8, 6 -4 -9")
  }; 
  static const ColourLines diag46[1] = { 
    ColourLines("1 2 -3, -1 -7, 3 4 9, 6 -8, 8 5 -9")
  }; 
  static const ColourLines diag47[1] = { 
    ColourLines("1 2 9, -1 -7, 4 5 8, -4 3 -9, 6 -8")
  }; 
  static const ColourLines diag48[1] = { 
    ColourLines("1 2 -3, -1 -7, 3 4 8, 6 -9, -8 5 9")
  }; 
  static const ColourLines diag49[1] = { 
    ColourLines("1 2 8, -1 -7, 4 5 9, -4 3 -8, 6 -9")
  }; 
  static const ColourLines diag50[2] = { 
    ColourLines("1 2 -3, -1 -7, 3 4 5 9, 6 -5 -8, 8 -9"), 
    ColourLines("1 2 -3, -1 -7, 3 4 5 8, 6 -5 -9, -8 9")
  }; 
  static const ColourLines diag51[4] = { 
    ColourLines("1 2 3 5 9, -1 -7, 4 -3 6, -4 -5 -8, 8 -9"), 
    ColourLines("1 2 3 5 8, -1 -7, 4 -3 6, -4 -5 -9, -8 9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 9, 6 -3 -5 -8, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -7, 4 5 8, 6 -3 -5 -9, -8 9")
  }; 
  static const ColourLines diag52[2] = { 
    ColourLines("1 2 3 9, -1 -8, 5 6, -5 4 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 6, -5 4 -9, -8 2 3 9")
  }; 
  static const ColourLines diag53[4] = { 
    ColourLines("1 2 9, -1 -8, 5 6, -5 4 3 -9, -7 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -7, 5 6, -5 4 3 -9, -8 2 9"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 6, -7 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 6, -5 4 3 2 -8, -7 -3 9")
  }; 
  static const ColourLines diag54[2] = { 
    ColourLines("1 2 6, -1 -8, 5 -4 9, -5 -7, 8 -2 -3 -9"), 
    ColourLines("1 8, -1 -2 -3 -9, 5 -4 9, -5 -7, 6 2 -8")
  }; 
  static const ColourLines diag55[4] = { 
    ColourLines("1 2 9, -1 -8, 5 -4 -3 -2 8, -5 -7, 6 3 -9"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -7, 6 3 -9, -8 2 9"), 
    ColourLines("1 2 3 6, -1 -8, 5 -4 -3 9, -5 -7, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 5 -4 -3 9, -5 -7, 6 3 2 -8")
  }; 
  static const ColourLines diag56[4] = { 
    ColourLines("1 2 6, -1 -8, 5 -4 -3 -2 8, -5 -9, -7 4 9"), 
    ColourLines("1 8, -1 -2 -3 -4 5, -5 -9, 6 2 -8, -7 4 9"), 
    ColourLines("1 2 6, -1 -8, 5 9, -5 4 -7, 8 -2 -3 -4 -9"), 
    ColourLines("1 8, -1 -2 -3 -4 -9, 5 9, -5 4 -7, 6 2 -8")
  }; 
  static const ColourLines diag57[4] = { 
    ColourLines("1 2 3 4 9, -1 -8, 5 -4 6, -5 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 -4 6, -5 -9, -8 2 3 4 9"), 
    ColourLines("1 2 3 4 -5, -1 -8, 5 9, 6 -4 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 5 9, -5 4 3 2 -8, 6 -4 -9")
  }; 
  static const ColourLines diag58[8] = { 
    ColourLines("1 2 -3, -1 -8, 3 4 5 6, -7 -5 9, 8 -2 -4 -9"), 
    ColourLines("1 8, -1 -2 -4 -9, 3 4 5 6, -3 2 -8, -7 -5 9"), 
    ColourLines("1 2 4 5 6, -1 -8, 3 -2 8, -3 -4 -9, -7 -5 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -9, 6 5 4 2 -8, -7 -5 9"), 
    ColourLines("1 2 -3, -1 -8, 3 4 9, 6 5 -9, -7 -5 -4 -2 8"), 
    ColourLines("1 8, -1 -2 -4 -5 -7, 3 4 9, -3 2 -8, 6 5 -9"), 
    ColourLines("1 2 4 9, -1 -8, 3 -2 8, -3 -4 -5 -7, 6 5 -9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -7, 6 5 -9, -8 2 4 9")
  }; 
  static const ColourLines diag59[8] = { 
    ColourLines("1 2 9, -1 -8, 4 -3 -2 8, -4 -5 -7, 6 5 3 -9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -7, 6 5 3 -9, -8 2 9"), 
    ColourLines("1 2 3 5 6, -1 -8, 4 -3 9, -4 -5 -7, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 4 -3 9, -4 -5 -7, 6 5 3 2 -8"), 
    ColourLines("1 2 9, -1 -8, 4 5 6, -4 3 -9, -7 -5 -3 -2 8"), 
    ColourLines("1 8, -1 -2 -3 -5 -7, 4 5 6, -4 3 -9, -8 2 9"), 
    ColourLines("1 2 3 -4, -1 -8, 4 5 6, -7 -5 -3 9, 8 -2 -9"), 
    ColourLines("1 8, -1 -2 -9, 4 5 6, -4 3 2 -8, -7 -5 -3 9")
  }; 
  static const ColourLines diag60[4] = { 
    ColourLines("1 2 -3, -1 -8, 3 4 5 9, 6 -9, -7 -4 -2 8"), 
    ColourLines("1 8, -1 -2 -4 -7, 3 4 5 9, -3 2 -8, 6 -9"), 
    ColourLines("1 2 4 5 9, -1 -8, 3 -2 8, -3 -4 -7, 6 -9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -7, 6 -9, -8 2 4 5 9")
  }; 
  static const ColourLines diag61[2] = { 
    ColourLines("1 2 3 -4, -1 -8, 4 5 9, 6 -9, -7 -2 8"), 
    ColourLines("1 8, -1 -2 -7, 4 5 9, -4 3 2 -8, 6 -9")
  }; 
  static const ColourLines diag62[4] = { 
    ColourLines("1 2 -3, -1 -8, 3 4 6, -7 9, 8 -2 -4 -5 -9"), 
    ColourLines("1 8, -1 -2 -4 -5 -9, 3 4 6, -3 2 -8, -7 9"), 
    ColourLines("1 2 4 6, -1 -8, 3 -2 8, -3 -4 -5 -9, -7 9"), 
    ColourLines("1 8, -1 -2 3, -3 -4 -5 -9, 6 4 2 -8, -7 9")
  }; 
  static const ColourLines diag63[2] = { 
    ColourLines("1 2 6, -1 -8, 4 -3 -2 8, -4 -5 -9, -7 9"), 
    ColourLines("1 8, -1 -2 -3 4, -4 -5 -9, 6 2 -8, -7 9")
  }; 
  static const ColourLines diag64[2] = { 
    ColourLines("1 2 3 8, -1 -9, 5 6, -5 4 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 6, -5 4 -8, 8 3 2 -9")
  }; 
  static const ColourLines diag65[4] = { 
    ColourLines("1 2 8, -1 -9, 5 6, -5 4 3 -8, -7 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -7, 5 6, -5 4 3 -8, 8 2 -9"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 6, -7 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 6, -5 4 3 2 -9, -7 -3 8")
  }; 
  static const ColourLines diag66[2] = { 
    ColourLines("1 2 6, -1 -9, 5 -4 8, -5 -7, -8 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -8, 5 -4 8, -5 -7, 6 2 -9")
  }; 
  static const ColourLines diag67[4] = { 
    ColourLines("1 2 8, -1 -9, 5 -4 -3 -2 9, -5 -7, 6 3 -8"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -7, 6 3 -8, 8 2 -9"), 
    ColourLines("1 2 3 6, -1 -9, 5 -4 -3 8, -5 -7, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 5 -4 -3 8, -5 -7, 6 3 2 -9")
  }; 
  static const ColourLines diag68[4] = { 
    ColourLines("1 2 6, -1 -9, 5 -4 -3 -2 9, -5 -8, -7 4 8"), 
    ColourLines("1 9, -1 -2 -3 -4 5, -5 -8, 6 2 -9, -7 4 8"), 
    ColourLines("1 2 6, -1 -9, 5 8, -5 4 -7, -8 -4 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -4 -8, 5 8, -5 4 -7, 6 2 -9")
  }; 
  static const ColourLines diag69[4] = { 
    ColourLines("1 2 3 4 8, -1 -9, 5 -4 6, -5 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 -4 6, -5 -8, 8 4 3 2 -9"), 
    ColourLines("1 2 3 4 -5, -1 -9, 5 8, 6 -4 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 5 8, -5 4 3 2 -9, 6 -4 -8")
  }; 
  static const ColourLines diag70[8] = { 
    ColourLines("1 2 -3, -1 -9, 3 4 5 6, -7 -5 8, -8 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -8, 3 4 5 6, -3 2 -9, -7 -5 8"), 
    ColourLines("1 2 4 5 6, -1 -9, 3 -2 9, -3 -4 -8, -7 -5 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -8, 6 5 4 2 -9, -7 -5 8"), 
    ColourLines("1 2 -3, -1 -9, 3 4 8, 6 5 -8, -7 -5 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -5 -7, 3 4 8, -3 2 -9, 6 5 -8"), 
    ColourLines("1 2 4 8, -1 -9, 3 -2 9, -3 -4 -5 -7, 6 5 -8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -7, 6 5 -8, 8 4 2 -9")
  }; 
  static const ColourLines diag71[8] = { 
    ColourLines("1 2 8, -1 -9, 4 -3 -2 9, -4 -5 -7, 6 5 3 -8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -7, 6 5 3 -8, 8 2 -9"), 
    ColourLines("1 2 3 5 6, -1 -9, 4 -3 8, -4 -5 -7, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 4 -3 8, -4 -5 -7, 6 5 3 2 -9"), 
    ColourLines("1 2 8, -1 -9, 4 5 6, -4 3 -8, -7 -5 -3 -2 9"), 
    ColourLines("1 9, -1 -2 -3 -5 -7, 4 5 6, -4 3 -8, 8 2 -9"), 
    ColourLines("1 2 3 -4, -1 -9, 4 5 6, -7 -5 -3 8, -8 -2 9"), 
    ColourLines("1 9, -1 -2 -8, 4 5 6, -4 3 2 -9, -7 -5 -3 8")
  }; 
  static const ColourLines diag72[4] = { 
    ColourLines("1 2 -3, -1 -9, 3 4 5 8, 6 -8, -7 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -7, 3 4 5 8, -3 2 -9, 6 -8"), 
    ColourLines("1 2 4 5 8, -1 -9, 3 -2 9, -3 -4 -7, 6 -8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -7, 6 -8, 8 5 4 2 -9")
  }; 
  static const ColourLines diag73[2] = { 
    ColourLines("1 2 3 -4, -1 -9, 4 5 8, 6 -8, -7 -2 9"), 
    ColourLines("1 9, -1 -2 -7, 4 5 8, -4 3 2 -9, 6 -8")
  }; 
  static const ColourLines diag74[4] = { 
    ColourLines("1 2 -3, -1 -9, 3 4 6, -7 8, -8 -5 -4 -2 9"), 
    ColourLines("1 9, -1 -2 -4 -5 -8, 3 4 6, -3 2 -9, -7 8"), 
    ColourLines("1 2 4 6, -1 -9, 3 -2 9, -3 -4 -5 -8, -7 8"), 
    ColourLines("1 9, -1 -2 3, -3 -4 -5 -8, 6 4 2 -9, -7 8")
  }; 
  static const ColourLines diag75[2] = { 
    ColourLines("1 2 6, -1 -9, 4 -3 -2 9, -4 -5 -8, -7 8"), 
    ColourLines("1 9, -1 -2 -3 4, -4 -5 -8, 6 2 -9, -7 8")
  }; 
  static const ColourLines diag76[1] = { 
    ColourLines("1 2 -3, -1 -4 -9, 3 6, -7 8, -8 -5 9")
  }; 
  static const ColourLines diag77[1] = { 
    ColourLines("1 2 9, -1 -5 -8, 4 6, -4 3 -9, -7 8")
  }; 
  static const ColourLines diag78[1] = { 
    ColourLines("1 2 -3, -1 -4 -8, 3 6, -7 9, 8 -5 -9")
  }; 
  static const ColourLines diag79[1] = { 
    ColourLines("1 2 8, -1 -5 -9, 4 6, -4 3 -8, -7 9")
  }; 
  static const ColourLines diag80[2] = { 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 6, -7 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 6, -7 5 8, -8 9")
  }; 
  static const ColourLines diag81[4] = { 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 6, -7 -2 5 9, 8 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 6, -7 -2 5 8, -8 9"), 
    ColourLines("1 5 9, -1 -2 -7, 4 6, -4 3 2 -5 -8, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -7, 4 6, -4 3 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag82[1] = { 
    ColourLines("1 4 9, -1 -2 3, -3 -7, 6 -8, 8 5 -9")
  }; 
  static const ColourLines diag83[1] = { 
    ColourLines("1 5 8, -1 -2 -9, 4 -3 9, -4 -7, 6 -8")
  }; 
  static const ColourLines diag84[1] = { 
    ColourLines("1 4 8, -1 -2 3, -3 -7, 6 -9, -8 5 9")
  }; 
  static const ColourLines diag85[1] = { 
    ColourLines("1 5 9, -1 -2 -8, 4 -3 8, -4 -7, 6 -9")
  }; 
  static const ColourLines diag86[2] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -7, 6 -5 -8, 8 -9"), 
    ColourLines("1 4 5 8, -1 -2 3, -3 -7, 6 -5 -9, -8 9")
  }; 
  static const ColourLines diag87[4] = { 
    ColourLines("1 2 6, -1 -5 -8, 4 -3 -2 5 9, -4 -7, 8 -9"), 
    ColourLines("1 2 6, -1 -5 -9, 4 -3 -2 5 8, -4 -7, -8 9"), 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -7, 6 2 -5 -8, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -7, 6 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag88[8] = { 
    ColourLines("1 4 5 6, -1 -2 3, -3 -8, -7 -5 9, 8 2 -4 -9"), 
    ColourLines("1 4 5 6, -1 -2 -8, 3 8, -3 2 -4 -9, -7 -5 9"), 
    ColourLines("1 2 8, -1 -4 -9, 3 -2 4 5 6, -3 -8, -7 -5 9"), 
    ColourLines("1 2 -3, -1 -4 -9, 3 8, 6 5 4 -2 -8, -7 -5 9"), 
    ColourLines("1 4 9, -1 -2 3, -3 -8, 6 5 -9, -7 -5 -4 2 8"), 
    ColourLines("1 4 9, -1 -2 -8, 3 8, -3 2 -4 -5 -7, 6 5 -9"), 
    ColourLines("1 2 8, -1 -4 -5 -7, 3 -2 4 9, -3 -8, 6 5 -9"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 8, 6 5 -9, -8 -2 4 9")
  }; 
  static const ColourLines diag89[8] = { 
    ColourLines("1 2 3 8, -1 -5 -7, 4 -3 9, -4 -8, 6 5 -2 -9"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 8, 6 5 -2 -9, -8 -3 9"), 
    ColourLines("1 2 9, -1 -5 -7, 4 -3 -2 5 6, -4 -8, 8 3 -9"), 
    ColourLines("1 2 9, -1 -5 -7, 4 8, -4 3 -9, 6 5 -2 -3 -8"), 
    ColourLines("1 5 6, -1 -2 -9, 4 -3 9, -4 -8, -7 -5 2 3 8"), 
    ColourLines("1 5 6, -1 -2 -9, 4 8, -4 3 2 -5 -7, -8 -3 9"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -8, -7 -5 2 9, 8 3 -9"), 
    ColourLines("1 5 6, -1 -2 -3 -8, 4 8, -4 3 -9, -7 -5 2 9")
  }; 
  static const ColourLines diag90[4] = { 
    ColourLines("1 4 5 9, -1 -2 3, -3 -8, 6 -9, -7 -4 2 8"), 
    ColourLines("1 4 5 9, -1 -2 -8, 3 8, -3 2 -4 -7, 6 -9"), 
    ColourLines("1 2 8, -1 -4 -7, 3 -2 4 5 9, -3 -8, 6 -9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 8, 6 -9, -8 -2 4 5 9")
  }; 
  static const ColourLines diag91[2] = { 
    ColourLines("1 5 9, -1 -2 -3 4, -4 -8, 6 -9, -7 3 8"), 
    ColourLines("1 5 9, -1 -2 -3 -8, 4 8, -4 3 -7, 6 -9")
  }; 
  static const ColourLines diag92[4] = { 
    ColourLines("1 4 6, -1 -2 3, -3 -8, -7 9, 8 2 -4 -5 -9"), 
    ColourLines("1 4 6, -1 -2 -8, 3 8, -3 2 -4 -5 -9, -7 9"), 
    ColourLines("1 2 8, -1 -4 -5 -9, 3 -2 4 6, -3 -8, -7 9"), 
    ColourLines("1 2 -3, -1 -4 -5 -9, 3 8, 6 4 -2 -8, -7 9")
  }; 
  static const ColourLines diag93[2] = { 
    ColourLines("1 2 3 8, -1 -5 -9, 4 -3 6, -4 -8, -7 9"), 
    ColourLines("1 2 3 -4, -1 -5 -9, 4 8, 6 -3 -8, -7 9")
  }; 
  static const ColourLines diag94[8] = { 
    ColourLines("1 4 5 6, -1 -2 3, -3 -9, -7 -5 8, -8 -4 2 9"), 
    ColourLines("1 4 5 6, -1 -2 -9, 3 9, -3 2 -4 -8, -7 -5 8"), 
    ColourLines("1 2 9, -1 -4 -8, 3 -2 4 5 6, -3 -9, -7 -5 8"), 
    ColourLines("1 2 -3, -1 -4 -8, 3 9, 6 5 4 -2 -9, -7 -5 8"), 
    ColourLines("1 4 8, -1 -2 3, -3 -9, 6 5 -8, -7 -5 -4 2 9"), 
    ColourLines("1 4 8, -1 -2 -9, 3 9, -3 2 -4 -5 -7, 6 5 -8"), 
    ColourLines("1 2 9, -1 -4 -5 -7, 3 -2 4 8, -3 -9, 6 5 -8"), 
    ColourLines("1 2 -3, -1 -4 -5 -7, 3 9, 6 5 -8, 8 4 -2 -9")
  }; 
  static const ColourLines diag95[8] = { 
    ColourLines("1 2 3 9, -1 -5 -7, 4 -3 8, -4 -9, 6 5 -2 -8"), 
    ColourLines("1 2 3 -4, -1 -5 -7, 4 9, 6 5 -2 -8, 8 -3 -9"), 
    ColourLines("1 2 8, -1 -5 -7, 4 -3 -2 5 6, -4 -9, -8 3 9"), 
    ColourLines("1 2 8, -1 -5 -7, 4 9, -4 3 -8, 6 5 -2 -3 -9"), 
    ColourLines("1 5 6, -1 -2 -8, 4 -3 8, -4 -9, -7 -5 2 3 9"), 
    ColourLines("1 5 6, -1 -2 -8, 4 9, -4 3 2 -5 -7, 8 -3 -9"), 
    ColourLines("1 5 6, -1 -2 -3 4, -4 -9, -7 -5 2 8, -8 3 9"), 
    ColourLines("1 5 6, -1 -2 -3 -9, 4 9, -4 3 -8, -7 -5 2 8")
  }; 
  static const ColourLines diag96[4] = { 
    ColourLines("1 4 5 8, -1 -2 3, -3 -9, 6 -8, -7 -4 2 9"), 
    ColourLines("1 4 5 8, -1 -2 -9, 3 9, -3 2 -4 -7, 6 -8"), 
    ColourLines("1 2 9, -1 -4 -7, 3 -2 4 5 8, -3 -9, 6 -8"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 9, 6 -8, 8 5 4 -2 -9")
  }; 
  static const ColourLines diag97[2] = { 
    ColourLines("1 5 8, -1 -2 -3 4, -4 -9, 6 -8, -7 3 9"), 
    ColourLines("1 5 8, -1 -2 -3 -9, 4 9, -4 3 -7, 6 -8")
  }; 
  static const ColourLines diag98[4] = { 
    ColourLines("1 4 6, -1 -2 3, -3 -9, -7 8, -8 -5 -4 2 9"), 
    ColourLines("1 4 6, -1 -2 -9, 3 9, -3 2 -4 -5 -8, -7 8"), 
    ColourLines("1 2 9, -1 -4 -5 -8, 3 -2 4 6, -3 -9, -7 8"), 
    ColourLines("1 2 -3, -1 -4 -5 -8, 3 9, 6 4 -2 -9, -7 8")
  }; 
  static const ColourLines diag99[2] = { 
    ColourLines("1 2 3 9, -1 -5 -8, 4 -3 6, -4 -9, -7 8"), 
    ColourLines("1 2 3 -4, -1 -5 -8, 4 9, 6 -3 -9, -7 8")
  }; 
  static const ColourLines diag100[8] = { 
    ColourLines("1 2 5 9, -1 -4 -7, 3 -2 4 6, -3 -5 -8, 8 -9"), 
    ColourLines("1 2 5 8, -1 -4 -7, 3 -2 4 6, -3 -5 -9, -8 9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -8, -7 -4 2 5 9, 8 -9"), 
    ColourLines("1 4 6, -1 -2 3, -3 -5 -9, -7 -4 2 5 8, -8 9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 5 9, 6 4 -2 -5 -8, 8 -9"), 
    ColourLines("1 2 -3, -1 -4 -7, 3 5 8, 6 4 -2 -5 -9, -8 9"), 
    ColourLines("1 4 6, -1 -2 -5 -8, 3 5 9, -3 2 -4 -7, 8 -9"), 
    ColourLines("1 4 6, -1 -2 -5 -9, 3 5 8, -3 2 -4 -7, -8 9")
  }; 
  static const ColourLines diag101[8] = { 
    ColourLines("1 2 4 6, -1 -5 -8, 3 -2 5 9, -3 -4 -7, 8 -9"), 
    ColourLines("1 2 4 6, -1 -5 -9, 3 -2 5 8, -3 -4 -7, -8 9"), 
    ColourLines("1 2 -3, -1 -5 -8, 3 4 6, -7 -4 -2 5 9, 8 -9"), 
    ColourLines("1 2 -3, -1 -5 -9, 3 4 6, -7 -4 -2 5 8, -8 9"), 
    ColourLines("1 5 9, -1 -2 3, -3 -4 -7, 6 4 2 -5 -8, 8 -9"), 
    ColourLines("1 5 8, -1 -2 3, -3 -4 -7, 6 4 2 -5 -9, -8 9"), 
    ColourLines("1 5 9, -1 -2 -4 -7, 3 4 6, -3 2 -5 -8, 8 -9"), 
    ColourLines("1 5 8, -1 -2 -4 -7, 3 4 6, -3 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag102[1] = { 
    ColourLines("1 4 8, -1 -2 3, -3 -5 -9, 6 -8, -7 9")
  }; 
  static const ColourLines diag103[1] = { 
    ColourLines("1 2 -3, -1 -5 -9, 3 4 8, 6 -8, -7 9")
  }; 
  static const ColourLines diag104[1] = { 
    ColourLines("1 4 9, -1 -2 3, -3 -5 -8, 6 -9, -7 8")
  }; 
  static const ColourLines diag105[1] = { 
    ColourLines("1 2 -3, -1 -5 -8, 3 4 9, 6 -9, -7 8")
  }; 

  static const int cv1at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv1at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv1at2[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv1at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv1at4[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv1at5[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv1at6[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv1at7[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv2at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv2at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv3at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv3at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv4at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv4at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv5at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv5at1[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv6at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv6at1[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv7at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv7at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv8at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv8at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv9at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv9at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv10at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv10at1[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv11at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv11at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv12at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv12at1[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv12at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv12at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv12at4[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv12at5[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv12at6[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv12at7[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv13at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv13at1[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv14at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv14at1[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv15at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv15at1[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv15at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv15at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv15at4[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv15at5[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv15at6[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv15at7[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv16at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv16at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv16at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv16at3[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv16at4[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv16at5[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv16at6[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv16at7[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv17at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv17at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv17at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv17at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv17at4[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv17at5[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv17at6[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv17at7[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv18at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv18at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv19at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv19at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv19at2[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv19at3[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv20at0[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv20at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv21at0[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv21at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv21at2[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv21at3[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv22at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv22at1[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv23at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv23at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv23at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv23at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv24at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv24at1[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv25at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv25at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv25at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv25at3[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv26at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv26at1[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv26at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv26at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv27at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv27at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv27at2[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv27at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv28at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv29at0[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv30at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv30at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv30at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv30at3[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv31at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv31at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv32at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv32at1[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv32at2[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv32at3[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv33at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv33at1[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv34at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv35at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv36at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv37at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv38at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv38at1[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv39at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv39at1[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv39at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv39at3[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv40at0[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv41at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv42at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv42at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv42at2[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv42at3[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv43at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv43at1[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv44at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv44at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv44at2[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv44at3[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv45at0[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv45at1[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv46at0[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv47at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv48at0[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv49at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv50at0[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv50at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv51at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv51at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv51at2[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv51at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv52at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv52at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv53at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv53at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv53at2[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv53at3[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv54at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv54at1[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv55at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv55at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv55at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv55at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv56at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv56at1[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv56at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv56at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv57at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv57at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv57at2[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv57at3[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv58at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv58at1[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv58at2[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv58at3[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv58at4[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv58at5[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv58at6[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv58at7[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv59at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv59at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv59at2[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv59at3[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv59at4[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv59at5[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv59at6[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv59at7[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv60at0[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv60at1[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv60at2[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv60at3[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv61at0[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv61at1[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv62at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv62at1[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv62at2[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv62at3[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv63at0[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv63at1[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv64at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv64at1[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv65at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv65at1[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv65at2[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv65at3[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv66at0[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv66at1[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv67at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv67at1[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv67at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv67at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv68at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv68at1[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv68at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv68at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv69at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv69at1[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv69at2[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv69at3[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv70at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv70at1[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv70at2[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv70at3[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv70at4[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv70at5[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv70at6[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv70at7[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv71at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv71at1[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv71at2[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv71at3[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv71at4[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv71at5[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv71at6[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv71at7[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv72at0[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv72at1[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv72at2[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv72at3[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv73at0[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv73at1[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv74at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv74at1[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv74at2[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv74at3[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv75at0[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv75at1[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv76at0[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv77at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv78at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv79at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv80at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv80at1[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv81at0[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv81at1[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv81at2[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv81at3[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv82at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv83at0[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv84at0[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv85at0[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv86at0[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv86at1[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv87at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv87at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv87at2[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv87at3[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv88at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv88at1[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv88at2[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv88at3[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv88at4[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv88at5[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv88at6[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv88at7[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv89at0[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv89at1[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv89at2[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv89at3[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv89at4[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv89at5[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv89at6[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv89at7[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv90at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv90at1[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv90at2[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv90at3[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv91at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv91at1[]  = { 1, 4, -1, 3, 0, 2, -999};
  static const int cv92at0[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv92at1[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv92at2[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv92at3[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv93at0[]  = { 1, 0, 3, -1, 4, 2, -999};
  static const int cv93at1[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv94at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv94at1[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv94at2[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv94at3[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv94at4[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv94at5[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv94at6[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv94at7[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv95at0[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv95at1[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv95at2[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv95at3[]  = { 1, 4, 0, 3, -1, 2, -999};
  static const int cv95at4[]  = { 1, -1, 3, 0, 4, 2, -999};
  static const int cv95at5[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv95at6[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv95at7[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv96at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv96at1[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv96at2[]  = { 1, 3, 0, 4, -1, 2, -999};
  static const int cv96at3[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv97at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv97at1[]  = { 1, 3, -1, 4, 0, 2, -999};
  static const int cv98at0[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv98at1[]  = { 1, -1, 4, 0, 3, 2, -999};
  static const int cv98at2[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv98at3[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv99at0[]  = { 1, 0, 4, -1, 3, 2, -999};
  static const int cv99at1[]  = { 1, 4, 0, -1, 3, 2, -999};
  static const int cv100at0[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv100at1[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv100at2[]  = { 1, -1, 0, 3, 4, 2, -999};
  static const int cv100at3[]  = { 1, -1, 0, 4, 3, 2, -999};
  static const int cv100at4[]  = { 1, 3, 4, 0, -1, 2, -999};
  static const int cv100at5[]  = { 1, 4, 3, 0, -1, 2, -999};
  static const int cv100at6[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv100at7[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv101at0[]  = { 1, -1, 3, 4, 0, 2, -999};
  static const int cv101at1[]  = { 1, -1, 4, 3, 0, 2, -999};
  static const int cv101at2[]  = { 1, 0, -1, 3, 4, 2, -999};
  static const int cv101at3[]  = { 1, 0, -1, 4, 3, 2, -999};
  static const int cv101at4[]  = { 1, 3, 4, -1, 0, 2, -999};
  static const int cv101at5[]  = { 1, 4, 3, -1, 0, 2, -999};
  static const int cv101at6[]  = { 1, 0, 3, 4, -1, 2, -999};
  static const int cv101at7[]  = { 1, 0, 4, 3, -1, 2, -999};
  static const int cv102at0[]  = { 1, 3, -1, 0, 4, 2, -999};
  static const int cv103at0[]  = { 1, 3, 0, -1, 4, 2, -999};
  static const int cv104at0[]  = { 1, 4, -1, 0, 3, 2, -999};
  static const int cv105at0[]  = { 1, 4, 0, -1, 3, 2, -999};

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
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int)),  &(diag3[1]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at2, sizeof(cv12at2)/sizeof(int)),  &(diag12[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at3, sizeof(cv12at3)/sizeof(int)),  &(diag12[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at4, sizeof(cv12at4)/sizeof(int)),  &(diag12[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at5, sizeof(cv12at5)/sizeof(int)),  &(diag12[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at6, sizeof(cv12at6)/sizeof(int)),  &(diag12[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at7, sizeof(cv12at7)/sizeof(int)),  &(diag12[7]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at2, sizeof(cv15at2)/sizeof(int)),  &(diag15[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at3, sizeof(cv15at3)/sizeof(int)),  &(diag15[3]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at4, sizeof(cv15at4)/sizeof(int)),  &(diag15[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at5, sizeof(cv15at5)/sizeof(int)),  &(diag15[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at6, sizeof(cv15at6)/sizeof(int)),  &(diag15[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at7, sizeof(cv15at7)/sizeof(int)),  &(diag15[7]) );
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
  } 
  else if( diag->id() == -19 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int)),  &(diag19[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at1, sizeof(cv19at1)/sizeof(int)),  &(diag19[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at2, sizeof(cv19at2)/sizeof(int)),  &(diag19[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at3, sizeof(cv19at3)/sizeof(int)),  &(diag19[3]) );
  } 
  else if( diag->id() == -20 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int)),  &(diag20[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at1, sizeof(cv20at1)/sizeof(int)),  &(diag20[1]) );
  } 
  else if( diag->id() == -21 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int)),  &(diag21[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at1, sizeof(cv21at1)/sizeof(int)),  &(diag21[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at2, sizeof(cv21at2)/sizeof(int)),  &(diag21[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at3, sizeof(cv21at3)/sizeof(int)),  &(diag21[3]) );
  } 
  else if( diag->id() == -22 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int)),  &(diag22[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at1, sizeof(cv22at1)/sizeof(int)),  &(diag22[1]) );
  } 
  else if( diag->id() == -23 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int)),  &(diag23[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at1, sizeof(cv23at1)/sizeof(int)),  &(diag23[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at2, sizeof(cv23at2)/sizeof(int)),  &(diag23[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at3, sizeof(cv23at3)/sizeof(int)),  &(diag23[3]) );
  } 
  else if( diag->id() == -24 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int)),  &(diag24[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at1, sizeof(cv24at1)/sizeof(int)),  &(diag24[1]) );
  } 
  else if( diag->id() == -25 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int)),  &(diag25[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at1, sizeof(cv25at1)/sizeof(int)),  &(diag25[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at2, sizeof(cv25at2)/sizeof(int)),  &(diag25[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at3, sizeof(cv25at3)/sizeof(int)),  &(diag25[3]) );
  } 
  else if( diag->id() == -26 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int)),  &(diag26[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at1, sizeof(cv26at1)/sizeof(int)),  &(diag26[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at2, sizeof(cv26at2)/sizeof(int)),  &(diag26[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at3, sizeof(cv26at3)/sizeof(int)),  &(diag26[3]) );
  } 
  else if( diag->id() == -27 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int)),  &(diag27[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at1, sizeof(cv27at1)/sizeof(int)),  &(diag27[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at2, sizeof(cv27at2)/sizeof(int)),  &(diag27[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at3, sizeof(cv27at3)/sizeof(int)),  &(diag27[3]) );
  } 
  else if( diag->id() == -28 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)),  &(diag28[0]) );
  } 
  else if( diag->id() == -29 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int)),  &(diag29[0]) );
  } 
  else if( diag->id() == -30 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int)),  &(diag30[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at1, sizeof(cv30at1)/sizeof(int)),  &(diag30[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at2, sizeof(cv30at2)/sizeof(int)),  &(diag30[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at3, sizeof(cv30at3)/sizeof(int)),  &(diag30[3]) );
  } 
  else if( diag->id() == -31 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int)),  &(diag31[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at1, sizeof(cv31at1)/sizeof(int)),  &(diag31[1]) );
  } 
  else if( diag->id() == -32 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int)),  &(diag32[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at1, sizeof(cv32at1)/sizeof(int)),  &(diag32[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at2, sizeof(cv32at2)/sizeof(int)),  &(diag32[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at3, sizeof(cv32at3)/sizeof(int)),  &(diag32[3]) );
  } 
  else if( diag->id() == -33 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int)),  &(diag33[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at1, sizeof(cv33at1)/sizeof(int)),  &(diag33[1]) );
  } 
  else if( diag->id() == -34 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv34at0, sizeof(cv34at0)/sizeof(int)),  &(diag34[0]) );
  } 
  else if( diag->id() == -35 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at0, sizeof(cv35at0)/sizeof(int)),  &(diag35[0]) );
  } 
  else if( diag->id() == -36 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int)),  &(diag36[0]) );
  } 
  else if( diag->id() == -37 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int)),  &(diag37[0]) );
  } 
  else if( diag->id() == -38 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int)),  &(diag38[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at1, sizeof(cv38at1)/sizeof(int)),  &(diag38[1]) );
  } 
  else if( diag->id() == -39 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int)),  &(diag39[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at1, sizeof(cv39at1)/sizeof(int)),  &(diag39[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at2, sizeof(cv39at2)/sizeof(int)),  &(diag39[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at3, sizeof(cv39at3)/sizeof(int)),  &(diag39[3]) );
  } 
  else if( diag->id() == -40 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int)),  &(diag40[0]) );
  } 
  else if( diag->id() == -41 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int)),  &(diag41[0]) );
  } 
  else if( diag->id() == -42 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int)),  &(diag42[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at1, sizeof(cv42at1)/sizeof(int)),  &(diag42[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at2, sizeof(cv42at2)/sizeof(int)),  &(diag42[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at3, sizeof(cv42at3)/sizeof(int)),  &(diag42[3]) );
  } 
  else if( diag->id() == -43 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at0, sizeof(cv43at0)/sizeof(int)),  &(diag43[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv43at1, sizeof(cv43at1)/sizeof(int)),  &(diag43[1]) );
  } 
  else if( diag->id() == -44 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at0, sizeof(cv44at0)/sizeof(int)),  &(diag44[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at1, sizeof(cv44at1)/sizeof(int)),  &(diag44[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at2, sizeof(cv44at2)/sizeof(int)),  &(diag44[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv44at3, sizeof(cv44at3)/sizeof(int)),  &(diag44[3]) );
  } 
  else if( diag->id() == -45 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at0, sizeof(cv45at0)/sizeof(int)),  &(diag45[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv45at1, sizeof(cv45at1)/sizeof(int)),  &(diag45[1]) );
  } 
  else if( diag->id() == -46 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv46at0, sizeof(cv46at0)/sizeof(int)),  &(diag46[0]) );
  } 
  else if( diag->id() == -47 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int)),  &(diag47[0]) );
  } 
  else if( diag->id() == -48 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int)),  &(diag48[0]) );
  } 
  else if( diag->id() == -49 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int)),  &(diag49[0]) );
  } 
  else if( diag->id() == -50 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at0, sizeof(cv50at0)/sizeof(int)),  &(diag50[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv50at1, sizeof(cv50at1)/sizeof(int)),  &(diag50[1]) );
  } 
  else if( diag->id() == -51 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at0, sizeof(cv51at0)/sizeof(int)),  &(diag51[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at1, sizeof(cv51at1)/sizeof(int)),  &(diag51[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at2, sizeof(cv51at2)/sizeof(int)),  &(diag51[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv51at3, sizeof(cv51at3)/sizeof(int)),  &(diag51[3]) );
  } 
  else if( diag->id() == -52 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at0, sizeof(cv52at0)/sizeof(int)),  &(diag52[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at1, sizeof(cv52at1)/sizeof(int)),  &(diag52[1]) );
  } 
  else if( diag->id() == -53 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int)),  &(diag53[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at1, sizeof(cv53at1)/sizeof(int)),  &(diag53[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at2, sizeof(cv53at2)/sizeof(int)),  &(diag53[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at3, sizeof(cv53at3)/sizeof(int)),  &(diag53[3]) );
  } 
  else if( diag->id() == -54 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int)),  &(diag54[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at1, sizeof(cv54at1)/sizeof(int)),  &(diag54[1]) );
  } 
  else if( diag->id() == -55 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int)),  &(diag55[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at1, sizeof(cv55at1)/sizeof(int)),  &(diag55[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at2, sizeof(cv55at2)/sizeof(int)),  &(diag55[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at3, sizeof(cv55at3)/sizeof(int)),  &(diag55[3]) );
  } 
  else if( diag->id() == -56 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int)),  &(diag56[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at1, sizeof(cv56at1)/sizeof(int)),  &(diag56[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at2, sizeof(cv56at2)/sizeof(int)),  &(diag56[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at3, sizeof(cv56at3)/sizeof(int)),  &(diag56[3]) );
  } 
  else if( diag->id() == -57 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int)),  &(diag57[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at1, sizeof(cv57at1)/sizeof(int)),  &(diag57[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at2, sizeof(cv57at2)/sizeof(int)),  &(diag57[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at3, sizeof(cv57at3)/sizeof(int)),  &(diag57[3]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv64at1, sizeof(cv64at1)/sizeof(int)),  &(diag64[1]) );
  } 
  else if( diag->id() == -65 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at0, sizeof(cv65at0)/sizeof(int)),  &(diag65[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at1, sizeof(cv65at1)/sizeof(int)),  &(diag65[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at2, sizeof(cv65at2)/sizeof(int)),  &(diag65[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv65at3, sizeof(cv65at3)/sizeof(int)),  &(diag65[3]) );
  } 
  else if( diag->id() == -66 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at0, sizeof(cv66at0)/sizeof(int)),  &(diag66[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv66at1, sizeof(cv66at1)/sizeof(int)),  &(diag66[1]) );
  } 
  else if( diag->id() == -67 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at0, sizeof(cv67at0)/sizeof(int)),  &(diag67[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at1, sizeof(cv67at1)/sizeof(int)),  &(diag67[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at2, sizeof(cv67at2)/sizeof(int)),  &(diag67[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at3, sizeof(cv67at3)/sizeof(int)),  &(diag67[3]) );
  } 
  else if( diag->id() == -68 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int)),  &(diag68[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at1, sizeof(cv68at1)/sizeof(int)),  &(diag68[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at2, sizeof(cv68at2)/sizeof(int)),  &(diag68[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at3, sizeof(cv68at3)/sizeof(int)),  &(diag68[3]) );
  } 
  else if( diag->id() == -69 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int)),  &(diag69[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at1, sizeof(cv69at1)/sizeof(int)),  &(diag69[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at2, sizeof(cv69at2)/sizeof(int)),  &(diag69[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at3, sizeof(cv69at3)/sizeof(int)),  &(diag69[3]) );
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
  } 
  else if( diag->id() == -77 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv77at0, sizeof(cv77at0)/sizeof(int)),  &(diag77[0]) );
  } 
  else if( diag->id() == -78 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv78at0, sizeof(cv78at0)/sizeof(int)),  &(diag78[0]) );
  } 
  else if( diag->id() == -79 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv79at0, sizeof(cv79at0)/sizeof(int)),  &(diag79[0]) );
  } 
  else if( diag->id() == -80 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at0, sizeof(cv80at0)/sizeof(int)),  &(diag80[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv80at1, sizeof(cv80at1)/sizeof(int)),  &(diag80[1]) );
  } 
  else if( diag->id() == -81 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at0, sizeof(cv81at0)/sizeof(int)),  &(diag81[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at1, sizeof(cv81at1)/sizeof(int)),  &(diag81[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at2, sizeof(cv81at2)/sizeof(int)),  &(diag81[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv81at3, sizeof(cv81at3)/sizeof(int)),  &(diag81[3]) );
  } 
  else if( diag->id() == -82 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv82at0, sizeof(cv82at0)/sizeof(int)),  &(diag82[0]) );
  } 
  else if( diag->id() == -83 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv83at0, sizeof(cv83at0)/sizeof(int)),  &(diag83[0]) );
  } 
  else if( diag->id() == -84 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv84at0, sizeof(cv84at0)/sizeof(int)),  &(diag84[0]) );
  } 
  else if( diag->id() == -85 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv85at0, sizeof(cv85at0)/sizeof(int)),  &(diag85[0]) );
  } 
  else if( diag->id() == -86 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at0, sizeof(cv86at0)/sizeof(int)),  &(diag86[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv86at1, sizeof(cv86at1)/sizeof(int)),  &(diag86[1]) );
  } 
  else if( diag->id() == -87 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at0, sizeof(cv87at0)/sizeof(int)),  &(diag87[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at1, sizeof(cv87at1)/sizeof(int)),  &(diag87[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at2, sizeof(cv87at2)/sizeof(int)),  &(diag87[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv87at3, sizeof(cv87at3)/sizeof(int)),  &(diag87[3]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at4, sizeof(cv94at4)/sizeof(int)),  &(diag94[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at5, sizeof(cv94at5)/sizeof(int)),  &(diag94[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at6, sizeof(cv94at6)/sizeof(int)),  &(diag94[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv94at7, sizeof(cv94at7)/sizeof(int)),  &(diag94[7]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at4, sizeof(cv100at4)/sizeof(int)),  &(diag100[4]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at5, sizeof(cv100at5)/sizeof(int)),  &(diag100[5]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at6, sizeof(cv100at6)/sizeof(int)),  &(diag100[6]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv100at7, sizeof(cv100at7)/sizeof(int)),  &(diag100[7]) );
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
  } 
  else if( diag->id() == -102 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv102at0, sizeof(cv102at0)/sizeof(int)),  &(diag102[0]) );
  } 
  else if( diag->id() == -103 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv103at0, sizeof(cv103at0)/sizeof(int)),  &(diag103[0]) );
  } 
  else if( diag->id() == -104 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv104at0, sizeof(cv104at0)/sizeof(int)),  &(diag104[0]) );
  } 
  else if( diag->id() == -105 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv105at0, sizeof(cv105at0)/sizeof(int)),  &(diag105[0]) );
  } 
  return sel;
}
