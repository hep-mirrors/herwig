// -*- C++ -*-
//
// NLOJetMEqbg2kkbqbg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEqbg2kkbqbg class.
//

#include "NLOJetMEqbg2kkbqbg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEqbg2kkbqbg::NLOJetMEqbg2kkbqbg() {}

NLOJetMEqbg2kkbqbg::~NLOJetMEqbg2kkbqbg() {}

IBPtr NLOJetMEqbg2kkbqbg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEqbg2kkbqbg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEqbg2kkbqbg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEqbg2kkbqbg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEqbg2kkbqbg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEqbg2kkbqbg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEqbg2kkbqbg("Herwig::NLOJetMEqbg2kkbqbg", "HwMatchboxNLOJet.so");

void NLOJetMEqbg2kkbqbg::Init() {

  static ClassDocumentation<NLOJetMEqbg2kkbqbg> documentation
    ("NLOJetMEqbg2kkbqbg");

}


void NLOJetMEqbg2kkbqbg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr k = quark[xi1]; 
    tcPDPtr kb = antiquark[xi1]; 
    for(unsigned xi2 = 0; xi2 != quark.size(); xi2++){ 
      tcPDPtr q = quark[xi2]; 
      tcPDPtr qb = antiquark[xi2]; 
      addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 3, qb, 4, k, 4, kb, 5, qb, 5, g, -1)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, k, g, 2, kb, 3, k, 5, kb, 1, qb, 5, g, -2)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, kb, g, 2, k, 5, k, 3, kb, 1, qb, 5, g, -3)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, g, g, 2, g, 5, k, 5, kb, 1, qb, 3, g, -4)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, g, 2, g, 5, k, 5, kb, 3, qb, 1, g, -5)));
      addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 4, g, 5, k, 5, kb, 3, qb, 4, g, -6)));
      addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, qb, 4, g, 5, k, 5, kb, 4, qb, 3, g, -7)));
      addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 4, k, 5, k, 4, kb, 3, qb, 5, g, -8)));
      addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 4, kb, 4, k, 5, kb, 3, qb, 5, g, -9)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, g, k, k, g, 4, k, 2, kb, 1, qb, 3, g, -10)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, g, g, k, g, 4, k, 3, kb, 1, qb, 2, g, -11)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, g, kb, kb, g, 2, k, 4, kb, 1, qb, 3, g, -12)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, g, g, kb, g, 3, k, 4, kb, 1, qb, 2, g, -13)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, g, kb, g, g, 2, k, 3, kb, 1, qb, 4, g, -14)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, g, k, g, g, 3, k, 2, kb, 1, qb, 4, g, -15)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 4, g, 5, k, 5, kb, 1, qb, 4, g, -16)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, g, g, 3, g, 5, k, 5, kb, 1, qb, 2, g, -17)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 4, k, 5, k, 4, kb, 1, qb, 5, g, -18)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, k, g, 3, k, 5, k, 2, kb, 1, qb, 5, g, -19)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 4, kb, 4, k, 5, kb, 1, qb, 5, g, -20)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, kb, g, 3, kb, 2, k, 5, kb, 1, qb, 5, g, -21)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, qb, g, k, g, 4, k, 3, kb, 2, qb, 1, g, -22)));
      addSafe(new_ptr((Tree2toNDiagram(5), qb, qb, g, kb, g, 3, k, 4, kb, 2, qb, 1, g, -23)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 2, qb, 4, g, 5, k, 5, kb, 4, qb, 1, g, -24)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, g, g, 3, g, 5, k, 5, kb, 2, qb, 1, g, -25)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, k, g, 1, qb, 3, k, 2, kb, 5, qb, 5, g, -26)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, g, kb, g, 1, qb, 2, k, 3, kb, 5, qb, 5, g, -27)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 4, g, 5, k, 5, kb, 2, qb, 4, g, -28)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, g, 1, g, 5, k, 5, kb, 3, qb, 2, g, -29)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 4, k, 5, k, 4, kb, 2, qb, 5, g, -30)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 4, kb, 4, k, 5, kb, 2, qb, 5, g, -31)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 1, qb, 4, g, 5, k, 5, kb, 4, qb, 2, g, -32)));
      addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, g, g, 1, g, 5, k, 5, kb, 2, qb, 3, g, -33)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 2, qb, 4, k, 4, kb, 5, qb, 5, g, -34)));
      addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 1, qb, 4, k, 4, kb, 5, qb, 5, g, -35)));
      if( ( xi1 ==  xi2 ) ){ 
        addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 3, qb, 4, q, 5, qb, 4, qb, 5, g, -36)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, q, g, 2, qb, 3, q, 1, qb, 5, qb, 5, g, -37)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, qb, g, 2, q, 5, q, 1, qb, 3, qb, 5, g, -38)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, g, g, 2, g, 5, q, 1, qb, 5, qb, 3, g, -39)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, g, 2, g, 5, q, 3, qb, 5, qb, 1, g, -40)));
        addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 4, g, 5, q, 3, qb, 5, qb, 4, g, -41)));
        addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, qb, 4, g, 5, q, 4, qb, 5, qb, 3, g, -42)));
        addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 4, q, 5, q, 3, qb, 4, qb, 5, g, -43)));
        addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, g, 4, qb, 4, q, 3, qb, 5, qb, 5, g, -44)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, g, q, q, g, 4, q, 1, qb, 2, qb, 3, g, -45)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, g, g, q, g, 4, q, 1, qb, 3, qb, 2, g, -46)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, g, qb, qb, g, 2, q, 1, qb, 4, qb, 3, g, -47)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, g, g, qb, g, 3, q, 1, qb, 4, qb, 2, g, -48)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, g, qb, g, g, 2, q, 1, qb, 3, qb, 4, g, -49)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, g, q, g, g, 3, q, 1, qb, 2, qb, 4, g, -50)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 4, g, 5, q, 1, qb, 5, qb, 4, g, -51)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, g, g, 3, g, 5, q, 1, qb, 5, qb, 2, g, -52)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 4, q, 5, q, 1, qb, 4, qb, 5, g, -53)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, q, g, 3, q, 5, q, 1, qb, 2, qb, 5, g, -54)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 4, qb, 4, q, 1, qb, 5, qb, 5, g, -55)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, qb, g, 3, qb, 2, q, 1, qb, 5, qb, 5, g, -56)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, qb, g, q, g, 4, q, 2, qb, 3, qb, 1, g, -57)));
        addSafe(new_ptr((Tree2toNDiagram(5), qb, qb, g, qb, g, 3, q, 2, qb, 4, qb, 1, g, -58)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 2, qb, 4, g, 5, q, 4, qb, 5, qb, 1, g, -59)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, g, g, 3, g, 5, q, 2, qb, 5, qb, 1, g, -60)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, q, g, 1, qb, 3, q, 5, qb, 2, qb, 5, g, -61)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 4, g, 5, q, 2, qb, 5, qb, 4, g, -62)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, g, 1, g, 5, q, 3, qb, 5, qb, 2, g, -63)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 4, q, 5, q, 2, qb, 4, qb, 5, g, -64)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 4, qb, 4, q, 2, qb, 5, qb, 5, g, -65)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, g, qb, g, 1, qb, 2, q, 5, qb, 3, qb, 5, g, -66)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 1, qb, 4, g, 5, q, 4, qb, 5, qb, 2, g, -67)));
        addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, g, g, 1, g, 5, q, 2, qb, 5, qb, 3, g, -68)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, g, 2, qb, 4, q, 5, qb, 4, qb, 5, g, -69)));
        addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 2, g, 1, qb, 4, q, 5, qb, 4, qb, 5, g, -70)));
      }
    }  
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEqbg2kkbqbg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv2at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv3at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv4at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv4at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv4at2[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv4at3[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv5at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv6at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv6at1[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv7at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv8at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv9at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv10at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv11at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv11at1[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv12at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv13at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv13at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv14at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv14at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv15at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv15at1[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv16at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv16at1[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv16at2[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv16at3[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv17at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv17at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv17at2[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv17at3[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv18at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv18at1[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv19at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv20at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv20at1[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv21at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv22at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv23at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv24at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv25at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv25at1[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv26at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv27at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv28at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv28at1[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv29at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv30at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv31at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv32at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv32at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv33at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv33at1[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv34at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv35at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv35at1[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv36at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv37at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv38at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv39at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv39at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv39at2[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv39at3[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv40at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv41at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv41at1[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv42at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv43at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv44at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv45at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv46at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv46at1[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv47at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv48at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv48at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv49at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv49at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv50at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv50at1[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv51at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv51at1[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv51at2[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv51at3[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv52at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv52at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv52at2[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv52at3[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv53at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv53at1[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv54at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv55at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv55at1[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv56at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv57at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv58at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv59at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv60at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv60at1[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv61at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv62at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv62at1[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv63at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv64at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv65at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv66at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv67at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv67at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv68at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv68at1[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv69at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv70at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv70at1[]  = { -1, 3, -999, 1, 0, 4, 2, -999};

  Selector <MEBase::DiagramIndex> sel;
  for(MEBase::DiagramIndex i=0; i < diags.size(); ++i){
    if( diags[i]->id() == -1 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -4 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at2, sizeof(cv4at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv4at3, sizeof(cv4at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -5 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -6 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -7 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -8 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -9 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -10 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -11 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -12 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv16at3, sizeof(cv16at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -17 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv17at0, sizeof(cv17at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at1, sizeof(cv17at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at2, sizeof(cv17at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv17at3, sizeof(cv17at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -18 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv18at0, sizeof(cv18at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv18at1, sizeof(cv18at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -19 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -20 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv20at1, sizeof(cv20at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -21 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -22 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -23 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -24 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -25 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv25at1, sizeof(cv25at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -26 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -27 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -28 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv28at1, sizeof(cv28at1)/sizeof(int)), i);   
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
          nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv32at1, sizeof(cv32at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -33 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv33at1, sizeof(cv33at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -34 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv34at0, sizeof(cv34at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -35 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv35at0, sizeof(cv35at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv35at1, sizeof(cv35at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -36 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -37 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -38 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int)), i);   
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
          nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv41at1, sizeof(cv41at1)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv46at1, sizeof(cv46at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -47 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -48 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv48at1, sizeof(cv48at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -49 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv49at1, sizeof(cv49at1)/sizeof(int)), i);   
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
        + nloJetAmplitude()->colourOrdered2(cv52at1, sizeof(cv52at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at2, sizeof(cv52at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv52at3, sizeof(cv52at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -53 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv53at1, sizeof(cv53at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -54 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -55 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv55at1, sizeof(cv55at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -56 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -57 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -58 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv58at0, sizeof(cv58at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -59 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv59at0, sizeof(cv59at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -60 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv60at0, sizeof(cv60at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv60at1, sizeof(cv60at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -61 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv61at0, sizeof(cv61at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -62 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv62at0, sizeof(cv62at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv62at1, sizeof(cv62at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -63 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv63at0, sizeof(cv63at0)/sizeof(int)), i);   
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
          nloJetAmplitude()->colourOrdered2(cv67at0, sizeof(cv67at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv67at1, sizeof(cv67at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -68 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv68at1, sizeof(cv68at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -69 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -70 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv70at0, sizeof(cv70at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv70at1, sizeof(cv70at1)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEqbg2kkbqbg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEqbg2kkbqbg

  static const ColourLines diag1[1] = { 
    ColourLines("-1 2, -2 -3 -4 -7, 6 4 -5 -9, -8 9")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("-1 -2 -5 -9, 4 6, -4 3 2 -8, -7 9")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("-1 -2 -3 4, -4 -7, 6 -9, -8 2 5 9")
  }; 
  static const ColourLines diag4[4] = { 
    ColourLines("-1 -2 -3 4, -4 -9, 6 5 2 -8, -7 -5 3 9"), 
    ColourLines("-1 -2 -3 -9, 4 9, -4 3 -5 -7, 6 5 2 -8"), 
    ColourLines("-1 -2 -5 -7, 4 -3 5 6, -4 -9, -8 2 3 9"), 
    ColourLines("-1 -2 -5 -7, 4 9, -4 3 2 -8, 6 5 -3 -9")
  }; 
  static const ColourLines diag5[1] = { 
    ColourLines("-1 -9, 4 -3 5 6, -4 -8, -7 -5 -2 9")
  }; 
  static const ColourLines diag6[2] = { 
    ColourLines("-1 2, -2 -3 -4 -9, 6 5 4 -8, -7 -5 9"), 
    ColourLines("-1 2, -2 -3 -4 -5 -7, 6 5 -9, -8 4 9")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("-1 2, -2 -3 -9, 6 5 -8, -7 -5 -4 9")
  }; 
  static const ColourLines diag8[1] = { 
    ColourLines("-1 2, -2 -3 -4 -7, 6 -9, -8 4 5 9")
  }; 
  static const ColourLines diag9[1] = { 
    ColourLines("-1 2, -2 -3 -4 -5 -9, 6 4 -8, -7 9")
  }; 
  static const ColourLines diag10[1] = { 
    ColourLines("-1 -2 -7, 5 6, -5 4 -9, -8 2 3 9")
  }; 
  static const ColourLines diag11[2] = { 
    ColourLines("-1 -2 -3 -7, 5 6, -5 4 3 -9, -8 2 9"), 
    ColourLines("-1 -2 -9, 5 6, -5 4 3 2 -8, -7 -3 9")
  }; 
  static const ColourLines diag12[1] = { 
    ColourLines("-1 -2 -3 -9, 5 -4 9, -5 -7, 6 2 -8")
  }; 
  static const ColourLines diag13[2] = { 
    ColourLines("-1 -2 -3 -4 5, -5 -7, 6 3 -9, -8 2 9"), 
    ColourLines("-1 -2 -9, 5 -4 -3 9, -5 -7, 6 3 2 -8")
  }; 
  static const ColourLines diag14[2] = { 
    ColourLines("-1 -2 -3 -4 5, -5 -9, 6 2 -8, -7 4 9"), 
    ColourLines("-1 -2 -3 -4 -9, 5 9, -5 4 -7, 6 2 -8")
  }; 
  static const ColourLines diag15[2] = { 
    ColourLines("-1 -2 -7, 5 -4 6, -5 -9, -8 2 3 4 9"), 
    ColourLines("-1 -2 -7, 5 9, -5 4 3 2 -8, 6 -4 -9")
  }; 
  static const ColourLines diag16[4] = { 
    ColourLines("-1 -2 -4 -9, 3 4 5 6, -3 2 -8, -7 -5 9"), 
    ColourLines("-1 -2 3, -3 -4 -9, 6 5 4 2 -8, -7 -5 9"), 
    ColourLines("-1 -2 -4 -5 -7, 3 4 9, -3 2 -8, 6 5 -9"), 
    ColourLines("-1 -2 3, -3 -4 -5 -7, 6 5 -9, -8 2 4 9")
  }; 
  static const ColourLines diag17[4] = { 
    ColourLines("-1 -2 -3 4, -4 -5 -7, 6 5 3 -9, -8 2 9"), 
    ColourLines("-1 -2 -9, 4 -3 9, -4 -5 -7, 6 5 3 2 -8"), 
    ColourLines("-1 -2 -3 -5 -7, 4 5 6, -4 3 -9, -8 2 9"), 
    ColourLines("-1 -2 -9, 4 5 6, -4 3 2 -8, -7 -5 -3 9")
  }; 
  static const ColourLines diag18[2] = { 
    ColourLines("-1 -2 -4 -7, 3 4 5 9, -3 2 -8, 6 -9"), 
    ColourLines("-1 -2 3, -3 -4 -7, 6 -9, -8 2 4 5 9")
  }; 
  static const ColourLines diag19[1] = { 
    ColourLines("-1 -2 -7, 4 5 9, -4 3 2 -8, 6 -9")
  }; 
  static const ColourLines diag20[2] = { 
    ColourLines("-1 -2 -4 -5 -9, 3 4 6, -3 2 -8, -7 9"), 
    ColourLines("-1 -2 3, -3 -4 -5 -9, 6 4 2 -8, -7 9")
  }; 
  static const ColourLines diag21[1] = { 
    ColourLines("-1 -2 -3 4, -4 -5 -9, 6 2 -8, -7 9")
  }; 
  static const ColourLines diag22[1] = { 
    ColourLines("-1 -9, 5 6, -5 4 3 -8, -7 -3 -2 9")
  }; 
  static const ColourLines diag23[1] = { 
    ColourLines("-1 -9, 5 -4 -3 -2 9, -5 -7, 6 3 -8")
  }; 
  static const ColourLines diag24[1] = { 
    ColourLines("-1 -9, 3 -2 9, -3 -4 -5 -7, 6 5 -8")
  }; 
  static const ColourLines diag25[2] = { 
    ColourLines("-1 -9, 4 -3 -2 9, -4 -5 -7, 6 5 3 -8"), 
    ColourLines("-1 -9, 4 5 6, -4 3 -8, -7 -5 -3 -2 9")
  }; 
  static const ColourLines diag26[1] = { 
    ColourLines("-1 -2 -7, 4 6, -4 3 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag27[1] = { 
    ColourLines("-1 -2 -3 4, -4 -7, 6 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag28[2] = { 
    ColourLines("-1 -4 -9, 3 -2 4 5 6, -3 -8, -7 -5 9"), 
    ColourLines("-1 -4 -5 -7, 3 -2 4 9, -3 -8, 6 5 -9")
  }; 
  static const ColourLines diag29[1] = { 
    ColourLines("-1 -5 -7, 4 -3 9, -4 -8, 6 5 -2 -9")
  }; 
  static const ColourLines diag30[1] = { 
    ColourLines("-1 -4 -7, 3 -2 4 5 9, -3 -8, 6 -9")
  }; 
  static const ColourLines diag31[1] = { 
    ColourLines("-1 -4 -5 -9, 3 -2 4 6, -3 -8, -7 9")
  }; 
  static const ColourLines diag32[2] = { 
    ColourLines("-1 -2 3, -3 -9, 6 5 -8, -7 -5 -4 2 9"), 
    ColourLines("-1 -2 -9, 3 9, -3 2 -4 -5 -7, 6 5 -8")
  }; 
  static const ColourLines diag33[2] = { 
    ColourLines("-1 -5 -7, 4 -3 -2 5 6, -4 -9, -8 3 9"), 
    ColourLines("-1 -5 -7, 4 9, -4 3 -8, 6 5 -2 -3 -9")
  }; 
  static const ColourLines diag34[1] = { 
    ColourLines("-1 -4 -7, 3 -2 4 6, -3 -5 -9, -8 9")
  }; 
  static const ColourLines diag35[2] = { 
    ColourLines("-1 -2 3, -3 -4 -7, 6 4 2 -5 -9, -8 9"), 
    ColourLines("-1 -2 -4 -7, 3 4 6, -3 2 -5 -9, -8 9")
  }; 
  static const ColourLines diag36[1] = { 
    ColourLines("-1 2, -2 -3 -4 -8, 6 4 -5 -9, -7 9")
  }; 
  static const ColourLines diag37[1] = { 
    ColourLines("-1 -2 -5 -9, 4 6, -4 3 2 -7, -8 9")
  }; 
  static const ColourLines diag38[1] = { 
    ColourLines("-1 -2 -3 4, -4 -8, 6 -9, -7 2 5 9")
  }; 
  static const ColourLines diag39[4] = { 
    ColourLines("-1 -2 -3 4, -4 -9, 6 5 2 -7, -8 -5 3 9"), 
    ColourLines("-1 -2 -3 -9, 4 9, -4 3 -5 -8, 6 5 2 -7"), 
    ColourLines("-1 -2 -5 -8, 4 -3 5 6, -4 -9, -7 2 3 9"), 
    ColourLines("-1 -2 -5 -8, 4 9, -4 3 2 -7, 6 5 -3 -9")
  }; 
  static const ColourLines diag40[1] = { 
    ColourLines("-1 -9, 4 -3 5 6, -4 -7, -8 -5 -2 9")
  }; 
  static const ColourLines diag41[2] = { 
    ColourLines("-1 2, -2 -3 -4 -9, 6 5 4 -7, -8 -5 9"), 
    ColourLines("-1 2, -2 -3 -4 -5 -8, 6 5 -9, -7 4 9")
  }; 
  static const ColourLines diag42[1] = { 
    ColourLines("-1 2, -2 -3 -9, 6 5 -7, -8 -5 -4 9")
  }; 
  static const ColourLines diag43[1] = { 
    ColourLines("-1 2, -2 -3 -4 -8, 6 -9, -7 4 5 9")
  }; 
  static const ColourLines diag44[1] = { 
    ColourLines("-1 2, -2 -3 -4 -5 -9, 6 4 -7, -8 9")
  }; 
  static const ColourLines diag45[1] = { 
    ColourLines("-1 -2 -8, 5 6, -5 4 -9, -7 2 3 9")
  }; 
  static const ColourLines diag46[2] = { 
    ColourLines("-1 -2 -3 -8, 5 6, -5 4 3 -9, -7 2 9"), 
    ColourLines("-1 -2 -9, 5 6, -5 4 3 2 -7, -8 -3 9")
  }; 
  static const ColourLines diag47[1] = { 
    ColourLines("-1 -2 -3 -9, 5 -4 9, -5 -8, 6 2 -7")
  }; 
  static const ColourLines diag48[2] = { 
    ColourLines("-1 -2 -3 -4 5, -5 -8, 6 3 -9, -7 2 9"), 
    ColourLines("-1 -2 -9, 5 -4 -3 9, -5 -8, 6 3 2 -7")
  }; 
  static const ColourLines diag49[2] = { 
    ColourLines("-1 -2 -3 -4 5, -5 -9, 6 2 -7, -8 4 9"), 
    ColourLines("-1 -2 -3 -4 -9, 5 9, -5 4 -8, 6 2 -7")
  }; 
  static const ColourLines diag50[2] = { 
    ColourLines("-1 -2 -8, 5 -4 6, -5 -9, -7 2 3 4 9"), 
    ColourLines("-1 -2 -8, 5 9, -5 4 3 2 -7, 6 -4 -9")
  }; 
  static const ColourLines diag51[4] = { 
    ColourLines("-1 -2 -4 -9, 3 4 5 6, -3 2 -7, -8 -5 9"), 
    ColourLines("-1 -2 3, -3 -4 -9, 6 5 4 2 -7, -8 -5 9"), 
    ColourLines("-1 -2 -4 -5 -8, 3 4 9, -3 2 -7, 6 5 -9"), 
    ColourLines("-1 -2 3, -3 -4 -5 -8, 6 5 -9, -7 2 4 9")
  }; 
  static const ColourLines diag52[4] = { 
    ColourLines("-1 -2 -3 4, -4 -5 -8, 6 5 3 -9, -7 2 9"), 
    ColourLines("-1 -2 -9, 4 -3 9, -4 -5 -8, 6 5 3 2 -7"), 
    ColourLines("-1 -2 -3 -5 -8, 4 5 6, -4 3 -9, -7 2 9"), 
    ColourLines("-1 -2 -9, 4 5 6, -4 3 2 -7, -8 -5 -3 9")
  }; 
  static const ColourLines diag53[2] = { 
    ColourLines("-1 -2 -4 -8, 3 4 5 9, -3 2 -7, 6 -9"), 
    ColourLines("-1 -2 3, -3 -4 -8, 6 -9, -7 2 4 5 9")
  }; 
  static const ColourLines diag54[1] = { 
    ColourLines("-1 -2 -8, 4 5 9, -4 3 2 -7, 6 -9")
  }; 
  static const ColourLines diag55[2] = { 
    ColourLines("-1 -2 -4 -5 -9, 3 4 6, -3 2 -7, -8 9"), 
    ColourLines("-1 -2 3, -3 -4 -5 -9, 6 4 2 -7, -8 9")
  }; 
  static const ColourLines diag56[1] = { 
    ColourLines("-1 -2 -3 4, -4 -5 -9, 6 2 -7, -8 9")
  }; 
  static const ColourLines diag57[1] = { 
    ColourLines("-1 -9, 5 6, -5 4 3 -7, -8 -3 -2 9")
  }; 
  static const ColourLines diag58[1] = { 
    ColourLines("-1 -9, 5 -4 -3 -2 9, -5 -8, 6 3 -7")
  }; 
  static const ColourLines diag59[1] = { 
    ColourLines("-1 -9, 3 -2 9, -3 -4 -5 -8, 6 5 -7")
  }; 
  static const ColourLines diag60[2] = { 
    ColourLines("-1 -9, 4 -3 -2 9, -4 -5 -8, 6 5 3 -7"), 
    ColourLines("-1 -9, 4 5 6, -4 3 -7, -8 -5 -3 -2 9")
  }; 
  static const ColourLines diag61[1] = { 
    ColourLines("-1 -2 -8, 4 6, -4 3 2 -5 -9, -7 9")
  }; 
  static const ColourLines diag62[2] = { 
    ColourLines("-1 -4 -9, 3 -2 4 5 6, -3 -7, -8 -5 9"), 
    ColourLines("-1 -4 -5 -8, 3 -2 4 9, -3 -7, 6 5 -9")
  }; 
  static const ColourLines diag63[1] = { 
    ColourLines("-1 -5 -8, 4 -3 9, -4 -7, 6 5 -2 -9")
  }; 
  static const ColourLines diag64[1] = { 
    ColourLines("-1 -4 -8, 3 -2 4 5 9, -3 -7, 6 -9")
  }; 
  static const ColourLines diag65[1] = { 
    ColourLines("-1 -4 -5 -9, 3 -2 4 6, -3 -7, -8 9")
  }; 
  static const ColourLines diag66[1] = { 
    ColourLines("-1 -2 -3 4, -4 -8, 6 2 -5 -9, -7 9")
  }; 
  static const ColourLines diag67[2] = { 
    ColourLines("-1 -2 3, -3 -9, 6 5 -7, -8 -5 -4 2 9"), 
    ColourLines("-1 -2 -9, 3 9, -3 2 -4 -5 -8, 6 5 -7")
  }; 
  static const ColourLines diag68[2] = { 
    ColourLines("-1 -5 -8, 4 -3 -2 5 6, -4 -9, -7 3 9"), 
    ColourLines("-1 -5 -8, 4 9, -4 3 -7, 6 5 -2 -3 -9")
  }; 
  static const ColourLines diag69[1] = { 
    ColourLines("-1 -4 -8, 3 -2 4 6, -3 -5 -9, -7 9")
  }; 
  static const ColourLines diag70[2] = { 
    ColourLines("-1 -2 3, -3 -4 -8, 6 4 2 -5 -9, -7 9"), 
    ColourLines("-1 -2 -4 -8, 3 4 6, -3 2 -5 -9, -7 9")
  }; 

  static const int cv1at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv2at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv3at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv4at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv4at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv4at2[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv4at3[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv5at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv6at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv6at1[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv7at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv8at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv9at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv10at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv11at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv11at1[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv12at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv13at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv13at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv14at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv14at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv15at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv15at1[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv16at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv16at1[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv16at2[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv16at3[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv17at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv17at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv17at2[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv17at3[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv18at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv18at1[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv19at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv20at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv20at1[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv21at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv22at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv23at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv24at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv25at0[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv25at1[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv26at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv27at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv28at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv28at1[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv29at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv30at0[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv31at0[]  = { -1, 4, 2, -999, 1, 0, 3, -999};
  static const int cv32at0[]  = { -1, 0, 4, 2, -999, 1, 3, -999};
  static const int cv32at1[]  = { -1, 4, 0, 2, -999, 1, 3, -999};
  static const int cv33at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv33at1[]  = { -1, 2, -999, 1, 4, 0, 3, -999};
  static const int cv34at0[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv35at0[]  = { -1, 0, 2, -999, 1, 4, 3, -999};
  static const int cv35at1[]  = { -1, 2, -999, 1, 0, 4, 3, -999};
  static const int cv36at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv37at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv38at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv39at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv39at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv39at2[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv39at3[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv40at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv41at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv41at1[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv42at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv43at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv44at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv45at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv46at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv46at1[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv47at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv48at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv48at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv49at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv49at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv50at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv50at1[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv51at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv51at1[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv51at2[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv51at3[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv52at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv52at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv52at2[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv52at3[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv53at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv53at1[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv54at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv55at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv55at1[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv56at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv57at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv58at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv59at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv60at0[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv60at1[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv61at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv62at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv62at1[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv63at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv64at0[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv65at0[]  = { -1, 4, 3, -999, 1, 0, 2, -999};
  static const int cv66at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv67at0[]  = { -1, 0, 4, 3, -999, 1, 2, -999};
  static const int cv67at1[]  = { -1, 4, 0, 3, -999, 1, 2, -999};
  static const int cv68at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv68at1[]  = { -1, 3, -999, 1, 4, 0, 2, -999};
  static const int cv69at0[]  = { -1, 3, -999, 1, 0, 4, 2, -999};
  static const int cv70at0[]  = { -1, 0, 3, -999, 1, 4, 2, -999};
  static const int cv70at1[]  = { -1, 3, -999, 1, 0, 4, 2, -999};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
  } 
  else if( diag->id() == -4 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int)),  &(diag4[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at1, sizeof(cv4at1)/sizeof(int)),  &(diag4[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at2, sizeof(cv4at2)/sizeof(int)),  &(diag4[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at3, sizeof(cv4at3)/sizeof(int)),  &(diag4[3]) );
  } 
  else if( diag->id() == -5 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int)),  &(diag5[0]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int)),  &(diag6[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at1, sizeof(cv6at1)/sizeof(int)),  &(diag6[1]) );
  } 
  else if( diag->id() == -7 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int)),  &(diag7[0]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int)),  &(diag8[0]) );
  } 
  else if( diag->id() == -9 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int)),  &(diag9[0]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)),  &(diag10[0]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)),  &(diag11[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at1, sizeof(cv11at1)/sizeof(int)),  &(diag11[1]) );
  } 
  else if( diag->id() == -12 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int)),  &(diag12[0]) );
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
  } 
  else if( diag->id() == -17 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at0, sizeof(cv17at0)/sizeof(int)),  &(diag17[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at1, sizeof(cv17at1)/sizeof(int)),  &(diag17[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at2, sizeof(cv17at2)/sizeof(int)),  &(diag17[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at3, sizeof(cv17at3)/sizeof(int)),  &(diag17[3]) );
  } 
  else if( diag->id() == -18 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at0, sizeof(cv18at0)/sizeof(int)),  &(diag18[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at1, sizeof(cv18at1)/sizeof(int)),  &(diag18[1]) );
  } 
  else if( diag->id() == -19 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int)),  &(diag19[0]) );
  } 
  else if( diag->id() == -20 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int)),  &(diag20[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at1, sizeof(cv20at1)/sizeof(int)),  &(diag20[1]) );
  } 
  else if( diag->id() == -21 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int)),  &(diag21[0]) );
  } 
  else if( diag->id() == -22 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int)),  &(diag22[0]) );
  } 
  else if( diag->id() == -23 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int)),  &(diag23[0]) );
  } 
  else if( diag->id() == -24 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int)),  &(diag24[0]) );
  } 
  else if( diag->id() == -25 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int)),  &(diag25[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at1, sizeof(cv25at1)/sizeof(int)),  &(diag25[1]) );
  } 
  else if( diag->id() == -26 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int)),  &(diag26[0]) );
  } 
  else if( diag->id() == -27 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int)),  &(diag27[0]) );
  } 
  else if( diag->id() == -28 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)),  &(diag28[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at1, sizeof(cv28at1)/sizeof(int)),  &(diag28[1]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at1, sizeof(cv32at1)/sizeof(int)),  &(diag32[1]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv35at1, sizeof(cv35at1)/sizeof(int)),  &(diag35[1]) );
  } 
  else if( diag->id() == -36 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv36at0, sizeof(cv36at0)/sizeof(int)),  &(diag36[0]) );
  } 
  else if( diag->id() == -37 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv37at0, sizeof(cv37at0)/sizeof(int)),  &(diag37[0]) );
  } 
  else if( diag->id() == -38 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int)),  &(diag38[0]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at1, sizeof(cv41at1)/sizeof(int)),  &(diag41[1]) );
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
  } 
  else if( diag->id() == -47 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv47at0, sizeof(cv47at0)/sizeof(int)),  &(diag47[0]) );
  } 
  else if( diag->id() == -48 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at0, sizeof(cv48at0)/sizeof(int)),  &(diag48[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv48at1, sizeof(cv48at1)/sizeof(int)),  &(diag48[1]) );
  } 
  else if( diag->id() == -49 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at0, sizeof(cv49at0)/sizeof(int)),  &(diag49[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv49at1, sizeof(cv49at1)/sizeof(int)),  &(diag49[1]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at2, sizeof(cv52at2)/sizeof(int)),  &(diag52[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv52at3, sizeof(cv52at3)/sizeof(int)),  &(diag52[3]) );
  } 
  else if( diag->id() == -53 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at0, sizeof(cv53at0)/sizeof(int)),  &(diag53[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv53at1, sizeof(cv53at1)/sizeof(int)),  &(diag53[1]) );
  } 
  else if( diag->id() == -54 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv54at0, sizeof(cv54at0)/sizeof(int)),  &(diag54[0]) );
  } 
  else if( diag->id() == -55 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at0, sizeof(cv55at0)/sizeof(int)),  &(diag55[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv55at1, sizeof(cv55at1)/sizeof(int)),  &(diag55[1]) );
  } 
  else if( diag->id() == -56 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv56at0, sizeof(cv56at0)/sizeof(int)),  &(diag56[0]) );
  } 
  else if( diag->id() == -57 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv57at0, sizeof(cv57at0)/sizeof(int)),  &(diag57[0]) );
  } 
  else if( diag->id() == -58 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv58at0, sizeof(cv58at0)/sizeof(int)),  &(diag58[0]) );
  } 
  else if( diag->id() == -59 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv59at0, sizeof(cv59at0)/sizeof(int)),  &(diag59[0]) );
  } 
  else if( diag->id() == -60 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at0, sizeof(cv60at0)/sizeof(int)),  &(diag60[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv60at1, sizeof(cv60at1)/sizeof(int)),  &(diag60[1]) );
  } 
  else if( diag->id() == -61 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv61at0, sizeof(cv61at0)/sizeof(int)),  &(diag61[0]) );
  } 
  else if( diag->id() == -62 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at0, sizeof(cv62at0)/sizeof(int)),  &(diag62[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv62at1, sizeof(cv62at1)/sizeof(int)),  &(diag62[1]) );
  } 
  else if( diag->id() == -63 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv63at0, sizeof(cv63at0)/sizeof(int)),  &(diag63[0]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv67at1, sizeof(cv67at1)/sizeof(int)),  &(diag67[1]) );
  } 
  else if( diag->id() == -68 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at0, sizeof(cv68at0)/sizeof(int)),  &(diag68[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv68at1, sizeof(cv68at1)/sizeof(int)),  &(diag68[1]) );
  } 
  else if( diag->id() == -69 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv69at0, sizeof(cv69at0)/sizeof(int)),  &(diag69[0]) );
  } 
  else if( diag->id() == -70 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at0, sizeof(cv70at0)/sizeof(int)),  &(diag70[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv70at1, sizeof(cv70at1)/sizeof(int)),  &(diag70[1]) );
  } 
  return sel;
}
