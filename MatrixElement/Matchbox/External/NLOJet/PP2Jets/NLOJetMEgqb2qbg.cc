// -*- C++ -*-
//
// NLOJetMEgqb2qbg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEgqb2qbg class.
//

#include "NLOJetMEgqb2qbg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEgqb2qbg::NLOJetMEgqb2qbg() {}

NLOJetMEgqb2qbg::~NLOJetMEgqb2qbg() {}

IBPtr NLOJetMEgqb2qbg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEgqb2qbg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEgqb2qbg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEgqb2qbg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEgqb2qbg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEgqb2qbg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEgqb2qbg("Herwig::NLOJetMEgqb2qbg", "HwMatchboxNLOJet.so");

void NLOJetMEgqb2qbg::Init() {

  static ClassDocumentation<NLOJetMEgqb2qbg> documentation
    ("NLOJetMEgqb2qbg");

}


void NLOJetMEgqb2qbg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr q = quark[xi1]; 
    tcPDPtr qb = antiquark[xi1]; 
    addSafe(new_ptr((Tree2toNDiagram(2), g, qb, 1, qb, 3, qb, 3, g, -1)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, q, qb, 1, qb, 2, g, -2)));
    addSafe(new_ptr((Tree2toNDiagram(3), g, g, qb, 2, qb, 1, g, -3)));
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEgqb2qbg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 0, -1, 2, 1, -999};
  static const int cv2at0[]  = { 0, 2, -1, 1, -999};
  static const int cv3at0[]  = { 0, -1, 2, 1, -999};
  static const int cv3at1[]  = { 0, 2, -1, 1, -999};

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
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEgqb2qbg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEgqb2qbg

  static const ColourLines diag1[1] = { 
    ColourLines("1 -2, -1 -3 -5, -4 5")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 2 5, -1 -4, -3 -5")
  }; 
  static const ColourLines diag3[2] = { 
    ColourLines("1 2 -3, -1 -5, -4 -2 5"), 
    ColourLines("1 5, -1 -2 -4, -3 2 -5")
  }; 

  static const int cv1at0[]  = { 0, -1, 2, 1, -999};
  static const int cv2at0[]  = { 0, 2, -1, 1, -999};
  static const int cv3at0[]  = { 0, -1, 2, 1, -999};
  static const int cv3at1[]  = { 0, 2, -1, 1, -999};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int)),  &(diag3[1]) );
  } 
  return sel;
}
