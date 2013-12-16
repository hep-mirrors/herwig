// -*- C++ -*-
//
// NLOJetMEqbg2qbg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEqbg2qbg class.
//

#include "NLOJetMEqbg2qbg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEqbg2qbg::NLOJetMEqbg2qbg() {}

NLOJetMEqbg2qbg::~NLOJetMEqbg2qbg() {}

IBPtr NLOJetMEqbg2qbg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEqbg2qbg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEqbg2qbg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEqbg2qbg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEqbg2qbg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEqbg2qbg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEqbg2qbg("Herwig::NLOJetMEqbg2qbg", "HwMatchboxNLOJet.so");

void NLOJetMEqbg2qbg::Init() {

  static ClassDocumentation<NLOJetMEqbg2qbg> documentation
    ("NLOJetMEqbg2qbg");

}


void NLOJetMEqbg2qbg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr q = quark[xi1]; 
    tcPDPtr qb = antiquark[xi1]; 
    addSafe(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb, 3, qb, 3, g, -1)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, g, g, 1, qb, 2, g, -2)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, g, 2, qb, 1, g, -3)));
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEqbg2qbg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 0, 2, 1, -999};
  static const int cv2at0[]  = { -1, 0, 2, 1, -999};
  static const int cv2at1[]  = { -1, 2, 0, 1, -999};
  static const int cv3at0[]  = { -1, 2, 0, 1, -999};

  Selector <MEBase::DiagramIndex> sel;
  for(MEBase::DiagramIndex i=0; i < diags.size(); ++i){
    if( diags[i]->id() == -1 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEqbg2qbg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEqbg2qbg

  static const ColourLines diag1[1] = { 
    ColourLines("-1 2, -2 -3 -5, -4 5")
  }; 
  static const ColourLines diag2[2] = { 
    ColourLines("-1 -2 3, -3 -5, -4 2 5"), 
    ColourLines("-1 -2 -5, 3 5, -3 2 -4")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("-1 -5, 3 -2 5, -3 -4")
  }; 

  static const int cv1at0[]  = { -1, 0, 2, 1, -999};
  static const int cv2at0[]  = { -1, 0, 2, 1, -999};
  static const int cv2at1[]  = { -1, 2, 0, 1, -999};
  static const int cv3at0[]  = { -1, 2, 0, 1, -999};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int)),  &(diag2[1]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
  } 
  return sel;
}
