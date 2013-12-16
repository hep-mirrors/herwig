// -*- C++ -*-
//
// NLOJetMEqbq2gg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEqbq2gg class.
//

#include "NLOJetMEqbq2gg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEqbq2gg::NLOJetMEqbq2gg() {}

NLOJetMEqbq2gg::~NLOJetMEqbq2gg() {}

IBPtr NLOJetMEqbq2gg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEqbq2gg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEqbq2gg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEqbq2gg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEqbq2gg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEqbq2gg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEqbq2gg("Herwig::NLOJetMEqbq2gg", "HwMatchboxNLOJet.so");

void NLOJetMEqbq2gg::Init() {

  static ClassDocumentation<NLOJetMEqbq2gg> documentation
    ("NLOJetMEqbq2gg");

}


void NLOJetMEqbq2gg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr q = quark[xi1]; 
    tcPDPtr qb = antiquark[xi1]; 
    addSafe(new_ptr((Tree2toNDiagram(2), qb, q, 1, g, 3, g, 3, g, -1)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 1, g, 2, g, -2)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 2, g, 1, g, -3)));
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEqbq2gg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 1, 2, 0, -999};
  static const int cv1at1[]  = { -1, 2, 1, 0, -999};
  static const int cv2at0[]  = { -1, 1, 2, 0, -999};
  static const int cv3at0[]  = { -1, 2, 1, 0, -999};

  Selector <MEBase::DiagramIndex> sel;
  for(MEBase::DiagramIndex i=0; i < diags.size(); ++i){
    if( diags[i]->id() == -1 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEqbq2gg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEqbq2gg

  static const ColourLines diag1[2] = { 
    ColourLines("-1 -3 -4, 2 3 5, 4 -5"), 
    ColourLines("-1 -3 -5, 2 3 4, -4 5")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("-1 -4, 3 5, 4 -2 -5")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("-1 -5, 3 4, -4 -2 5")
  }; 

  static const int cv1at0[]  = { -1, 1, 2, 0, -999};
  static const int cv1at1[]  = { -1, 2, 1, 0, -999};
  static const int cv2at0[]  = { -1, 1, 2, 0, -999};
  static const int cv3at0[]  = { -1, 2, 1, 0, -999};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int)),  &(diag1[1]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
  } 
  return sel;
}
