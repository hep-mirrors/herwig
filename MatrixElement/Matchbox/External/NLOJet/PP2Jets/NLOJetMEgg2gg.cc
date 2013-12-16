// -*- C++ -*-
//
// NLOJetMEgg2gg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEgg2gg class.
//

#include "NLOJetMEgg2gg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEgg2gg::NLOJetMEgg2gg() {}

NLOJetMEgg2gg::~NLOJetMEgg2gg() {}

IBPtr NLOJetMEgg2gg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEgg2gg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEgg2gg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEgg2gg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEgg2gg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEgg2gg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEgg2gg("Herwig::NLOJetMEgg2gg", "HwMatchboxNLOJet.so");

void NLOJetMEgg2gg::Init() {

  static ClassDocumentation<NLOJetMEgg2gg> documentation
    ("NLOJetMEgg2gg");

}


void NLOJetMEgg2gg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  addSafe(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 3, g, -1)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 2, g, -2)));
  addSafe(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 1, g, -3)));
}


Selector <MEBase::DiagramIndex> NLOJetMEgg2gg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 2, 1, 0, -998};
  static const int cv1at1[]  = { -1, 0, 2, 1, -998};
  static const int cv1at2[]  = { -1, 1, 2, 0, -998};
  static const int cv1at3[]  = { -1, 0, 1, 2, -998};
  static const int cv2at0[]  = { -1, 2, 0, 1, -998};
  static const int cv2at1[]  = { -1, 1, 2, 0, -998};
  static const int cv2at2[]  = { -1, 0, 2, 1, -998};
  static const int cv2at3[]  = { -1, 1, 0, 2, -998};
  static const int cv3at0[]  = { -1, 1, 0, 2, -998};
  static const int cv3at1[]  = { -1, 2, 1, 0, -998};
  static const int cv3at2[]  = { -1, 0, 1, 2, -998};
  static const int cv3at3[]  = { -1, 2, 0, 1, -998};

  Selector <MEBase::DiagramIndex> sel;
  for(MEBase::DiagramIndex i=0; i < diags.size(); ++i){
    if( diags[i]->id() == -1 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at2, sizeof(cv1at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv1at3, sizeof(cv1at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at2, sizeof(cv2at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv2at3, sizeof(cv2at3)/sizeof(int)), i);   
    else if( diags[i]->id() == -3 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at2, sizeof(cv3at2)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv3at3, sizeof(cv3at3)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEgg2gg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEgg2gg

  static const ColourLines diag1[4] = { 
    ColourLines("1 3 5, -1 2, -2 -3 -4, 4 -5"), 
    ColourLines("1 -2, -1 -3 -4, 2 3 5, 4 -5"), 
    ColourLines("1 3 4, -1 2, -2 -3 -5, -4 5"), 
    ColourLines("1 -2, -1 -3 -5, 2 3 4, -4 5")
  }; 
  static const ColourLines diag2[4] = { 
    ColourLines("1 2 5, -1 -4, 3 -2 4, -3 -5"), 
    ColourLines("1 4, -1 -2 3, -3 -5, -4 2 5"), 
    ColourLines("1 2 -3, -1 -4, 3 5, 4 -2 -5"), 
    ColourLines("1 4, -1 -2 -5, 3 5, -3 2 -4")
  }; 
  static const ColourLines diag3[4] = { 
    ColourLines("1 2 4, -1 -5, 3 -2 5, -3 -4"), 
    ColourLines("1 5, -1 -2 3, -3 -4, 4 2 -5"), 
    ColourLines("1 2 -3, -1 -5, 3 4, -4 -2 5"), 
    ColourLines("1 5, -1 -2 -4, 3 4, -3 2 -5")
  }; 

  static const int cv1at0[]  = { -1, 2, 1, 0, -998};
  static const int cv1at1[]  = { -1, 0, 2, 1, -998};
  static const int cv1at2[]  = { -1, 1, 2, 0, -998};
  static const int cv1at3[]  = { -1, 0, 1, 2, -998};
  static const int cv2at0[]  = { -1, 2, 0, 1, -998};
  static const int cv2at1[]  = { -1, 1, 2, 0, -998};
  static const int cv2at2[]  = { -1, 0, 2, 1, -998};
  static const int cv2at3[]  = { -1, 1, 0, 2, -998};
  static const int cv3at0[]  = { -1, 1, 0, 2, -998};
  static const int cv3at1[]  = { -1, 2, 1, 0, -998};
  static const int cv3at2[]  = { -1, 0, 1, 2, -998};
  static const int cv3at3[]  = { -1, 2, 0, 1, -998};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at1, sizeof(cv1at1)/sizeof(int)),  &(diag1[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at2, sizeof(cv1at2)/sizeof(int)),  &(diag1[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at3, sizeof(cv1at3)/sizeof(int)),  &(diag1[3]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at1, sizeof(cv2at1)/sizeof(int)),  &(diag2[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at2, sizeof(cv2at2)/sizeof(int)),  &(diag2[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at3, sizeof(cv2at3)/sizeof(int)),  &(diag2[3]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at0, sizeof(cv3at0)/sizeof(int)),  &(diag3[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at1, sizeof(cv3at1)/sizeof(int)),  &(diag3[1]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at2, sizeof(cv3at2)/sizeof(int)),  &(diag3[2]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv3at3, sizeof(cv3at3)/sizeof(int)),  &(diag3[3]) );
  } 
  return sel;
}
