// -*- C++ -*-
//
// NLOJetMEgq2kqkb.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEgq2kqkb class.
//

#include "NLOJetMEgq2kqkb.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEgq2kqkb::NLOJetMEgq2kqkb() {}

NLOJetMEgq2kqkb::~NLOJetMEgq2kqkb() {}

IBPtr NLOJetMEgq2kqkb::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEgq2kqkb::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEgq2kqkb::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEgq2kqkb::persistentOutput(PersistentOStream &) const {}

void NLOJetMEgq2kqkb::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEgq2kqkb,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEgq2kqkb("Herwig::NLOJetMEgq2kqkb", "HwMatchboxNLOJet.so");

void NLOJetMEgq2kqkb::Init() {

  static ClassDocumentation<NLOJetMEgq2kqkb> documentation
    ("NLOJetMEgq2kqkb");

}


void NLOJetMEgq2kqkb::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr k = quark[xi1]; 
    tcPDPtr kb = antiquark[xi1]; 
    for(unsigned xi2 = 0; xi2 != quark.size(); xi2++){ 
      tcPDPtr q = quark[xi2]; 
      tcPDPtr qb = antiquark[xi2]; 
      addSafe(new_ptr((Tree2toNDiagram(2), g, q, 1, q, 3, g, 4, k, 3, q, 4, kb, -1)));
      addSafe(new_ptr((Tree2toNDiagram(4), g, kb, g, q, 1, k, 3, q, 2, kb, -2)));
      addSafe(new_ptr((Tree2toNDiagram(3), g, qb, q, 2, g, 4, k, 1, q, 4, kb, -3)));
      addSafe(new_ptr((Tree2toNDiagram(4), g, k, g, q, 2, k, 3, q, 1, kb, -4)));
      addSafe(new_ptr((Tree2toNDiagram(3), g, g, q, 1, g, 4, k, 2, q, 4, kb, -5)));
      if( ( xi1 ==  xi2 ) ){ 
        addSafe(new_ptr((Tree2toNDiagram(2), g, q, 1, q, 3, g, 3, q, 4, q, 4, qb, -6)));
        addSafe(new_ptr((Tree2toNDiagram(3), g, qb, q, 2, g, 1, q, 4, q, 4, qb, -7)));
        addSafe(new_ptr((Tree2toNDiagram(4), g, qb, g, q, 3, q, 1, q, 2, qb, -8)));
        addSafe(new_ptr((Tree2toNDiagram(4), g, q, g, q, 3, q, 2, q, 1, qb, -9)));
        addSafe(new_ptr((Tree2toNDiagram(3), g, g, q, 1, g, 2, q, 4, q, 4, qb, -10)));
      }
    }  
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEgq2kqkb::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 1, -1, 0, -999, 2, 3, -999};
  static const int cv2at0[]  = { 1, -1, 0, -999, 2, 3, -999};
  static const int cv3at0[]  = { 2, -1, 3, -999, 1, 0, -999};
  static const int cv4at0[]  = { 2, -1, 3, -999, 1, 0, -999};
  static const int cv5at0[]  = { 1, -1, 0, -999, 2, 3, -999};
  static const int cv5at1[]  = { 2, -1, 3, -999, 1, 0, -999};
  static const int cv6at0[]  = { 2, -1, 0, -999, 1, 3, -999};
  static const int cv7at0[]  = { 1, -1, 3, -999, 2, 0, -999};
  static const int cv8at0[]  = { 2, -1, 0, -999, 1, 3, -999};
  static const int cv9at0[]  = { 1, -1, 3, -999, 2, 0, -999};
  static const int cv10at0[]  = { 2, -1, 0, -999, 1, 3, -999};
  static const int cv10at1[]  = { 1, -1, 3, -999, 2, 0, -999};

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
          nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -5 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -6 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int)), i);   
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
          nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEgq2kqkb::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEgq2kqkb

  static const ColourLines diag1[1] = { 
    ColourLines("1 3 4 5, -1 2, 6 -4 -7")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 5, -1 -2 -3 4, 6 3 -7")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 6, -1 -2 -4 -7, 3 4 5")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("1 2 3 6, -1 -7, 4 -3 5")
  }; 
  static const ColourLines diag5[2] = { 
    ColourLines("1 4 5, -1 -2 3, 6 2 -4 -7"), 
    ColourLines("1 2 6, -1 -4 -7, 3 -2 4 5")
  }; 
  static const ColourLines diag6[1] = { 
    ColourLines("1 3 4 6, -1 2, 5 -4 -7")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("1 5, -1 -2 -4 -7, 3 4 6")
  }; 
  static const ColourLines diag8[1] = { 
    ColourLines("1 6, -1 -2 -3 4, 5 3 -7")
  }; 
  static const ColourLines diag9[1] = { 
    ColourLines("1 2 3 5, -1 -7, 4 -3 6")
  }; 
  static const ColourLines diag10[2] = { 
    ColourLines("1 4 6, -1 -2 3, 5 2 -4 -7"), 
    ColourLines("1 2 5, -1 -4 -7, 3 -2 4 6")
  }; 

  static const int cv1at0[]  = { 1, -1, 0, -999, 2, 3, -999};
  static const int cv2at0[]  = { 1, -1, 0, -999, 2, 3, -999};
  static const int cv3at0[]  = { 2, -1, 3, -999, 1, 0, -999};
  static const int cv4at0[]  = { 2, -1, 3, -999, 1, 0, -999};
  static const int cv5at0[]  = { 1, -1, 0, -999, 2, 3, -999};
  static const int cv5at1[]  = { 2, -1, 3, -999, 1, 0, -999};
  static const int cv6at0[]  = { 2, -1, 0, -999, 1, 3, -999};
  static const int cv7at0[]  = { 1, -1, 3, -999, 2, 0, -999};
  static const int cv8at0[]  = { 2, -1, 0, -999, 1, 3, -999};
  static const int cv9at0[]  = { 1, -1, 3, -999, 2, 0, -999};
  static const int cv10at0[]  = { 2, -1, 0, -999, 1, 3, -999};
  static const int cv10at1[]  = { 1, -1, 3, -999, 2, 0, -999};

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
  } 
  else if( diag->id() == -5 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at0, sizeof(cv5at0)/sizeof(int)),  &(diag5[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv5at1, sizeof(cv5at1)/sizeof(int)),  &(diag5[1]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int)),  &(diag6[0]) );
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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at1, sizeof(cv10at1)/sizeof(int)),  &(diag10[1]) );
  } 
  return sel;
}
