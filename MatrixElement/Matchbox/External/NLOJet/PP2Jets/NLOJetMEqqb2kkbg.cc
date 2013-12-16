// -*- C++ -*-
//
// NLOJetMEqqb2kkbg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEqqb2kkbg class.
//

#include "NLOJetMEqqb2kkbg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEqqb2kkbg::NLOJetMEqqb2kkbg() {}

NLOJetMEqqb2kkbg::~NLOJetMEqqb2kkbg() {}

IBPtr NLOJetMEqqb2kkbg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEqqb2kkbg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEqqb2kkbg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEqqb2kkbg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEqqb2kkbg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEqqb2kkbg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEqqb2kkbg("Herwig::NLOJetMEqqb2kkbg", "HwMatchboxNLOJet.so");

void NLOJetMEqqb2kkbg::Init() {

  static ClassDocumentation<NLOJetMEqqb2kkbg> documentation
    ("NLOJetMEqqb2kkbg");

}


void NLOJetMEqqb2kkbg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr k = quark[xi1]; 
    tcPDPtr kb = antiquark[xi1]; 
    for(unsigned xi2 = 0; xi2 != quark.size(); xi2++){ 
      tcPDPtr q = quark[xi2]; 
      tcPDPtr qb = antiquark[xi2]; 
      addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 4, k, 4, kb, 3, g, -1)));
      addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, k, 4, k, 3, kb, 4, g, -2)));
      addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, kb, 3, k, 4, kb, 4, g, -3)));
      addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 4, k, 4, kb, 1, g, -4)));
      addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 4, k, 4, kb, 2, g, -5)));
      if( ( xi1 ==  xi2 ) ){ 
        addSafe(new_ptr((Tree2toNDiagram(4), q, g, g, qb, 1, q, 3, qb, 2, g, -6)));
        addSafe(new_ptr((Tree2toNDiagram(4), q, g, q, qb, 1, q, 2, qb, 3, g, -7)));
        addSafe(new_ptr((Tree2toNDiagram(3), q, g, qb, 2, qb, 1, q, 4, qb, 4, g, -8)));
        addSafe(new_ptr((Tree2toNDiagram(4), q, q, g, qb, 2, q, 3, qb, 1, g, -9)));
        addSafe(new_ptr((Tree2toNDiagram(3), q, g, qb, 1, q, 4, q, 2, qb, 4, g, -10)));
      }
    }  
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEqqb2kkbg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 0, 2, -999, 1, 3, -1, -999};
  static const int cv1at1[]  = { 1, -1, -999, 0, 3, 2, -999};
  static const int cv2at0[]  = { 0, 2, -999, 1, 3, -1, -999};
  static const int cv3at0[]  = { 1, -1, -999, 0, 3, 2, -999};
  static const int cv4at0[]  = { 0, 2, -999, 1, 3, -1, -999};
  static const int cv5at0[]  = { 1, -1, -999, 0, 3, 2, -999};
  static const int cv6at0[]  = { 0, 3, -1, -999, 1, 2, -999};
  static const int cv6at1[]  = { 0, -1, -999, 1, 3, 2, -999};
  static const int cv7at0[]  = { 0, 3, -1, -999, 1, 2, -999};
  static const int cv8at0[]  = { 0, -1, -999, 1, 3, 2, -999};
  static const int cv9at0[]  = { 0, 3, -1, -999, 1, 2, -999};
  static const int cv10at0[]  = { 0, -1, -999, 1, 3, 2, -999};

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
    else if( diags[i]->id() == -4 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int)), i);   
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
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEqqb2kkbg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEqqb2kkbg

  static const ColourLines diag1[2] = { 
    ColourLines("1 3 7, -2 -3 -4 -6, 5 4 -7"), 
    ColourLines("1 3 4 5, -2 -3 -7, -6 -4 7")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 3 4 7, -2 -3 -6, 5 -7")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 3 5, -2 -3 -4 -7, -6 7")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("1 7, -3 -4 -6, 5 4 2 -7")
  }; 
  static const ColourLines diag5[1] = { 
    ColourLines("1 4 5, -3 -7, -6 -4 2 7")
  }; 
  static const ColourLines diag6[2] = { 
    ColourLines("1 2 7, -4 3 -7, 5 -2 -3 -6"), 
    ColourLines("1 2 3 -4, 5 -2 -7, -6 -3 7")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("1 2 3 7, -4 -7, 5 -2 -6")
  }; 
  static const ColourLines diag8[1] = { 
    ColourLines("1 2 -3, 5 -2 -4 -7, -6 7")
  }; 
  static const ColourLines diag9[1] = { 
    ColourLines("1 7, -4 3 2 -7, 5 -3 -6")
  }; 
  static const ColourLines diag10[1] = { 
    ColourLines("1 2 -3, 5 -7, -6 -2 4 7")
  }; 

  static const int cv1at0[]  = { 0, 2, -999, 1, 3, -1, -999};
  static const int cv1at1[]  = { 1, -1, -999, 0, 3, 2, -999};
  static const int cv2at0[]  = { 0, 2, -999, 1, 3, -1, -999};
  static const int cv3at0[]  = { 1, -1, -999, 0, 3, 2, -999};
  static const int cv4at0[]  = { 0, 2, -999, 1, 3, -1, -999};
  static const int cv5at0[]  = { 1, -1, -999, 0, 3, 2, -999};
  static const int cv6at0[]  = { 0, 3, -1, -999, 1, 2, -999};
  static const int cv6at1[]  = { 0, -1, -999, 1, 3, 2, -999};
  static const int cv7at0[]  = { 0, 3, -1, -999, 1, 2, -999};
  static const int cv8at0[]  = { 0, -1, -999, 1, 3, 2, -999};
  static const int cv9at0[]  = { 0, 3, -1, -999, 1, 2, -999};
  static const int cv10at0[]  = { 0, -1, -999, 1, 3, 2, -999};

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
  else if( diag->id() == -4 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv4at0, sizeof(cv4at0)/sizeof(int)),  &(diag4[0]) );
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
  return sel;
}
