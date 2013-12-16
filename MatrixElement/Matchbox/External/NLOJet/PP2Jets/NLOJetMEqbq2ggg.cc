// -*- C++ -*-
//
// NLOJetMEqbq2ggg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEqbq2ggg class.
//

#include "NLOJetMEqbq2ggg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEqbq2ggg::NLOJetMEqbq2ggg() {}

NLOJetMEqbq2ggg::~NLOJetMEqbq2ggg() {}

IBPtr NLOJetMEqbq2ggg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEqbq2ggg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEqbq2ggg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEqbq2ggg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEqbq2ggg::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEqbq2ggg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEqbq2ggg("Herwig::NLOJetMEqbq2ggg", "HwMatchboxNLOJet.so");

void NLOJetMEqbq2ggg::Init() {

  static ClassDocumentation<NLOJetMEqbq2ggg> documentation
    ("NLOJetMEqbq2ggg");

}


void NLOJetMEqbq2ggg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr q = quark[xi1]; 
    tcPDPtr qb = antiquark[xi1]; 
    addSafe(new_ptr((Tree2toNDiagram(2), qb, q, 1, g, 3, g, 4, g, 4, g, 3, g, -1)));
    addSafe(new_ptr((Tree2toNDiagram(2), qb, q, 1, g, 3, g, 4, g, 3, g, 4, g, -2)));
    addSafe(new_ptr((Tree2toNDiagram(2), qb, q, 1, g, 3, g, 3, g, 4, g, 4, g, -3)));
    addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, q, 1, g, 3, g, 2, g, -4)));
    addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, q, 1, g, 2, g, 3, g, -5)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 2, g, 1, g, 4, g, 4, g, -6)));
    addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, q, 3, g, 1, g, 2, g, -7)));
    addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, q, 2, g, 1, g, 3, g, -8)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 2, g, 4, g, 1, g, 4, g, -9)));
    addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, q, 3, g, 2, g, 1, g, -10)));
    addSafe(new_ptr((Tree2toNDiagram(4), qb, qb, qb, q, 2, g, 3, g, 1, g, -11)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 2, g, 4, g, 4, g, 1, g, -12)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 1, g, 2, g, 4, g, 4, g, -13)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 1, g, 4, g, 2, g, 4, g, -14)));
    addSafe(new_ptr((Tree2toNDiagram(3), qb, qb, q, 1, g, 4, g, 4, g, 2, g, -15)));
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEqbq2ggg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv1at1[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv1at2[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv1at3[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv2at0[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv2at1[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv2at2[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv2at3[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv3at0[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv3at1[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv3at2[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv3at3[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv4at0[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv5at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv6at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv6at1[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv7at0[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv8at0[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv9at0[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv9at1[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv10at0[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv11at0[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv12at0[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv12at1[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv13at0[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv13at1[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv14at0[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv14at1[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv15at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv15at1[]  = { -1, 2, 1, 3, 0, -999};

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
          nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -10 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -11 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)), i);   
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
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEqbq2ggg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEqbq2ggg

  static const ColourLines diag1[4] = { 
    ColourLines("-1 -3 -4 -5, 2 3 7, 5 -6, 6 4 -7"), 
    ColourLines("-1 -3 -4 -6, 2 3 7, 5 4 -7, -5 6"), 
    ColourLines("-1 -3 -7, 2 3 4 6, 5 -6, -5 -4 7"), 
    ColourLines("-1 -3 -7, 2 3 4 5, -5 6, -6 -4 7")
  }; 
  static const ColourLines diag2[4] = { 
    ColourLines("-1 -3 -4 -5, 2 3 6, 5 -7, -6 4 7"), 
    ColourLines("-1 -3 -4 -7, 2 3 6, 5 4 -6, -5 7"), 
    ColourLines("-1 -3 -6, 2 3 4 7, 5 -7, -5 -4 6"), 
    ColourLines("-1 -3 -6, 2 3 4 5, -5 7, 6 -4 -7")
  }; 
  static const ColourLines diag3[4] = { 
    ColourLines("-1 -3 -4 -6, 2 3 5, -5 4 7, 6 -7"), 
    ColourLines("-1 -3 -4 -7, 2 3 5, -5 4 6, -6 7"), 
    ColourLines("-1 -3 -5, 2 3 4 7, 5 -4 -6, 6 -7"), 
    ColourLines("-1 -3 -5, 2 3 4 6, 5 -4 -7, -6 7")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("-1 -5, 4 6, 5 -2 -7, -6 -3 7")
  }; 
  static const ColourLines diag5[1] = { 
    ColourLines("-1 -5, 4 7, 5 -2 -6, 6 -3 -7")
  }; 
  static const ColourLines diag6[2] = { 
    ColourLines("-1 -5, 3 4 7, 5 -2 -4 -6, 6 -7"), 
    ColourLines("-1 -5, 3 4 6, 5 -2 -4 -7, -6 7")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("-1 -6, 4 5, -5 -3 7, 6 -2 -7")
  }; 
  static const ColourLines diag8[1] = { 
    ColourLines("-1 -6, 4 7, 5 -3 -7, -5 -2 6")
  }; 
  static const ColourLines diag9[2] = { 
    ColourLines("-1 -6, 3 4 7, 5 -7, -5 -4 -2 6"), 
    ColourLines("-1 -6, 3 4 5, -5 7, 6 -2 -4 -7")
  }; 
  static const ColourLines diag10[1] = { 
    ColourLines("-1 -7, 4 5, -5 -3 6, -6 -2 7")
  }; 
  static const ColourLines diag11[1] = { 
    ColourLines("-1 -7, 4 6, 5 -3 -6, -5 -2 7")
  }; 
  static const ColourLines diag12[2] = { 
    ColourLines("-1 -7, 3 4 6, 5 -6, -5 -4 -2 7"), 
    ColourLines("-1 -7, 3 4 5, -5 6, -6 -4 -2 7")
  }; 
  static const ColourLines diag13[2] = { 
    ColourLines("-1 -4 -6, 3 5, -5 -2 4 7, 6 -7"), 
    ColourLines("-1 -4 -7, 3 5, -5 -2 4 6, -6 7")
  }; 
  static const ColourLines diag14[2] = { 
    ColourLines("-1 -4 -5, 3 6, 5 -7, -6 -2 4 7"), 
    ColourLines("-1 -4 -7, 3 6, 5 4 -2 -6, -5 7")
  }; 
  static const ColourLines diag15[2] = { 
    ColourLines("-1 -4 -5, 3 7, 5 -6, 6 4 -2 -7"), 
    ColourLines("-1 -4 -6, 3 7, 5 4 -2 -7, -5 6")
  }; 

  static const int cv1at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv1at1[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv1at2[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv1at3[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv2at0[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv2at1[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv2at2[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv2at3[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv3at0[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv3at1[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv3at2[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv3at3[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv4at0[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv5at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv6at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv6at1[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv7at0[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv8at0[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv9at0[]  = { -1, 2, 1, 3, 0, -999};
  static const int cv9at1[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv10at0[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv11at0[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv12at0[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv12at1[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv13at0[]  = { -1, 2, 3, 1, 0, -999};
  static const int cv13at1[]  = { -1, 3, 2, 1, 0, -999};
  static const int cv14at0[]  = { -1, 1, 3, 2, 0, -999};
  static const int cv14at1[]  = { -1, 3, 1, 2, 0, -999};
  static const int cv15at0[]  = { -1, 1, 2, 3, 0, -999};
  static const int cv15at1[]  = { -1, 2, 1, 3, 0, -999};

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
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at1, sizeof(cv9at1)/sizeof(int)),  &(diag9[1]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)),  &(diag10[0]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)),  &(diag11[0]) );
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
  return sel;
}
