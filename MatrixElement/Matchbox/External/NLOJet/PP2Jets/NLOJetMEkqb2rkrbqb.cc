// -*- C++ -*-
//
// NLOJetMEkqb2rkrbqb.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEkqb2rkrbqb class.
//

#include "NLOJetMEkqb2rkrbqb.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEkqb2rkrbqb::NLOJetMEkqb2rkrbqb() {}

NLOJetMEkqb2rkrbqb::~NLOJetMEkqb2rkrbqb() {}

IBPtr NLOJetMEkqb2rkrbqb::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEkqb2rkrbqb::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEkqb2rkrbqb::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEkqb2rkrbqb::persistentOutput(PersistentOStream &) const {}

void NLOJetMEkqb2rkrbqb::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEkqb2rkrbqb,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEkqb2rkrbqb("Herwig::NLOJetMEkqb2rkrbqb", "HwMatchboxNLOJet.so");

void NLOJetMEkqb2rkrbqb::Init() {

  static ClassDocumentation<NLOJetMEkqb2rkrbqb> documentation
    ("NLOJetMEkqb2rkrbqb");

}


void NLOJetMEkqb2rkrbqb::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr k = quark[xi1]; 
    tcPDPtr kb = antiquark[xi1]; 
    for(unsigned xi2 = 0; xi2 != quark.size(); xi2++){ 
      if ( xi1 == xi2 )
	continue;
      tcPDPtr q = quark[xi2]; 
      tcPDPtr qb = antiquark[xi2]; 
      for(unsigned xi3 = 0; xi3 != quark.size(); xi3++){ 
        tcPDPtr r = quark[xi3]; 
        tcPDPtr rb = antiquark[xi3]; 
        addSafe(new_ptr((Tree2toNDiagram(4), k, g, g, qb, 2, g, 5, r, 1, k, 5, rb, 3, qb, -1)));
        addSafe(new_ptr((Tree2toNDiagram(5), k, g, rb, g, qb, 2, r, 1, k, 3, rb, 4, qb, -2)));
        addSafe(new_ptr((Tree2toNDiagram(5), k, g, r, g, qb, 3, r, 1, k, 2, rb, 4, qb, -3)));
        addSafe(new_ptr((Tree2toNDiagram(3), k, g, qb, 2, qb, 4, g, 5, r, 1, k, 5, rb, 4, qb, -4)));
        addSafe(new_ptr((Tree2toNDiagram(4), k, g, q, qb, 3, g, 5, r, 1, k, 5, rb, 2, qb, -5)));
        addSafe(new_ptr((Tree2toNDiagram(3), k, g, qb, 1, k, 4, g, 5, r, 4, k, 5, rb, 2, qb, -6)));
        addSafe(new_ptr((Tree2toNDiagram(4), k, k, g, qb, 1, g, 5, r, 2, k, 5, rb, 3, qb, -7)));
        if( ( xi3 ==  xi2 ) ){ 
          addSafe(new_ptr((Tree2toNDiagram(4), k, g, g, qb, 2, g, 5, q, 1, k, 3, qb, 5, qb, -15)));
          addSafe(new_ptr((Tree2toNDiagram(5), k, g, qb, g, qb, 2, q, 1, k, 4, qb, 3, qb, -16)));
          addSafe(new_ptr((Tree2toNDiagram(5), k, g, q, g, qb, 3, q, 1, k, 4, qb, 2, qb, -17)));
          addSafe(new_ptr((Tree2toNDiagram(3), k, g, qb, 2, qb, 4, g, 5, q, 1, k, 4, qb, 5, qb, -18)));
          addSafe(new_ptr((Tree2toNDiagram(4), k, g, q, qb, 3, g, 5, q, 1, k, 2, qb, 5, qb, -19)));
          addSafe(new_ptr((Tree2toNDiagram(3), k, g, qb, 1, k, 4, g, 5, q, 4, k, 2, qb, 5, qb, -20)));
          addSafe(new_ptr((Tree2toNDiagram(4), k, k, g, qb, 1, g, 5, q, 2, k, 3, qb, 5, qb, -21)));
        }
        if( ( xi3 ==  xi1 ) ){ 
          addSafe(new_ptr((Tree2toNDiagram(4), k, g, g, qb, 2, g, 1, k, 5, k, 5, kb, 3, qb, -22)));
          addSafe(new_ptr((Tree2toNDiagram(5), k, g, kb, g, qb, 1, k, 2, k, 3, kb, 4, qb, -23)));
          addSafe(new_ptr((Tree2toNDiagram(5), k, g, k, g, qb, 1, k, 3, k, 2, kb, 4, qb, -24)));
          addSafe(new_ptr((Tree2toNDiagram(3), k, g, qb, 2, qb, 4, g, 1, k, 5, k, 5, kb, 4, qb, -25)));
          addSafe(new_ptr((Tree2toNDiagram(4), k, g, q, qb, 3, g, 1, k, 5, k, 5, kb, 2, qb, -26)));
          addSafe(new_ptr((Tree2toNDiagram(3), k, g, qb, 1, k, 4, g, 4, k, 5, k, 5, kb, 2, qb, -27)));
          addSafe(new_ptr((Tree2toNDiagram(4), k, k, g, qb, 1, g, 2, k, 5, k, 5, kb, 3, qb, -28)));
        }
      }  
    }  
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEkqb2rkrbqb::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv1at1[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv2at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv3at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv4at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv5at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv6at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv7at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv8at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv8at1[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv9at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv10at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv11at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv12at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv13at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv14at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv15at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv15at1[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv16at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv17at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv18at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv19at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv20at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv21at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv22at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv22at1[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv23at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv24at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv25at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv26at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv27at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv28at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv29at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv29at1[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv30at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv30at1[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv31at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv32at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv33at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv34at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv35at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv36at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv37at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv38at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv39at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv40at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv41at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv42at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};

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
          nloJetAmplitude()->colourOrdered2(cv6at0, sizeof(cv6at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -7 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -8 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -9 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -10 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -11 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -12 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -13 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -14 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -15 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -16 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv16at0, sizeof(cv16at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -17 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv17at0, sizeof(cv17at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -18 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv18at0, sizeof(cv18at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -19 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -20 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -21 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -22 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv22at1, sizeof(cv22at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -23 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -24 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -25 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -26 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -27 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -28 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -29 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv29at1, sizeof(cv29at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -30 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int))
        + nloJetAmplitude()->colourOrdered2(cv30at1, sizeof(cv30at1)/sizeof(int)), i);   
    else if( diags[i]->id() == -31 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -32 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -33 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int)), i);   
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
          nloJetAmplitude()->colourOrdered2(cv38at0, sizeof(cv38at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -39 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -40 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -41 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -42 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEkqb2rkrbqb::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEkqb2rkrbqb

  static const ColourLines diag1[2] = { 
    ColourLines("1 2 5 6, -4 3 -5 -8, 7 -2 -3 -9"), 
    ColourLines("1 2 3 -4, 6 5 -3 -9, 7 -2 -5 -8")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 2 6, -5 4 -8, 7 -2 -3 -4 -9")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 2 3 4 -5, 6 -4 -9, 7 -2 -8")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("1 2 -3, 6 5 -9, 7 -2 -4 -5 -8")
  }; 
  static const ColourLines diag5[1] = { 
    ColourLines("1 2 3 5 6, -4 -5 -8, 7 -2 -9")
  }; 
  static const ColourLines diag6[1] = { 
    ColourLines("1 2 -3, 6 5 4 -2 -9, 7 -5 -8")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("1 5 6, -4 3 2 -5 -8, 7 -3 -9")
  }; 
  static const ColourLines diag8[2] = { 
    ColourLines("1 3 5 7, -2 -3 -4 -8, 6 4 -5 -9"), 
    ColourLines("1 3 4 6, -2 -3 -5 -9, 7 5 -4 -8")
  }; 
  static const ColourLines diag9[1] = { 
    ColourLines("1 3 7, -2 -3 -4 -5 -8, 6 5 -9")
  }; 
  static const ColourLines diag10[1] = { 
    ColourLines("1 3 4 5 6, -2 -3 -9, 7 -5 -8")
  }; 
  static const ColourLines diag11[1] = { 
    ColourLines("1 3 6, -2 -3 -4 -5 -9, 7 5 -8")
  }; 
  static const ColourLines diag12[1] = { 
    ColourLines("1 3 4 5 7, -2 -3 -8, 6 -5 -9")
  }; 
  static const ColourLines diag13[1] = { 
    ColourLines("1 4 6, -3 -5 -9, 7 5 2 -4 -8")
  }; 
  static const ColourLines diag14[1] = { 
    ColourLines("1 5 7, -3 -4 -8, 6 4 2 -5 -9")
  }; 
  static const ColourLines diag15[2] = { 
    ColourLines("1 2 5 6, -4 3 -5 -9, 7 -2 -3 -8"), 
    ColourLines("1 2 3 -4, 6 5 -3 -8, 7 -2 -5 -9")
  }; 
  static const ColourLines diag16[1] = { 
    ColourLines("1 2 6, -5 4 -9, 7 -2 -3 -4 -8")
  }; 
  static const ColourLines diag17[1] = { 
    ColourLines("1 2 3 4 -5, 6 -4 -8, 7 -2 -9")
  }; 
  static const ColourLines diag18[1] = { 
    ColourLines("1 2 -3, 6 5 -8, 7 -2 -4 -5 -9")
  }; 
  static const ColourLines diag19[1] = { 
    ColourLines("1 2 3 5 6, -4 -5 -9, 7 -2 -8")
  }; 
  static const ColourLines diag20[1] = { 
    ColourLines("1 2 -3, 6 5 4 -2 -8, 7 -5 -9")
  }; 
  static const ColourLines diag21[1] = { 
    ColourLines("1 5 6, -4 3 2 -5 -9, 7 -3 -8")
  }; 
  static const ColourLines diag22[2] = { 
    ColourLines("1 2 5 7, -4 3 -5 -8, 6 -2 -3 -9"), 
    ColourLines("1 2 3 -4, 6 -2 -5 -8, 7 5 -3 -9")
  }; 
  static const ColourLines diag23[1] = { 
    ColourLines("1 2 7, -5 4 -8, 6 -2 -3 -4 -9")
  }; 
  static const ColourLines diag24[1] = { 
    ColourLines("1 2 3 4 -5, 6 -2 -8, 7 -4 -9")
  }; 
  static const ColourLines diag25[1] = { 
    ColourLines("1 2 -3, 6 -2 -4 -5 -8, 7 5 -9")
  }; 
  static const ColourLines diag26[1] = { 
    ColourLines("1 2 3 5 7, -4 -5 -8, 6 -2 -9")
  }; 
  static const ColourLines diag27[1] = { 
    ColourLines("1 2 -3, 6 -5 -8, 7 5 4 -2 -9")
  }; 
  static const ColourLines diag28[1] = { 
    ColourLines("1 5 7, -4 3 2 -5 -8, 6 -3 -9")
  }; 
  static const ColourLines diag29[2] = { 
    ColourLines("1 3 5 7, -2 -3 -4 -9, 6 4 -5 -8"), 
    ColourLines("1 3 4 6, -2 -3 -5 -8, 7 5 -4 -9")
  }; 
  static const ColourLines diag30[2] = { 
    ColourLines("1 2 5 7, -4 3 -5 -9, 6 -2 -3 -8"), 
    ColourLines("1 2 3 -4, 6 -2 -5 -9, 7 5 -3 -8")
  }; 
  static const ColourLines diag31[1] = { 
    ColourLines("1 3 7, -2 -3 -4 -5 -9, 6 5 -8")
  }; 
  static const ColourLines diag32[1] = { 
    ColourLines("1 3 4 5 6, -2 -3 -8, 7 -5 -9")
  }; 
  static const ColourLines diag33[1] = { 
    ColourLines("1 3 6, -2 -3 -4 -5 -8, 7 5 -9")
  }; 
  static const ColourLines diag34[1] = { 
    ColourLines("1 3 4 5 7, -2 -3 -9, 6 -5 -8")
  }; 
  static const ColourLines diag35[1] = { 
    ColourLines("1 2 7, -5 4 -9, 6 -2 -3 -4 -8")
  }; 
  static const ColourLines diag36[1] = { 
    ColourLines("1 2 3 4 -5, 6 -2 -9, 7 -4 -8")
  }; 
  static const ColourLines diag37[1] = { 
    ColourLines("1 2 -3, 6 -2 -4 -5 -9, 7 5 -8")
  }; 
  static const ColourLines diag38[1] = { 
    ColourLines("1 2 3 5 7, -4 -5 -9, 6 -2 -8")
  }; 
  static const ColourLines diag39[1] = { 
    ColourLines("1 2 -3, 6 -5 -9, 7 5 4 -2 -8")
  }; 
  static const ColourLines diag40[1] = { 
    ColourLines("1 5 7, -4 3 2 -5 -9, 6 -3 -8")
  }; 
  static const ColourLines diag41[1] = { 
    ColourLines("1 4 6, -3 -5 -8, 7 5 2 -4 -9")
  }; 
  static const ColourLines diag42[1] = { 
    ColourLines("1 5 7, -3 -4 -9, 6 4 2 -5 -8")
  }; 

  static const int cv1at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv1at1[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv2at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv3at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv4at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv5at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv6at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv7at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv8at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv8at1[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv9at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv10at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv11at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv12at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv13at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv14at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv15at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv15at1[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv16at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv17at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv18at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv19at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv20at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv21at0[]  = { 1, -1, -999, 0, 4, -999, 2, 3, -999};
  static const int cv22at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv22at1[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv23at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv24at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv25at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv26at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv27at0[]  = { 0, -1, -999, 1, 3, -999, 2, 4, -999};
  static const int cv28at0[]  = { 2, -1, -999, 0, 3, -999, 1, 4, -999};
  static const int cv29at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv29at1[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv30at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv30at1[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv31at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv32at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv33at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv34at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv35at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv36at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv37at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv38at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv39at0[]  = { 0, -1, -999, 1, 4, -999, 2, 3, -999};
  static const int cv40at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};
  static const int cv41at0[]  = { 1, -1, -999, 0, 3, -999, 2, 4, -999};
  static const int cv42at0[]  = { 2, -1, -999, 0, 4, -999, 1, 3, -999};

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
  } 
  else if( diag->id() == -7 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv7at0, sizeof(cv7at0)/sizeof(int)),  &(diag7[0]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at0, sizeof(cv8at0)/sizeof(int)),  &(diag8[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv8at1, sizeof(cv8at1)/sizeof(int)),  &(diag8[1]) );
  } 
  else if( diag->id() == -9 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv9at0, sizeof(cv9at0)/sizeof(int)),  &(diag9[0]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv10at0, sizeof(cv10at0)/sizeof(int)),  &(diag10[0]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv11at0, sizeof(cv11at0)/sizeof(int)),  &(diag11[0]) );
  } 
  else if( diag->id() == -12 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv12at0, sizeof(cv12at0)/sizeof(int)),  &(diag12[0]) );
  } 
  else if( diag->id() == -13 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv13at0, sizeof(cv13at0)/sizeof(int)),  &(diag13[0]) );
  } 
  else if( diag->id() == -14 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv14at0, sizeof(cv14at0)/sizeof(int)),  &(diag14[0]) );
  } 
  else if( diag->id() == -15 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at0, sizeof(cv15at0)/sizeof(int)),  &(diag15[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv15at1, sizeof(cv15at1)/sizeof(int)),  &(diag15[1]) );
  } 
  else if( diag->id() == -16 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv16at0, sizeof(cv16at0)/sizeof(int)),  &(diag16[0]) );
  } 
  else if( diag->id() == -17 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv17at0, sizeof(cv17at0)/sizeof(int)),  &(diag17[0]) );
  } 
  else if( diag->id() == -18 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv18at0, sizeof(cv18at0)/sizeof(int)),  &(diag18[0]) );
  } 
  else if( diag->id() == -19 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv19at0, sizeof(cv19at0)/sizeof(int)),  &(diag19[0]) );
  } 
  else if( diag->id() == -20 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv20at0, sizeof(cv20at0)/sizeof(int)),  &(diag20[0]) );
  } 
  else if( diag->id() == -21 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv21at0, sizeof(cv21at0)/sizeof(int)),  &(diag21[0]) );
  } 
  else if( diag->id() == -22 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at0, sizeof(cv22at0)/sizeof(int)),  &(diag22[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv22at1, sizeof(cv22at1)/sizeof(int)),  &(diag22[1]) );
  } 
  else if( diag->id() == -23 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv23at0, sizeof(cv23at0)/sizeof(int)),  &(diag23[0]) );
  } 
  else if( diag->id() == -24 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv24at0, sizeof(cv24at0)/sizeof(int)),  &(diag24[0]) );
  } 
  else if( diag->id() == -25 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv25at0, sizeof(cv25at0)/sizeof(int)),  &(diag25[0]) );
  } 
  else if( diag->id() == -26 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv26at0, sizeof(cv26at0)/sizeof(int)),  &(diag26[0]) );
  } 
  else if( diag->id() == -27 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv27at0, sizeof(cv27at0)/sizeof(int)),  &(diag27[0]) );
  } 
  else if( diag->id() == -28 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv28at0, sizeof(cv28at0)/sizeof(int)),  &(diag28[0]) );
  } 
  else if( diag->id() == -29 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at0, sizeof(cv29at0)/sizeof(int)),  &(diag29[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv29at1, sizeof(cv29at1)/sizeof(int)),  &(diag29[1]) );
  } 
  else if( diag->id() == -30 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at0, sizeof(cv30at0)/sizeof(int)),  &(diag30[0]) );
   sel.insert( nloJetAmplitude()->colourOrdered2(cv30at1, sizeof(cv30at1)/sizeof(int)),  &(diag30[1]) );
  } 
  else if( diag->id() == -31 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv31at0, sizeof(cv31at0)/sizeof(int)),  &(diag31[0]) );
  } 
  else if( diag->id() == -32 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv32at0, sizeof(cv32at0)/sizeof(int)),  &(diag32[0]) );
  } 
  else if( diag->id() == -33 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv33at0, sizeof(cv33at0)/sizeof(int)),  &(diag33[0]) );
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
  } 
  else if( diag->id() == -39 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv39at0, sizeof(cv39at0)/sizeof(int)),  &(diag39[0]) );
  } 
  else if( diag->id() == -40 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv40at0, sizeof(cv40at0)/sizeof(int)),  &(diag40[0]) );
  } 
  else if( diag->id() == -41 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv41at0, sizeof(cv41at0)/sizeof(int)),  &(diag41[0]) );
  } 
  else if( diag->id() == -42 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv42at0, sizeof(cv42at0)/sizeof(int)),  &(diag42[0]) );
  } 
  return sel;
}
