// -*- C++ -*-
//
// NLOJetMEkbq2qkb.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEkbq2qkb class.
//

#include "NLOJetMEkbq2qkb.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEkbq2qkb::NLOJetMEkbq2qkb() {}

NLOJetMEkbq2qkb::~NLOJetMEkbq2qkb() {}

IBPtr NLOJetMEkbq2qkb::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEkbq2qkb::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEkbq2qkb::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEkbq2qkb::persistentOutput(PersistentOStream &) const {}

void NLOJetMEkbq2qkb::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEkbq2qkb,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEkbq2qkb("Herwig::NLOJetMEkbq2qkb", "HwMatchboxNLOJet.so");

void NLOJetMEkbq2qkb::Init() {

  static ClassDocumentation<NLOJetMEkbq2qkb> documentation
    ("NLOJetMEkbq2qkb");

}


void NLOJetMEkbq2qkb::doGetDiagrams() const {
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
      addSafe(new_ptr((Tree2toNDiagram(3), kb, g, q, 2, q, 1, kb, -1)));
    }  
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEkbq2qkb::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { -1, 0, -999, 1, 2, -999};
  static const int cv2at0[]  = { -1, 2, -999, 1, 0, -999};

  Selector <MEBase::DiagramIndex> sel;
  for(MEBase::DiagramIndex i=0; i < diags.size(); ++i){
    if( diags[i]->id() == -1 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)), i);   
    else if( diags[i]->id() == -2 )
      sel.insert( 
          nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)), i);   
  }
  return sel;
}


Selector<const ColourLines *> NLOJetMEkbq2qkb::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEkbq2qkb

  static const ColourLines diag1[1] = { 
    ColourLines("-1 -2 3, 4 2 -5")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("-1 -3 -5, 2 3 4")
  }; 

  static const int cv1at0[]  = { -1, 0, -999, 1, 2, -999};
  static const int cv2at0[]  = { -1, 2, -999, 1, 0, -999};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
  } 
  return sel;
}
