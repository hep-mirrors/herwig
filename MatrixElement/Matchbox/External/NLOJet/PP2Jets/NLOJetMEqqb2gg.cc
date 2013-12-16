// -*- C++ -*-
//
// NLOJetMEqqb2gg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEqqb2gg class.
//

#include "NLOJetMEqqb2gg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEqqb2gg::NLOJetMEqqb2gg() {}

NLOJetMEqqb2gg::~NLOJetMEqqb2gg() {}

IBPtr NLOJetMEqqb2gg::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEqqb2gg::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEqqb2gg::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEqqb2gg::persistentOutput(PersistentOStream &) const {}

void NLOJetMEqqb2gg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEqqb2gg,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEqqb2gg("Herwig::NLOJetMEqqb2gg", "HwMatchboxNLOJet.so");

void NLOJetMEqqb2gg::Init() {

  static ClassDocumentation<NLOJetMEqqb2gg> documentation
    ("NLOJetMEqqb2gg");

}


void NLOJetMEqqb2gg::doGetDiagrams() const {
  // get the particle data objects

  PDPtr g = getParticleData(ParticleID::g); 
  for(unsigned xi1 = 0; xi1 != quark.size(); xi1++){ 
    tcPDPtr q = quark[xi1]; 
    tcPDPtr qb = antiquark[xi1]; 
    addSafe(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, g, 3, g, -1)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 2, g, -2)));
    addSafe(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 1, g, -3)));
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEqqb2gg::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 0, 1, 2, -1, -999};
  static const int cv1at1[]  = { 0, 2, 1, -1, -999};
  static const int cv2at0[]  = { 0, 2, 1, -1, -999};
  static const int cv3at0[]  = { 0, 1, 2, -1, -999};

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


Selector<const ColourLines *> NLOJetMEqqb2gg::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEqqb2gg

  static const ColourLines diag1[2] = { 
    ColourLines("1 3 5, -2 -3 -4, 4 -5"), 
    ColourLines("1 3 4, -2 -3 -5, -4 5")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 4, -3 -5, -4 2 5")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 5, -3 -4, 4 2 -5")
  }; 

  static const int cv1at0[]  = { 0, 1, 2, -1, -999};
  static const int cv1at1[]  = { 0, 2, 1, -1, -999};
  static const int cv2at0[]  = { 0, 2, 1, -1, -999};
  static const int cv3at0[]  = { 0, 1, 2, -1, -999};

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
