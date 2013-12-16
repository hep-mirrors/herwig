// -*- C++ -*-
//
// NLOJetMEkqb2kqb.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NLOJetMEkqb2kqb class.
//

#include "NLOJetMEkqb2kqb.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/External/NLOJet/NLOJetPhasespace.h"

using namespace Herwig;

NLOJetMEkqb2kqb::NLOJetMEkqb2kqb() {}

NLOJetMEkqb2kqb::~NLOJetMEkqb2kqb() {}

IBPtr NLOJetMEkqb2kqb::clone() const {
  return new_ptr(*this);
}

IBPtr NLOJetMEkqb2kqb::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void NLOJetMEkqb2kqb::doinit() {
  NLOJetMEBase<0,2,0>::doinit();
}

void NLOJetMEkqb2kqb::persistentOutput(PersistentOStream &) const {}

void NLOJetMEkqb2kqb::persistentInput(PersistentIStream &, int) {}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Herwig::NLOJetMEkqb2kqb,Herwig::NLOJetMEBase<0,2,0> >
  describeNLOJetMEkqb2kqb("Herwig::NLOJetMEkqb2kqb", "HwMatchboxNLOJet.so");

void NLOJetMEkqb2kqb::Init() {

  static ClassDocumentation<NLOJetMEkqb2kqb> documentation
    ("NLOJetMEkqb2kqb");

}


void NLOJetMEkqb2kqb::doGetDiagrams() const {
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
      addSafe(new_ptr((Tree2toNDiagram(3), k, g, qb, 1, k, 2, qb, -1)));
    }  
  }  
}


Selector <MEBase::DiagramIndex> NLOJetMEkqb2kqb::diagrams(const DiagramVector & diags) const {
  // select the diagram
  matchboxAmplitude()->prepareAmplitudes(this);

  static const int cv1at0[]  = { 0, -1, -999, 1, 2, -999};
  static const int cv2at0[]  = { 1, -1, -999, 0, 2, -999};

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


Selector<const ColourLines *> NLOJetMEkqb2kqb::colourGeometries(tcDiagPtr diag) const {
  // colour lines for NLOJetMEkqb2kqb

  static const ColourLines diag1[1] = { 
    ColourLines("1 2 -3, 4 -2 -5")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 3 4, -2 -3 -5")
  }; 

  static const int cv1at0[]  = { 0, -1, -999, 1, 2, -999};
  static const int cv2at0[]  = { 1, -1, -999, 0, 2, -999};

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv1at0, sizeof(cv1at0)/sizeof(int)),  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( nloJetAmplitude()->colourOrdered2(cv2at0, sizeof(cv2at0)/sizeof(int)),  &(diag2[0]) );
  } 
  return sel;
}
